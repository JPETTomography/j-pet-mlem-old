#ifndef RECONSTRUCTOR
#define RECONSTRUCTOR

#include <vector>
#include <algorithm>

#include "2d/barrel/lor_info.h"
#include "2d/strip/event.h"
#include "3d/geometry/point.h"
#include "3d/geometry/voxel_grid.h"

#include "util/grapher.h"

template <typename Scanner, typename Kernel2D> class Reconstructor {
 public:
  using F = typename Scanner::F;
  using S = typename Scanner::S;
  using LorPixelInfo = PET2D::Barrel::LorPixelnfo<F, S>;
  using Response = typename Scanner::Response;
  using LOR = PET2D::Barrel::LOR<S>;
  using StripEvent = PET2D::Strip::Event<F>;
  using PixelInfo = typename LorPixelInfo::PixelInfo;
  using Pixel = typename LorPixelInfo::Pixel;
  using Point2D = PET2D::Point<F>;
  using Point = PET3D::Point<F>;
  using Vector2D = PET2D::Vector<F>;

  struct FrameEvent {
    LOR lor;
    F up;
    F right;
    F tan;
    F sec;
    F half_box_up;
    F half_box_right;
    typename LorPixelInfo::PixelInfoContainer::const_iterator first_pixel;
    typename LorPixelInfo::PixelInfoContainer::const_iterator last_pixel;
    S first_plane;
    S last_plane;
    Response response;
  };

  Reconstructor(const Scanner& scanner,
                const LorPixelInfo& lor_pixel_info,
                F z_left,
                int n_planes)
      : scanner_(scanner),
        lor_pixel_info_(lor_pixel_info),
        z_left(z_left),
        n_planes(n_planes),
        v_grid(lor_pixel_info_.grid, z_left, n_planes),
        n_voxels(lor_pixel_info_.grid.n_columns * lor_pixel_info_.grid.n_rows *
                 n_planes),
        kernel_(scanner.sigma_z(), scanner.sigma_dl()),
        rho_new_(n_voxels, F(0.0)),
        rho_(n_voxels, F(1.0)),
        kernel_cache_(n_voxels, F(0.0)) {}

  Point translate_to_point(const Response& response) {

    auto segment = lor_pixel_info_[response.lor].segment;
    F t = 0.5 - response.dl / (2 * segment->length);
    return PET3D::interpolate(
        t,
        Point(segment->start.x, segment->start.y, response.z_dn),
        Point(segment->end.x, segment->end.y, response.z_up));
  }

  FrameEvent translate_to_frame(const Response& response) {
    FrameEvent event;
    event.response = response;
    event.lor = response.lor;

    auto R = lor_pixel_info_[event.lor].segment->length / 2;
    StripEvent strip_event(response.z_up, response.z_dn, response.dl);

    strip_event.transform(R, event.tan, event.up, event.right);
    F A, B, C;
    kernel_.ellipse_bb(
        event.tan, event.sec, A, B, C, event.half_box_up, event.half_box_right);

    auto ev_z_left = event.right - event.half_box_right;
    auto ev_z_right = event.right + event.half_box_right;
    event.first_plane = plane(ev_z_left);
    event.last_plane = plane(ev_z_right) + 1;

    auto y_up = event.up + event.half_box_up;
    auto y_dn = event.up - event.half_box_up;
    auto t_up = (y_up + R) / (2 * R);
    auto t_dn = (y_dn + R) / (2 * R);

    PixelInfo pix_info_up, pix_info_dn;
    pix_info_up.t = t_up;
    pix_info_dn.t = t_dn;
    event.last_pixel =
        std::upper_bound(lor_pixel_info_[event.lor].pixels.begin(),
                         lor_pixel_info_[event.lor].pixels.end(),
                         pix_info_up,
                         [](const PixelInfo& a, const PixelInfo& b)
                             -> bool { return a.t < b.t; });

    event.first_pixel =
        std::lower_bound(lor_pixel_info_[event.lor].pixels.begin(),
                         lor_pixel_info_[event.lor].pixels.end(),
                         pix_info_dn,
                         [](const PixelInfo& a, const PixelInfo& b)
                             -> bool { return a.t < b.t; });

    return event;
  }

  S plane(F z) {
    return S(std::floor((z - z_left) / lor_pixel_info_.grid.pixel_size));
  }

  int fscanf_responses(std::istream& in) {
    int count = 0;
    auto response = fscanf_response(in);
    while (in) {
      count++;
      events_.push_back(translate_to_frame(response));
      response = fscanf_response(in);
    }
    return count;
  }

  int n_events() const { return events_.size(); }
  FrameEvent frame_event(int i) const { return events_[i]; }

  int iterate() {
    event_count_ = 0;
    voxel_count_ = 0;
    pixel_count_ = 0;
    auto grid = lor_pixel_info_.grid;
    std::fill(rho_new_.begin(), rho_new_.end(), 0);
    int i;

    /* ------- Event loop ------*/
    for (i = 0; i < n_events(); ++i) {
      event_count_++;
      auto event = frame_event(i);
      auto lor = event.lor;
      auto segment = *lor_pixel_info_[lor].segment;
      auto R = segment.length / 2;
      auto width = lor_pixel_info_[lor].width;

      //      std::cout << "event " << event.up << " " << event.right << " "
      //                << event.tan << "\n";

      /* ---------  Voxel loop  - denominator ----------- */
      double denominator = 0;
      for (auto it = event.first_pixel; it != event.last_pixel; ++it) {
        pixel_count_++;
        auto pix = it->pixel;
        auto ix = pix.x;
        auto iy = pix.y;

        auto center = grid.center_at(ix, iy);
        auto distance = segment.distance_from(center);
        auto up = segment.projection_relative_middle(center);

        for (int iz = event.first_plane; iz < event.last_plane; ++iz) {
          voxel_count_++;
          auto z = v_grid.center_z_at(ix, iy, iz);
          int index = v_grid.index(ix, iy, iz);

          auto diff = Point2D(up, z) - Point2D(event.up, event.right);
          auto weight =
              kernel_(
                  event.up, event.tan, event.sec, R, Vector2D(diff.y, diff.x)) *
              rho_[index];
          kernel_cache_[index] = weight;
          denominator += weight;
        }
      }  // Voxel loop - denominator

      F inv_denominator;
      if (denominator > 0) {
        inv_denominator = 1 / denominator;
      } else {
        std::cerr << "denminator == 0 !";
        return i;
      }

      /* ---------  Voxel loop ------------ */
      for (auto it = event.first_pixel; it != event.last_pixel; ++it) {
        auto pix = it->pixel;
        auto ix = pix.x;
        auto iy = pix.y;
        for (int iz = event.first_plane; iz < event.last_plane; ++iz) {
          int index = v_grid.index(ix, iy, iz);
          rho_new_[index] += kernel_cache_[index] * inv_denominator;
        }
      }  // Voxel loop
    }
    std::swap(rho_, rho_new_);
    return i;
  }

  typename std::vector<F>::const_iterator rho_begin() const {
    return rho_.begin();
  }
  typename std::vector<F>::const_iterator rho_end() const { return rho_.end(); }
  typename std::vector<F>::iterator rho_begin() { return rho_.begin(); }

  int voxel_count() const { return voxel_count_; }
  int pixel_count() const { return pixel_count_; }
  int event_count() const { return event_count_; }

  void graph_frame_event(Graphics<F>& graphics, int event_index) {
    auto event = events_[event_index];
    auto lor = event.lor;
    graphics.add(scanner_.barrel, lor);
    graphics.add(*lor_pixel_info_[lor].segment);
    for (auto pix = event.first_pixel; pix != event.last_pixel; ++pix) {
      graphics.addPixel(lor_pixel_info_.grid, pix->pixel);
    }

    auto p = translate_to_point(event.response);

    graphics.add(Point2D(p.x, p.y));
  }

 private:
  Response fscanf_response(std::istream& in) {
    S d1, d2;
    F z_up, z_dn, dl;
    Response response;
    in >> d1 >> d2 >> response.z_up >> response.z_dn >> response.dl;
    response.lor = LOR(d1, d2);
    return response;
  }

 public:
  const Scanner& scanner_;
  const LorPixelInfo& lor_pixel_info_;
  const F z_left;
  const S n_planes;
  const PET3D::VoxelGrid<F, S> v_grid;
  const int n_voxels;

 private:
  std::vector<FrameEvent> events_;
  Kernel2D kernel_;
  std::vector<F> rho_new_;
  std::vector<F> rho_;
  int event_count_;
  int voxel_count_;
  int pixel_count_;
  std::vector<F> kernel_cache_;
};

#endif  // RECONSTRUCTOR
