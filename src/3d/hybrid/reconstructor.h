#ifndef RECONSTRUCTOR
#define RECONSTRUCTOR

#include <vector>
#include <algorithm>

#include "2d/barrel/lor_info.h"
#include "2d/strip/event.h"

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
  };

  Reconstructor(const Scanner& scanner,
                const LorPixelInfo& lor_pixel_info,
                F fov_length,
                F z_left)
      : scanner_(scanner),
        lor_pixel_info_(lor_pixel_info),
        fov_length(fov_length),
        z_left(z_left),
        n_planes((int)std::ceil(fov_length / lor_pixel_info_.grid.pixel_size)),
        kernel_(scanner.sigma_z(), scanner.sigma_dl()),
        rho_new_(lor_pixel_info_.grid.n_columns * lor_pixel_info_.grid.n_rows *
                     n_planes,
                 F(0.0)),
        rho_(lor_pixel_info_.grid.n_columns * lor_pixel_info_.grid.n_rows *
                 n_planes,
             F(1.0)) {}

  FrameEvent translate_to_frame(const Response& response) {
    FrameEvent event;
    event.lor = response.lor;

    auto R = lor_pixel_info_[event.lor].segment.length;
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

  Response fscanf_responses(std::istream& in) {
    auto response = fscanf_response(in);
    while (in) {
      events_.push_back(translate_to_frame(response));
      response = fscanf_response(in);
    }
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
      auto segment = lor_pixel_info_[lor].segment;
      auto R = segment.length / 2;
      auto width = lor_pixel_info_[lor].width;

      /* ---------  Voxel loop  - denominator ----------- */
      F denominator = 0;
      for (auto it = event.first_pixel; it != event.last_pixel; ++it) {
        pixel_count_++;
        auto pix = it->pixel;
        auto ix = pix.x;
        auto iy = pix.y;
        auto center = grid.center_at(ix, iy);
        auto distance = segment.distance_from(center);
        auto up = segment.projection_relative_middle(center);

        for (int plane = event.first_plane; plane < event.last_plane; ++plane) {
          voxel_count_++;
          auto z = (plane + 0.5) * grid.pixel_size + z_left;
          int index =
              ix + iy * grid.n_columns + plane * grid.n_columns * grid.n_rows;
          auto weight =
              kernel_(event.up,
                      event.tan,
                      event.sec,
                      R,
                      Point2D(up, z) - Point2D(event.up, event.right));
          std::cerr << event.up << " " << event.right << " " << up << " " << z
                    << " ";
          std::cerr << weight << " " << rho_[index] << "\n";
          denominator += weight * rho_[index];
        }
      }  // Voxel loop - denominator

      F inv_denominator;
      if (denominator > 0)
        inv_denominator = 1 / denominator;
      else {
        std::cerr << "denminator == 0 !";
        return i;
      }

      /* ---------  Voxel loop ------------ */
      for (auto it = event.first_pixel; it != event.last_pixel; ++it) {
        auto pix = it->pixel;
        auto ix = pix.x;
        auto iy = pix.y;
        auto center = grid.center_at(ix, iy);
        auto distance = segment.distance_from(center);
        auto up = segment.projection_relative_middle(center);

        for (int plane = event.first_plane; plane < event.last_plane; ++plane) {
          auto z = (plane + 0.5) * grid.pixel_size - z_left;
          int index =
              ix + iy * grid.n_columns + plane * grid.n_columns * grid.n_rows;
          auto weight =
              kernel_(event.up,
                      event.tan,
                      event.sec,
                      R,
                      Point2D(up, z) - Point2D(event.up, event.right));
          rho_new_[index] += weight * rho_[index] * inv_denominator;
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

  int n_voxels() const {
    return lor_pixel_info_.grid.n_columns * lor_pixel_info_.grid.n_rows *
           n_planes;
  }

  void graph_frame_event(std::ostream& out, int event_index) {}

 private:
  Response fscanf_response(std::istream& in) {
    S d1, d2;
    F z_up, z_dn, dl;
    Response response;
    in >> d1 >> d2 >> response.z_up >> response.z_dn >> response.dl;
    response.lor = LOR(d1, d2);
    return response;
  }

  const Scanner& scanner_;
  const LorPixelInfo& lor_pixel_info_;
  const F fov_length;
  const F z_left;
  const S n_planes;
  std::vector<FrameEvent> events_;
  Kernel2D kernel_;
  std::vector<F> rho_new_;
  std::vector<F> rho_;
  int event_count_;
  int voxel_count_;
  int pixel_count_;
};

#endif  // RECONSTRUCTOR
