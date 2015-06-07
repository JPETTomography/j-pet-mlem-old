#ifndef VOXEL_SET_BUILDER
#define VOXEL_SET_BUILDER

#include "3d/geometry/voxel_set.h"

namespace PET3D {
template <typename FType, typename SType> class VoxelSetBuilder {
  using F = FType;
  using S = SType;
  using VoxelSet = PET3D::VoxelSet<F,S>;

  public:
  /* Builds a traingular cut of a Z slice at the given z-plane
   * assuming that the pixel grid is centered on (0,0).
   *
   */
  static void BuildTriagularZSlice(VoxelSet& set, S iz, F fov_radius) {
    auto& v_grid = set.v_grid;
    auto& p_grid = v_grid.pixel_grid;
    auto cix=p_grid.n_columns/2;
    auto ciy=p_grid.n_rows/2;
    for(S ix = cix; ix<p_grid.n_columns;ix++)
          for(S iy = ciy; iy<=ix;iy++) {
              auto p = p_grid.center_at(ix,iy);
              if(p.x*p.x+p.y*p.y<=fov_radius*fov_radius) {
                  set.push_back(Voxel<S>(ix, iy, iz));
              }
          }

  }

  static void BuildYSlice(VoxelSet& set, S iy, F fov_radius) {
    auto& v_grid = set.v_grid;
    auto& p_grid = v_grid.pixel_grid;
    auto cix=p_grid.n_columns/2;
    auto ciz=v_grid.n_planes/2;
    for(S ix = cix; ix<p_grid.n_columns;ix++)
          for(S iz = ciz; iz<=v_grid.n_planes;iz++) {
              auto p = p_grid.center_at(ix, iy);
              if(p.x*p.x+p.y*p.y<=fov_radius*fov_radius) {
                  set.push_back(Voxel<S>(ix, iy, iz));
              }
          }
  }
};


}
#endif  // VOXEL_SET_BUILDER
