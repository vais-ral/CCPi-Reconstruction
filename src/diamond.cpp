
#include <iostream>
#include "base_types.hpp"
#include "instruments.hpp"

bool CCPi::Diamond::setup_experimental_geometry(const std::string path,
						const std::string file,
						const bool phantom)
{
  std::cerr << "Todo\n";
  return false;
}

bool CCPi::Diamond::read_scans(const std::string path, const bool phantom)
{
  std::cerr << "Todo\n";
  return false;
}

bool CCPi::Diamond::finish_voxel_geometry(real voxel_origin[3],
					  real voxel_size[3],
					  const voxel_data &voxels) const
{
  std::cerr << "Todo\n";
  return false;
}

void CCPi::Diamond::apply_beam_hardening()
{
  std::cerr << "Todo\n";
}

long CCPi::Diamond::get_data_size() const
{
  std::cerr << "Todo\n";
  return 0;
}

pixel_type *const CCPi::Diamond::get_pixel_data() const
{
  std::cerr << "Todo\n";
  return 0;
}

void CCPi::Diamond::forward_project(pixel_type *pixels,
				    voxel_type *const voxels,
				    const real origin[3],
				    const real width[3], const int nx,
				    const int ny, const int nz) const
{
  std::cerr << "Todo\n";
}

void CCPi::Diamond::backward_project(pixel_type *pixels,
				     voxel_type *const voxels,
				     const real origin[3],
				     const real width[3], const int nx,
				     const int ny, const int nz) const
{
  std::cerr << "Todo\n";
}

void CCPi::Diamond::backward_project(voxel_type *const voxels,
				     const real origin[3],
				     const real width[3], const int nx,
				     const int ny, const int nz) const
{
  std::cerr << "Todo\n";
}
