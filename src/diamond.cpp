
#include <iostream>
#include "base_types.hpp"
#include "instruments.hpp"

bool CCPi::Diamond::setup_experimental_geometry(const std::string file,
						const bool phantom)
{
  std::cerr << "Todo\n";
  return false;
}

bool CCPi::Diamond::read_scans(const bool phantom)
{
  std::cerr << "Todo\n";
  return false;
}

void CCPi::Diamond::forward_project(voxel_data &voxels, const real origin[3],
				    const real width[3])
{
  std::cerr << "Todo\n";
}

void CCPi::Diamond::backward_project(voxel_data &voxels, const real origin[3],
				     const real width[3])
{
  std::cerr << "Todo\n";
}
