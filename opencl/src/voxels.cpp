
#include "base_types.hpp"
#include "utils.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "results.hpp"
#include "voxels.hpp"
#include "blas.hpp"
#include "mpi.hpp"
#include "ui_calls.hpp"

void calculate_block_sizes(int &nx_voxels, int &ny_voxels, int &nz_voxels,
			   int &maxz_voxels, int &block_size, int &block_step,
			   const int num_processors, const int blocking_factor,
			   const int pixels_per_voxel,
			   CCPi::instrument *instrument,
			   const bool recon_blocks)
{
  // calculate blocks
  instrument->get_xy_size(nx_voxels, ny_voxels, pixels_per_voxel);
  maxz_voxels = instrument->get_z_size(std::max(nx_voxels, ny_voxels),
				       pixels_per_voxel);
  add_output("voxel sizes ");
  add_output(nx_voxels);
  add_output(" ");
  add_output(ny_voxels);
  add_output(" ");
  add_output(maxz_voxels);
  send_output();
  nz_voxels = 0;
  block_size = 0;
  block_step = 0;
  if (blocking_factor == 0 and num_processors == 1) {
    nz_voxels = maxz_voxels;
    block_size = nz_voxels;
    block_step = nz_voxels;
  } else if (instrument->supports_blocks() and recon_blocks) {
    int sz = 1;
    if (blocking_factor > 0)
      sz = blocking_factor;
    if (maxz_voxels / (sz * num_processors) < 1)
      report_error("Reduce blocking factor or number of processors");
    block_size = sz;
    block_step = block_size * num_processors;
    nz_voxels = block_size;
  } else if (num_processors == 1) {
    report_error("Ignoring blocking factor - not supported by device");
    nz_voxels = maxz_voxels;
    block_size = nz_voxels;
    block_step = nz_voxels;
  }
}

voxel_data *reconstruct(CCPi::instrument *device,
			CCPi::reconstruction_alg *algorithm,
			const std::string data_file,
			const std::string output_name, real full_vox_origin[3],
			real voxel_size[3], const real rotation_centre,
			const int pixels_per_voxel, const int blocking_factor,
			const bool beam_harden,
			const CCPi::output_format write_format,
			const bool clamp_output, const bool phantom)
{
  int num_processors = machine::get_number_of_processors();
  voxel_data *voxels = 0;
  std::string path;
  std::string filename;
  CCPi::split_path_and_name(data_file, path, filename);
  if (device->setup_experimental_geometry(path, filename, rotation_centre,
					  pixels_per_voxel, phantom)) {
    int nx_voxels = 0;
    int ny_voxels = 0;
    int maxz_voxels = 0;
    int nz_voxels = 0;
    int block_size = 0;
    int block_step = 0;
    calculate_block_sizes(nx_voxels, ny_voxels, nz_voxels, maxz_voxels,
			  block_size, block_step, num_processors,
			  blocking_factor, pixels_per_voxel,
			  device, algorithm->supports_blocks());
    int z_data_size = block_size * pixels_per_voxel;
    int z_data_step = block_step * pixels_per_voxel;
    device->set_v_block(z_data_size);
    int block_offset = machine::get_processor_id() * block_size;
    int z_data_offset = block_offset * pixels_per_voxel;
    //real full_vox_origin[3];
    //real voxel_size[3];
    if (device->finish_voxel_geometry(full_vox_origin, voxel_size,
				      nx_voxels, ny_voxels, maxz_voxels)) {
      // can modify offsets and end if parallel beam to solve subregion
      int end_value = device->total_num_v_pixels();
      bool ok = false;
      bool first = true;
      do {
	ok = false;
	if (block_offset + block_size > maxz_voxels)
	  block_size = maxz_voxels - block_offset;
	if (z_data_offset + z_data_size > end_value) {
	  z_data_size = end_value - z_data_offset;
	  device->set_v_block(z_data_size);
	}
	nz_voxels = block_size;
	real voxel_origin[3];
	voxel_origin[0] = full_vox_origin[0];
	voxel_origin[1] = full_vox_origin[1];
	voxel_origin[2] = full_vox_origin[2]
	  + block_offset * voxel_size[2];
	if (device->read_scans(path, z_data_offset,
			       z_data_size, first, phantom)) {
	  voxels =
	    new voxel_data(boost::extents[nx_voxels][ny_voxels][nz_voxels]);
	  if (beam_harden)
	    device->apply_beam_hardening();
	  ok = algorithm->reconstruct(device, *voxels,
				      voxel_origin, voxel_size);
	  if (ok) {
	    clamp_min(*voxels, 0.0, nx_voxels, ny_voxels, nz_voxels);
	    CCPi::write_results(output_name, *voxels, full_vox_origin,
				voxel_size, block_offset, maxz_voxels,
				write_format, clamp_output);
	  }
	} else
	  ok = false;
	first = false;
	block_offset += block_step;
	z_data_offset += z_data_step;
      } while (ok and z_data_offset < end_value);
    }
  }
  return voxels;
}

// Todo - alot of this is copied from the above one, which isn't ideal
// but simplify since we will reconstruct all the pixel data provided as
// a single block probably. At least with DLS SAVU
voxel_data *reconstruct(CCPi::instrument *device,
			CCPi::reconstruction_alg *algorithm,
			const numpy_3d &pixels, const numpy_1d &angles,
			const real rotation_centre, const int pixels_per_voxel,
			const int blocking_factor, const bool beam_harden)
{
  int num_processors = machine::get_number_of_processors();
  voxel_data *voxels = 0;
  if (device->setup_experimental_geometry(pixels, angles, rotation_centre,
					  pixels_per_voxel)) {
    int nx_voxels = 0;
    int ny_voxels = 0;
    int maxz_voxels = 0;
    int nz_voxels = 0;
    int block_size = 0;
    int block_step = 0;
    calculate_block_sizes(nx_voxels, ny_voxels, nz_voxels, maxz_voxels,
			  block_size, block_step, num_processors,
			  blocking_factor, pixels_per_voxel,
			  device, algorithm->supports_blocks());
    int z_data_size = block_size * pixels_per_voxel;
    int z_data_step = block_step * pixels_per_voxel;
    device->set_v_block(z_data_size);
    int block_offset = machine::get_processor_id() * block_size;
    int z_data_offset = block_offset * pixels_per_voxel;
    real full_vox_origin[3];
    real voxel_size[3];
    if (device->finish_voxel_geometry(full_vox_origin, voxel_size,
				      nx_voxels, ny_voxels, maxz_voxels)) {
      // can modify offsets and end if parallel beam to solve subregion
      int end_value = device->total_num_v_pixels();
      bool ok = false;
      //bool first = true;
      do {
	ok = false;
	if (block_offset + block_size > maxz_voxels)
	  block_size = maxz_voxels - block_offset;
	if (z_data_offset + z_data_size > end_value) {
	  z_data_size = end_value - z_data_offset;
	  device->set_v_block(z_data_size);
	}
	nz_voxels = block_size;
	real voxel_origin[3];
	voxel_origin[0] = full_vox_origin[0];
	voxel_origin[1] = full_vox_origin[1];
	voxel_origin[2] = full_vox_origin[2]
	  + block_offset * voxel_size[2];
	if (device->read_scans(pixels, z_data_offset, z_data_size)) {
	  voxels =
	    new voxel_data(boost::extents[nx_voxels][ny_voxels][nz_voxels]);
	  if (beam_harden)
	    device->apply_beam_hardening();
	  ok = algorithm->reconstruct(device, *voxels,
				      voxel_origin, voxel_size);
	  if (ok) {
	    // truncate negative values
	    clamp_min(*voxels, 0.0, nx_voxels, ny_voxels, nz_voxels);
	  }
	} else
	  ok = false;
	//first = false;
	block_offset += block_step;
	z_data_offset += z_data_step;
      } while (ok and z_data_offset < end_value);
    }
  }
  /*
  voxel_data *voxels = new voxel_data(boost::extents[2][2][2]);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	(*voxels)[i][j][k] = i + j + k;
  */
  return voxels;
}

// Todo - another copy of the savu one - can we merge in some way?
voxel_data *reconstruct(CCPi::instrument *device,
			CCPi::reconstruction_alg *algorithm,
			const numpy_3d &pixels, const numpy_1d &angles,
			const numpy_1d &h_offsets, const numpy_1d &v_offsets,
			const int pixels_per_voxel, const real source_x,
			const real detector_x, const real pixel_h_size,
			const real pixel_v_size, const real mask_radius,
			const bool beam_harden, real full_vox_origin[3],
			real voxel_size[3], const bool has_offsets)
{
  int num_processors = machine::get_number_of_processors();
  const int blocking_factor = 0;
  voxel_data *voxels = 0;
  if (device->setup_experimental_geometry(pixels, angles, h_offsets, v_offsets,
					  pixels_per_voxel, source_x,
					  detector_x, pixel_h_size,
					  pixel_v_size, mask_radius,
					  has_offsets)) {
    int nx_voxels = 0;
    int ny_voxels = 0;
    int maxz_voxels = 0;
    int nz_voxels = 0;
    int block_size = 0;
    int block_step = 0;
    calculate_block_sizes(nx_voxels, ny_voxels, nz_voxels, maxz_voxels,
			  block_size, block_step, num_processors,
			  blocking_factor, pixels_per_voxel,
			  device, algorithm->supports_blocks());
    int z_data_size = block_size * pixels_per_voxel;
    int z_data_step = block_step * pixels_per_voxel;
    device->set_v_block(z_data_size);
    int block_offset = machine::get_processor_id() * block_size;
    int z_data_offset = block_offset * pixels_per_voxel;
    //real full_vox_origin[3];
    //real voxel_size[3];
    if (device->finish_voxel_geometry(full_vox_origin, voxel_size,
				      nx_voxels, ny_voxels, maxz_voxels)) {
      // can modify offsets and end if parallel beam to solve subregion
      int end_value = device->total_num_v_pixels();
      bool ok = false;
      //bool first = true;
      do {
	ok = false;
	if (block_offset + block_size > maxz_voxels)
	  block_size = maxz_voxels - block_offset;
	if (z_data_offset + z_data_size > end_value) {
	  z_data_size = end_value - z_data_offset;
	  device->set_v_block(z_data_size);
	}
	nz_voxels = block_size;
	real voxel_origin[3];
	voxel_origin[0] = full_vox_origin[0];
	voxel_origin[1] = full_vox_origin[1];
	voxel_origin[2] = full_vox_origin[2]
	  + block_offset * voxel_size[2];
	if (device->read_scans(pixels, z_data_offset, z_data_size)) {
	  voxels =
	    new voxel_data(boost::extents[nx_voxels][ny_voxels][nz_voxels]);
	  if (beam_harden)
	    device->apply_beam_hardening();
	  ok = algorithm->reconstruct(device, *voxels,
				      voxel_origin, voxel_size);
	  if (ok) {
	    // truncate negative values
	    clamp_min(*voxels, 0.0, nx_voxels, ny_voxels, nz_voxels);
	  }
	} else
	  ok = false;
	//first = false;
	block_offset += block_step;
	z_data_offset += z_data_step;
      } while (ok and z_data_offset < end_value);
    }
  }
  /*
  voxel_data *voxels = new voxel_data(boost::extents[2][2][2]);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	(*voxels)[i][j][k] = i + j + k;
  */
  return voxels;
}
