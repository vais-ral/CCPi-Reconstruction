#include "ReconstructionImpl.h"
#include "omp.h"
#include "cgls.hpp"
#include "sirt.hpp"
#include "mlem.hpp"
#include "ui_calls.hpp"
#include "utils.hpp"
#include "voxels.hpp"
#include "blas.hpp"
#include <boost/filesystem.hpp>

CCPi::ReconstructionImpl::ReconstructionImpl(void)
{
	deviceId = CCPi::dev_Nikon_XTek;
	numberOfProcessors = omp_get_max_threads();
	algorithmId = CCPi::alg_CGLS;
	bHyperThreads = true;
	numberOfIterations = 20;
	resolution = 1;
	regularise = 0.01;
}


CCPi::ReconstructionImpl::~ReconstructionImpl(void)
{
	if(voxels!=NULL)
		delete voxels;
}

void CCPi::ReconstructionImpl::setDeviceId(devices id)
{
	deviceId = id;
}

void CCPi::ReconstructionImpl::setAlgorithmId(algorithms id)
{
	algorithmId = id;
}

void CCPi::ReconstructionImpl::setNumberOfProcessors(int noProcs)
{
	numberOfProcessors = noProcs;
}

void CCPi::ReconstructionImpl::enableHyperThreads()
{
	bHyperThreads = true;
}
void CCPi::ReconstructionImpl::disableHyperThreads()
{
	bHyperThreads = false;
}

void CCPi::ReconstructionImpl::setResolution(int resolution)
{
	this->resolution = resolution;
}
void CCPi::ReconstructionImpl::setNumberOfIterations(int iterations)
{
	numberOfIterations = iterations;
}

void CCPi::ReconstructionImpl::setRegularisation(double value)
{
	regularise = value;
}

void CCPi::ReconstructionImpl::enableBeamHardening()
{
	bBeamHardening = true;
}

void CCPi::ReconstructionImpl::disableBeamHardening()
{
	bBeamHardening = false;
}

void CCPi::ReconstructionImpl::setBeamHardening(bool value)
{
	bBeamHardening = value;
}

void CCPi::ReconstructionImpl::setFilename(std::string name)
{
	filename = name;
}

bool CCPi::ReconstructionImpl::run()
{
	if (!bHyperThreads) {
		int num_processors = numberOfProcessors / 2;
		omp_set_num_threads(num_processors);
	}else{
		omp_set_num_threads(numberOfProcessors);
	}
	bool noerror = true;
	
	CCPi::instrument *device = 0;
	
	switch (deviceId) {
	case CCPi::dev_Nikon_XTek:
		device = new CCPi::Nikon_XTek;
		break;
	default:
		noerror = false;
		break;
	}

	CCPi::reconstruction_alg *recon_algorithm = 0;
	switch (algorithmId) {
	case CCPi::alg_CGLS:
	  recon_algorithm = new CCPi::cgls_3d(numberOfIterations);
	  break;
	case CCPi::alg_SIRT:
	  recon_algorithm = new CCPi::sirt(numberOfIterations);
	  break;
	case CCPi::alg_MLEM:
	  recon_algorithm = new CCPi::mlem(numberOfIterations);
	  break;
	case CCPi::alg_CGLS_Tikhonov:
	  recon_algorithm = new CCPi::cgls_tikhonov(numberOfIterations, regularise);
	  break;
	case CCPi::alg_CGLS_TVreg:
	  recon_algorithm = new CCPi::cgls_tv_reg(numberOfIterations, regularise);
	  break;
	default:
		noerror = false;
	}

	if (noerror) {
		noerror = false;
		bool phantom = false;
		real rotation_centre = -1.0;
		boost::filesystem::path p(filename);
		if (device->setup_experimental_geometry(p.parent_path().generic_string(), p.filename().generic_string(), rotation_centre, resolution, phantom)) {
		  int nx_voxels = 0;
		  int ny_voxels = 0;
		  int maxz_voxels = 0;
		  int nz_voxels = 0;
		  int block_size = 0;
		  int block_step = 0;
		  calculate_block_sizes(nx_voxels, ny_voxels, nz_voxels, maxz_voxels, block_size, block_step, 1,
		    0, resolution, device, recon_algorithm->supports_blocks());
		  int z_data_size = block_size * resolution;
		  int z_data_step = block_step * resolution;
		  device->set_v_block(z_data_size);
		  int block_offset = 0;
		  int z_data_offset = block_offset * resolution;
				if (device->finish_voxel_geometry(voxel_origin, voxel_size, nx_voxels, ny_voxels, nz_voxels)) {
					if (device->read_scans(p.parent_path().generic_string(), 0, z_data_size, true, phantom)) {
						voxels = new voxel_data(boost::extents[nx_voxels][ny_voxels][nz_voxels], boost::c_storage_order());
						init_data(*voxels, nx_voxels, ny_voxels, nz_voxels);
						if (bBeamHardening)
							device->apply_beam_hardening();
						noerror = recon_algorithm->reconstruct(device, *voxels, voxel_origin, voxel_size);
					}
				}
		}
	}
	if (!noerror)
		delete voxels;
	return noerror;
}

bool CCPi::ReconstructionImpl::saveResults(std::string outputFilename, CCPi::output_format outputFormat)
{
	bool clamp_output = true;
	const voxel_data::size_type *s = voxels->shape();
	clamp_min(*voxels, 0.0, s[0], s[1], s[2]);
	CCPi::write_results(outputFilename, *voxels, voxel_origin, voxel_size, 0, (int)s[2], outputFormat, clamp_output);
	return true;
}