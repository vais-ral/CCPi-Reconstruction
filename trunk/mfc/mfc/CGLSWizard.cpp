// CGLSWizard.cpp : implementation file
//

#include "stdafx.h"
#include "mfc.h"
#include "omp.h"
#include "CGLSWizard.h"
#include "Progress_page.h"
#include "src/ui_calls.hpp"
#include "src/utils.hpp"
#include "src/voxels.hpp"
#include "src/blas.hpp"

CGLSWizard *my_sheet = 0;
static voxel_data *voxels = 0;
static real voxel_origin[3];
static real voxel_size[3];
//extern CWinThread *compute_thread;

// CGLSWizard

IMPLEMENT_DYNAMIC(CGLSWizard, CPropertySheet)

CGLSWizard::CGLSWizard(UINT nIDCaption, CWnd* pParentWnd, UINT iSelectPage)
	:CPropertySheet(nIDCaption, pParentWnd, iSelectPage)
	, resolution(1)
	, niterations(20)
	, instrument(CCPi::dev_Nikon_XTek)
	, algorithm(CCPi::alg_CGLS)
	, out_format(CCPi::unsigned_short_tiff)
	, beam_harden(true), hyper_threads(true), progress(0)
{

}

CGLSWizard::CGLSWizard(LPCTSTR pszCaption, CWnd* pParentWnd, UINT iSelectPage)
	:CPropertySheet(pszCaption, pParentWnd, iSelectPage)
	, resolution(1)
	, niterations(20)
	, instrument(CCPi::dev_Nikon_XTek)
	, algorithm(CCPi::alg_CGLS)
	, out_format(CCPi::unsigned_short_tiff)
	, beam_harden(true), hyper_threads(true), progress(0)
{

}

CGLSWizard::~CGLSWizard()
{
	if (progress != 0)
		delete progress;
}


BEGIN_MESSAGE_MAP(CGLSWizard, CPropertySheet)
END_MESSAGE_MAP()


// CGLSWizard message handlers

UINT __cdecl main_loop(LPVOID param)
{
	bool ok = my_sheet->main_loop();
	//PostMessage((HWND)param, WM_THREADFINISHED, 0, 0);
	if (ok) {
		// Todo - is this safe?
		my_sheet->SetWizardButtons(PSWIZB_NEXT);
		return 0;
	} else {
		return 1;
	}
}

bool CGLSWizard::main_loop()
{
	if (!use_hyper_threads()) {
		int num_processors = omp_get_max_threads() / 2;
		omp_set_num_threads(num_processors);
	}
	bool ok = true;
	CCPi::instrument *device = 0;
	switch (instrument) {
	case CCPi::dev_Nikon_XTek:
		device = new CCPi::Nikon_XTek;
		break;
	default:
		//std::cerr << "ERROR: Unknown device type\n";
		ok = false;
		break;
	}
	CCPi::reconstruction_alg *recon_algorithm = 0;
	switch (algorithm) {
	case CCPi::alg_CGLS:
		recon_algorithm = new CCPi::cgls_3d(niterations);
		break;
	default:
		//std::cerr << "ERROR: Unknown algorithm\n";
		ok = false;
	}
	if (ok) {
		ok = false;
		bool phantom = false;
		real rotation_centre = -1.0;
		if (device->setup_experimental_geometry(get_data_path(), get_data_name(), rotation_centre, phantom)) {
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
					if (device->read_scans(get_data_path(), 0, z_data_size, true, phantom)) {
						voxels = new voxel_data(boost::extents[nx_voxels][ny_voxels][nz_voxels], boost::c_storage_order());
						init_data(*voxels, nx_voxels, ny_voxels, nz_voxels);
						if (beam_harden)
							device->apply_beam_hardening();
						ok = recon_algorithm->reconstruct(device, *voxels, voxel_origin, voxel_size);
					}
				}
			//}
		}
	}
	if (!ok)
		delete voxels;
	return ok;
}

void CGLSWizard::save_results()
{
	bool clamp_output = true;
	std::string output_name;
	CCPi::combine_path_and_name(get_output_path(), get_output_name(), output_name);
	const voxel_data::size_type *s = voxels->shape();
	CCPi::write_results(output_name, *voxels, voxel_origin, voxel_size, 0, (int)s[2], out_format, clamp_output);
	delete voxels;
}

void CGLSWizard::add_progress()
{
	progress = new Progress_page;
	AddPage(progress);
}

void CGLSWizard::initialise_progress(const int length, const char label[])
{
	progress->initialise_progress(length, label);
}

void CGLSWizard::update_progress(const int value)
{
	progress->update_progress(value);
}

void CGLSWizard::send_output(const std::string str)
{
	progress->send_output(str);
}