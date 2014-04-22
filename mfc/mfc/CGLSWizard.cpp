// CGLSWizard.cpp : implementation file
//

#include "stdafx.h"
#include "mfc.h"
#include "omp.h"
#include "CGLSWizard.h"
#include "Progress_page.h"
#include "src/ui_calls.hpp"
#include "src/utils.hpp"

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
		recon_algorithm = new CCPi::cgls_base(niterations);
		break;
	default:
		//std::cerr << "ERROR: Unknown algorithm\n";
		ok = false;
	}
	if (ok) {
		ok = false;
		bool phantom = false;
		//std::string path;
		//std::string name;
		//CCPi::split_path_and_name(filename, path, name);
		if (device->setup_experimental_geometry(get_data_path(), get_data_name(), phantom)) {
			if (device->read_data_size(get_data_path(), phantom)) {
				int nx_voxels = device->get_num_h_pixels() / resolution;
				if (device->get_num_h_pixels() % resolution != 0)
					nx_voxels++;
				int ny_voxels = nx_voxels;
				int nz_voxels = device->get_num_v_pixels() / resolution;
				if (device->get_num_v_pixels() % resolution != 0)
					nz_voxels++;
				int z_data_size = nz_voxels * resolution;
				device->set_v_block(z_data_size);
				if (device->finish_voxel_geometry(voxel_origin, voxel_size, nx_voxels, ny_voxels, nz_voxels)) {
					if (device->read_scans(get_data_path(), 0, z_data_size, true, phantom)) {
						voxels = new voxel_data(boost::extents[nx_voxels][ny_voxels][nz_voxels], boost::fortran_storage_order());
						for (int i = 0; i < nz_voxels; i++) {
							for (int j = 0; j < ny_voxels; j++) {
								for (int k = 0; k < nx_voxels; k++) {
									(*voxels)[k][j][i] = 0.0;
								}
							}
						}
						if (beam_harden)
							device->apply_beam_hardening();
						//if (fast_projection and first)
						//	device->setup_projection_matrix(voxel_origin, voxel_size, nx_voxels, ny_voxels, nz_voxels);
						ok = recon_algorithm->reconstruct(device, *voxels, voxel_origin, voxel_size);
					}
				}
			}
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