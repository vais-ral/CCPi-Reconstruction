// CGLS_params.cpp : implementation file
//

#include "stdafx.h"
#include "mfc.h"
#include "CGLS_params.h"
#include "afxdialogex.h"
#include "CGLSWizard.h"


// CGLS_params dialog

IMPLEMENT_DYNAMIC(CGLS_params, CPropertyPage)

CGLS_params::CGLS_params()
	: CPropertyPage(CGLS_params::IDD)
{

}

CGLS_params::~CGLS_params()
{
}

void CGLS_params::DoDataExchange(CDataExchange* pDX)
{
	CPropertyPage::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_ITER_SLIDER, iter_slider);
	DDX_Control(pDX, IDC_RES_SLIDER, pix_slider);
}


BEGIN_MESSAGE_MAP(CGLS_params, CPropertyPage)
	ON_WM_HSCROLL()
	ON_BN_CLICKED(IDC_BEAM_HARDEN, &CGLS_params::OnBnClickedBeamHarden)
END_MESSAGE_MAP()


// CGLS_params message handlers


void CGLS_params::OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar)
{
	if (pScrollBar == (CScrollBar *)&pix_slider) {
		CGLSWizard *s = (CGLSWizard *)GetParent();
		int resolution = pix_slider.GetPos();
		s->set_resolution(resolution);
		CString text;
		text.Format(_T("%d"), resolution);
		SetDlgItemText(IDC_RES_VALUE, text);
	}
	else if (pScrollBar == (CScrollBar *)&iter_slider) {
		int niterations = iter_slider.GetPos();
		CGLSWizard *s = (CGLSWizard *)GetParent();
		s->set_iterations(niterations);
		CString text;
		text.Format(_T("%d"), niterations);
		SetDlgItemText(IDC_ITER_VALUE, text);
	}

	CPropertyPage::OnHScroll(nSBCode, nPos, pScrollBar);
}


void CGLS_params::OnBnClickedBeamHarden()
{
	CGLSWizard *s = (CGLSWizard *)GetParent();
	s->toggle_beam_harden();
}


BOOL CGLS_params::OnSetActive()
{
	my_sheet->SetWizardButtons(PSWIZB_BACK | PSWIZB_NEXT);
	pix_slider.SetRangeMin(1);
	pix_slider.SetRangeMax(16);
	pix_slider.SetPos(my_sheet->get_resolution());
	CString text;
	text.Format(_T("%d"), my_sheet->get_resolution());
	SetDlgItemText(IDC_RES_VALUE, text);
	iter_slider.SetRangeMin(5);
	iter_slider.SetRangeMax(30);
	iter_slider.SetPos(my_sheet->get_iterations());
	text.Format(_T("%d"), my_sheet->get_iterations());
	SetDlgItemText(IDC_ITER_VALUE, text);
	CheckDlgButton(IDC_BEAM_HARDEN, (int)my_sheet->get_beam_harden());
	return CPropertyPage::OnSetActive();
}
