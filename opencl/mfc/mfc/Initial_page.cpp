// Initial_page.cpp : implementation file
//

#include "stdafx.h"
#include "mfc.h"
#include "Initial_page.h"
#include "afxdialogex.h"
#include "CGLSWizard.h"

// Initial_page dialog

IMPLEMENT_DYNAMIC(Initial_page, CPropertyPage)

Initial_page::Initial_page()
	: CPropertyPage(Initial_page::IDD)
{

}

Initial_page::~Initial_page()
{
}

void Initial_page::DoDataExchange(CDataExchange* pDX)
{
	CPropertyPage::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(Initial_page, CPropertyPage)
	ON_BN_CLICKED(IDC_RADIO_XTEK, &Initial_page::OnBnClickedRadioXtek)
	ON_BN_CLICKED(IDC_CGLS, &Initial_page::OnBnClickedCgls)
	ON_BN_CLICKED(IDC_HT_CHECK, &Initial_page::OnBnClickedHtCheck)
	ON_BN_CLICKED(IDC_MLEM, &Initial_page::OnBnClickedMLEM)
	ON_BN_CLICKED(IDC_SIRT, &Initial_page::OnBnClickedSIRT)
	ON_BN_CLICKED(IDC_CGLS_TIK, &Initial_page::OnBnClickedCGLSTik)
	ON_BN_CLICKED(IDC_RADIO1, &Initial_page::OnBTNClickedCGLS_TV)
END_MESSAGE_MAP()


// Initial_page message handlers


void Initial_page::OnBnClickedRadioXtek()
{
	my_sheet->set_instrument(CCPi::dev_Nikon_XTek);
}

void Initial_page::OnBnClickedCgls()
{
	my_sheet->set_algorithm(CCPi::alg_CGLS);
}

BOOL Initial_page::OnSetActive()
{
	CPropertySheet *parent = (CPropertySheet *)GetParent();
	parent->SetWizardButtons(PSWIZB_NEXT);
	CheckDlgButton(IDC_RADIO_XTEK, 1);
	CheckDlgButton(IDC_CGLS, 1);

	return CPropertyPage::OnSetActive();
}

void Initial_page::OnBnClickedHtCheck()
{
	my_sheet->toggle_hyper_threads();
}

void Initial_page::OnBnClickedMLEM()
{
  my_sheet->set_algorithm(CCPi::alg_MLEM);
}

void Initial_page::OnBnClickedSIRT()
{
  my_sheet->set_algorithm(CCPi::alg_SIRT);
}

void Initial_page::OnBnClickedCGLSTik()
{
  my_sheet->set_algorithm(CCPi::alg_CGLS_Tikhonov);
}


void Initial_page::OnBTNClickedCGLS_TV()
{
  my_sheet->set_algorithm(CCPi::alg_CGLS_TVreg);
}
