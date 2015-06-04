// Results_page.cpp : implementation file
//

#include "stdafx.h"
#include "mfc.h"
#include "Results_page.h"
#include "afxdialogex.h"
#include "file_browser.h"
#include "CGLSWizard.h"


// Results_page dialog

IMPLEMENT_DYNAMIC(Results_page, CPropertyPage)

Results_page::Results_page()
: CPropertyPage(Results_page::IDD), browser(0)
{

}

Results_page::~Results_page()
{
	if (browser != 0)
		delete browser;
}

void Results_page::DoDataExchange(CDataExchange* pDX)
{
	CPropertyPage::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(Results_page, CPropertyPage)
	ON_BN_CLICKED(IDC_TIFF16_OUT, &Results_page::OnBnClickedTiff16Out)
	ON_BN_CLICKED(IDC_FLOAT_OUT, &Results_page::OnBnClickedFloatOut)
	ON_BN_CLICKED(IDC_SAVE_BUTTON, &Results_page::OnBnClickedSaveButton)
END_MESSAGE_MAP()


// Results_page message handlers


void Results_page::OnBnClickedTiff16Out()
{
	my_sheet->set_output(CCPi::unsigned_short_tiff);
}


void Results_page::OnBnClickedFloatOut()
{
	my_sheet->set_output(CCPi::bgs_float_dump);
	//my_sheet->set_output(CCPi::native_dump);
}


void Results_page::OnBnClickedSaveButton()
{
	browser = new file_browser(FALSE, NULL, NULL, OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, NULL, my_sheet);
	if (browser->DoModal() == IDOK) {
		my_sheet->SetWizardButtons(PSWIZB_BACK | PSWIZB_FINISH);
		CString name = browser->GetFileName();
		if (name.GetLength() < 1024) {
			char data_file[1024];
			// its w_char so might be a problem here
			for (int i = 0; i < name.GetLength(); i++)
				data_file[i] = (char)name[i];
			data_file[name.GetLength()] = '\0';
			my_sheet->set_output_name(data_file);
		}
		name = browser->GetFolderPath();
		if (name.GetLength() < 1024) {
			char data_file[1024];
			// its w_char so might be a problem here
			for (int i = 0; i < name.GetLength(); i++)
				data_file[i] = (char)name[i];
			data_file[name.GetLength()] = '\0';
			my_sheet->set_output_path(data_file);
		}
	}
}


BOOL Results_page::OnSetActive()
{
	CPropertySheet *s = (CPropertySheet *)GetParent();
	CheckDlgButton(IDC_TIFF16_OUT, 1);
	s->SetWizardButtons(PSWIZB_BACK | PSWIZB_DISABLEDFINISH);

	return CPropertyPage::OnSetActive();
}
