// Load_page.cpp : implementation file
//

#include "stdafx.h"
#include "mfc.h"
#include "Load_page.h"
#include "afxdialogex.h"
#include "CGLSWizard.h"
#include "file_browser.h"


// Load_page dialog

IMPLEMENT_DYNAMIC(Load_page, CPropertyPage)

Load_page::Load_page()
: CPropertyPage(Load_page::IDD), browser(0)
{

}

Load_page::~Load_page()
{
	if (browser != 0)
		delete browser;
}

void Load_page::DoDataExchange(CDataExchange* pDX)
{
	CPropertyPage::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(Load_page, CPropertyPage)
	ON_BN_CLICKED(IDC_LOAD_BROWSE, &Load_page::OnBnClickedLoadBrowse)
END_MESSAGE_MAP()


// Load_page message handlers


void Load_page::OnBnClickedLoadBrowse()
{
	CPropertySheet *s = (CPropertySheet *)GetParent();
	//LPCTSTR filter = _T("XTek file (*.xtekct)|*.xtekct|");
	//file_browser fb(TRUE, NULL, NULL, OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, filter, s);
	if (browser->DoModal() == IDOK) {
		s->SetWizardButtons(PSWIZB_BACK | PSWIZB_NEXT);
		set_filename();
	}
}


BOOL Load_page::OnSetActive()
{
	CPropertySheet *s = (CPropertySheet *)GetParent();
	s->SetWizardButtons(PSWIZB_BACK);
	BOOL result = CPropertyPage::OnSetActive();
	LPCTSTR filter = _T("XTek file (*.xtekct)|*.xtekct|");
	browser = new file_browser(TRUE, NULL, NULL, OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, filter, s);
	if (browser->DoModal() == IDOK) {
		s->SetWizardButtons(PSWIZB_BACK | PSWIZB_NEXT);
		set_filename();
	}
	return result;
}

void Load_page::set_filename()
{
	CString name = browser->GetFileName();
	if (name.GetLength() < 1024) {
		char data_file[1024];
		// its w_char so might be a problem here
		for (int i = 0; i < name.GetLength(); i++)
			data_file[i] = (char)name[i];
		data_file[name.GetLength()] = '\0';
		my_sheet->set_data_name(data_file);
	}
	name = browser->GetFolderPath();
	if (name.GetLength() < 1024) {
		char data_file[1024];
		// its w_char so might be a problem here
		for (int i = 0; i < name.GetLength(); i++)
			data_file[i] = (char)name[i];
		data_file[name.GetLength()] = '\0';
		my_sheet->set_data_path(data_file);
	}
}
