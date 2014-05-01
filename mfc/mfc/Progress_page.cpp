// Progress_page.cpp : implementation file
//

#include "stdafx.h"
#include "mfc.h"
#include "Progress_page.h"
#include "afxdialogex.h"
#include "CGLSWizard.h"

//CWinThread *compute_thread = 0;

// Progress_page dialog

IMPLEMENT_DYNAMIC(Progress_page, CPropertyPage)

Progress_page::Progress_page()
	: CPropertyPage(Progress_page::IDD)
{

}

Progress_page::~Progress_page()
{
}

void Progress_page::DoDataExchange(CDataExchange* pDX)
{
	CPropertyPage::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_PROGRESS_BAR, progress_bar);
}


BEGIN_MESSAGE_MAP(Progress_page, CPropertyPage)
END_MESSAGE_MAP()


// Progress_page message handlers

extern UINT __cdecl main_loop(LPVOID param);

BOOL Progress_page::OnSetActive()
{
	CPropertySheet *s = (CPropertySheet *)GetParent();
	s->SetWizardButtons(0);
	CWinThread *compute_thread = AfxBeginThread(&main_loop, my_sheet->GetSafeHwnd(), THREAD_PRIORITY_NORMAL, 0, 0, 0);
	//compute_thread->m_bAutoDelete = FALSE;
	return CPropertyPage::OnSetActive();
}

void Progress_page::initialise_progress(const int length, const char label[])
{
	SetDlgItemText(IDC_PROGRESS_LABEL, CA2T(label));
	progress_bar.SetRange(0, length);
	progress_bar.SetPos(0);
}

void Progress_page::update_progress(const int value)
{
	progress_bar.SetPos(value);
}

void Progress_page::send_output(const std::string str)
{
	char buff[1024];
	strncpy(buff, str.c_str(), 1024);
	buff[1023] = '\0';
	SetDlgItemText(IDC_MESSAGES, CA2T(buff));
}
