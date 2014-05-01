// file_browser.cpp : implementation file
//

#include "stdafx.h"
#include "mfc.h"
#include "file_browser.h"


// file_browser

IMPLEMENT_DYNAMIC(file_browser, CFileDialog)

file_browser::file_browser(BOOL bOpenFileDialog, LPCTSTR lpszDefExt, LPCTSTR lpszFileName,
		DWORD dwFlags, LPCTSTR lpszFilter, CWnd* pParentWnd) :
		CFileDialog(bOpenFileDialog, lpszDefExt, lpszFileName, dwFlags, lpszFilter, pParentWnd)
{

}

file_browser::~file_browser()
{
}


BEGIN_MESSAGE_MAP(file_browser, CFileDialog)
END_MESSAGE_MAP()



// file_browser message handlers


