#pragma once


// file_browser

class file_browser : public CFileDialog
{
	DECLARE_DYNAMIC(file_browser)

public:
	file_browser(BOOL bOpenFileDialog, // TRUE for FileOpen, FALSE for FileSaveAs
		LPCTSTR lpszDefExt = NULL,
		LPCTSTR lpszFileName = NULL,
		DWORD dwFlags = OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
		LPCTSTR lpszFilter = NULL,
		CWnd* pParentWnd = NULL);
	virtual ~file_browser();

protected:
	DECLARE_MESSAGE_MAP()
};


