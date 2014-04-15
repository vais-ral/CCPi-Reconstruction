#pragma once
#include "file_browser.h"


// Load_page dialog

class Load_page : public CPropertyPage
{
	DECLARE_DYNAMIC(Load_page)

public:
	Load_page();
	virtual ~Load_page();

// Dialog Data
	enum { IDD = IDD_LOAD_PAGE };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedLoadBrowse();
	virtual BOOL OnSetActive();
private:
	class file_browser *browser;

	void set_filename();
};
