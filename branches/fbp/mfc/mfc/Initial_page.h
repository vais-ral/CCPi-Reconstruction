#pragma once


// Initial_page dialog

class Initial_page : public CPropertyPage
{
	DECLARE_DYNAMIC(Initial_page)

public:
	Initial_page();
	virtual ~Initial_page();

// Dialog Data
	enum { IDD = IDD_DEV_PAGE };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedRadioXtek();
	afx_msg void OnBnClickedCgls();
	virtual BOOL OnSetActive();
	afx_msg void OnBnClickedHtCheck();
	afx_msg void OnBnClickedMLEM();
	afx_msg void OnBnClickedSIRT();
};
