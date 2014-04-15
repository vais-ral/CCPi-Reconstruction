#pragma once


// Results_page dialog

class Results_page : public CPropertyPage
{
	DECLARE_DYNAMIC(Results_page)

public:
	Results_page();
	virtual ~Results_page();

// Dialog Data
	enum { IDD = IDD_RESULTS_PAGE };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedTiff16Out();
	afx_msg void OnBnClickedFloatOut();
	afx_msg void OnBnClickedSaveButton();
	virtual BOOL OnSetActive();

private:
	class file_browser *browser;
};
