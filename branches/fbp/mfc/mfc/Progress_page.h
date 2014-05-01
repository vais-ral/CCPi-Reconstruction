#pragma once
#include "afxcmn.h"

#include <string>
// Progress_page dialog

class Progress_page : public CPropertyPage
{
	DECLARE_DYNAMIC(Progress_page)

public:
	Progress_page();
	virtual ~Progress_page();
	void initialise_progress(const int length, const char label[]);
	void update_progress(const int value);
	void send_output(const std::string str);

// Dialog Data
	enum { IDD = IDD_PROGRESS_PAGE };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	virtual BOOL OnSetActive();
private:
	CProgressCtrl progress_bar;
};
