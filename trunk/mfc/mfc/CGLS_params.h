#pragma once
#include "afxcmn.h"


// CGLS_params dialog

class CGLS_params : public CPropertyPage
{
	DECLARE_DYNAMIC(CGLS_params)

public:
	CGLS_params();
	virtual ~CGLS_params();

// Dialog Data
	enum { IDD = IDD_CGLS_PAGE };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	afx_msg void OnBnClickedBeamHarden();
	virtual BOOL OnSetActive();
private:
	CSliderCtrl iter_slider;
	CSliderCtrl pix_slider;
};
