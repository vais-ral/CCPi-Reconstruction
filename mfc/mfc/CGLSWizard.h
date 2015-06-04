#pragma once

#include <string>
#include "src/base_types.hpp"
#include "src/algorithms.hpp"
#include "src/instruments.hpp"
#include "src/results.hpp"

// CGLSWizard

class CGLSWizard : public CPropertySheet
{
	DECLARE_DYNAMIC(CGLSWizard)

public:
	CGLSWizard(UINT nIDCaption, CWnd* pParentWnd = NULL, UINT iSelectPage = 0);
	CGLSWizard(LPCTSTR pszCaption, CWnd* pParentWnd = NULL, UINT iSelectPage = 0);
	virtual ~CGLSWizard();
	bool main_loop();
	void save_results();
	void initialise_progress(const int length, const char label[]);
	void update_progress(const int value);
	void add_progress();
	void send_output(const std::string str);

protected:
	DECLARE_MESSAGE_MAP()

	friend class Initial_page;
	friend class Load_page;
	friend class CGLS_params;
	friend class Results_page;

	void set_resolution(const int r);
	void set_iterations(const int n);
	void set_instrument(const CCPi::devices d);
	void set_algorithm(const CCPi::algorithms a);
	void toggle_beam_harden();
	void set_data_name(const std::string f);
	void set_data_path(const std::string f);
	void set_output(const CCPi::output_format f);
	void set_output_name(const std::string name);
	void set_output_path(const std::string name);
	void toggle_hyper_threads();

	int get_resolution() const;
	int get_iterations() const;
	CCPi::devices get_instrument() const;
	CCPi::algorithms get_algorithm() const;
	bool get_beam_harden() const;
	std::string get_data_name() const;
	std::string get_data_path() const;
	std::string get_output_name() const;
	std::string get_output_path() const;
	bool use_hyper_threads() const;
	double get_regularisation() const;

private:
	int resolution;
	int niterations;
	CCPi::devices instrument;
	CCPi::algorithms algorithm;
	CCPi::output_format out_format;
	std::string data_file;
	std::string data_path;
	std::string out_file;
	std::string out_path;
	bool beam_harden;
	bool hyper_threads;
	class Progress_page *progress;
public:
	double regularise;
};

extern CGLSWizard *my_sheet;

inline void CGLSWizard::set_resolution(const int r)
{
	resolution = r;
}

inline void CGLSWizard::set_iterations(const int n)
{
	niterations = n;
}

inline void CGLSWizard::set_instrument(const CCPi::devices d)
{
	instrument = d;
}

inline void CGLSWizard::set_algorithm(const CCPi::algorithms a)
{
	algorithm = a;
}

inline void CGLSWizard::toggle_beam_harden()
{
	beam_harden = !beam_harden;
}

inline void CGLSWizard::set_data_name(const std::string f)
{
	data_file = f;
}

inline void CGLSWizard::set_data_path(const std::string f)
{
	data_path = f;
}

inline void CGLSWizard::set_output(const CCPi::output_format f)
{
	out_format = f;
}

inline void CGLSWizard::set_output_name(const std::string name)
{
	out_file = name;
}

inline void CGLSWizard::set_output_path(const std::string name)
{
	out_path = name;
}

inline void CGLSWizard::toggle_hyper_threads()
{
	hyper_threads = !hyper_threads;
}

inline int CGLSWizard::get_resolution() const
{
	return resolution;
}

inline int CGLSWizard::get_iterations() const
{
	return niterations;
}

inline CCPi::devices CGLSWizard::get_instrument() const
{
	return instrument;
}

inline CCPi::algorithms CGLSWizard::get_algorithm() const
{
	return algorithm;
}

inline bool CGLSWizard::get_beam_harden() const
{
	return beam_harden;
}

inline double CGLSWizard::get_regularisation() const
{
  return regularise;
}

inline std::string CGLSWizard::get_data_name() const
{
	return data_file;
}

inline std::string CGLSWizard::get_data_path() const
{
	return data_path;
}

inline std::string CGLSWizard::get_output_name() const
{
	return out_file;
}

inline std::string CGLSWizard::get_output_path() const
{
	return out_path;
}

inline bool CGLSWizard::use_hyper_threads() const
{
	return hyper_threads;
}