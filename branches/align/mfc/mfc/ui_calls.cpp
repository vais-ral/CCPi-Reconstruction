
#include "stdafx.h"
#include "mfc.h"
#include "CGLSWizard.h"
#include "src/base_types.hpp"
#include "src/ui_calls.hpp"

static std::string buffer;

void initialise_progress(const int length, const char label[])
{
	my_sheet->initialise_progress(length, label);
}

void update_progress(const int value)
{
	my_sheet->update_progress(value);
}

void report_error(const std::string message)
{
	char buff[1024];
	strncpy(buff, message.c_str(), 1024);
	buff[1023] = '\0';
	my_sheet->MessageBox(CA2T(buff), _T("CGLS Error"), MB_ICONERROR | MB_OK);
}

void report_error(const std::string message, const std::string arg)
{
	char buff[1024];
	strncpy(buff, message.c_str(), 1024);
	buff[1023] = '\0';
	int l = strlen(buff);
	strncat(buff, arg.c_str(), 1024 - l);
	buff[1023] = '\0';
	my_sheet->MessageBox(CA2T(buff), _T("CGLS Error"), MB_ICONERROR | MB_OK);
}

void report_error(const std::string message, const std::string arg1, const std::string arg2)
{
	char buff[1024];
	strncpy(buff, message.c_str(), 1024);
	buff[1023] = '\0';
	int l = strlen(buff);
	strncat(buff, arg1.c_str(), 1024 - l);
	buff[1023] = '\0';
	l = strlen(buff);
	strncat(buff, arg2.c_str(), 1024 - l);
	buff[1023] = '\0';
	my_sheet->MessageBox(CA2T(buff), _T("CGLS Error"), MB_ICONERROR | MB_OK);
}

void add_output(const std::string str)
{
	buffer += str;
}

void add_output(const char c)
{
	buffer += c;
}

void add_output(const int i)
{
	char buff[32];
	snprintf(buff, 32, "%1d", i);
	buffer += buff;
}

void add_output(const sl_int i)
{
	char buff[32];
	snprintf(buff, 32, "%1ld", i);
	buffer += buff;
}

void add_output(const int i, const int w, const bool fill)
{
	char buff[32];
	if (fill)
		snprintf(buff, 32, "%0*d", w, i);
	else
		snprintf(buff, 32, "%*d", w, i);
	buffer += buff;
}

void add_output(const real r)
{
	char buff[32];
	snprintf(buff, 32, "%f", r);
	buffer += buff;
}

void send_output()
{
	my_sheet->send_output(buffer);
	buffer = "";
}