
#include <iostream>
#include <cstdio>
#ifdef MATLAB_MEX_FILE
#  ifdef WIN32
#    define snprintf _snprintf
#  endif // WIN32
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // MEX_FILE
#include "ui_calls.hpp"

static std::string buffer;
#include "UserInterfaceQT.h"
UserInterfaceQT ui_calls;

void initialise_progress(const int length, const char label[])
{
	ui_calls.set_initialise_progress(length, label);
//  std::cout << label << '\n';
}

void update_progress(const int value)
{
	ui_calls.set_update_progress(value);
}

void report_error(const std::string message)
{
	ui_calls.set_report_error(message);
//  std::cerr << message << '\n';
}

void report_error(const std::string message, const std::string arg)
{
	ui_calls.set_report_error(message,arg);
//  std::cerr << message << arg << '\n';
}

void report_error(const std::string message, const std::string arg1,
		  const std::string arg2)
{
	ui_calls.set_report_error(message,arg1,arg2);
//  std::cerr << message << arg1 << arg2 << '\n';
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
  snprintf(buff, 32, "%g", r);
  buffer += buff;
}

void send_output()
{
  ui_calls.set_send_output(QString(buffer.c_str()));
  buffer = "";
}
