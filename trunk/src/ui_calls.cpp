
#include <iostream>
#include <cstdio>
#include "base_types.hpp"
#include "ui_calls.hpp"

static std::string buffer;

void initialise_progress(const int length, const char label[])
{
  std::cout << label << '\n';
}

void update_progress(const int value)
{
}

void report_error(const std::string message)
{
  std::cerr << message << '\n';
}

void report_error(const std::string message, const std::string arg)
{
  std::cerr << message << arg << '\n';
}

void report_error(const std::string message, const std::string arg1,
		  const std::string arg2)
{
  std::cerr << message << arg1 << arg2 << '\n';
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
  std::cout << buffer << '\n';
  buffer = "";
}
