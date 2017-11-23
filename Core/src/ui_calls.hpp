
#ifndef CCPI_UI_CALLS
#define CCPI_UI_CALLS

#include <string>

void initialise_progress(const int length, const char label[]);
void update_progress(const int value);
void report_error(const std::string message);
void report_error(const std::string message, const std::string arg);
void report_error(const std::string message, const std::string arg1,
		  const std::string arg2);
void add_output(const std::string str);
void add_output(const char c);
void add_output(const int i);
void add_output(const int i, const int w, const bool fill);
void add_output(const sl_int i);
void add_output(const real r);
void send_output();

#endif // CCPI_UI_CALLS
