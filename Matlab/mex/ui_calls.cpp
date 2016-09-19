
#include <string>
#include "mex.h"

void report_error(const std::string message)
{
  mexPrintf("%s\n", message.c_str());
}
