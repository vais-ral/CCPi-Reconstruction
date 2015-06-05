
#include "base_types.hpp"
#include "utils.hpp"

#ifdef WIN32
#  include <io.h>
#  define DIR_SEPARATOR '\\'
#else
#  define DIR_SEPARATOR '/'
#endif // WIN32

void CCPi::split_path_and_name(const std::string fullname, std::string &path,
			       std::string &name)
{
  int i;
  for (i = fullname.length() - 1; i >= 0; i--) {
    if (fullname[i] == DIR_SEPARATOR)
      break;
  }
  if (i < 0) {
    path = "";
    name = fullname;
  } else {
    path = fullname.substr(0, i + 1);
    name = fullname.substr(i + 1);
  }
}

void CCPi::combine_path_and_name(const std::string path, const std::string name,
				 std::string &fullname)
{
  fullname = path;
  if (path.length() > 0)
	if (path[path.length() - 1] != DIR_SEPARATOR)
		fullname += DIR_SEPARATOR;
  fullname += name;
}

#ifdef WIN32

bool CCPi::access(const char name[])
{
  return _access(name, 0) == 0;
}

#else

bool CCPi::access(const char name[])
{
  return ::access(name, R_OK) == 0;
}

#endif // WIN32
