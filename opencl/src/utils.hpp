
#ifndef CCPI_UTILITY_ROUTINES
#define CCPI_UTILITY_ROUTINES

namespace CCPi {

  void split_path_and_name(const std::string fullname, std::string &path,
			   std::string &name);
  void combine_path_and_name(const std::string path, const std::string name,
			     std::string &fullname);
  bool access(const char name[]);

}

#endif // CCPI_UTILITY_ROUTINES
