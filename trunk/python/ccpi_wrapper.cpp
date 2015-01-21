
#include <boost/python.hpp>
#include "reconstruct.hpp"

BOOST_PYTHON_MODULE(ccpi)
{
  using namespace boost::python;
  def ("cgls", reconstruct_cgls);
  def ("tvreg", reconstruct_tvreg);
}
