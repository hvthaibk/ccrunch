#include <boost/python.hpp>

char const* hello()
{
   return "Hello, world! This is the Python wrapper for the GEM project.";
}

BOOST_PYTHON_MODULE(libGemPython)
{
    using namespace boost::python;

    def("hello", hello);
}
