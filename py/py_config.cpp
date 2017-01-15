#include "config.h"
#include "msg.h"

#include "config.h"
#include "py_config.h"

PyObject* py_config_sizeof_float(PyObject* self, PyObject* args)
{
  // Return sizeof(Float) in the library and the Python module
  return Py_BuildValue("nn", config_sizeof_float(), sizeof(Float));
}
