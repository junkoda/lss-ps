#include "msg.h"
#include "py_msg.h"

PyObject* py_msg_get_loglevel(PyObject* self, PyObject* args)
{
  return Py_BuildValue("i", (int) msg_get_loglevel());
}

PyObject* py_msg_set_loglevel(PyObject* self, PyObject* args)
{
  int lv;
  if(!PyArg_ParseTuple(args, "i", &lv)) {
    return NULL;
  }    
  
  msg_set_loglevel((LogLevel) lv);
  
  Py_RETURN_NONE;
}


