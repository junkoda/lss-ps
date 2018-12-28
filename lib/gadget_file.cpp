#include <iostream> // debug
#include <string>
#include <cstdio>
#include <cassert>
#include "msg.h"
#include "error.h"
#include "gadget_file.h"

using namespace std;

GadgetFile::GadgetFile(const char filename_[])
{
  filename = string(filename_);
}

void GadgetFile::open()
{
  fp= fopen(filename.c_str(), "r");
  if(fp == 0) {
    msg_printf(msg_error, "Error: Gadget file not found: %s",
	       filename.c_str());
    throw FileNotFoundError();
  }

  int sep_begin, sep_end;
  int ret= fread(&sep_begin, sizeof(int), 1, fp); assert(ret == 1);
  if(sep_begin != size_header) {
    msg_printf(msg_error, "Error: Unexpeced format for Gadget binary file: %s",
	       filename.c_str());
    throw IOError();
  }

  assert(sizeof(GadgetFileHeader) == size_header);
  ret= fread(&h, sizeof(GadgetFileHeader), 1, fp); assert(ret == 1);
  
  ret= fread(&sep_end, sizeof(int), 1, fp); assert(ret == 1);
  assert(sep_begin == sep_end);
}

void GadgetFile::read(const char component,
		      const size_t ibegin, const size_t iend,
		      const size_t stride,
		      float* buf) {
  assert(ibegin <= iend);
  assert(ibegin <= static_cast<size_t>(h.np[1]));
  assert(iend <= static_cast<size_t>(h.np[1]));

  const size_t np_all_type =
    (size_t) h.np[0] + h.np[1] + h.np[2] + h.np[3] + h.np[4] + h.np[5];

  if(component == 'x') {
    fseek(fp, size_header + 8, SEEK_SET);
  }
  else if(component == 'v') {
    fseek(fp, size_header + 8 + 3*sizeof(float)*np_all_type + 8, SEEK_SET);
  }
  else {
    msg_printf(msg_error, "Error: unknown component to read: %c", component);
    throw ValueError();
  }

  int sep_begin, sep_end;
  int ret= fread(&sep_begin, sizeof(int), 1, fp); assert(ret == 1);

  assert(static_cast<size_t>(sep_begin) == 3*h.np[1]*sizeof(float));

  fseek(fp, 3*sizeof(float)*ibegin, SEEK_CUR);

  const size_t np= iend - ibegin;
  for(size_t i=0; i<np; ++i) {
    ret= fread(buf, sizeof(float), 3, fp); assert(ret == 3);
    buf = (float*)(reinterpret_cast<char*>(buf) + stride);
  }

  fseek(fp, 3*sizeof(float)*(h.np[1] - iend), SEEK_CUR);

  ret= fread(&sep_end, sizeof(int), 1, fp); assert(ret == 1);

  assert(sep_end == sep_begin);
}

void GadgetFile::close()
{
  fclose(fp);
}
