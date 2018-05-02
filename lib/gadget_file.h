#ifndef GADGET_FILE_H
#define GADGET_FILE_H 1

#include <string>
#include <cstdio>

struct GadgetFileHeader {
  int      np[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned int np_total[6];
  int      flag_cooling;
  int      num_files;
  double   boxsize;
  double   omega0;
  double   omega_lambda;
  double   hubble_param; 
  int flag_stellarage;
  int flag_metals;
  unsigned int np_total_highword[6];
  int  flag_entropy_instead_u;
  char fill[60];
  //char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];
           /* fills to 256 Bytes */
};

class GadgetFile {
 public:
  GadgetFile(const char filename_[]);
  void open();
  void close();
  void read(const char component,
	    const size_t ibegin, const size_t iend,
	    const size_t stride,
	    float* const buf);
 private:
  std::string filename;
  FILE* fp;
  GadgetFileHeader h;
  static const int size_header= 256;
};

#endif
