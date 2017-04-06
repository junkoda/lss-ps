#ifndef _CATALOGUE_H
#define _CATALOGUE_H 1

#include <iostream>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "config.h"
#include "error.h"
#include "msg.h"

struct Particle {
  Float x[3], w, nbar;
};

class Catalogue : public std::vector<Particle> {
 public:
  Float boxsize;
};

class CatalogueFile {
 public:
  CatalogueFile() = default;
  virtual ~CatalogueFile() = default;
  virtual void open() = 0;
  virtual void read(const size_t nbuf, std::vector<Particle>& v) = 0;
  virtual void close() = 0;
};

class CatalogueFileAscii : public CatalogueFile {
 public:
  CatalogueFileAscii(const char filename_[],
		     const std::vector<int> ipos_,
		     const std::vector<int> iweights_,
		     const int inbar_,
		     const double Pest_);
  virtual ~CatalogueFileAscii() = default;

  virtual void open();
  virtual void read(const size_t nbuf, std::vector<Particle>& v);
  virtual void close();

 private:
  std::string filename;
  std::vector<int> ipos, iweights;
  int inbar;
  int imax;
  double Pest;
  FILE* fp;
};



void catalogue_read_text(Catalogue* const cat, const char filename[]);


//
// Convertor
//
enum class AngleUnit {rad, deg};

struct XYZ {
  // When the input file is xyz, nothing is needed to do
  void operator()(Particle& p) const {}
};

class Spherical {
 public:
  Spherical(const AngleUnit angle_unit) {
    if(angle_unit == AngleUnit::rad)
      unit_angle= 1.0;
    else if(angle_unit == AngleUnit::deg)
      unit_angle= M_PI/180.0;
    else
      assert(false);

    msg_printf(msg_info, "Using sperical coordinate position = (RA, Dec, r)\n");
  }
  
  void operator()(Particle& p) const {
    // (RA, Dec, r) -> (x, y, z)
    double r= p.x[2];
    double theta= unit_angle*p.x[1]; // measured from xy plane
    double phi= -unit_angle*p.x[0];  // minus for definition of RA
    double rcosO= r*cos(theta);
    p.x[0]= rcosO*cos(phi);
    p.x[1]= rcosO*sin(phi);
    p.x[2]= r*sin(theta);
  }
 private:
  double unit_angle;
};

/*
template<typename Converter>
void catalogue_read(Catalogue* const cat, const char filename[],
		    const std::vector<int>& ipos,
		    const std::vector<int>& iweights, const int inbar1,
		    const bool Pest,
		    Converter convert)
{
  // Read catalogue from an ascii file
  //
  // ipos:     column numbers for positions; the column index starts from 1
  // iweights: column numbers for weights
  // inbar1:   column number for nbar
  //
  // Pest:  1/(1 + nbar*Pest) in FKP weight.
  //        set Pest to 0 when FKP weight is not used
  //
  // The position is converted to xyz by the Convertor
  //
  
  auto ts = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Reading catalogue %s\n", filename);

  const int ipos1= ipos[0] - 1; assert(ipos1 >= 0);
  const int ipos2= ipos[1] - 1; assert(ipos2 >= 0);
  const int ipos3= ipos[2] - 1; assert(ipos3 >= 0);
  std::vector<int> weight_cols;
  
  for(auto i : iweights)
    weight_cols.push_back(i - 1);
  const int inbar= inbar1 - 1;

  const bool use_fkp_weight= Pest > 0.0;

  msg_printf(msg_info, "Columns for positions: %d %d %d\n",
	     ipos1 + 1, ipos2 + 1, ipos3 + 1);
  msg_printf(msg_info, "Columns for weights:");
  for(auto i : weight_cols)
    msg_printf(msg_info, " %d", i + 1);
  if(weight_cols.empty())
    msg_printf(msg_info, "none");
  msg_printf(msg_info, "\n");


  if(inbar > 0)
    msg_printf(msg_info, "Column for mean number density: %d", inbar + 1);
  else
    msg_printf(msg_info, "Column for mean number density: none");

  
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    msg_printf(msg_fatal,
	       "Error: unable to open catalogue file %s\n", filename);
    throw IOError();
  }
  
  constexpr size_t nbuf= 1024;
  char buf[nbuf];

  Particle p;
  std::vector<double> v;
  p.w = 1;
  p.nbar= 1.0;
  
  while(fgets(buf, nbuf - 1, fp)) {
    if(buf[0] == '#')
      continue;

    split(buf, v);

#ifdef DEBUG
    assert(ix < v.size());
    assert(iy < v.size());
    assert(iz < v.size());
    assert(inbar < (int) v.size());
#endif

    p.x[0]= v[ipos1];
    p.x[1]= v[ipos2];
    p.x[2]= v[ipos3];
    p.w= 1.0;

    for(const int i : weight_cols) {
      p.w *= v[i];
    }
    if(inbar >= 0) {
      p.nbar = v[inbar];
      //std::cerr << p.nbar << std::endl;
    }

    if(use_fkp_weight)
      p.w /=  (1.0 + p.nbar*Pest);


    convert(p);
    cat->push_back(p);
  }

  int ret = fclose(fp);
  if(ret != 0) {
    msg_printf(msg_error,
	       "Error: error upon closing the catalogue file %s\n", filename);
    throw IOError();
  }

  msg_printf(msg_verbose, "Read %lu particles\n", cat->size());
  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Time read catalogue template %le\n",
	     std::chrono::duration<double>(te - ts).count());
}
*/

void catalogue_compute_range(const Catalogue& cat,
			     double x0_out[], double& boxsize_out);

#endif
