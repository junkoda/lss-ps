//
// Catalogue is an array of particles
//

#include <vector>
#include <algorithm> // std::max
#include <chrono>
#include <cstdio>

#include "config.h"
#include "msg.h"
#include "catalogue.h"

using namespace std;

//
//
//
static inline int skip_space(const char line[], int i)
{
  // return the first position after i that is not white space
  while(line[i] != '\0') {
    if(line[i] == ' ' || line[i] == '\t')
      i++;
    else
      break;
  }
  return i;
}

static inline int end_of_number(const char line[], int i)
{
  // return the first position after i that is not number (not space)
  while(line[i] != '\0') {
    if(line[i] == ' ' || line[i] == '\t' || line[i] == '\n')
      break;
    else
      i++;
  }
  return i;
}

static inline void split(char line[], std::vector<double>& v)
{
  // split a space spearated char line into a vector of double
  int i=0;
  int ibegin= 0;

  v.clear();
  
  while(line[i] != '\0') {
    ibegin= skip_space(line, i);
    if(line[i] == '\n' || line[i] == '\0')
      break;
    
    i= end_of_number(line, ibegin+1);
    line[i]= '\0';
    double x= atof(line + ibegin);
    v.push_back(x);
    i++;
  }
}

//
// CatalogueFileAscii
//

CatalogueFileAscii::CatalogueFileAscii(const char filename_[],
				       const std::vector<int> ipos_,
				       const std::vector<int> iweights_,
				       const int inbar_,
				       const double Pest_) :

  fp(0)
{
  filename= string(filename_);
  ipos= ipos_; 
  iweights= iweights_;
  inbar= inbar_;
  Pest= Pest_;

  if(ipos.size() != 3) {
    msg_printf(msg_fatal, "Error: ipos must have 3 numbers\n");
    throw ValError();
  }

  // compute imax
  // imax is the maximum column number among ipos, iweights, and inbar
  imax= inbar;
  
  for(int k=0; k<3; ++k) {
    imax= max(imax, ipos[k]);
    assert(ipos[k] >= 0);
  }
  
  imax= max(imax, inbar);
  for(int i : iweights) {
    imax = max(imax, i);
    assert(i >= 0);
  }
}

void CatalogueFileAscii::open()
{
  fp= fopen(filename.c_str(), "r");
  
  if(fp == 0) {
    msg_printf(msg_fatal,
	       "Error: unable to open catalogue file %s\n", filename.c_str());
    throw IOError();
  }
}

void CatalogueFileAscii::close()
{
  fclose(fp);
  fp= 0;
}
  
void CatalogueFileAscii::read(const size_t nbatch, vector<Particle>& cat)
{
  assert(fp); // File not opened; called open()?
  cat.clear();
  
  const bool use_fkp_weight= Pest > 0.0;
  constexpr size_t nbuf= 1024; // maximum length of the line

  char buf[nbuf];
  vector<double> v;

  Particle p;
  p.nbar= 1.0;

  size_t n= 0;

  #pragma omp critical (__CATALOG_READ_BUF__)
  {
    while(n < nbatch) {
      char* ret= fgets(buf, nbuf - 1, fp);
      if(ret == NULL) break;

      if(buf[0] == '#')
	continue;

      split(buf, v); assert(imax < v.size());

      p.x[0]= v[ipos[0]];
      p.x[1]= v[ipos[1]];
      p.x[2]= v[ipos[2]];

      p.w= 1.0;
      for(const int i : iweights)
	p.w *= v[i];
    
      if(inbar >= 0)
	p.nbar = v[inbar];

      if(use_fkp_weight)
	p.w /=  (1.0 + p.nbar*Pest);

      cat.push_back(p);
      n++;
    }
  }
}

// Read ascii file to Catalogue
// File format:
//   x y z -- white space separated cartisian corrdinate in 1/h Mpc
//   Lines starting with # is a commnet


void catalogue_read_text(Catalogue* const cat, const char filename[])
{
  auto ts = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Reading catalogue %s\n", filename);
    
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    msg_printf(msg_fatal,
	       "Error: unable to open catalogue file %s\n", filename);
    throw IOError();
  }

  
  const size_t nbuf= 1024;
  char buf[nbuf];

  Particle p;
  p.w = 1;
  p.nbar = 0.0;
  
  while(fgets(buf, nbuf - 1, fp)) {
    if(buf[0] == '#')
      continue;

#ifdef DOUBLEPRECISION
    int ret= sscanf(buf, "%le %le %le", p.x, p.x + 1, p.x + 2);
#else
    int ret= sscanf(buf, "%e %e %e", p.x, p.x + 1, p.x + 2);
#endif
    
    if(ret != 3) {
      msg_printf(msg_error, "Error: unable to read 3 numbers, %s\n", buf);
      throw IOError();
    }

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
  msg_printf(msg_verbose, "Time read catalogue txt %le\n",
	     std::chrono::duration<double>(te - ts).count());

}

void catalogue_compute_range(const Catalogue& cat,
			     double x0_out[], double& boxsize_out)
{
  // Input:
  //   catalogue cat
  // Output:
  //   x0: minimum corrdinate
  //   boxsize:

  double left[3], right[3];

  for(int k=0; k<3; ++k)
    left[k]= right[k]= cat.front().x[k];
    
  for(auto&& p : cat) {
    for(int k=0; k<3; ++k) {
      if(p.x[k] < left[k]) left[k]= p.x[k];
      if(p.x[k] > right[k]) right[k]= p.x[k];
    }
  }

  boxsize_out= 0.0;
  for(int k=0; k<3; ++k) {
    x0_out[k]= left[k];
    double diff= right[k] - left[k];
    if(diff > boxsize_out) boxsize_out= diff;
  }
}
