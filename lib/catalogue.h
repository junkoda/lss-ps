#ifndef _CATALOGUE_H
#define _CATALOGUE_H 1

#include <iostream>
#include <vector>
#include <cassert>
#include <cstdlib>

#include "config.h"
#include "error.h"
#include "msg.h"

//#define DEBUG

struct Particle {
  Float x[3], w, nbar;
};

class Catalogue : public std::vector<Particle> {
 public:
  Float boxsize;
};

void catalogue_read_text(Catalogue* const cat, const char filename[]);

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

static inline void read_line(char line[], std::vector<double>& v)
{
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

class XYZ {
 public:
  XYZ(const std::vector<int>& xyz, const std::vector<int>& weights, const int inbar_) {
    assert(xyz.size() == 3);
    ix= xyz[0] - 1;
    iy= xyz[1] - 1;
    iz= xyz[2] - 1;
    for(auto i : weights)
      weight_cols.push_back(i - 1);
    inbar= inbar_ - 1;

    msg_printf(msg_info, "Columns for xyz: %d %d %d\n",
	       ix + 1, iy + 1, iz + 1);
    msg_printf(msg_info, "Columns for weights:");
    for(auto i : weight_cols)
      msg_printf(msg_info, " %d", i + 1);
    if(weight_cols.empty())
      msg_printf(msg_info, "none");
    msg_printf(msg_info, "\n");

    msg_printf(msg_info, "Column for mean number density: ");
    if(inbar > 0)
      msg_printf(msg_info, "%d\n", inbar + 1);
    else
      msg_printf(msg_info, "none\n");
  }
  
  
  void operator()(const std::vector<double>& v, Particle& p) const {
#ifdef DEBUG
    assert(ix < v.size());
    assert(iy < v.size());
    assert(iz < v.size());
    assert(inbar < (int) v.size());
#endif
    
    p.x[0]= v[ix];
    p.x[1]= v[iy];
    p.x[2]= v[iz];
    p.w= 1.0;
    p.nbar= 0.0;
    for(const int i : weight_cols) {
#ifdef DEBUG
      assert(i < v.size());
#endif
      p.w *= v[i];
    }
    if(inbar >= 0)
      p.nbar = v[inbar];
  }

  void push_back(const int icol_weight);
 private:
  int ix, iy, iz, inbar;
  std::vector<int> weight_cols;	
};

template<typename Converter>
void catalogue_read(Catalogue* const cat, const char filename[],
		    Converter f)
{
  auto ts = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Reading catalogue %s\n", filename);
    
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
  
  while(fgets(buf, nbuf - 1, fp)) {
    if(buf[0] == '#')
      continue;

    read_line(buf, v);
    f(v, p);
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


#endif
