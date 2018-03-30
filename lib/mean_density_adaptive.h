#ifndef MEAN_DENSITY_ADAPTIVE_H
#define MEAN_DENSITY_ADAPTIVE_H 1

#include <vector>
#include "config.h"    // -> Float
#include "catalogue.h" // -> Paritcle

struct KDPoint {
  Float x[3]; // cartisian corrdinate of this point
  Float w;    // weight of this point
  Float n_local, n_average;
};

struct Node {
  int k;                   // axis 0,1,2 which this node is dived into two    
  //Float left[3], right[3]; // min & max coordinate of the particles
  Float left, right;       // min & max coordinate of the particles
  index_t ibegin, iend;    // index range of particles in this node
};

class KDTree {
 public:
  //KDTree();
  KDTree(std::vector<KDPoint>& v_, const int quota_);
  ~KDTree();
  //~KDtree();

  void compute_bounding_box(const size_t ibegin, const size_t iend,
			    Float left[], Float right[]);
 private:
  std::vector<KDPoint>& v;
  index_t quota;
  //int nlevel;

  Node* nodes;
  size_t n_nodes;

  void construct_balanced();
  void construct_balanced_recursive(const size_t inode, 
			    const index_t ibegin, const index_t iend,
			    Float left[], Float right[], Float boxsize3[]);
};


//
// Templates
//

#endif
