#ifndef MEAN_DENSITY_ADAPTIVE_H
#define MEAN_DENSITY_ADAPTIVE_H 1

#include <vector>
#include "config.h"    // -> Float
#include "catalogue.h" // -> Paritcle

struct KDPoint {
  index_t idx; // original index
  Float x[3];  // cartisian corrdinate of this point
  Float w;     // weight of this point
  Float rk;    // distance to the kth nearest neighbor
  Float n_local;
};

struct Node {
  int k;                   // axis 0,1,2 which this node is dived into two    
  Float left, right;       // min & max coordinate of the particles
  index_t ibegin, iend;    // index range of particles in this node
  Float density_sum;
  Float left_h, right_h;   // min & max of density support (left - rk, right+rk)
};

//
// Class KNeighbors
//   store k smallest radius^2
class KNeighbors {
 public:
  KNeighbors(const int knbr);

  void clear();
  
  void push_back(const Float r2, const index_t idx) {
    // Insert r2 to the orderd vector v
    // This simple insersion algorithm assumes knbr is small 8 ~ 32
    // For large knbar, there are better way to organise the data
    // such as a binary tree
    if(r2 > v_r2.back()) return;
    
    int i= static_cast<int>(v_r2.size()) - 2;
    assert(i > 0);
    
    while(r2 < v_r2[i]) {
      assert(0 < i);
      assert(static_cast<size_t>(i + 1) < v_r2.size());
      v_r2[i + 1]= v_r2[i];
      v_idx[i + 1]= v_idx[i];
      --i;
    }

    assert(i >= 0);
    assert(static_cast<size_t>(i + 1) < v_r2.size());
    v_r2[i + 1]= r2;
    v_idx[i + 1]= idx;
  }

  Float r2(const size_t i) const {
    return v_r2[i];
  }

  index_t idx(const size_t i) const {
    return v_idx[i];
  }
  
  Float r2_max() const {
    return v_r2.back();
  }

 private:
  std::vector<Float> v_r2;
  std::vector<index_t> v_idx;
  // v_r2[0] = 0
  // v_r2[1] ... v[knbr]: k smallerst values so far
  // v_idx[1] ... v_idx[knbr]: index of k nearest neighbors
};


class KDTree {
 public:
  KDTree(std::vector<KDPoint>& v_, const int quota_);
  ~KDTree();

  void compute_rk(const int knbr);

  void compute_bounding_box(const size_t ibegin, const size_t iend,
			    Float left[], Float right[]);
  void update_node_statistics();
  //void adaptive_kernel_density(std::vector<KDPoint>& v, const int knbr);
  //Float adaptive_kernel_density(const Float x[]);
  Float adaptive_kernel_density(const Float x[], size_t inode=0);
  Float adaptive_average_density(const Float x[], const int knbr, Float r_guess);
 private:
  std::vector<KDPoint>& v;
  index_t quota; // number of particles in the leaf (bottom node)
  int knbr;      // number of neighbors for density estimation

  Node* nodes;
  size_t n_nodes;
  int height_max;

  Float left3_root[3], right3_root[3];

  void construct_balanced();
  void construct_balanced_recursive(const size_t inode, const int height,
				 const index_t ibegin, const index_t iend,
				 Float left[], Float right[], Float boxsize3[]);
  void collect_k_nearest_neighbors_recursive(const size_t inode,
					     const Float x[],
					     KNeighbors& nbr);
    
  //Float compute_density_sum_recursive(const size_t inode);
  Float update_node_statistice_recursive(const size_t inode,
					 Float left3[], Float right3[]);

  void average_density_recursive(const size_t inode,
				 const Float x[], const Float r, Float dx[],
				 size_t* const n, Float* sum);

  Float estimate_approx_density(const Float x[], const int knbr);
  Float estimate_approx_density_recursive(const size_t inode,
				   const Float x[],
				   const int knbr,
				   Float boxsize3[]);


};


//
// Templates
//

#endif
