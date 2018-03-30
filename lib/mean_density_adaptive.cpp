#include <iostream>
#include <limits>
#include <algorithm>
#include <cassert>
#include "config.h"
#include "error.h"
#include "mean_density_adaptive.h"

using namespace std;

//
// local functions
//

// Indexing of a binary tree
static inline size_t left_child(const size_t i)
{
  return (i << 1) + 1;
}

static inline size_t right_child(const size_t i)
{
  return (i << 1) + 2;
}

/*
static inline bool is_leaf(KdTree const * const tree)
{
  return (tree->iend - tree->ibegin) < KdTree::quota;
}
*/

// Comparison function object
template<typename T> class CompPoints {
 public:
  CompPoints(const int direction) : k(direction) {}
  bool operator()(const T& p, const T& q) const {
    // compare the corrdinate of k'th axis
    return p.x[k] < q.x[k];
  }
 private:
  int k;
};


template <typename float_type>
void construct_vector(float_type const * const xyz,
		      const size_t xyz_stride,
		      float_type const * weight,
		      const size_t weight_stride,
		      const size_t np)
{
  vector<KDPoint> v;
  v.reserve(np);

  KDPoint p;
  p.w = 0;
  p.n_local = 0;
  p.n_average = 0;
  
  for(size_t i=0; i<np; ++i) {
    p.x[0] = xyz[0];
    p.x[1] = xyz[1];
    p.x[2] = xyz[2];
    xyz= (float_type*) ((char*) xyz + xyz_stride);
    
    if(weight) {
      p.w = *weight;
      weight = (float_type*) ((char*) weight + weight_stride);
    }
    v.push_back(p);
  }
}


void KDTree::compute_bounding_box(const size_t ibegin,
				  const size_t iend,
				  Float left[], Float right[])
{
  // Compute the bounding box, the cuboid that encloses all particles
  // v[ibegin] to v[iend - 1]
  assert(v.size() > 0);
  
  left[0]= right[0]= v[ibegin].x[0];
  left[1]= right[1]= v[ibegin].x[1];
  left[2]= right[2]= v[ibegin].x[2];

  for(size_t i=ibegin+1; i<iend; ++i) {
    left[0]= min(left[0], v[i].x[0]);
    left[1]= min(left[1], v[i].x[1]);
    left[2]= min(left[2], v[i].x[2]);

    right[0]= max(right[0], v[i].x[0]);
    right[1]= max(right[1], v[i].x[1]);
    right[2]= max(right[2], v[i].x[2]);
  }
}

static inline int cut_direction(const Float boxsize3[])
{
  // Find the longest edge of a cuboid
  int k=1;
  if(boxsize3[2] > boxsize3[1])
    k=2;
  if(boxsize3[0] > boxsize3[k])
    k=0;

  return k;
}



//
// KDTree class
//

KDTree::KDTree(vector<KDPoint>& v_, const int quota_) :
  v(v_), quota(quota_), nodes(nullptr), n_nodes(0)
{

}

KDTree::~KDTree()
{

}

//
// KDTree private member functions
//
void KDTree::construct_balanced()
{
  // Construct a balanced kd tree calling construct_balanced_recursive
  assert(quota > 0);
  // Compute the size of the kdtree and allocate the memory for nodes
  const size_t np= v.size();
  size_t nleaf= 1;
  size_t height_new= 0;
  while(quota*nleaf < np) {
    nleaf= nleaf << 1;
    height_new++;
  }

  assert(nleaf == (static_cast<size_t>(1) << height_new));

  n_nodes= 2*nleaf - 1;

  nodes= (Node*) malloc(sizeof(Node)*n_nodes);
  size_t mem_size= sizeof(Node)*n_nodes/(1000*1000); // in Mbytes
  
  if(nodes) {
    msg_printf(msg_info, "Allocated %lu nodes, %lu Mbytes for kdtree\n",
	       n_nodes, mem_size);
    //msg_printf(msg_verbose, "Tree depth= %d, number of nodes= %lu\n",
    //height_new, n_nodes);
  }
  else {	       
    msg_printf(msg_fatal, "Error: unable to allocate memory for %lu nodes, "
	       "%lu mbytes required\n", n_nodes, mem_size);
    throw MemoryError();
  }

  // Compute the bounding box for particles
  Float left[3], right[3];
  Float boxsize3[3];
  
  compute_bounding_box(0, v.size(), left, right);
  boxsize3[0]= right[0] - left[0];
  boxsize3[1]= right[1] - left[1];
  boxsize3[2]= right[2] - left[2];

  construct_balanced_recursive(0, 0, v.size(), left, right, boxsize3);

  nodes[0].k= 0;
  nodes[0].left= left[0];
  nodes[0].right= right[0];
}


void KDTree::construct_balanced_recursive(const size_t inode, 
				  const index_t ibegin, const index_t iend,
				  Float left[], Float right[], Float boxsize3[])
{
  // Construct a balanced kd tree recursively.
  // First called from construct_balanced()
  //
  // Args:
  //   inode: index of the node constucting
  //   ibegin-iend: particle index range within this node
  // Output:
  //   left[]: left corner of this node
  //   right[]: right corner of this node
  assert(0 <= inode && inode < n_nodes);
  Node* const node= nodes + inode;

  node->ibegin= ibegin;
  node->iend= iend;
  
  if(iend - ibegin <= quota) {
    compute_bounding_box(ibegin, iend, left, right);
    node->left= left[node->k];
    node->right= right[node->k];
    return;
  }

  // Divide particles into two along k axis
  const int k= cut_direction(boxsize3);

  boxsize3[k] /= 2;

  const size_t imid= ibegin + (iend - ibegin)/2;
  nth_element(v.begin() + ibegin, v.begin() + imid, v.begin() + iend,
	      CompPoints<KDPoint>(k));

  // construct left sub nodes (contains particles with smaller x[k])
  const size_t ileft= left_child(inode);
  assert(0 <= ileft && ileft < n_nodes);
  node[ileft].k= k;
  construct_balanced_recursive(ileft, ibegin, imid, left, right, boxsize3);
  node[ileft].left= left[k];
  node[ileft].right= right[k];

  // construct right sub nodes (contains particles with larger x[k])
  Float left1[3], right1[3];
  const size_t iright= right_child(inode);
  assert(0 <= iright && ileft < iright);
  node[iright].k= k;
  construct_balanced_recursive(iright, imid, iend, left1, right1, boxsize3);
  node[iright].left= left1[k];
  node[iright].right= right1[k];

  // Compute left and right with the values returned from children nodes
  left[0]= min(left[0], left1[0]);
  left[1]= min(left[1], left1[1]);
  left[2]= min(left[2], left1[2]);

  right[0]= max(right[0], right1[0]);
  right[1]= max(right[1], right1[1]);
  right[2]= max(right[2], right1[2]);
}


