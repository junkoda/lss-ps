#include <iostream>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cassert>
#include "config.h"
#include "error.h"
#include "mean_density_adaptive.h"

using namespace std;

static size_t count_debug= 0;

//
// local functions
//

// Cubic spline kernel: [Monaghan1992]
//   J. Monaghan, Smoothed Particle Hydrodynamics,
//   Annual Review of Astronomy and Astrophysics, 30 (1992), pp. 543-574.
static inline Float kernel(const Float r, const Float rk)
{
  const Float h= 0.5*rk;
  const Float fac= 1.0/(M_PI*h*h*h);
  const Float q= r/h;

  if(q > 2.0)
    return 0.0;
  else if(q > 1.0)
    return 0.25*fac*pow(2.0 - q, 3);

  return fac*(1.0 - 1.5*q*q*(1.0 - 0.5*q));
}


// Node index of the left child node
static inline size_t left_child(const size_t i)
{
  return (i << 1) + 1;
}

// Index of the right child node
static inline size_t right_child(const size_t i)
{
  return (i << 1) + 2;
}

// Index of a child (left for ichild=0, right for ichild=1
static inline size_t child(const size_t i, const size_t ichild)
{
  return (i << 1) + (1 + ichild);
}

// Comparison function object
// Compares the coordinate of kth axis of two points p.x[k] < q.x[k]
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


//
// Class KNeighbors
//
KNeighbors::KNeighbors(const int knbr)
{
  v_r2.resize(knbr + 1, numeric_limits<Float>::max());
  v_r2[0]= 0;
  
  v_idx.resize(knbr + 1, -1);
}

void KNeighbors::clear()
{
  v_r2.assign(v_r2.size(), numeric_limits<Float>::max());
  v_idx.assign(v_idx.size(), -1);

  v_r2[0]= 0;
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

static inline Float dist2(const Float x[], const Float y[])
{
  // 3D Euclidian distance with periodic boundary condition
  Float dx= x[0] - y[0];
  Float dy= x[1] - y[1];
  Float dz= x[2] - y[2];

  return dx*dx + dy*dy + dz*dz;
}

//
// KDTree class
//

KDTree::KDTree(vector<KDPoint>& v_, const int quota_) :
  v(v_), quota(quota_), nodes(nullptr), n_nodes(0), height_max(0)
{
  construct_balanced();
}

KDTree::~KDTree()
{
  free(nodes);
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

//
// KDTree private member functions
//
void KDTree::construct_balanced()
{
  // Construct a balanced kd tree calling construct_balanced_recursive
  assert(quota > 0);
  // Compute the size of the kdtree and allocate the memory for nodes
  const size_t np= v.size(); assert(np > 0);
  
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
  }
  else {	       
    msg_printf(msg_fatal, "Error: unable to allocate memory for %lu nodes, "
	       "%lu mbytes required\n", n_nodes, mem_size);
    throw MemoryError();
  }

  height_max= height_new;

  // Compute the bounding box for particles
  Float left[3], right[3];
  Float boxsize3[3];
  
  compute_bounding_box(0, v.size(), left, right);

  for(int k=0; k<3; ++k) {
    boxsize3[k]= right[k] - left[k];
    left3_root[k]= left[k];
    right3_root[k]= right[k];
  }

  construct_balanced_recursive(0, 0, 0, v.size(), left, right, boxsize3);

  nodes[0].k= 0;
  nodes[0].left= left[0];
  nodes[0].right= right[0];

  msg_printf(msg_debug, "KDtree constructed kdtree\n");
}


void KDTree::construct_balanced_recursive(const size_t inode, const int height,
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
    // Update bounding box based on the particles in this node
    compute_bounding_box(ibegin, iend, left, right);
    node->left= left[node->k];
    node->right= right[node->k];
    return;
  }

  // Divide particles into two along k axis
  const int k= cut_direction(boxsize3);

  boxsize3[k] /= 2;

  const index_t imid= ibegin + (iend - ibegin)/2;
  assert(ibegin < imid && imid < iend - 1);
  nth_element(v.begin() + ibegin, v.begin() + imid, v.begin() + iend,
	      CompPoints<KDPoint>(k));

  // construct left sub nodes (contains particles with smaller x[k])
  const index_t ileft= left_child(inode);
  assert(0 <= ileft && ileft < static_cast<index_t>(n_nodes));
  nodes[ileft].k= k;

  construct_balanced_recursive(ileft, height + 1, ibegin, imid, left, right, boxsize3);

  nodes[ileft].left= left[k];
  nodes[ileft].right= right[k];

  // construct right sub nodes (contains particles with larger x[k])
  Float left1[3], right1[3];
  const index_t iright= right_child(inode);
  assert(0 <= iright && ileft < iright);
  nodes[iright].k= k;
  construct_balanced_recursive(iright, height + 1, imid, iend, left1, right1, boxsize3);
  nodes[iright].left= left1[k];
  nodes[iright].right= right1[k];

  // Update left and right
  // Original left/right was half of the parent node, this is updated to a smaller bounding box
  // based on the children nodes
  left[0]= min(left[0], left1[0]);
  left[1]= min(left[1], left1[1]);
  left[2]= min(left[2], left1[2]);

  right[0]= max(right[0], right1[0]);
  right[1]= max(right[1], right1[1]);
  right[2]= max(right[2], right1[2]);
}

//
// Find neighbors recurlively
//
void KDTree::collect_k_nearest_neighbors_recursive(const size_t inode,
						   const Float x[],
						   KNeighbors& nbr)
{
  // Find k nearest neighbor particles near position x[]
  // and push it to the Nbr object

  count_debug++;

  assert(0 <= inode && inode < n_nodes);
  Node const * const node= nodes + inode;
  const int k= node->k;
  const Float eps2= nbr.r2_max();
  const Float eps= sqrt(eps2);

  // If this node is more than eps away from x, no k-nearest neighbors
  // are in this node
  if(x[k] < node->left  - eps || x[k] > node->right + eps) {
    return; // This node is far enough from position x
  }

  // Add neightbor particles if this node is a leaf
  if(node->iend - node->ibegin <= quota) {
    for(index_t j=node->ibegin; j<node->iend; ++j) {
      assert(0 <= j && j < static_cast<index_t>(v.size()));

      Float r2= dist2(x, v[j].x);
      if(r2 > 0) { // Adding the density of the particle itself makes
	           // the density systematically higher
	nbr.push_back(r2, j);
      }
    }
    return;
  }

  // Search subtrees recursively
  size_t iright= right_child(inode);
  int inear= x[nodes[iright].k] > nodes[iright].left;
  // inear = 0: search left node first
  // inear = 1: search right node first
  
  collect_k_nearest_neighbors_recursive(child(inode, inear), x, nbr);
  collect_k_nearest_neighbors_recursive(child(inode, !inear), x, nbr);
}

void KDTree::compute_rk(const int knbr)
{
  // Compute rk and n_local
  //
  //
  // Args:
  //   knbr (int): number of neighbors used for adaptive kernel
  //               density estimation
  //
  // Input:
  //   Particles positions v[i].x
  // Output:
  //   Distance to kth neighbor v[i].rk
  //   Estimated density v[i].n_local
  //
  const size_t n= v.size();

  // Reset n_local
  for(size_t i=0; i<n; ++i)
    v[i].n_local= 0;  

  // For each particle in v, search for k nearest neighbors
  // and assign the density of particle i to the neighbors
  KNeighbors nbr(knbr);
  for(size_t i=0; i<n; ++i) {
    nbr.clear();
    count_debug= 0;
    collect_k_nearest_neighbors_recursive(0, v[i].x, nbr);
    //cerr << "count: " << count_debug << endl;

    const Float rk= sqrt(nbr.r2_max());
    const Float w= v[i].w;
    v[i].rk= rk;

    for(int j=1; j<=knbr; ++j) {
      index_t inbr = nbr.idx(j);
      assert(0 <= inbr && static_cast<size_t>(inbr) < v.size());

      v[inbr].n_local += w*kernel(sqrt(nbr.r2(j)), rk);
    }
  }
  msg_printf(msg_debug, "KDTree rk and local density computed\n");
}


Float KDTree::adaptive_kernel_density(const Float x[], size_t inode)
{
  //
  // Args:
  //   x[3]:  poisiton for the density estiamtion
  //   inode: compute density within node i (0 for first call)
  //
  Node const * const node= nodes + inode;
  const int k= node->k;

  // The resulting densityg
  Float dens= 0;
  
  if(x[k] < node->left_h || x[k] > node->right_h) {
    // This node does not contribute to the density at x
    return 0;
  }

  if(node->iend - node->ibegin <= quota) {
    // This is a leaf; compute density from particles in this node
    for(index_t j=node->ibegin; j<node->iend; ++j) {
      assert(0 <= j && j < static_cast<index_t>(v.size()));
      Float r2= dist2(x, v[j].x);
      if(r2 > 0) {
	dens += v[j].w*kernel(sqrt(r2), v[j].rk);
      }
    }
    return dens;
  }

  dens += adaptive_kernel_density(x, left_child(inode));
  dens += adaptive_kernel_density(x, right_child(inode));

  return dens;
}


void KDTree::update_node_statistics()
{
  Float left3[3], right3[3];
  for(int k=0; k<3; ++k) {
    left3[k]= left3_root[k];
    right3[k]= right3_root[k];
  }
  
  update_node_statistice_recursive(0, left3, right3);
}

Float KDTree::update_node_statistice_recursive(const size_t inode,
				     Float left3[], Float right3[])
{
  //
  // Output
  //   density_sum: sum of v[i].n_local in this node
  //   left_h:      minimum of v[i].x[k] - v[i].rk in this node
  //   right_h:     maximum of v[i].x[k] + v[i].rk in this node
  //                where k= node->k
  //
  Float sum= 0;
  
  assert(0 <= inode && inode < n_nodes);
  Node* const node= nodes + inode;

  const int k= node->k;
  left3[k]= node->left;
  right3[k]= node->right;

  // Sum n_local  if this node is a leaf
  if(node->iend - node->ibegin <= quota) {
    node->left_h= left3[k];
    node->right_h= right3[k];
    
    for(index_t j=node->ibegin; j<node->iend; ++j) {
      assert(0 <= j && j < static_cast<index_t>(v.size()));
      node->left_h= min(node->left_h, node->left - v[j].rk);
      node->right_h= max(node->right_h, node->right + v[j].rk);

      sum += v[j].n_local;
    }
    node->density_sum= sum;
    return sum;
  }

  // Traverse subnodes recursively
  sum += update_node_statistice_recursive(left_child(inode), left3, right3);

  Float left3b[]= {left3[0], left3[1], left3[2]};
  Float right3b[]= {right3[0], right3[1], right3[2]};
  sum += update_node_statistice_recursive(right_child(inode), left3b, right3b);

  // Update using subnode information
  Node* const left_node= nodes + left_child(inode);
  Node* const right_node= nodes + right_child(inode);

  node->left_h= min(left_node->left_h, right_node->left_h);
  node->right_h= min(left_node->right_h, right_node->right_h);

  node->density_sum= sum;
  
  return sum;
}

void KDTree::average_density_recursive(const size_t inode,
				       const Float x[], const Float r,
				       Float dx[],
				       size_t* const n, Float* sum)
{
  // Compute the average density of neighbors within radius R about x
  // Args:
  //   left3: minum of particle positions
  //   right3: maximum of particle positions
  // Output:
  //   *n: number of particles within the sphere
  //   *sum: sum of local densities (n_local) within the sphere
  assert(0 <= inode && inode < n_nodes);
  Node* const node= nodes + inode;
  const int k= node->k;

  count_debug++;


  // If this node is more than eps away from x, no k-nearest neighbors
  // are in this node
  if(x[k] < node->left  - r || x[k] > node->right + r) {
    return; // This node is far enough from position x
  }

  dx[k]= max(abs(x[k] - node->left), abs(x[k] - node->right));
  if(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] < r*r) {
    // distance to the furthest corner is less than r
    // => node is completely within radius r

    *sum += node->density_sum;
    *n   += node->iend - node->ibegin;

    return;
  }

  Float density_sum= 0;
  size_t num= 0;
  
  // Check all particles if this node is a leaf
  if(node->iend - node->ibegin <= quota) {
    for(index_t j=node->ibegin; j<node->iend; ++j) {
      assert(0 <= j && j < static_cast<index_t>(v.size()));

      Float rj2= dist2(x, v[j].x);

      
      if(rj2 <= r*r) {
	density_sum += v[j].n_local;
	num++;
      }
    }
    
    *sum += density_sum;
    *n   += num;
    return;
  }

  // Traverse subtrees recursively
  Float dxa[]= {dx[0], dx[1], dx[2]};
  average_density_recursive(left_child(inode), x, r, dxa,
			    n, sum);
  
  average_density_recursive(right_child(inode), x, r, dx,
			    n, sum);
}

Float KDTree::estimate_approx_density(const Float x[], const int knbr)
{
  // Roughly estimate the density around x containing knbr particles
  Float boxsize3[]= {right3_root[0] - left3_root[0],
		     right3_root[1] - left3_root[1],
		     right3_root[2] - left3_root[2]};


  return estimate_approx_density_recursive(0, x, knbr, boxsize3);		     
}

Float KDTree::estimate_approx_density_recursive(const size_t inode,
				      const Float x[],
				      const int knbr,
				      Float boxsize3[])
{
  // Estimate density on scales containing knbr particles
  // Returns node cuboid volume / number of particles
  // when number of particles first gets below knbr
  
  Node const * const node= nodes + inode;
  boxsize3[node->k]= node->right - node->left;
  const Float vol= boxsize3[0]*boxsize3[1]*boxsize3[2];
  
  if(node->iend - node->ibegin < knbr) {
    // This is the scale containg about knbr particles
    return (node->iend - node->ibegin)/vol;
  }

  const size_t ileft= left_child(inode);
  const int k= nodes[ileft].k;
  if(nodes[ileft].left <= x[k] && x[k] <= nodes[ileft].right)
    return estimate_approx_density_recursive(ileft, x, knbr, boxsize3);

  const size_t iright= left_child(inode);
  assert(k == nodes[iright].k);
  if(nodes[iright].left <= x[k] && x[k] <= nodes[iright].right)
    return estimate_approx_density_recursive(iright, x, knbr, boxsize3);

  return (node->iend - node->ibegin)/vol;
}

Float KDTree::adaptive_average_density(const Float x[], const int knbr,
				       Float r_guess)
{
  // Average n_local of k neighbor particles around position x
  //
  // Args:
  //    x[3]: position
  //    knbr: approximate number of neighbor, k, that the local density will be averaged
  //    r_guess: initial guess of the radius containing k neighbors (0 for automatic guess)
  count_debug= 0;

  if(r_guess == 0) {
    // estiamte r_guess from density
    Float dens_guess= estimate_approx_density(x, knbr);
    r_guess = pow(knbr/(4.0/3.0*M_PI*dens_guess), 1.0/3.0);
  }
  
  Float dx3[3];
  for(int k=0; k<3; ++k) {
    dx3[k]= max(abs(right3_root[k] - x[k]), abs(x[k] - left3_root[k]));
  }

  size_t n_particles;
  Float density_sum;

  const int iter_max= 2;

  for(int iter=0; iter<iter_max; ++iter) {
    n_particles= 0;
    density_sum= 0;
    average_density_recursive(0, x, r_guess, dx3,
			      &n_particles, &density_sum);

    if(n_particles == 0) {
      r_guess= 10*r_guess;
    }
    else {
      Float r_guess_new=
	r_guess*pow(static_cast<double>(knbr)/n_particles, 1.0/3);
      //cerr << r_guess << " " << density_sum << " / " << n_particles << endl;
      //cerr << iter << "/" << iter_max << " " <<
      //n_particles << " " << r_guess << " -> " << r_guess_new << endl;

      r_guess= r_guess_new;
    }
  }

  //cerr << n_particles << endl;
  //abort();

  //cerr << "node visit " << count_debug << endl;

  return density_sum/n_particles;

}

