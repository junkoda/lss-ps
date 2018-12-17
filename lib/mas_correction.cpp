#include "mas_correction.h"
#include "msg.h"

using namespace std;

static vector<Float> mas_correction_array;
static int mas_correction_n;

vector<Float>&
mas_correction_get_array(const int nc, const int n_mas)
{
  if(mas_correction_array.size() == static_cast<size_t>(nc)
     && n_mas == mas_correction_n)
    return mas_correction_array;

  mas_correction_array.clear();
  mas_correction_array.assign(nc, static_cast<Float>(1));

  if(n_mas == 0) // mas_correction = false
    return mas_correction_array;
  
  const int knq = nc/2;
  const Float fac= M_PI/nc;

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for(int i=1; i<nc; ++i) {
    int k= i <= knq ? i : i - nc;
    Float sinc = sin(fac*k)/(fac*k);
    mas_correction_array[i] = 1.0/pow(sinc, 2*n_mas);
  }

  msg_printf(msg_info, "MAS correction array initialised\n");
  
  return mas_correction_array;
}
