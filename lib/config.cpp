#include <assert.h>
#include "config.h"
#include "msg.h"


void config_assert(void)
{
  assert(ALGN % sizeof(Float) == 0);
  assert(ALGN % sizeof(double) == 0);
}

size_t size_align(size_t size)
{
  if(size % ALGN != 0)
    size += ALGN - (size % ALGN);

  assert(size % ALGN == 0);
  
  return size;
}


