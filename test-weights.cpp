#include "weights.h"


/////////////////////hypothetical test code////////////////
weights w(Nx*Ny*Nz*23);
w.update(0, 0.3);
w.update(1, 0.6);
assert(w.get_total() == 0.9);
assert(w.lookup(0.2) == 0);
//////////////////////////////////////////////////////
