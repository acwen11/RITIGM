#ifndef SEPARATION_H
#define SEPARATION_H

#include "pizzacactus.h"

namespace Pizza {
namespace BNSTrackerGen {

pz_real proper_distance(const cactus_glop& glop,
                        const pz_real x0, const pz_real y0, 
                        const pz_real x1, const pz_real y1, int num_points);

                        
}
}

#endif
 
