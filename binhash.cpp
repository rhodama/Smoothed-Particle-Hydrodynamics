#include <string.h>

#include "zmorton.hpp"
#include "binhash.hpp"

/*@q
 * ====================================================================
 */

/*@T
 * \subsection{Spatial hashing implementation}
 * 
 * In the current implementation, we assume [[HASH_DIM]] is $2^b$,
 * so that computing a bitwise of an integer with [[HASH_DIM]] extracts
 * the $b$ lowest-order bits.  We could make [[HASH_DIM]] be something
 * other than a power of two, but we would then need to compute an integer
 * modulus or something of that sort.
 * 
 *@c*/

#define HASH_MASK (HASH_DIM-1)

unsigned particle_bucket(particle_t* p, float h)
{
    unsigned ix = p->x[0]/h;
    unsigned iy = p->x[1]/h;
    unsigned iz = p->x[2]/h;
    return zm_encode(ix & HASH_MASK, iy & HASH_MASK, iz & HASH_MASK);
}


    unsigned particle_neighborhood(unsigned int * buckets, particle_t* p, float h)
{
    /* BEGIN TASK */
    unsigned ix = p->x[0] / h;
    unsigned iy = p->x[1] / h;
    unsigned iz = p->x[2] / h;
    unsigned count = 0;

    unsigned x, y, z;

    // dx = -1, dy = -1, dz = -1
    x = (ix - 1) & HASH_MASK;
    y = (iy - 1) & HASH_MASK;
    z = (iz - 1) & HASH_MASK;
    BUCKET(buckets, 0) = zm_encode(x, y, z);

    // dx = -1, dy = -1, dz = 0
    x = (ix - 1) & HASH_MASK;
    y = (iy - 1) & HASH_MASK;
    z = iz & HASH_MASK;
    BUCKET(buckets, 1) = zm_encode(x, y, z);

    // dx = -1, dy = -1, dz = 1
    x = (ix - 1) & HASH_MASK;
    y = (iy - 1) & HASH_MASK;
    z = (iz + 1) & HASH_MASK;
    BUCKET(buckets, 2) = zm_encode(x, y, z);

    // dx = -1, dy = 0, dz = -1
    x = (ix - 1) & HASH_MASK;
    y = iy & HASH_MASK;
    z = (iz - 1) & HASH_MASK;
    BUCKET(buckets, 3) = zm_encode(x, y, z);

    // dx = -1, dy = 0, dz = 0
    x = (ix - 1) & HASH_MASK;
    y = iy & HASH_MASK;
    z = iz & HASH_MASK;
    BUCKET(buckets, 4) = zm_encode(x, y, z);

    // dx = -1, dy = 0, dz = 1
    x = (ix - 1) & HASH_MASK;
    y = iy & HASH_MASK;
    z = (iz + 1) & HASH_MASK;
    BUCKET(buckets, 5) = zm_encode(x, y, z);

    // dx = -1, dy = 1, dz = -1
    x = (ix - 1) & HASH_MASK;
    y = (iy + 1) & HASH_MASK;
    z = (iz - 1) & HASH_MASK;
    BUCKET(buckets, 6) = zm_encode(x, y, z);

    // dx = -1, dy = 1, dz = 0
    x = (ix - 1) & HASH_MASK;
    y = (iy + 1) & HASH_MASK;
    z = iz & HASH_MASK;
    BUCKET(buckets, 7) = zm_encode(x, y, z);

    // dx = -1, dy = 1, dz = 1
    x = (ix - 1) & HASH_MASK;
    y = (iy + 1) & HASH_MASK;
    z = (iz + 1) & HASH_MASK;
    BUCKET(buckets, 8) = zm_encode(x, y, z);

    // dx = 0, dy = -1, dz = -1
    x = ix & HASH_MASK;
    y = (iy - 1) & HASH_MASK;
    z = (iz - 1) & HASH_MASK;
    BUCKET(buckets, 9) = zm_encode(x, y, z);

    // dx = 0, dy = -1, dz = 0
    x = ix & HASH_MASK;
    y = (iy - 1) & HASH_MASK;
    z = iz & HASH_MASK;
    BUCKET(buckets, 10) = zm_encode(x, y, z);

    // dx = 0, dy = -1, dz = 1
    x = ix & HASH_MASK;
    y = (iy - 1) & HASH_MASK;
    z = (iz + 1) & HASH_MASK;
    BUCKET(buckets, 11) = zm_encode(x, y, z);

    // dx = 0, dy = 0, dz = -1
    x = ix & HASH_MASK;
    y = iy & HASH_MASK;
    z = (iz - 1) & HASH_MASK;
    BUCKET(buckets, 12) = zm_encode(x, y, z);

    // dx = 0, dy = 0, dz = 0
    x = ix & HASH_MASK;
    y = iy & HASH_MASK;
    z = iz & HASH_MASK;
    BUCKET(buckets, 13) = zm_encode(x, y, z);

    // dx = 0, dy = 0, dz = 1
    x = ix & HASH_MASK;
    y = iy & HASH_MASK;
    z = (iz + 1) & HASH_MASK;
    BUCKET(buckets, 14) = zm_encode(x, y, z);

    // dx = 0, dy = 1, dz = -1
    x = ix & HASH_MASK;
    y = (iy + 1) & HASH_MASK;
    z = (iz - 1) & HASH_MASK;
    BUCKET(buckets, 15) = zm_encode(x, y, z);

    // dx = 0, dy = 1, dz = 0
    x = ix & HASH_MASK;
    y = (iy + 1) & HASH_MASK;
    z = iz & HASH_MASK;
    BUCKET(buckets, 16) = zm_encode(x, y, z);

    // dx = 0, dy = 1, dz = 1
    x = ix & HASH_MASK;
    y = (iy + 1) & HASH_MASK;
    z = (iz + 1) & HASH_MASK;
    BUCKET(buckets, 17) = zm_encode(x, y, z);

    // dx = 1, dy = -1, dz = -1
    x = (ix + 1) & HASH_MASK;
    y = (iy - 1) & HASH_MASK;
    z = (iz - 1) & HASH_MASK;
    BUCKET(buckets, 18) = zm_encode(x, y, z);

    // dx = 1, dy = -1, dz = 0
    x = (ix + 1) & HASH_MASK;
    y = (iy - 1) & HASH_MASK;
    z = iz & HASH_MASK;
    BUCKET(buckets, 19) = zm_encode(x, y, z);

    // dx = 1, dy = -1, dz = 1
    x = (ix + 1) & HASH_MASK;
    y = (iy - 1) & HASH_MASK;
    z = (iz + 1) & HASH_MASK;
   BUCKET(buckets, 20) = zm_encode(x, y, z);

    // dx = 1, dy = 0, dz = -1
    x = (ix + 1) & HASH_MASK;
    y = iy & HASH_MASK;
    z = (iz - 1) & HASH_MASK;
    BUCKET(buckets, 21)= zm_encode(x, y, z);

    // dx = 1, dy = 0, dz = 0
    x = (ix + 1) & HASH_MASK;
    y = iy & HASH_MASK;
    z = iz & HASH_MASK;
    BUCKET(buckets, 22) = zm_encode(x, y, z);

    // dx = 1, dy = 0, dz = 1
    x = (ix + 1) & HASH_MASK;
    y = iy & HASH_MASK;
    z = (iz + 1) & HASH_MASK;
    BUCKET(buckets, 23) = zm_encode(x, y, z);

    // dx = 1, dy = 1, dz = -1
    x = (ix + 1) & HASH_MASK;
    y = (iy + 1) & HASH_MASK;
    z = (iz - 1) & HASH_MASK;
    BUCKET(buckets, 24) = zm_encode(x, y, z);

    // dx = 1, dy = 1, dz = 0
    x = (ix + 1) & HASH_MASK;
    y = (iy + 1) & HASH_MASK;
    z = iz & HASH_MASK;
    BUCKET(buckets, 25) = zm_encode(x, y, z);

    // dx = 1, dy = 1, dz = 1
    x = (ix + 1) & HASH_MASK;
    y = (iy + 1) & HASH_MASK;
    z = (iz + 1) & HASH_MASK;
    BUCKET(buckets, 26) = zm_encode(x, y, z);
    count = 27;  // Since we're always processing all 27 neighbors

    return count;
    /* END TASK */

   
}

    



void hash_particles(sim_state_t* s, float h)
{
    /* BEGIN TASK */
    int n = s->n;
    particle_t* part = s->part;
    particle_t** hash = s->hash;

    // Clear the hash table
    for (int i = 0; i < HASH_SIZE; i++) {
        hash[i] = nullptr;
    }
 

#pragma omp parallel
{
    #pragma omp for nowait
    for (int i = 0; i < n; i++) {
        particle_t* p = &part[i];
        unsigned b = particle_bucket(p, h);
        p->next = hash[b];

        #pragma omp atomic write
        hash[b] = p;
    }
}



}
