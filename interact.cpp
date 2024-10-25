#include <omp.h>  
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "vec3.hpp"
#include "zmorton.hpp"

#include "params.hpp"
#include "state.hpp"
#include "interact.hpp"
#include "binhash.hpp"
#include <immintrin.h>
#define HASH_MASK (HASH_DIM-1)

/* Define this to use the bucketing version of the code */
#define USE_BUCKETING


#define CACHE_LINE_SIZE 64

/*@T
 * \subsection{Density computations}
 * 
 * The formula for density is
 * \[
 *   \rho_i = \sum_j m_j W_{p6}(r_i-r_j,h)
 *          = \frac{315 m}{64 \pi h^9} \sum_{j \in N_i} (h^2 - r^2)^3.
 * \]
 * We search for neighbors of node $i$ by checking every particle,
 * which is not very efficient.  We do at least take advange of
 * the symmetry of the update ($i$ contributes to $j$ in the same
 * way that $j$ contributes to $i$).
 *@c*/

typedef struct {
    unsigned value;
    char padding[CACHE_LINE_SIZE - sizeof(unsigned)];
} aligned_unsigned;
static_assert(sizeof(aligned_unsigned) == CACHE_LINE_SIZE, "aligned_unsigned size is not equal to cache line size");


inline
void update_density(particle_t* pi, particle_t* pj, float h2, float C)
{
    float r2 = vec3_dist2(pi->x, pj->x);
    float z  = h2-r2;
    if (z > 0) {
        float rho_ij = C*z*z*z;
        #pragma omp atomic
        pi->rho += rho_ij;
        #pragma omp atomic
        pj->rho += rho_ij;
    }
}

void compute_density(sim_state_t* s, sim_param_t* params)
{
    int n = s->n;
    particle_t* p = s->part;
    particle_t** hash = s->hash;

    float h  = params->h;
    float h2 = h*h;
    float h3 = h2*h;
    float h9 = h3*h3*h3;
    float mass = s->mass;
    float C  = ( 315.0/64.0/M_PI ) * mass / h9;
    

    // Clear densities

#pragma omp parallel
{
    #pragma omp for schedule(static) 
    for (int i = 0; i < n; ++i) {
        p[i].rho = 0;
    }
}


    // Accumulate density info
#ifdef USE_BUCKETING
    /* BEGIN TASK */
    

float C2=(315.0/64.0/M_PI) * mass / h3;
#pragma omp parallel
{   
    float local_rho[n] = {0};
    unsigned int buckets[27 * INTS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    #pragma omp for schedule (dynamic, 64) 
    for (int i = 0; i < n; ++i) {
        particle_t* pi = p + i;
        pi->rho += C2;
        
        int np = particle_neighborhood(buckets, pi, h);
        
        for (int k = 0; k < 27; ++k) {
            particle_t* pj = s->hash[BUCKET(buckets, k)];
            if (pj) {
                __builtin_prefetch(pj->next, 0, 1);  
            for (; pj; pj = pj->next) {
                if (pj > pi) {
                    float r2 = vec3_dist2(pi->x, pj->x);
                    float z = h2 - r2;
                    if (z > 0) {
                        float rho_ij = C * z * z * z;

                        local_rho[i] += rho_ij;
                        local_rho[pj - p] += rho_ij;  
                    }
                }
                if (pj->next) {
                    __builtin_prefetch(pj->next->next, 0, 1);  
            }
        }
    }
}
}


    #pragma omp critical
        {
            for (int i = 0; i < n; ++i) {
                s->part[i].rho += local_rho[i];
            }
        }


}




    /* END TASK */
#else
    for (int i = 0; i < n; ++i) {
        particle_t* pi = s->part+i;
        pi->rho += ( 315.0/64.0/M_PI ) * s->mass / h3;
        for (int j = i+1; j < n; ++j) {
            particle_t* pj = s->part+j;
            update_density(pi, pj, h2, C);
        }
    }
#endif
}


/*@T
 * \subsection{Computing forces}
 * 
 * The acceleration is computed by the rule
 * \[
 *   \bfa_i = \frac{1}{\rho_i} \sum_{j \in N_i} 
 *     \bff_{ij}^{\mathrm{interact}} + \bfg,
 * \]
 * where the pair interaction formula is as previously described.
 * Like [[compute_density]], the [[compute_accel]] routine takes
 * advantage of the symmetry of the interaction forces
 * ($\bff_{ij}^{\mathrm{interact}} = -\bff_{ji}^{\mathrm{interact}}$)
 * but it does a very expensive brute force search for neighbors.
 *@c*/

inline
void update_forces(particle_t* pi, particle_t* pj, float h2,
                   float rho0, float C0, float Cp, float Cv)
{
    float dx[3];
    vec3_diff(dx, pi->x, pj->x);
    float r2 = vec3_len2(dx);
    if (r2 < h2) {
        const float rhoi = pi->rho;
        const float rhoj = pj->rho;
        float q = sqrt(r2/h2);
        float u = 1-q;
        float w0 = C0 * u/rhoi/rhoj;
        float wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q;
        float wv = w0 * Cv;
        float dv[3];
        vec3_diff(dv, pi->v, pj->v);

        // Equal and opposite pressure forces
        vec3_saxpy(pi->a,  wp, dx);
        vec3_saxpy(pj->a, -wp, dx);
        
        // Equal and opposite viscosity forces
        vec3_saxpy(pi->a,  wv, dv);
        vec3_saxpy(pj->a, -wv, dv);
    }
}

void compute_accel(sim_state_t* state, sim_param_t* params)
{
    // Unpack basic parameters
    const float h    = params->h;
    const float rho0 = params->rho0;
    const float k    = params->k;
    const float mu   = params->mu;
    const float g    = params->g;
    const float mass = state->mass;
    const float h2   = h*h;

    double t1, t2;
    // Unpack system state
    particle_t* p = state->part;
    particle_t** hash = state->hash;
    int n = state->n;

    // 1. Rehash particles
   
    hash_particles(state, h);
  
    

    // 2. Compute density
 
    compute_density(state, params);
    
    

    // Start with gravity and surface forces
    // 3. Apply gravity
   
    for (int i = 0; i < n; ++i) {
        vec3_set(p[i].a, 0, -g, 0);
    }


    // Constants for interaction term
    float C0 = 45 * mass / M_PI / ( (h2)*(h2)*h );
    float Cp = k/2;
    float Cv = -mu;
    // 4. Compute interaction forces

    // Accumulate forces
#ifdef USE_BUCKETING
     /* BEGIN TASK */
    
    
 #pragma omp parallel
{
      
      unsigned int buckets[27 * INTS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));

    float local_a[n][3] = {0};  
 

    #pragma omp for schedule(dynamic, 64) 
    for (int i = 0; i < n; ++i) {
        particle_t* pi = p + i;
        int np = particle_neighborhood(buckets, pi, h);

 
        float pi_x[3] = {pi->x[0], pi->x[1], pi->x[2]};
        float pi_v[3] = {pi->v[0], pi->v[1], pi->v[2]};
        float pi_rho = pi->rho;

        for (int k = 0; k < 27; ++k) {
            for (particle_t* pj = hash[BUCKET(buckets, k)]; pj; pj = pj->next) {
                __builtin_prefetch(pj->next, 0, 1);
                
                if (pj > pi) {
                    float dx[3], dv[3];
                    float r2 = 0;

                    
                    dx[0] = pi_x[0] - pj->x[0];
                    dx[1] = pi_x[1] - pj->x[1];
                    dx[2] = pi_x[2] - pj->x[2];
                    dv[0] = pi_v[0] - pj->v[0];
                    dv[1] = pi_v[1] - pj->v[1];
                    dv[2] = pi_v[2] - pj->v[2];
                    r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

                    if (r2 < h2) {
                        float rhoj = pj->rho;
                        float q = sqrt(r2/h2);
                        float u = 1-q;
                        float w0 = C0 * u / (pi_rho * rhoj);
                        float wp = w0 * Cp * (pi_rho + rhoj - 2*rho0) * u/q;
                        float wv = w0 * Cv;

                       
                        float force0 = wp * dx[0] + wv * dv[0];
                        float force1 = wp * dx[1] + wv * dv[1];
                        float force2 = wp * dx[2] + wv * dv[2];
                        int j = pj - p;
                        if (j >= 0 && j < n) {  
                    
                        local_a[i][0] += force0;
                   
                        local_a[i][1] += force1;
                 
                        local_a[i][2] += force2;
                     
                        local_a[j][0] -= force0;
                     
                        local_a[j][1] -= force1;
                       
                        local_a[j][2] -= force2;
                    }
                    }
                    
                }
                if (pj->next) {
                    __builtin_prefetch(pj->next->next, 0, 1);
                }
            }
        }

    }   


#pragma omp critical
          {
    for (int i = 0; i < n; ++i) {
                p[i].a[0] += local_a[i][0];
                p[i].a[1] += local_a[i][1];
                p[i].a[2] += local_a[i][2];
    }
          }
          
}


    /* END TASK */
#else
    for (int i = 0; i < n; ++i) {
        particle_t* pi = p+i;
        for (int j = i+1; j < n; ++j) {
            particle_t* pj = p+j;
            update_forces(pi, pj, h2, rho0, C0, Cp, Cv);
        }
    }
#endif

}

