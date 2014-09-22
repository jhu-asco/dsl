#ifndef DSL_PARTICLEPRM_H
#define DSL_PARTICLEPRM_H

#include "state.h"
#include "trajectory.h"

namespace dsl {

  class ParticlePrm {
  public:
    /**
     * Initialize a particle prm using its system and sampling
     * bounds defined by the states slb and sub
     * @param sys the system
     * @param slb state sampling lower bound
     * @param sub state sampling upper bound
     */
    ParticlePrm(const State &slb,
                const State &sub);

    virtual ~ParticlePrm();
    
    bool Compute(Trajectory &traj, const State &si, const State &sf);

    double ds;   ///< arc-length step size

    double eps;  ///< min distance between different vertices    
   
    double tf;

    double h;
 
    /**
     * Helper method for optimal planning
     */
    static bool sint_opt(double ts[2], double as[2], 
                         double xi, double vi, 
                         double xf, double vf,
                         double *am, int n, double tf);
    
    /**
     * Helper method for optimal planning
     */
    static bool sint_a(double a[2], 
                       double xi, double vi, 
                       double xf, double vf, double tf);


    /**
     * Helper method for optimal planning
     */
    static double sint_profile(double tss[][2], double ass[][2], 
                               const State &sa, const State &sb);
      

  };
};

#endif
