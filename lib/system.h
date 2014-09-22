#ifndef DSL_SYSTEM_H
#define DSL_SYSTEM_H

namespace dsl {

  class State;

  /**
   * Basic control system with state manifold \f$X\f$,
   * \f$\operatorname{dim}(X)=n \f$,
   * and control manifold \f$U\f$, 
   * \f$\operatorname{dim}(U)=c \f$. The system should contain the basic
   * system model parameters such as dimensions, mass porperties, control
   * input vector fields, etc... 
   * This description would be useful in using systems in the form of an 
   * ODE \f$\dot x = f(x,u,t) \f$, or for discrete integrators such as
   * \f$ x_{k+1} = f(x_k, u_{0:k+1}, t_{0:k+1}) \f$, or implicit ones.
   *
   * The class serves as a base class for implementing specific mechanical systems.
   *
   * Author: Marin Kobilarov -- Copyright (C) 2006
   */
  class System {
  public:
    /**
     * Initialize a basic DMOC system
     * @param n state dimension
     * @param c control input dimension
     */
    System(int n,
           int c);

    virtual ~System();

    /**
     * Set bounds on the state and control variables,
     * e.g. box constraints, joint limits, control input bounds.
     * Pass NULL to any of the parameters to ignore those bounds.
     * @param xlb (n-dim. array) lower bound on the state
     * @param xub (n-dim. array) upper bound on the state
     * @param ulb (c-dim. array) lower bound on the controls
     * @param uub (c-dim. array) upper bound on the controls  
     */
    void SetBounds(const double *xlb,
                   const double *xub,
                   const double *ulb = 0,
                   const double *uub = 0);


    /**
     * Create a new state - subclasses can choose to overwrite this method
     * in order to create the appropriate state for their system. 
     * @param x state vector
     * @param u controls
     * @param t time
     */
    virtual State* Create(const double *x = 0,
                          const double *u = 0,
                          double t = 0) const;
    

    /**
     * A basic method to interpolate a state tuple
     * \f$s=(x,u,t)\in X\times U\times \mathbb{R}_{+}\f$ between
     * two given state tuples \f$s_a=(x_a,u_a,t_a)\f$
     * and \f$s_b=(x_b,u_b,t_b)\f$, for \f$t_a < t_b \f$.
     * For instance, in case of a vector (or coordinate)
     * space one can set \f$ s=(1-\alpha)s_a + \alpha s_b \f$, 
     * for \f$ \alpha \in [0,1]\f$. The method is virtual since each
     * system can provide its own way of interpolation. Often coordinate 
     * spaces such as ones composed of the Euclidean group require 
     * careful interpolation that respects the charts chosen. The base
     * implementation is for linear interpolation in vector spaces. Subclasses
     * that requires special or higher odrer interpolation can choose
     * to override this method.
     * @param s new state 
     * @param sa first state
     * @param sb second state
     * @param a an interpolating number in [0,1]     
     * @param time should time be interpolated too? (true by default)
     */
    virtual void Get(State &s, 
                     const State &sa, 
                     const State &sb, 
                     double a, 
                     bool time = true) const;
       
    int n;         ///< \f$ n=\operatorname{dim}(X) \f$ state space dimension
    int c;         ///< \f$ c=\operatorname{dim}(U) \f$ control input space dimension

    double *xlb;         ///< state lower bound
    double *xub;         ///< state upper bound
    double *ulb;         ///< controls lower bound
    double *uub;         ///< controls upper bound  
  };

};

#endif
