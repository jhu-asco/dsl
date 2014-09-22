#ifndef DSL_STATE_H
#define DSL_STATE_H

#include <iostream>
#include <iomanip>
#include <string>
#include "system.h"


namespace dsl {
  class System;

  /**
   * Basic control system state tuple \f$(x,u,t)\f$ where 
   * \f$x \in X \f$ is the state of the mechanical system with state manifold \f$X\f$,
   * \f$\operatorname{dim}(X)=n \f$,
   * \f$u \in U \f$ are the controls with control manifold \f$U\f$, 
   * \f$\operatorname{dim}(U)=c \f$, 
   * and \f$t\f$ is the time.
   * Note that technically the system state is \f$x\f$ and a more 
   * appropriate name for this class would be "state+control tuple" since it has
   * extra information. But for conciseness we call this class simply State.
   *              
   * In addition to state, controls, and time, there is an extra field
   * of type dgc::Data that gives the ability to store additional abstract
   * data.
   *
   * Author: Marin Kobilarov -- Copyright (C) 2006
   */
  class State {
  public:
    /**
     * Initialize a basic control system state
     * @param sys system
     * @param x state vector (optional)
     * @param u control vector (optional)
     * @param t time (optional)
     * @param data additional abstract data (optional - empty by default, if provided the contents of data will be _copied_ into this object's data field)
     */
    State(const System &sys,
          const double *x = 0,
          const double *u = 0,
          double t = 0);
    
    /**
     * Copy constructor
     * @param s state 
     */
    State(const State &s);

    /**
     * Initialize a basic control system state from an input stream
     * @param sys system
     * @param istr input stream
    */
    State(const System &sys, std::istream &istr);

    virtual ~State();

    /**
     * Copy assignement
     * @param s state
     * @return a copy of this state
     */
    State& operator=(const State &s);

    /**
     * Produce a copy of the object
     * @return copy of this state
     */
    virtual State* Clone() const;

    const System& sys;   ///< the control system

    double *x;           ///< state vector
    double *u;           ///< control vector
    double t;            ///< time

  private:
    friend std::ostream& operator<<(std::ostream &os, const State &s);
    friend std::istream& operator>>(std::istream &is, State &s);
  };

  /**
   * Output the state to a stream
   * @param os output stream
   * @param s state
   * @return the output stream
   */
  std::ostream& operator<<(std::ostream &os, const dsl::State &s);

  /**
   * Input the state from a stream
   * @param is input stream
   * @param s state
   * @return the input stream
   */
  std::istream& operator>>(std::istream &is, dsl::State &s);
};

#endif
