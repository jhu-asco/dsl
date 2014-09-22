#ifndef DSL_TRAJECTORY_H
#define DSL_TRAJECTORY_H

#include <iostream>
#include <cstring>
#include <map>
#include "system.h"
#include "state.h"

namespace dsl {

  //  class Trajectory;

  /**
   * Basic control system trajectory represented by a set of discrete nodes
   * 
   * This should serve as a base class for implementing trajectories for a particular
   * mechanical system described by the System class
   *
   * Author: Marin Kobilarov -- Copyright (C) 2006
   */
  class Trajectory {
  public:
    /**
     * Initialize a trajectory for the system sys
     * Once a trajectory is initialized with this constructor
     * one should call Init(sn,...) with the appropriate number of
     * discrete points desired. Otherwise the trajectory would contain
     * only one point (usually the origin)
     * @param sys mechanical system
     * @param sn number of discrete segments (optional)
     */
    Trajectory(const System& sys, int sn = 0);

    /**
     * Copy constructor
     * @param traj a trajectory
     */
    Trajectory(const Trajectory &traj);


    /**
     * Copy assignment
     * @param traj a trajectory 
     */
    Trajectory& operator=(const Trajectory &traj);

    
    /**
     * Initialize a trajectory from a serialized stream
     * @param sys mechanical system
     * @param istr stream
     */    
    Trajectory(const System& sys, std::istream& istr);

    
    virtual ~Trajectory();

    /**
     * Clone the trajectory
     * @return a copy of this trajectory
     */
    virtual Trajectory* Clone() const;


    /**
     * Resize the trajectory keeping any
     * previous states and if necessary adding empty states
     * at the end.
     * @param sn new number of discrete segments
     */ 
    void Resize(int sn);

    /**
     * Reverse the states+control alongs trajectory. keep same times
     */ 
    void Reverse();


    /**
     * Initialize the trajectory using sn segments
     * and optionally (if both si and sf are povided)
     * interpolate the states between
     * states si and sf. In fact, Init(sn) is equivalent to
     * Resize(sn).
     * @param sn number of segments
     * @param si start state (optional)
     * @param sf final state (optional)
     */
    void Init(int sn,
              const State *si = 0,
              const State *sf = 0);
      

    /**
     * Attach trajectory traj to this trajectory. This trajectory
     * is modified and its size becomes as explained below. This operation
     * is fairly efficient since it is based on raw memory manipulation.
     * @param traj trajectory to be attached
     * @param back if true attach to back of this trajectory, else to the front
     * @param time if true then adjust time along the newly added states by
     *        incrementing/decrementing each state's time based on the 
     *        last/first state's time in the current trajectory depending
     *        on whether attaching to back/front, respectively.
     * @param js if true then treat the end point of this trajectory and 
     *           the first point of the given trajectory traj as the same 
     *           points (with same time)--in this case the total
     *           number of segments would be this->sn + traj.sn,
     *           otherwise it is this->sn + traj.sn + 1
     */
    void Attach(const Trajectory &traj, 
                bool back = true,
                bool time = false,
                bool js = true);

    /**
     * Append a state to the end of this trajectory
     * @param s state to be attached
     */
    void Append(const State &s);



    /**
     * Clear all states along this trajectory and free
     * associated memory. 
     */
    void Clear(); 

    /**
     * Set the path times starting from t0 with timestep h
     * 
     * @param t0 start time
     * @param h time step
     */
    void SetTime(double t0, double h);

    /**
     * Set any additional trajectory parameters (usually these are 
     * useful for optimization purposes and 
     * are added at the end of the optimization vector, e.g. the
     * homotopy continuation method uses these as homotopy parameters
     * for deforming obstacles -- see Obstacle)
     * @param pn number of additional params
     * @param ps parameter values
     * @param pmin lower bound
     * @param pmax upper bound
     */
    void SetParams(int pn, 
                   const double *ps = 0, 
                   const double *pmin = 0, 
                   const double *pmax = 0);
    
    
    /**
     * Get the state at time t (t should already be set in the state s)
     * t should be in the valid range of trajectory times. This uses
     * the base interpolation define in System::Get. For now this method
     * assumes that all segments have equal time duration. 
     *
     * @param s state
     */
    void Get(State &s) const;

    
    /**
     * Get the index of trajectory segment in which time t falls. 
     * For now this method 
     * assumes that all segments have equal time duration.
     * @param t time
     * @return discrete trajectory segment where t falls (returns -1 if t
     * is out of bounds)
     */
    int Get(double t) const;


    /**
     * Get a subtrajectory given two starting times. Utilizes the System::Get
     * implementation.
     * @param traj to be filled in (the number of discrete segments traj.sn is used to determine the time step, so the discrete length of traj is not modified). Hence, one should pass traj with already existing states, and with traj->sn >0. If traj->sn=0 then just use GetState.
     * @param ti initial time (must be within this trajectory)
     * @param tf final time (must be within this trajectory)
     */
    void Get(Trajectory &traj, double ti, double tf) const;
    
    
    /**
     * Add a point interpolated at s inside every trajectory segment
     * new trajectory has 2*sn+1 points
     * @param a number in [0,1] indicating where inside the segment to put the new point
     */
    void Refine(double a = .5);
    
    /**
     * Resample the trajectory using sn new equaly spaced segments
     * @param sn new number of discrete segments
     */    
    void Resample(int sn);


    /**
     * Add mn new states to the trajectory after points 
     * with indices in mi and interpolated between states at
     * indices mi[i] and mi[i+1] using the number mu[i] in [0,1].
     * @param mn number of new states to insert in the trajectory
     * @param mi (mn-array) indices of old states after which the new point will be inserted
     * @param mu (mn-array) numbers in the range [0,1] indicating where in the
     * mi[i]-th segment the new state will be interpolated
     */
    void Modify(int mn, const int *mi, const double *mu);


    /**
     * Checks whether this trajectory contains time t
     * @param t given time
     * @return true if t is between the start and end times of the trajectory
     */
    bool IsValidTime(double t) const;
    

    const System &sys;    ///< mechanical system
    
    int sn;               ///< sn number of discrete segments, (sn+1) points
    State **states;       ///< array of trajectory state pointers

    bool ext;             ///< extended formulation (used for optimization purposes -- this flag is internally set when one requires variations in the exteneded stapce: space + time)
    
    int pn;               ///< number of additional parameters
    double *ps;           ///< other parameters
    double *pmin;         ///< ps lower bound
    double *pmax;         ///< ps upper bound

    /**
     * Log the trajectory to a file.
     * @param logName log filename
     * @param di log every di-th state (optional, default is 1)
     */
    virtual void Log(const char *logName, int di = 1);    

    /**
     * Log a state in simple human readable format
     * @param os output stream
     * @param di log every di-th state (optional, default is 1)
     */
    virtual void Log(std::ofstream& os, int di = 1);

    /**
     * Write a formatted number
     * @param os output stream
     * @param fw field width
     * @param a the number of write
     */
    static void Log(std::ofstream& os, int fw, double a);       


    double hmin;          ///< minimum timestep size (used internally for checking trajectory consistency) it is 1e-10 by default
    
  private:
    friend std::ostream& operator<<(std::ostream &os, const Trajectory &traj);
    friend std::istream& operator>>(std::istream &is, Trajectory &traj);
  };

  /**
   * Output the trajectory to a stream
   * @param os output stream
   * @param traj trajectory 
   * @return the output stream
   */
  std::ostream& operator<<(std::ostream &os, const dsl::Trajectory &traj);

  /**
   * Input the trajectory from a stream
   * @param is input stream
   * @param traj trajectory
   * @return the input stream
   */
  std::istream& operator>>(std::istream &is, dsl::Trajectory &traj);
};

#endif
