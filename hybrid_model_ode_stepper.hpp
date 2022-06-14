#ifndef HYBRID_MODEL_ODE_STEPPER_HPP
#define HYBRID_MODEL_ODE_STEPPER_HPP


#include<iostream>
#include<cmath>
#include<string>
#include<boost/numeric/odeint.hpp>
#include<boost/random.hpp>
#include<cstdlib>
#include<ctime>

#include "hybrid_model_ode_system.hpp"

namespace hybrid_model {

using namespace std;
using namespace boost::numeric::odeint;

template < typename State >
class HybridModelOdeStepper {
 public:
	// Renames types
	typedef State state_type;
	typedef state_type::value_type value_type;
	typedef value_type time_type;
	typedef unsigned short order_type;

	// Defines the stepper_tag so that integrate_time() can choose which
  // integration method to use, e.g. `stepper_tag` is for simple stepper,
  // 'controlled_stepper_tag` is for controlled stepper.
	typedef boost::numeric::odeint::stepper_tag stepper_category;

	// Constructor
	HybridModelOdeStepper(int max_interation, double time_error_tolerance, double x_error_tolerance, const string& time_estimate_method): max_interation_for_binary_search_(max_interation),  time_error_tolerance_(time_error_tolerance),  x_error_tolerance_(x_error_tolerance), time_estimate_method_(time_estimate_method) {}

  // The method is called by integrate-function to do each ODE step estimate.
  template<Class System>
  void do_step(System system, state_type &x, time_type t, time_type dt) const;

 private:
  // The method to estimate each step's change.
	template<class System>
	void ApplyRungeKutta4(System system, state_type &x, time_type t, time_type dt) const;

	// Get the estimated time for selected @x from given bounds.
  // @estimate_method can be "binary_search" or "linear_interpolation".
	time_type EstimateTime(double x, time_type t_lower, time_type t_upper, double x_lower, double x_upper, const string& time_estimate_method) const; 

	// class members
	int max_interation_for_binary_search_;
	double time_error_tolerance_;
	double x_error_tolerance_; 
  sting time_estimate_method_;
}


template < typename State >
template < class System >
void HybridModelOdeStepper<State>::do_step(System system, state_type &x, time_type t, time_type dt) const {
	//========================
  //  Initialize variables
  //========================
  state_type x0, x1;                      // temporary state variables.
	time_type t_start = t, t_end = t + dt;  // start and end time point
	time_type t_cur = t;                    // current time point
			    	
  //===========================================
  // Case 0. if all cells are dead, do nothing
  //===========================================
  if (x[m*N+2] == 0.0) {
    return;
  }
        
  //=================================================================
	// This method helps find the exact location of flipping 
	// the genes and also consider the case that multiple flippings 
	// which happens in the same time block exist.
	//=================================================================
	while (t_cur < t_end) {
		// Updates dt, which is the time left from current time point to the end time of this step.
		dt = t_end - t_cur;
	  
    // Updates the system's status by checking the integral value.
		state_type x_temp;
    system.second(x, x_temp, t_cur);
		x = x_temp;
    
    // if all cells are dead
    if (x[m*N+2] == 0.0) {
      return;
    }

					// Processes for getting new x states
					x0 = x;	
					fTarg<System> (system, x0, t_cur, dt);

					// This case means no more flips in this time block now, 
					// so continue to next one
					if (x0[N*m] <= x0[N*m+1]) {
						x = x0;
						break;
					} 

					//========================================================
					//  Using bracketing method to find the exact location 
					//  of flipping
					//=========================================================

					// time and y setup
					time_type t_left = t_cur, t_right = t_cur+dt, tNew;
					double y_bar = x[N*m+1], y_left = x[N*m], y_right = x0[N*m], yNew;
					// t- and y-tolerances and rounding error limits
					double tTol = rtTol * fabs(t_right-t_left);
					double yTol = ryTol * fabs(y_right-y_left);
					double tEps = rtEps * fabs(t_right-t_left);
					double yEps = ryEps * fabs(y_right-y_left);
                    
					// Variables
					double dtHiLo = t_right - t_left;  // time difference of lower and upper bound
					double dtUpDn = 0;                 // time bound change
					double dyHiLo = y_right - y_left;  // y difference of lower and upper bound
					double dyHiTarg = 0;
					double dyLoTarg = 0;

					// Get initial 3rd grid point for quadratic interpolation by
					// linear interpolation of t between initial t_left and t_right
					time_type t3rd = 0;
					double y3rd = 0;
					if (jIpol == "quadratic") {
						t3rd = tIpol(y_bar, t_left, t_right, t3rd, y_left, y_right, y3rd, "linear");
						time_type temp_dt = t_cur - t3rd;
						state_type x2 = x;
						fTarg<System> (system, x2, t_cur, temp_dt);
						y3rd = x2[N*m];
					}

					int jConverg = 0;    // 0: iteration not yet converged
					int cnt = 0;         // Counter
					int jTypeIt = 1;     // Iteration step type
					int jUpDn = 1;       // 1: t_left moved to tNew

					// Loop for finding the t location
					while (jConverg == 0 && cnt < nItMax) {
						cnt++;           // Increase counter

						// Type 1 iteration step
						if (jTypeIt == 1) {
							tNew = tIpol(y_bar, t_left, t_right, t3rd, y_left, y_right, y3rd, jIpol);
						}
						// Type 2 iteration step
						if (jTypeIt == 2) {
							tNew = 0.5 * (t_left + t_right);
							if (dtUpDn > 0) tNew = t_left+(1+mUpDn)*dtUpDn;
							else if (dtUpDn < 0) tNew = t_right+(1+mUpDn)*dtUpDn;
						}
						// Type 3 iteration step
						if (jTypeIt == 3) {
							tNew = 0.5 * (t_left + t_right);
						}

						// Get yNew
						dt = tNew - t_cur;
						x0 = x;
						fTarg<System> (system, x0, t_cur, dt);
						yNew = x0[N*m];

						// tNew is the new lower bound or "exact" solution
						if ((yNew-y_bar)*(y_left-y_bar) >= 0) {
							t3rd = t_left;
							y3rd = y_left;
							dtUpDn = tNew - t_left;
							
							// Set jUpDn=1 or 2 if t_left is moved up to tNew
							jUpDn = 1;
							if ((yNew - y_bar)*(y_left-y_bar) == 0) {
								jUpDn = 2;
							}
							t_left = tNew;
							y_left = yNew;

							dtHiLo = t_right - t_left;
							dyHiLo = y_right - y_left;

						} else if ((yNew-y_bar)*(y_right-y_bar) > 0) {
							// tNew is the new upper bound
							t3rd = t_right;
							y3rd = y_right;
							dtUpDn = tNew - t_right;
							
							t_right = tNew;
							y_right = yNew;

							dtHiLo = t_right - t_left;
							dyHiLo = y_right - y_left;
						}

						// Switch to next iteration step type by cycling thru 
						// 1->2->3->..., except if jIpol=="bisect".
						// Switch 1->3 if 1->2 is not allowed.
						if (jIpol == "bisect") {
							jTypeIt = 1;
						} else {
							if ((jTypeIt==2) || (jTypeIt==3) || ((jTypeIt==1) && (dtHiLo>((2+fabs(mUpDn))*fabs(dtUpDn))))) {
								jTypeIt = 1 + jTypeIt % 3;
							} else if ((jTypeIt == 1) && (dtHiLo <= ((2+fabs(mUpDn))*fabs(dtUpDn)))) {
								jTypeIt = 3;
							}

							if (fabs(dtUpDn) < tEps) {
								jTypeIt = 3;
							}
						}

						// Check convergence of the iteration
						dyLoTarg = y_left - y_bar;
						dyHiTarg = y_right - y_bar;
						if ((fabs(dyHiTarg) < yEps) || (fabs(dyLoTarg) < yEps)) {
							jConverg = 1;
						}
						if (fabs(dyHiLo) < yEps) {
							jConverg = 2;
						}
						if (dtHiLo < tEps) {
							jConverg = 3;
						}
						if (dtHiLo < tTol && fabs(dyHiLo) < yTol) {
							jConverg = 4;
						}
						if (jUpDn == 2) {
							jConverg = 5;
						}
						
						// Get iteration results
						if (jConverg==1) {
							if (fabs(dyLoTarg) < fabs(dyHiTarg)) {
								t_cur = t_left;
							} else if (fabs(dyLoTarg) > fabs(dyHiTarg)) {
								t_cur = t_right;
							} else if (fabs(dyLoTarg) == fabs(dyHiTarg)) {
								t_cur = 0.5 * (t_left + t_right);
							}
						}
						if (jConverg==2 || jConverg==3 || jConverg==4) {
							t_cur = 0.5 * (t_left+t_right);
						}
						if (jConverg==5) {
							t_cur = t_left;
						}
					}
					x = x0;                     // Renew x states
					x[N*m] = y_bar+0.1;             // Set y_bar

}


} // namespace hybrid_model

#endif
