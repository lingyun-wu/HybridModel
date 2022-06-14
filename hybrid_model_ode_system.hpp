#ifndef HYBRID_MODEL_ODE_SYSTEM_HPP
#define HYBRID_MODEL_ODE_SYSTEM_HPP

#include<iostream>
#include<cmath>
#include<string>
#include<boost/numeric/odeint.hpp>
#include<boost/random.hpp>
#include<cstdlib>
#include<ctime>

namespace hybrid_model {

enum CellSpecies {
	wc_1_r0 = 0,
	wc_1_r1 = 1,
	wc_1 = 2,
	frq_r = 3,
	frq = 4,
	wcc = 5,
	ccg_r1 = 6,
	ccg = 7,
	s_i = 8,
	g_f = 9,
	g_c = 10,
	life_state = 11
}

enum CellGroupSpecies {
	integral = 0,
	integral_bar = 1,
	alive_cell_num = 2,
	inactive_cell_num = 3
}

enum CellParam {
	A_f = 0,
  B_f = 1,
  S_1 = 2,
  S_3 = 3,
  S_4 = 4,
  D_1 = 5,
  D_3 = 6,
  C_1 = 7,
  L_1 = 8,
  L_3 = 9,
  D_4 = 10,
  D_6 = 11,
  D_7 = 12,
  D_8 = 13,
  C_2 = 14,
  P = 15,
  A_c = 16,
  B_c = 17,
  S_c = 18,
  L_c = 19,
  D_cr = 20,
  D_cp = 21,
  C_3 = 22,
  eta = 23,
  D_9 = 24,
  C_4 = 25,
  K = 26,
  death_active_coeff = 27,
  birth_coeff = 28,
  deactive_coeff = 29,
  active_coeff = 30,
  death_inactive_coeff = 31,
	Q = 32,
	light_interval = 33,
	dark_interval = 34,
	mf = 35,
	c_eff_exp = 36
}

struct SystemDeterministicPart {
	typedef vector<double> state_type;
	
	vector<double>& paramter_;
	vector<double>& Se_;

	// Constructor
	SystemDeterministicPart(vector<double>& parameter, vector<double>& Se): parameter_(parameter), Se_(Se) {}

	// Overloads () operator
	void operator() (const state_type& x, state_type& dxdt, double t, double dt) const {
		
	}
}
	




} // namespace hybrid_model

#endif
