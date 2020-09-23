// Version 4.3

#pragma once
#include <math.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_trig.h>


const double R = 8.314510;           //In J / mol K from Woan
const double kB = 1.380658E-23;      //In J / K from Woan
const double hbar = 1.05457266E-34;  //Also from Woan
const double N_A = 6.0221367E23;    //Also from Woan

struct A_params { double T; double k_D; double G_inf; double density; double viscosity; };
struct C_params { double T;  double k_D; double sound_speed; };

double omega_T_func(double k, double G_inf, double density, double visc ) {
	// Is mass density, not number density
	// visc currently redundant but implemented for future use
	return k * sqrt(G_inf / density);
}

double omega_L_func(double sound_speed, double k) {
	return sound_speed * k;
}

double k_prop_func(double density, double G_inf, double viscosity) {
	// This is mass density, not number density
	return sqrt(G_inf*density)/(viscosity*M_SQRT2);
}

double g_k(double k_D, double k) {
	return 3.0 * N_A * pow(k, 2.0) * pow(k_D, -3.0);
}

double full_integrand_A(double k, void* p) {
	//For U_T
	struct A_params* A_p = (struct A_params*)p;
	double T = (A_p->T);
	double k_D = (A_p->k_D);
	double G_inf = (A_p->G_inf);
	double density = (A_p->density);
	double viscosity = (A_p->viscosity);
	return hbar * g_k(k_D, k) * omega_T_func(k, G_inf, density, viscosity) / (exp(hbar * omega_T_func(k, G_inf, density, viscosity) / (kB * T)) - 1.0);
}

double full_integrand_C(double k, void* p) {
	//For U_L
	//Named C for consistency with version 1
	struct C_params* C_p = (struct C_params*)p;
	double T = (C_p->T);
	double k_D = (C_p->k_D);				
	double sound_speed = (C_p->sound_speed);

	return hbar * omega_L_func(sound_speed, k) * g_k(k_D, k) / (exp(hbar* omega_L_func(sound_speed, k) / (kB*T)) - 1.0);
}

double U_T(double anharm_param, double beta, double T, double k_lower, double k_D, double G_inf, double density, double viscosity) {
	gsl_integration_workspace* w
		= gsl_integration_workspace_alloc(10000);

	double error = 0.0;
	double int_result = 0.0;

	gsl_function F;
	struct A_params A_p = { T, k_D, G_inf, density, viscosity };
	F.function = &full_integrand_A;
	F.params = &A_p;
	gsl_integration_qags(&F, k_lower, k_D, 0, 1E-5, 10000, w, &int_result, &error);
	gsl_integration_workspace_free(w);
	return 2.0 * (1.0 + (anharm_param * beta * T)) * int_result;
}

double U_L(double anharm_param, double beta, double T, double k_D, double sound_speed) {
	gsl_integration_workspace* w
		= gsl_integration_workspace_alloc(10000);

	double error = 0.0;
	double int_result = 0.0;

	gsl_function F;
	struct C_params C_p = { T, k_D, sound_speed };
	F.function = &full_integrand_C;
	F.params = &C_p;

	gsl_integration_qags(&F, 0, k_D, 0, 1E-5, 10000, w, &int_result, &error);
	gsl_integration_workspace_free(w);
	return (1.0 + (anharm_param * beta * T)) * int_result;
}