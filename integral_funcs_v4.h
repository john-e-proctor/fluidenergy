// Version 4.3
const int inputdatarows = 100;
const int anharm_max_vals = 2;     // max numbers of values we are allowed to try for each of these variables
const int k_D_max_vals = 100;   // to define size of the array
const int G_inf_max_vals = 100;

#ifndef INTEGRALFUNCS_V4_H
#define INTEGRALFUNCS_V4
double omega_T_func(double k, double G_inf, double density, double visc);
double omega_L_func(double sound_speed, double k);
double k_prop_func(double density, double G_inf, double viscosity);
double full_integrand_A(double k, void* p);
double full_integrand_C(double k, void* p);
double U_T(double gamma, double beta, double T, double k_lower, double k_D, double G_inf, double density, double viscosity);
double U_L(double gamma, double beta, double T, double k_D, double sound_speed);
int parameter_write(double fit_parameters[3][3],
    double best_gamma, double best_k_D, double best_G_inf, double fluid_molar_mass, double visc_coeff, double reference_num,
    string P_filename, string input_filename);
int U_write(double U_values[inputdatarows][2], string U_filename);
int k_write(double k_values[inputdatarows][4], string k_filename);
int fit_quality_write(int p_best, int q_best, int r_best, double parameter_vals[4][anharm_max_vals][k_D_max_vals][G_inf_max_vals],
    int num_k_D_vals, int num_G_inf_vals, string k_D_fit_q_filename, string G_inf_fit_q_filename);

#endif // !INTEGRALFUNCS_V4