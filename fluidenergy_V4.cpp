// Version 4.3
const double version_num = 4.3;

using namespace std;

#include <iostream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <fstream>
#include <integral_funcs_v4.h>


const double k_gap_initial_ref = 0.1;  // Initial refinement proportion for k_gap
const double R = 8.314510;          //In J / mol K from Woan
const double kB = 1.380658E-23;     //In J / K from Woan
const double hbar = 1.05457266E-34; //Also from Woan
const double N_A = 6.0221367E23;    //Also from Woan
const double molar_mass = 0.02018;  // for Ne, in kg
const double ref_num = 1.0;      // To be written to .par file for reference

double experimentaldata[inputdatarows][9]{};    //T, P, V, U_exp, sound speed, visc, beta, U_th (buffered values), density in kg/m3;
double lower_isochore[inputdatarows][6]{};        //T, P, V, U, sound speed, visc
double upper_isochore[inputdatarows][6]{};        //T, P, V, U, sound speed, visc
double parameter_values[4][anharm_max_vals][k_D_max_vals][G_inf_max_vals]{};
//                      ^ in order anharm. parameter, k_D, G_inf, Error in U
double bestfit_U[inputdatarows][2]{}; //T, U
double bestfit_U_L[inputdatarows]{};
double refined_values[inputdatarows][4]{};  // T, k_prop (init), k_prop (fitted), k_prop refinement increment
double refined_U[inputdatarows][2]{};

int main()
{
    double input_buffer = 0.0;        //input buffer
    double k_prop = 0.0;           //For buffering k_gap values
    double exp_U_mean = 0.0;        // To hold mean value of experimental U from NIST
    double th_U_mean = 0.0;         // To hold mean theoretical value of U

    int p_best = 0;
    int q_best = 0;
    int r_best = 0;

    double parameter_space[3][3]{}; //Parameter space limits loaded from file

    //Load NIST data
#pragma region Region A Load NIST data, calculate beta, mean value of U, densities;
    ifstream in_more_data;
    string inputstem = "hello";
    string inputfilename = "inputdata.txt";
    cout << "Enter input data file name STEM only:\n" << endl;
    cin >> inputstem;
    cout << "\n";
    inputfilename = inputstem + "_input.txt";
    in_more_data.open(inputfilename);
    if (!in_more_data) {
        cerr << "Input data could not be loaded" << endl;
        exit(1); //Exits the program
    }
    int IP_s = 0;
    while (!in_more_data.eof()) {
        in_more_data >> experimentaldata[IP_s][0] >> experimentaldata[IP_s][1] >> experimentaldata[IP_s][2] >> experimentaldata[IP_s][3] >> experimentaldata[IP_s][4] >> experimentaldata[IP_s][5];
        IP_s++;
    }
    in_more_data.close();

    if (IP_s == inputdatarows) { //this IS correct
        cout << "Correct number of input rows" << endl;
        cout << "\n";
    }
    else {
        cout << "***** WRONG number of input rows *****" << endl;
        cout << "\n";
        exit(1);

    }

    // Load the upper and lower isochores for calculating thermal expansion coefficient
    // Upper first
    string upperisocfilename = "hello";
    upperisocfilename = inputstem + "_upperisochore.txt";
    in_more_data.open(upperisocfilename);
    if (!in_more_data) {
        cerr << "Upper isochore could not be loaded" << endl;
        exit(1); //Exits the program
    }
    IP_s = 0;
    while (!in_more_data.eof()) {
        in_more_data >> upper_isochore[IP_s][0] >> upper_isochore[IP_s][1] >> upper_isochore[IP_s][2] >> upper_isochore[IP_s][3] >> upper_isochore[IP_s][4] >> upper_isochore[IP_s][5];
        IP_s++;
    }
    in_more_data.close();

    if (IP_s == inputdatarows) { //this IS correct
        cout << "Correct number of input rows in upper isochore" << endl;
        cout << "\n";
    }
    else {
        cout << "***** WRONG number of input rows in upper isochore *****" << endl;
        cout << "\n";
        exit(1);
    }

    // Then lower
    string lowerisocfilename = "hello";
    lowerisocfilename = inputstem + "_lowerisochore.txt";

    in_more_data.open(lowerisocfilename);
    if (!in_more_data) {
        cerr << "Lower isochore could not be loaded" << endl;
        exit(1); //Exits the program
    }
    IP_s = 0;
    while (!in_more_data.eof()) {
        in_more_data >> lower_isochore[IP_s][0] >> lower_isochore[IP_s][1] >> lower_isochore[IP_s][2] >> lower_isochore[IP_s][3] >> lower_isochore[IP_s][4] >> lower_isochore[IP_s][5];
        IP_s++;
    }
    in_more_data.close();

    if (IP_s == inputdatarows) { //this IS correct
        cout << "Correct number of input rows in lower isochore" << endl;
        cout << "\n";
    }
    else {
        cout << "***** WRONG number of input rows in lower isochore *****" << endl;
        cout << "\n";
        exit(1);
    }

    //Subtract gas-like viscosity component if desired
    double gas_visc_coeff = 0.0;
    cout << "Enter parameter for gas-like viscosity to subtract\n" << endl;
    cout << "In microPascal.seconds, 0 if no subtraction required\n" << endl;
    cin >> gas_visc_coeff;
    cout << "\n";
    for (int s = 0; s < inputdatarows; s++) {
        experimentaldata[s][5] = experimentaldata[s][5] - (gas_visc_coeff * pow(experimentaldata[s][0], 0.5));
    }


    //Calculate beta at each point, calculate densities and tot up for mean U, convert viscosities
    exp_U_mean = 0.0;
    for (int s = 1; s < (inputdatarows - 1); s++) {
        //Thermal expansion first
        experimentaldata[s][6] = -(1.0 / experimentaldata[s][2]) * ((experimentaldata[s + 1][1] - experimentaldata[s - 1][1]) / (experimentaldata[s + 1][0] - experimentaldata[s - 1][0])) * ((upper_isochore[s][2] - lower_isochore[s][2]) / (upper_isochore[s][1] - lower_isochore[s][1]));
        exp_U_mean = exp_U_mean + experimentaldata[s][3];                   // Adding to mean value of U
        experimentaldata[s][8] = molar_mass * 1000 / experimentaldata[s][2]; //Density in kg / m3
        experimentaldata[s][5] = experimentaldata[s][5] * 1E-6; //to convert viscosity to Pa.s
    }
    exp_U_mean = exp_U_mean / (double(inputdatarows) - 2.0);

#pragma endregion

#pragma region Region B Get parameter space limits from file, check with user
    ifstream indata;
    indata.open("parameter_space.txt");
    if (!indata) {
        cerr << "Parameter space could not be loaded" << endl;
        exit(1); //Exits the program
    }
    int param_r = 0;
    while (!indata.eof()) {
        indata >> parameter_space[param_r][0] >> parameter_space[param_r][1] >> parameter_space[param_r][2];
        param_r++;
    }
    indata.close();

    cout << "Lambda min = " << parameter_space[0][0] << endl;
    cout << "Enter new value, or 0 to accept existing value:" << endl;
    cin >> input_buffer;
    if (input_buffer != 0.0) parameter_space[0][0] = input_buffer;
    cout << "Lambda min = " << parameter_space[0][0] << endl;
    cout << "\n";

    cout << "Lambda max = " << parameter_space[1][0] << endl;
    cout << "Enter new value, or 0 to accept existing value:" << endl;
    cin >> input_buffer;
    if (input_buffer != 0.0) parameter_space[1][0] = input_buffer;
    cout << "Lambda max = " << parameter_space[1][0] << endl;
    cout << "\n";

    cout << "Number of Lambda values to search = " << parameter_space[2][0] << endl;
    cout << "Enter new value, or 0 to accept existing value:" << endl;
    cin >> input_buffer;
    if (input_buffer != 0.0) parameter_space[2][0] = input_buffer;
    if (parameter_space[2][0] > ((double)anharm_max_vals)) parameter_space[2][0] = ((double)anharm_max_vals);
    cout << "Number of Lambda to search = " << parameter_space[2][0] << endl;
    cout << "\n";

    cout << "k_D_min = " << parameter_space[0][1] << endl;
    cout << "Enter new value, or 0 to accept existing value:" << endl;
    cin >> input_buffer;
    if (input_buffer != 0.0) parameter_space[0][1] = input_buffer;
    cout << "k_D_min = " << parameter_space[0][1] << endl;
    cout << "\n";

    cout << "k_D_max = " << parameter_space[1][1] << endl;
    cout << "Enter new value, or 0 to accept existing value:" << endl;
    cin >> input_buffer;
    if (input_buffer != 0.0) parameter_space[1][1] = input_buffer;
    cout << "k_D_max = " << parameter_space[1][1] << endl;
    cout << "\n";

    cout << "Number of k_D values to search = " << parameter_space[2][1] << endl;
    cout << "Enter new value, or 0 to accept existing value:" << endl;
    cin >> input_buffer;
    if (input_buffer != 0.0) parameter_space[2][1] = input_buffer;
    if (parameter_space[2][1] > ((double)k_D_max_vals)) parameter_space[2][1] = ((double)k_D_max_vals);
    cout << "Number of k_D values to search = " << parameter_space[2][1] << endl;
    cout << "\n";

    cout << "G_inf_min = " << parameter_space[0][2] << endl;
    cout << "Enter new value, or 0 to accept existing value:" << endl;
    cin >> input_buffer;
    if (input_buffer != 0.0) parameter_space[0][2] = input_buffer;
    cout << "G_inf_min = " << parameter_space[0][2] << endl;
    cout << "\n";

    cout << "G_inf_max = " << parameter_space[1][2] << endl;
    cout << "Enter new value, or 0 to accept existing value:" << endl;
    cin >> input_buffer;
    if (input_buffer != 0.0) parameter_space[1][2] = input_buffer;
    cout << "G_inf_max = " << parameter_space[1][2] << endl;
    cout << "\n";

    cout << "Number of G_inf values to search = " << parameter_space[2][2] << endl;
    cout << "Enter new value, or 0 to accept existing value:" << endl;
    cin >> input_buffer;
    if (input_buffer != 0.0) parameter_space[2][2] = input_buffer;
    if (parameter_space[2][2] > ((double)G_inf_max_vals)) parameter_space[2][2] = ((double)G_inf_max_vals);
    cout << "Number of G_inf values to search = " << parameter_space[2][2] << endl;
    cout << "\n";

#pragma endregion

#pragma region Region C Populate array of actual parameter space values
    for (int p = 0; p < ((int)parameter_space[2][0]); p++) {           //cycle through lambda values
        for (int q = 0; q < ((int)parameter_space[2][1]); q++) {       //cycle through k_D values
            for (int r = 0; r < ((int)parameter_space[2][2]); r++) {   //cycle through G_inf values
                //lambda
                if (parameter_space[2][0] < 1.5) {  // This is to check if it is equal to 1, but checking == with doubles potentially not robust
                    parameter_values[0][p][q][r] = (parameter_space[0][0] + parameter_space[1][0]) / 2.0;
                }
                else {
                    parameter_values[0][p][q][r]
                        = parameter_space[0][0] + ((double)p) * (parameter_space[1][0] - parameter_space[0][0]) / (parameter_space[2][0] - 1.0);

                }

                //k_D
                if (parameter_space[2][1] < 1.5) {  // This is to check if it is equal to 1, but checking == with doubles potentially not robust
                    parameter_values[1][p][q][r] = (parameter_space[0][1] + parameter_space[1][1]) / 2.0;
                }
                else {
                    parameter_values[1][p][q][r]
                        = parameter_space[0][1] + ((double)q) * (parameter_space[1][1] - parameter_space[0][1]) / (parameter_space[2][1] - 1.0);
                }

                //G_inf
                if (parameter_space[2][2] < 1.5) {  // This is to check if it is equal to 1, but checking == with doubles potentially not robust
                    parameter_values[2][p][q][r] = (parameter_space[0][2] + parameter_space[1][2]) / 2.0;
                }
                else {
                    parameter_values[2][p][q][r]
                        = parameter_space[0][2] + ((double)r) * (parameter_space[1][2] - parameter_space[0][2]) / (parameter_space[2][2] - 1.0);
                }
            }
        }
    }
#pragma endregion 
    
#pragma region Region D Calculate values of U_th, and least sq error for each point in parameter space
    double U_adjust = 0.0;
    for (int p = 0; p < (int(parameter_space[2][0])); p++) {           //lambda
        for (int q = 0; q < (int(parameter_space[2][1])); q++) {       //Omega_D
            for (int r = 0; r < (int(parameter_space[2][2])); r++) {   //G_inf
                th_U_mean = 0.0;
                for (int s = 1; s < (inputdatarows - 1); s++) { //Cycle through input data rows
                    // Viscosity already converted to Pa.s
                    k_prop = k_prop_func(experimentaldata[s][8], parameter_values[2][p][q][r], experimentaldata[s][5]);
                    experimentaldata[s][7] = U_L(parameter_values[0][p][q][r],
                        experimentaldata[s][6],
                        experimentaldata[s][0],
                        parameter_values[1][p][q][r],
                        experimentaldata[s][4]
                    ) + R*experimentaldata[s][0]; // To add the kinetic energy for transverse dof
                    if (k_prop < parameter_values[1][p][q][r]) {
                        experimentaldata[s][7] = experimentaldata[s][7] +
                            0.5*U_T(parameter_values[0][p][q][r],
                                experimentaldata[s][6],
                                experimentaldata[s][0],
                                k_prop,
                                parameter_values[1][p][q][r],
                                parameter_values[2][p][q][r],
                                experimentaldata[s][8],
                                experimentaldata[s][5]
                            );
                    }
                    experimentaldata[s][7] = experimentaldata[s][7] * 1.0E-3; //Convert to kJ/mol to compare with NIST data
                    th_U_mean = th_U_mean + experimentaldata[s][7];
                }
                th_U_mean = th_U_mean / (double(inputdatarows) - 2.0);
                //Adjust by U_0 and calculate least sq. error
                U_adjust = exp_U_mean - th_U_mean; //U_adjust needs to be kept for later also!
                parameter_values[3][p][q][r] = 0;
                for (int s = 1; s < (inputdatarows - 1); s++) {
                    experimentaldata[s][7] = experimentaldata[s][7] + U_adjust;
                    parameter_values[3][p][q][r] = parameter_values[3][p][q][r] + pow((experimentaldata[s][3] - experimentaldata[s][7]), 2.0);
                }
                parameter_values[3][p][q][r] = parameter_values[3][p][q][r] / (double(inputdatarows) - 2.0);

            }
        }
    }
#pragma endregion

#pragma region Region E Recalculate U with best fit parameters and  subtract for fit
    double fit_error = DBL_MAX;
    for (int p = 0; p < (int(parameter_space[2][0])); p++) {           //lambda
        for (int q = 0; q < (int(parameter_space[2][1])); q++) {       //Omega_D
            for (int r = 0; r < (int(parameter_space[2][2])); r++) {   //G_inf
                if (parameter_values[3][p][q][r] < fit_error) {
                    p_best = p;
                    q_best = q;
                    r_best = r;
                    fit_error = parameter_values[3][p][q][r];
                }
            }
        }
    }


    th_U_mean = 0;
    for (int s = 1; s < (inputdatarows - 1); s++) { //Cycle through input data rows
        bestfit_U[s][0] = experimentaldata[s][0]; //Copying T across
        refined_U[s][0] = experimentaldata[s][0]; // And also for the array to store U following refinement of k_prop
        refined_values[s][0] = experimentaldata[s][0]; //And to array for saving k_prop values
        refined_values[s][3] = k_gap_initial_ref;
        k_prop = k_prop_func(experimentaldata[s][8], parameter_values[2][p_best][q_best][r_best], experimentaldata[s][5]);
        refined_values[s][1] = k_prop;
        refined_values[s][2] = k_prop;
        bestfit_U_L[s] = U_L(parameter_values[0][p_best][q_best][r_best],
            experimentaldata[s][6],
            experimentaldata[s][0],
            parameter_values[1][p_best][q_best][r_best],
            experimentaldata[s][4]
        );
        bestfit_U[s][1] = bestfit_U_L[s] + R * experimentaldata[s][0];

        if (k_prop < parameter_values[1][p_best][q_best][r_best]) {
            bestfit_U[s][1] = bestfit_U[s][1] +
                0.5*U_T(parameter_values[0][p_best][q_best][r_best],
                    experimentaldata[s][6],
                    experimentaldata[s][0],
                    k_prop,
                    parameter_values[1][p_best][q_best][r_best],
                    parameter_values[2][p_best][q_best][r_best],
                    experimentaldata[s][8],
                    experimentaldata[s][5]
                );
        }

        bestfit_U[s][1] = bestfit_U[s][1] * 1.0E-3;
        th_U_mean = th_U_mean + bestfit_U[s][1];
    }
    th_U_mean = th_U_mean / (double(inputdatarows) - 2.0);
    U_adjust = exp_U_mean - th_U_mean;
    for (int s = 1; s < (inputdatarows - 1); s++) {
        bestfit_U[s][1] = bestfit_U[s][1] + U_adjust;
    }
#pragma endregion

#pragma region Region F direct refinement of k_prop
    double k_prop_plus = 0;
    double k_prop_minus = 0;
    double U_plus = 0;
    double U_minus = 0;
    double U_adjust_JperMol = 0;
    U_adjust_JperMol = U_adjust*1E3;

    for (int v = 0; v < 50; v++) { // 50 is number of rounds of refinement (for now...)
        for (int s = 1; s < (inputdatarows - 1); s++) {
            k_prop_plus = refined_values[s][2] * (1 + refined_values[s][3]);
            k_prop_minus = refined_values[s][2] * (1 - refined_values[s][3]);

            refined_U[s][1] = bestfit_U_L[s] + R * experimentaldata[s][0] + U_adjust_JperMol; // Fill in longitudinal contribution, RT term and zero adjust
            U_plus = refined_U[s][1];
            U_minus = refined_U[s][1];

            if (refined_values[s][2] < parameter_values[1][p_best][q_best][r_best]) {
                refined_U[s][1] = refined_U[s][1] +
                    0.5*U_T(parameter_values[0][p_best][q_best][r_best],
                        experimentaldata[s][6],
                        experimentaldata[s][0],
                        refined_values[s][2], //k prop
                        parameter_values[1][p_best][q_best][r_best],
                        parameter_values[2][p_best][q_best][r_best],
                        experimentaldata[s][8],
                        experimentaldata[s][5]
                    );
            }

            //k prop plus
            if (k_prop_plus < parameter_values[1][p_best][q_best][r_best]) {
                U_plus = U_plus +
                    0.5*U_T(parameter_values[0][p_best][q_best][r_best],
                        experimentaldata[s][6],
                        experimentaldata[s][0],
                        k_prop_plus, //k prop
                        parameter_values[1][p_best][q_best][r_best],
                        parameter_values[2][p_best][q_best][r_best],
                        experimentaldata[s][8],
                        experimentaldata[s][5]
                    );
            }

            //k prop minus
            if (k_prop_minus < parameter_values[1][p_best][q_best][r_best]) {
                U_minus = U_minus +
                    0.5*U_T(parameter_values[0][p_best][q_best][r_best],
                        experimentaldata[s][6],
                        experimentaldata[s][0],
                        k_prop_minus, //k prop
                        parameter_values[1][p_best][q_best][r_best],
                        parameter_values[2][p_best][q_best][r_best],
                        experimentaldata[s][8],
                        experimentaldata[s][5]
                    );
            }

            double U_exp_JperMol = 0;
            U_exp_JperMol = experimentaldata[s][3]*1E3;
            // Is existing k prop best fit?
            if (pow((U_exp_JperMol - refined_U[s][1]), 2.0) < pow((U_exp_JperMol - U_plus), 2.0)) {
                if (pow((U_exp_JperMol - refined_U[s][1]), 2.0) < pow((U_exp_JperMol - U_minus), 2.0)) {
                    refined_values[s][3] = refined_values[s][3] * 0.3; // Decrease increment
                }
            }

            // Is k prop plus best fit?
            if (pow((U_exp_JperMol - U_plus), 2.0) < pow((U_exp_JperMol - refined_U[s][1]), 2.0)) {
                if (pow((U_exp_JperMol - U_plus), 2.0) < pow((U_exp_JperMol - U_minus), 2.0)) {
                    refined_values[s][2] = k_prop_plus;
                }
            }

            // Is k prop minus best fit?
            if (pow((U_exp_JperMol - U_minus), 2.0) < pow((U_exp_JperMol - U_plus), 2.0)) {
                if (pow((U_exp_JperMol - U_minus), 2.0) < pow((U_exp_JperMol - refined_U[s][1]), 2.0)) {
                    refined_values[s][2] = k_prop_minus;
                }
            }

        }
    }

#pragma endregion

#pragma region Region G Recalculate U with final k-prop values and save to file
    for (int s = 1; s < (inputdatarows - 1); s++) {
        refined_U[s][1] = bestfit_U_L[s] + R * experimentaldata[s][0] + U_adjust_JperMol;

        if (refined_values[s][1] < parameter_values[1][p_best][q_best][r_best]) {
            refined_U[s][1] = refined_U[s][1] +
                0.5*U_T(parameter_values[0][p_best][q_best][r_best],
                    experimentaldata[s][6],
                    experimentaldata[s][0],
                    refined_values[s][2], //k prop
                    parameter_values[1][p_best][q_best][r_best],
                    parameter_values[2][p_best][q_best][r_best],
                    experimentaldata[s][8],
                    experimentaldata[s][5]
                );
        }

        refined_U[s][1] = refined_U[s][1] * 1E-3;
    }


    //Print parameters to screen and save best fit U to file

    string filename = "hello.txt";
    string U_filename = "hello.txt";
    string U_ref_filename = "hello.txt";
    string P_filename = "hello.txt;";
    string k_filename = "hello.txt";
    string k_D_fit_q_filename = "hello.txt";
    string G_inf_fit_q_filename = "hello.txt";

    cout << "Please enter a file name to write (NOT including extension): ";
    cin >> filename;
    U_filename = filename + "_search" + ".dat";
    U_ref_filename = filename + "_refined" + ".dat";
    P_filename = filename + ".par";
    k_filename = filename + "_kgaps" + ".dat";
    k_D_fit_q_filename = filename + "k_D_fitq" + ".dat";
    G_inf_fit_q_filename = filename + "G_inf_fitq" + ".dat";

    //U  values after initial parameter search and refinement
    int init_U_write = 0;
    init_U_write = U_write(bestfit_U, U_filename);
    init_U_write = U_write(refined_U, U_ref_filename);

    // Parameters
    int param_write = 0;
    param_write = parameter_write(parameter_space,
        parameter_values[0][p_best][q_best][r_best], parameter_values[1][p_best][q_best][r_best], parameter_values[2][p_best][q_best][r_best],
        molar_mass, gas_visc_coeff, ref_num, P_filename, inputfilename);

    // k prop values
    int kgapwrite = 0;
    kgapwrite = k_write(refined_values, k_filename);

    // fit quality
    int fit_write = 0;
    fit_write = fit_quality_write(p_best, q_best, r_best, parameter_values,
        (int)parameter_space[2][1], (int)parameter_space[2][2], k_D_fit_q_filename, G_inf_fit_q_filename);

#pragma endregion

    cout << "Finished!\n";
    return 0;
}

//Write U values
int U_write(double U_values[inputdatarows][2], string U_filename) {
    FILE* excitingresults;
    errno_t err;
    err = fopen_s(&excitingresults, U_filename.c_str(), "w");
    if (err || !excitingresults) {
        cout << "Unable to save U values\n";
    }
    else {
        for (int s = 1; s < (inputdatarows - 1); s++) {
            fprintf(excitingresults, "%E\t", U_values[s][0]);
            fprintf(excitingresults, "%E\n", U_values[s][1]);
        }
    }
    fclose(excitingresults);
    return 0;
}

//Write k values
int k_write(double k_values[inputdatarows][4], string k_filename) {
    FILE* excitingresults;
    errno_t err;
    err = fopen_s(&excitingresults, k_filename.c_str(), "w");
    if (err || !excitingresults) {
        cout << "Unable to save k_prop values\n";
    }
    else {
        for (int s = 1; s < (inputdatarows - 1); s++) {
            fprintf(excitingresults, "%E\t", k_values[s][0]);
            fprintf(excitingresults, "%E\t", k_values[s][1]);
            fprintf(excitingresults, "%E\t", k_values[s][2]);
            fprintf(excitingresults, "%E\n", k_values[s][3]);
        }
    }
    fclose(excitingresults);
    return 0;
}

// Write parameters to file
int parameter_write(double fit_parameters[3][3],
    double best_anharm_param, double best_k_D, double best_G_inf, double fluid_molar_mass, double visc_coeff, double reference_num,
    string P_filename, string input_filename) {
    // Parameters
    errno_t err;
    FILE* excitingparameters;
    err = fopen_s(&excitingparameters, P_filename.c_str(), "w");
    if (err || !excitingparameters) {
        cout << "Unable to save parameters\n";
    }
    else {
        string input_fileinfo = "hello";
        input_fileinfo = "Input file: " + input_filename + "- IS THIS CORRECT?";
        cout << input_fileinfo << endl;
        fprintf(excitingparameters, "Anharm. parameter\t");
        fprintf(excitingparameters, "k_D\t");
        fprintf(excitingparameters, "G_inf\t");
        fprintf(excitingparameters, "Other\n");

        fprintf(excitingparameters, "%E\t", best_anharm_param);
        fprintf(excitingparameters, "%E\t", best_k_D);
        fprintf(excitingparameters, "%E\t", best_G_inf);
        fprintf(excitingparameters, "%E\n", fluid_molar_mass);

        fprintf(excitingparameters, "%E\t", fit_parameters[0][0]);
        fprintf(excitingparameters, "%E\t", fit_parameters[0][1]);
        fprintf(excitingparameters, "%E\t", fit_parameters[0][2]);
        fprintf(excitingparameters, "%E\n", reference_num);

        fprintf(excitingparameters, "%E\t", fit_parameters[1][0]);
        fprintf(excitingparameters, "%E\t", fit_parameters[1][1]);
        fprintf(excitingparameters, "%E\t", fit_parameters[1][2]);
        fprintf(excitingparameters, "%E\n", visc_coeff);

        fprintf(excitingparameters, "%E\t", fit_parameters[2][0]);
        fprintf(excitingparameters, "%E\t", fit_parameters[2][1]);
        fprintf(excitingparameters, "%E\t", fit_parameters[2][2]);
        fprintf(excitingparameters, "%E\n", version_num);

    }
    fclose(excitingparameters);

    cout << "Best fit anharmonicity parameter = " << best_anharm_param << endl;
    cout << "Best fit k_D = " << best_k_D << endl;
    cout << "Best fit G_inf = " << best_G_inf << endl;
    cout << "\n";
    return 0;
}

// Write trends in fit quality to file
int fit_quality_write(int p_best, int q_best, int r_best, double parameter_vals[4][anharm_max_vals][k_D_max_vals][G_inf_max_vals],
    int num_k_D_vals, int num_G_inf_vals, string k_D_fit_q_filename, string G_inf_fit_q_filename) {
    FILE* excitingresults;
    errno_t err;
    
    // first write k_D values at best fit G_inf
    err = fopen_s(&excitingresults, k_D_fit_q_filename.c_str(), "w");
    if (err || !excitingresults) {
        cout << "Unable to save fit quality values\n";
    }
    else {
        for (int s = 1; s < (num_k_D_vals - 1); s++) {
            fprintf(excitingresults, "%E\t", parameter_vals[1][p_best][s][r_best]);   // k_D
            fprintf(excitingresults, "%E\n", parameter_vals[3][p_best][s][r_best]);   // fit quality
        }
    }
    fclose(excitingresults);

    // Then G_inf values at best fit k_D
    err = fopen_s(&excitingresults, G_inf_fit_q_filename.c_str(), "w");
    if (err || !excitingresults) {
        cout << "Unable to save fit quality values\n";
    }
    else {
        for (int s = 1; s < (num_G_inf_vals - 1); s++) {
            fprintf(excitingresults, "%E\t", parameter_vals[2][p_best][q_best][s]);   // G_inf
            fprintf(excitingresults, "%E\n", parameter_vals[3][p_best][q_best][s]);   // fit quality
        }
    }
    fclose(excitingresults);

    return 0;
}