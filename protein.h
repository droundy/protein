#pragma once
using namespace std;
#include <stdio.h>
#include <string>
#include <math.h>
#include <cassert>
#include "weights.h"

extern double dx; //grid spacing
//number of gridpoints in each direction
extern int Nx, Ny, Nz;

const double difD = 2.5; // (um)^2 s^- 1
const double difE = 2.5; // (um)^2 s^-1
const double rate_ADP_ATP = 1; // s^-1
const double rate_D = .025; // um s^-1
const double rate_dD = .0015; // (um)^3 s^-1
const double rate_de = .7; // s^-1
const double rate_E = .093; // (um)^3 s^-1

enum reaction {ADP_to_ATP, DE_to_ADP_E, ATP_to_D, E_D_to_DE, X_ADP_pos, X_ADP_neg, Y_ADP_pos, Y_ADP_neg, Z_ADP_pos,
               Z_ADP_neg, X_ATP_pos, X_ATP_neg, Y_ATP_pos, Y_ATP_neg, Z_ATP_pos, Z_ATP_neg, X_E_pos, X_E_neg,
               Y_E_pos, Y_E_neg, Z_E_pos, Z_E_neg};
const int num_pos_reactions = Z_E_neg+1;
const int d[6][3] =
  {
    {1,0,0},
    {-1,0,0},
    {0,1,0},
    {0,-1,0},
    {0,0,1},
    {0,0,-1},
  };

extern int min_xi, min_yi, min_zi;

extern double A, B, C, D;
extern std::string mem_f_shape; //cell shape argument
extern std::string sim_type; //type of simulation type - exact, stochastic full_array, or stochastic half_array

struct stoch_params {
  int xi;
  int yi;
  int zi;
  int reaction;
};

//extern double *mem_A;

const int starting_num_guassians=20;
const int random_num_guassians=5;
const double Norm = 15.0; //This is the height of the guassians that make the cell wall
extern double guass[3*starting_num_guassians]; //stores y,z, and sigma for each guassian when creating random cell wall
extern int rand_seed; //=14; at this point I have this passed in from the command line as the D argument



//protein_utils:
void trim_grid(double **pointer_to_mem_A, double *first_mem_A);
void sym_check (double * mem_A);
std::string triangle_section(double y, double z);
void test_the_amount_of_area(double *first_mem_A, std::string mem_f_shape);
void set_insideArr(bool *insideArr);
bool inside(int xi, int yi, int zi);



//protein_weights:
void update_densities_and_weighting_for_reaction(stoch_params p, weights *ws, bool *insideArr, int *N_ATP, int *N_ADP,
                                    int *N_E, int *ND_st, int *NDE_st, double *mem_A);
void update_densities_and_weighting_for_diffusion(stoch_params p, weights *ws, bool *insideArr, int *N_ATP, int *N_ADP,
                                    int *N_E, int *ND_st, int *NDE_st, double *mem_A);
void initialize_densities_and_weighting(weights *ws, int *N_ATP, int *N_ADP, int *N_E,
                                        int *ND_st, int *NDE_st, double *mem_A);



//protein_membrane:
void randomize_cell_wall(double guass[]);
double f_2D_TIE_fighter(double y, double z);
double f_2D_triangle(double y, double z);
double f_2D_randst(double y, double z);
double f_2D_stad(double y, double z);
double mem_f(double x, double y, double z);
double find_intersection(const double fXYZ, const double fXYz, const double fXyZ, const double fxYZ,
                         const double fxyZ, const double fxYz, const double fXyz, const double fxyz,
                         const double f_minus_C, bool debugme = false);
void set_membrane(double mem_A[]);



//protein_microscopy
int set_density(double *nATP, double *nADP, double *nE, double *ND, double *NDE, int *ND_st, int *NDE_st, int *N_ATP, int *N_ADP, int *N_E, double *mem_A);
int get_next_density(double *mem_A, bool *insideArr, double *nATP, double *nADP,
                     double *nE, double *Nd, double *Nde, double *NflD, double *NflE,
                     double *JxATP, double *JyATP, double *JzATP,
                     double *JxADP, double *JyADP, double *JzADP,
                     double *JxE, double *JyE, double *JzE);
int get_J(double difD, double *nATP, double *nADP, double *nE,
          double *JxATP, double *JyATP, double *JzATP,
          double *JxADP, double *JyADP, double *JzADP,
          double *JxE, double *JyE, double *JzE);

