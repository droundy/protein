#include <stdio.h>
#include "weights.h"
struct stoch_params;


int set_density(double *nATP, double *nADP, double *nE, double *ND, double *NDE, int *ND_st, int *NDE_st, int *N_ATP, int *N_ADP, int *N_E, double *mem_A);
int get_next_stochastic_state(double *mem_A, bool *insideArr, double *nATP, double *nADP,
                   double *nE, double *Nd, double *Nde, double *NflD, double *NflE,
                   double *JxATP, double *JyATP, double *JzATP,
                   double *JxADP, double *JyADP, double *JzADP,
                              double *JxE, double *JyE, double *JzE);
int get_next_density(double *mem_A, bool *insideArr, double *nATP, double *nADP,
                     double *nE, double *Nd, double *Nde, double *NflD, double *NflE,
                     double *JxATP, double *JyATP, double *JzATP,
                     double *JxADP, double *JyADP, double *JzADP,
                     double *JxE, double *JyE, double *JzE);
int get_J(double difD, double *nATP, double *nADP, double *nE,
          double *JxATP, double *JyATP, double *JzATP,
          double *JxADP, double *JyADP, double *JzADP,
          double *JxE, double *JyE, double *JzE);
void update_densities_and_weighting_for_reaction(stoch_params p, weights *ws, int *N_ATP, int *N_ADP,
                                    int *N_E, int *ND_st, int *NDE_st, double *mem_A);
void update_densities_and_weighting_for_diffusion(stoch_params p, weights *ws, int *N_ATP, int *N_ADP,
                                    int *N_E, int *ND_st, int *NDE_st, double *mem_A);
void initialize_densities_and_weighting(weights *ws, int *N_ATP, int *N_ADP, int *N_E,
                                        int *ND_st, int *NDE_st, double *mem_A);
void set_membrane(double mem_A[]);
void set_curvature(double mem_A[], double curvature[]);
void set_insideArr(bool *insideArr);
bool inside(int xi, int yi, int zi);
double find_intersection(const double fXYZ, const double fXYz, const double fXyZ, const double fxYZ,
                         const double fxyZ, const double fxYz, const double fXyz, const double fxyz,
                         const double f_minus_C, bool debugme = false);
