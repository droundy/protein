#include <iostream>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "protein.h"
#include "MersenneTwister.h"
#include "weights.h"
#include <unistd.h>

const bool stochastic = true;

const int stoch_iter = 1000; //this is just a placeholder, not sure how we want to do time yet

const double nATP_starting_density = 1000.0/(M_PI*0.5*0.5); //proteins per micrometer^3 (values from paper)
const double nE_starting_density = 350.0/(M_PI*0.5*0.5); // proteins per micrometer
double density_factor;

int area_rating_flag = 0;
int slice_flag = 0;
int dump_flag = 0;
int debug_flag = 0;

char* hires_flag_str = new char[1024];
char* slice_flag_str = new char[1024];
char* debug_flag_str = new char[1024];


double tot_time; //total simulation time
double time_step; //simulation time step
int iter; //total # of iterations
int printout_iterations; //iterations between file printout

int box_divider_left;
int box_divider_right;

int Nx, Ny, Nz;
int min_xi, min_yi, min_zi;
double dx;

double x, y, z;

string mem_f_shape; //cell shape argument
string sim_type; //type of simulation type - exact, stochastic full_array, or stochastic half_array
double A, B, C, D; //specific shape parameters, set by command line args

//N denotes protein number, n denotes protein number density
double *nATP; //min D bound to an ATP
double *nADP; //min D bound to an ADP
double *nE; //loose min E in cytoplasm
double *ND; //min D bound to ATP on the wall
double *NDE; //min D bound to ATP and min E on the wall
int *s_N_ATP; //min D bound to an ATP integers for stochastic
int *s_N_ADP; //min D bound to an ADP integers for stochastic
int *s_N_E; //loose min E in cytoplasm integers for stochastic
int *s_ND; //min D bound to ATP on the wall
int *s_NDE; //min D bound to ATP and min E on the wall
double *NflD;
double *NflE;
double *f_mem;


double guass[3*starting_num_guassians]; //stores y,z, and sigma for each guassian when creating random cell wall
int rand_seed = 0; //=14; at this point I have this passed in from the command line as the D argument


 //struct for storing plot information
struct protein {
  char* name;

  //time_map
  double* sum;

  //box_plot
  double* numLeft;
  double* numMid;
  double* numRight;

  double* numRightUp;
  double* numRightDown;
  double* numLeftUp;
  double* numLeftDown;

  //arrow_plot
  double* maxval;
  int* ymax;
  int* zmax;
};

char* print_filename(const char* plotname, const char* proteinname) {
  char* filename = new char[1024];
  sprintf(filename,"data/shape-%s/%s%s%s%s-%s-%s-%1.2f-%1.2f-%1.2f-%1.2f-%1.2f-%s.dat",mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,plotname,proteinname,mem_f_shape.c_str(),A,B,C,D,density_factor,sim_type.c_str());
  return filename;
}

void count_and_print_proteins(int iteration, int *s_N_ATP, int *s_N_ADP, int *s_N_E, int *s_ND,
                              int *s_NDEs_N_ATP, double *nATP, double *nADP, double *nE, double *ND,
                              double *NDE, double *NflD, double *NflE, double *mem_A, bool *insideArr);


stoch_params index_to_parameters(int index) {
  stoch_params p;
  //should I worry about these floors when the integer is right on?
  p.reaction = floor( index/(Nx*Ny*Nz) );
  p.xi = floor( (index - p.reaction*Nx*Ny*Nz) / (Ny*Nz) );
  p.yi = floor( (index - p.reaction*Nx*Ny*Nz - p.xi*Ny*Nz) / Nz);
  p.zi = index - p.reaction*Nx*Ny*Nz - p.xi*Ny*Nz - p.yi*Nz;
  return p;
}


const char *reaction_name(int reaction) {
  switch (reaction) {
    case ADP_to_ATP:
      return "ADP_to_ATP";
    case DE_to_ADP_E:
      return "DE_to_ADP_E";
    case ATP_to_D:
      return "ATP_to_D";
    case E_D_to_DE:
      return "E_D_to_DE";

    case X_ADP_pos:
      return "X_ADP_pos";
    case Y_ADP_pos:
      return "Y_ADP_pos";
    case Z_ADP_pos:
      return "Z_ADP_pos";
    case X_ADP_neg:
      return "X_ADP_neg";
    case Y_ADP_neg:
      return "Y_ADP_neg";
    case Z_ADP_neg:
      return "Z_ADP_neg";

    case X_ATP_pos:
      return "X_ATP_pos";
    case Y_ATP_pos:
      return "Y_ATP_pos";
    case Z_ATP_pos:
      return "Z_ATP_pos";
    case X_ATP_neg:
      return "X_ATP_neg";
    case Y_ATP_neg:
      return "Y_ATP_neg";
    case Z_ATP_neg:
      return "Z_ATP_neg";

    case X_E_pos:
      return "X_E_pos";
    case Y_E_pos:
      return "Y_E_pos";
    case Z_E_pos:
      return "Z_E_pos";
    case X_E_neg:
      return "X_E_neg";
    case Y_E_neg:
      return "Y_E_neg";
    case Z_E_neg:
      return "Z_E_neg";
    default:
      printf("I am crashing syou silly\n");
      exit(1);
  }
}


int main (int argc, char *argv[]) {
  //command line parameters
  mem_f_shape = argv[1];
  A = atof(argv[2]);
  B = atof(argv[3]);
  C = atof(argv[4]);
  D = atof(argv[5]);
  density_factor = atof(argv[6]);
  sim_type = argv[8];
  dx=.05;

  memset(hires_flag_str,0,1024*sizeof(char));
  memset(slice_flag_str,0,1024*sizeof(char));
  memset(debug_flag_str,0,1024*sizeof(char));

  //flag checking
  for (int i=0; i<argc; i++) {
    if (strcmp(argv[i],"-area")==0) {
      area_rating_flag = 1;
      printf("Area rating printout.\n");
    }
    if (strcmp(argv[i],"-hires")==0) {
      dx=.025;
      printf("Using high resolution.\n");
      sprintf(hires_flag_str,"hires-");
    }
    if (strcmp(argv[i],"-slice")==0) {
      slice_flag = 1;
      printf("Printing middle slice data.\n");
      sprintf(slice_flag_str,"slice-");
    }
    if (strcmp(argv[i],"-dump")==0) {
      dump_flag = 1;
      printf("Printing all 501 data files.\n");
    }
    if (strcmp(argv[i],"-debug")==0) {
      dx = .05;
      debug_flag = 1;
      sprintf(debug_flag_str,"debug-");
      printf("=============================================\nDebug mode. dx=.15 um^3, tot_time=250s\n=============================================\n");
    }
  }

  //compute grid size based on cell parameters
  if (mem_f_shape=="p") {
    Nx = round(2*B/dx) + 4;
    Ny = round(2*B/dx) + 4;
    Nz = round((A + 2*B)/dx) + 4;
    printf("Nx = %d Ny = %d Nz = %d\n",Nx,Ny,Nz);
    box_divider_left = int(Nz/3);
    box_divider_right = int(2*Nz/3);
  }
  if (mem_f_shape=="stad") {
    Nx = round(A/dx) + 3;
    Ny = round((2*B+2*A)/dx) + 3;
    Nz = round((C+2*D+2*A)/dx) + 3;
  }
  if (mem_f_shape=="randst") {
    Nx = round(A/dx) + 5;
    Ny = round(B/dx) + 5;
    Nz = round(C/dx) + 5;
    rand_seed = int(D);
    if (D != round(D)) {
      printf("WARNING!!! When using randst the last argument, the rand_seed, should be an integer!  For now I've truncated it!!!\n");
    }
  }
  if (mem_f_shape=="TIE_fighter") {
    Nx = ceil(A/dx) + 4;
    Ny = ceil(B/dx) + 4;
    Nz = ceil(C/dx) + 4;
  }
  if (mem_f_shape=="triangle") {
    Nx = ceil(A/dx) + 5;
    Nz = ceil((A+B)/dx) + 5;
    //Using law of cosines we get height of triangle:
    double theta = acos((B*B+D*D-C*C)/(2*B*D));
    Ny = ceil((A+D*sin(theta))/dx) + 4;
  }
  if (mem_f_shape=="sp") {
    Nx = ceil(2*A/dx) + 5;
    Ny = ceil(2*A/dx) + 5;
    Nz = ceil(2*A/dx) + 5;
  }
  if (mem_f_shape=="e") {
    Nx = ceil(1/dx) + 4;
    Ny = ceil(2*A/dx) + 4;
    Nz = ceil(2*B/dx) + 4;
  }

  //fixed simulation parameters
  tot_time = 4000; //sec
  if (debug_flag==1) {
    tot_time = 15;
  }
  time_step = .1*dx*dx/difD;//sec
  iter = int(tot_time/time_step);//this will give us triangle data in about two days and randst in four days?
  printout_iterations = int(0.5/time_step);
  int total_arrow_time = 1000;
  int arrow_iter = int(total_arrow_time/time_step);//this tells us how many seconds we'll do the arrow plots for
  int print_denominator = 1000; //This dictates how often (in iterations) you add a line to the box, time, and ave plots
  printf("printout iterations = %d\n",printout_iterations);//approximately 5 seconds between each printout
  double dV = dx*dx*dx;

  printf("Simulation arguments: %s %g %g %g %g %g\n",mem_f_shape.c_str(),A,B,C,D,density_factor);

  //In the following, for every set of three numbers, the 1st is y and he 2nd is z and the 3rd is quassian width
  double guass92[] = {2.0,2.0,.3,2.0,2.5,.15,2.0,1.5,.15,2.5,2.0,.15,1.5,2.0,.15};
  double guass93[] = {1.0,1.0,.4,1.0,1.8,.2,1.0,2.6,.4,1.8,2.6,.2,2.6,2.6,.4,2.6,1.8,.2,2.6,1.0,.4,1.8,1.0,.2};
  double guass94[] = {2.6,3.2,.3,2.28,2.75,.25,3.0,3.9,.6,3.1,3.6,.4,3.3,3.9,.4,3.5,4.7,.5,2.9,5.6,.5,3.1,5.2,.4,2.4,5.8,.2,3.6,5.1,.3};
  double guass95[] = {2.2,2.4,.3,2.5,3.2,.6,2.7,3.5,.4,2.9,3.5,.4,3.5,4.2,.5,3.8,4.1,.8,3.1,4.6,.6,3.15,4.3,.5};
  double guass96[] = {1.3,1.3,.7,2.1,2,.7,3,2,.7,3.9,2,.7,4.7,1.3,.7,4,2.1,.7,4,3,.7,4,3.9,.7,4.7,4.7,.7,3.9,4,.7,3,4,.7,2.3,4,.7,1.3,4.7,.7,2.1,3.9,.7,3,3.9,.7,2.1,3.9,.7};
  double guass97[] = {1.4,3,.4,1.8,3,.4,2.2,3,.4,2.6,3,.4,3,3,.4,3.4,3,.4,3.8,3,.4,4.2,3,.4,4.6,3,.4,5,3,.4,5.4,3,.4,3.4,2.4,.6};
  double guass98[] = {2.0,2.0,.3,3,3,.6,4.2,3.4,.3,4.6,4.6,.6,3.4,5.6,.6};
  double guass99[] = {2.0,2.2,.5,3,3,.50,4.0,3.6,.50,3,4.2,.50,2.0,5,.5};

  double size_modifier_92 = 1.5;
  double size_modifier_93 = 1.5;
  double size_modifier_94 = 0.3;
  double size_modifier_95 = 0.3;
  double size_modifier_96 = 1.0;
  double size_modifier_97 = 1.0;
  double size_modifier_98 = 1.0;
  double size_modifier_99 = 1.3;

  {
    char *hn = (char *)malloc(1024);
    gethostname(hn, 1023);
    printf("Hostname:  %s\n", hn);
    free(hn);
  }
  printf("\nsize_modifier 92 is = %g",size_modifier_92);
  printf("\nsize_modifier 93 is = %g",size_modifier_93);
  printf("\nsize_modifier 94 is = %g",size_modifier_94);
  printf("\nsize_modifier 95 is = %g\n",size_modifier_95);
  printf("\nsize_modifier 99 is = %g\n",size_modifier_99);
  fflush(stdout);

  bzero(guass,int(sizeof(guass)/sizeof(double)));
  if (rand_seed == 92) {
    printf("Hello? rand_seed = %d\n",rand_seed);
    for (int i=0;i<int(sizeof(guass92)/sizeof(double));i++){
      guass[i] = size_modifier_92*guass92[i];
    }
  }
  else if (rand_seed == 93) {
    for (int i=0;i<int(sizeof(guass93)/sizeof(double));i++){
      guass[i] = size_modifier_93*guass93[i];
    }
  }
  else if (rand_seed == 94) {
    for (int i=0;i<int(sizeof(guass94)/sizeof(double));i++){
      guass[i] = size_modifier_94*guass94[i];
    }
  }
  else if (rand_seed == 95) {
    for (int i=0;i<int(sizeof(guass95)/sizeof(double));i++){
      guass[i] = size_modifier_95*guass95[i];
    }
  }
  else if (rand_seed == 96) {
    for (int i=0;i<int(sizeof(guass96)/sizeof(double));i++){
      guass[i] = size_modifier_96*guass96[i];
    }
  }
  else if (rand_seed == 97) {
    for (int i=0;i<int(sizeof(guass97)/sizeof(double));i++){
      guass[i] = size_modifier_97*guass97[i];
    }
  }
  else if (rand_seed == 98) {
    for (int i=0;i<int(sizeof(guass98)/sizeof(double));i++){
      guass[i] = size_modifier_98*guass98[i];
    }
  }
  else if (rand_seed == 99) {
    for (int i=0;i<int(sizeof(guass99)/sizeof(double));i++){
      guass[i] = size_modifier_99*guass99[i];
    }
  }
  else {
    if (mem_f_shape == "randst"){
      randomize_cell_wall(guass);
    }
  }
  printf("Those are all the guassians!\n");
  //end random stuff

  double *first_mem_A = new double[Nx*Ny*Nz];
  set_membrane(first_mem_A);

  bool *first_insideArr = new bool[Nx*Ny*Nz];
  set_insideArr(first_insideArr);

  if (mem_f_shape == "p" || mem_f_shape == "stad" || mem_f_shape == "sp"){
    test_the_amount_of_area(first_mem_A,mem_f_shape);
  }

  ///////////////?
  if (mem_f_shape == "stad"){
    printf("Trying to do stadium shape.  Not coded propperly to do that yet\n");
    exit(1);
  }
  ///////////////


  //Trimming the grid
  double *mem_A;
  double **pointer_to_mem_A;
  pointer_to_mem_A = &mem_A;
  bool *insideArr;
  bool **pointer_to_insideArr;
  pointer_to_insideArr = &insideArr;
  printf("\nBefore Trim Nx = %d Ny = %d Nz = %d\n",Nx,Ny,Nz);
  trim_grid(pointer_to_mem_A, first_mem_A, pointer_to_insideArr, first_insideArr);
  printf("\nAfter trim, min_xi = %d min_yi = %d min_zi = %d\n",min_xi,min_yi,min_zi);
  fflush(stdout);

  ///////////////////////////////////////////////
  // //Print outs comparing mem_A and insideArr:
  // for (int i=0;i<Nx;i++){
  //   printf("\ninsideArr:\n");
  //   for (int j=0;j<Ny;j++){
  //     for (int k=0;k<Nz;k++){
  //       if (insideArr[i*Ny*Nz+j*Nz+k]){
  //         printf("%d ",1);
  //       }
  //       else {
  //         printf("%d ",0);
  //       }
  //     }
  //     printf("\n");
  //   }
  //   printf("Membrane:\n");
  //   for (int j=0;j<Ny;j++){
  //     for (int k=0;k<Nz;k++){
  //       if (mem_A[i*Ny*Nz+j*Nz+k] != 0.0){
  //         printf("%d ",1);
  //       }
  //       else {
  //         printf("%d ",0);
  //       }
  //     }
  //     printf("\n");
  //   }
  // }
  // fflush(stdout);
  // exit(1);
  ///////////////////////////////////////////////

  //global arrays for storing simulation data
  nATP = new double[Nx*Ny*Nz];
  nADP = new double[Nx*Ny*Nz];
  nE = new double[Nx*Ny*Nz];
  s_N_ATP = new int[Nx*Ny*Nz];
  s_N_ADP = new int[Nx*Ny*Nz];
  s_N_E = new int[Nx*Ny*Nz];
  s_ND = new int[Nx*Ny*Nz];
  s_NDE = new int[Nx*Ny*Nz];
  ND = new double[Nx*Ny*Nz];
  NDE = new double[Nx*Ny*Nz];
  NflD = new double[Nx*Ny*Nz];
  NflE = new double[Nx*Ny*Nz];
  f_mem = new double[Nx*Ny*Nz];
  double *JxATP = new double[Nx*Ny*Nz];
  double *JyATP = new double[Nx*Ny*Nz];
  double *JzATP = new double[Nx*Ny*Nz];
  double *JxADP = new double[Nx*Ny*Nz];
  double *JyADP = new double[Nx*Ny*Nz];
  double *JzADP = new double[Nx*Ny*Nz];
  double *JxE = new double[Nx*Ny*Nz];
  double *JyE = new double[Nx*Ny*Nz];
  double *JzE = new double[Nx*Ny*Nz];

  const int numProteins = 7;

  protein* nATP_plot = new protein;
  protein* nE_plot = new protein;
  protein* nADP_plot = new protein;
  protein* NDE_plot = new protein;
  protein* ND_plot = new protein;
  protein* NflD_plot = new protein;
  protein* NflE_plot = new protein;

  protein* proteinList[numProteins] = { nATP_plot, nADP_plot, nE_plot, ND_plot, NDE_plot, NflD_plot, NflE_plot };
  double* protein_arrs[numProteins] = { nATP, nADP, nE, ND, NDE, NflD, NflE };

  //initialize things
  for (int pNum=0; pNum<numProteins; pNum++) {
    int total_print_iter = iter/print_denominator+2;
    proteinList[pNum]->sum = new double[Ny*Nz];
    proteinList[pNum]->name = new char[1024];

    proteinList[pNum]->numLeft = new double[total_print_iter];
    proteinList[pNum]->numMid = new double[total_print_iter];
    proteinList[pNum]->numRight = new double[total_print_iter];

    proteinList[pNum]->numRightUp = new double[total_print_iter];
    proteinList[pNum]->numRightDown = new double[total_print_iter];
    proteinList[pNum]->numLeftUp = new double[total_print_iter];
    proteinList[pNum]->numLeftDown = new double[total_print_iter];

    proteinList[pNum]->maxval = new double[arrow_iter];
    proteinList[pNum]->ymax = new int[arrow_iter];
    proteinList[pNum]->zmax = new int[arrow_iter];

    bzero(proteinList[pNum]->sum,Ny*Nz*sizeof(double));
    bzero(proteinList[pNum]->name,1024*sizeof(char));

    bzero(proteinList[pNum]->numLeft,total_print_iter*sizeof(double));
    bzero(proteinList[pNum]->numMid,total_print_iter*sizeof(double));
    bzero(proteinList[pNum]->numRight,total_print_iter*sizeof(double));

    bzero(proteinList[pNum]->numRightUp,total_print_iter*sizeof(double));
    bzero(proteinList[pNum]->numRightDown,total_print_iter*sizeof(double));
    bzero(proteinList[pNum]->numLeftUp,total_print_iter*sizeof(double));
    bzero(proteinList[pNum]->numLeftDown,total_print_iter*sizeof(double));

    bzero(proteinList[pNum]->maxval,arrow_iter*sizeof(double));
    bzero(proteinList[pNum]->ymax,arrow_iter*sizeof(int));
    bzero(proteinList[pNum]->zmax,arrow_iter*sizeof(int));
  }

  sprintf(proteinList[0]->name,"D_nATP");
  sprintf(proteinList[1]->name,"D_nADP");
  sprintf(proteinList[2]->name,"E_nE");
  sprintf(proteinList[3]->name,"D_ND");
  sprintf(proteinList[4]->name,"D_E_NDE");
  sprintf(proteinList[5]->name,"NflD");
  sprintf(proteinList[6]->name,"NflE");


  set_density(nATP, nADP, nE, ND, NDE, s_ND, s_NDE, s_N_ATP, s_N_ADP, s_N_E, mem_A, insideArr);

  //Starting the Sections file set up
  double left_area_total = 0;
  double middle_area_total = 0;
  double right_area_total = 0;
  double right_up_area_total = 0;
  double right_down_area_total = 0;
  double left_up_area_total = 0;
  double left_down_area_total = 0;


  if (mem_f_shape=="triangle") {
    char* outfilename_sections = new char[1024];
    sprintf(outfilename_sections, "data/shape-%s/membrane_files/%s%s%ssections-%s-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),
            debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor);
    FILE *outfile_sections = fopen((const char*)outfilename_sections,"w");
    for (int j=0;j<Ny;j++){
      for (int i=0;i<Nz;i++) {
        double marker = 0;
        if (triangle_section(j*dx,i*dx)=="Right"){
          marker = 1;
          for (int a=0;a<Nx;a++){
            right_area_total += mem_A[a*Ny*Nz+j*Nz+i];
          }
        }
        if (triangle_section(j*dx,i*dx)=="Mid"){
          marker = 2;
          for (int a=0;a<Nx;a++){
            middle_area_total += mem_A[a*Ny*Nz+j*Nz+i];
          }
        }
        if (triangle_section(j*dx,i*dx)=="Left"){
          marker = 3;
          for (int a=0;a<Nx;a++){
            left_area_total += mem_A[a*Ny*Nz+j*Nz+i];
          }
        }
        if (inside(int(Nx/2),j,i)==false) {
          marker = 0;
        }
        fprintf(outfile_sections,"%g ",marker);
      }
      fprintf(outfile_sections,"\n");
    }
    fclose(outfile_sections);
    delete[] outfilename_sections;
  }

  if (mem_f_shape=="p" || mem_f_shape=="stad"){
    box_divider_left = int(Nz/3);
    box_divider_right = int(2*Nz/3);
    char* outfilename_sections = new char[1024];
    sprintf(outfilename_sections, "data/shape-%s/membrane_files/%s%s%ssections-%s-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),
            debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor);
    FILE *outfile_sections = fopen((const char*)outfilename_sections,"w");
    for (int a=0; a<Ny; a++) {
      for (int b=0; b<Nz; b++) {
        double marker = 0;
        if (b < box_divider_left) {
          marker = 3;
          for (int c=0;c<Nx;c++){
            left_area_total += mem_A[c*Ny*Nz+a*Nz+b];
          }
        }
        if (b > box_divider_right) {
          marker = 1;
          for (int c=0;c<Nx;c++){
            right_area_total += mem_A[c*Ny*Nz+a*Nz+b];
          }
        }
        if ((b <= box_divider_right) && (b >= box_divider_left)) {
          marker = 2;
          for (int c=0;c<Nx;c++){
            middle_area_total += mem_A[c*Ny*Nz+a*Nz+b];
          }
        }
        if (inside(int(Nx/2),a,b)==false) {
          marker = 0;
        }
        fprintf(outfile_sections, "%g ",marker);
      }
      fprintf(outfile_sections,"\n");
    }
    fclose(outfile_sections);
    delete[] outfilename_sections;
  }

  double vert_div = 0;
  double vert_div_two = 0;
  double hor_div = 0;
  double hor_div_two = 0;

  if (mem_f_shape=="randst") {
    if (rand_seed == 92) {
      vert_div = size_modifier_92*(1.2/dx)-min_zi+1;
      vert_div_two = size_modifier_92*(2.7/dx)-min_zi+1;
    }
    else if (rand_seed == 93) {
      vert_div = size_modifier_93*(1.2/dx)-min_zi+1;
      vert_div_two = size_modifier_93*(2.7/dx)-min_zi+1;
    }
    else if (rand_seed == 94) {
      vert_div = size_modifier_94*(3.6/dx)-min_zi+1;
      vert_div_two = size_modifier_94*(5.0/dx)-min_zi+1;
    }
    else if (rand_seed == 95) {
      vert_div = size_modifier_95*(3.3/dx)-min_zi+1;
      vert_div_two = size_modifier_95*(4.4/dx)-min_zi+1;
    }
    else if (rand_seed == 96) {
      vert_div = size_modifier_96*(2.8/dx)-min_zi+1;
      hor_div = size_modifier_96*(3.0/dx)-min_yi+1;
    }
    else if (rand_seed == 97) {
      vert_div = size_modifier_97*(2.4/dx)-min_zi+1;
      hor_div = size_modifier_97*(2.45/dx)-min_yi+1;
      hor_div_two = size_modifier_97*(4.6/dx)-min_yi+1;
    }
    else if (rand_seed == 98) {
      vert_div = size_modifier_98*(3.0/dx)-min_zi+1;
      vert_div_two = size_modifier_98*(4.8/dx)-min_zi+1;
    }
    else if (rand_seed == 99) {
      vert_div = size_modifier_99*(2.6/dx)-min_zi+1;
      vert_div_two = size_modifier_99*(4.6/dx)-min_zi+1;
    }
  }
  if (mem_f_shape=="randst") {
    char* outfilename_sections = new char[1024];
    sprintf(outfilename_sections, "data/shape-%s/membrane_files/%s%s%ssections-%s-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),
            debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor);
    FILE *outfile_sections = fopen((const char*)outfilename_sections,"w");
    for (int a=0; a<Ny; a++) {
      for (int b=0; b<Nz; b++) {
        double marker = 0;
        if (rand_seed == 99 || rand_seed == 98 || rand_seed == 95 || rand_seed == 94 || rand_seed == 93 || rand_seed == 92) {
          if (b < vert_div) {
            marker = 3;
            for (int c=0;c<Nx;c++){
              left_area_total += mem_A[c*Ny*Nz+a*Nz+b];
            }
          }
          else if (b > vert_div_two) {
            marker = 1;
            for (int c=0;c<Nx;c++){
              right_area_total += mem_A[c*Ny*Nz+a*Nz+b];
            }
          }
          else {
            marker = 2;
            for (int c=0;c<Nx;c++){
              middle_area_total += mem_A[c*Ny*Nz+a*Nz+b];
            }
          }
        }
        if (rand_seed == 97) {
          if (b < vert_div) {
            marker = 4;
            for (int c=0;c<Nx;c++){
              left_area_total += mem_A[c*Ny*Nz+a*Nz+b];
            }
          }
          else if (a > hor_div_two) {
            marker = 1;
            for (int c=0;c<Nx;c++){
              right_up_area_total += mem_A[c*Ny*Nz+a*Nz+b];
            }
          }
          else if (a < hor_div) {
            marker = 3;
            for (int c=0;c<Nx;c++){
              right_down_area_total += mem_A[c*Ny*Nz+a*Nz+b];
            }
          }
          else {
            marker = 2;
            for (int c=0;c<Nx;c++){
              middle_area_total += mem_A[c*Ny*Nz+a*Nz+b];
            }
          }
        }
        if (rand_seed == 96) {
          if (b < vert_div && a < hor_div) {
            marker = 3;
            for (int c=0;c<Nx;c++){
              left_down_area_total += mem_A[c*Ny*Nz+a*Nz+b];
            }
          }
          else if (b < vert_div && a >= hor_div) {
            marker = 4;
            for (int c=0;c<Nx;c++){
              left_up_area_total += mem_A[c*Ny*Nz+a*Nz+b];
            }
          }
          else if (b >= vert_div && a < hor_div) {
            marker = 1;
            for (int c=0;c<Nx;c++){
              right_down_area_total += mem_A[c*Ny*Nz+a*Nz+b];
            }
          }
          else {
            marker = 2;
            for (int c=0;c<Nx;c++){
              right_up_area_total += mem_A[c*Ny*Nz+a*Nz+b];
            }
          }
        }
        if (inside(int(Nx/2),a,b)==false) {
          marker = 0;
        }
        fprintf(outfile_sections, "%g ",marker);
      }
      fprintf(outfile_sections, "\n");
    }
    fclose(outfile_sections);
    delete[] outfilename_sections;
  }
  printf("\nSections file has printed\n\n");
  fflush(stdout);

  weights *ws;
  if (sim_type == "full_array") {
    ws = new weights(num_pos_reactions*Nx*Ny*Nz,'f');
    printf("\nSimulation type = full_array\n");
  }
  else if (sim_type == "half_array") {
    ws = new weights(num_pos_reactions*Nx*Ny*Nz,'h');
    printf("\nSimulation type = half_array\n");
  }
  else {
    ws = NULL;
    printf("\nSimulation type = exact\n");
  }
  printf("\nHave initialized weights object\n");
  fflush(stdout);

  if (sim_type != "exact") {
    initialize_densities_and_weighting(ws, insideArr, s_N_ATP, s_N_ADP, s_N_E, s_ND, s_NDE, mem_A);
    printf("\nHave initialized the values of the weights array\n");
    fflush(stdout);
  }

  /////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////

  double start_time = double(clock()/CLOCKS_PER_SEC);
  double time = start_time;

  //starting simuation
  double spill_over_time = 0;
  for (int i=0;i<iter;i++){
    int debug = 0;
    if (sim_type != "exact") {
      double elapsed_time = spill_over_time;
      while (elapsed_time < time_step) {
        int index = ws->lookup( (double)rand()/(RAND_MAX) );//random number is from 0 to 1
        stoch_params p = index_to_parameters(index);
        debug += 1;
        if (p.reaction <= E_D_to_DE) {
          switch (p.reaction) {
            case ADP_to_ATP:
              s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] -= 1;
              s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] += 1;
              break;
            case DE_to_ADP_E:
              s_NDE[p.xi*Ny*Nz+p.yi*Nz+p.zi] -= 1;
              s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] += 1;
              s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] += 1;
              break;
            case ATP_to_D:
              s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] -= 1;
              s_ND[p.xi*Ny*Nz+p.yi*Nz+p.zi] += 1;
              break;
            case E_D_to_DE:
              s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] -= 1;
              s_ND[p.xi*Ny*Nz+p.yi*Nz+p.zi] -= 1;
              s_NDE[p.xi*Ny*Nz+p.yi*Nz+p.zi] += 1;
              break;
            default:
              printf("reaction should be a reaction, but not any of the specific cases\n!");
              fflush(stdout);
              exit(1);
          }
          update_all_densities_and_weighting_for_changing_gridpt(ws, insideArr, s_N_ATP, s_N_ADP, s_N_E, s_ND,
                                                        s_NDE, mem_A, p.xi, p.yi, p.zi);
        }
        else if (p.reaction >= X_ADP_pos && p.reaction <= Z_ADP_neg) {
          s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] -= 1;
          update_all_densities_and_weighting_for_changing_gridpt(ws, insideArr, s_N_ATP, s_N_ADP, s_N_E, s_ND,
                                                                 s_NDE, mem_A, p.xi, p.yi, p.zi);
          int q = p.reaction - X_ADP_pos;
          s_N_ADP[ (p.xi+d[q][0])*Ny*Nz + (p.yi+d[q][1])*Nz + (p.zi+d[q][2]) ] += 1;
          update_all_densities_and_weighting_for_changing_gridpt(ws, insideArr, s_N_ATP, s_N_ADP, s_N_E, s_ND,
                                                        s_NDE, mem_A, p.xi+d[q][0], p.yi+d[q][1], p.zi+d[q][2]);
        }
        else if (p.reaction >= X_ATP_pos && p.reaction <= Z_ATP_neg) {
          s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] -= 1;
          update_all_densities_and_weighting_for_changing_gridpt(ws, insideArr, s_N_ATP, s_N_ADP, s_N_E, s_ND,
                                                                 s_NDE, mem_A, p.xi, p.yi, p.zi);
          int q = p.reaction - X_ATP_pos;
          s_N_ATP[ (p.xi+d[q][0])*Ny*Nz + (p.yi+d[q][1])*Nz + (p.zi+d[q][2]) ] += 1;
          update_all_densities_and_weighting_for_changing_gridpt(ws, insideArr, s_N_ATP, s_N_ADP, s_N_E, s_ND,
                                                        s_NDE, mem_A, p.xi+d[q][0], p.yi+d[q][1], p.zi+d[q][2]);
        }
        else if (p.reaction >= X_E_pos && p.reaction <= Z_E_neg) {
          s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] -= 1;
          update_all_densities_and_weighting_for_changing_gridpt(ws, insideArr, s_N_ATP, s_N_ADP, s_N_E, s_ND,
                                                                 s_NDE, mem_A, p.xi, p.yi, p.zi);
          int q = p.reaction - X_E_pos;
          s_N_E[ (p.xi+d[q][0])*Ny*Nz + (p.yi+d[q][1])*Nz + (p.zi+d[q][2]) ] += 1;
          update_all_densities_and_weighting_for_changing_gridpt(ws, insideArr, s_N_ATP, s_N_ADP, s_N_E, s_ND,
                                                        s_NDE, mem_A, p.xi+d[q][0], p.yi+d[q][1], p.zi+d[q][2]);
        }
        if (s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] < 0 || s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] < 0 || s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] < 0) {
          printf("There was a negative protein!! x %d,y %d,z %d,reaction %s, weight = %g, ADP = %d, ATP = %d, E = %d, wall = %g\n",
                 p.xi,p.yi,p.zi,reaction_name(p.reaction+i),ws->lookup_prob_for_specific_index((p.reaction+i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi),
                 s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi],s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi],s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi],mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] );
          fflush(stdout);
          exit(1);
        }
        elapsed_time -= log( (double)rand()/(RAND_MAX) ) / ws->get_total();
      }
      spill_over_time = elapsed_time - time_step;
      for (int j=0;j<Nx*Ny*Nz;j++) {
        nATP[j] = double(s_N_ATP[j])/(dx*dx*dx);
        nADP[j] = double(s_N_ADP[j])/(dx*dx*dx);
        nE[j] = double(s_N_E[j])/(dx*dx*dx);
        ND[j] = double(s_ND[j]);
        NDE[j] = double(s_NDE[j]);
        NflD[j] = (nATP[j] + nADP[j])*(dx*dx*dx) + ND[j] + NDE[j];
        NflE[j] = nE[j]*dx*dx*dx + NDE[j];
      }
    }
    else {
      get_J(difD, nATP, nADP, nE, JxATP, JyATP,
            JzATP, JxADP, JyADP, JzADP, JxE, JyE, JzE);
      get_next_density(mem_A, insideArr, nATP, nADP, nE, ND, NDE, NflD, NflE, JxATP, JyATP, JzATP,
                       JxADP, JyADP, JzADP, JxE, JyE, JzE);
    }
    if ( (sim_type == "exact" && i%100000 == 0) || (sim_type != "exact" && i%100000 == 0) ) {
      bool negative_protein = false;
      for(int xi=0;xi<Nx;xi++){
        for(int yi=0;yi<Ny;yi++){
          for(int zi=0;zi<Nz;zi++){
            if (s_N_E[xi*Ny*Nz+yi*Nz+zi] < 0 || s_N_ATP[xi*Ny*Nz+yi*Nz+zi] < 0 || s_N_ADP[xi*Ny*Nz+yi*Nz+zi] < 0) {
              printf("There is a negative protein at x %d, y %d, z %d! and its iteration %d\n",xi,yi,zi,i);
              negative_protein = true;
            }
          }
        }
      }
      if (negative_protein) {
        exit(1);
      }
      count_and_print_proteins(i,s_N_ATP, s_N_ADP, s_N_E, s_ND, s_NDE, nATP, nADP, nE, ND, NDE, NflD, NflE, mem_A, insideArr);
      time = double(clock()/CLOCKS_PER_SEC) - start_time;
      printf("This is a %s simulation and it took %g seconds to get to interation %d\n",
             sim_type.c_str(), time, i);
      printf("Which means it would take %g seconds to get 1000 simulation seconds\n",
             (double(time)/double(i))*(1000.0/time_step));
      fflush(stdout);
    }
    if (i%(int(100*print_denominator))==0) {
      printf("Finished sim loop # i=%d, We're %1.2f percent done\n",i,100*double(i)/iter);
    }
    //time map
    for (int pNum=0; pNum<numProteins; pNum++) {
      for (int a=0; a<Ny; a++) {
        for (int b=0; b<Nz; b++) {
          if (slice_flag==0) {
            for (int c=0; c<Nx; c++) {
              proteinList[pNum]->sum[a*Nz+b] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b];
            }
          }
          else {
            proteinList[pNum]->sum[a*Nz+b] += protein_arrs[pNum][int(Nx/2)*Ny*Nz+a*Nz+b];
          }
        }
      }

      //box plot ...
      if (i%print_denominator==0) {
        if ((strcmp(proteinList[pNum]->name,"D_nATP")==0) || (strcmp(proteinList[pNum]->name,"E_nE")==0) || (strcmp(proteinList[pNum]->name,"D_nADP")==0)) {
          dV = dx*dx*dx;
        }
        else {
          dV = 1;
        }
        int i_dat = i/print_denominator;
        for (int a=0; a<Ny; a++) {
          for (int b=0; b<Nz; b++) {
            for (int c=0; c<Nx; c++) {
              if (mem_f_shape == "p" || mem_f_shape == "stad") {
                if (b < box_divider_left) {
                  proteinList[pNum]->numLeft[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                }
                if (b > box_divider_right) {
                  proteinList[pNum]->numRight[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                }
                if ((b <= box_divider_right) && (b >= box_divider_left)) {
                  proteinList[pNum]->numMid[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                }
              }
              if (mem_f_shape == "randst") {
                if (rand_seed == 99 || rand_seed == 98 || rand_seed == 95 || rand_seed == 94 || rand_seed == 93 || rand_seed == 92) {
                  if (b < vert_div) {
                    proteinList[pNum]->numLeft[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else if (b > vert_div_two) {
                    proteinList[pNum]->numRight[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else {
                    proteinList[pNum]->numMid[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                }
                if (rand_seed == 97) {
                  if (b < vert_div) {
                    proteinList[pNum]->numLeft[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else if (a > hor_div_two) {
                    proteinList[pNum]->numRightUp[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else if (a < hor_div) {
                    proteinList[pNum]->numRightDown[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else {
                    proteinList[pNum]->numMid[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                }
                if (rand_seed == 96) {
                  if (b < vert_div && a < hor_div) {
                    proteinList[pNum]->numLeftDown[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else if (b < vert_div && a >= hor_div) {
                    proteinList[pNum]->numLeftUp[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else if (b >= vert_div && a < hor_div) {
                    proteinList[pNum]->numRightDown[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else {
                    proteinList[pNum]->numRightUp[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                }
              }
              if (mem_f_shape == "triangle") {
                if (triangle_section(a*dx,b*dx) == "Left") {
                  proteinList[pNum]->numLeft[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                }
                else if (triangle_section(a*dx,b*dx) == "Right") {
                  proteinList[pNum]->numRight[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                }
                else {
                  proteinList[pNum]->numMid[i_dat] += protein_arrs[pNum][c*Ny*Nz+a*Nz+b]*dV;
                }
              }
            }
          }
        }
      }

      //arrow plot
      if (i<arrow_iter) {
        double storemaxval = 0;
        double currentval;
        for (int a=0; a<Ny; a++) {
          for (int b=0; b<Nz; b++) {
            if (slice_flag==0) {
              currentval=0;
              for (int c=0; c<Nx; c++) {
                currentval += protein_arrs[pNum][c*Ny*Nz+a*Nz+b];
              }
              if (currentval > storemaxval) {
                storemaxval = currentval;
                proteinList[pNum]->maxval[i] = storemaxval;
                proteinList[pNum]->ymax[i] = a;
                proteinList[pNum]->zmax[i] = b;
              }
            }
            else {
              currentval = protein_arrs[pNum][int(Nx/2)*Ny*Nz+a*Nz+b];
              if (currentval > storemaxval) {
                storemaxval = currentval;
                proteinList[pNum]->maxval[i] = storemaxval;
                proteinList[pNum]->ymax[i] = a;
                proteinList[pNum]->zmax[i] = b;
              }
            }
          }
        }
      }
    }

    //begin file printing
    if ((dump_flag == 1) && (i%printout_iterations == 0)) {
      dV = dx*dx*dx;
      int k = i/printout_iterations;

      //begin nATP printing.
      char *outfilenameATP = new char[1024];
      sprintf(outfilenameATP, "data/shape-%s/%s%s%snATP-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d-%s.dat", mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor,k,sim_type.c_str());
      FILE *nATPfile = fopen((const char *)outfilenameATP,"w");
      delete[] outfilenameATP;

      if (slice_flag==1) {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            fprintf(nATPfile, "%1.2f ", nATP[(int(Nx/2))*Ny*Nz+a*Nz+b]);
          }
          fprintf(nATPfile, "\n");
        }
        fclose(nATPfile);
      }

      else {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            double nATPsum = 0;
            for (int c=0;c<Nx;c++){
              nATPsum += nATP[c*Ny*Nz+a*Nz+b];
            }
            fprintf(nATPfile, "%1.2f ", nATPsum);
          }
          fprintf(nATPfile, "\n");
        }
        fclose(nATPfile);
      }
      //end nATP printing

      //nE printing
      char *outfilenameE = new char[1000];
      sprintf(outfilenameE, "data/shape-%s/%s%s%snE-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d-%s.dat", mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor,k,sim_type.c_str());
      FILE *nEfile = fopen((const char *)outfilenameE,"w");
      delete[] outfilenameE;

      if (slice_flag==1) {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            fprintf(nEfile, "%1.2f ", nE[(int(Nx/2))*Ny*Nz+a*Nz+b]);
          }
          fprintf(nEfile, "\n");
        }
        fclose(nEfile);
      }

      else {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            double nEsum = 0;
            for (int c=0;c<Nx;c++){
              nEsum += nE[c*Ny*Nz+a*Nz+b];
            }
            fprintf(nEfile, "%1.2f ", nEsum);
          }
          fprintf(nEfile, "\n");
        }
        fclose(nEfile);
      }
      //end nE printing

      //nADP printing
      char *outfilenameADP = new char[1000];
      sprintf(outfilenameADP, "data/shape-%s/%s%s%snADP-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d-%s.dat", mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor,k,sim_type.c_str());
      FILE *nADPfile = fopen((const char *)outfilenameADP,"w");
      delete[] outfilenameADP;

      if (slice_flag==1) {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            fprintf(nADPfile, "%1.2f ", nADP[(int(Nx/2))*Ny*Nz+a*Nz+b]);
          }
          fprintf(nADPfile, "\n");
        }
        fclose(nADPfile);
      }

      else {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            double nADPsum = 0;
            for (int c=0;c<Nx;c++){
              nADPsum += nADP[c*Ny*Nz+a*Nz+b];
            }
            fprintf(nADPfile, "%1.2f ", nADPsum);
          }
          fprintf(nADPfile, "\n");
        }
        fclose(nADPfile);
      }
      //end nADP printing

      //begin ND printing
      char *outfilenameD = new char[1000];
      sprintf(outfilenameD, "data/shape-%s/%s%s%sND-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d-%s.dat", mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor,k,sim_type.c_str());
      FILE *NDfile = fopen((const char *)outfilenameD,"w");
      delete[] outfilenameD;

      if (slice_flag==1) {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            fprintf(NDfile, "%1.2f ", ND[(int(Nx/2))*Ny*Nz+a*Nz+b]);
          }
          fprintf(NDfile, "\n");
        }
        fclose(NDfile);
      }

      else {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            double NDsum = 0;
            for (int c=0;c<Nx;c++){
              NDsum += ND[c*Ny*Nz+a*Nz+b];
            }
            fprintf(NDfile, "%1.2f ", NDsum);
          }
          fprintf(NDfile, "\n");
        }
        fclose(NDfile);
      }
      //end ND printing

      //begin NDE printing
      char *outfilenameDE = new char[1000];
      sprintf(outfilenameDE, "data/shape-%s/%s%s%sNDE-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d-%s.dat", mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor,k,sim_type.c_str());
      FILE *NDEfile = fopen((const char *)outfilenameDE,"w");
      delete[] outfilenameDE;

      if (slice_flag==1) {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            fprintf(NDEfile, "%1.2f ", NDE[(int(Nx/2))*Ny*Nz+a*Nz+b]);
          }
          fprintf(NDEfile, "\n");
        }
        fclose(NDEfile);
      }

      else {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            double NDEsum = 0;
            for (int c=0;c<Nx;c++){
              NDEsum += NDE[c*Ny*Nz+a*Nz+b];
            }
            fprintf(NDEfile, "%1.2f ", NDEsum);
          }
          fprintf(NDEfile, "\n");
        }
        fclose(NDEfile);
      }
      //end NDE printing

      //begin NflE printing
      char *outfilenameflE = new char[1000];
      sprintf(outfilenameflE, "data/shape-%s/%s%s%sNflE-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d-%s.dat", mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor,k,sim_type.c_str());
      FILE *NflEfile = fopen((const char *)outfilenameflE,"w");
      delete[] outfilenameflE;

      if (slice_flag==1) {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            fprintf(NflEfile, "%1.2f ", nE[(int(Nx/2))*Ny*Nz+a*Nz+b]*dV + NDE[(int(Nx/2))*Ny*Nz+a*Nz+b]);
          }
          fprintf(NflEfile, "\n");
        }
        fclose(NflEfile);
      }

      else {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            double NflEsum = 0;
            for (int c=0;c<Nx;c++){
              NflEsum += nE[c*Ny*Nz+a*Nz+b]*dV + NDE[c*Ny*Nz+a*Nz+b];
            }
            fprintf(NflEfile, "%1.2f ", NflEsum);
          }
          fprintf(NflEfile, "\n");
        }
        fclose(NflEfile);
      }
      //end NflE printing

      //begin NflD printing
      char *outfilenameflD = new char[1000];
      sprintf(outfilenameflD, "data/shape-%s/%s%s%sNflD-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d-%s.dat", mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor,k,sim_type.c_str());
      FILE *NflDfile = fopen((const char *)outfilenameflD,"w");
      delete[] outfilenameflD;

      if (slice_flag==1) {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            fprintf(NflDfile, "%1.2f ", NDE[(int(Nx/2))*Ny*Nz+a*Nz+b] + nADP[(int(Nx/2))*Ny*Nz+a*Nz+b]*dV + nATP[(int(Nx/2))*Ny*Nz+a*Nz+b]*dV + ND[(int(Nx/2))*Ny*Nz+a*Nz+b]);
          }
          fprintf(NflDfile, "\n");
        }
        fclose(NflDfile);
      }

      else {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            double NflDsum = 0;
            for (int c=0;c<Nx;c++){
              NflDsum += NDE[c*Ny*Nz+a*Nz+b] + nADP[c*Ny*Nz+a*Nz+b]*dV + nATP[c*Ny*Nz+a*Nz+b]*dV + ND[c*Ny*Nz+a*Nz+b];
            }
            fprintf(NflDfile, "%1.2f ", NflDsum);
          }
          fprintf(NflDfile, "\n");
        }
        fclose(NflDfile);
      }
      //end NflD printing
      k++;
    }

    time_step = .1*dx*dx/difD;//sec
    int plot_denominator = 100000;
    int i_dat = i/print_denominator;
    if (i%plot_denominator==0){
      //boxplot
      char *boxname = print_filename("box-plot","");
      printf("\nPrinting the box plots %s\n",boxname);
      fflush(stdout);
      FILE* box_plot = fopen(boxname,"w");
      delete[] boxname;
      for (int pNum=0; pNum<numProteins; pNum++) {

        if (mem_f_shape == "p" || rand_seed == 99 || rand_seed == 98 || rand_seed == 95 || rand_seed == 94 || rand_seed == 93 || rand_seed == 92 ||
            mem_f_shape == "triangle" || mem_f_shape == "stad") {
          fprintf(box_plot,"%s\tleft\t",proteinList[pNum]->name);
          for (int i_plot_dat=0; i_plot_dat<i_dat+1; i_plot_dat++) {
            fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numLeft[i_plot_dat]));
          }
          fprintf(box_plot,"\n");
          fprintf(box_plot,"%s\tmid\t",proteinList[pNum]->name);
          for (int i_plot_dat=0; i_plot_dat<i_dat+1; i_plot_dat++) {
            fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numMid[i_plot_dat]));
          }
          fprintf(box_plot,"\n");
          fprintf(box_plot,"%s\tright\t",proteinList[pNum]->name);
          for (int i_plot_dat=0; i_plot_dat<i_dat+1; i_plot_dat++) {
            fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numRight[i_plot_dat]));
          }
          fprintf(box_plot,"\n");
          fprintf(box_plot,"\n");
        }
        if (rand_seed == 97) {
          fprintf(box_plot,"%s\tleft\t",proteinList[pNum]->name);
          for (int i_plot_dat=0; i_plot_dat<i_dat+1; i_plot_dat++) {
            fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numLeft[i_plot_dat]));
          }
          fprintf(box_plot,"\n");
          fprintf(box_plot,"%s\trightup\t",proteinList[pNum]->name);
          for (int i_plot_dat=0; i_plot_dat<i_dat+1; i_plot_dat++) {
            fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numRightUp[i_plot_dat]));
          }
          fprintf(box_plot,"\n");
          fprintf(box_plot,"%s\tmid\t",proteinList[pNum]->name);
          for (int i_plot_dat=0; i_plot_dat<i_dat+1; i_plot_dat++) {
            fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numMid[i_plot_dat]));
          }
          fprintf(box_plot,"\n");
          fprintf(box_plot,"%s\trightdown\t",proteinList[pNum]->name);
          for (int i_plot_dat=0; i_plot_dat<i_dat+1; i_plot_dat++) {
            fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numRightDown[i_plot_dat]));
          }
          fprintf(box_plot,"\n");
          fprintf(box_plot,"\n");
        }
        if (rand_seed == 96) {
          fprintf(box_plot,"%s\trightup\t",proteinList[pNum]->name);
          for (int i_plot_dat=0; i_plot_dat<i_dat+1; i_plot_dat++) {
            fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numRightUp[i_plot_dat]));
          }
          fprintf(box_plot,"\n");
          fprintf(box_plot,"%s\tleftup\t",proteinList[pNum]->name);
          for (int i_plot_dat=0; i_plot_dat<i_dat+1; i_plot_dat++) {
            fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numLeftUp[i_plot_dat]));
          }
          fprintf(box_plot,"\n");
          fprintf(box_plot,"%s\tleftdown\t",proteinList[pNum]->name);
          for (int i_plot_dat=0; i_plot_dat<i_dat+1; i_plot_dat++) {
            fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numLeftDown[i_plot_dat]));
          }
          fprintf(box_plot,"\n");
          fprintf(box_plot,"%s\trightdown\t",proteinList[pNum]->name);
          for (int i_plot_dat=0; i_plot_dat<i_dat+1; i_plot_dat++) {
            fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numRightDown[i_plot_dat]));
          }
          fprintf(box_plot,"\n");
          fprintf(box_plot,"\n");
        }
      }
      fclose(box_plot);

      for (int pNum=0; pNum<numProteins; pNum++) {
        char *avename = new char[1024];
        sprintf(avename,"%s",print_filename("ave_plot",""));
        FILE* ave_plot = fopen(avename,"w");

        if (mem_f_shape == "p" || mem_f_shape == "triangle" || rand_seed == 98 || rand_seed == 95 || rand_seed == 94 || rand_seed == 93 || rand_seed == 92 ||
            rand_seed == 99 || mem_f_shape == "stad"){
          for (int pNum=3; pNum<numProteins; pNum++) {
            fprintf(ave_plot,"%s\tleft\t",proteinList[pNum]->name);
            for (int i_plot_dat=0; i_plot_dat<i_dat; i_plot_dat++) {
              fprintf(ave_plot,"%1.2f\t",(proteinList[pNum]->numLeft[i_plot_dat]/left_area_total));
            }
            fprintf(ave_plot,"\n");
            fprintf(ave_plot,"%s\tmid\t",proteinList[pNum]->name);
            for (int i_plot_dat=0; i_plot_dat<i_dat; i_plot_dat++) {
              fprintf(ave_plot,"%1.2f\t",(proteinList[pNum]->numMid[i_plot_dat]/middle_area_total));
            }
            fprintf(ave_plot,"\n");
            fprintf(ave_plot,"%s\tright\t",proteinList[pNum]->name);
            for (int i_plot_dat=0; i_plot_dat<i_dat; i_plot_dat++) {
              fprintf(ave_plot,"%1.2f\t",(proteinList[pNum]->numRight[i_plot_dat]/right_area_total));
            }
            fprintf(ave_plot,"\n");
            fprintf(ave_plot,"\n");
          }
        }
        if (rand_seed == 97){
          for (int pNum=3; pNum<numProteins; pNum++) {
            fprintf(ave_plot,"%s\tleft\t",proteinList[pNum]->name);
            for (int i_plot_dat=0; i_plot_dat<i_dat; i_plot_dat++) {
              fprintf(ave_plot,"%1.2f\t",(proteinList[pNum]->numLeft[i_plot_dat]/left_area_total));
            }
            fprintf(ave_plot,"\n");
            fprintf(ave_plot,"%s\trightup\t",proteinList[pNum]->name);
            for (int i_plot_dat=0; i_plot_dat<i_dat; i_plot_dat++) {
              fprintf(ave_plot,"%1.2f\t",(proteinList[pNum]->numRightUp[i_plot_dat]/right_up_area_total));
            }
            fprintf(ave_plot,"\n");
            fprintf(ave_plot,"%s\tmid\t",proteinList[pNum]->name);
            for (int i_plot_dat=0; i_plot_dat<i_dat; i_plot_dat++) {
              fprintf(ave_plot,"%1.2f\t",(proteinList[pNum]->numMid[i_plot_dat]/middle_area_total));
            }
            fprintf(ave_plot,"\n");
            fprintf(ave_plot,"%s\trightdown\t",proteinList[pNum]->name);
            for (int i_plot_dat=0; i_plot_dat<i_dat; i_plot_dat++) {
              fprintf(ave_plot,"%1.2f\t",(proteinList[pNum]->numRightDown[i_plot_dat]/right_down_area_total));
            }
            fprintf(ave_plot,"\n");
            fprintf(ave_plot,"\n");
          }
        }
        if (rand_seed == 96){
          for (int pNum=3; pNum<numProteins; pNum++) {
            fprintf(ave_plot,"%s\trightup\t",proteinList[pNum]->name);
            for (int i_plot_dat=0; i_plot_dat<i_dat; i_plot_dat++) {
              fprintf(ave_plot,"%1.2f\t",(proteinList[pNum]->numRightUp[i_plot_dat]/right_up_area_total));
            }
            fprintf(ave_plot,"\n");
            fprintf(ave_plot,"%s\tleftup\t",proteinList[pNum]->name);
            for (int i_plot_dat=0; i_plot_dat<i_dat; i_plot_dat++) {
              fprintf(ave_plot,"%1.2f\t",(proteinList[pNum]->numLeftUp[i_plot_dat]/left_up_area_total));
            }
            fprintf(ave_plot,"\n");
            fprintf(ave_plot,"%s\tleftdown\t",proteinList[pNum]->name);
            for (int i_plot_dat=0; i_plot_dat<i_dat; i_plot_dat++) {
              fprintf(ave_plot,"%1.2f\t",(proteinList[pNum]->numLeftDown[i_plot_dat]/left_down_area_total));
            }
            fprintf(ave_plot,"\n");
            fprintf(ave_plot,"%s\trightdown\t",proteinList[pNum]->name);
            for (int i_plot_dat=0; i_plot_dat<i_dat; i_plot_dat++) {
              fprintf(ave_plot,"%1.2f\t",(proteinList[pNum]->numRightDown[i_plot_dat]/right_down_area_total));
            }
            fprintf(ave_plot,"\n");
            fprintf(ave_plot,"\n");
          }
        }
        fclose(ave_plot);
        delete[] avename;
      }

      for (int pNum=0; pNum<numProteins; pNum++) {
        //time map
        char *timename = print_filename("time-map",proteinList[pNum]->name);
        FILE* time_map = fopen(timename,"w");
        delete[] timename;
        for (int a=0; a<Ny; a++) {
          for (int b=0; b<Nz; b++) {
            fprintf(time_map,"%1.2f\t",(proteinList[pNum]->sum[a*Nz+b])/((double)iter));
          }
          fprintf(time_map,"\n");
        }
        fclose(time_map);
      }
    }

    for (int pNum=0; pNum<numProteins; pNum++) {
      if (i<arrow_iter) {
        //would like to know how long this algorithm takes and change it if should
        if (i%printout_iterations == 0) {
          int* time_maxima_y = new int[arrow_iter];
          int* time_maxima_z = new int[arrow_iter];
          double* time_maxima_value = new double[arrow_iter];
          //printf("We're in the arrow printout loop now!!!\n");
          for (int p=1; p<(i-(.5/time_step)); p++) {
            double max_value = 0;
            int max_k = 0;
            for (int k = int(p-(.5/time_step)); k<(p+(.5/time_step)); k++){
              if (proteinList[pNum]->maxval[k] > max_value) {
                max_value = proteinList[pNum]->maxval[k];
                max_k = k;
              }
            }
            if( max_k == p) {
              time_maxima_y[p] = proteinList[pNum]->ymax[p];
              time_maxima_z[p] = proteinList[pNum]->zmax[p];
              time_maxima_value[p] = proteinList[pNum]->maxval[p];
            }
            else {
              time_maxima_y[p] = 0;
              time_maxima_z[p] = 0;
              time_maxima_value[p] = 0;
            }
          }
          //print to file
          char *arrowname = print_filename("arrow-plot",proteinList[pNum]->name);
          FILE* arrowfile = fopen(arrowname,"w");
          delete[] arrowname;
          for (int p=1; p<(i-1); p++) {
            if ((time_maxima_y[p] != 0) && (time_maxima_z[p] != 0)) {
              fprintf(arrowfile,"%d\t%d\t%g\t%g\n",time_maxima_y[p],time_maxima_z[p],p*time_step,time_maxima_value[p]);
            }
          }
          fclose(arrowfile);
          delete[] time_maxima_y;
          delete[] time_maxima_z;
          delete[] time_maxima_value;
        }
      }
    }
  }
  //end file printing
  //end simulation

  for (int pNum=0; pNum<numProteins; pNum++) {
    delete[] proteinList[pNum]->sum;
    delete[] proteinList[pNum]->maxval;
    delete[] proteinList[pNum]->ymax;
    delete[] proteinList[pNum]->zmax;
    delete[] proteinList[pNum]->numLeft;
    delete[] proteinList[pNum]->numMid;
    delete[] proteinList[pNum]->numRight;
    delete[] proteinList[pNum]->numRightUp;
    delete[] proteinList[pNum]->numRightDown;
    delete[] proteinList[pNum]->numLeftUp;
    delete[] proteinList[pNum]->numLeftDown;
  }

  delete[] JxATP;
  delete[] JyATP;
  delete[] JzATP;
  delete[] JxADP;
  delete[] JyADP;
  delete[] JzADP;
  delete[] JxE;
  delete[] JyE;
  delete[] JzE;

  //printing plot information

  //printing to the project directory so we have a shortlist of what we've done.
  char *fname = new char[1024];
  sprintf(fname,"catalog.txt");
  FILE * catalog;
  int catalog_exists;
  catalog = fopen(fname,"r");
  if (catalog==NULL) {
    catalog_exists=0;
  }
  else {
    catalog_exists=1;
    fclose(catalog);
  }
  if (catalog_exists==1) {
    catalog=fopen(fname,"a+b");
  }
  else {
    catalog=fopen(fname,"w+b");
  }
  if (catalog!=NULL) {
    fprintf(catalog,"%s %1.2f %1.2f %1.2f %1.2f %1.2f", mem_f_shape.c_str(),A,B,C,D,density_factor);
    if (dx==.05) {
      fprintf(catalog," -hires\n");
    }
    else {
      fprintf(catalog,"\n");
    }
    fclose(catalog);
  }
  delete[] fname;
  delete[] debug_flag_str;
  delete[] hires_flag_str;
  delete[] slice_flag_str;

  for (int pNum=0; pNum<numProteins; pNum++) {
    delete[] proteinList[pNum]->name;
    delete proteinList[pNum];
  }
  return 0;
}



void count_and_print_proteins(int iteration, int *s_N_ATP, int *s_N_ADP, int *s_N_E, int *s_ND, int *s_NDEs_N_ATP, double *nATP, double *nADP,
                              double *nE, double *ND, double *NDE, double *NflD, double *NflE, double *mem_A, bool *insideArr) {
  double dV = dx*dx*dx;
  double total_D_n = 0;
  double total_E_n = 0;
  double total_D_N = 0;
  double total_E_N = 0;
  double total_NflD = 0;
  double total_NflE = 0;
  double total_NE = 0;
  double total_NDE = 0;
  double total_NATP = 0;
  double total_NADP = 0;
  double total_ND = 0;
  for (int h=0;h<Nx*Ny*Nz;h++) {
    total_D_n += (nATP[h] + nADP[h])*dx*dx*dx + ND[h] + NDE[h];
    total_E_n += nE[h]*dx*dx*dx + NDE[h];
    total_D_N += s_N_ATP[h] + s_N_ADP[h] + s_ND[h] + s_NDE[h];
    total_E_N += s_N_E[h] + s_NDE[h];
    total_NflD += NflD[h];
    total_NflE += NflE[h];
    total_NE += nE[h]*dx*dx*dx;
    total_NDE += NDE[h];
    total_NATP += nATP[h]*dx*dx*dx;
    total_NADP += nADP[h]*dx*dx*dx;
    total_ND += ND[h];
  }
  double gridpoints_inside=0;
  for (int h=0;h<Nx*Ny*Nz;h++) {
    if (insideArr[h]){
      gridpoints_inside += 1;
    }
  }
  double number_of_membrane_gridpts=0;
  double tot_mem_A=0;
  for (int h=0;h<Nx*Ny*Nz;h++) {
    if (mem_A[h] != 0.0) {
      number_of_membrane_gridpts += 1;
      tot_mem_A += mem_A[h];
    }
  }
  double ave_mem_A = tot_mem_A/number_of_membrane_gridpts;
  printf("\nFor %s simulation:\ntot_mem_A=%g number_of_membrane_gridpts=%g gridpoints_inside=%g ave_mem_A=%g",
         sim_type.c_str(), tot_mem_A, number_of_membrane_gridpts, gridpoints_inside, ave_mem_A);
  // printf("\nAt iteration %d:\ntotal_D_n = %g\ntotal_E_n = %g\ntotal_D_N = %g\ntotal_E_N = %g\ntotal_NflD = %g\ntotal_NflE = %g\n"
  //        ,iteration,total_D_n,total_E_n,total_D_N,total_E_N,total_NflD,total_NflE);
  printf("\ntotal_NE = %g\ntotal_NDE = %g\ntotal_NATP = %g\ntotal_NADP = %g\ntotal_ND = %g\n"
         ,total_NE,total_NDE,total_NATP,total_NADP,total_ND);
  printf("\nAverage for one gridpt:\nnE = %g\nsigmaDE = %g\nnATP = %g\nnADP = %g\nsigmaND = %g\n",
         total_NE/(gridpoints_inside*dV),total_NDE/(number_of_membrane_gridpts*ave_mem_A),total_NATP/(gridpoints_inside*dV),
         total_NADP/(gridpoints_inside*dV),total_ND/(number_of_membrane_gridpts*ave_mem_A));
  printf("\nMore Average for one gridpt:\nNE = %g\nNDE = %g\nNATP = %g\nNADP = %g\nND = %g\n",
         total_NE/gridpoints_inside,total_NDE/number_of_membrane_gridpts,total_NATP/gridpoints_inside,
         total_NADP/gridpoints_inside,total_ND/number_of_membrane_gridpts);
  fflush(stdout);
  return;
}



int get_J(double difD, double *nATP, double *nADP, double *nE,
          double *JxATP, double *JyATP, double *JzATP,
          double *JxADP, double *JyADP, double *JzADP,
          double *JxE, double *JyE, double *JzE){
  for(int xi=0;xi<Nx-1;xi++){
    for(int yi=0;yi<Ny-1;yi++){
      for(int zi=0;zi<Nz-1;zi++){
        JzATP[xi*Ny*Nz+yi*Nz+zi] = -difD*(nATP[xi*Ny*Nz+yi*Nz+zi+1]-nATP[xi*Ny*Nz+yi*Nz+zi])/dx;
        JyATP[xi*Ny*Nz+yi*Nz+zi] = -difD*(nATP[xi*Ny*Nz+(yi+1)*Nz+zi]-nATP[xi*Ny*Nz+yi*Nz+zi])/dx;
        JxATP[xi*Ny*Nz+yi*Nz+zi] = -difD*(nATP[(xi+1)*Ny*Nz+yi*Nz+zi]-nATP[xi*Ny*Nz+yi*Nz+zi])/dx;
        JzADP[xi*Ny*Nz+yi*Nz+zi] = -difD*(nADP[xi*Ny*Nz+yi*Nz+zi+1]-nADP[xi*Ny*Nz+yi*Nz+zi])/dx;
        JyADP[xi*Ny*Nz+yi*Nz+zi] = -difD*(nADP[xi*Ny*Nz+(yi+1)*Nz+zi]-nADP[xi*Ny*Nz+yi*Nz+zi])/dx;
        JxADP[xi*Ny*Nz+yi*Nz+zi] = -difD*(nADP[(xi+1)*Ny*Nz+yi*Nz+zi]-nADP[xi*Ny*Nz+yi*Nz+zi])/dx;
        JzE[xi*Ny*Nz+yi*Nz+zi] = -difD*(nE[xi*Ny*Nz+yi*Nz+zi+1]-nE[xi*Ny*Nz+yi*Nz+zi])/dx;
        JyE[xi*Ny*Nz+yi*Nz+zi] = -difD*(nE[xi*Ny*Nz+(yi+1)*Nz+zi]-nE[xi*Ny*Nz+yi*Nz+zi])/dx;
        JxE[xi*Ny*Nz+yi*Nz+zi] = -difD*(nE[(xi+1)*Ny*Nz+yi*Nz+zi]-nE[xi*Ny*Nz+yi*Nz+zi])/dx;
      }
    }
  }
  return 0;
}



int get_next_density(double *mem_A, bool *insideArr, double *nATP, double *nADP,
                     double *nE, double *ND, double *NDE, double *NflD, double *NflE,
                     double *JxATP, double *JyATP, double *JzATP,
                     double *JxADP, double *JyADP, double *JzADP,
                     double *JxE, double *JyE, double *JzE){
  for(int xi=0;xi<Nx-1;xi++){
    for(int yi=0;yi<Ny-1;yi++){
      for(int zi=0;zi<Nz-1;zi++){
        //for the diffusion terms, we use dn/dt = dA*J/dV thinking in terms of tot #, but dA/dV=1/dx
        if (insideArr[(xi+1)*Ny*Nz+yi*Nz+zi] && insideArr[xi*Ny*Nz+yi*Nz+zi]){
          nADP[(xi+1)*Ny*Nz+yi*Nz+zi] += JxADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nADP[xi*Ny*Nz+yi*Nz+zi] -= JxADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nATP[(xi+1)*Ny*Nz+yi*Nz+zi] += JxATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nATP[xi*Ny*Nz+yi*Nz+zi] -= JxATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nE[(xi+1)*Ny*Nz+yi*Nz+zi] += JxE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nE[xi*Ny*Nz+yi*Nz+zi] -= JxE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
        }
        if (insideArr[xi*Ny*Nz+(yi+1)*Nz+zi] && insideArr[xi*Ny*Nz+yi*Nz+zi]){
          nADP[xi*Ny*Nz+(yi+1)*Nz+zi] += JyADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nADP[xi*Ny*Nz+yi*Nz+zi] -= JyADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nATP[xi*Ny*Nz+(yi+1)*Nz+zi] += JyATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nATP[xi*Ny*Nz+yi*Nz+zi] -= JyATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nE[xi*Ny*Nz+(yi+1)*Nz+zi] += JyE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nE[xi*Ny*Nz+yi*Nz+zi] -= JyE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
        }
        if (insideArr[xi*Ny*Nz+yi*Nz+(zi+1)] && insideArr[xi*Ny*Nz+yi*Nz+zi]){
          nADP[xi*Ny*Nz+yi*Nz+(zi+1)] += JzADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nADP[xi*Ny*Nz+yi*Nz+zi] -= JzADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nATP[xi*Ny*Nz+yi*Nz+(zi+1)] += JzATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nATP[xi*Ny*Nz+yi*Nz+zi] -= JzATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nE[xi*Ny*Nz+yi*Nz+(zi+1)] += JzE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nE[xi*Ny*Nz+yi*Nz+zi] -= JzE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
        }
      }
    }
  }
  double ADP_to_ATP;
  double de_to_ADP_E;
  double ATP_to_d;
  double E_d_to_de;
  for(int xi=0;xi<Nx;xi++){
    for(int yi=0;yi<Ny;yi++){
      for(int zi=0;zi<Nz;zi++){
        ADP_to_ATP = rate_ADP_ATP*nADP[xi*Ny*Nz+yi*Nz+zi]*time_step;
        de_to_ADP_E = rate_de*NDE[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*time_step;
        ATP_to_d = (rate_D + rate_dD*(ND[xi*Ny*Nz+yi*Nz+zi] + NDE[xi*Ny*Nz+yi*Nz+zi])/mem_A[xi*Ny*Nz+yi*Nz+zi])
          *nATP[xi*Ny*Nz+yi*Nz+zi]*time_step;
        E_d_to_de = rate_E*ND[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*nE[xi*Ny*Nz+yi*Nz+zi]*time_step;
        //Jeff!  remember that when you gain cyto density and lose the same amount of wall density,
        //the numbers of proteins gained/lost will be different, and it's the numbers that you want to be the same!!
        //also, keep thinking about the issue below where all additions are divided and then mult by the same mem_A fun
        nADP[xi*Ny*Nz+yi*Nz+zi] -= ADP_to_ATP;
        nATP[xi*Ny*Nz+yi*Nz+zi] += ADP_to_ATP;
        if (mem_A[xi*Ny*Nz+yi*Nz+zi] != 0){
          NDE[xi*Ny*Nz+yi*Nz+zi] += -de_to_ADP_E*mem_A[xi*Ny*Nz+yi*Nz+zi];
          nADP[xi*Ny*Nz+yi*Nz+zi] += de_to_ADP_E*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);
          nE[xi*Ny*Nz+yi*Nz+zi] += de_to_ADP_E*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);

          nATP[xi*Ny*Nz+yi*Nz+zi] -= ATP_to_d*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);
          ND[xi*Ny*Nz+yi*Nz+zi] += ATP_to_d*mem_A[xi*Ny*Nz+yi*Nz+zi];

          nE[xi*Ny*Nz+yi*Nz+zi] -= E_d_to_de*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);
          ND[xi*Ny*Nz+yi*Nz+zi] -= E_d_to_de*mem_A[xi*Ny*Nz+yi*Nz+zi];
          NDE[xi*Ny*Nz+yi*Nz+zi] += E_d_to_de*mem_A[xi*Ny*Nz+yi*Nz+zi];
        }
        NflD[xi*Ny*Nz+yi*Nz+zi] = (nATP[xi*Ny*Nz+yi*Nz+zi] + nADP[xi*Ny*Nz+yi*Nz+zi])*(dx*dx*dx) + ND[xi*Ny*Nz+yi*Nz+zi] + NDE[xi*Ny*Nz+yi*Nz+zi];
        NflE[xi*Ny*Nz+yi*Nz+zi] = nE[xi*Ny*Nz+yi*Nz+zi]*dx*dx*dx + NDE[xi*Ny*Nz+yi*Nz+zi];
      }
    }
  }
  return 0;
}



int set_density(double *nATP, double *nADP, double *nE, double *ND, double *NDE,
                int *s_ND, int *s_NDE, int *s_N_ATP, int *s_N_ADP, int *s_N_E, double *mem_A, bool *insideArr){
  double dV = dx*dx*dx;
  printf("In set_density function, Nx = %d Ny = %d Nz = %d\n",Nx,Ny,Nz);
  int right_most_point_z=0; //left and right most points for z
  int left_most_point_z=Nz;
  int right_most_point_y=0; //"left" and "right" most points for y, in terms of magnitude (right = larger y value)
  int left_most_point_y=Ny;
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (insideArr[i*Ny*Nz+j*Nz+k]){
          if (k>right_most_point_z){
            right_most_point_z = k;
          }
          if(k<left_most_point_z){
            left_most_point_z = k;
          }
          if (j>right_most_point_y){
            right_most_point_y = j;
          }
          if (j<left_most_point_y){
            left_most_point_y = j;
          }
        }
      }
    }
  }
  int density_divider_right = int(right_most_point_z - (right_most_point_z - left_most_point_z)/3);
  //int density_divider_left = int(right_most_point_z - 2*(right_most_point_z - left_most_point_z)/3);
  int gridpoints_low_dens = 0;
  int gridpoints_high_dens = 0;
  int gridpoints_total = 0;
  double wall_area_high = 0;
  double wall_area_low = 0;
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (insideArr[i*Ny*Nz+j*Nz+k]){
          gridpoints_total++;
          if (k>density_divider_right) {
            gridpoints_high_dens++;
            wall_area_high += mem_A[i*Ny*Nz+j*Nz+k];
          }
          else {
            gridpoints_low_dens++;
            wall_area_low += mem_A[i*Ny*Nz+j*Nz+k];
          }
        }
      }
    }
  }
  double density_factor_low_dens = gridpoints_total/(gridpoints_low_dens + density_factor*gridpoints_high_dens);
  double density_factor_high_dens = density_factor*gridpoints_total/(gridpoints_low_dens + density_factor*gridpoints_high_dens);
  //begin setting density at each gridpoint:
  if (sim_type == "exact"){
    for (int i=0;i<Nx;i++){
      for (int j=0;j<Ny;j++){
        for (int k=0;k<Nz;k++){
          if (insideArr[i*Ny*Nz+j*Nz+k]){
            if (k>density_divider_right) {
              nE[i*Ny*Nz+j*Nz+k] = nE_starting_density*density_factor_high_dens;
            }
            else {
              nE[i*Ny*Nz+j*Nz+k] = nE_starting_density*density_factor_low_dens;
            }
            if (k>density_divider_right) {
              nATP[i*Ny*Nz+j*Nz+k] = nATP_starting_density*density_factor_high_dens;
            }
            else {
              nATP[i*Ny*Nz+j*Nz+k] = nATP_starting_density*density_factor_low_dens;
            }
          }
          else {
            nATP[i*Ny*Nz+j*Nz+k] = 0;
            nE[i*Ny*Nz+j*Nz+k] = 0;
          }
          nADP[i*Ny*Nz+j*Nz+k] =0;
          ND[i*Ny*Nz+j*Nz+k] =0;
          NDE[i*Ny*Nz+j*Nz+k] = 0;
          NflD[i*Ny*Nz+j*Nz+k] = ND[i*Ny*Nz+j*Nz+k] + NDE[i*Ny*Nz+j*Nz+k] + (nATP[i*Ny*Nz+j*Nz+k] + nADP[i*Ny*Nz+j*Nz+k])*dV;
          NflE[i*Ny*Nz+j*Nz+k] = NDE[i*Ny*Nz+j*Nz+k] + nE[i*Ny*Nz+j*Nz+k]*dV;
        }
      }
    }
    printf("This is for exact!:\n\n");
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        printf("%g ",NflD[(Nx/2)*Ny*Nz+j*Nz+k]);
      }
      printf("\n");
    }
    fflush(stdout);
  }
  else {
    int E_total = nE_starting_density*gridpoints_total*dV;
    int ATP_total = nATP_starting_density*gridpoints_total*dV;
    printf("E_total = %d ATP_total = %d\n",E_total,ATP_total);
    while (E_total > 0) {
      for (int i=0;i<Nx;i++){
        for (int j=0;j<Ny;j++){
          for (int k=0;k<Nz;k++){
            if (insideArr[i*Ny*Nz+j*Nz+k]){
              if (k>density_divider_right) {
                if ( (double)rand()/(RAND_MAX) < density_factor_high_dens/10000.0) {
                  s_N_E[i*Ny*Nz+j*Nz+k] += 1;
                  E_total -= 1;
                }
              }
              else {
                if ( (double)rand()/(RAND_MAX) < density_factor_low_dens/10000.0) {
                  s_N_E[i*Ny*Nz+j*Nz+k] += 1;
                  E_total -= 1;
                }
              }
            }
            else {
              s_N_E[i*Ny*Nz+j*Nz+k] = 0.0;
            }
          }
        }
      }
    }
    while (ATP_total > 0) {
      for (int i=0;i<Nx;i++){
        for (int j=0;j<Ny;j++){
          for (int k=0;k<Nz;k++){
            if (insideArr[i*Ny*Nz+j*Nz+k]){
              if (k>density_divider_right) {
                if ( (double)rand()/(RAND_MAX) < density_factor_high_dens/10000.0) {
                  s_N_ATP[i*Ny*Nz+j*Nz+k] += 1;
                  ATP_total -= 1;
                }
              }
              else {
                if ( (double)rand()/(RAND_MAX) < density_factor_low_dens/10000.0) {
                  s_N_ATP[i*Ny*Nz+j*Nz+k] += 1;
                  ATP_total -= 1;
                }
              }
            }
            else {
              s_N_ATP[i*Ny*Nz+j*Nz+k] = 0.0;
            }
          }
        }
      }
    }
    for (int i=0;i<Nx*Ny*Nz;i++){
      s_N_ADP[i] = 0;
      s_ND[i] = 0;
      s_NDE[i] = 0;
      NflD[i] = s_ND[i] + s_NDE[i] + s_N_ATP[i] + s_N_ADP[i];
      NflE[i] = s_NDE[i] + s_N_E[i];
    }
    printf("This is for %s simulation!:\n\n",sim_type.c_str());
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        printf("%g ",NflD[(Nx/2)*Ny*Nz+j*Nz+k]);
      }
      printf("\n");
    }
  fflush(stdout);
  }
  return 0;
}
