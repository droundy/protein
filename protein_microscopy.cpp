#include <iostream>
using namespace std;
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "protein.h"
#include "MersenneTwister.h"
#include "weights.h"
#include <cassert>
#include <unistd.h>

const bool stochastic = true;

const double difD = 2.5; // (um)^2 s^- 1
const double difE = 2.5; // (um)^2 s^-1
const double rate_ADP_ATP = 1; // s^-1
const double rate_D = .025; // um s^-1
const double rate_dD = .0015; // (um)^3 s^-1
const double rate_de = .7; // s^-1
const double rate_E = .093; // (um)^3 s^-1

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


double dx; //grid spacing
double tot_time; //total simulation time
double time_step; //simulation time step
int iter; //total # of iterations
int printout_iterations; //iterations between file printout

int box_divider_left;
int box_divider_right;

double x, y, z;
int Nx, Ny, Nz; //number of gridpoints in each direction
int min_xi;
int min_yi;
int min_zi;

bool once_inside = true;

string mem_f_shape; //cell shape argument
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

enum reaction {ADP_to_ATP, DE_to_ADP_E, ATP_to_D, E_D_to_DE, X_ADP_pos, X_ADP_neg, Y_ADP_pos, Y_ADP_neg, Z_ADP_pos,
               Z_ADP_neg, X_ATP_pos, X_ATP_neg, Y_ATP_pos, Y_ATP_neg, Z_ATP_pos, Z_ATP_neg, X_E_pos, X_E_neg,
               Y_E_pos, Y_E_neg, Z_E_pos, Z_E_neg};
const int num_pos_reactions = Z_E_neg+1;
const int d [6][3] =
  {
    {1,0,0},
    {-1,0,0},
    {0,1,0},
    {0,-1,0},
    {0,0,1},
    {0,0,-1},
  };

const int starting_num_guassians=20;
double guass[3*starting_num_guassians]; //stores y,z, and sigma for each guassian when creating random cell wall
const int random_num_guassians=5;
int rand_seed=0; //=14; at this point I have this passed in from the command line as the D argument
double Norm = 15.0; //This is the height of the guassians that make the cell wall

//begin randst shape functions
double rand_dis(double d0,double d_fac,int i) {
  int fac = 1;
  srand(i+rand_seed);
  double x = (rand()%1000);
  x = x/1000.0;
  if ( (rand()%1000) < 500) fac = -1;
  return d0*(1 + fac*(-d_fac*log(1-x)));
}

void randomize_cell_wall(double guass[]){
  //double X = Nx*dx; unused
  double Y = Ny*dx;
  double Z = Nz*dx;
  guass[0]=Y/2.0;
  guass[1]=Z/2.0;
  guass[2] = rand_dis(0.2,.13,2);
  for (int i = 1; random_num_guassians; i++){
    srand(i+rand_seed);
    double sigma = rand_dis(0.2,.13,i);
    guass[3*i+2]=sigma;
    double d = rand_dis(0.4,.13,i);
    double theta = fmod(rand(),2*M_PI);
    double y_change = d*sin(theta);
    double z_change = d*cos(theta);
    guass[i*3] = guass[(i-1)*3]+y_change;
    guass[i*3+1] = guass[(i-1)*3+1]+z_change;
  }
}

double f_2D_TIE_fighter(double y, double z){
  double f = 0;double f1 = 0;double f2 = 0;double f3 = 0;
  double Y = Ny*dx;
  double Z = Nz*dx;
  if (y<2.6){
    f1 = (z-Z/2.0)*(z-Z/2.0)/4 + (y-1.5)*(y-1.5)/0.6 - 1.0;
  } else {f1 = 0.1;}
  if (Y-y < 2.6){
    f2 = (z-Z/2.0)*(z-Z/2.0)/4 + (y-(Y-1.5))*(y-(Y-1.5))/0.6 - 1.0;
  } else {f2 = 0.1;}
  if (abs(z-Z/2.0) < 0.8 && abs(y-Y/2.0) < 2.0){
    f3 = abs(z-Z/2.0)-.5;
  } else {f3 = 0.1;}
  f = f1+f2+f3;
  return f;
}

bool only_print_once = true;
double f_2D_triangle(double y, double z){
  double Y = Ny*dx; double Z = Nz*dx; // total height and width of grid
  double z1 = A/2.0+2*dx; double y1 = A/2.0+2*dx; // lower left corner of triangle
  if (y < y1) {
    return 0.1; // it's too low to be in the triangle
  }

  double z2 = Z-A/2.0-2*dx; double y2 = y1; // lower right corner of triangle
  //Using law of cosines from lengths of sides we get top corner:
  double cos_theta = (B*B+D*D-C*C)/(2*B*D);
  double z3 = A/2.0+2*dx+D*cos_theta; double y3 = Y-A/2.0-2*dx; // top corner of triangle
  if (only_print_once==true){
    printf("z1 = %g y1 = %g\nz2 = %g y2 = %g\nz3 = %g y3 = %g\n",z1,y1,z2,y2,z3,y3);
    printf("cos_theta = %g\nZ = %g Y = %g A = %g",cos_theta,Z,Y,A/2.0);
    only_print_once = false;
  }

  if (z >= z3) {
    double fac = (z-z3)/(z2-z3); //how far along the line the z coordinate is
    double y_line = fac*(y2-y3)+y3; //the y coordinate of the line at the z point
    if (y > y_line) {
      return 0.1; // it's outside the triangle on the right side
    }
  }
  if (z < z3) {
    double fac = (z-z1)/(z3-z1); //how far along the line the z coordinate is
    double y_line = fac*(y3-y1)+y1; //the y coordinate of the line at the z point
    if (y > y_line) {
      return 0.1; // its outside the triangle on the left side
    }
  }
  //double rad = 1.75*(z2-z1)*sqrt(3.0)/6.0;
  //double y_circle = Y/2.0; double z_circle = z1 + sqrt(3)*(y2-y1)/6.0;
  //if ((z < zl1) && (z < zl3) && (z > z1) && ((z-z_circle)*(z-z_circle) + (y-y_circle)*(y-y_circle)) < rad*rad){
  //  return -0.1;
  //} else {
  return -0.1;
}

double f_2D_randst(double y, double z){
  int num_guassians = 0;
  for(int i=0;i<starting_num_guassians;i++){
    if (guass[i*3]!=0.0 || guass[i*3+1]!=0.0 || guass[i*3+2]!=0.0){
      num_guassians++;
    }
  }
  double f=0;
  for (int i = 0; i<num_guassians; i++){
    double arg = ((y-guass[3*i])*(y-guass[3*i]) + (z-guass[3*i+1])*(z-guass[3*i+1]))/(guass[3*i+2]*guass[3*i+2]);
    double fi = -((Norm*exp(-arg))-exp(-1.0));
    f+=fi;
  }
  return f;
}

bool only = true;
double f_2D_stad(double y, double z){
  double f;
  double Y = Ny*dx;
  double Z = Nz*dx;
  double z1 = (Z-C)/2;
  if(only){
    printf("z1 = %g Z = %g C = %g\n",z1,Z,C);
    fflush(stdout);
    only = false;
  }
  double z2 = (C +(Z-C)/2);
  double y1 = Y/2;
  if (z < z1) {
    f = (y-y1)*(y-y1)/(B*B)+(z-z1)*(z-z1)/(D*D)-1;
    if (f > 0.0){
      return 0.1;
    }
  }
  if (z > z2) {
    f = (y-y1)*(y-y1)/(B*B)+(z-z2)*(z-z2)/(D*D)-1;
    if (f > 0.0){
      return 0.1;
    }
  }
  if (abs(y-y1) > B){
    return 0.1;
  }
  return -0.1;
}
//end randst shape functions

//mem_f produces a scalar field on the grid. points where mem_f = 0 are cell wall.
double mem_f(double x, double y, double z) {
  if(mem_f_shape=="randst" || mem_f_shape=="TIE_fighter" || mem_f_shape=="triangle" || mem_f_shape=="stad"){
    double X = Nx*dx;
    double f_2d = 0;
    if(mem_f_shape=="randst") f_2d = f_2D_randst(y,z);
    else if(mem_f_shape=="TIE_fighter") f_2d = f_2D_TIE_fighter(y,z);
    else if(mem_f_shape=="triangle") f_2d = f_2D_triangle(y,z);
    else if(mem_f_shape=="stad") f_2d = f_2D_stad(y,z);
    else {
      printf("somethings wrong with the shape argument!!!");
      exit(1);
    }
    if (f_2d<=0) {
      // if (abs(2*(x-(X/2.0))/A - 0.99999) < .1){
      //   printf("f_2d = %g and X/2.0 = %g and 2*(x-(X/2.0))/A - 0.99999 = %g and x = %g\n",f_2d,X/2.0,2*(x-(X/2.0))/A - 0.99999,x);
      //   fflush(stdout);
      // }
      return abs(2*(x-(X/2.0))/A) - 0.99999;//1.00001;
    }
    double closest_y0 = -100.0;
    double closest_z0 = -100.0;
    //bool there_is_closest_point=0; unused
    for (double y0 = y-(A/2.0+2.0*dx); y0<y+(A/2.0+2.0*dx); y0+=dx/2.0) {
      for (double z0 = z-(A/2.0+2.0*dx); z0<z+(A/2.0+2.0*dx); z0+=dx/2.0) {
        if ( (y-y0)*(y-y0)+(z-z0)*(z-z0) < (y-closest_y0)*(y-closest_y0)+(z-closest_z0)*(z-closest_z0) ) {
          double f0 = 0;
          if(mem_f_shape=="randst") f0 = f_2D_randst(y0,z0);
          if(mem_f_shape=="TIE_fighter") f0 = f_2D_TIE_fighter(y0,z0);
          if(mem_f_shape=="triangle") f0 = f_2D_triangle(y0,z0);
          if(mem_f_shape=="stad") f0 = f_2D_stad(y0,z0);
          if (f0 <= 0) {
            closest_y0 = y0;
            closest_z0 = z0;
          }
        }
      }
    }
    double dis = sqrt((y-closest_y0)*(y-closest_y0) + (z-closest_z0)*(z-closest_z0) + (x-X/2.0)*(x-X/2.0));
    if (dis < A) {
      return 2.0*(dis/A-.49999);
    } else {
      return 1.0;
    }
  }

  if (mem_f_shape=="p"){
    //A = length, B = radius of endcap and cylinder
    double f;
    double X = Nx*dx;
    double Y = Ny*dx;
    double Z = Nz*dx;
    double z1 = (Z-A)/2.0;
    double z2 = (A+(Z-A)/2.0);
    double x1 = X/2.0;
    double y1 = Y/2.0;
    f = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1))-B;
    if (z < z1) {
      f = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1))-B;
    }
    if (z > z2) {
      f = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z2)*(z-z2))-B;
    }
    return f;
  }
  if (mem_f_shape=="sp"){
    // A = radius
    double X = Nx*dx;
    double Y = Ny*dx;
    double Z = Nz*dx;
    double x1 = X/2.0;
    double y1 = Y/2.0;
    double z1 = Z/2.0;
    double f;
    f = sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1)) - A;
    return f;
  }
  if (mem_f_shape=="e"){
    //B = x axis radius radius, C = y axis radius radius, A = z axis radius radius
    double X = Nx*dx;
    double Y = Ny*dx;
    double Z = Nz*dx;
    double x1 = X/2;
    double y1 = Y/2;
    double z1 = Z/2;
    double f = sqrt( (x-x1)*(x-x1)/(B*B) + (y-y1)*(y-y1)/(C*C)+ (z-z1)*(z-z1)/(A*A) ) - 1;
    return f;
  }

   else {
     double f = 1;
     return f;
   }
 }


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
   sprintf(filename,"data/shape-%s/%s%s%s%s-%s-%s-%1.2f-%1.2f-%1.2f-%1.2f-%1.2f.dat",mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,plotname,proteinname,mem_f_shape.c_str(),A,B,C,D,density_factor);
   return filename;
 }

void sym_check (double * mem_A){
  char *out_filename = new char[1024];
  sprintf(out_filename,"testingsym.txt");
  FILE *out_file = fopen((const char *)out_filename, "w");
  delete[] out_filename;
  fprintf(out_file, "does this show up?\n");

  FILE * out_files[10];
  for (int count=0;count<10;count++){
    char* filename = new char[1024];
    sprintf(filename, "testingsym.txt%d",count);
    out_files[count] = fopen((const char *)filename,"w");
    delete[] filename;
  }
  //fprintf(out_file,"Nx = %d Ny = %d Nz = %d\n",Nx,Ny,Nz);
  for (int count=0;count<10;count++){
    int xi = count+10;
    for(int zi=0;zi<Nz;zi++){
      for(int yi=0;yi<Ny;yi++){
        fprintf(out_files[count], "%g\t",mem_A[xi*Ny*Nz+yi*Nz+zi]);
      }
    fprintf(out_files[count],"\n");
    }
  }
  fprintf(out_file,"Done with mem_A print outs\n");
  for(int zi=0;zi<Nz;zi++){
    for(int yi=-Ny/2;yi<Ny/2;yi++){
      //for(int yi=0;yi<Ny;yi++){
      int xi = 7;
      //int xa = Nx/2 + xi;
      int ya = Ny/2 + yi;
      //int xb = Nx/2 - xi;
      int yb = Ny/2 - yi;
      if (abs(mem_A[xi*Ny*Nz+ya*Nz+zi] - mem_A[xi*Ny*Nz+yb*Nz+zi]) > .000001) {
        fprintf(out_file,"xi = %d ya = %d zi = %d mem(xa) = %g\n",xi,yi,zi,mem_A[xi*Ny*Nz+ya*Nz+zi]);
        //fprintf(out_file,"\n");
      }
      //}
      //fprintf(out_file, "%d\t",int(inside(xi,yi,Nz/2)));
    }
    fprintf(out_file,"\n");
  }
  //fprintf(out_file,"Done with mem_A symmetry comparisons\n");
  // for(int zi=0;zi<Nz;zi++){
  //   for(int xi=0;xi<Nx/2;xi++){
  //     for(int yi=0;yi<Ny;yi++){
  //       int xa = Nx/2 + xi;
  //       //int ya = Ny/2 + yi;
  //       int xb = Nx/2 - xi;
  //       //int yb = Ny/2 - yi;
  //       if (mem_f(xa*dx,yi*dx,zi*dx) - mem_f(xb*dx,yi*dx,zi*dx) > .000001) {
  //         fprintf(out_file,"xa = %g ya = %g z = %g mem_A = %g\n",xa*dx,yi*dx,zi*dx,mem_f(xa*dx,yi*dx,zi*dx));//[xi*Ny*Nz+ya*Nz+zi]);
  //         fprintf(out_file,"xa = %g ya = %g z = %g mem_A = %g\n",xb*dx,yi*dx,zi*dx,mem_f(xb*dx,yi*dx,zi*dx));//[xi*Ny*Nz+ya*Nz+zi]);
  //       }
  //     }
  //     //fprintf(out_file, "%d\t",int(inside(xi,yi,Nz/2)));
  //   }
  //   fprintf(out_file,"\n");
  // }
  //fprintf(out_file,"Starting mem_f check\n!!!");
  // for (int i=0;i<100000000;i++){
  //   double xa = Nx*dx*(double(rand())/RAND_MAX);
  //   double xb = Nx*dx -xa;
  //   double y = Ny*dx*(double(rand())/RAND_MAX);
  //   double z = Nz*dx*(double(rand())/RAND_MAX);
  // for(double z=0;z<Nz*dx;z += dx/4.0){
  //   for(double xa=0;xa<Nx*dx/2.0;xa += dx/4.0){
  //     for(double y=0;y<Ny*dx;y += dx/4.0){
  //       double xb = Nx*dx - xa;
  //       if (abs(mem_f(xa,y,z) - mem_f(xb,y,z)) > .00000001) {
  //         fprintf(out_file,"Nx*dx = %g\n",Nx*dx);
  //         fprintf(out_file,"xa = %g yi = %g z = %g mem_A = %g\n",xa,y,z,mem_f(xa,y,z));//[xi*Ny*Nz+ya*Nz+zi]);
  //         fprintf(out_file,"xb = %g yi = %g z = %g mem_A = %g\n",xb,y,z,mem_f(xb,y,z));//[xi*Ny*Nz+ya*Nz+zi]);
  //       }
  //       //fprintf(out_file,"No! ");
  //     }
  //   }
  // }
  // fprintf(out_file,"Done with mem_f symmetry comparisons\n");
  for (int i=0;i<10;i++){
    fclose(out_files[i]);
  }
  fclose(out_file);
  // printf("Done with sym check!!!!!");
  fflush(stdout);
}


bool only_once = true;
string triangle_section (double y, double z) {
  //needs editing!!!!!
  //double Y = Ny*dx; double Z = Nz*dx; // total height and width of grid
  double z1 = A/2.0+2*dx; double y1 = A/2.0+2*dx; // lower left corner of triangle
  double z2 = B+A/2.0+2*dx; double y2 = y1; // lower right corner of triangle
  //Using law of cosines from lengths of sides we get top corner:
  double theta = acos((B*B+D*D-C*C)/(2*B*D));
  double z3 = A/2.0+2*dx+D*cos(theta); double y3 = y1+D*sin(theta); // top corner of triangle
  if (only_once == true) {
    printf("z1 = %g y1 = %g\nz2 = %g y2 = %g\nz3 %g y3 = %g\n",z1,y1,z2,y2,z3,y3);
  }

  //get bisecting points and lines:
  double y_21 = (y1 + y2)/2.0; double z_21 = (z2 + z1)/2.0;
  double y_32 = (y3 + y2)/2.0; double z_32 = (z3 + z2)/2.0;
  double y_13 = (y1 + y3)/2.0; double z_13 = (z1 + z3)/2.0;
  double slope1 = (y_32-y1)/(z_32-z1); // from left corner to right line
  double slope2 = (y2-y_13)/(z2-z_13); //from right corner to left line
  double slope3 = (y_21-y3)/(z_21-z3); //from top corner to bottom
  //running into nan issues when z3 is same as z_21, so I brute force a large negative slope:
  if (abs(z_21-z3) < .000001){
    slope3 = 1000000*(y_21-y3);
  }
  if (abs(z2-z_13) < .000001){
    slope2 = 1000000*(y2-y_13);
  }
  if (only_once==true){
    printf("slope1 = %g slope2 = %g slope3 = %g\n",slope1,slope2,slope3);
  }
  //find centroid, which is where all three lines above intersect:
  double z_cen = (y3 - y1 + slope1*z1 - slope3*z3)/(slope1 -slope3);
  double y_cen = slope1*(z_cen-z1) + y1;
  if (only_once ==true){
    printf("z_cen = %g y_cen = %g\n",z_cen,y_cen);
    only_once = false;
  }
  //The density will start higher in the Right section, although there won't be a
  //lot of symmetry (will also have some in the middle section).
  if (z > z_cen) {
    if (y > slope1*(z-z_cen)+y_cen) {
      return "Mid";
    }
  } else {
    if (y > slope2*(z-z_13)+y_13) {
      return "Mid";
    }
  }
  if (y > slope3*(z-z_cen)+y_cen) {
    return "Right";
  }
  return "Left";
}

void trim_grid(double **pointer_to_mem_A, double *first_mem_A){//, int min_xi, int min_yi, int min_zi){
  int max_xi = 0;
  min_xi = Nx;
  int max_yi = 0;
  min_yi = Ny;
  int max_zi = 0;
  min_zi = Nz;
  for(int xi=0;xi<Nx;xi++){
    for(int yi=0;yi<Ny;yi++){
      for(int zi=0;zi<Nz;zi++){
        if (first_mem_A[xi*Ny*Nz+yi*Nz+zi] != 0) {
          if (xi > max_xi){
            max_xi = xi;
          }
          if (yi > max_yi){
            max_yi = yi;
          }
          if (zi > max_zi){
            max_zi = zi;
          }
          if (xi < min_xi){
            min_xi = xi;
          }
          if (yi < min_yi){
            min_yi = yi;
          }
          if (zi < min_zi){
            min_zi = zi;
          }
        }
      }
    }
  }
  int new_Nx = max_xi-min_xi+3;
  int new_Ny = max_yi-min_yi+3;
  int new_Nz = max_zi-min_zi+3;
  *pointer_to_mem_A = new double[new_Nx*new_Ny*new_Nz];
  for(int xi=0;xi<new_Nx;xi++){
    for(int yi=0;yi<new_Ny;yi++){
      for(int zi=0;zi<new_Nz;zi++){
        (*pointer_to_mem_A)[xi*new_Ny*new_Nz+yi*new_Nz+zi] = first_mem_A[(xi+min_xi-1)*Ny*Nz+(yi+min_yi-1)*Nz+(zi+min_zi-1)];
      }
    }
  }
  Nx = new_Nx;
  Ny = new_Ny;
  Nz = new_Nz;
}

void test_the_amount_of_area(double *first_mem_A, string mem_f_shape){
  printf("\n\n\nDid this print ok %s???\n",(char *)mem_f_shape.c_str());
  double A_ours = 0;
  for (int i =0;i<Nx*Ny*Nz;i++){
    A_ours += first_mem_A[i];
  }
  double A_ideal=0;
  if (mem_f_shape=="stad"){
    printf("This in here eh?\n");
    double d_phi = 0.0001;
    double Integral = 0;
    for (double phi = 0; phi<M_PI; phi += d_phi){
      double denominator = sqrt( (B*B*cos(phi)*cos(phi)) + (D*D*sin(phi)*sin(phi)) );
      Integral += A*D*B*M_PI/2/denominator*d_phi;
    }
    A_ideal = 4*C*B + 2*M_PI*B*D + C*A*M_PI + M_PI*A*A/2 + 2*Integral;
  }
  else if (mem_f_shape=="p"){
    A_ideal = 2*M_PI*B*A + 4*M_PI*B*B;
  }
  else if (mem_f_shape=="sp"){
    A_ideal = 4*M_PI*A*A;
  }
  double percent_off = 100.0*abs(A_ideal-A_ours)/A_ideal;
  if (A_ideal > A_ours) {
    printf("\nOur mem_A area is too small and is %g%% off\n\n",percent_off);
  }
  else {
    printf("\nOur mem_A area is too great and is %g%% off\n\n",percent_off);
  }
}



struct stoch_params {
  int xi;
  int yi;
  int zi;
  int reaction;
};


stoch_params index_to_perameters(int index) {
  stoch_params p;
  //should I worry about these floors when the integer is right on?
  p.reaction = floor( index/(Nx*Ny*Nz) );
  p.xi = floor( (index - p.reaction*Nx*Ny*Nz) / (Ny*Nz) );
  p.yi = floor( (index - p.reaction*Nx*Ny*Nz - p.xi*Ny*Nz) / Nz);
  p.zi = index - p.reaction*Nx*Ny*Nz - p.xi*Ny*Nz - p.yi*Nz;
  return p;
}



int main (int argc, char *argv[]) {
  //command line parameters
  mem_f_shape = argv[1];
  A = atof(argv[2]);
  B = atof(argv[3]);
  C = atof(argv[4]);
  D = atof(argv[5]);
  density_factor = atof(argv[6]);
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
  double guass94[] = {2.6,3.2,.3,2.28,2.75,.25,3.0,3.9,.6,3.1,3.6,.4,3.3,3.9,.4,3.5,4.7,.5,2.9,5.6,.5,3.1,5.2,.4,2.4,5.8,.2,3.6,5.1,.3};
  double guass95[] = {2.2,2.4,.3,2.5,3.2,.6,2.7,3.5,.4,2.9,3.5,.4,3.5,4.2,.5,3.8,4.1,.8,3.1,4.6,.6,3.15,4.3,.5};
  double guass96[] = {1.3,1.3,.7,2.1,2,.7,3,2,.7,3.9,2,.7,4.7,1.3,.7,4,2.1,.7,4,3,.7,4,3.9,.7,4.7,4.7,.7,3.9,4,.7,3,4,.7,2.3,4,.7,1.3,4.7,.7,2.1,3.9,.7,3,3.9,.7,2.1,3.9,.7};
  double guass97[] = {1.4,3,.4,1.8,3,.4,2.2,3,.4,2.6,3,.4,3,3,.4,3.4,3,.4,3.8,3,.4,4.2,3,.4,4.6,3,.4,5,3,.4,5.4,3,.4,3.4,2.4,.6};
  double guass98[] = {2.0,2.0,.3,3,3,.6,4.2,3.4,.3,4.6,4.6,.6,3.4,5.6,.6};
  double guass99[] = {2.0,2.2,.5,3,3,.50,4.0,3.6,.50,3,4.2,.50,2.0,5,.5};

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
  printf("\nsize_modifier 94 is = %g",size_modifier_94);
  printf("\nsize_modifier 95 is = %g\n",size_modifier_95);
  printf("\nsize_modifier 99 is = %g\n",size_modifier_99);
  fflush(stdout);

  bzero(guass,int(sizeof(guass)/sizeof(double)));
  if (rand_seed == 94) {
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

  if (mem_f_shape == "p" || mem_f_shape == "stad" || mem_f_shape == "sp"){
    test_the_amount_of_area(first_mem_A,mem_f_shape);
  }

  //Trimming the grid
  double *mem_A;
  if (mem_f_shape!="stad"){
    double **pointer_to_mem_A;
    pointer_to_mem_A = &mem_A;
    printf("\nBefore Trim Nx = %d Ny = %d Nz = %d\n",Nx,Ny,Nz);
    trim_grid(pointer_to_mem_A,first_mem_A);//,min_xi,min_yi,min_zi);
    printf("After Trim Nx = %d Ny = %d Nz = %d\n",Nx,Ny,Nz);
    fflush(stdout);
  }
  else {
    mem_A = new double[Nx*Ny*Nz];
    for (int i=0;i<Nx*Ny*Nz;i++){
      mem_A[i] = first_mem_A[i];
    }
  }

  printf("\nAfter trim, min_xi = %d min_yi = %d min_zi = %d\n",min_xi,min_yi,min_zi);
  fflush(stdout);

  bool *insideArr = new bool[Nx*Ny*Nz];
  set_insideArr(insideArr);
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


  set_density(nATP, nADP, nE, ND, NDE, s_ND, s_NDE, s_N_ATP, s_N_ADP, s_N_E, mem_A);

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
    if (rand_seed == 94) {
      vert_div = size_modifier_94*(3.6/dx)-min_zi+1;
      vert_div_two = size_modifier_94*(5.0/dx)-min_zi+1;
    }
    if (rand_seed == 95) {
      vert_div = size_modifier_95*(3.3/dx)-min_zi+1;
      vert_div_two = size_modifier_95*(4.4/dx)-min_zi+1;
    }
    if (rand_seed == 96) {
      vert_div = size_modifier_96*(2.8/dx)-min_zi+1;
      hor_div = size_modifier_96*(3.0/dx)-min_yi+1;
    }
    if (rand_seed == 97) {
      vert_div = size_modifier_97*(2.4/dx)-min_zi+1;
      hor_div = size_modifier_97*(2.45/dx)-min_yi+1;
      hor_div_two = size_modifier_97*(4.6/dx)-min_yi+1;
    }
    if (rand_seed == 98) {
      vert_div = size_modifier_98*(3.0/dx)-min_zi+1;
      vert_div_two = size_modifier_98*(4.8/dx)-min_zi+1;
    }
    if (rand_seed == 99) {
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
        if (rand_seed == 99 || rand_seed == 98 || rand_seed == 95 || rand_seed == 94) {
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

  weights ws(num_pos_reactions*Nx*Ny*Nz);
  initialize_densities_and_weighting(ws, s_N_ATP, s_N_ADP, s_N_E, s_ND, s_NDE, mem_A);

  ///////////////////////////////////////
  //I need to look into whether I rally need to use s_ND and not just the same ND?
  ////////////////////////////////////

  //starting simuation
  double spill_over_time = 0;
  for (int i=0;i<iter;i++){
    if (stochastic) {
      double elapsed_time = spill_over_time;
      while (elapsed_time < time_step) {
        int index = ws.lookup( (double)rand()/(RAND_MAX) );//random number is from 0 to 1
        stoch_params p = index_to_perameters(index);
        if (p.reaction <= E_D_to_DE) {
          update_densities_and_weighting_for_reaction(p, ws, s_N_ATP, s_N_ADP, s_N_E, s_ND, s_NDE, mem_A);
        }
        else {
          update_densities_and_weighting_for_diffusion(p, ws, s_N_ATP, s_N_ADP, s_N_E, s_ND, s_NDE, mem_A);
        }
        elapsed_time += log( (double)rand()/(RAND_MAX) ) / ws.get_total();
      }
      spill_over_time = elapsed_time - time_step;
      for (int i=0;i<Nx*Ny*Nz;i++) {
        nATP[i] = double(s_N_ATP[i])/dx*dx*dx;
        nADP[i] = double(s_N_ADP[i])/dx*dx*dx;
        nE[i] = double(s_N_E[i])/dx*dx*dx;
        ND[i] = double(s_ND[i]);
        NDE[i] = double(s_NDE[i]);
      }
    }
    else {
      get_J(difD, nATP, nADP, nE, JxATP, JyATP,
            JzATP, JxADP, JyADP, JzADP, JxE, JyE, JzE);
      get_next_density(mem_A, insideArr, nATP, nADP, nE, ND, NDE, NflD, NflE, JxATP, JyATP, JzATP,
                       JxADP, JyADP, JzADP, JxE, JyE, JzE);
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

      // printf("Nine Got Here!\n");
      // fflush(stdout);

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
                if (rand_seed == 99 || rand_seed == 98 || rand_seed == 95 || rand_seed == 94) {
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
      sprintf(outfilenameATP, "data/shape-%s/%s%s%snATP-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor,k);
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
      sprintf(outfilenameE, "data/shape-%s/%s%s%snE-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor,k);
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
      sprintf(outfilenameADP, "data/shape-%s/%s%s%snADP-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor,k);
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
      sprintf(outfilenameD, "data/shape-%s/%s%s%sND-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor,k);
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
      sprintf(outfilenameDE, "data/shape-%s/%s%s%sNDE-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor,k);
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

      // printf("Sixteen Got Here!\n");
      // fflush(stdout);

      //begin NflE printing
      char *outfilenameflE = new char[1000];
      sprintf(outfilenameflE, "data/shape-%s/%s%s%sNflE-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor,k);
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
      sprintf(outfilenameflD, "data/shape-%s/%s%s%sNflD-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", mem_f_shape.c_str(),debug_flag_str,hires_flag_str,slice_flag_str,mem_f_shape.c_str(),A,B,C,D,density_factor,k);
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
      FILE* box_plot = fopen(boxname,"w");
      delete[] boxname;
      for (int pNum=0; pNum<numProteins; pNum++) {

        if (mem_f_shape == "p" || rand_seed == 99 || rand_seed == 98 || rand_seed == 95 || rand_seed == 94 ||
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

        if (mem_f_shape == "p" || mem_f_shape == "triangle" || rand_seed == 98 || rand_seed == 95 || rand_seed == 94 ||
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
          clock_t start_time = clock();
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
          clock_t time_took = clock()-start_time;
          printf("\n\nThe arrow algorithm at iteration time %d took %f seconds to complete\n\n\n",int(i/time_step),double(time_took/CLOCKS_PER_SEC));
        }
      }
    }
  }
  //end file printing
  //end simulation

  printf("Twelve Got Here!\n");
  fflush(stdout);


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

  printf("Thirteen Got Here!\n");
  fflush(stdout);

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
  printf("Fourteen Got Here!\n");
  fflush(stdout);

  return 0;
}


void initialize_densities_and_weighting(weights ws, int *s_N_ATP, int *s_N_ADP, int *s_N_E, int *s_ND, int *s_NDE, double *mem_A){
  /* reminder of the order of the reaction enums and what d is:
     enum reaction {ADP_to_ATP, DE_to_ADP_E, ATP_to_D, E_D_to_DE,
                X_ADP_pos, X_ADP_neg, Y_ADP_pos, Y_ADP_neg, Z_ADP_pos, Z_ADP_neg,
                X_ATP_pos, X_ATP_neg, Y_ATP_pos, Y_ATP_neg, Z_ATP_pos, Z_ATP_neg,
                X_E_pos, X_E_neg, Y_E_pos, Y_E_neg, Z_E_pos, Z_E_neg};
     const int d [6][3] = { {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1} };
     Also jz = -difD*(n[z+1]-n[z])/dx and jz is #/(dA*dt) and I want dN/dt units for our probabilities
     so dN/dt = jz*dA = difD*(n[z+1]-n[z])*dA/dx = difD*(N[z+1]-N[z])/(dA)
  */
  double dV = dx*dx*dx;
  double dA = dx*dx;
  for(int xi=0;xi<Nx;xi++){
    for(int yi=0;yi<Ny;yi++){
      for(int zi=0;zi<Nz;zi++){
        ws.update(rate_ADP_ATP*s_N_ADP[xi*Ny*Nz+yi*Nz+zi], ADP_to_ATP*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
        if (mem_A[xi*Ny*Nz+yi*Nz+zi] != 0.0) {
          ws.update((rate_D + rate_dD*(ND[xi*Ny*Nz+yi*Nz+zi] + NDE[xi*Ny*Nz+yi*Nz+zi])/mem_A[xi*Ny*Nz+yi*Nz+zi])
                    *s_N_ATP[xi*Ny*Nz+yi*Nz+zi], ATP_to_D*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
          ws.update(rate_E*ND[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*s_N_E[xi*Ny*Nz+yi*Nz+zi], E_D_to_DE*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
          ws.update(rate_de*NDE[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*dV, DE_to_ADP_E*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
        }
        else {
          ws.update( 0.0, ATP_to_D*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
          ws.update( 0.0, E_D_to_DE*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
          ws.update( 0.0, DE_to_ADP_E*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
        }
        for (int i=0;i<6;i++){
          if ( s_N_ATP[xi*Ny*Nz+yi*Nz+zi] - s_N_ATP[(xi+d[i][0])*Ny*Nz + (yi+d[i][1])*Nz + (zi+d[i][2])] >= 0) {
            ws.update( difD*( s_N_ATP[xi*Ny*Nz+yi*Nz+zi] - s_N_ATP[ (xi+d[i][0])*Ny*Nz + (yi+d[i][1])*Nz + (zi+d[i][2]) ] ) / dA,
                       (X_ATP_pos + i)*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
          }
          else {
            ws.update( 0.0, (X_ATP_pos + i)*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
          }
          if (s_N_ADP[xi*Ny*Nz+yi*Nz+zi] - s_N_ADP[(xi+d[i][0])*Ny*Nz + (yi+d[i][1])*Nz + (zi+d[i][2])] >= 0) {
            ws.update( difD*( s_N_ADP[xi*Ny*Nz+yi*Nz+zi] - s_N_ADP[(xi+d[i][0])*Ny*Nz + (yi+d[i][1])*Nz + (zi+d[i][2])] ) / dA,
                       (X_ADP_pos + i)*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
          }
          else {
            ws.update( 0.0, (X_ADP_pos + i)*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
          }
          if ( s_N_E[xi*Ny*Nz+yi*Nz+zi] - s_N_E[(xi+d[i][0])*Ny*Nz + (yi+d[i][1])*Nz + (zi+d[i][2])] >= 0) {
            ws.update( difD*( s_N_E[xi*Ny*Nz+yi*Nz+zi] - s_N_E[ (xi+d[i][0])*Ny*Nz + (yi+d[i][1])*Nz + (zi+d[i][2]) ] ) / dA,
                       (X_E_pos + i)*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
          }
          else {
            ws.update( 0.0, (X_E_pos + i)*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
          }
        }
      }
    }
  }
  printf("We have allocated and initialized the probability weighting memory\n");
}



void update_densities_and_weighting_for_reaction(stoch_params p, weights ws, int *s_N_ATP, int *s_N_ADP, int *s_N_E, int *s_ND, int *s_NDE, double *mem_A){
  /* reminder of the order of the reaction enums and what d is:
     enum reaction {ADP_to_ATP, DE_to_ADP_E, ATP_to_D, E_D_to_DE,
                X_ADP_pos, X_ADP_neg, Y_ADP_pos, Y_ADP_neg, Z_ADP_pos, Z_ADP_neg,
                X_ATP_pos, X_ATP_neg, Y_ATP_pos, Y_ATP_neg, Z_ATP_pos, Z_ATP_neg,
                X_E_pos, X_E_neg, Y_E_pos, Y_E_neg, Z_E_pos, Z_E_neg};
     const int d [6][3] = { {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1} };
     Also jz = -difD*(n[z+1]-n[z])/dx and jz is #/(dA*dt) and I want dN/dt units for our probabilities
     so dN/dt = difD*(n[z+1]-n[z])*dA/dx = difD*(N[z+1]-N[z])/(dA)
  */
  double dV = dx*dx*dx;
  double dA = dx*dx;
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
  }
  if (p.reaction == ADP_to_ATP || p.reaction == DE_to_ADP_E) { //s_N_ADP has changed, need to change probs effected by this
    ws.update(rate_ADP_ATP*s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi], ADP_to_ATP*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    //This loop goes around the changed lattice point and updates the diffusion propabilities for going to each adjacent lattice point
    for (int i=0;i<6;i++){
      if (s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ADP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
        ws.update( difD*( s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ADP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] ) / dA,
                   (X_ADP_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
      else {
        ws.update( 0.0, (X_ADP_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
    }
  }
  if (p.reaction == ADP_to_ATP || p.reaction == ATP_to_D) { //s_N_ATP has changed
    if (mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] != 0.0) {
      ws.update((rate_D + rate_dD*(ND[p.xi*Ny*Nz+p.yi*Nz+p.zi] + NDE[p.xi*Ny*Nz+p.yi*Nz+p.zi])/mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi])
                *s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi], ATP_to_D*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
    for (int i=0;i<6;i++){
      if ( s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ATP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
        ws.update( difD*( s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ATP[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ] ) / dA,
                   (X_ATP_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
      else {
        ws.update( 0.0, (X_ATP_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
    }
  }
  if (p.reaction == DE_to_ADP_E || p.reaction == E_D_to_DE) { //s_N_E has changed
    if (mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] != 0.0) {
      ws.update(rate_E*ND[p.xi*Ny*Nz+p.yi*Nz+p.zi]/mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi]*s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi], E_D_to_DE*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
    for (int i=0;i<6;i++){
      if ( s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_E[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
        ws.update( difD*( s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_E[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ] ) / dA,
                   (X_E_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
      else {
        ws.update( 0.0, (X_E_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
    }
  }
  if (p.reaction == ATP_to_D || p.reaction == E_D_to_DE) { //s_ND has changed
    if (mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] != 0.0) {
      ws.update((rate_D + rate_dD*(ND[p.xi*Ny*Nz+p.yi*Nz+p.zi] + NDE[p.xi*Ny*Nz+p.yi*Nz+p.zi])/mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi])
                *s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi], ATP_to_D*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      ws.update(rate_E*ND[p.xi*Ny*Nz+p.yi*Nz+p.zi]/mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi]*s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi], E_D_to_DE*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
  }
  if (p.reaction == DE_to_ADP_E || p.reaction == E_D_to_DE) { //s_NDE has changed
    if (mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] != 0.0) {
      ws.update(rate_de*NDE[p.xi*Ny*Nz+p.yi*Nz+p.zi]/mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi]*dV, DE_to_ADP_E*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      ws.update((rate_D + rate_dD*(ND[p.xi*Ny*Nz+p.yi*Nz+p.zi] + NDE[p.xi*Ny*Nz+p.yi*Nz+p.zi])/mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi])
                *s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi], ATP_to_D*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
  }
  return;
}




void update_densities_and_weighting_for_diffusion(stoch_params p, weights ws, int *s_N_ATP, int *s_N_ADP, int *s_N_E, int *s_ND, int *s_NDE, double *mem_A){
  /* reminder of the order of the reaction enums and what d is:
     enum reaction {ADP_to_ATP, DE_to_ADP_E, ATP_to_D, E_D_to_DE,
                X_ADP_pos, X_ADP_neg, Y_ADP_pos, Y_ADP_neg, Z_ADP_pos, Z_ADP_neg,
                X_ATP_pos, X_ATP_neg, Y_ATP_pos, Y_ATP_neg, Z_ATP_pos, Z_ATP_neg,
                X_E_pos, X_E_neg, Y_E_pos, Y_E_neg, Z_E_pos, Z_E_neg};
     const int d [6][3] = { {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1} };
     Also jz = -difD*(n[z+1]-n[z])/dx and jz is #/(dA*dt) and I want dN/dt units for our probabilities
     so dN/dt = difD*(n[z+1]-n[z])*dx  */
  double dA = dx*dx;
  if (p.reaction >= X_ADP_pos && p.reaction <= Z_ADP_neg) { //s_N_ADP has diffused
    s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] -= 1;
    ws.update(rate_ADP_ATP*s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi], ADP_to_ATP*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    //This loop goes around the changed lattice point and updates the diffusion propabilities for going to each adjacent lattice point
    for (int i=0;i<6;i++){
      if ( s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ADP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
        ws.update( difD*( s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ADP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] ) / dA,
                   (X_ADP_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
      else {
        ws.update( 0.0, (X_ADP_pos + 1)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
    }
    switch (p.reaction) {
      case X_ADP_pos:
        p.xi += 1; break;
      case Y_ADP_pos:
        p.yi += 1; break;
      case Z_ADP_pos:
        p.zi += 1; break;
      case X_ADP_neg:
        p.xi -= 1; break;
      case Y_ADP_neg:
        p.yi -= 1; break;
      case Z_ADP_neg:
        p.zi -= 1; break;
    }
    s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] += 1;
    ws.update(rate_ADP_ATP*s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi]*dx, ADP_to_ATP*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    for (int i=0;i<6;i++){
      if ( s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ADP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
        ws.update( difD*( s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ADP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] ) / dA,
                   (X_ADP_pos + 1)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
      else {
        ws.update( 0.0, (X_ADP_pos + 1)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
    }
    return;
  }
  else if (p.reaction >= X_ATP_pos && p.reaction <= Z_ATP_neg) { //s_N_ATP has diffused
    s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] -= 1;
    if (mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] != 0.0) {
      ws.update((rate_D + rate_dD*(ND[p.xi*Ny*Nz+p.yi*Nz+p.zi] + NDE[p.xi*Ny*Nz+p.yi*Nz+p.zi])/mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi])
                *s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi], ATP_to_D*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
    for (int i=0;i<6;i++){
      if ( s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ATP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
        ws.update( difD*( s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ATP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] ) / dA,
                   (X_ATP_pos + 1)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
      else {
        ws.update( 0.0, (X_ATP_pos + 1)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
    }
    switch (p.reaction) {
      case X_ATP_pos:
        p.xi += 1; break;
      case Y_ATP_pos:
        p.yi += 1; break;
      case Z_ATP_pos:
        p.zi += 1; break;
      case X_ATP_neg:
        p.xi -= 1; break;
      case Y_ATP_neg:
        p.yi -= 1; break;
      case Z_ATP_neg:
        p.zi -= 1; break;
    }
    s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] += 1;
    if (mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] != 0.0) {
      ws.update((rate_D + rate_dD*(ND[p.xi*Ny*Nz+p.yi*Nz+p.zi] + NDE[p.xi*Ny*Nz+p.yi*Nz+p.zi])/mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi])
                *s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi], ATP_to_D*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
    for (int i=0;i<6;i++){
      if ( s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ATP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
        ws.update( difD*( s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ATP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] ) / dA,
                   (X_ATP_pos + 1)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
      else {
        ws.update( 0.0, (X_ATP_pos + 1)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
    }
    return;
  }
  else if (p.reaction >= X_E_pos && p.reaction <= Z_E_neg) { //s_N_E has diffused
    s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] -= 1;
    if (mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] != 0.0) {
      ws.update(rate_E*ND[p.xi*Ny*Nz+p.yi*Nz+p.zi]/mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi]*s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi], E_D_to_DE*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
    for (int i=0;i<6;i++){
      if ( s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_E[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
        ws.update( difD*( s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_E[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] ) / dA,
                   (X_E_pos + 1)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
      else {
        ws.update( 0.0, (X_E_pos + 1)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
    }
    switch (p.reaction) {
      case X_E_pos:
        p.xi += 1; break;
      case Y_E_pos:
        p.yi += 1; break;
      case Z_E_pos:
        p.zi += 1; break;
      case X_E_neg:
        p.xi -= 1; break;
      case Y_E_neg:
        p.yi -= 1; break;
      case Z_E_neg:
        p.zi -= 1; break;
    }
    s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] += 1;
    if (mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] != 0.0) {
      ws.update(rate_E*ND[p.xi*Ny*Nz+p.yi*Nz+p.zi]/mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi]*s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi], E_D_to_DE*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
    for (int i=0;i<6;i++){
      if ( s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_E[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
        ws.update( difD*( s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_E[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] ) / dA,
                   (X_E_pos + 1)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
      else {
        ws.update( 0.0, (X_E_pos + 1)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      }
    }
    return;
  }
}





void set_membrane(double mem_A[]) {
  //int count = 0;
  for(int xi=0;xi<Nx;xi++){
    for(int yi=0;yi<Ny;yi++){
      for(int zi=0;zi<Nz;zi++){
        //count++;
        double fXYZ = mem_f((xi+0.5)*dx, (yi+0.5)*dx, (zi+0.5)*dx);
        double fXYz = mem_f((xi+0.5)*dx, (yi+0.5)*dx, (zi-0.5)*dx);
        double fXyZ = mem_f((xi+0.5)*dx, (yi-0.5)*dx, (zi+0.5)*dx);
        double fXyz = mem_f((xi+0.5)*dx, (yi-0.5)*dx, (zi-0.5)*dx);
        double fxYZ = mem_f((xi-0.5)*dx, (yi+0.5)*dx, (zi+0.5)*dx);
        double fxYz = mem_f((xi-0.5)*dx, (yi+0.5)*dx, (zi-0.5)*dx);
        double fxyZ = mem_f((xi-0.5)*dx, (yi-0.5)*dx, (zi+0.5)*dx);
        double fxyz = mem_f((xi-0.5)*dx, (yi-0.5)*dx, (zi-0.5)*dx);
        double f = mem_f(xi*dx, yi*dx, zi*dx);
        //printf("In set_membrane xi=%d yi=%d zi=%d and count = and Nx*Ny*Nz=%d\n",xi,yi,zi,Nx*Ny*Nz);
        //fflush(stdout);
        mem_A[xi*Ny*Nz+yi*Nz+zi] = find_intersection(fXYZ, fXYz, fXyZ, fXyz, fxYZ, fxYz, fxyZ, fxyz, f, false);
      }
    }
  }
}


double find_intersection(const double fXYZ, const double fXYz, const double fXyZ, const double fXyz,
                         const double fxYZ, const double fxYz, const double fxyZ, const double fxyz,
                         double f, bool debugme) {
  //Instead of evaluating mem_f in the middle of each cube we will average the corners for greater accuracy:
  f =  (fXYZ + fXYz + fXyZ + fXyz + fxYZ + fxYz + fxyZ + fxyz + f)/9;
  double dA = 0;
  if (debugme) printf("I am debugging\n");
  // "pts" is a set of points in 3D space (in units of distance) where
  // the plane of f=0 intersects with the *edges* of the cube.
  double *ptsx = new double[8];
  double *ptsy = new double[8];
  double *ptsz = new double[8];
  for (int i=0;i<8;i++) ptsx[i] = 0;
  for (int i=0;i<8;i++) ptsy[i] = 0;
  for (int i=0;i<8;i++) ptsz[i] = 0;
  int np = 0; // np is the number of intersections between edges of the cube and the plane f=0
  double df_dx = (fXYZ + fXYz + fXyZ + fXyz - fxYZ - fxyZ - fxYz - fxyz)/(4*dx);
  double df_dy = (fXYZ + fXYz + fxYZ + fxYz - fXyZ - fxyZ - fXyz - fxyz)/(4*dx);
  double df_dz = (fXYZ + fXyZ + fxYZ + fxyZ - fXYz - fXyz - fxYz - fxyz)/(4*dx);
  if (debugme) {
    printf("df_dx = %g\n", df_dx);
    printf("df_dy = %g\n", df_dy);
    printf("df_dz = %g\n", df_dz);
  }
  for (double j=-0.5; j<1.0; j++){
    for (double k=-0.5; k<1.0; k++){
      if ((-f - j*dx*df_dy - k*dx*df_dz)/(df_dx*dx) < 0.5 && (-f - j*dx*df_dy - k*dx*df_dz)/(df_dx*dx) > -0.5){
        ptsx[np] = (-f - j*dx*df_dy - k*dx*df_dz)/df_dx;
        ptsy[np] = j*dx;
        ptsz[np] = k*dx;
        np++;
      }
    }
  }
  //printf("df_dy = %g and f = %g and -f/(df_dy*dx) = %g\n",df_dy,f,-f/(df_dy*dx));
  for (double i=-0.5; i<1.0; i++){
    for (double k=-0.5; k<1.0; k++){
      if ((-f - i*dx*df_dx - k*dx*df_dz)/(df_dy*dx) < 0.5 && (-f - i*dx*df_dx - k*dx*df_dz)/(df_dy*dx) > -0.5){
        ptsy[np] = (-f - i*dx*df_dx - k*dx*df_dz)/df_dy;
        ptsx[np] = i*dx;
        ptsz[np] = k*dx;
        np++;
      }
    }
  }
  for (double j=-0.5; j<1.0; j++){
    for (double i=-0.5; i<1.0; i++){
      if ((-f - j*dx*df_dy - i*dx*df_dx)/(df_dz*dx) < 0.5 && (-f - j*dx*df_dy - i*dx*df_dx)/(df_dz*dx) > -0.5){
        ptsz[np] = (-f - j*dx*df_dy - i*dx*df_dx)/df_dz;
        ptsy[np] = j*dx;
        ptsx[np] = i*dx;
        np++;
      }
    }
  }
  //printf("np = %d\n",np);
  if (np == 0) return 0.0; // no intersections ===> no area!
  if (debugme) {
    //printf("np = %d\n", np);
    for (int i=0;i<np;i++) {
      printf("\t%g %g %g\n", ptsx[i], ptsy[i], ptsz[i]);
    }
  }
  for (int i=0; i<np;i++){
    //printf("ptsx[%d]=%g, ptsy[%d]=%g, ptsy[%d]=%g\n",i,ptsx[i],i,ptsy[i],i,ptsz[i]);
  }
  int nz0 = 0;  // how many of our pts are on the negative z side
  int ny0 = 0;
  int nx0 = 0;
  // "pt" is the same set of points in 3D space as pts, but reordered
  // such that they go around the polygon of intersection in order.
  double *ptx = new double[8];
  double *pty = new double[8];
  double *ptz = new double[8];
  for (int i=0;i<8;i++) ptx[i] = 0;
  for (int i=0;i<8;i++) pty[i] = 0;
  for (int i=0;i<8;i++) ptz[i] = 0;
  double *line = new double[8];
  for (int i=0;i<8;i++) line[i] = 0;
  double *cos = new double[8];
  double cos_max = -2;
  int as=0;
  int bs=0;
  int cs=0;
  int ds=0;
  int s = 1;
  double eline = 0;
  double p = 0;
  for (int i=0;i<np;i++){
    if (ptsz[i] == -0.5*dx) {ptz[nz0]=ptsz[i]; ptx[nz0]=ptsx[i]; pty[nz0]=ptsy[i]; nz0++;}
  }
  // at this point, we know how many times the f=0 plane intersects
  // the edgest of the negative z square, which will be either 0, 2,
  // or 4.
  if (debugme) {
    printf("nz0 = %d\n", nz0);
  }
  if (nz0 == 4) {
    // the plane is coplanar with one edge, so our area is very easy.
    return dx*dx;
  }
  if (nz0 == 2){
    line[0] = sqrt((ptx[1]-ptx[0])*(ptx[1]-ptx[0]) + (pty[1]-pty[0])*(pty[1]-pty[0])
                   + (ptz[1]-ptz[0])*(ptz[1]-ptz[0]));
    for (int i=0;i<np;i++){
      if (ptsz[i] != -0.5*dx){
        line[s] = sqrt((ptsx[i]-ptx[0])*(ptsx[i]-ptx[0]) + (ptsy[i]-pty[0])*(ptsy[i]-pty[0])
                       + (ptsz[i]-ptz[0])*(ptsz[i]-ptz[0]));
        cos[s] = ((ptsx[i]-ptx[0])*(ptx[1]-ptx[0]) + (ptsy[i]-pty[0])*(pty[1]-pty[0])
                  + (ptsz[i]-ptz[0])*(ptz[1]-ptz[0])) / (line[s]*line[0]);
        ptz[s+1] = ptsz[i]; ptx[s+1] = ptsx[i]; pty[s+1] = ptsy[i];
        s++;
      }
    }
  } else {
    // The plane doesn't intersect with the negative z side of the
    // cube, so let's look at the negative *x* side...
    for (int i=0;i<np;i++){
      if (ptsx[i] == -0.5*dx) {ptz[nx0]=ptsz[i]; ptx[nx0]=ptsx[i]; pty[nx0]=ptsy[i]; nx0++;}
    }
    if (debugme) {
      printf("nx0 = %d\n", nx0);
    }
    if (nx0 == 4) {
      return dx*dx; // coplanar, as before
    }
    if (nx0 == 2){
      line[0] = sqrt((ptx[1]-ptx[0])*(ptx[1]-ptx[0]) + (pty[1]-pty[0])*(pty[1]-pty[0])
                     + (ptz[1]-ptz[0])*(ptz[1]-ptz[0]));
      for (int i=0;i<np;i++){
        if (ptsx[i] != -0.5*dx){
          line[s] = sqrt((ptsx[i]-ptx[0])*(ptsx[i]-ptx[0]) + (ptsy[i]-pty[0])*(ptsy[i]-pty[0])
                         + (ptsz[i]-ptz[0])*(ptsz[i]-ptz[0]));
          cos[s] = ((ptsx[i]-ptx[0])*(ptx[1]-ptx[0]) + (ptsy[i]-pty[0])*(pty[1]-pty[0])
                    + (ptsz[i]-ptz[0])*(ptz[1]-ptz[0])) / (line[s]*line[0]);
          ptz[s+1] = ptsz[i]; ptx[s+1] = ptsx[i]; pty[s+1] = ptsy[i];
          s++;
        }
      }
    } else {
      // No intersection with negative z or negative x side of the cube...
      for (int i=0;i<np;i++){
        if (ptsy[i] == -0.5*dx) {ptz[ny0]=ptsz[i]; ptx[ny0]=ptsx[i]; pty[ny0]=ptsy[i]; ny0++;}
      }
      if (debugme) {
        printf("ny0 = %d\n", ny0);
      }
      if (ny0 == 4) {
        return dx*dx; // coplanar, as before
      }
      if (ny0 == 2){
        line[0] = sqrt((ptx[1]-ptx[0])*(ptx[1]-ptx[0]) + (pty[1]-pty[0])*(pty[1]-pty[0])
                       + (ptz[1]-ptz[0])*(ptz[1]-ptz[0]));
        for (int i=0;i<np;i++){
          if (ptsy[i] != -0.5*dx){
            line[s] = sqrt((ptsx[i]-ptx[0])*(ptsx[i]-ptx[0]) + (ptsy[i]-pty[0])*(ptsy[i]-pty[0])
                           + (ptsz[i]-ptz[0])*(ptsz[i]-ptz[0]));
            cos[s] = ((ptsx[i]-ptx[0])*(ptx[1]-ptx[0]) + (ptsy[i]-pty[0])*(pty[1]-pty[0])
                      + (ptsz[i]-ptz[0])*(ptz[1]-ptz[0])) / (line[s]*line[0]);
            ptz[s+1] = ptsz[i]; ptx[s+1] = ptsx[i]; pty[s+1] = ptsy[i];
            s++;
          }
        }
      } else {
        // We now know that the plane must be going through the +x, +y, +z corner!
        assert(np == 3);
        //return 0.0;
        const double dist01 = sqrt((ptsx[1]-ptsx[0])*(ptsx[1]-ptsx[0])
                                   + (ptsy[1]-ptsy[0])*(ptsy[1]-ptsy[0])
                                   + (ptsz[1]-ptsz[0])*(ptsz[1]-ptsz[0]));
        const double dist02 = sqrt((ptsx[2]-ptsx[0])*(ptsx[2]-ptsx[0])
                                   + (ptsy[2]-ptsy[0])*(ptsy[2]-ptsy[0])
                                   + (ptsz[2]-ptsz[0])*(ptsz[2]-ptsz[0]));
        const double dot012 = (ptsx[1]-ptsx[0])*(ptsx[2]-ptsx[0]) + (ptsy[1]-ptsy[0])*(ptsy[2]-ptsy[0]) + (ptsz[1]-ptsz[0])*(ptsz[2]-ptsz[0]);
        return 0.5*sqrt(dist01*dist01*dist02*dist02 - dot012*dot012);
      }
    }
  }
  assert(s == np-1);
  for (int i=1;i<s;i++){
    if (cos[i] > cos_max){cos_max = cos[i]; as=i;}
  }
  if (debugme) {
    printf("s = %d\n", s);
  }
  eline = sqrt((ptx[as+1]-ptx[1])*(ptx[as+1]-ptx[1]) + (pty[as+1]-pty[1])*(pty[as+1]-pty[1])
               + (ptz[as+1]-ptz[1])*(ptz[as+1]-ptz[1]));
  p = (line[as]+line[0]+eline)/2;
  dA = sqrt(p*(p-line[as])*(p-line[0])*(p-eline));
  //printf("First dA = %g\n",dA);
  if (debugme) {
    printf("dA = %g\n", dA);
    printf("p = %g\n", p);
    printf("eline = %g\n", eline);
    printf("line[as] = %g\n", line[as]);
    printf("line[0] = %g\n", line[0]);
    printf("as = %d\n", as);
  }
  cos[as] = -2; cos_max = -2;

  if (np==4 || np==5 || np==6){
    for (int i=1;i<s;i++){
      if (cos[i] > cos_max){cos_max = cos[i]; bs=i;}
    }
    eline = sqrt((ptx[bs+1]-ptx[as+1])*(ptx[bs+1]-ptx[as+1]) + (pty[bs+1]-pty[as+1])*(pty[bs+1]-pty[as+1])
                 + (ptz[bs+1]-ptz[as+1])*(ptz[bs+1]-ptz[as+1]));
    p = (line[bs]+line[as]+eline)/2;
    dA += sqrt(p*(p-line[bs])*(p-line[as])*(p-eline));
    //printf("Second dA = %g\n",sqrt(p*(p-line[bs])*(p-line[as])*(p-eline)));
    cos[bs] = -2; cos_max = -2;
  }

  if (np==5 || np==6){
    for (int i=1;i<s;i++){
      if (cos[i] > cos_max){cos_max = cos[i]; cs=i;}
    }
    eline = sqrt((ptx[cs+1]-ptx[bs+1])*(ptx[cs+1]-ptx[bs+1]) + (pty[cs+1]-pty[bs+1])*(pty[cs+1]-pty[bs+1])
                 + (ptz[cs+1]-ptz[bs+1])*(ptz[cs+1]-ptz[bs+1]));
    p = (line[cs]+line[bs]+eline)/2;
    dA += sqrt(p*(p-line[cs])*(p-line[bs])*(p-eline));
    //printf("Third dA = %g\n",sqrt(p*(p-line[cs])*(p-line[bs])*(p-eline)));
    cos[cs] = -2; cos_max = -2;
  }

  if (np == 6) {
    for (int i=1;i<s;i++){
      if (cos[i] > cos_max){cos_max = cos[i]; ds=i;}
    }
    eline = sqrt((ptx[ds+1]-ptx[cs+1])*(ptx[ds+1]-ptx[cs+1]) + (pty[ds+1]-pty[cs+1])*(pty[ds+1]-pty[cs+1])
                 + (ptz[ds+1]-ptz[cs+1])*(ptz[ds+1]-ptz[cs+1]));
    p = (line[ds]+line[cs]+eline)/2;
    dA += sqrt(p*(p-line[ds])*(p-line[cs])*(p-eline));
    //printf("Fourth dA = %g\n",sqrt(p*(p-line[ds])*(p-line[cs])*(p-eline)));
  }
  delete[] ptsx;
  delete[] ptsy;
  delete[] ptsz;
  delete[] ptx;
  delete[] pty;
  delete[] ptz;
  delete[] line;
  delete[] cos;
  return dA;
}

void set_insideArr(bool *insideArr){
  for(int xi=0;xi<Nx;xi++){
    for(int yi=0;yi<Ny;yi++){
      for(int zi=0;zi<Nz;zi++){
        insideArr[xi*Ny*Nz+yi*Nz+zi] = inside(xi,yi,zi);
      }
    }
  }
}

bool inside(int xi, int yi, int zi){
  if (once_inside==true){
    printf("\nxi = %d yi = %d zi = %d\n",xi,yi,zi);
  }
  xi = xi + min_xi-1;
  yi = yi + min_yi-1;
  zi = zi + min_zi-1;
  if (once_inside==true){
    printf("xi = %d yi = %d zi = %d\n\n",xi,yi,zi);
    once_inside = false;
  }
  fflush(stdout);
  if (mem_f((xi-0.5)*dx,(yi-0.5)*dx,(zi-0.5)*dx) <= 0) {return true;}
  if (mem_f((xi+0.5)*dx,(yi-0.5)*dx,(zi-0.5)*dx) <= 0) {return true;}
  if (mem_f((xi-0.5)*dx,(yi+0.5)*dx,(zi-0.5)*dx) <= 0) {return true;}
  if (mem_f((xi-0.5)*dx,(yi-0.5)*dx,(zi+0.5)*dx) <= 0) {return true;}
  if (mem_f((xi-0.5)*dx,(yi+0.5)*dx,(zi+0.5)*dx) <= 0) {return true;}
  if (mem_f((xi+0.5)*dx,(yi-0.5)*dx,(zi+0.5)*dx) <= 0) {return true;}
  if (mem_f((xi+0.5)*dx,(yi+0.5)*dx,(zi-0.5)*dx) <= 0) {return true;}
  if (mem_f((xi+0.5)*dx,(yi+0.5)*dx,(zi+0.5)*dx) <= 0) {return true;}
  return false;
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

int set_density(double *nATP, double *nADP, double *nE, double *ND, double *NDE, int *s_ND, int *s_NDE, int *s_N_ATP, int *s_N_ADP, int *s_N_E, double *mem_A){
  double dV = dx*dx*dx;
  printf("In set_density function, Nx = %d Ny = %d Nz = %d\n",Nx,Ny,Nz);
  int right_most_point_z=0; //left and right most points for z
  int left_most_point_z=Nz;
  int right_most_point_y=0; //"left" and "right" most points for y, in terms of magnitude (right = larger y value)
  int left_most_point_y=Ny;
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (inside(i,j,k)){
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
  //get total gridpoints, gridpoints left of divide, gridpoints right of divide for protein count
  double vert_div = 2.8/dx-min_zi+1;
  double hor_div = 3.0/dx-min_yi+1;
  int density_divider_right = int(right_most_point_z - (right_most_point_z - left_most_point_z)/3);
  // printf("right_most_point_z = %d and left_most_point_z = %d and density_divider_right = %d and Nz = %d\n",
  //       right_most_point_z,left_most_point_z,density_divider_right,Nz);
  int density_divider_left = int(right_most_point_z - 2*(right_most_point_z - left_most_point_z)/3);
  int gridpoints_low_dens = 0;
  int gridpoints_high_dens = 0;
  int gridpoints_total = 0;
  double wall_area_high = 0;
  double wall_area_low = 0;
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (inside(i,j,k)){
          gridpoints_total++;
          if ((rand_seed == 96 && k>vert_div && j<hor_div) || (rand_seed != 96 && k>density_divider_right)) {
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
  //compute density scale factors low_dens and high_dens of divide (to ensure correct protein #)
  double density_factor_low_dens = gridpoints_total/(gridpoints_low_dens + density_factor*gridpoints_high_dens);
  double density_factor_high_dens = density_factor*gridpoints_total/(gridpoints_low_dens + density_factor*gridpoints_high_dens);

  //double MinD_proteins_per_area_high = dx*dx*dx*nATP_starting_density*gridpoints_high_dens*density_factor_high_dens/wall_area_high;
  //double MinD_proteins_per_area_low = dx*dx*dx*nATP_starting_density*gridpoints_low_dens*density_factor_low_dens/wall_area_low;
  // printf("wall_area_high = %g and proteins_per_area_high = %g\n",wall_area_high,MinD_proteins_per_area_high);
  // printf("Density factors: low_dens: %f, high_dens: %f, ratio: %f\n", density_factor_low_dens, density_factor_high_dens, density_factor_hig
  //      h_dens/density_factor_low_dens);

  // printf("Gridpoints low_dens of the divider: %d\n",gridpoints_low_dens);
  // printf("Gridpoints high_dens of the divider:%d\n",gridpoints_high_dens);
  // printf("Gridpoints total: %d\n",gridpoints_total);

  //begin setting density at each gridpoint:
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (inside(i,j,k)){
          if ((rand_seed == 96 && k>vert_div && j<hor_div) || (rand_seed != 96 && k>density_divider_right)) {
            ND[i*Ny*Nz+j*Nz+k] = 0;//MinD_proteins_per_area_high*mem_A[i*Ny*Nz+j*Nz+k]/2.0;
            nE[i*Ny*Nz+j*Nz+k] = nE_starting_density*density_factor_high_dens;
            NflE[i*Ny*Nz+j*Nz+k] = NDE[i*Ny*Nz+j*Nz+k] + nE[i*Ny*Nz+j*Nz+k]*dV;
          }
          else {
            ND[i*Ny*Nz+j*Nz+k] = 0;//MinD_proteins_per_area_low*mem_A[i*Ny*Nz+j*Nz+k]/2.0;
            nE[i*Ny*Nz+j*Nz+k] = nE_starting_density*density_factor_low_dens;
            NflE[i*Ny*Nz+j*Nz+k] = NDE[i*Ny*Nz+j*Nz+k] + nE[i*Ny*Nz+j*Nz+k]*dV;
          }
          // if (mem_f_shape=="triangle" && j > (B*sqrt(3)/2.0)/2.0){
          //   nATP[i*Ny*Nz+j*Nz+k] = nATP_starting_density*density_factor_high_dens;
          //   NflD[i*Ny*Nz+j*Nz+k] = ND[i*Ny*Nz+j*Nz+k] + NDE[i*Ny*Nz+j*Nz+k] + (nATP[i*Ny*Nz+j*Nz+k] + nADP[i*Ny*Nz+j*Nz+k])*dV;
          // }
        }
        else {
          nATP[i*Ny*Nz+j*Nz+k] = 0;
          nE[i*Ny*Nz+j*Nz+k] = 0;
          nADP[i*Ny*Nz+j*Nz+k] =0;
          ND[i*Ny*Nz+j*Nz+k] =0;
          NDE[i*Ny*Nz+j*Nz+k] = 0;
          NflD[i*Ny*Nz+j*Nz+k] = ND[i*Ny*Nz+j*Nz+k] + NDE[i*Ny*Nz+j*Nz+k] + (nATP[i*Ny*Nz+j*Nz+k] + nADP[i*Ny*Nz+j*Nz+k])*dV;
          NflE[i*Ny*Nz+j*Nz+k] = NDE[i*Ny*Nz+j*Nz+k] + nE[i*Ny*Nz+j*Nz+k]*dV;
        }
      }
    }
  }

  gridpoints_low_dens = 0;
  gridpoints_high_dens = 0;
  gridpoints_total = 0;
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (inside(i,j,k)){
          gridpoints_total++;
          if ((rand_seed == 96 && k<vert_div && j<hor_div) || (rand_seed != 96 && k<density_divider_left)) {
            gridpoints_high_dens++;
          }
          else {
            gridpoints_low_dens++;
          }
        }
      }
    }
  }

  //compute density scale factors low_dens and high_dens of divide (to ensure correct protein #)
  density_factor_low_dens = gridpoints_total/(gridpoints_low_dens + density_factor*gridpoints_high_dens);
  density_factor_high_dens = density_factor*gridpoints_total/(gridpoints_low_dens + density_factor*gridpoints_high_dens);

  //  printf("Density factors: low_dens: %f, high_dens: %f, ratio: %f\n", density_factor_low_dens, density_factor_high_dens, density_factor_high_dens/density_factor_low_dens);

  // printf("Gridpoints low_dens of the divider: %d\n",gridpoints_low_dens);
  // printf("Gridpoints high_dens of the divider:%d\n",gridpoints_high_dens);
  // printf("Gridpoints total: %d\n",gridpoints_total);

  //begin setting density at each gridpoint:
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (inside(i,j,k)){
          if ((rand_seed == 96 && k<vert_div && j<hor_div) || (rand_seed != 96 && k<density_divider_left)) {
            if (mem_f_shape=="stad"){
              if (j > Ny/2) {
                nATP[i*Ny*Nz+j*Nz+k] = 2*nATP_starting_density*density_factor_high_dens;
              }
              else {
                nATP[i*Ny*Nz+j*Nz+k] = 0;
              }
            }
            else {
              nATP[i*Ny*Nz+j*Nz+k] = nATP_starting_density*density_factor_high_dens;
            }
          }
          else {
            nATP[i*Ny*Nz+j*Nz+k] = nATP_starting_density*density_factor_low_dens;
          }
          nADP[i*Ny*Nz+j*Nz+k] =0;
          NDE[i*Ny*Nz+j*Nz+k] = 0;
          NflD[i*Ny*Nz+j*Nz+k] = ND[i*Ny*Nz+j*Nz+k] + NDE[i*Ny*Nz+j*Nz+k] + (nATP[i*Ny*Nz+j*Nz+k] + nADP[i*Ny*Nz+j*Nz+k])*dV;
        }
      }
    }
  }
  if (stochastic) {
    for (int i=0; i<Nx*Ny*Nz;i++) {
      s_N_ATP[i] = round(nATP[i]*dV);
      s_N_ADP[i] = round(nADP[i]*dV);
      s_N_E[i] = round(nATP[i]*dV);
      s_ND[i] = round(ND[i]);
      s_NDE[i] = round(NDE[i]);
    }
  }
  return 0;
}






// int get_next_stochastic_state(double *mem_A, bool *insideArr, double *nATP, double *nADP,
//                               double *nE, double *ND, double *NDE, double *NflD, double *NflE,
//                               double *JxATP, double *JyATP, double *JzATP,
//                               double *JxADP, double *JyADP, double *JzADP,
//                               double *JxE, double *JyE, double *JzE) {
//   double dV = dx*dx*dx;
//   double ADP_to_ATP;
//   double de_to_ADP_E;
//   double ATP_to_d;
//   double E_d_to_de;
//   double *W_arr = new double[4*Nx*Ny*Nz+1];
//   W_arr[0] = 0;
//   double dist_value = 0;
//   for(int xi=0;xi<Nx;xi++){
//     for(int yi=0;yi<Ny;yi++){
//       for(int zi=0;zi<Nz;zi++){
//         ADP_to_ATP = rate_ADP_ATP*nADP[xi*Ny*Nz+yi*Nz+zi]*dV;
//         W_arr[4*(xi*Ny*Nz+yi*Nz+zi)+1] = dist_value + ADP_to_ATP;
//         de_to_ADP_E = rate_de*NDE[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*dV;
//         W_arr[4*(xi*Ny*Nz+yi*Nz+zi)+2] = dist_value + ADP_to_ATP + de_to_ADP_E;
//         ATP_to_d = (rate_D + rate_dD*(ND[xi*Ny*Nz+yi*Nz+zi] + NDE[xi*Ny*Nz+yi*Nz+zi])/mem_A[xi*Ny*Nz+yi*Nz+zi])
//           *nATP[xi*Ny*Nz+yi*Nz+zi]*dV;
//         W_arr[4*(xi*Ny*Nz+yi*Nz+zi)+3] = dist_value + ADP_to_ATP + de_to_ADP_E + ATP_to_d;
//         E_d_to_de = rate_E*ND[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*nE[xi*Ny*Nz+yi*Nz+zi]*dV;
//         W_arr[4*(xi*Ny*Nz+yi*Nz+zi)+4] = dist_value + ADP_to_ATP + de_to_ADP_E + ATP_to_d + E_d_to_de;
//         dist_value += ADP_to_ATP + de_to_ADP_E + ATP_to_d + E_d_to_de;
//       }
//     }
//   }
//   double C_fun = dist_value;
//   for (int i =0;i<Nx*Ny*Nz-1;i++){
//     if (W_arr[i+1] < W_arr[i]) {
//       printf("Theres a problem in the W_arr assignment loops/n");
//       exit(1);
//     }
//     if (W_arr[4*Nx*Ny*Nz] != C_fun) {
//       printf("Theres a problem with the C_fun assignment\n");
//       exit(1);
//     }
//   }
//   double max_W = 0;
//   for (int i=0; i<4*Nx*Ny*Nz+1; i++) {
//     if (W_arr[i] > max_W) {
//       max_W = W_arr[i];
//     }
//   }
//   //starting the loop that I think hopefully we can run a number of times before updating W_arr?
//   double elapsed_time=0;
//   while (stoch_iter*time_step > elapsed_time) {
//     double ran = (double)rand()/(RAND_MAX)*C_fun; //between 0 and C_fun
//     int interval = int(ran/C_fun*4*Nx*Ny*Nz);
//     while (ran < W_arr[interval] or ran >= W_arr[interval+1]) {
//       if ( abs(ran - W_arr[interval]) < max_W) {
//         if (ran - W_arr[interval] > 0) {
//           interval++;
//         }
//         else {
//           interval--;
//         }
//       }
//       else {
//         interval = interval + (ran - W_arr[interval])/C_fun*4*Nx*Ny*Nz;
//         if (interval < 0) {
//           interval = -interval;
//         }
//       }
//     }
//     printf("ran = %f\nW_arr[interval] = %f\nW_arr[interval+1] = %f\n",ran,W_arr[interval],W_arr[interval+1]);
//     printf("Ran should be in between the other two!\n");
//     if (ran >= W_arr[interval+1] or ran < W_arr[interval]) {
//       printf("And its not so Im quitting :(");
//       exit(1);
//     }
//     int xi = floor( interval/(4*Ny*Nz) );
//     int yi = floor( (interval - xi*Ny*Nz)/(4*Nz) );
//     int zi = floor( (interval - xi*Ny*Nz - yi*Nz)/4 );
//     int n = interval -xi*Ny*Nz - yi*Nz - zi;
//     //should I worry about these floors when the integer is right on?
//     if (n == 0) {
//       nADP[xi*Ny*Nz+yi*Nz+zi] -= 1/dV;
//       nATP[xi*Ny*Nz+yi*Nz+zi] += 1/dV;
//     }
//     if (n == 1) {
//       NDE[xi*Ny*Nz+yi*Nz+zi] -= 1;
//       nADP[xi*Ny*Nz+yi*Nz+zi] += 1/dV;
//       nE[xi*Ny*Nz+yi*Nz+zi] += 1/dV;
//     }
//     if (n == 2) {
//       nATP[xi*Ny*Nz+yi*Nz+zi] -= 1/dV;
//       ND[xi*Ny*Nz+yi*Nz+zi] += 1;
//     }
//     if (n == 3) {
//       nE[xi*Ny*Nz+yi*Nz+zi] -= 1/dV;
//       ND[xi*Ny*Nz+yi*Nz+zi] -= 1;
//       NDE[xi*Ny*Nz+yi*Nz+zi] += 1;
//     }
//     double delta_t = log( (double)rand()/(RAND_MAX) ) / C_fun;
//     elapsed_time += delta_t;
//   }
//   delete[]  W_arr;
//   return 0;
// }


// void set_curvature(double mem_A[], double curvature[]){
//   printf("doing set curvature!!!\n");
//   //fflush(stdout);
//   double X = Nx*dx;
//   double x1 = (X-A)/2.0;
//   double x2 = (X+A)/2.0;
//   for(int xi=0;xi<Nx;xi++){
//     for(int yi=0;yi<Ny;yi++){
//       for(int zi=0;zi<Nz;zi++){
//         if (mem_A[xi*Ny*Nz+yi*Nz+zi]==0 || (xi*dx+0.05)>x2 || (xi*dx-0.05)<x1 ){
//           curvature[xi*Ny*Nz+yi*Nz+zi]=0;
//         } else {
//           double fX = mem_f((xi+0.5)*dx, yi*dx, zi*dx);
//           double fx = mem_f((xi-0.5)*dx, yi*dx, zi*dx);
//           double fY = mem_f(xi*dx, (yi+0.5)*dx, zi*dx);
//           double fy = mem_f(xi*dx, (yi-0.5)*dx, zi*dx);
//           double fZ = mem_f(xi*dx, yi*dx, (zi+0.5)*dx);
//           double fz = mem_f(xi*dx, yi*dx, (zi-0.5)*dx);
//           double f = mem_f(xi*dx, yi*dx, zi*dx);

//           double df_dx = (fX-fx)/dx;
//           double df_dy = (fY-fy)/dx;
//           double df_dz = (fZ-fz)/dx;
//           double constant = sqrt((df_dx)*(df_dx) + (df_dy)*(df_dy) + (df_dz)*(df_dz));

//           curvature[xi*Ny*Nz+yi*Nz+zi] = 4*(fX + fx + fY + fy + fZ + fz - 6*f)/dx/dx/constant;
//         }
//       }
//     }
//   }
//   //fflush(stdout);
// }
