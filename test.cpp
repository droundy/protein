#include <iostream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct things {
  int mynum;
  int * pnum;
  char * mystring;
  string otherstring;
};

int main() {
  printf("Hello\n");
//   // things people;
//   // things * objects;
//   // people.mynum = 5;
//   // objects->mynum = 6;
//   // people.mystring = "this is my string";
//   // printf("%s",people.mystring);
//   int Nx = 5;
//   int Ny = 5;
//   int Nz = 5;
//   double num = 0;
//   double ran = (double)rand()/(RAND_MAX); //between 0 and 1
//   for(int xi=0;xi<Nx;xi++){
//     for(int yi=0;yi<Ny;yi++){
//       for(int zi=0;zi<Nz;zi++){
//         num += ran*.0001;
//         printf("%d %d %d %f  ",xi,yi,zi,num);
//         if (num < 1.0) {
//           break;
//         }
//       }
//     }
//   }
//   printf("All done!!\n");
  enum MyEnum {zero, one, two, three};
  printf("Does this work? %d %d %d and size of = %f\n",one*2,two*2,zero*2,sizeof(zero));


  int number = 2;
  printf("Hello? %d %d\n",number, two);
  if (two <= two) {
    printf("the boolean worked so Im here! %d %d\n",number,two);
  }
  if (two >= two) {
    printf("the other boolean worked so Im here! %d %d\n",number,two);
  }

  int d [6][3] =
    {
      {1,0,0},
      {-1,0,0},
      {0,1,0},
      {0,-1,0},
      {0,0,1},
      {0,0,-1},
    };
  for (int i=0;i<6;i++){
    for (int j=0;j<3;j++) {
      printf("i=%d j=%d is %d\n",i,j,d[i][j]);
    }
  }
  int const jeff = 5;
  int dan = 4;
  const int * krebs = &jeff;
  krebs = &dan;
  printf("jeff and krebs are %d and %d\n",jeff, *krebs);

  return 0;
}



// int get_next_stochastic_state(double *mem_A, bool *insideArr, double *nATP, double *nADP,
// 			      double *nE, double *Nd, double *Nde, double *NflD, double *NflE,
// 			      double *JxATP, double *JyATP, double *JzATP,
// 			      double *JxADP, double *JyADP, double *JzADP,
//                               double *JxE, double *JyE, double *JzE, double elapsed_time) {
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
// 	de_to_ADP_E = rate_de*NDE[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*dV;
// 	W_arr[4*(xi*Ny*Nz+yi*Nz+zi)+2] = dist_value + ADP_to_ATP + de_to_ADP_E;
//         ATP_to_d = (rate_D + rate_dD*(ND[xi*Ny*Nz+yi*Nz+zi] + NDE[xi*Ny*Nz+yi*Nz+zi])/mem_A[xi*Ny*Nz+yi*Nz+zi])
//           *nATP[xi*Ny*Nz+yi*Nz+zi]*dV;
// 	W_arr[4*(xi*Ny*Nz+yi*Nz+zi)+3] = dist_value + ADP_to_ATP + de_to_ADP_E + ATP_to_d;
//         E_d_to_de = rate_E*ND[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*nE[xi*Ny*Nz+yi*Nz+zi]*dV;
// 	W_arr[4*(xi*Ny*Nz+yi*Nz+zi)+4] = dist_value + ADP_to_ATP + de_to_ADP_E + ATP_to_d + E_d_to_de;
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
//   double ran = (double)rand()/(RAND_MAX)*C_fun; //between 0 and C_fun
//   int interval = int(ran/C_fun*4*Nx*Ny*Nz);
//   while (ran < W_arr[interval] or ran >= W_arr[interval+1]) {
//     if ( abs(ran - W_arr[interval]) < C_fun/Nx*Ny*Nz) { //should figure out the max in the array above and use that
//       if (ran - W_arr[interval] > 0) {
// 	interval++;
//       }
//       else {
// 	interval--;
//       }
//     }
//     else {
//       interval = interval + (ran - W_arr[interval])/C_fun*4*Nx*Ny*Nz;
//       if (interval < 0) {
// 	interval = -interval;
//       }
//     }
//   }
//   printf("ran = %f\nW_arr[interval] = %f\nW_arr[interval+1] = %f\n",ran,W_arr[interval],W_arr[interval+1]);
//   printf("Ran should be in between the other two!\n");
//   if (ran >= W_arr[interval+1] or ran < W_arr[interval]) {
//     printf("And its not so Im quitting :(");
//     exit(1);
//   }
//   int xi = floor( interval/(4*Ny*Nz) );
//   int yi = floor( (interval - xi*Ny*Nz)/(4*Nz) );
//   int zi = floor( (interval - xi*Ny*Nz - yi*Nz)/4 );
//   int n = interval -xi*Ny*Nz - yi*Nz - zi;
//   //should I worry about these floors when the integer is right on?

//   if (n == 0) {
//     nADP[xi*Ny*Nz+yi*Nz+zi] -= 1/dV;
//     nATP[xi*Ny*Nz+yi*Nz+zi] += 1/dV;
//   }
//     if (n == 1) {
//     Nde[xi*Ny*Nz+yi*Nz+zi] -= 1;
//     nADP[xi*Ny*Nz+yi*Nz+zi] += 1/dV;
//     nE[xi*Ny*Nz+yi*Nz+zi] += 1/dV;
//   }
//   if (n == 2) {
//     nATP[xi*Ny*Nz+yi*Nz+zi] -= 1/dV;
//     Nd[xi*Ny*Nz+yi*Nz+zi] += 1;
//   }
//   if (n == 3) {
//     nE[xi*Ny*Nz+yi*Nz+zi] -= 1/dV;
//     Nd[xi*Ny*Nz+yi*Nz+zi] -= 1;
//     Nde[xi*Ny*Nz+yi*Nz+zi] += 1;
//   }
//   double delta_t = log( (double)rand()/(RAND_MAX) ) / C_fun;
//   elapsed_time = += delta_t;
//   delete[]  W_arr;
//   return 0;
// }





// w_ADP_ATP  = rate_ADP_ATP*nADP[xi*Ny*Nz+yi*Nz+zi]

// for(int xi=0;xi<Nx;xi++){
//     for(int yi=0;yi<Ny;yi++){
//       for(int zi=0;zi<Nz;zi++){
//         ADP_to_ATP = rate_ADP_ATP*nADP[xi*Ny*Nz+yi*Nz+zi]*time_step;
//         de_to_ADP_E = rate_de*NDE[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*time_step;
//         ATP_to_d = (rate_D + rate_dD*(ND[xi*Ny*Nz+yi*Nz+zi] + NDE[xi*Ny*Nz+yi*Nz+zi])/mem_A[xi*Ny*Nz+yi*Nz+zi])
//           *nATP[xi*Ny*Nz+yi*Nz+zi]*time_step;
//         E_d_to_de = rate_E*ND[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*nE[xi*Ny*Nz+yi*Nz+zi]*time_step;
//         //Jeff!  remember that when you gain cyto density and lose the same amount of wall density,
//         //the numbers of proteins gained/lost will be different, and it's the numbers that you want to be the same!!
//         //also, keep thinking about the issue below where all additions are divided and then mult by the same mem_A fun
//         nADP[xi*Ny*Nz+yi*Nz+zi] -= ADP_to_ATP;
//         nATP[xi*Ny*Nz+yi*Nz+zi] += ADP_to_ATP;
//         if (mem_A[xi*Ny*Nz+yi*Nz+zi] != 0){
//           NDE[xi*Ny*Nz+yi*Nz+zi] += -de_to_ADP_E*mem_A[xi*Ny*Nz+yi*Nz+zi];
//           nADP[xi*Ny*Nz+yi*Nz+zi] += de_to_ADP_E*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);
//           nE[xi*Ny*Nz+yi*Nz+zi] += de_to_ADP_E*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);

//           nATP[xi*Ny*Nz+yi*Nz+zi] -= ATP_to_d*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);
//           ND[xi*Ny*Nz+yi*Nz+zi] += ATP_to_d*mem_A[xi*Ny*Nz+yi*Nz+zi];

//           nE[xi*Ny*Nz+yi*Nz+zi] -= E_d_to_de*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);
//           ND[xi*Ny*Nz+yi*Nz+zi] -= E_d_to_de*mem_A[xi*Ny*Nz+yi*Nz+zi];
//           NDE[xi*Ny*Nz+yi*Nz+zi] += E_d_to_de*mem_A[xi*Ny*Nz+yi*Nz+zi];
//         }
//         NflD[xi*Ny*Nz+yi*Nz+zi] = (nATP[xi*Ny*Nz+yi*Nz+zi] + nADP[xi*Ny*Nz+yi*Nz+zi])*(dx*dx*dx) + ND[xi*Ny*Nz+yi*Nz+zi] + NDE[xi*Ny*Nz+yi*Nz+zi];
//         NflE[xi*Ny*Nz+yi*Nz+zi] = nE[xi*Ny*Nz+yi*Nz+zi]*dx*dx*dx + NDE[xi*Ny*Nz+yi*Nz+zi];
//       }
//     }




// int get_next_density(double *mem_A, bool *insideArr, double *nATP, double *nADP,
//                      double *nE, double *Nd, double *Nde, double *NflD, double *NflE,
//                      double *JxATP, double *JyATP, double *JzATP,
//                      double *JxADP, double *JyADP, double *JzADP,
//                      double *JxE, double *JyE, double *JzE){
//   for(int xi=0;xi<Nx-1;xi++){
//     for(int yi=0;yi<Ny-1;yi++){
//       for(int zi=0;zi<Nz-1;zi++){
//         //for the diffusion terms, we use dn/dt = dA*J/dV thinking in terms of tot #, but dA/dV=1/dx
//         if (insideArr[(xi+1)*Ny*Nz+yi*Nz+zi] && insideArr[xi*Ny*Nz+yi*Nz+zi]){
//           nADP[(xi+1)*Ny*Nz+yi*Nz+zi] += JxADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//           nADP[xi*Ny*Nz+yi*Nz+zi] -= JxADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//           nATP[(xi+1)*Ny*Nz+yi*Nz+zi] += JxATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//           nATP[xi*Ny*Nz+yi*Nz+zi] -= JxATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//           nE[(xi+1)*Ny*Nz+yi*Nz+zi] += JxE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//           nE[xi*Ny*Nz+yi*Nz+zi] -= JxE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//         }
//         if (insideArr[xi*Ny*Nz+(yi+1)*Nz+zi] && insideArr[xi*Ny*Nz+yi*Nz+zi]){
//           nADP[xi*Ny*Nz+(yi+1)*Nz+zi] += JyADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//           nADP[xi*Ny*Nz+yi*Nz+zi] -= JyADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//           nATP[xi*Ny*Nz+(yi+1)*Nz+zi] += JyATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//           nATP[xi*Ny*Nz+yi*Nz+zi] -= JyATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//           nE[xi*Ny*Nz+(yi+1)*Nz+zi] += JyE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//           nE[xi*Ny*Nz+yi*Nz+zi] -= JyE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//         }
//         if (insideArr[xi*Ny*Nz+yi*Nz+(zi+1)] && insideArr[xi*Ny*Nz+yi*Nz+zi]){
//           nADP[xi*Ny*Nz+yi*Nz+(zi+1)] += JzADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//           nADP[xi*Ny*Nz+yi*Nz+zi] -= JzADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//           nATP[xi*Ny*Nz+yi*Nz+(zi+1)] += JzATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//           nATP[xi*Ny*Nz+yi*Nz+zi] -= JzATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//           nE[xi*Ny*Nz+yi*Nz+(zi+1)] += JzE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//           nE[xi*Ny*Nz+yi*Nz+zi] -= JzE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
//         }
//       }
//     }
//   }
//   double ADP_to_ATP;
//   double de_to_ADP_E;
//   double ATP_to_d;
//   double E_d_to_de;
//   for(int xi=0;xi<Nx;xi++){
//     for(int yi=0;yi<Ny;yi++){
//       for(int zi=0;zi<Nz;zi++){
//         ADP_to_ATP = rate_ADP_ATP*nADP[xi*Ny*Nz+yi*Nz+zi]*time_step;
//         de_to_ADP_E = rate_de*NDE[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*time_step;
//         ATP_to_d = (rate_D + rate_dD*(ND[xi*Ny*Nz+yi*Nz+zi] + NDE[xi*Ny*Nz+yi*Nz+zi])/mem_A[xi*Ny*Nz+yi*Nz+zi])
//           *nATP[xi*Ny*Nz+yi*Nz+zi]*time_step;
//         E_d_to_de = rate_E*ND[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*nE[xi*Ny*Nz+yi*Nz+zi]*time_step;
//         //Jeff!  remember that when you gain cyto density and lose the same amount of wall density,
//         //the numbers of proteins gained/lost will be different, and it's the numbers that you want to be the same!!
//         //also, keep thinking about the issue below where all additions are divided and then mult by the same mem_A fun
//         nADP[xi*Ny*Nz+yi*Nz+zi] -= ADP_to_ATP;
//         nATP[xi*Ny*Nz+yi*Nz+zi] += ADP_to_ATP;
//         if (mem_A[xi*Ny*Nz+yi*Nz+zi] != 0){
//           NDE[xi*Ny*Nz+yi*Nz+zi] += -de_to_ADP_E*mem_A[xi*Ny*Nz+yi*Nz+zi];
//           nADP[xi*Ny*Nz+yi*Nz+zi] += de_to_ADP_E*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);
//           nE[xi*Ny*Nz+yi*Nz+zi] += de_to_ADP_E*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);

//           nATP[xi*Ny*Nz+yi*Nz+zi] -= ATP_to_d*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);
//           ND[xi*Ny*Nz+yi*Nz+zi] += ATP_to_d*mem_A[xi*Ny*Nz+yi*Nz+zi];

//           nE[xi*Ny*Nz+yi*Nz+zi] -= E_d_to_de*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);
//           ND[xi*Ny*Nz+yi*Nz+zi] -= E_d_to_de*mem_A[xi*Ny*Nz+yi*Nz+zi];
//           NDE[xi*Ny*Nz+yi*Nz+zi] += E_d_to_de*mem_A[xi*Ny*Nz+yi*Nz+zi];
//         }
//         NflD[xi*Ny*Nz+yi*Nz+zi] = (nATP[xi*Ny*Nz+yi*Nz+zi] + nADP[xi*Ny*Nz+yi*Nz+zi])*(dx*dx*dx) + ND[xi*Ny*Nz+yi*Nz+zi] + NDE[xi*Ny*Nz+yi*Nz+zi];
//         NflE[xi*Ny*Nz+yi*Nz+zi] = nE[xi*Ny*Nz+yi*Nz+zi]*dx*dx*dx + NDE[xi*Ny*Nz+yi*Nz+zi];
//       }
//     }
//   }
//   return 0;
// }

// double delta_t = find_delta_t_for_next_event;
//         get_next_state(mem_A, insideArr, nATP, nADP, nE, ND, NDE, NflD, NflE, JxATP, JyATP, JzATP,
//                        JxADP, JyADP, JzADP, JxE, JyE, JzE);

// get_next_density(mem_A, insideArr, nATP, nADP, nE, ND, NDE, NflD, NflE, JxATP, JyATP, JzATP,
//                        JxADP, JyADP, JzADP, JxE, JyE, JzE);

