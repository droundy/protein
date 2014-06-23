#include "protein.h"
#include <cassert>


void initialize_densities_and_weighting(weights *ws,  const bool *insideArr, int *s_N_ATP, int *s_N_ADP, int *s_N_E,
                                        int *s_ND, int *s_NDE, const double *mem_A){
  /* reminder of the order of the reaction enums and what d is:
     enum reaction {ADP_to_ATP, DE_to_ADP_E, ATP_to_D, E_D_to_DE,
                X_ADP_pos, X_ADP_neg, Y_ADP_pos, Y_ADP_neg, Z_ADP_pos, Z_ADP_neg,
                X_ATP_pos, X_ATP_neg, Y_ATP_pos, Y_ATP_neg, Z_ATP_pos, Z_ATP_neg,
                X_E_pos, X_E_neg, Y_E_pos, Y_E_neg, Z_E_pos, Z_E_neg};
     const int d [6][3] = { {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1} };
     Also jz = -difD*(n[z+1]-n[z])/dx and jz is #/(dA*dt) and I want dN/dt units for our probabilities
     so dN/dt = jz*dA = difD*(n[z+1]-n[z])*dA/dx = difD*(N[z+1]-N[z])/(dA)
  */
  // for (int i=0; i<num_pos_reactions*Nx*Ny*Nz;i++) {
  //   ws->update(0.0,i);
  // }
  double dV = dx*dx*dx;
  double dA = dx*dx;
  for(int xi=0;xi<Nx;xi++){
    for(int yi=0;yi<Ny;yi++){
      for(int zi=0;zi<Nz;zi++){
        if (insideArr[xi*Ny*Nz+yi*Nz+zi]) {
          ws->update(rate_ADP_ATP*s_N_ADP[xi*Ny*Nz+yi*Nz+zi], ADP_to_ATP*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);//1
          if (mem_A[xi*Ny*Nz+yi*Nz+zi] != 0.0) {
            ws->update((rate_D*mem_A[xi*Ny*Nz+yi*Nz+zi] + rate_dD*(s_ND[xi*Ny*Nz+yi*Nz+zi] + s_NDE[xi*Ny*Nz+yi*Nz+zi]))
                       *s_N_ATP[xi*Ny*Nz+yi*Nz+zi]/dV, ATP_to_D*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);//74
            ws->update(rate_E*s_ND[xi*Ny*Nz+yi*Nz+zi]*s_N_E[xi*Ny*Nz+yi*Nz+zi]/dV, E_D_to_DE*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);//74
            ws->update(rate_de*s_NDE[xi*Ny*Nz+yi*Nz+zi], DE_to_ADP_E*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);//.7
          }
          else {
            ws->update( 0.0, ATP_to_D*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
            ws->update( 0.0, E_D_to_DE*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
            ws->update( 0.0, DE_to_ADP_E*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
          }
          for (int i=0;i<6;i++){
            if (insideArr[ (xi+d[i][0])*Ny*Nz + (yi+d[i][1])*Nz + (zi+d[i][2]) ]) {
              ws->update( difD*s_N_ATP[xi*Ny*Nz+yi*Nz+zi] / dA,
                          (X_ATP_pos + i)*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);//1000
              ws->update( difD*s_N_ADP[xi*Ny*Nz+yi*Nz+zi] / dA,
                          (X_ADP_pos + i)*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
              ws->update( difD*s_N_E[xi*Ny*Nz+yi*Nz+zi] / dA,
                          (X_E_pos + i)*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
            }
          }
        }
      }
    }
  }
  printf("We have allocated and initialized the probability weighting memory\n");
  fflush(stdout);
}









void update_all_densities_and_weighting_for_changing_gridpt(weights *ws, const bool *insideArr, int *s_N_ATP,
                                                   int *s_N_ADP, int *s_N_E, int *s_ND, int *s_NDE, const double *mem_A,
                                                            int xi, int yi, int zi, char type){
 /* reminder of the order of the reaction enums and what d is:
     enum reaction {ADP_to_ATP, DE_to_ADP_E, ATP_to_D, E_D_to_DE,
                X_ADP_pos, X_ADP_neg, Y_ADP_pos, Y_ADP_neg, Z_ADP_pos, Z_ADP_neg,
                X_ATP_pos, X_ATP_neg, Y_ATP_pos, Y_ATP_neg, Z_ATP_pos, Z_ATP_neg,
                X_E_pos, X_E_neg, Y_E_pos, Y_E_neg, Z_E_pos, Z_E_neg};
     const int d [6][3] = { {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1} };
     Also jz = -difD*(n[z+1]-n[z])/dx and jz is #/(dA*dt) and I want dN/dt units for our probabilities
     so dN/dt = difD*(n[z+1]-n[z])*dA/dx = difD*(N[z+1]-N[z])/(dA)
 */
  assert(insideArr[xi*Ny*Nz + yi*Nz + zi]);
  // if (!insideArr[xi*Ny*Nz + yi*Nz + zi]) {
  //   printf("Error!  You're trying to update probs for a grid point that is completely outside the cell.  There shouldn't be a reason to do that.\n");
  //   printf("Grid pt is at x %d, y %d, z %d and type = %c\n",xi,yi,zi,type);
  //   exit(1);
  // }
  double dV = dx*dx*dx;
  double dA = dx*dx;
  ws->update(rate_ADP_ATP*s_N_ADP[xi*Ny*Nz+yi*Nz+zi], ADP_to_ATP*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
  if (mem_A[xi*Ny*Nz+yi*Nz+zi] != 0.0) {
    ws->update((rate_D*mem_A[xi*Ny*Nz+yi*Nz+zi] + rate_dD*(s_ND[xi*Ny*Nz+yi*Nz+zi] + s_NDE[xi*Ny*Nz+yi*Nz+zi]))
               *s_N_ATP[xi*Ny*Nz+yi*Nz+zi]/dV, ATP_to_D*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
    ws->update(rate_E*s_ND[xi*Ny*Nz+yi*Nz+zi]*s_N_E[xi*Ny*Nz+yi*Nz+zi]/dV, E_D_to_DE*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
    ws->update(rate_de*s_NDE[xi*Ny*Nz+yi*Nz+zi], DE_to_ADP_E*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
  }
  //These loops go around the changed lattice point and update the diffusion propabilities for going to each adjacent lattice point
  if (type == 'r'){
    for (int i=0;i<6;i++){
      if (insideArr[ (xi+d[i][0])*Ny*Nz + (yi+d[i][1])*Nz + (zi+d[i][2]) ]) {
        ws->update( difD*s_N_ADP[xi*Ny*Nz+yi*Nz+zi] / dA,
                    (X_ADP_pos + i)*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
        ws->update( difD*s_N_ATP[xi*Ny*Nz+yi*Nz+zi] / dA,
                    (X_ATP_pos + i)*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
        ws->update( difD*s_N_E[xi*Ny*Nz+yi*Nz+zi] / dA,
                    (X_E_pos + i)*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
      }
    }
  }
  else if (type == 'd'){
    for (int i=0;i<6;i++){
      if (insideArr[ (xi+d[i][0])*Ny*Nz + (yi+d[i][1])*Nz + (zi+d[i][2]) ]) {
        ws->update( difD*s_N_ADP[xi*Ny*Nz+yi*Nz+zi] / dA,
                    (X_ADP_pos + i)*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
      }
    }
  }
  else if (type == 't'){
    for (int i=0;i<6;i++){
      if (insideArr[ (xi+d[i][0])*Ny*Nz + (yi+d[i][1])*Nz + (zi+d[i][2]) ]) {
        ws->update( difD*s_N_ATP[xi*Ny*Nz+yi*Nz+zi] / dA,
                    (X_ATP_pos + i)*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
      }
    }
  }
  else if (type == 'e'){
    for (int i=0;i<6;i++){
      if (insideArr[ (xi+d[i][0])*Ny*Nz + (yi+d[i][1])*Nz + (zi+d[i][2]) ]) {
        ws->update( difD*s_N_E[xi*Ny*Nz+yi*Nz+zi] / dA,
                    (X_E_pos + i)*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
      }
    }
  }
  return;
}
