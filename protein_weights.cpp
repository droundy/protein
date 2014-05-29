#include "protein.h"




void update_densities_and_weighting_for_reaction(stoch_params p, weights *ws, bool *insideArr, int *s_N_ATP,
                                                        int *s_N_ADP, int *s_N_E, int *s_ND, int *s_NDE, double *mem_A){
  /* reminder of the order of the reaction enums and what d is:
     enum reaction {ADP_to_ATP, DE_to_ADP_E, ATP_to_D, E_D_to_DE,
                X_ADP_pos, X_ADP_neg, Y_ADP_pos, Y_ADP_neg, Z_ADP_pos, Z_ADP_neg,
                X_ATP_pos, X_ATP_neg, Y_ATP_pos, Y_ATP_neg, Z_ATP_pos, Z_ATP_neg,
                X_E_pos, X_E_neg, Y_E_pos, Y_E_neg, Z_E_pos, Z_E_neg};
     const int d [6][3] = { {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1} };
     Also jz = -difD*(n[z+1]-n[z])/dx and jz is #/(dA*dt) and I want dN/dt units for our probabilities
     so dN/dt = difD*(n[z+1]-n[z])*dA/dx = difD*(N[z+1]-N[z])/(dA)
  */
  if (!insideArr[p.xi*Ny*Nz + p.yi*Nz + p.zi] || p.xi == 0 || p.yi == 0 || p.zi == 0) {
    return;
  }
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
    ws->update(rate_ADP_ATP*s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi], ADP_to_ATP*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    //This loop goes around the changed lattice point and updates the diffusion propabilities for going to each adjacent lattice point
    for (int i=0;i<6;i++){
      if (insideArr[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ]) {
        if (s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ADP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
          ws->update( difD*( s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ADP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] ) / dA,
                      (X_ADP_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
        else {
          ws->update( 0.0, (X_ADP_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
      }
    }
  }
  if (p.reaction == ADP_to_ATP || p.reaction == ATP_to_D) { //s_N_ATP has changed
    if (mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] != 0.0) {
      ws->update((rate_D*mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] + rate_dD*(s_ND[p.xi*Ny*Nz+p.yi*Nz+p.zi] + s_NDE[p.xi*Ny*Nz+p.yi*Nz+p.zi]))
                *s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi]/dV, ATP_to_D*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
    for (int i=0;i<6;i++){
      if (insideArr[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ]) {
        if ( s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ATP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
          ws->update( difD*( s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ATP[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ] ) / dA,
                      (X_ATP_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
        else {
          ws->update( 0.0, (X_ATP_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
      }
    }
  }
  if (p.reaction == DE_to_ADP_E || p.reaction == E_D_to_DE) { //s_N_E has changed
    if (mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] != 0.0) {
      ws->update(rate_E*s_ND[p.xi*Ny*Nz+p.yi*Nz+p.zi]*s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi]/dV, E_D_to_DE*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
    for (int i=0;i<6;i++){
      if (insideArr[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ]) {
        if ( s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_E[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
          ws->update( difD*( s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_E[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ] ) / dA,
                      (X_E_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
        else {
          ws->update( 0.0, (X_E_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
      }
    }
  }
  if (p.reaction == ATP_to_D || p.reaction == E_D_to_DE) { //s_ND has changed
    if (mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] != 0.0) {
      ws->update((rate_D*mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] + rate_dD*(s_ND[p.xi*Ny*Nz+p.yi*Nz+p.zi] + s_NDE[p.xi*Ny*Nz+p.yi*Nz+p.zi]))
                *s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi]/dV, ATP_to_D*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      ws->update(rate_E*s_ND[p.xi*Ny*Nz+p.yi*Nz+p.zi]*s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi]/dV, E_D_to_DE*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
  }
  if (p.reaction == DE_to_ADP_E || p.reaction == E_D_to_DE) { //s_NDE has changed
    if (mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] != 0.0) {
      ws->update(rate_de*s_NDE[p.xi*Ny*Nz+p.yi*Nz+p.zi], DE_to_ADP_E*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
      ws->update((rate_D*mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] + rate_dD*(s_ND[p.xi*Ny*Nz+p.yi*Nz+p.zi] + s_NDE[p.xi*Ny*Nz+p.yi*Nz+p.zi]))
                *s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi]/dV, ATP_to_D*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
  }
  return;
}





void update_densities_and_weighting_for_diffusion(stoch_params p, weights *ws, bool *insideArr, int *s_N_ATP,
                                                         int *s_N_ADP, int *s_N_E, int *s_ND, int *s_NDE, double *mem_A){
  /* reminder of the order of the reaction enums and what d is:
     enum reaction {ADP_to_ATP, DE_to_ADP_E, ATP_to_D, E_D_to_DE,
                X_ADP_pos, X_ADP_neg, Y_ADP_pos, Y_ADP_neg, Z_ADP_pos, Z_ADP_neg,
                X_ATP_pos, X_ATP_neg, Y_ATP_pos, Y_ATP_neg, Z_ATP_pos, Z_ATP_neg,
                X_E_pos, X_E_neg, Y_E_pos, Y_E_neg, Z_E_pos, Z_E_neg};
     const int d [6][3] = { {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1} };
     Also jz = -difD*(n[z+1]-n[z])/dx and jz is #/(dA*dt) and I want dN/dt units for our probabilities
     so dN/dt = difD*(n[z+1]-n[z])*dx  */
  if (!insideArr[p.xi*Ny*Nz + p.yi*Nz + p.zi] || p.xi < 0.5 || p.yi < 0.5 || p.zi < 0.5) {
    return;
  }
  double dA = dx*dx;
  double dV = dx*dx*dx;
  if (p.reaction >= X_ADP_pos && p.reaction <= Z_ADP_neg) { //s_N_ADP has diffused
    for (int i=0;i<6;i++){
      if (p.reaction == (X_ADP_pos + i) && !insideArr[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ] ) {
        return;
      }
    }
    s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] -= 1;
    ws->update(rate_ADP_ATP*s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi], ADP_to_ATP*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    //This loop goes around the changed lattice point and updates the diffusion propabilities for going to each adjacent lattice point
    for (int i=0;i<6;i++){
      if (insideArr[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ]) {
        if ( s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ADP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
          ws->update( difD*( s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ADP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] ) / dA,
                      (X_ADP_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
        else {
          ws->update( 0.0, (X_ADP_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
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
    ws->update(rate_ADP_ATP*s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi], ADP_to_ATP*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    for (int i=0;i<6;i++){
      if (insideArr[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ]) {
        if ( s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ADP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
          ws->update( difD*( s_N_ADP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ADP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] ) / dA,
                      (X_ADP_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
        else {
          ws->update( 0.0, (X_ADP_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
      }
    }
    return;
  }
  else if (p.reaction >= X_ATP_pos && p.reaction <= Z_ATP_neg) { //s_N_ATP has diffused
    for (int i=0;i<6;i++){
      if (p.reaction == (X_ATP_pos + i) && !insideArr[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ]) {
        return;
      }
    }
    s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] -= 1;
    if (mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] != 0.0) {
      ws->update((rate_D*mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] + rate_dD*(s_ND[p.xi*Ny*Nz+p.yi*Nz+p.zi] + s_NDE[p.xi*Ny*Nz+p.yi*Nz+p.zi]))
                *s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi]/dV, ATP_to_D*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
    for (int i=0;i<6;i++){
      if (insideArr[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ]) {
        if ( s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ATP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
          ws->update( difD*( s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ATP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] ) / dA,
                      (X_ATP_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
        else {
          ws->update( 0.0, (X_ATP_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
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
      ws->update((rate_D*mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] + rate_dD*(s_ND[p.xi*Ny*Nz+p.yi*Nz+p.zi] + s_NDE[p.xi*Ny*Nz+p.yi*Nz+p.zi]))
                *s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi]/dV, ATP_to_D*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
    for (int i=0;i<6;i++){
      if (insideArr[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ]) {
        if ( s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ATP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
          ws->update( difD*( s_N_ATP[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_ATP[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] ) / dA,
                      (X_ATP_pos + 1)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
        else {
          ws->update( 0.0, (X_ATP_pos + 1)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
      }
    }
    return;
  }
  else if (p.reaction >= X_E_pos && p.reaction <= Z_E_neg) { //s_N_E has diffused
    for (int i=0;i<6;i++){
      if (p.reaction == (X_E_pos + i) && !insideArr[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ]) {
        return;
      }
    }
    s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] -= 1;
    if (mem_A[p.xi*Ny*Nz+p.yi*Nz+p.zi] != 0.0) {
      ws->update(rate_E*s_ND[p.xi*Ny*Nz+p.yi*Nz+p.zi]*s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi]/dV, E_D_to_DE*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
    for (int i=0;i<6;i++){
      if (insideArr[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ]) {
        if ( s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_E[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
          ws->update( difD*( s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_E[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] ) / dA,
                      (X_E_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
        else {
          ws->update( 0.0, (X_E_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
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
      ws->update(rate_E*s_ND[p.xi*Ny*Nz+p.yi*Nz+p.zi]*s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi]/dV, E_D_to_DE*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
    }
    for (int i=0;i<6;i++){
      // printf("Third p.reaction = %d p.xi = %d p.yi = %d p.zi = %d i = %d\n",p.reaction,p.xi,p.yi,p.zi,i);
      // fflush(stdout);
      if (insideArr[ (p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2]) ]) {
        if ( s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_E[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] >= 0) {
          ws->update( difD*( s_N_E[p.xi*Ny*Nz+p.yi*Nz+p.zi] - s_N_E[(p.xi+d[i][0])*Ny*Nz + (p.yi+d[i][1])*Nz + (p.zi+d[i][2])] ) / dA,
                      (X_E_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
        else {
          ws->update( 0.0, (X_E_pos + i)*Nx*Ny*Nz + p.xi*Ny*Nz+p.yi*Nz+p.zi);
        }
      }
      // printf("Fourth p.xi = %d p.yi = %d p.zi = %d\n",p.xi,p.yi,p.zi);
      // fflush(stdout);
    }
    return;
  }
}





void initialize_densities_and_weighting(weights *ws,  bool *insideArr, int *s_N_ATP, int *s_N_ADP, int *s_N_E,
                                        int *s_ND, int *s_NDE, double *mem_A){
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
  printf("We have allocated and initialized the probability weighting memory\n");
}







void update_all_densities_and_weighting_for_changing_gridpt(weights *ws, bool *insideArr, int *s_N_ATP,
                                                   int *s_N_ADP, int *s_N_E, int *s_ND, int *s_NDE, double *mem_A,
                                                   int xi, int yi, int zi){
 /* reminder of the order of the reaction enums and what d is:
     enum reaction {ADP_to_ATP, DE_to_ADP_E, ATP_to_D, E_D_to_DE,
                X_ADP_pos, X_ADP_neg, Y_ADP_pos, Y_ADP_neg, Z_ADP_pos, Z_ADP_neg,
                X_ATP_pos, X_ATP_neg, Y_ATP_pos, Y_ATP_neg, Z_ATP_pos, Z_ATP_neg,
                X_E_pos, X_E_neg, Y_E_pos, Y_E_neg, Z_E_pos, Z_E_neg};
     const int d [6][3] = { {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1} };
     Also jz = -difD*(n[z+1]-n[z])/dx and jz is #/(dA*dt) and I want dN/dt units for our probabilities
     so dN/dt = difD*(n[z+1]-n[z])*dA/dx = difD*(N[z+1]-N[z])/(dA)
  */
  if (!insideArr[xi*Ny*Nz + yi*Nz + zi]) {
    printf("Error!  You're trying to update probs for a gird point that is completely outside the cell.  There shouldn't be a reason to do that.\n");
    exit(1);
  }
  double dV = dx*dx*dx;
  double dA = dx*dx;
  ws->update(rate_ADP_ATP*s_N_ADP[xi*Ny*Nz+yi*Nz+zi], ADP_to_ATP*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
  //This loop goes around the changed lattice point and updates the diffusion propabilities for going to each adjacent lattice point
  if (mem_A[xi*Ny*Nz+yi*Nz+zi] != 0.0) {
    ws->update((rate_D*mem_A[xi*Ny*Nz+yi*Nz+zi] + rate_dD*(s_ND[xi*Ny*Nz+yi*Nz+zi] + s_NDE[xi*Ny*Nz+yi*Nz+zi]))
               *s_N_ATP[xi*Ny*Nz+yi*Nz+zi]/dV, ATP_to_D*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
    ws->update(rate_E*s_ND[xi*Ny*Nz+yi*Nz+zi]*s_N_E[xi*Ny*Nz+yi*Nz+zi]/dV, E_D_to_DE*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
    ws->update(rate_de*s_NDE[xi*Ny*Nz+yi*Nz+zi], DE_to_ADP_E*Nx*Ny*Nz + xi*Ny*Nz+yi*Nz+zi);
  }
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
  return;
}

