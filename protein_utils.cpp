#include "protein.h"


void trim_grid(double **pointer_to_mem_A, double *first_mem_A){
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
  printf("Here first\nNx=%d,Ny=%d,Nz=%d",Nx,Ny,Nz);
  printf("Here first\nnew_Nx=%d,newNy=%d,newNz=%d",new_Nx,new_Ny,new_Nz);
  printf("Here first\nmax_Nx=%d,maxNy=%d,maxNz=%d",max_xi,max_yi,max_zi);
  printf("Here first\nmin_Nx=%d,minNy=%d,minNz=%d",min_xi,min_yi,min_zi);
  fflush(stdout);
  *pointer_to_mem_A = new double[new_Nx*new_Ny*new_Nz];
  printf("Here next\n");
  fflush(stdout);
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
      if (fabs(mem_A[xi*Ny*Nz+ya*Nz+zi] - mem_A[xi*Ny*Nz+yb*Nz+zi]) > .000001) {
        fprintf(out_file,"xi = %d ya = %d zi = %d mem(xa) = %g\n",xi,yi,zi,mem_A[xi*Ny*Nz+ya*Nz+zi]);
        //fprintf(out_file,"\n");
      }
      //}
      //fprintf(out_file, "%d\t",int(inside(xi,yi,Nz/2)));
    }
    fprintf(out_file,"\n");
  }
  for (int i=0;i<10;i++){
    fclose(out_files[i]);
  }
  fclose(out_file);
  // printf("Done with sym check!!!!!");
  fflush(stdout);
}



bool only_once = true;
std::string triangle_section (double y, double z) {
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
  if (fabs(z_21-z3) < .000001){
    slope3 = 1000000*(y_21-y3);
  }
  if (fabs(z2-z_13) < .000001){
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


void test_the_amount_of_area(double *first_mem_A, std::string mem_f_shape){
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
  double percent_off = 100.0*fabs(A_ideal-A_ours)/A_ideal;
  if (A_ideal > A_ours) {
    printf("\nOur mem_A area is too small and is %g%% off\n\n",percent_off);
  }
  else {
    printf("\nOur mem_A area is too great and is %g%% off\n\n",percent_off);
  }
}




void set_insideArr(bool *insideArr){
  for(int xi=0;xi<Nx;xi++){
    for(int yi=0;yi<Ny;yi++){
      for(int zi=0;zi<Nz;zi++){
        insideArr[xi*Ny*Nz+yi*Nz+zi] = inside(xi,yi,zi);
        if (xi < 0.5 || yi < 0.5 || zi < 0.5){
          insideArr[xi*Ny*Nz+yi*Nz+zi] = false;
        }
      }
    }
  }
}

bool inside(int xi, int yi, int zi){
  xi = xi + min_xi-1;
  yi = yi + min_yi-1;
  zi = zi + min_zi-1;
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
