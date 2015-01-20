#include "protein.h"
#include <stdlib.h>


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

double f_2D_stad(double y, double z){
  double f;
  double Y = Ny*dx;
  double Z = Nz*dx;
  double z1 = (Z-B)/2.0;
  double z2 = z1 + B;
  double y1 = Y/2.0;
  f = sqrt((y-y1)*(y-y1))-C;
  if (z < z1) {
    f = sqrt((y-y1)*(y-y1)+(z-z1)*(z-z1))-C;
  }
  if (z > z2) {
    f = sqrt((y-y1)*(y-y1)+(z-z2)*(z-z2))-C;
  }
  return f;
}




//mem_f produces a scalar field on the grid. points where mem_f = 0 are cell wall.
double mem_f(double x, double y, double z) {
  if(mem_f_shape=="randst" || mem_f_shape=="TIE_fighter" || mem_f_shape=="triangle" || mem_f_shape=="stad"){
    double X = Nx*dx;
    double f_2d = 0;
    if(mem_f_shape=="triangle") f_2d = f_2D_triangle(y,z);
    else if(mem_f_shape=="TIE_fighter") f_2d = f_2D_TIE_fighter(y,z);
    else if(mem_f_shape=="randst") {
      f_2d = f_2D_randst(y,z);
      for (double yd = y-((A-.25)/2.0); yd <= y+((A-.25)/2.0); yd+=dx/4.0) {
        for (double zd = z-((A-.25)/2.0); zd <= z+((A-.25)/2.0); zd+=dx/4.0) {
          if (f_2D_randst(yd,zd) <= 0) f_2d = -1.0;
        }
      }
    }
    else if(mem_f_shape=="stad") {
      f_2d = f_2D_stad(y,z);
      for (double yd = y-((A-.25)/2.0); yd <= y+((A-.25)/2.0); yd+=dx/4.0) {
        for (double zd = z-((A-.25)/2.0); zd <= z+((A-.25)/2.0); zd+=dx/4.0) {
          if (f_2D_stad(yd,zd) <= 0) f_2d = -1.0;
        }
      }
    }
    else {
      printf("somethings wrong with the shape argument!!!");
      exit(1);
    }
    //    printf("\nX/2.0 = %g\n",X/2.0);
    if (f_2d<=0) {
      return fabs(2*(x-(X/2.0))/A) - 0.99999;//1.00001;
    }
    double closest_y0 = -100.0;
    double closest_z0 = -100.0;
    //closest pt is closest that is still inside 2d guassian boundaries of cell
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


void set_membrane(double mem_A[]) {
  for(int xi=0;xi<Nx;xi++){
    for(int yi=0;yi<Ny;yi++){
      for(int zi=0;zi<Nz;zi++){
        double fXYZ = mem_f((xi+0.5)*dx, (yi+0.5)*dx, (zi+0.5)*dx);
        double fXYz = mem_f((xi+0.5)*dx, (yi+0.5)*dx, (zi-0.5)*dx);
        double fXyZ = mem_f((xi+0.5)*dx, (yi-0.5)*dx, (zi+0.5)*dx);
        double fXyz = mem_f((xi+0.5)*dx, (yi-0.5)*dx, (zi-0.5)*dx);
        double fxYZ = mem_f((xi-0.5)*dx, (yi+0.5)*dx, (zi+0.5)*dx);
        double fxYz = mem_f((xi-0.5)*dx, (yi+0.5)*dx, (zi-0.5)*dx);
        double fxyZ = mem_f((xi-0.5)*dx, (yi-0.5)*dx, (zi+0.5)*dx);
        double fxyz = mem_f((xi-0.5)*dx, (yi-0.5)*dx, (zi-0.5)*dx);
        double f = mem_f(xi*dx, yi*dx, zi*dx);
        // if (xi == Nx-1 && (f < 0 || fXYZ < 0 || fXYz < 0 || fXyZ < 0 || fXyz < 0 || fxYZ < 0 || fxYz < 0 || fxyZ < 0 || fxyz < 0)) {
        //   printf("\nx %d,y %d,z %d and fXYZ = %g fxyz = %g f = %g",xi,yi,zi,fXYZ,fxyz,f);
        //   fflush(stdout);
        // }
        mem_A[xi*Ny*Nz+yi*Nz+zi] = find_intersection(fXYZ, fXYz, fXyZ, fXyz, fxYZ, fxYz, fxyZ, fxyz, f, false);
      }
    }
  }
}

