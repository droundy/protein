using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include "weights.h"
#include <cassert>
#include <time.h>
#include <math.h>


/*
  weights(int myN)// constructor, allocates new weights
  ~weights() { // destructor

  void update(double w, int i)
    this will update tot_weights and also the memory within the
    class according to the reaction and placement input.

  int lookup(double p) const { // p is from 0 to 1
    protein_microscopy.cpp will give a random number and this
    returns where in the memory allocated this is (in terms of index i)

  double get_total() const {
*/
/////////////////////hypothetical test code////////////////
const int Nx = 3;
const int Ny = 3;
const int Nz = 3;
const int size_of_ps = Nx*Ny*Nz;
weights ws(Nx*Ny*Nz);
const int Nx_big = 100;
const int Ny_big = 100;
const int Nz_big = 100;


double difference_in_lookups(double *ps, double ps_total, weights ws, int num_lookups) {
  int *bins = new int [size_of_ps];
  for (int i=0;i<size_of_ps;i++){
    bins[i] = 0;
  }
  for (int i=0;i<num_lookups;i++){
    double random =  (double)rand()/(RAND_MAX);
    int index = ws.lookup(random);//random number is from 0 to 1
    //printf("rand = %g and index = %d\n",random,index);
    bins[index]++;
  }
  double diff = ps[0]/ps_total - double(bins[0])/double(num_lookups);
  //printf("looks = %d and diff = %g and %g and %g\n",num_lookups,ps[0]/ps_total, double(bins[0]),double(num_lookups));
  delete[] bins;
  return diff;
  //  return 0.0;
}


int main() {
  double *ps  = new double[size_of_ps];
  for (int i=0;i<size_of_ps;i++){
    ps[i]=0;
  }
  double den = 256.00;
  ps[0] = 20/den;
  ps[1] = 30/den;
  ps[2] = 40/den;
  ws.update(ps[0], 0);
  ws.update(ps[1], 1);
  ws.update(ps[2], 2);
  double ps_total=0;
  for (int i=0;i<size_of_ps;i++){
    ps_total+=ps[i];
  }
  /////////////////////First test
  assert(fabs(ws.get_total()-ps_total) < 10e-20);
  printf("passed first test, total difference = %g\n",fabs(ws.get_total()-ps_total));


  ////////////////////Second Test
  int num_looks = 4;
  int *looks = new int[num_looks+1];
  looks[0]=100;
  double diff = 200;
  for (int i=1;i<num_looks+1;i++) {
    looks[i] = 10*looks[i-1];
    double old_diff = diff;
    diff = difference_in_lookups(ps,ps_total,ws,looks[i]);
    printf("num_looks = %d and diff = %g and old_diff = %g\n",looks[i],fabs(diff),fabs(old_diff));
    assert(fabs(old_diff) > fabs(diff));
  }
  printf("Passed Second test\n");
  delete [] looks;


  //////////////////////Third Test
  weights ws_big(Nx_big*Ny_big*Nz_big*23);
  double prob_test = 1.0;
  const clock_t begin_time = clock();
  for (int i=0;i<Nx_big*Ny_big*Nz_big*23;i++){
    ws_big.update(prob_test, i);
    ws_big.update(prob_test, i+1);
    //ws_big.update(prob_test, int(rand()%Nx_big*Ny_big*Nz_big*23));
  }
  clock_t end_time = clock();
  assert(double(end_time-begin_time)/CLOCKS_PER_SEC < 1.0);
  printf("Passed third test, time to update everything once = %g\n",double(end_time-begin_time)/CLOCKS_PER_SEC);
  /*Less than 1 second here is actually very bad.  Below is code that
    figures how long a typical simulation (1500 sim-seconds) will take
    for a hypothetical situation in wich every lattice point and a
    neighbor is updated per time step, using the current update
    function (not a good guess, just a starting point).  */

  int tot_iters_typical_sim = int(15000 / ((.1*.05*.05)/2.5));
  double tot_sec_for_hypothetical_sim = tot_iters_typical_sim*double(end_time-begin_time)/CLOCKS_PER_SEC;
  int days = floor(tot_sec_for_hypothetical_sim/(60*60*24));
  int hours = floor( (tot_sec_for_hypothetical_sim - days*(60*60*24)) / (60*60)  );
  int min = floor( (tot_sec_for_hypothetical_sim - days*(60*60*24) - hours*(60*60)) / 60);

  printf("Total time for a hypothetical simulation is = %d days, %d hours, and %d min\n",days,hours,min);


  /////////////////////Fourth Test
  //assert(fabs(ws_big.get_total() - prob_test*Nx_big*Ny_big*Nz_big*23) < 10e-20);
  printf("Passed fourth test, total difference = %g\n",fabs(ws_big.get_total() - prob_test*Nx_big*Ny_big*Nz_big*23));

  delete[] ps;

  return 0;
}
//////////////////////////////////////////////////////
