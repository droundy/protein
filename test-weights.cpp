using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include "weights.h"
#include <cassert>



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
weights ws(Nx*Ny*Nz*23);

double difference_in_probs(double *ps, double ps_total, weights ws, int num_lookups) {
  int bins[9];
  for (int i;i<9;i++){
    bins[i] = 0;
  }
  for (int i=0;i<num_lookups;i++){
    int index = ws.lookup( (double)rand()/(RAND_MAX) );//random number is from 0 to 1
    bins[index]++;
  }
  double diff = ps[0]/ps_total - double(bins[0])/double(num_lookups);
  return diff;
}


int main() {
  double *ps  = new double[9];
  double den = 256.00;
  ps[0] = 20/den;
  ps[1] = 30/den;
  ps[2] = 40/den;
  double ps_total=0;
  for (int i=3;i<9;i++){
    ps[i]=0;
  }
  for (int i=0;i<9;i++){
    ps_total+=ps[i];
  }

  ws.update(ps[0], 0);
  ws.update(ps[1], 1);
  ws.update(ps[2], 2);

  int times=900000;
  double result = difference_in_probs(ps,ps_total,ws,times);
  printf("result = %g\n",result);


  printf("get_total diff = %g\n",ws.get_total()-ps[0]-ps[1]);
  assert(ws.get_total() == ps_total);

  printf("lookup(0.8) = %d\n",ws.lookup(0.8));
  assert(ws.lookup(0.8) == 2);

  return 0;
}
//////////////////////////////////////////////////////
