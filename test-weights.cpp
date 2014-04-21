#include <cfloat>
#include <stdio.h>
#include <stdlib.h>
#include "weights.h"
#include <cassert>
#include <time.h>
#include <math.h>

const int Nx_big = 100;
const int Ny_big = 100;
const int Nz_big = 100;

int main() {
    /*At this point the two main issues are how long things take and
      how much is lost in terms of updating the very large arrays
      (little errors adding up to big ones)*/

  int num_possibilities = 16;
  weights ws1 = weights(num_possibilities);
  for (int i=0;i<num_possibilities;i++) {
    ws1.update(1.0,i);
  }
  ws1.update_one_spot_and_look_at_entire_array(1.0, 12);

  int test_pass[4];
  test_pass[0] = ws1.test_lookup_from_zero_to_one(num_possibilities);
  printf("\ntesting lookup from zero to one = %s\n",test_pass[0]? "PASS" : "FAIL");

  //The following test should be tested on much larger arrays
  test_pass[1] = ws1.test_get_total_matches_summing_up(num_possibilities);
  printf("\ntesting get_total matches with summing up probs = %s\n",test_pass[1]? "PASS" : "FAIL");

  test_pass[2] = ws1.test_looking_up_shows_correct_probs(num_possibilities);
  printf("\ntesting whether looking up many times shows the right probabilities = %s\n\n",
         test_pass[2]? "PASS" : "FAIL");

  int big_num_possibilities = Nx_big*Ny_big*Nz_big*23;
  weights ws_big = weights(big_num_possibilities);
  test_pass[3] = ws_big.test_time_of_simulation(big_num_possibilities);
  printf("\ntesting the time it takes to simulate = %s\n\n",
         test_pass[3]? "PASS" : "FAIL");

  return 0;
}

//////////////////////////////////////////////////////


void weights::update_one_spot_and_look_at_entire_array(double w, int index) {
  for (int i=0;i<size_of_array;i++) {
    printf("ws[%d] = %g\n",i,ws[i]);
  }
  update(w, index);
  for (int i=0;i<size_of_array;i++) {
    printf("ws[%d] = %g\n",i,ws[i]);
  }
  return;
}


int weights::test_lookup_from_zero_to_one(int num_possibilities) {
  for (int i=0;i<num_possibilities;i++) {
    update(1.0,i);
  }
  int index = 0;
  double lookup_inc = 1.0/num_possibilities;
  double lookup_val = 0.5*(lookup_inc);
  for (int i=0;i<num_possibilities;i++) {
    //printf("lookup_inc = %g lookup_val = %g index = %d\n",lookup_inc,lookup_val,index);
    if (index != lookup(lookup_val)){
      return 0;
    }
    lookup_val += lookup_inc;
    index++;
  }
  return 1;
}


int weights::test_get_total_matches_summing_up(int num_possibilities) {
  if (num_possibilities < 8) {
    printf("\nnum_possibilities needs to be >= 8 for the get total matches test!!\n");
    return 0;
  }
  double total = 0;
  for (int i=0;i<num_possibilities;i++) {
    update(1.0,i);
    total += 1.0;
  }
  if (fabs(get_total() - total) > DBL_EPSILON){
    return 0;
  }
  update(1.0,1); total += 0.0;
  update(2.2,12); total += 1.2;
  update(3.6,14); total += 2.6;
  update(4.9,9); total += 3.9;
  update(5.2,7); total += 4.2;
  update(6.0,3); total += 5.0;
  if (fabs(get_total() - total) > 10e-15){//smaller than this it doesn't pass. Is that ok?
    return 0;
  }
  return 1;
}


int weights::test_looking_up_shows_correct_probs(int num_possibilities){
  if (num_possibilities < 8) {
    printf("\nnum_possibilities needs to be >= 8 for the shows correct probs test!!\n");
    return 0;
  }
  for (int i=0;i<num_possibilities;i++) {
    update(1.0,i);
  }
  int update_indexes[5] = {1,7,9,12,14}; //picked randomly from my head
  double update_values[5] = {1.0,2.2,3.6,4.9,5.2};//so are these
  int num_updates = 5;
  for (int i=0;i<num_updates;i++) {
    update(update_values[i],update_indexes[i]);
  }
  int num_lookup_tests = 3;
  int num_lookups = 1;
  int *bins = new int[num_possibilities];
  double* old_diffs = new double[num_updates];
  for (int i=0;i<num_updates;i++){
    old_diffs[i] = 123456.0; //just needed big number
  }
  for (int i=0;i<num_lookup_tests;i++) {
    printf("\nNew lookup test\n");
    num_lookups = 150*num_lookups;
    for (int j=0;j<num_possibilities;j++){
      bins[j] = 0;
    }
    for (int j=0;j<num_lookups;j++){
      bins[lookup((double)rand()/(RAND_MAX))]++; //random number is from 0 to
    }
    for (int k=0;k<num_updates;k++) {
      double diff = fabs(update_values[k]/get_total()
                         - double(bins[update_indexes[k]])/double(num_lookups));
      printf("update_value = %g, update_index = %d\n",update_values[k],update_indexes[k]);
      printf("looks = %d and diff = %g, old_diff = %g and %g and %g num_lookups = %d\n",num_lookups, diff,old_diffs[k],
             update_values[k]/get_total(),double(bins[update_indexes[k]])/double(num_lookups),num_lookups);
      if (fabs(old_diffs[k]) < fabs(diff)) {
        printf("test failed on the last printed\n\n");
        return 0;
      }
      old_diffs[k] = diff;
    }
  }
  delete[] bins;
  return 1;
}


int weights::test_time_of_simulation(int num_possibilities) {
    const clock_t begin_time = clock();
    for (int i=0;i<Nx_big*Ny_big*Nz_big*23-1;i++){
        update(1.0, i);
        update(1.0, i+1);
    }
    clock_t end_time = clock();
    printf("Time it update everything once = %g\n",
           double(end_time-begin_time)/CLOCKS_PER_SEC);
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
    if (double(end_time-begin_time)/CLOCKS_PER_SEC > 1.0) {
        return 0;
    }
    return 1;
}
