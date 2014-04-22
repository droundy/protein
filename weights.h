#ifndef __WEIGHTS_H_INCLUDED__
#define __WEIGHTS_H_INCLUDED__

#include <stdio.h>
#include <math.h>

class weights {
 public:
  weights(int num_of_possibilities);
  ~weights() { // destructor
    delete[] ws;
  }
  void update(double w, int i);
  int lookup(double p) const;
  double get_total() const {
    return ws[0];
  }
  double lookup_prob_for_specific_index(int index);

  /*Testing Functions*/
  void update_one_spot_and_look_at_entire_array(double w, int index);
  int test_lookup_from_zero_to_one(int num_possibilities);
  int test_get_total_matches_summing_up(int num_possibilities);
  int test_looking_up_shows_correct_probs(int num_possibilities);
  int test_time_of_simulation(int num_possibilities);
  int test_of_lookup_prob_for_specific_index(int num_possibilities);

 private:
  int size_of_array;
  double *ws;
  int num_on_top_row;
  int num_below;
  double tot_weights;
};

#endif
