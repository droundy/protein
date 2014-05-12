#pragma once

#include <stdio.h>
#include <math.h>

class weights {
 public:
  weights(int num_of_possibilities, char sim_type); //sim_type either exact, 'f' for full_array, or 'h' for half_array
  ~weights() { // destructor
    delete[] ws;
  }
  void update(double w, int i);
  int lookup(double p) const;
  double get_total();
  double lookup_prob_for_specific_index(int index);

  /*Testing Functions*/
  void update_one_spot_and_look_at_entire_array(double w, int index);
  int test_lookup_from_zero_to_one(int num_possibilities);
  int test_get_total_matches_summing_up(int num_possibilities);
  int test_looking_up_shows_correct_probs(int num_possibilities);
  int test_time_of_simulation(int num_possibilities);
  int test_of_lookup_prob_for_specific_index(int num_possibilities);

 private:
  char sim_type;
  int size_of_array;
  double *ws;
  int num_on_top_row;
  int num_below;
  double tot_weights;
};

