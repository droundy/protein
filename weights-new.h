#ifndef __WEIGHTS_H_INCLUDED__
#define __WEIGHTS_H_INCLUDED__

#include <stdio.h>
#include <math.h>

class weights {
 public:
  weights()= {}
  virtual void update(double w, int i) = 0;
  int lookup(double p) const = 0;
  double get_total() = 0;
  double lookup_prob_for_specific_index(int index) = 0;

  /*Testing Functions*/
  void update_one_spot_and_look_at_entire_array(double w, int index) = 0;
  int test_lookup_from_zero_to_one(int num_possibilities) = 0;
  int test_get_total_matches_summing_up(int num_possibilities) = 0;
  int test_looking_up_shows_correct_probs(int num_possibilities) = 0;
  int test_time_of_simulation(int num_possibilities) = 0;
  int test_of_lookup_prob_for_specific_index(int num_possibilities) = 0;
};

class fullweights {
 public:
  fullweights()= {}
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
  data foo;
};

class halfweights {
 public:
  halfweights()= {}
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
};

#endif
