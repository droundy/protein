#include "weights.h"
#include <math.h>


weights::weights(int num_of_possibilities, char input_sim_type) { // constructor, allocates new weights
  //will want to allocate and initialize all the memory here
  sim_type = input_sim_type;
  if (sim_type == 'h') {
    num_on_top_row = 1;
    size_of_array = 2;
    while (num_of_possibilities > num_on_top_row*2) {
      num_on_top_row *=2;
      size_of_array += num_on_top_row;
    }
    ws = new double[size_of_array];
    printf("size of array is = %d\n",size_of_array);
    for (int i=0;i<size_of_array;i++){
      ws[i] = 0.0;
    }
    num_below = size_of_array/2;
  }
  else {//if (sim_type == 'f') {
    num_on_top_row = 2;
    size_of_array = 2;
    while (num_of_possibilities > num_on_top_row) {
      num_on_top_row *=2;
      size_of_array += num_on_top_row;
    }
    ws = new double[size_of_array];
    num_below = size_of_array/2 - 1;
    for (int i=0;i<size_of_array;i++){
      ws[i] = 0.0;
    }
  }
}

void weights::update(double w, int i) {
  /*this will update tot_weights and also the memory within the
    class according to the reaction and placement input.*/
  if (sim_type == 'h') {
    int index = i;
    double w_old;
    if (index%2 == 0) {
      w_old = ws[index/2 + num_below];
    }
    else {
      index = index/2 + num_below;
      while (index%2 != 0 && index>0) {
        index = index/2;
      }
      index = index/2;
      if (index==0) {
        w_old = ws[0] - ws[1];
        index++;
      }
      else {
        w_old = ws[index] - ws[index*2];
        index *= 2;
      }
      while (index<num_below) {
        index = index*2 + 1;
        w_old = w_old - ws[index];
      }
    }
    double delta_w = w - w_old;
    ws[0] += delta_w;
    index = 1;
    int N_row = num_on_top_row;
    while (index < size_of_array) {
      if (i >= N_row) {
        index = index*2 + 1;
        i -= N_row;
      }
      else {
        ws[index] += delta_w;
        index *= 2;
      }
      N_row = N_row/2;
    }
    tot_weights += delta_w;
    return;
  }
  else {//if (sim_type == 'f') {
    int index = i + num_below;
    double delta_w = w - ws[index];
    ws[index] += delta_w;
    while (index > 1) {
      index = index/2 - 1;
      ws[index] += delta_w;
    }
  }
}




int weights::lookup(double p) const { // p is from 0 to 1
  /*protein_microscopy.cpp will give a random number and this
    returns where in the memory allocated this is (in terms of index i)*/
  if (sim_type == 'h') {
    p *= ws[0];
    int index = 1;
    while (index < num_below) {
      if (p > ws[index]) {
        p -= ws[index];
        index = index*2 + 1;
      }
      else {
        index *= 2;
      }
    }
    if (p > ws[index]) {
      return (index - num_below)*2 + 1;
    }
    else {
      return (index - num_below)*2;
    }
  }
  else {//if (sim_type == 'f') {
    p *= (ws[0] + ws[1]);
    int index = 0;
    while (index < num_below) {
      if (p > ws[index]) {
        p -= ws[index];
        index = (index+2)*2;
      }
      else {
        index = (index+1)*2;
      }
    }
    if (p > ws[index]) {
      return index + 1 - num_below;
    }
    else {
      return index - num_below;
    }
  }
}


double weights::get_total() {
  if (sim_type == 'h') {
    return ws[0];
  }
  else {//if (sim_type == 'f') {
    return (ws[0] + ws[1]);
  }
}

double weights::lookup_prob_for_specific_index(int index) {
  if (sim_type == 'h') {
    if (index%2 == 0) {
      return ws[index/2 + num_below];
    }
    else {
      index = index/2 + num_below;
      while (index%2 != 0 && index>0) {
        index = index/2;
      }
      index = index/2;
      double prob = ws[index];
      if (index != 0) {
        index *= 2;
        prob -= ws[index];
      }
      while (index < num_below) {
        index = index*2 + 1;
        prob -= ws[index];
      }
      return prob;
    }
  }
  else {//if (sim_type == 'f') {
    return ws[index + num_below];
  }
}
