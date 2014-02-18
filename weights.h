

class weights {
 public:
  weights(int myN) { // constructor, allocates new weights
    //will want to allocate and initialize all the memory here
    N = myN;
    tot_weights = 0;
    ws = new double[N];
  }
  ~weights() { // destructor
    delete[] ws;
  }
  void update(double w, int i) {
    /*this will update tot_weights and also the memory within the
      class according to the reaction and placement input.*/
    double tot_weights += w - ws[i];
    ws[i] = w;
  }
  int lookup(double p) const { // p is from 0 to 1
    /*protein_microscopy.cpp will give a random number and this
    returns where in the memory allocated this is (in terms of index i)*/
    p *= tot_weights;
    for (int i=0; i<N; i++) {
      if (p < ws[i]) {
        return i;
      }
      p -= ws[i];
    }
    return N-1;
  }
  double get_total() const {
    return tot_weights;
  }
 private:
  int N;
  double tot_weights;
  double *ws; // or something smarter
};
