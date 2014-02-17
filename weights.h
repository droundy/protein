
class weights {
 public:
  weights(int myN) { // constructor, allocates new weights
    N = myN;
    tot_weights = 0;
    ws = new double[N];
  }
  ~weights() { // destructor
    delete[] ws;
  }
  void update(int i, double w) {
    tot_weights += w - ws[i];
    ws[i] = w;
  }
  int lookup(double p) const { // p is from 0 to 1
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
