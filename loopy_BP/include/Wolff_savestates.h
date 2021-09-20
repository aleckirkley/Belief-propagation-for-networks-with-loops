#ifndef WOLFF_HEADER 
#define WOLFF_HEADER

#include "Graph.h"
#include <math.h>
#include <random>
#include <deque>

std::random_device RD;
std::uniform_real_distribution<double> unif(0.0,1.0);
std::mt19937 RNG_Wolff(RD());

class WolffSampler {
  public:
    std::vector<int> S;
    std::vector<std::vector<int>> states;
    std::vector<int> H;
    Graph G;
    double p;
    void n_cluster_flips(int);
    void sample_states(int, double, int);
    void reorient();
    WolffSampler(Graph);
    WolffSampler(){};
};


WolffSampler::WolffSampler(Graph G_in):
  S(G_in.number_of_nodes(),1)
{
  G = G_in;
  p = -1;  //initialize with invalid value 
}


void WolffSampler::n_cluster_flips(int n) {
  std::uniform_int_distribution<int> unif_node(0, G.number_of_nodes()-1);
  for (int k=0; k<G.number_of_nodes(); ++k) {
      S[k] = (unif(RNG_Wolff) < 1.0) ? 1.0 : -1.0;
  } 
  for (int iteration=0; iteration<n; ++iteration) {
    int k = unif_node(RNG_Wolff);
    int spin = S[k];
    S[k] = -spin;
    std::deque<int> to_explore = {k};
    while (to_explore.size()>0) {
      int j = to_explore.front();
      to_explore.pop_front();
      for (int l : G.neighbors(j)){
        if (S[l]==spin) {
            if (unif(RNG_Wolff)<p){
                to_explore.push_back(l);
                S[l] = -spin;
            }
          
        }
      }
    }
    reorient(); // make sure node 0 is always S[0]=1;
  }
}


// state sampler -- only sample if value of p has changed
void WolffSampler::sample_states(int n, double beta, int period) {
  if (p!=1-exp(-beta)){
    p = 1 - exp(-beta);
    n_cluster_flips(n*period);
    states.resize(n);
    H.resize(n);
    for (int x=0; x<n; x++){
      n_cluster_flips(period);
      states[x] = S;
      H[x] = 0;
      for (int k : G.nodes())
	for (int j : G.neighbors(k))
	  if (k<j) H[x] += (S[k]!=S[j]);
    }
  }
}


void WolffSampler::reorient() {
  if (S[0] != 1)
    for ( int i=0; i<S.size(); ++i )
      S[i] = -S[i];
}


#endif
