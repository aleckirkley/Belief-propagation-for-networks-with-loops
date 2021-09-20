#include <iostream>
#include <boost/math/differentiation/autodiff.hpp>
#include "Graph.h"
#include "mp_auto_MC.h"
#include "free_energy_approximations.h"

Graph G;
BeliefPropagation BP;

template <typename T>
std::vector<T> linspace(T start, T stop, int num) {
	std::vector<T> X(num);
	for (int i=0; i<num; ++i)
		X[i] = start + i*(stop-start)/(num-1);
	return X;
}

template <typename T>
bool n_iterations(T beta, double h, int n, double tol, int MC) {
	double current_C = (double)BP.C_total;
	for (int s=0; s<n; ++s){
		for (int i=0; i<5; ++i){
			BP.diff_iteration(beta,h,MC);
		}
		BP.iteration(beta,h,MC);
		if ( abs((double)BP.C_total-current_C)<tol )
			return true;
		current_C = (double)BP.C_total;

	}
	return false;
}

int main(int argc, const char *argv[]) {

	std::string fname = argv[1];
	int r = std::stoi(argv[2]);
	double T_min = std::stof(argv[3]);
	double T_max = std::stof(argv[4]);
	int n_steps = std::stoi(argv[5]);
	double h = std::stof(argv[6]);
	bool MC = std::stoi(argv[7]);
    int max_iters = std::stoi(argv[8]);
    double tol = std::stof(argv[9]);
    std::string outfname = argv[10];
    std::ofstream outfile;
	outfile.open(outfname);

	G = read_edgelist(fname);
	BP.init(G,r);
    
	for ( double T : linspace(T_min, T_max, n_steps) ) {
        
		a_Double beta=boost::math::differentiation::make_fvar<double, 4>(1./T);
        
        std::cout << n_iterations(beta,h,max_iters,tol,MC) << std::endl;
        double M = (double)BP.M_total;
        double U = (double)BP.U_total;
        double C = (double)BP.C_total;
        double S = (double)entropy(G, BP, 2.0*beta, h);
        
		std::cout << T << " " << M << " " << S << " " << C << std::endl;
		outfile << T << " " << M << " " << S << " " << C << std::endl;
        
		std::cout << std::endl;
	}

	return 0;
}

