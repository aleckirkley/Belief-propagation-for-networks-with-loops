#pragma once
#include "mp_auto_MC.h"

template <typename T>
T distribution_entropy(const std::vector<T>& P) {
	T ans = 0;
	for ( int i=0; i<P.size(); ++i )
		ans -= P[i]*log0(P[i]);
	return ans;
}

template <typename T>
inline T mu_log(T mu_i, T mu_ij) {
	return mu_i*log0(mu_ij) + (1-mu_i)*log0(1-mu_ij);
}

template <typename T>
T energy(Graph& G, BeliefPropagation & BP, const T& beta, double h){
	T U = 0.0;
	for ( int i : G.nodes() ){
		U += h*BP.N_nodes[i].value;
		for (int j : G.neighbors(i)){
			auto P = BP.N_nodes[i].compute_two_point_distribution(beta, BP.N_diff, h, j);
			U += 0.5*(P[1]+P[2]);
		}
	}
	return U;
}

template <typename T>
T energy_I(Graph& G, BeliefPropagation & BP, const T& beta, double h){
	T U = 0.0;
	for ( int i : G.nodes() ){
		U += h*BP.N_nodes[i].value;
		for (int j : G.neighbors(i)){
			auto P = BP.N_intersect[i][j].compute_two_point_distribution(beta, BP.N_diff, 0, i, BP.N_diff[j][i].value); 
			U += 0.5*(P[1]+P[2]);
		}
	}
	return U;
}

template <typename T>
T Bethe_free_energy(Graph &G, BeliefPropagation & BP, const T& beta, double h) {
	T F = -beta*energy(G, BP, beta, h);
	for (int i : G.nodes()){
		F += (G.degree(i)-1) * ( BP.N_nodes[i].value * log0( BP.N_nodes[i].value ) );
		F += (G.degree(i)-1) * ((1- BP.N_nodes[i].value) * log0(1- BP.N_nodes[i].value ) );
		for (int j : G.neighbors(i)){
			auto P = BP.N_nodes[i].compute_two_point_distribution(beta, BP.N_diff, h, j); 
			F += 0.5*distribution_entropy(P);
		}
	}
	return F;
}


template <typename T, typename U>
void compute_counting_numbers(Graph &G, BeliefPropagation &BP, T& K, U& w){
	for (int i : BP.nodes){
		w[i];
		for (int j : BP.N_nodes[i].nodes){
			w[i][j] = 1;
		}
	}
	for (int i : BP.nodes){
		for (int j : BP.N_nodes[i].nodes){
			double cap = (double)BP.N_intersect[i][j].size;
			for (auto edge : BP.N_intersect[i][j].edges){
				w[edge.first][edge.second] -= 1./(cap * (cap-1));
				w[edge.second][edge.first] -= 1./(cap * (cap-1));
			}
		}
	}

	for (int i : BP.nodes) {
		K[i] = 1;
	}

	for (int i : BP.nodes) {
		for (int j : BP.N_nodes[i].nodes) {
			double cap = (double)BP.N_intersect[i][j].size;
			K[j] -= 1.0/ ((cap)*(cap-1.0));
			for (int k : BP.N_intersect[i][j].nodes){
				K[k] -= 1.0/ ((cap)*(cap-1.0));
			}
		}
		for (int j : G.neighbors(i)){
			K[i] -= w[i][j];
		}
	}

}

template <typename T>
T energy_neighborhoods(Graph &G, BeliefPropagation & BP, const T& beta, double h) {

	std::unordered_map<int, std::unordered_map<int, double>> w;
	std::unordered_map<int, T> K;
	compute_counting_numbers(G, BP, K, w);

	T U = 0.0;
	for (int i : BP.nodes){
		for (int j : BP.N_nodes[i].nodes){
			T cap = (T)BP.N_intersect[i][j].size;
			//U += BP.N_intersect[i][j].energy(beta, BP.N_diff, 0.0, BP.N_diff[j][i].value)/ ((cap)*(cap-1.0));
			U += BP.N_intersect[i][j].energy(beta, BP.N_diff, 0.0, BP.N_diff[j][i].value)/ ((cap)*(cap-1.0));

		}
	}

	for (int i : G.nodes()){
		for (int j : G.neighbors(i)){
			auto P = BP.N_intersect[i][j].compute_two_point_distribution(beta, BP.N_diff, 0, i, BP.N_diff[j][i].value); 
			//U -= 0.5*w[i][j]*(P[1]+P[2]);
			U += 0.5*w[i][j]*(P[1]+P[2]);
		}
	}

	for (int i : G.nodes()){
		T p_i = BP.N_nodes[i].value;
		//U += (K[i]-1.0)*p_i*h;
		U += p_i*h;
		//U += (K[i])*p_i*h;
	}
	return U;
}

template <typename T>
T entropy(Graph &G, BeliefPropagation & BP, const T& beta, double h) {

	std::unordered_map<int, std::unordered_map<int, double>> w;
	std::unordered_map<int, T> K;
	compute_counting_numbers(G, BP, K, w);
    
	T S = 0.0;

	// intersection entropy
	for (int i : BP.nodes){
		for (int j : BP.N_nodes[i].nodes){
			T cap = (T)BP.N_intersect[i][j].size;
			S += BP.N_intersect[i][j].entropy(beta, BP.N_diff, 0.0, BP.N_diff[j][i].value)/ ((cap)*(cap-1.0));
		}
	}
	// edge entropy 
	for (int i : G.nodes()){
		for (int j : G.neighbors(i)){
			//auto P = BP.N_nodes[i].compute_two_point_distribution(beta, BP.N_diff, h, j);
			auto P = BP.N_intersect[i][j].compute_two_point_distribution(beta, BP.N_diff, 0, i, BP.N_diff[j][i].value); 
			//S -= 0.5*w[i][j]*distribution_entropy(P);
			S += 0.5*w[i][j]*distribution_entropy(P);
		}
	}
    
	// node entropy
	for (int i : G.nodes()){
		T p_i = BP.N_nodes[i].value;
		//S -= (K[i]-1.0)* (p_i*log0(p_i) + (1-p_i)*log0(1-p_i));
		S += K[i] * (-(p_i*log0(p_i) + (1-p_i)*log0(1-p_i)));
		//S += (p_i*log0(p_i) + (1-p_i)*log0(1-p_i));
	}

	return S;
}

