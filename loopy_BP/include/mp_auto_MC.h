#ifndef MP_HEADER 
#define MP_HEADER

#include "Graph.h"
#include "find_neighborhoods.h"
#include "neighborhoods.h"
#include <math.h>
#include <random>
#include <boost/math/differentiation/autodiff.hpp>

class BeliefPropagation {
    public:
        std::unordered_map<int, Neighborhood> N_nodes;
        std::unordered_map<int, std::unordered_map<int, Neighborhood>> N_diff;
        std::unordered_map<int, std::unordered_map<int, Neighborhood>> N_intersect;
        a_Double sum_marginal;
        a_Double M_total;
        a_Double U_total;
        a_Double C_total;
        a_Double total_num_edges;
        void diff_iteration(a_Double const&, double const&, int);
        void iteration(a_Double const&, double const&, int);
        void init(Graph &, int);
        std::vector<int> nodes;
};


void BeliefPropagation::init(Graph &G, int r){
    std::unordered_set<int> node_set = G.nodes();
    nodes.insert(nodes.end(), node_set.begin(), node_set.end());
    total_num_edges = G.number_of_edges();

    for (int i : nodes){
        N_nodes[i];
        N_nodes[i].init( find_neighborhood_edges(G,i,r) ,i);
    }

    for (int i : nodes){
        N_diff[i];
        for ( int j : N_nodes[i].nodes ){
            N_diff[i][j];
            N_diff[i][j].init(difference(G, N_nodes[i].edges, N_nodes[j].edges), i);
        }
    }

    for (int j : nodes){
        N_intersect[j];
        for (int i : N_nodes[j].nodes){
            N_intersect[j][i];
            N_intersect[j][i].init(intersection(G, N_nodes[i].edges, N_nodes[j].edges), i);
        }
    }
}

void BeliefPropagation::diff_iteration( a_Double const &beta, double const &h, int MC=0) {

    std::shuffle(nodes.begin(), nodes.end(), RNG);
    #pragma omp parallel for
    for (int i=0; i<nodes.size(); ++i){
        for (int j : N_nodes.at(i).nodes){
            if (MC)
                N_diff[i][j].compute_value_MC(beta, N_diff, h);
            else
                N_diff[i][j].compute_value(beta, N_diff, h);
        }
    }
    
    std::shuffle(nodes.begin(), nodes.end(), RNG);
    #pragma omp parallel for
    for (int i=0; i<nodes.size(); ++i){
        for (int j : N_nodes.at(i).nodes){
            if (MC)
                N_diff[i][j].compute_value_der_MC(beta, N_diff, h);
            else
                N_diff[i][j].compute_value_der(beta, N_diff, h);
        }
    }

}

void BeliefPropagation::iteration( a_Double const &beta, double const &h, int MC=0) {

    std::shuffle(nodes.begin(), nodes.end(), RNG);
    #pragma omp parallel for
    for (int i=0; i<nodes.size(); ++i){
        for (int j : N_nodes.at(i).nodes){
            if (MC)
                N_diff[i][j].compute_value_MC(beta, N_diff, h);
            else
                N_diff[i][j].compute_value(beta, N_diff, h);
        }
    }
    
    std::shuffle(nodes.begin(), nodes.end(), RNG);
    #pragma omp parallel for
    for (int i=0; i<nodes.size(); ++i){
        for (int j : N_nodes.at(i).nodes){
            if (MC)
                N_diff[i][j].compute_value_der_MC(beta, N_diff, h);
            else
                N_diff[i][j].compute_value_der(beta, N_diff, h);
        }
    }

    std::shuffle(nodes.begin(), nodes.end(), RNG);
    #pragma omp parallel for
    for (int i=0; i<nodes.size(); ++i){
        if (MC)
            M_total += N_nodes.at(i).compute_value_MC(beta, N_diff, h);
        else
            M_total += N_nodes.at(i).compute_value(beta, N_diff, h);
        if (MC)
            U_total += N_nodes.at(i).compute_local_U_MC(beta, N_diff, h);
        else
            U_total += N_nodes.at(i).compute_local_U(beta, N_diff, h);
        if (MC)
            C_total += N_nodes.at(i).compute_local_C_MC(beta, N_diff, h);
        else
            C_total += N_nodes.at(i).compute_local_C(beta, N_diff, h);
         
    }

    M_total = 0.0;
    U_total = 0.0;
    C_total = 0.0;
    for (int i : nodes){
        M_total += N_nodes.at(i).value;
        U_total += N_nodes.at(i).local_U/2.0;
        C_total += -beta*beta*N_nodes.at(i).local_C;
    }
    M_total = 2.0*(double)M_total/nodes.size()-1.0;
    U_total = 2.0*(double)U_total - total_num_edges;
    C_total *= 4.0;

}

#endif
