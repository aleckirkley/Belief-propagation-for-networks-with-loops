#ifndef NEIGHBORHOOD_HEADER 
#define NEIGHBORHOOD_HEADER

#include "Graph.h"
#include "find_neighborhoods.h"
#include "Wolff_savestates.h"
#include <math.h>
#include <random>
#include <boost/math/differentiation/autodiff.hpp>

std::mt19937 RNG(RD());

const int max_neig_size_exact = 20;
const int samples_per_spin = 1000;
const int max_samples = 100000;
const int wolff_interval = 1;

typedef boost::math::differentiation::detail::fvar<double,4> a_Double; 

typedef std::pair<int,int> EDGE;
typedef std::vector<EDGE> EDGE_LIST;

class Neighborhood;

typedef std::unordered_map<int,std::unordered_map<int,Neighborhood>> N_map;


template <typename T>
inline T log0(T const x) {
    return (x<=0) ? 0 : log(x);
}

bool iterate_spin_configuration( std::vector<int> &S ) {
    for (int i=0; i<S.size(); ++i)
        if (1==(S[i]=(1-S[i]))) return true;
    return false;
}

class Neighborhood {
    public:
        int i;
        EDGE_LIST edges;
        std::unordered_set<int> nodes;
        a_Double value;
        a_Double value_der;
        a_Double local_U;
        a_Double local_C;
        int size;
        std::unordered_map<int,int> local_name;
        void init(EDGE_LIST, int);
        std::pair<a_Double,a_Double> compute_Z0_Z1(a_Double const&, N_map&, double const&);
        a_Double compute_value(a_Double const&, N_map&, double const&);
        a_Double compute_value_MC(a_Double const&, N_map&, double const&);
        a_Double compute_value_der(a_Double const&, N_map&, double const&);
        a_Double compute_value_der_MC(a_Double const&, N_map&, double const&);
        a_Double compute_local_U(a_Double const&, N_map&, double const&);
        a_Double compute_local_U_MC(a_Double const&, N_map&, double const&);
        a_Double compute_local_C(a_Double const&, N_map&, double const&);
        a_Double compute_local_C_MC(a_Double const&, N_map&, double const&);
        std::vector<a_Double> compute_Z00_Z11(a_Double const&, N_map&, double const&, int);
        a_Double compute_two_point(a_Double const&, N_map&, double const&, int);
	std::vector<a_Double> compute_two_point_distribution(a_Double const&, N_map&,
			double const&, int, a_Double const&);
	std::vector<a_Double> compute_two_point_distribution(a_Double const&, N_map&,
			double const&, int);
        WolffSampler WS;
        a_Double energy(a_Double const&, N_map&, double const&);
        a_Double energy(a_Double const&, N_map&, double const&, a_Double const& );
        a_Double entropy(a_Double const&, N_map&, double const&, a_Double const& );
};

void Neighborhood::init(EDGE_LIST neighborhood_edges, int node) {
    i = node;
    edges = neighborhood_edges;
    nodes.insert(i);
    for (auto edge : edges){
        nodes.insert(edge.first);
        nodes.insert(edge.second);
    }
    size = nodes.size();
    nodes.erase(i);
    value = 0.5 + 0.0000001*(unif(RNG));
    value_der = 0.1*(unif(RNG));
    local_name[i] = 0;
    int temp = 1;
    for ( auto j : nodes ){
        local_name[j] = temp;
        ++temp;
    }

    Graph G_local;
    G_local.add_node(0);
    for (int node : nodes)
        G_local.add_node(local_name[node]);
    for (auto edge : edges){
        G_local.add_edge(local_name[edge.first],local_name[edge.second]);
    }
    WolffSampler WS_temp(G_local);
    WS = WS_temp;
}

std::pair<a_Double, a_Double> Neighborhood::compute_Z0_Z1( a_Double const& beta, N_map &N_diff, double const& h ) {
    a_Double Z0=0;
    a_Double Z1=0;
    std::unordered_map<int, a_Double> q;
    for (auto k : nodes) 
        q[k] = N_diff.at(k).at(i).value;
    std::vector<int> S(size,0);

    do {
        a_Double H = (S[0]==1) ? h : 0;
        for (auto edge : edges)
            H += (S[local_name[edge.first]]!=S[local_name[edge.second]]);
        a_Double P = exp(-beta*H);
        for ( int j : nodes ) 
            P *= (S[local_name[j]]==1) ? q[j] : 1-q[j];
        if (S[0]==0)
            Z0 += P;
        else
            Z1 += P;
    } while( iterate_spin_configuration(S) );

    return {Z0, Z1};
}

a_Double Neighborhood::energy( a_Double const& beta, N_map &N_diff, double const& h ) {
    a_Double Z=0;
    a_Double U=0;
    std::unordered_map<int, a_Double> q;
    for (auto k : nodes) 
        q[k] = N_diff.at(k).at(i).value;
    std::vector<int> S(size,0);

    do {
        a_Double H = (S[0]==1) ? h : 0;
        for (auto edge : edges)
            H += (S[local_name[edge.first]]!=S[local_name[edge.second]]);
        a_Double P = exp(-beta*H);
        for ( int j : nodes ) 
            P *= (S[local_name[j]]==1) ? q[j] : 1-q[j];
        Z += P;
        U += H*P;
    } while( iterate_spin_configuration(S) );

    return U/Z;
}

a_Double Neighborhood::energy( a_Double const& beta, N_map &N_diff, double const& h, a_Double const &qq ) {
    a_Double Z=0;
    a_Double U=0;
    std::unordered_map<int, a_Double> q;
    for (auto k : nodes) 
        q[k] = N_diff.at(k).at(i).value;
    std::vector<int> S(size,0);

    do {
        a_Double H = (S[0]==1) ? h : 0;
        for (auto edge : edges)
            H += (S[local_name[edge.first]]!=S[local_name[edge.second]]);
        a_Double P = exp(-beta*H);
        P *= (S[0]==1) ? qq : 1-qq;
        for ( int j : nodes ) 
            P *= (S[local_name[j]]==1) ? q[j] : 1-q[j];
        Z += P;
        U += H*P;
    } while( iterate_spin_configuration(S) );

    return U/Z;
}

a_Double Neighborhood::entropy( a_Double const& beta, N_map &N_diff,
		double const& h, a_Double const &qq ) {
    a_Double Z=0;
    a_Double QlnQ=0;
    std::unordered_map<int, a_Double> q;
    for (auto k : nodes) 
        q[k] = N_diff.at(k).at(i).value;
    std::vector<int> S(size,0);

    do {
        a_Double H = (S[0]==1) ? h : 0;
        for (auto edge : edges)
            H += (S[local_name[edge.first]]!=S[local_name[edge.second]]);
        a_Double P = exp(-beta*H);
        P *= (S[0]==1) ? qq : 1-qq;
        for ( int j : nodes ) 
            P *= (S[local_name[j]]==1) ? q[j] : 1-q[j];
        Z += P;
        QlnQ += P*log0(P);
    } while( iterate_spin_configuration(S) );

    return (-QlnQ/Z + log(Z));
}

std::vector<a_Double> Neighborhood::compute_two_point_distribution( a_Double const& beta,
		N_map &N_diff, double const& h, int j2, a_Double const& qq ) {
    std::vector<a_Double> Z = {0,0,0,0};
    std::unordered_map<int, a_Double> q;
    for (auto k : nodes) 
        q[k] = N_diff.at(k).at(i).value;
    std::vector<int> S(size,0);

    do {
        a_Double H = (S[0]==1) ? h : 0;
        for (auto edge : edges)
            H += (S[local_name[edge.first]]!=S[local_name[edge.second]]);
        a_Double P = exp(-beta*H);
        P *= (S[0]==1) ? qq : 1-qq;
        for ( int j : nodes ) 
            P *= (S[local_name[j]]==1) ? q[j] : 1-q[j];
        Z[(2*(S[0]==1))+(S[local_name[j2]]==1)] += P;
    } while( iterate_spin_configuration(S) );

    auto sum_Z = Z[0]+Z[1]+Z[2]+Z[3];

    return {Z[0]/sum_Z, Z[1]/sum_Z, Z[2]/sum_Z, Z[3]/sum_Z};
}

std::vector<a_Double> Neighborhood::compute_two_point_distribution( a_Double const& beta,
		N_map &N_diff, double const& h, int j2) {
    std::vector<a_Double> Z = {0,0,0,0};
    std::unordered_map<int, a_Double> q;
    for (auto k : nodes) 
        q[k] = N_diff.at(k).at(i).value;
    std::vector<int> S(size,0);

    do {
        a_Double H = (S[0]==1) ? h : 0;
        for (auto edge : edges)
            H += (S[local_name[edge.first]]!=S[local_name[edge.second]]);
        a_Double P = exp(-beta*H);
        for ( int j : nodes ) 
            P *= (S[local_name[j]]==1) ? q[j] : 1-q[j];
        Z[(2*(S[0]==1))+(S[local_name[j2]]==1)] += P;
    } while( iterate_spin_configuration(S) );

    auto sum_Z = Z[0]+Z[1]+Z[2]+Z[3];

    return {Z[0]/sum_Z, Z[1]/sum_Z, Z[2]/sum_Z, Z[3]/sum_Z};
}


a_Double Neighborhood::compute_value( a_Double const& beta, N_map &N_diff, double const& h){
    
    a_Double Z=0;
    a_Double numerator=0;
    std::unordered_map<int, a_Double> q;
    for (auto k : nodes) 
        q[k] = N_diff.at(k).at(i).value;
    std::vector<int> S(size,0);

    do {
        a_Double H = 0;
        for (auto edge : edges)
            H += (S[local_name[edge.first]]!=S[local_name[edge.second]]);
        
        a_Double P = exp(-2.0*beta*H);
        for ( int j : nodes ) 
            P *= (S[local_name[j]]==1) ? q[j] : 1-q[j];
        Z += P;
        if (S[0] == 1){
            numerator += P;
        }
        
    } while( iterate_spin_configuration(S) );

    value = numerator/Z;
    return value;
}

a_Double Neighborhood::compute_value_MC( a_Double const& beta, N_map &N_diff, double const& h ) {
    if (size<=max_neig_size_exact){
        return compute_value(beta, N_diff, h);
    }

    a_Double numerator=0;
    a_Double Z=0;

    int num_samples = (size+1)*samples_per_spin;
    if (num_samples > max_samples){
        num_samples = max_samples;
    }
 
    WS.sample_states(num_samples, (double)(2.0*beta), wolff_interval); 
    double beta_W = 0;//-log(1-WS.p);

    for (int x=0; x<num_samples; ++x) {
        a_Double H = WS.H[x];
        a_Double p_up = 1;//exp(-2.0*beta*H);
        a_Double p_down = 1;//exp(-2.0*beta*H);
        for ( int j : nodes ) {
            a_Double qj = N_diff.at(j).at(i).value;
            if (WS.states[x][local_name[j]]==1) {
                p_up *= qj;
                p_down *= 1-qj;
            }
            else {
                p_up *= 1-qj;
                p_down *= qj;
            }
        
        }
        numerator += p_up;
        Z += p_up + p_down;
    }
    value = numerator/Z ;
    if (value < 0.00001) value = 0.00001;
    if (value > 0.99999) value = 0.99999;
    return value ;
}


a_Double Neighborhood::compute_value_der( a_Double const& beta, N_map &N_diff, double const& h){
    
    a_Double Z=0;
    a_Double numerator = 0;
    std::unordered_map<int, a_Double> q;
    std::unordered_map<int, a_Double> eta;
    for (auto k : nodes) 
        q[k] = N_diff.at(k).at(i).value;
    for (auto k : nodes) 
        eta[k] = N_diff.at(k).at(i).value_der;
    std::vector<int> S(size,0);
    
    do {
        a_Double H = 0;
        for (auto edge : edges)
            H += (S[local_name[edge.first]]!=S[local_name[edge.second]]);
        
        a_Double eta_sum = 0;
        for ( int k : nodes ) {
            eta_sum += (S[local_name[k]]==1) ? eta[k]/(q[k] + 0.0000000001) : -eta[k]/(1.0-q[k] + 0.0000000001); 
        }
        
        a_Double prefactor = (S[0]==1) ? value-1.0 : value;
        
        a_Double P = exp(-2.0*beta*H);
        for ( int j : nodes ) 
            P *= (S[local_name[j]]==1) ? q[j] : 1-q[j];
        Z += P;
        numerator += prefactor*(H - eta_sum)*P;
    } while( iterate_spin_configuration(S) );
    
    value_der = numerator/Z;
    return value_der;
}

a_Double Neighborhood::compute_value_der_MC( a_Double const& beta, N_map &N_diff, double const& h ) {
    if (size<=max_neig_size_exact){
        return compute_value_der(beta, N_diff, h);
    }
    
    a_Double Z_up=0;
    a_Double Z_down=0;
    a_Double numerator_up=0;
    a_Double numerator_down=0;

    int num_samples = (size+1)*samples_per_spin;
    if (num_samples > max_samples){
        num_samples = max_samples;
    }
    WS.sample_states(num_samples, (double)(2.0*beta), wolff_interval); 
    double beta_W = 0;//-log(1-WS.p);
    double num_edges = (double)edges.size();
    a_Double q_val = 0;
    a_Double eta_val = 0;
    
    for (int x=0; x<num_samples; ++x) {
        a_Double H_up = WS.H[x];
        a_Double H_down = WS.H[x];
        a_Double p_up = 1;//exp(-2.0*beta*H_up);
        a_Double p_down = 1;//exp(-2.0*beta*H_down);
        a_Double eta_sum_up = 0;
        a_Double eta_sum_down = 0;
        
        for ( int k : nodes ) {
            q_val = N_diff.at(k).at(i).value;
            eta_val = N_diff.at(k).at(i).value_der;
            eta_sum_up += (WS.states[x][local_name[k]]==1) ? eta_val/(q_val + 0.0000000001) : -eta_val/(1.0-q_val + 0.0000000001); 
            eta_sum_down += (WS.states[x][local_name[k]]==1) ? -eta_val/(1.0-q_val + 0.0000000001) : eta_val/(q_val + 0.0000000001); 
        }
    
        for ( int j : nodes ) {
            a_Double qj = N_diff.at(j).at(i).value;
            if (WS.states[x][local_name[j]]==1) {
                p_up *= qj;
                p_down *= 1-qj;
            }
            else {
                p_up *= 1-qj;
                p_down *= qj;
            }

        }
        
        a_Double prefactor_up = value - 1.0;
        a_Double prefactor_down = value;
        numerator_up += prefactor_up*(H_up - eta_sum_up)*p_up;
        numerator_down += prefactor_down*(H_down - eta_sum_down)*p_down;
        Z_up += p_up;
        Z_down += p_down;
    }
    value_der = (numerator_up+numerator_down)/(Z_up+Z_down);

    return value_der ;
}

a_Double Neighborhood::compute_local_U( a_Double const& beta, N_map &N_diff, double const& h){
    a_Double Z=0;
    a_Double Upartial=0;
    std::unordered_map<int, a_Double> q;
    for (auto k : nodes) 
        q[k] = N_diff.at(k).at(i).value;
    std::vector<int> S(size,0);

    do {
        a_Double Hpartial = 0;
        a_Double Hneig = 0;
        a_Double delta_H = 0;
        for (EDGE edge : edges){
            delta_H = (S[local_name[edge.first]]!=S[local_name[edge.second]]);
            if (local_name[edge.first]*local_name[edge.second] == 0){
                Hpartial += delta_H;
            }
            Hneig += delta_H;    
        }    
        a_Double P = exp(-2.0*beta*Hneig);
        for ( int j : nodes ) 
            P *= (S[local_name[j]]==1) ? q[j] : 1-q[j];
        Z += P;
        Upartial += Hpartial*P;
    } while( iterate_spin_configuration(S) );

    local_U = Upartial/Z;
    return local_U;
}


a_Double Neighborhood::compute_local_U_MC( a_Double const& beta, N_map &N_diff, double const& h ) {
    if (size<=max_neig_size_exact){
        return compute_local_U(beta, N_diff, h);
    }
    
    a_Double numerator_up=0;
    a_Double numerator_down=0;
    a_Double Z_up=0;
    a_Double Z_down=0;

    int num_samples = (size+1)*samples_per_spin;
    if (num_samples > max_samples){
        num_samples = max_samples;
    }
    WS.sample_states(num_samples, (double)(2.0*beta), wolff_interval); 
    double beta_W = 0; //-log(1-WS.p);
    
    for (int x=0; x<num_samples; ++x) {
        a_Double Hneig_up = WS.H[x];
        a_Double Hneig_down = WS.H[x];
        a_Double p_up = 1;//exp(-2.0*beta*Hneig_up);
        a_Double p_down = 1;//exp(-2.0*beta*Hneig_down);
        
        a_Double Hpartial_up = 0.0;
        a_Double Hpartial_down = 0.0;
        a_Double delta_H = 0;
        for (EDGE edge : edges){
            if (local_name[edge.first]*local_name[edge.second] == 0){
                delta_H = (WS.states[x][local_name[edge.first]]!=WS.states[x][local_name[edge.second]]);
                Hpartial_up += delta_H;
                Hpartial_down += delta_H;
            }
        }
    
        for ( int j : nodes ) {
            a_Double qj = N_diff.at(j).at(i).value;
            if (WS.states[x][local_name[j]]==1) {
                p_up *= qj;
                p_down *= 1-qj;
            }
            else {
                p_up *= 1-qj;
                p_down *= qj;
            }
            //p_up *= (WS.states[x][local_name[j]]==1) ? qj : 1-qj;
            //p_down *= (WS.states[x][local_name[j]]==1) ? 1-qj : qj;
        }
        
        numerator_up += Hpartial_up*p_up;
        numerator_down += Hpartial_down*p_down;
        Z_up += p_up;
        Z_down += p_down;
    }
    local_U = (numerator_up+numerator_down)/(Z_up+Z_down);
    return local_U;
}

a_Double Neighborhood::compute_local_C( a_Double const& beta, N_map &N_diff, double const& h){
    a_Double Z=0;
    a_Double Uneig=0;
    a_Double Upartial=0;
    a_Double avg_eta_sum=0;
    a_Double UU=0;
    a_Double Ueta_sum=0;
    std::unordered_map<int, a_Double> q;
    std::unordered_map<int, a_Double> eta;
    for (auto k : nodes) 
        q[k] = N_diff.at(k).at(i).value;
    for (auto k : nodes) 
        eta[k] = N_diff.at(k).at(i).value_der;
    std::vector<int> S(size,0);

    do {
        a_Double Hpartial = 0.0;
        a_Double Hneig = 0.0;
        a_Double delta_H = 0.0;
        for (EDGE edge : edges){
            delta_H = (S[local_name[edge.first]]!=S[local_name[edge.second]]);
            if (local_name[edge.first]*local_name[edge.second] == 0){
                Hpartial += delta_H;
            }
            Hneig += delta_H;
        }
        a_Double eta_sum = 0;
        for ( int k : nodes ) {
            eta_sum += (S[local_name[k]]==1) ? eta[k]/(q[k] + 0.0000000001) : -eta[k]/(1.0-q[k] + 0.0000000001); 
        }
    
        a_Double P = exp(-2.0*beta*Hneig);
        for ( int j : nodes ) 
            P *= (S[local_name[j]]==1) ? q[j] : 1-q[j];
        Z += P;
        Uneig += Hneig*P;
        Upartial += Hpartial*P;
        avg_eta_sum += eta_sum*P;
        UU += Hpartial*Hneig*P;
        Ueta_sum += Hpartial*eta_sum*P;
    } while( iterate_spin_configuration(S) );

    local_C = Uneig*Upartial/2.0/Z/Z - UU/2.0/Z + Ueta_sum/2.0/Z - Upartial*avg_eta_sum/2.0/Z/Z;
    return local_C;
}

a_Double Neighborhood::compute_local_C_MC( a_Double const& beta, N_map &N_diff, double const& h ) {
    if (size<=max_neig_size_exact){
        return compute_local_C(beta, N_diff, h);
    }
    
    a_Double Uneig_up=0;
    a_Double Upartial_up=0;
    a_Double avg_eta_sum_up=0;
    a_Double UU_up=0;
    a_Double Ueta_sum_up=0;
    a_Double Uneig_down=0;
    a_Double Upartial_down=0;
    a_Double avg_eta_sum_down=0;
    a_Double UU_down=0;
    a_Double Ueta_sum_down=0;
    a_Double Z_up=0;
    a_Double Z_down=0;
    a_Double q_val=0;
    a_Double eta_val=0;

    int num_samples = (size+1)*samples_per_spin;
    if (num_samples > max_samples){
        num_samples = max_samples;
    }
    WS.sample_states(num_samples, (double)(2.0*beta), wolff_interval); 
    double beta_W = 0; //-log(1-WS.p);
    
    for (int x=0; x<num_samples; ++x) {
        a_Double Hneig_up = WS.H[x];
        a_Double Hneig_down = WS.H[x];
        a_Double p_up = 1;//exp(-2.0*beta*Hneig_up);
        a_Double p_down = 1;//exp(-2.0*beta*Hneig_down);
        
        a_Double Hpartial_up = 0.0;
        a_Double Hpartial_down = 0.0;
        a_Double delta_H = 0;
        for (EDGE edge : edges){
            if (local_name[edge.first]*local_name[edge.second] == 0){
                delta_H = (WS.states[x][local_name[edge.first]]!=WS.states[x][local_name[edge.second]]);
                Hpartial_up += delta_H;
                Hpartial_down += delta_H;
            }
        }
        a_Double eta_sum_up = 0;
        a_Double eta_sum_down = 0;
        for ( int k : nodes ) {
            q_val = N_diff.at(k).at(i).value;
            eta_val = N_diff.at(k).at(i).value_der;
            eta_sum_up += (WS.states[x][local_name[k]]==1) ? eta_val/(q_val + 0.0000000001) : -eta_val/(1.0-q_val + 0.0000000001); 
            eta_sum_down += (WS.states[x][local_name[k]]==1) ? -eta_val/(1.0-q_val + 0.0000000001) : eta_val/(q_val + 0.0000000001); 
        }
        
        for ( int j : nodes ) {
            a_Double qj = N_diff.at(j).at(i).value;
            if (WS.states[x][local_name[j]]==1) {
                p_up *= qj;
                p_down *= 1-qj;
            }
            else {
                p_up *= 1-qj;
                p_down *= qj;
            }

        }
        
        Uneig_up += Hneig_up*p_up;
        Uneig_down += Hneig_down*p_down;
        Upartial_up += Hpartial_up*p_up;
        Upartial_down += Hpartial_down*p_down;
        avg_eta_sum_up += eta_sum_up*p_up;
        avg_eta_sum_down += eta_sum_down*p_down;
        UU_up += Hneig_up*Hpartial_up*p_up;
        UU_down += Hneig_down*Hpartial_down*p_down;
        Ueta_sum_up += eta_sum_up*Hpartial_up*p_up;
        Ueta_sum_down += eta_sum_down*Hpartial_down*p_down;
        Z_up += p_up;
        Z_down += p_down;
    }
    a_Double Uneig = (Uneig_up+Uneig_down)/(Z_up+Z_down);
    a_Double Upartial = (Upartial_up+Upartial_down)/(Z_up+Z_down);
    a_Double avg_eta_sum = (avg_eta_sum_up+avg_eta_sum_down)/(Z_up+Z_down);
    a_Double UU = (UU_up+UU_down)/(Z_up+Z_down);
    a_Double Ueta_sum = (Ueta_sum_up+Ueta_sum_down)/(Z_up+Z_down);

    local_C = Uneig*Upartial/2.0 - UU/2.0 + Ueta_sum/2.0 - Upartial*avg_eta_sum/2.0;
    return local_C;
}








#endif
