#pragma once
#include "Graph.h"
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <tuple>

std::vector<std::vector<int>> list_paths(Graph &G, int i, int j, int r) {
    if (r==1) {
        if (G.has_edge(i,j))
            return {{i,j}};
        else
            return {};
    }
    else if (r>1) {
        std::vector<std::vector<int>> ans;
        for (int k : G.neighbors(i)){
            if (k!=i){
                for (auto x : list_paths(G, k, j, r-1)){
                    if (std::find(x.begin(), x.end(), i) == x.end() ) {
                        x.insert(x.begin(),i);
                        ans.push_back(x);
                    }
                }
            }
        }
        return ans;
    }
    return {};
}

unsigned long long int eid(Graph &G, int i, int j){
    if (i<j) {
        return G.number_of_nodes()*i + j;
    }
    else {
        return G.number_of_nodes()*j + i;
    }
}

std::pair<int,int> eid(Graph &G, unsigned long long int e){
    return {e/G.number_of_nodes(), e%G.number_of_nodes()};
}

std::vector<std::pair<int,int>> find_neighborhood_edges(Graph &G, int i, int r){
    std::unordered_set<unsigned long long int> edges;
    for (int j : G.neighbors(i)){
        edges.insert( eid(G,i,j) );
    }
    for (int j : G.neighbors(i)){
        for (int k : G.neighbors(i)){
            if (j<k) {
                for (int rr=1; rr<=r; ++rr){
                    for (auto path : list_paths(G, j, k, rr)){
                        for (int y=1; y<path.size(); ++y){
                            edges.insert( eid(G, path[y-1], path[y]) );
                        }
                    }
                }
            }
        }
    }
    std::vector<std::pair<int,int>> ans;
    for (auto edge : edges){
        ans.push_back( eid(G, edge) );
    }
    return ans;
}

std::vector<std::pair<int,int>> intersection( Graph &G, std::vector<std::pair<int,int>> N1, std::vector<std::pair<int,int>> N2 ){
    std::unordered_set<int> edges1;
    std::vector<std::pair<int,int>> edges_intersect;
    for (auto x : N1){
        edges1.insert( eid(G, x.first, x.second) );
    }
    for (auto x : N2){
        if (edges1.count( eid(G, x.first, x.second) ))
            edges_intersect.push_back(x);
    }
    return edges_intersect;
}

std::vector<std::pair<int,int>> difference( Graph &G, std::vector<std::pair<int,int>> N1, std::vector<std::pair<int,int>> N2 ){
    std::unordered_set<int> edges2;
    std::vector<std::pair<int,int>> edges_difference;
    for (auto x : N2){
        edges2.insert( eid(G, x.first, x.second) );
    }
    for (auto x : N1){
        if (!(edges2.count( eid(G, x.first, x.second) )))
            edges_difference.push_back(x);
    }
    return edges_difference;
}
