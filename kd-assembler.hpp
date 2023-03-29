//
//  kd-assembler.hpp
//  kd-assembler
//
//  Created by Joe Song on 3/19/18.
//  Copyright © 2018 Joe Song. All rights reserved.


#ifndef kd_assembler_h
#define kd_assembler_h

#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <sstream>

using namespace std;

struct Node_kd {
    pair <string,string> m_label; // a unique node label

    // adjacency list: the same node index can
    //   show up mutiple times on the list to
    //   indicate multiple edges between the same
    //   pair of nodes
    list<size_t> m_outgoing;

    // the total number of incoming edges
    size_t m_num_of_incoming;
};

struct DiGraph_kd { // directed graph data structure
    vector<Node_kd> m_nodes;
};

vector<pair <string,string> > get_kdmers(string seq, size_t k, size_t d, bool randomized=true);
void create_deBruijn_graph_by_hashing_kd(const vector<pair<string,string> >& kdmers, DiGraph_kd& g);
string build_kd_sequence(const list<size_t>& path, const DiGraph_kd& g, const size_t k, const size_t d);
void printDOTFileKd(const DiGraph_kd& g, const string& file);
list<size_t> find_Eulerian_path(DiGraph_kd& g);
bool has_Eulerian_path(const DiGraph_kd& g);
string assemble_kdmers(const vector<pair<string, string> >& kdmers, const string& method, const string & dotfile, const size_t k, const size_t d);

#endif /* k_assembler_h */
 