//
//  k-assembler.hpp
//  k-assembler
//
//  Created by Joe Song on 3/19/18.
//  Copyright © 2018 Joe Song. All rights reserved.
//

#ifndef k_assembler_h
#define k_assembler_h

#include <list>
#include <vector>
#include <string>

using namespace std;

struct Node {
    string m_label; // a unique node label
    
    // adjacency list: the same node index can
    //   show up mutiple times on the list to
    //   indicate multiple edges between the same
    //   pair of nodes
    list<size_t> m_outgoing;
    
    // the total number of incoming edges
    size_t m_num_of_incoming;
};

struct DiGraph { // directed graph data structure
    vector<Node> m_nodes;
};

void create_deBruijn_graph_by_string_comp(const vector<string> & kmers, DiGraph & g);

void create_deBruijn_graph_by_hashing(const vector<string> & kmers, DiGraph & g);

list<size_t> find_Eulerian_path(DiGraph & g);

bool has_Eulerian_path(const DiGraph & g);

string build_sequence(const list<size_t> & path, const DiGraph & g);

string assemble_kmers(const vector<string> & kmers, const string & method,
                      const string & dotfile="");

void printDOTFile(const DiGraph & g, const string & file);

#endif /* k_assembler_h */
