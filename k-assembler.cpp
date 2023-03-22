//
//  k-assembler.cpp
//  k-assembler
//
//  Created by Joe Song on 3/19/18.
//  Copyright © 2018 Joe Song. All rights reserved.
//

#include "k-assembler.hpp"
#include <fstream>


string assemble_kmers(const vector<string> & kmers, const string & method,const string & dotfile)
// assembler the given k-mers using the specified method and
//   return the assembled sequence
{
    string seq;
    
    DiGraph g;
    
    if(method == "k-mer pairwise comparison") {
        create_deBruijn_graph_by_string_comp(kmers, g);
    } else if(method == "k-mer hashing") {
        create_deBruijn_graph_by_hashing(kmers, g);
    } else {
        throw "ERROR: unknown methods!";
    }
    
    if (!dotfile.empty()) {
        printDOTFile(g, dotfile);
    }
    
    if(! has_Eulerian_path(g)) {
        
        throw "ERROR: Eulerian path does not exist!";
        
    } else {
        
        list<size_t> path = find_Eulerian_path(g);
        
        seq = build_sequence(path, g);
    }
    
    return seq;
}

string build_sequence(const list<size_t> & path, const DiGraph & g)
// build a sequence by following a path in the input graph
{
    vector<Node> nodes = g.m_nodes;
    auto k = nodes[path.front()].m_label.size() + 1u;
    
    string seq;
    
    seq.resize(k - 1u + path.size() - 1u);
    
    seq.replace(0, k-1u, nodes[path.front()].m_label);
    
    auto pos = path.begin();
    pos ++;
    
    size_t i = k-1u;
    
    while (pos != path.end()) {
        
        seq[i] = nodes[*pos].m_label.back();
        i ++;
        pos ++;
    }
    
    return seq;
}

void printDOTFile(const DiGraph & g, const string & file)
{
    ofstream ofs(file);
    
    ofs << "digraph {" << endl;
    ofs << "label=\"de Bruijn graph\"" << endl;
    
    auto & nodes = g.m_nodes;
    
    for (auto & node: nodes) {
        for (auto to = node.m_outgoing.begin(); to != node.m_outgoing.end(); ++ to) {
            auto & prefix = node.m_label;
            auto & suffix = nodes[*to].m_label;
            ofs << prefix << "->" << suffix
            << "[label=" << prefix << suffix.back() << "];" << endl;
        }
    }
    
    ofs << "}" << endl;
}

