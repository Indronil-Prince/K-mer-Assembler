//
//  k-assembler.cpp
//  k-assembler
//
//  Created by Joe Song on 3/19/18.
//  Copyright � 2018 Joe Song. All rights reserved.
/*

#include "k-assembler.hpp"
#include "kd-assembler.hpp"
#include <unordered_map>
#include <vector>
#include <fstream>

void create_deBruijn_graph_by_kdmer(const vector<string>& reads, DiGraph& g, size_t k, size_t d) {
    unordered_map<string, size_t> node_index; // map node labels to node indices

    
        }
    }
}

vector<string> generate_kdmers(const string& read, size_t k, size_t d) {
    vector<string> kmers;
    size_t n = read.length();
    for (size_t i = 0; i+k+d <= n-k; ++i) {
        string kmer = read.substr(i, k);
        kmer+= read.substr(i+k+d, k);
        kmers.push_back(kmer);
    }
    return kmers;
}

string assemble_kdmers(const vector<string>& kdmers, size_t k, size_t d, const string& dotfile) {
    string seq;
    DiGraph g;

    // Create de Bruijn graph from k-dmers
    create_deBruijn_graph_by_kdmer(kdmers, g, k, d);

    if (!dotfile.empty()) {
        printDOTFile(g, dotfile);
    }

    // Check if the graph has an Eulerian path
    if (!has_Eulerian_path(g)) {
        throw "ERROR: Eulerian path does not exist!";
    }
    else {
        // Find Eulerian path and build sequence
        list<size_t> path = find_Eulerian_path(g);
        seq = build_kd_sequence(path, g, k, d);
    }

    return seq;
}

string build_kd_sequence(const list<size_t>& path, const DiGraph& g, const size_t k, const size_t d)
{
    vector<Node> nodes = g.m_nodes;

    size_t overlap_length = k - 1 - d;
    size_t seq_length = (path.size() - 1) * (k - d) + k;
    string seq(seq_length, 'N');

    auto pos = path.begin();
    size_t i = 0;

    // First k-mer
    string first_kmer = nodes[*pos].m_label;
    seq.replace(0, k, first_kmer);
    i += k;

    // Subsequent k-mers with overlaps
    pos++;
    while (pos != path.end()) {
        string kmer = nodes[*pos].m_label;
        string overlap = kmer.substr(0, overlap_length);
        seq.replace(i - d, d + overlap_length, overlap);
        seq.replace(i, k - overlap_length, kmer.substr(overlap_length));
        i += k - d;
        pos++;
    }

    return seq;
}
*/