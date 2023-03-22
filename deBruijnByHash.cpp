//
//  deBruijnByHash.cpp
//  k-assembler
//
//  Created by Joe Song on 11/24/15.
//  Copyright © 2015 Joe Song. All rights reserved.
//
//  Updated 3/19/2018

#include "k-assembler.hpp"

#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>

using namespace std;

struct DNAHasher
// Hash function used for DNA sequence
{
    std::size_t operator()(const string & seq) const
    {
        // TO DO: Write a DNA sequence hash function here
        
        // BEGIN your code here:
        std::size_t hash_value = 0;
        std::size_t max_width = 200000;

        // Create a lookup table for the nucleotide codes
        std::unordered_map<char, int> nucleotide_codes = {
            {'A', 0},
            {'C', 1},
            {'G', 2},
            {'T', 3}
        };
        // Iterate over the sequence and compute the hash value
        for (std::size_t i = 0; i < seq.size() && i < max_width; ++i) {
            char nucleotide = toupper(seq[i]);
            if (nucleotide_codes.count(nucleotide)) {
                hash_value = (hash_value << 2) | nucleotide_codes[nucleotide];
            }
        }
        return hash_value;
        // END your code above
    }
};

struct AlphabetHasher
// An example hash function used for the English alphabet
{
    std::size_t operator()(const string & seq) const
    {
        size_t val = 0;
        size_t max_width=20;
        for(size_t i=0; i<seq.size() && i<max_width; ++i) {
            val = val << 5;
            val += tolower(seq[i])-'a';
        }
        return val;
    }
};

// define the hash table class
// typedef unordered_multimap<string, size_t, AlphabetHasher> CSeqHash;

typedef unordered_multimap<string, size_t, DNAHasher> CSeqHash;

CSeqHash create_hash_table(const vector<string> & kmers)
// create one hash table by inserting both the prefix and suffix of each
//   k-mer. The prefix and suffix is the key. Associated with each element
//   in the hash table is the node id for that prefix or suffix in the
//   de Bruijn graph to be constructed.
{
    CSeqHash ht;
    size_t node_id=0; // the node id will be used in the de Bruijn graph
    for (auto i=0u; i<kmers.size(); ++i) {
        for(auto j=0u; j<2; ++j) { // j=0: prefix; j=1: suffix
            auto key = kmers[i].substr(j, kmers[i].length()-1);
            if (ht.find(key) == ht.end()) {
                ht.insert(make_pair(key, node_id ++));
            }
        }
    }
    return ht;
}

void create_deBruijn_graph_by_hashing(const vector<string> & kmers, DiGraph & g)
// create a deBruijn graph by inserting all k-mers into the graph by hashing
{
    // TO DO:
    
   // BEGIN your code below:
   // Create hash table for both the k-1 prefix and suffix of each k-mer
    CSeqHash node_ids = create_hash_table(kmers);

    // Initialize empty node vector for graph g
    g.m_nodes.clear();

    for (auto i = 0u; i < kmers.size(); ++i) {
        // find the prefix node id from_id from the hash table
        auto prefix = kmers[i].substr(0, kmers[i].length() - 1);
        auto range = node_ids.equal_range(prefix);
        auto from_id = range.first->second;

        // update node from_id's label to prefix if necessary
        if (g.m_nodes.size() <= from_id) {
            g.m_nodes.resize(from_id + 1);
        }
        auto& node_from = g.m_nodes[from_id];
        if (node_from.m_label.empty()) {
            node_from.m_label = prefix;
        }

        // find the suffix node id to_id from the hash table
        auto suffix = kmers[i].substr(1, kmers[i].length() - 1);
        range = node_ids.equal_range(suffix);
        auto to_id = range.first->second;

        // update node to_id's label to suffix if necessary
        if (g.m_nodes.size() <= to_id) {
            g.m_nodes.resize(to_id + 1);
        }
        auto& node_to = g.m_nodes[to_id];
        if (node_to.m_label.empty()) {
            node_to.m_label = suffix;
        }
        // Create a new edge (from_id, to_id) by inserting node to_id into the adjacency list of node from_id
        g.m_nodes[from_id].m_outgoing.push_back(to_id);

        // Update the number of incoming edges of node to_id
        g.m_nodes[to_id].m_num_of_incoming++;
    }

}
