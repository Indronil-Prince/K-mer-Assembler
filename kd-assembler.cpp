//
//  k-assembler.cpp
//  k-assembler
//
//  Created by Joe Song on 3/19/18.
//  Copyright � 2018 Joe Song. All rights reserved.


#include "kd-assembler.hpp"
#include <unordered_map>
#include <stack>
#include <vector>
#include <fstream>

vector< pair<string, string> > get_kdmers(string seq, size_t k, size_t d, bool randomized)
// obtain all kd-mers of a given sequence. The order of the kd-mers is randomized by default.
{
    vector< pair <string, string> > kdmers(seq.length() - (2 * k) - d + 1);

    for (size_t i = 0; i < kdmers.size(); ++i) {
        kdmers[i].first = seq.substr(i, k);
        kdmers[i].second = seq.substr(i + k + d, k);
    }

    if (randomized) {

        std::random_device rd;
        std::mt19937 gen(rd());

        size_t nkdmers = kdmers.size();

        for (size_t i = 0; i < nkdmers - 1; ++i) {
            std::uniform_int_distribution<size_t> dis(i, kdmers.size() - 1);
            size_t j = dis(gen);
            pair<string, string> kdmer = make_pair(kdmers[j].first, kdmers[j].second);
            kdmers[j] = kdmers[i];
            kdmers[i] = kdmer;
        }
    }

    return kdmers;
}
struct DNAHasherKd
    // Hash function used for DNA sequence
{
    std::size_t operator()(const pair<string, string>& seq) const
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
        for (std::size_t i = 0; i < seq.first.size() && i < max_width; ++i) {
            char nucleotide = toupper(seq.first[i]);
            if (nucleotide_codes.count(nucleotide)) {
                hash_value = (hash_value << 2) | nucleotide_codes[nucleotide];
            }
        }
        for (std::size_t j = 0; j < seq.second.size() && j < max_width; ++j) {
            char nucleotide = toupper(seq.second[j]);
            if (nucleotide_codes.count(nucleotide)) {
                hash_value = (hash_value << 2) | nucleotide_codes[nucleotide];
            }

        }
        return hash_value;
        // END your code above
    }
};

typedef unordered_multimap< pair<string,string>, size_t, DNAHasherKd> CSeqHashKd;
CSeqHashKd create_hash_table_kd(const vector< pair<string, string> >& kdmers)
// create one hash table by inserting both the prefix and suffix of each
//   k-mer. The prefix and suffix is the key. Associated with each element
//   in the hash table is the node id for that prefix or suffix in the
//   de Bruijn graph to be constructed.
{
    CSeqHashKd ht;
    size_t node_id = 0; // the node id will be used in the de Bruijn graph
    for (auto i = 0u; i < kdmers.size(); ++i) {
        for (auto j = 0u; j < 2; ++j) { // j=0: prefix; j=1: suffix
            auto key = kdmers[i].first.substr(j, kdmers[i].first.length() - 1);
            auto key2 = kdmers[i].second.substr(j, kdmers[i].second.length() - 1);
            auto keys = make_pair(key, key2);
            if (ht.find(keys) == ht.end()) {
                ht.insert(make_pair(keys, node_id++));
            }
        }
    }
    return ht;
}

void create_deBruijn_graph_by_hashing_kd(const vector<pair<string,string> >& kdmers, DiGraph_kd& g)
// create a deBruijn graph by inserting all k-mers into the graph by hashing
{
    // TO DO:

   // BEGIN your code below:
   // Create hash table for both the k-1 prefix and suffix of each k-mer
    CSeqHashKd node_ids = create_hash_table_kd(kdmers);

    // Initialize empty node vector for graph g
    g.m_nodes.clear();

    for (auto i = 0u; i < kdmers.size(); ++i) {
        // find the prefix node id from_id from the hash table
        auto prefix = make_pair((kdmers[i].first.substr(0, kdmers[i].first.length() - 1)), (kdmers[i].second.substr(0, kdmers[i].second.length() - 1)));
        auto range = node_ids.equal_range(prefix);
        auto from_id = range.first->second;

        // update node from_id's label to prefix if necessary
        if (g.m_nodes.size() <= from_id) {
            g.m_nodes.resize(from_id + 1);
        }
        auto& node_from = g.m_nodes[from_id];
        if (node_from.m_label.first.empty() || node_from.m_label.second.empty()) {
            node_from.m_label = prefix;
        }

        // find the suffix node id to_id from the hash table
        auto suffix = make_pair((kdmers[i].first.substr(1, kdmers[i].first.length() - 1)), (kdmers[i].second.substr(1, kdmers[i].second.length() - 1)));
        range = node_ids.equal_range(suffix);
        auto to_id = range.first->second;

        // update node to_id's label to suffix if necessary
        if (g.m_nodes.size() <= to_id) {
            g.m_nodes.resize(to_id + 1);
        }
        auto& node_to = g.m_nodes[to_id];
        if (node_from.m_label.first.empty() || node_from.m_label.second.empty()) {
            node_to.m_label = suffix;
        }
        // Create a new edge (from_id, to_id) by inserting node to_id into the adjacency list of node from_id
        g.m_nodes[from_id].m_outgoing.push_back(to_id);

        // Update the number of incoming edges of node to_id
        g.m_nodes[to_id].m_num_of_incoming++;
    }

}


string assemble_kdmers(const vector< pair<string,string> >& kdmers, const string& method, const string& dotfile, const size_t k, const size_t d) {
    string seq;
    DiGraph_kd g;

    // Create de Bruijn graph from k-dmers
    create_deBruijn_graph_by_hashing_kd(kdmers, g);

    if (!dotfile.empty()) {
        printDOTFileKd(g, dotfile);
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


string build_kd_sequence(const list<size_t>& path, const DiGraph_kd& g, const size_t k, const size_t d)
{
    vector<Node_kd> nodes = g.m_nodes;

    string seq,next;

    seq.resize(k - 1u + path.size() - 1u);
    next.resize(k - 1u + path.size() - 1u);

    seq.replace(0, k - 1u, nodes[path.front()].m_label.first);
    next.replace(0, k - 1u, nodes[path.front()].m_label.second);
    auto pos = path.begin();
    pos++;

    size_t i = k - 1;

    while (pos != path.end()){

        seq[i] = nodes[*pos].m_label.first.back();
        next[i] = nodes[*pos].m_label.second.back();
        i++;
        pos++;
    }
    seq += next.substr(next.length() - (k + d));
    cout << seq << endl;
    return seq;
}

void printDOTFileKd(const DiGraph_kd& g, const string& file)
{
    ofstream ofs(file);

    ofs << "digraph for (K,d)-mer{" << endl;
    ofs << "label=\"de Bruijn graph\"" << endl;

    auto& nodes = g.m_nodes;

    for (auto& node : nodes) {
        for (auto to = node.m_outgoing.begin(); to != node.m_outgoing.end(); ++to) {
            auto& prefix1 = node.m_label.first;
            auto& prefix2 = node.m_label.second;
            auto& suffix1 = nodes[*to].m_label.first;
            auto& suffix2 = nodes[*to].m_label.second;
            ofs << prefix1 << "->" << suffix1
                << "[label=" << prefix1 << suffix1.back() << "];" << "---" << prefix2 << "->" << suffix2
                << "[label=" << prefix2 << suffix2.back() << "];" << endl;
        }
    }

    ofs << "}" << endl;
}



size_t source(const DiGraph_kd& g)
// Find a source node from g: the node has one more
// outgoing edge than incoming edge. When such a node
// does not exist, a value of the total number of nodes
// is returned.
{
    size_t i;

    for (i = 0; i < g.m_nodes.size(); ++i) {
        if (g.m_nodes[i].m_outgoing.size()
            == g.m_nodes[i].m_num_of_incoming + 1) {
            break;
        }
    }
    return i;
}

size_t sink(const DiGraph_kd& g)
// Find a source node from g: the node has one more
// outgoing edge than incoming edge. When such a node
// does not exist, a value of the total number of nodes
// is returned.
{
    size_t i;
    for (i = 0; i < g.m_nodes.size(); ++i) {
        if (g.m_nodes[i].m_outgoing.size() + 1
            == g.m_nodes[i].m_num_of_incoming) {
            break;
        }
    }
    return i;
}


list<size_t> find_Eulerian_cycle(DiGraph_kd& g)
// find an Eulerian cycle from graph g
{
    list <size_t> cycle; // main cycle

    // BEGIN your code here:
    stack<size_t> stk;
    vector<bool> used(g.m_nodes.size(), false);
    size_t start = 0;

    // Find a starting node with outgoing edges
    for (size_t i = 0; i < g.m_nodes.size(); ++i) {
        if (!g.m_nodes[i].m_outgoing.empty()) {
            start = i;
            break;
        }
    }

    stk.push(start);

    while (!stk.empty()) {
        size_t node = stk.top();
        if (g.m_nodes[node].m_outgoing.empty()) {
            cycle.push_front(node);
            stk.pop();
        }
        else {
            size_t next = g.m_nodes[node].m_outgoing.back();
            g.m_nodes[node].m_outgoing.pop_back();
            stk.push(next);
        }
    }

    // Check if all edges have been used
    for (size_t i = 0; i < g.m_nodes.size(); ++i) {
        if (!g.m_nodes[i].m_outgoing.empty()) {
            throw "Graph is not Eulerian!";
        }
    }

    return cycle;
}

list<size_t> find_Eulerian_path(DiGraph_kd& g)
// find an Eulerian path from graph g, assuming g has such a path
{
    list <size_t> path, cycle;

    size_t src = source(g);    // find the source node
    size_t dest = sink(g);     // find the sink nod e

    // In the special case graph with only cycles
    //   and no source nor sink, we choose node 0 as the
    //   start and end of the Eulerian path:
    src = src >= g.m_nodes.size() ? 0 : src;
    dest = dest >= g.m_nodes.size() ? 0 : dest;

    vector<Node_kd>& nodes = g.m_nodes;

    // add an edge from the sink node to the source node
    nodes[dest].m_outgoing.push_back(src);

    // increase the incoming degree of the source node by one
    nodes[src].m_num_of_incoming++;

    cycle = find_Eulerian_cycle(g);

    list<size_t>::iterator pos_src, pos_dest;

    for (pos_dest = cycle.begin(); pos_dest != cycle.end(); ++pos_dest) {
        pos_src = pos_dest;
        pos_src++;
        if (pos_src != cycle.end()) {
            if (*pos_src == src && *pos_dest == dest) {
                break;
            }
        }
        else {
            break;
        }
    }

    if (pos_src != cycle.end() && pos_dest != cycle.end()) {

        auto pos = pos_src;
        do {
            path.push_back(*pos);
            pos++;
            if (pos == cycle.end()) {
                pos = cycle.begin();
                pos++; // skip the first node in the cycle
                // which is the same with the last
                // last node in the cycle.
            }
        } while (pos != pos_src);

    }
    else {
        throw "Searching for Eulerian path has failed!";
    }

    // return the path
    return path;

}

bool has_Eulerian_path(const DiGraph_kd& g)
// determine if graph g has an Eulerian path. This path could
//   be a cycle in special cases.
{
    bool exist = true;

    size_t numSources = 0;
    size_t numSinks = 0;

    for (auto& node : g.m_nodes) {
        size_t out = node.m_outgoing.size();
        size_t in = node.m_num_of_incoming;
        if (out == in) { // check for intermediate balanced node
            continue;
        }
        else if (out == in + 1) { // check for source node
            numSources++;
            if (numSources > 1) {
                exist = false;
                break;
            }
        }
        else if (out + 1 == in) { // check for sink node
            numSinks++;
            if (numSinks > 1) {
                exist = false;
                break;
            }
        }
        else {
            exist = false;
            break;
        }
    }

    return exist;
}


