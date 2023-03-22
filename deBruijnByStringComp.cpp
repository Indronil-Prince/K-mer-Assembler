//
//  deBruijnByStringComp.cpp
//  k-assembler
//
//  Created by Joe Song on 12/14/15.
//  Copyright © 2015 Joe Song. All rights reserved.
//
//  Updated 3/19/2018

#include "k-assembler.hpp"

void create_deBruijn_graph_by_string_comp(const vector<string> & kmers, DiGraph & g)
// Insert all k mers into graph g by pair-wise sequence comparison
{
    list<Node> nodes;
    
    for (auto & kmer: kmers) { // for each kmer of size k
        
        size_t k = kmer.size();
        
        // The "from" node:
        string prefix = kmer.substr(0, k-1);
        
        // find for a node that matches the k-1 prefix
        size_t i=0; // i is the node id
        
        auto itr=nodes.begin();
        for (;
             itr!=nodes.end(); ++itr, ++i) {
            if (itr->m_label == prefix) {
                break;
            }
        }
        
        if(i >= nodes.size()) {
            // if the k-1 prefix does not exist as a node
            //   in the graph, create a new node
            
            Node from;
            from.m_label = prefix;
            from.m_num_of_incoming = 0;
            nodes.push_back(from);
            itr = nodes.end();
            itr --;
        }
        
        auto from_itr = itr;
        
        // The "to" node:
        string suffix = kmer.substr(1, k-1);
        
        // find the k-1 suffix
        size_t j=0;
        itr = nodes.begin();
        for (; itr!=nodes.end(); ++itr, ++j) {
            if (itr->m_label == suffix) {
                break;
            }
        }
        
        if(j >= nodes.size()) {
            // if the k-1 suffix does not exist as a node
            //   in the graph, create a new node
            Node to;
            to.m_label = suffix;
            to.m_num_of_incoming = 0;
            
            // insert the new node to nodes
            nodes.push_back(to);
            
            // remember the new node position on the list
            itr = nodes.end();
            itr --;
        }
        
        auto to_itr = itr;
        
        from_itr -> m_outgoing.push_back(j);
        to_itr -> m_num_of_incoming ++;
    }
    
    // transfer the nodes from the list to a vector
    g.m_nodes.resize(nodes.size());
    auto litr=nodes.begin();
    auto vitr=g.m_nodes.begin();
    while(litr != nodes.end()) {
        * vitr ++ = * litr ++;
    }
}
