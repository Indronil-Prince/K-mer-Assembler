//
//  EulerPath.cpp
//  k-assembler
//
//  Created by Joe Song on 11/23/15.
//  Copyright © 2015-2021 Joe Song. All rights reserved.
//
//  Updated 3/19/2018
//  Updated 9/19/2021

#include "k-assembler.hpp"
#include <stack>
#include <unordered_map>

using namespace std;

size_t source(const DiGraph & g)
// Find a source node from g: the node has one more
// outgoing edge than incoming edge. When such a node
// does not exist, a value of the total number of nodes
// is returned.
{
    size_t i;
    
    for (i=0; i<g.m_nodes.size(); ++i) {
        if (g.m_nodes[i].m_outgoing.size()
            == g.m_nodes[i].m_num_of_incoming + 1) {
            break;
        }
    }
    return i;
}

size_t sink(const DiGraph & g)
// Find a source node from g: the node has one more
// outgoing edge than incoming edge. When such a node
// does not exist, a value of the total number of nodes
// is returned.
{
    size_t i;
    for (i=0; i<g.m_nodes.size(); ++i) {
        if (g.m_nodes[i].m_outgoing.size() + 1
            == g.m_nodes[i].m_num_of_incoming) {
            break;
        }
    }
    return i;
}


list<size_t> find_Eulerian_cycle(DiGraph & g)
// find an Eulerian cycle from graph g
{
    list <size_t> cycle; // main cycle
    
    // TO-DO: insert code to find Eulerian cycle represented
    //   as a list of node id.
    // E.g., 3 -> 2 -> 0 -> 3 -> 4 -> 1 -> 3 is a cycle on
    //   a graph with 5 nodes.
    
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
        } else {
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

list<size_t> find_Eulerian_path(DiGraph & g)
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

    vector<Node> & nodes = g.m_nodes;

    // add an edge from the sink node to the source node
    nodes[dest].m_outgoing.push_back(src);

    // increase the incoming degree of the source node by one
    nodes[src].m_num_of_incoming ++;

    cycle = find_Eulerian_cycle(g);

    list<size_t>::iterator pos_src, pos_dest;

    for (pos_dest = cycle.begin(); pos_dest != cycle.end(); ++ pos_dest) {
        pos_src = pos_dest;
        pos_src ++;
        if(pos_src != cycle.end()) {
            if (*pos_src == src && *pos_dest == dest) {
                break;
            }
        } else {
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

bool has_Eulerian_path(const DiGraph & g)
// determine if graph g has an Eulerian path. This path could
//   be a cycle in special cases.
{
    bool exist = true;
    
    size_t numSources=0;
    size_t numSinks=0;
    
    for (auto & node : g.m_nodes) {
        size_t out = node.m_outgoing.size();
        size_t in = node.m_num_of_incoming;
        if(out == in) { // check for intermediate balanced node
            continue;
        } else if(out == in + 1) { // check for source node
            numSources ++;
            if (numSources > 1) {
                exist = false;
                break;
            }
        } else if(out + 1 == in) { // check for sink node
            numSinks ++;
            if (numSinks > 1) {
                exist = false;
                break;
            }
        } else {
            exist = false;
            break;
        }
    }
    
    return exist;
}