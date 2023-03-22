//
//  kd-assembler.hpp
//  kd-assembler
//
//  Created by Joe Song on 3/19/18.
//  Copyright © 2018 Joe Song. All rights reserved.
/*

#ifndef kd_assembler_h
#define kd_assembler_h

#include <list>
#include <vector>
#include <string>

using namespace std;

void create_deBruijn_graph_by_string_comp(const vector<string> & kmers, DiGraph & g);

void create_deBruijn_graph_by_hashing(const vector<string> & kmers, DiGraph & g);

list<size_t> find_Eulerian_path(DiGraph & g);

bool has_Eulerian_path(const DiGraph & g);

vector<string> generate_kdmers(const string& read, size_t k, size_t d);

string build_kd_sequence(const list<size_t> & path, const DiGraph & g, const size_t k, const size_t d);

string assemble_kdmers(const vector<string> & kmers, size_t k, size_t d, const string & dotfile="");

void printDOTFile(const DiGraph & g, const string & file);

#endif /* k_assembler_h */
 