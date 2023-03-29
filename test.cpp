//
//  test.cpp
//  k-assembler
//
//  Created by Joe Song on 11/24/15.
//  Copyright © 2015 Joe Song. All rights reserved.
//
//  Updated 3/19/2018
//  Updated 9/19/2021

#include "sequence.hpp"
#include "k-assembler.hpp"
#include "kd-assembler.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

void test_and_print_message(const string & seq, const string & seq_truth,
                            size_t k, const string & message)
{
    //cout << seq << endl;
    if(seq == seq_truth) {
        cout << "Passed " << message <<
            " (assembled original sequence). Congratulations!" << endl;
    } else if(compare_composition(seq, seq_truth, k)) {
        cout << "Passed " << message <<
            " (assembled a sequence of the same composition with the original sequence). Congratulations!" << endl;
    } else {
        cerr << "FAILED " + message + "!" << endl;
    }
}

static void test_1(const string & method)
{
    cout << "Testing k-assembler by " << method << endl;
    
    // Testing sequences:
    vector<string> seqs_truth = {
        "aaaaaaaaaaa",
        "agcagctcagc",
        "agcagctcagg",
        random_DNA_sequence(10000, 20000)
    };
    
    // The value of k for k-mers to be used for each test sequence:
    vector<size_t> ks = {5, 3, 3, 20};
    
    for (size_t i=0; i<seqs_truth.size(); ++i) {
             
        cout << endl << "Example " << i << ":" << endl;
        
        string & seq_truth = seqs_truth[i];
        size_t & k = ks[i];
        
        auto kmers = get_kmers(seq_truth, k, true);
        
        DiGraph g;
        clock_t begin = clock();
        
        if(method == "k-mer pairwise comparison") {
            create_deBruijn_graph_by_string_comp(kmers, g);
        } else if(method == "k-mer hashing") {
            create_deBruijn_graph_by_hashing(kmers, g);
        } else {
            cerr << "ERROR: unknown methods!" << endl;
        }
        
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        
        cout << "Elapsed time for building de Bruijn graph: " << elapsed_secs << endl;
        
        if(has_Eulerian_path(g)) {
            std::cout << "Passed test for existence of Eulerian path. Congratulations!" << endl;
        } else {
            std::cout << "Failed test for existence of Eulerian path!" << endl;
        }
        
        try {
            list<size_t> path = find_Eulerian_path(g);
            /*
            for (auto const& i : path) {
                cout << i << endl;
            }
            */
            string seq = build_sequence(path, g);
            
            stringstream ss;
            ss << "Test 1 Example " << i;
            
            test_and_print_message(seq, seq_truth, k, ss.str());
            
        } catch (string err_message) {
            cerr << "ERROR: " << err_message << endl;
        }
    }
}

static void test_2(const string & method)
{
    string seq_truth = random_DNA_sequence();
    size_t k=10;
    
    vector<string> kmers = get_kmers(seq_truth, k);
    
    try {
        string seq = assemble_kmers(kmers, method);
    
        test_and_print_message(seq, seq_truth, k, "Test 2");

    } catch (char const * s) {
        cerr << s << endl;
    }
}

void test_3(const string & method)
{
    string seq_truth = random_DNA_sequence(15, 15);
    size_t k=4;

    cout << "Sequence: " << seq_truth << endl;
    
    vector<string> kmers = get_kmers(seq_truth, k);
    cout << "kmers:" << endl;
    for (auto & kmer : kmers) {
        cout << kmer << endl;
    }
    
    try {
        string seq = assemble_kmers(kmers, method, "deBruijn.dot");
        test_and_print_message(seq, seq_truth, k, "Test 3");
        
    } catch (char const * s) {
        cerr << s << endl;
    }
}


void test_4(const string& method)
{
    string seq_truth = "aaaaagcctttaatt";
    size_t k = 3;
    size_t d = 2;
    cout << "Sequence: " << seq_truth << endl;

    vector<pair <string,string> > kdmers = get_kdmers(seq_truth, k, d);
    cout << "kdmers:" << endl;
    for (auto& kdmer : kdmers) {
        cout << kdmer.first << "-" << kdmer.second << endl;
    }
    try {
        string seq2 = assemble_kdmers(kdmers, method, "deBruijn2.dot", k, d);
        test_and_print_message(seq2, seq_truth, k, "Test 4");

    }
    catch (char const* s) {
        cerr << s << endl;
    }
}


void test_seq_assembly()
{
    vector<string> methods = {
        "k-mer pairwise comparison",
        "k-mer hashing"
    };
    
    for (auto method: methods) {
        
        cout << "-----------" << endl;
        
        test_1(method);
        cout << endl;
        test_2(method);
        cout << endl;

    }

    cout << "-----------" << endl;

    test_3("k-mer hashing");
    cout << "-----------" << endl;
    test_4("k-mer hashing");
}
