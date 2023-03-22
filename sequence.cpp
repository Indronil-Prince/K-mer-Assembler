//
//  sequence.cpp
//  k-assembler
//
//  Created by Joe Song on 12/14/15.
//  Copyright © 2015 Joe Song. All rights reserved.
//
//  Updated 3/19/2018
//  Updated 9/19/2021

#include "k-assembler.hpp"

#include <iostream>
#include <random>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <sstream>

using namespace std;

string random_DNA_sequence(size_t min_length=10, size_t max_length=10000)
// generate a DNA sequence of length in the given range
{
    string DNA;
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> length_dis(min_length, max_length);
    
    // random generate a length
    DNA.resize(length_dis(gen));
    
    // random generate the DNA sequence using alphabet atgc
    std::uniform_int_distribution<unsigned char> nucleotide_dis(0, 3);
    
    string atgc = "atgc";
    
    for (int n=0; n<DNA.size(); ++n) {
        DNA[n] = atgc[nucleotide_dis(gen)];
    }
    
    return DNA;
}

vector<string> get_kmers(string seq, size_t k, bool randomized=true)
// obtain all k-mers of a given sequence. The order of the k-mers
//   is randomized by default.
{
    vector<string> kmers(seq.length() - k + 1);
    
    for (size_t i=0; i<kmers.size(); ++i) {
        kmers[i] = seq.substr(i, k);
    }
    
    if (randomized) {
        vector<string> randomized_kmers(kmers.size());
        
        std::random_device rd;
        std::mt19937 gen(rd());
        
        size_t nkmers = kmers.size();
        
        for (size_t i=0; i<nkmers-1; ++i) {
            std::uniform_int_distribution<size_t> dis(i, kmers.size()-1);
            size_t j = dis(gen);
            string kmer = kmers[j];
            kmers[j] = kmers[i];
            kmers[i] = kmer;
        }
    }
    
    return kmers;
}

bool compare_composition(const string & s1, const string & s2, size_t k)
// compare the canonical composition of two sequences. Canonical
//   composition is all k-mers arranged in dictionary order.
{
    bool same_composition = true;
    
    if (s1.length() != s2.length()) {
        same_composition = false;
    } else if(s1 == s2) {
        same_composition = true;
    } else {
        auto composition1 = get_kmers(s1, k);
        sort(composition1.begin(), composition1.end());
        auto composition2 = get_kmers(s2, k);
        sort(composition2.begin(), composition2.end());
        same_composition = composition1 == composition2;
    }
    return same_composition;
}
