#ifndef PLACEMENT_H
#define PLACEMENT_H

#include "tools.h"
#include "mexttree.h"
#include "phylotesting.h"
#include "nnisearch.h"
#include "mutation.h"
#include "matoptimize.h"
#include "fstream"
#include <filesystem>

const int INF = (int)1e9 + 7;

int readFile(ifstream &inFileStream, char *outFileName, int numRow);

// read VCF file per 8 lines
int readVCFFile(IQTree *tree, Alignment *alignment, Params &params);

// add more K row using IQTree
int addMoreRowIQTree(IQTree *tree, Alignment *alignment);

// add more K row using PLL core 
int addMoreRowPLL(IQTree *tree, Alignment *alignment, Params &params);

// hill climbing to update candidates like mpboot
int updatePermutation(IQTree *tree, Alignment *alignment, Params &params, const vector<int> permCol, vector<vector<int> > &candidates, vector<int> &candidateScore);

// compute parsimony score after add more k row
int computeParsimonyPermutation(IQTree *tree, Alignment *alignment, Params &params, const vector<int> &permCol, const vector<int> &perm);

// add more K row using mutation like usher
void addMoreRowMutation(Params &params);

// check if origin tree doesn't change.
void checkCorrectTree(char *originTreeFile, char *newTreeFile);
#endif
