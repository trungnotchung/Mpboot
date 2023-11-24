#ifndef MATOPTIMIZE_H
#define MATOPTIMIZE_H

#include <phylotree.h>
#include <iqtree.h>
#include <mutation.h>

// SPR with mutation
void matOptimize(IQTree *tree, Alignment *alignment, string name,int sprDist);
#endif