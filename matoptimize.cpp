#include <matoptimize.h>

void matOptimize(IQTree *tree, Alignment *alignment, string name,int sprDist)
{
    cout << "\n========== Start matOptimize ========== \n";
    PhyloNode *startNode = tree->findNode(name);
    assert(startNode != NULL);
    cout << "startNode: " << startNode->name << endl;
    
    vector<PhyloNode *> candidateNodes; 
    vector<PhyloNode *> sourceNodes;
    sourceNodes.push_back(startNode);
    while(true) {   
        // Init tree's infomartion.
        tree->depth_first_search();

        // Find candidate nodes.
        candidateNodes.clear();
        tree->initCandidate(sourceNodes, sprDist);

        cout << "init candidates done\n";
        tree->tryMatOptimize(sprDist);
        cout << "try mat optimize done\n";
        break;
    }
}