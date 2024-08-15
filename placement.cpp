#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iqtree_config.h>
#include "phylotree.h"
#include "phylosupertree.h"
#include "phylosupertreeplen.h"
#include "phyloanalysis.h"
#include "alignment.h"
#include "superalignment.h"
#include "iqtree.h"
#include "model/modelgtr.h"
#include "model/modeldna.h"
#include "myreader.h"
#include "model/rateheterogeneity.h"
#include "model/rategamma.h"
#include "model/rateinvar.h"
#include "model/rategammainvar.h"
// #include "modeltest_wrapper.h"
#include "model/modelprotein.h"
#include "model/modelbin.h"
#include "model/modelcodon.h"
#include "stoprule.h"

#include "mtreeset.h"
#include "mexttree.h"
#include "model/ratemeyerhaeseler.h"
#include "whtest_wrapper.h"
#include "model/partitionmodel.h"
#include "guidedbootstrap.h"
#include "model/modelset.h"
#include "timeutil.h"
#include "parstree.h"
#include "tinatree.h"
#include "sprparsimony.h"
#include "placement.h"
#include <algorithm>
#include <thread>

void checkCorectTree(char *originTreeFile, char *newTreeFile)
{
	cout << "================= Check correct tree ================\n";
	IQTree *originTree = new IQTree;
	bool originIsRooted = false;
	originTree->readTree(originTreeFile, originIsRooted);

	IQTree *newTree = new IQTree;
	bool newIsRooted = false;
	newTree->readTree(newTreeFile, newIsRooted);

	vector<string> originLeafName;
	originTree->getLeafName(originLeafName);

	newTree->assignRoot(originLeafName[0]);
	sort(originLeafName.begin(), originLeafName.end());
	newTree->initInfoNode(originLeafName);

	if (newTree->compareTree(originTree))
		cout << "Correct tree\n";
	else
		cout << "Wrong tree\n";

	delete originTree;
	delete newTree;
}

void configLeafNames(IQTree *tree, Node *node, Node *dad)
{
	if (node->isLeaf())
	{
		node->id = tree->aln->getSeqID(node->name);
	}
	FOR_NEIGHBOR_IT(node, dad, it)
	configLeafNames(tree, (*it)->node, node);
}

string ppRunOriginalTbr(Alignment *alignment, Params &params, string newickTree = "")
{
	cout << "\n========== Start TBR core ==========\n";
	IQTree *tree = new IQTree(alignment);
	tree->params = &params;
	if (newickTree == "")
	{
		bool isRooted = params.is_rooted;
		tree->readTree(params.mutation_tree_file, isRooted);
	}
	else
	{
		tree->readTreeString(newickTree);
	}

	ofstream fout1("tree_tbr_1.txt");
	tree->drawTree(fout1, WT_SORT_TAXA | WT_NEWLINE);

	// becasue the tree is read from file and file contains its name, not its id
	configLeafNames(tree, tree->root, NULL);

	tree->initializeAllPartialPars();
	tree->clearAllPartialLH();
	tree->curScore = tree->computeParsimony();
	cout << "Tree's score before running tbr: " << tree->curScore << '\n';

	double start_time = getCPUTime();
	string newTreeString = tree->ppRunOriginalTbr();
	ofstream fout2("tree_tbr_2.txt");
	tree->drawTree(fout2, WT_SORT_TAXA | WT_NEWLINE);
	double end_time = getCPUTime();
	cout << "Time running TBR: " << fixed << setprecision(3) << (double)(end_time - start_time) << " seconds\n";
	cout << "Memory: " << getMemory() << " KB\n";

	if (params.pp_test_optimize)
	{
		ofstream fout("newTree.txt");
		tree->printTree(fout, WT_SORT_TAXA | WT_NEWLINE);
		fout.close();
		checkCorectTree(params.original_tree_file, "newTree.txt");
	}

	return newTreeString;
}

string ppRunOriginalSpr(Alignment *alignment, Params &params, string newickTree = "")
{
	cout << "\n========== Start spr core ==========\n";
	IQTree *tree = new IQTree(alignment);
	tree->params = &params;
	if (newickTree == "")
	{
		bool isRooted = params.is_rooted;
		tree->readTree(params.mutation_tree_file, isRooted);
	}
	else
	{
		tree->readTreeString(newickTree);
	}

	ofstream fout1("tree1.txt");
	tree->drawTree(fout1, WT_SORT_TAXA | WT_NEWLINE);
	configLeafNames(tree, tree->root, NULL);

	tree->initializeAllPartialPars();
	tree->clearAllPartialLH();
	tree->curScore = tree->computeParsimony();
	cout << "Tree's score before running spr: " << tree->curScore << '\n';

	double startTime = getCPUTime();
	string newTreeString = tree->ppRunOriginalSpr();
	ofstream fout2("tree2.txt");
	tree->drawTree(fout2, WT_SORT_TAXA | WT_NEWLINE);
	double endTime = getCPUTime();
	cout << "Time running SPR: " << fixed << setprecision(3) << (double)(endTime - startTime) << " seconds\n";
	cout << "Memory: " << getMemory() << " KB\n";

	if (params.pp_test_optimize)
	{
		ofstream fout("newTree.txt");
		tree->printTree(fout, WT_SORT_TAXA | WT_NEWLINE);
		fout.close();
		checkCorectTree(params.original_tree_file, "newTree.txt");
	}
	return newTreeString;
}

void initialize(IQTree *tree, Alignment *alignment, vector<int> &savePermCol, vector<int> &permCol, vector<int> &compressedPermCol)
{
	permCol.resize(savePermCol.size());
	compressedPermCol.resize(savePermCol.size());
	if (alignment->existingSamples.size())
	{
		for (int j = 0; j < savePermCol.size(); ++j)
		{
			int p = savePermCol[j];
			compressedPermCol[j] = alignment->existingSamples[0][p].compressed_position;
			permCol[j] = alignment->existingSamples[0][p].position;
		}
	}
	alignment->ungroupSitePattern();
	tree->add_row = true;
	tree->save_branch_states_dad = new UINT[(alignment->size() + 7) / 8 + 1];
	tree->computeParsimony();
	tree->initMutation(permCol, compressedPermCol);
}

int readFile(ifstream &inFileStream, char *outFileName, int numRow)
{
	ofstream outFile(outFileName);
	if (!outFile.is_open())
	{
		cout << "Cannot open file " << outFileName << '\n';
		return 0;
	}
	string line;
	int curRow = 0;
	while (getline(inFileStream, line))
	{
		if (line == "")
		{
			continue;
		}
		outFile << line << '\n';
		++curRow;
		if (curRow >= numRow)
		{
			break;
		}
	}
	outFile.close();
	return curRow;
}

int readVCFFile(IQTree *tree, Alignment **alignment, Params &params)
{
	char *alnFile = params.aln_file;
	ifstream in;
	in.exceptions(ios::failbit | ios::badbit);
	in.open(alnFile);
	string line;
	in.exceptions(ios::badbit);

	int totalColumn = readFile(in, "temp.vcf", 12) - 1;
	*alignment = new Alignment("temp.vcf", params.sequence_type, params.intype, params.numStartRow);
	(*alignment)->ungroupSitePattern();
	std::remove("temp.vcf");

	tree->setAlignment(*alignment);
	tree->aln = *alignment;

	vector<int> permCol = (*alignment)->findPermCol();
	vector<int> savePermCol = permCol;
	vector<int> compressedPermCol = permCol;
	initialize(tree, *alignment, savePermCol, permCol, compressedPermCol);

	auto startTime = getCPUTime();
	while (true)
	{
		startTime = getCPUTime();
		int numColumn = (*alignment)->readPartialVCF(in, params.sequence_type, savePermCol, params.numStartRow, totalColumn, 8);
		if (numColumn == 0)
		{
			break;
		}

		tree->clearAllPartialLH();
		totalColumn += numColumn;
		startTime = getCPUTime();
		initialize(tree, *alignment, savePermCol, permCol, compressedPermCol);
	}

	in.close();
	return totalColumn;
}

CandidateNode calcBestBranchPlaceNode(vector<MutationNode> &missingSamples, IQTree *tree, vector<int> &visited_missing_sample_mutations, vector<Mutation> &cur_missing_sample_mutations, vector<int> &visited_ancestral_mutations, vector<Mutation> &cur_ancestral_mutations, vector<int> &visited_excess_mutations, vector<Mutation> &cur_excess_mutations, int i, int &timerRegular, int &timerOptimized)
{
	int totalNodes = (int)tree->bfs.size();
	CandidateNode inp;
	int bestSetDifference = INF;
	size_t bestNodeNumLeaves = INF;
	size_t bestDistance = INF;
	std::vector<Mutation> excessMutations;
	std::vector<bool> nodeHasUnique(totalNodes, false);
	bool bestNodeHasUnique = false;
	size_t bestJ = 0;

	inp.best_node = NULL;
	inp.best_node_branch = NULL;
	inp.best_set_difference = &bestSetDifference;
	inp.best_node_num_leaves = &bestNodeNumLeaves;
	inp.best_distance = &bestDistance;
	inp.node = (PhyloNode *)tree->root->neighbors[0]->node;
	inp.node_branch = (PhyloNeighbor *)inp.node->findNeighbor(tree->root);
	inp.missing_sample_mutations = &missingSamples[i].mutations;
	inp.excess_mutations = &excessMutations;
	inp.has_unique = &bestNodeHasUnique;
	inp.node_has_unique = &(nodeHasUnique);
	inp.best_j = &bestJ;

	tree->initDataCalculatePlacementMutation(inp, visited_missing_sample_mutations, cur_missing_sample_mutations, timerOptimized);
	tree->optimizedCalculatePlacementMutation(inp, 0, true, visited_missing_sample_mutations, cur_missing_sample_mutations, visited_ancestral_mutations, cur_ancestral_mutations, visited_excess_mutations, cur_excess_mutations, timerOptimized);

	for (int j = 0; j < totalNodes; ++j)
	{
		if (inp.best_node == tree->bfs[j].first)
		{
			bestJ = j;
			break;
		}
	}
	*inp.best_set_difference = INF;
	inp.j = bestJ;
	inp.node = tree->bfs[bestJ].first;
	inp.node_branch = tree->bfs[bestJ].second;
	return inp;
}

void placeMutation(vector<MutationNode> &missingSamples, IQTree *tree, vector<int> &visited_missing_sample_mutations, vector<Mutation> &cur_missing_sample_mutations, vector<int> &visited_ancestral_mutations, vector<Mutation> &cur_ancestral_mutations, vector<int> &visited_excess_mutations, vector<Mutation> &cur_excess_mutations, int i, int &timerRegular, int &timerOptimized, CandidateNode &inp)
{
	std::vector<Mutation> excessMutations;
	inp.excess_mutations = &excessMutations;
	tree->calculatePlacementMutation(inp, false, true, visited_missing_sample_mutations, cur_missing_sample_mutations, visited_ancestral_mutations, cur_ancestral_mutations, timerRegular);
	tree->addNewSample(inp.node, inp.node_branch, *inp.excess_mutations, i, missingSamples[i].name, visited_ancestral_mutations, cur_ancestral_mutations, visited_excess_mutations, cur_excess_mutations, timerRegular);
}

void addMoreRowMutation(Params &params)
{
	Alignment *alignment;

	IQTree *tree;
	tree = new IQTree;

	char *fileName = params.mutation_tree_file;
	bool isRooted = false;

	if (params.tree_zip_file != NULL)
	{
		tree->readTree(params.tree_zip_file, fileName, isRooted);
	}
	else
	{
		tree->readTree(fileName, isRooted);
	}

	int vecSize = readVCFFile(tree, &alignment, params) + 1;
	const int NUM_THREADS = 3;

	std::vector<std::vector<Mutation>> cur_excess_mutations, cur_missing_sample_mutations, cur_ancestral_mutations;
	std::vector<std::vector<int>> visited_missing_sample_mutations, visited_ancestral_mutations;
	std::vector<std::vector<int>> visited_excess_mutations;
	std::vector<int> timerRegular, timerOptimized;

	vector<std::future<CandidateNode>> result;

	// Init new tree's memory
	cur_missing_sample_mutations.resize(NUM_THREADS, std::vector<Mutation>(vecSize));
	cur_ancestral_mutations.resize(NUM_THREADS, std::vector<Mutation>(vecSize));
	visited_missing_sample_mutations.resize(NUM_THREADS, std::vector<int>(vecSize));
	visited_ancestral_mutations.resize(NUM_THREADS, std::vector<int>(vecSize));
	cur_excess_mutations.resize(NUM_THREADS, std::vector<Mutation>(vecSize));
	visited_excess_mutations.resize(NUM_THREADS, std::vector<int>(vecSize));
	timerRegular.resize(NUM_THREADS, 0);
	timerOptimized.resize(NUM_THREADS, 0);

	cout << "\n========== Start placement core ==========\n";

	// free memory
	delete[] tree->save_branch_states_dad;
	tree->add_row = false;

	cout << "Tree parsimony after init mutations: " << tree->computeParsimonyScoreMutation() << '\n';
	int numSample = (int)alignment->missingSamples.size();
	vector<MutationNode> missingSamples(numSample);
	for (int i = 0; i < (int)alignment->missingSamples.size(); ++i)
	{
		missingSamples[i].mutations = alignment->missingSamples[i];
		missingSamples[i].name = alignment->remainName[i];
	}
	numSample = min(numSample, params.numAddRow);

	auto startTime = getCPUTime();

	for (int k = 0; k < numSample; k += NUM_THREADS)
	{
		result.clear();
		tree->breadth_first_expansion();
		for (int i = k; i < min(numSample, k + NUM_THREADS); ++i)
		{
			result.push_back(std::async(std::launch::async, calcBestBranchPlaceNode, std::ref(missingSamples), std::ref(tree),
										std::ref(visited_missing_sample_mutations[i % NUM_THREADS]),
										std::ref(cur_missing_sample_mutations[i % NUM_THREADS]),
										std::ref(visited_ancestral_mutations[i % NUM_THREADS]),
										std::ref(cur_ancestral_mutations[i % NUM_THREADS]),
										std::ref(visited_excess_mutations[i % NUM_THREADS]),
										std::ref(cur_excess_mutations[i % NUM_THREADS]), i,
										std::ref(timerRegular[i % NUM_THREADS]),
										std::ref(timerOptimized[i % NUM_THREADS])));
		}
		for (int i = k; i < min(numSample, k + NUM_THREADS); ++i)
		{
			CandidateNode inp = result[i % NUM_THREADS].get();
			placeMutation(missingSamples, tree, visited_missing_sample_mutations[i % NUM_THREADS], cur_missing_sample_mutations[i % NUM_THREADS], visited_ancestral_mutations[i % NUM_THREADS], cur_ancestral_mutations[i % NUM_THREADS], visited_excess_mutations[i % NUM_THREADS], cur_excess_mutations[i % NUM_THREADS], i, timerRegular[i % NUM_THREADS], timerOptimized[i % NUM_THREADS], inp);
		}
	}
	cout << "New tree's parsimony score: " << tree->computeParsimonyScoreMutation() << '\n';
	cout << "Time: " << fixed << setprecision(3) << (double)(getCPUTime() - startTime) << " seconds\n";
	cout << "Memory: " << getMemory() << " KB\n";

	// free memory
	cur_missing_sample_mutations.clear();
	cur_ancestral_mutations.clear();
	visited_missing_sample_mutations.clear();
	visited_ancestral_mutations.clear();
	cur_excess_mutations.clear();
	visited_excess_mutations.clear();

	delete alignment;
	alignment = NULL;
	delete tree;
}