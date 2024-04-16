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

int addMoreRowIQTree(IQTree *tree, Alignment *alignment)
{
	auto startTime = getCPUTime();
	// add more K row to this tree
	int k = alignment->remainSeq.size();
	IQTree resTree;
	vector<int> perm;
	for (int i = 0; i < k; ++i)
		perm.push_back(i);
	vector<int> permCol = alignment->findPermCol();
	int bestScore = INF;
	int timer = 0;
	do
	{
		// my_random_shuffle(perm.begin(), perm.end());
		++timer;
		if (timer > 1)
			break;
		IQTree newTree;
		char *file_name = "tree.treefile";
		tree->printTree(file_name, WT_SORT_TAXA | WT_NEWLINE);
		bool is_rooted = false;
		(newTree).readTree(file_name, is_rooted);
		// newTree.copyPhyloTree(tree);
		newTree.setAlignment(alignment);
		newTree.aln = new Alignment;
		newTree.aln->copyAlignment(alignment);

		cout << "tree parsimony score before add more k rows: " << tree->computeParsimonyScore() << " " << newTree.computeParsimonyScore() << '\n';
		newTree.addRemainRow(alignment->remainName, alignment->remainSeq, perm, permCol);
		newTree.initializeAllPartialPars();
		newTree.clearAllPartialLH();
		cout << newTree.computeParsimonyScore() << " " << newTree.computeParsimony() << '\n';
		int curScore = newTree.computeParsimonyScore();
		if (curScore < bestScore)
		{
			bestScore = curScore;
		}
		delete newTree.aln;
		newTree.aln = NULL;
	} while (next_permutation(perm.begin(), perm.end()));
	cout << "\nBest tree parsimony found after add more k rows: ";
	cout << bestScore << '\n';
	cout << "Time: " << fixed << setprecision(3) << getCPUTime() - startTime << " seconds\n";
	return bestScore;
}

const double e = 2.7182818;
int addMoreRowPLL(IQTree *tree, Alignment *alignment, Params &params)
{
	auto startTime = getCPUTime();
	int k = alignment->remainSeq.size();
	IQTree resTree;
	vector<int> permCol = alignment->findPermCol();
	int bestScore;

	if (k <= 5)
	{
		int bestScore = INF;
		vector<int> bestPerm;

		vector<int> perm;
		for (int i = 0; i < k; ++i)
			perm.push_back(i);
		do
		{
			int newScore = computeParsimonyPermutation(tree, alignment, params, permCol, perm);
			if (bestScore > newScore)
				bestPerm = perm;
			bestScore = min(bestScore, newScore);
		} while (next_permutation(perm.begin(), perm.end()));
		cout << '\n';
		cout << "Best tree parsimony found after add more k rows using PLL core: ";
		cout << bestScore << '\n';
		cout << "Time: " << getCPUTime() - startTime << " seconds\n";
		return bestScore;
	}

	int numCandidate = 5;
	vector<vector<int>> candidates;
	vector<int> candidateScore;

	// init 100 permutation
	for (int id = 1; id <= 100; ++id)
	{
		vector<int> perm;
		for (int i = 0; i < k; ++i)
			perm.push_back(i);
		my_random_shuffle(perm.begin(), perm.end());
		int maxLoop = 10;
		int maxDist = max(1, k / 20);
		int maxLen = max(2, k / 20);
		int bestScoreHit = 1;
		int curScore = computeParsimonyPermutation(tree, alignment, params, permCol, perm);
		for (int loop = 1; loop <= maxLoop; ++loop)
		{
			// rand 2..maxLen
			int len = 2 + rand() % (maxLen - 1);
			// rand 0..k-len
			int l = rand() % (k - len);
			int r = l + len - 1;
			for (int i = 0; i < k; ++i)
			{
				if (l <= i && i <= r)
					continue;
				if (min(abs(i - l), abs(i - r)) > maxDist)
					continue;
				vector<int> newPerm;
				for (int j = 0; j <= i; ++j)
				{
					if (l <= j && j <= r)
						continue;
					newPerm.push_back(perm[j]);
				}
				for (int j = l; j <= r; ++j)
					newPerm.push_back(perm[j]);
				for (int j = i + 1; j < k; ++j)
				{
					if (l <= j && j <= r)
						continue;
					newPerm.push_back(perm[j]);
				}
				assert((int)newPerm.size() == k);
				int newScore = computeParsimonyPermutation(tree, alignment, params, permCol, newPerm);

				if (newScore < curScore)
				{
					curScore = newScore;
					perm = newPerm;
					bestScoreHit = 1;
				}
				else if (newScore == curScore)
				{
					if ((double)rand() / RAND_MAX >= 1.0 / bestScoreHit)
						perm = newPerm;
					++bestScoreHit;
				}
			}
		}

		candidates.push_back(perm);
		candidateScore.push_back(curScore);
		if ((int)candidates.size() > numCandidate)
		{
			int choice = 0;
			for (int i = 0; i < (int)candidates.size(); ++i)
			{
				if (candidateScore[i] > candidateScore[choice])
					choice = i;
				else if (candidateScore[i] == candidateScore[choice])
				{
					if (rand() % 2 == 1)
						choice = i;
				}
			}
			candidates.erase(candidates.begin() + choice);
			candidateScore.erase(candidateScore.begin() + choice);
		}
	}

	bestScore = updatePermutation(tree, alignment, params, permCol, candidates, candidateScore);
	cout << '\n';
	cout << "Best tree parsimony found after add more k rows using PLL core: ";
	cout << bestScore << '\n';
	cout << "Time: " << getCPUTime() - startTime << " seconds\n";
	return bestScore;
}

int updatePermutation(IQTree *tree, Alignment *alignment, Params &params, const vector<int> permCol, vector<vector<int>> &candidates, vector<int> &candidateScore)
{
	cout << "we are here\n";
	assert((int)candidates.size() > 0);
	// main loop
	int cnt = 0;
	while (true)
	{
		int id = rand() % ((int)candidates.size());
		vector<int> perm = candidates[id];
		int curScore = candidateScore[id];
		int k = (int)perm.size();
		int bestScoreHit = 1;
		int maxDist = max(1, k / 20);
		int maxLoop = 10;
		for (int loop = 1; loop <= maxLoop; ++loop)
		{
			for (int i = 0; i < k; ++i)
			{
				for (int j = i + 1; j <= min(k - 1, i + maxDist); ++j)
				{
					swap(perm[i], perm[j]);
					int newScore = computeParsimonyPermutation(tree, alignment, params, permCol, perm);
					if (newScore < curScore)
					{
						curScore = newScore;
						bestScoreHit = 1;
					}
					else if (newScore == curScore)
					{
						++bestScoreHit;
						if ((double)rand() / RAND_MAX < 1.0 / bestScoreHit) // not accept
							swap(perm[i], perm[j]);
					}
					else
						swap(perm[i], perm[j]);
				}
			}
		}

		if (curScore == candidateScore[id])
			++cnt;
		else
			cnt = 1;

		candidates.push_back(perm);
		candidateScore.push_back(curScore);
		// cout << curScore << " " << candidateScore[id] << " " << cnt << '\n';
		int choice = 0;
		for (int i = 0; i < (int)candidates.size(); ++i)
		{
			if (candidateScore[i] > candidateScore[choice])
				choice = i;
			else if (candidateScore[i] == candidateScore[choice])
			{
				if (rand() % 2 == 1)
					choice = i;
			}
		}
		candidates.erase(candidates.begin() + choice);
		candidateScore.erase(candidateScore.begin() + choice);

		if (cnt == 10)
			break;
	}

	int bestScore = candidateScore[0];
	for (int i = 0; i < (int)candidateScore.size(); ++i)
		bestScore = min(bestScore, candidateScore[i]);
	return bestScore;
}

int computeParsimonyPermutation(IQTree *tree, Alignment *alignment, Params &params, const vector<int> &permCol, const vector<int> &perm)
{
	IQTree newTree;
	(newTree).copyPhyloTree(tree);
	newTree.aln = new Alignment;
	newTree.aln->copyAlignment(alignment);
	int score = newTree.addRemainRowSPR(alignment->remainName, alignment->remainSeq, perm, permCol, params);
	delete newTree.aln;
	newTree.aln = NULL;
	return score;
}

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
		bool is_rooted = params.is_rooted;
		tree->readTree(params.mutation_tree_file, is_rooted);
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
		bool is_rooted = params.is_rooted;
		tree->readTree(params.mutation_tree_file, is_rooted);
	}
	else
	{
		tree->readTreeString(newickTree);
	}

	ofstream fout1("tree1.txt");
	tree->drawTree(fout1, WT_SORT_TAXA | WT_NEWLINE);

	// becasue the tree is read from file and file contains its name, not its id
	configLeafNames(tree, tree->root, NULL);

	tree->initializeAllPartialPars();
	tree->clearAllPartialLH();
	tree->curScore = tree->computeParsimony();
	cout << "Tree's score before running spr: " << tree->curScore << '\n';

	double start_time = getCPUTime();
	string newTreeString = tree->ppRunOriginalSpr();
	ofstream fout2("tree2.txt");
	tree->drawTree(fout2, WT_SORT_TAXA | WT_NEWLINE);
	double end_time = getCPUTime();
	cout << "Time running SPR: " << fixed << setprecision(3) << (double)(end_time - start_time) << " seconds\n";
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

void initialize(IQTree *tree, Alignment *alignment, vector<int> &savePermCol, vector<int> &permCol, vector<int> &compressedPermCol) {
	permCol.resize(savePermCol.size());
	compressedPermCol.resize(savePermCol.size());
	if (alignment->existingSamples.size())
	{
		for(int j = 0; j < savePermCol.size(); ++j) {
			int p = savePermCol[j];
			compressedPermCol[j] = alignment->existingSamples[0][p].compressed_position;
			permCol[j] = alignment->existingSamples[0][p].position;
		}
	}
	alignment->ungroupSitePattern();
	tree->add_row = true;
	tree->save_branch_states_dad = new UINT[(alignment->size() + 7) / 8 + 1];
	cout << "Tree parsimony before add k rows: " << tree->computeParsimony() << '\n';
	tree->initMutation(permCol, compressedPermCol);
}

int readFile(ifstream &inFileStream, char *outFileName, int numRow) {
	ofstream outFile(outFileName);
	if (!outFile.is_open())
	{
		cout << "Cannot open file " << outFileName << '\n';
		return 0;
	}
	string line;
	int cnt = 0;
	while (getline(inFileStream, line))
	{
		if(line == "") {
			continue;
		}
		outFile << line << '\n';
		++cnt;
		if (cnt >= numRow) {
			break;
		}
	}
	outFile.close();
	return cnt;
}


int readVCFFile(IQTree *tree, Alignment **alignment, Params &params) {
	char *aln_file = params.aln_file;
	ifstream in;
    in.exceptions(ios::failbit | ios::badbit);
    in.open(aln_file);
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
	double sumTimeReadPartialVCF = 0, sumTimeInitialize = 0;

	while (true) {
		startTime = getCPUTime();
		int cnt = (*alignment)->readPartialVCF(in, params.sequence_type, savePermCol, params.numStartRow, totalColumn, 8); 
		if(cnt == 0) {
			break;
		}
		sumTimeReadPartialVCF += (double)(getCPUTime() - startTime);

		tree->clearAllPartialLH();
		totalColumn += cnt;
		startTime = getCPUTime();
		initialize(tree, *alignment, savePermCol, permCol, compressedPermCol);
		sumTimeInitialize += (double)(getCPUTime() - startTime);
	}

	cout << "TimeReadPartialVCF: " << fixed << setprecision(3) << sumTimeReadPartialVCF << endl;
	cout << "TimeInitialize: " << fixed << setprecision(3) << sumTimeInitialize << endl;

	in.close();
	return totalColumn;
}

void addMoreRowMutation(Params &params)
{
	Alignment *alignment;

	IQTree *tree;
	tree = new IQTree;

	char *file_name = params.mutation_tree_file;
	bool is_rooted = false;

	if(params.tree_zip_file != NULL) {
		tree->readTree(params.tree_zip_file, file_name, is_rooted);
	} else {
		tree->readTree(file_name, is_rooted);
	}

	int sz = readVCFFile(tree, &alignment, params);
	// Init new tree's memory
	tree->cur_missing_sample_mutations.resize(sz + 1);
	tree->cur_ancestral_mutations.resize(sz + 1);
	tree->visited_missing_sample_mutations.resize(sz + 1);
	tree->visited_ancestral_mutations.resize(sz + 1);
	tree->cur_excess_mutations.resize(sz + 1);
	tree->visited_excess_mutations.resize(sz + 1);

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

	for (int i = 0; i < numSample; ++i)
	{
		vector<pair<PhyloNode *, PhyloNeighbor *>> bfs = tree->breadth_first_expansion();

		int total_nodes = (int)bfs.size();
		CandidateNode inp;
		int best_set_difference = INF;
		size_t best_node_num_leaves = INF;
		size_t best_distance = INF;
		std::vector<Mutation> excess_mutations;
		std::vector<bool> node_has_unique(total_nodes, false);
		bool best_node_has_unique = false;
		size_t best_j = 0;

		inp.best_set_difference = &best_set_difference;
		inp.best_node_num_leaves = &best_node_num_leaves;
		inp.best_distance = &best_distance;
		inp.node = (PhyloNode*)tree->root->neighbors[0]->node;
		inp.node_branch = (PhyloNeighbor*)inp.node->findNeighbor(tree->root);
		inp.missing_sample_mutations = &missingSamples[i].mutations;
		inp.excess_mutations = &excess_mutations;
		inp.has_unique = &best_node_has_unique;
		inp.node_has_unique = &(node_has_unique);
		inp.best_j = &best_j;

		tree->initDataCalculatePlacementMutation(inp);
		tree->optimizedCalculatePlacementMutation(inp, 0, true);

		for(int j = 0; j < total_nodes; ++j) {
			if(inp.best_node == bfs[j].first) {
				best_j = j;
				break;
			}
		}
		*inp.best_set_difference = INF;
		inp.j = best_j;
		inp.node = bfs[best_j].first;
		inp.node_branch = bfs[best_j].second;
		tree->calculatePlacementMutation(inp, false, true);
		tree->addNewSample(bfs[best_j].first, bfs[best_j].second, excess_mutations, i, missingSamples[i].name);
	}
	cout << "New tree's parsimony score: " << tree->computeParsimonyScoreMutation() << '\n';
	cout << "Time: " << fixed << setprecision(3) << (double)(getCPUTime() - startTime) << " seconds\n";
	cout << "Memory: " << getMemory() << " KB\n";
	
	// free memory
	tree->cur_missing_sample_mutations.clear();
	tree->cur_ancestral_mutations.clear();
	tree->visited_missing_sample_mutations.clear();
	tree->visited_ancestral_mutations.clear();
	tree->cur_excess_mutations.clear();
	tree->visited_excess_mutations.clear();

	// ofstream fout("addedTree.txt");
	// tree->printTree(fout, WT_SORT_TAXA | WT_NEWLINE);

	// stringstream ss;
	// tree->printTree(ss, WT_SORT_TAXA | WT_NEWLINE);
	// string treeAfterPhase1 = ss.str();

	// alignment->addToAlignmentNewSeq(alignment->remainName, alignment->remainSeq, savePermCol);
	// // tree->checkMutation(pos);
	// params.numStartRow = alignment->getNSeq();
	// // params.numAddRow = 0;

	// cout << "========= Starting optimization =========\n";
	// const int start_time = getCPUTime();
	// if (params.ppoptspr) {
	// 	treeAfterPhase1 = ppRunOriginalSpr(alignment, params, treeAfterPhase1);
	// }

	// if (params.ppopttbr) {
	// 	treeAfterPhase1 = ppRunOriginalTbr(alignment, params, treeAfterPhase1);
	// }
	// const int end_time = getCPUTime();
	// cout << "Time: " << fixed << setprecision(3) << (double)(end_time - start_time) << " seconds\n";
	
	
	delete alignment;
	alignment = NULL;
	delete tree;
}