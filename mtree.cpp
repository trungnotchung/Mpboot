/***************************************************************************
 *   Copyright (C) 2006 by BUI Quang Minh, Steffen Klaere, Arndt von Haeseler   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "mtree.h"
#include <iostream>
//#include <fstream>
#include <iterator>
#include "splitgraph.h"
using namespace std;

/*********************************************
	class MTree
*********************************************/

MTree::MTree() {
    root = NULL;
    leafNum = 0;
    nodeNum = 0;
    rooted = false;
    num_precision = 6;
    len_scale = 1.0;
	fig_char = "|-+++";
}

MTree::MTree(const char *userTreeFile, bool &is_rooted)
{
    init(userTreeFile, is_rooted);
}

void MTree::init(const char *userTreeFile, bool &is_rooted) {
    num_precision = 10;
    len_scale = 1.0;
    readTree(userTreeFile, is_rooted);
    //printInfo();
	fig_char = "|-+++";
}


/**
	constructor, get from another tree
*/
MTree::MTree(MTree &tree) {
    init(tree);
}

void MTree::init(MTree &tree) {
    root = tree.root;
    leafNum = tree.leafNum;
    nodeNum = tree.nodeNum;
    rooted = tree.rooted;
    //userFile = tree.userFile;
    // have to delete the root when exchange to another object
    tree.root = NULL;
    num_precision = tree.num_precision;
    len_scale = tree.len_scale;
    fig_char = tree.fig_char;
}

void MTree::copyTree(MTree *tree) {
    if (root) freeNode();
    stringstream ss;
    tree->printTree(ss);
    readTree(ss, tree->rooted);
}

void MTree::copyTree(MTree *tree, string &taxa_set) {
    if (tree->leafNum != taxa_set.length()) outError("#leaves and taxa_set do not match!");
    leafNum = nodeNum = branchNum = 0;
    for (string::iterator it = taxa_set.begin(); it != taxa_set.end(); it++)
        nodeNum += (*it);
    double new_len;
    if (root) freeNode();
    root = NULL;
    root = copyTree(tree, taxa_set, new_len);
}

Node* MTree::copyTree(MTree *tree, string &taxa_set, double &len, Node *node, Node *dad) {
    if (!node) {
        if (taxa_set[tree->root->id]) {
            node = tree->root;
        } else {
            for (int i = 0; i < tree->leafNum; i++)
                if (taxa_set[i]) {
                    node = tree->findNodeID(i);
                    break;
                }
        }
    }
    Node *new_node = NULL;
    NodeVector new_nodes;
    DoubleVector new_lens;
    if (node->isLeaf()) {
        len = 0.0;
        if (taxa_set[node->id]) {
            new_node = newNode(leafNum++, node->name.c_str());
        }
        if (dad) return new_node;
    }
    if (new_node) {
        new_nodes.push_back(new_node);
        new_lens.push_back(len);
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
        double new_len;
        new_node = copyTree(tree, taxa_set, new_len, (*it)->node, node);
        if (new_node) {
            new_nodes.push_back(new_node);
            new_lens.push_back((*it)->length + new_len);
        }
    }
    if (new_nodes.empty()) return NULL;
    if (new_nodes.size() == 1) {
        len = new_lens[0];
        return new_nodes[0];
    }
    if (!dad && new_nodes.size() == 2) {
        double sum_len = new_lens[0] + new_lens[1];
        new_nodes[0]->addNeighbor(new_nodes[1], sum_len, branchNum);
        new_nodes[1]->addNeighbor(new_nodes[0], sum_len, branchNum);
        branchNum++;
        return new_nodes[0];
    }
    Node* int_node = newNode(nodeNum++, node->name.c_str());
    len = 0.0;
    for (int i = 0; i < new_nodes.size(); i++) {
        int_node->addNeighbor(new_nodes[i], new_lens[i], branchNum);
        new_nodes[i]->addNeighbor(int_node, new_lens[i], branchNum);
        branchNum++;
    }
    return int_node;
}

Node* MTree::newNode(int node_id, const char* node_name) {
    return new Node(node_id, node_name);
}

Node* MTree::newNode(int node_id, int node_name) {
    return new Node(node_id, node_name);
}

bool MTree::isBifurcating(Node *node, Node *dad) {
	if (!node) node = root;
	if (!node->isLeaf() && node->degree() != 3) return false;
	FOR_NEIGHBOR_IT(node, dad, it) {
		if (!(*it)->node->isLeaf() && (*it)->node->degree() != 3) return false;
		if (!isBifurcating((*it)->node, node)) return false;
	}
	return true;
}

void MTree::printBranchLengths(ostream &out, Node *node, Node *dad)
{
    if (node == NULL) {
    	node = root;
    	sortTaxa();
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (node->name != "") out << node->name; else out << node->id;
        out << "\t";
        if ((*it)->node->name != "") out << (*it)->node->name; else out << (*it)->node->id;
        out << "\t" << (*it)->length << endl;
        printBranchLengths(out, (*it)->node, node);
    }
}

int MTree::countZeroBranches(Node *node, Node *dad, double epsilon) {
    int count = 0;
    if (node == NULL) node = root;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if ((*it)->length <= epsilon) count++;
        count += countZeroBranches((*it)->node, node, epsilon);
    }
    return count;
}

int MTree::countZeroInternalBranches(Node *node, Node *dad, double epsilon) {
    int count = 0;
    if (node == NULL) node = root;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if ((*it)->length <= epsilon && !(*it)->node->isLeaf() && !node->isLeaf()) count++;
        count += countZeroInternalBranches((*it)->node, node, epsilon);
    }
    return count;

}

int MTree::countLongBranches(Node *node, Node *dad, double upper_limit) {
    int count = 0;
    if (node == NULL) node = root;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if ((*it)->length >= upper_limit) count++;
        count += countLongBranches((*it)->node, node, upper_limit);
    }
    return count;
}


void MTree::printTree(const char *ofile, int brtype)
{
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        if (brtype & WT_APPEND)
            out.open(ofile, ios_base::out | ios_base::app);
        else
            out.open(ofile);
        printTree(out, brtype);
        out.close();
        if (verbose_mode >= VB_DEBUG)
            cout << "Tree was printed to " << ofile << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, ofile);
    }
}

string MTree::getTreeString() {
	stringstream tree_stream;
	printTree(tree_stream);
	return tree_stream.str();
}

void MTree::printTree(ostream &out, int brtype) {
    if (root->isLeaf()) {
        if (root->neighbors[0]->node->isLeaf()) {
            // tree has only 2 taxa!
            out << "(";
            printTree(out, brtype, root);
            out << ",";
            if (brtype & WT_TAXON_ID)
                out << root->neighbors[0]->node->id;
            else
                out << root->neighbors[0]->node->name;

            if (brtype & WT_BR_LEN)
                out << ":0";
            out << ")";
        } else
            // tree has more than 2 taxa
            printTree(out, brtype, root->neighbors[0]->node);
    } else
        printTree(out, brtype, root);

    out << ";";
    if (brtype & WT_NEWLINE) out << endl;
}

struct IntString {
    int id;
    string str;
};

/**
	nodecmp, for pruning algorithm
*/
struct IntStringCmp
{
    /**
    	nodecmp, for pruning algorithm
    */
    bool operator()(const IntString* s1, const IntString* s2) const
    {
        return (s1->id) < (s2->id);
    }
};

typedef set<IntString*, IntStringCmp> IntStringSet;

int MTree::printTree(ostream &out, int brtype, Node *node, Node *dad)
{
    int smallest_taxid = leafNum;
    out.precision(num_precision);
    if (!node) node = root;
    if (node->isLeaf()) {
        smallest_taxid = node->id;
        if (brtype & WT_TAXON_ID)
            out << node->id;
        else
            out << node->name;

        if (brtype & WT_BR_LEN) {
        	out.setf( std::ios::fixed, std:: ios::floatfield ); // some sofware does handle number format like '1.234e-6'
            out.precision(15); // increase precision to avoid zero branch (like in RAxML)
        	double len = node->neighbors[0]->length;
            if (brtype & WT_BR_SCALE) len *= len_scale;
            if (brtype & WT_BR_LEN_ROUNDING) len = round(len);
            if (brtype & WT_BR_LEN_FIXED_WIDTH)
                out << ":" << fixed << len;
            else
                out << ":" << len;
        }
    } else {
        // internal node
        out << "(";
        bool first = true;
        double length = 0.0;
        //for (int i = 0; i < node->neighbors.size(); i++)
        //if (node->neighbors[i]->node != dad)
        if (! (brtype & WT_SORT_TAXA)) {
            FOR_NEIGHBOR_IT(node, dad, it) {
                if ((*it)->node->name != ROOT_NAME) {
                    if (!first)
                        out << ",";
                    int taxid = printTree(out, brtype, (*it)->node, node);
                    if (taxid < smallest_taxid) smallest_taxid = taxid;
                    first = false;
                } else
                    length = (*it)->length;
            } else {
                length = (*it)->length;
            }
        } else {
            IntStringSet strout;
            FOR_NEIGHBOR_IT(node, dad, it) {
                if ((*it)->node->name != ROOT_NAME) {
                    ostringstream ss;
                    IntString *str = new IntString;
                    str->id = printTree(ss, brtype, (*it)->node, node);
                    //ss.flush();
                    str->str = ss.str();
                    strout.insert(str);
                } else
                    length = (*it)->length;
            } else {
                length = (*it)->length;
            }
            smallest_taxid = (*strout.begin())->id;
            IntStringSet::iterator iss;
            for (iss = strout.begin(); iss != strout.end(); iss++) {
                if (!first) out << ",";
                out << (*iss)->str;
                first = false;
            }
            for (iss = strout.begin(); iss != strout.end(); iss++)
                delete (*iss);
        }
        out << ")";
        if (!node->name.empty())
            out << node->name;
        else if (brtype & WT_INT_NODE)
            out << node->id;
        if (dad != NULL || length > 0.0) {
            if (brtype & WT_BR_SCALE) length *= len_scale;
            if (brtype & WT_BR_LEN_ROUNDING) length = round(length);
            if (brtype & WT_BR_LEN) {
                if (brtype & WT_BR_LEN_FIXED_WIDTH)
                    out << ":" << fixed << length;
                else
                    out << ":" << length;
            } else if (brtype & WT_BR_CLADE) {
                if (! node->name.empty()) out << "/";
                out << length;
            }
        }
    }
    return smallest_taxid;
}


void MTree::printSubTree(ostream &out, NodeVector &subtree) {
    if (root->isLeaf())
        printSubTree(out, subtree, root->neighbors[0]->node);
    else
        printSubTree(out, subtree, root);
    out << ";";
}

void MTree::printSubTree(ostream &out, NodeVector &subtree, Node *node, Node *dad) {
    if (!node) node = root;

    NeighborVec::iterator it;
    double length = 0.0, dad_length = 0.0;
    // go down if only 1 child available
    Node *child = NULL;
    int degree;
    do {
        degree = 0;
        FOR_NEIGHBOR(node, dad, it) {
            if (subtree[(*it)->node->id] != NULL) {
                degree++;
                child = (*it)->node;
            }
        } else dad_length = (*it)->length;

        if (degree == 1) {
            dad = node;
            node = child;
            length += dad_length;
        }
    } while (degree == 1 && !node->isLeaf());

    if (node->isLeaf())
        out << node->name << ":" << node->neighbors[0]->length + length;
    else
    {
        // internal node
        out << "(";
        bool first = true;

        FOR_NEIGHBOR(node, dad, it)	{
            if (subtree[(*it)->node->id] != NULL) {
                if ((*it)->node->name != ROOT_NAME) {
                    if (!first)
                        out << ",";
                    printSubTree(out, subtree, (*it)->node, node);
                    first = false;
                } else
                    length += (*it)->length;
            }
        } else {
            length += (*it)->length;
        }
        out << ")";
        if (!node->name.empty())
            out << node->name;
        if (dad != NULL || length > 1e-20)
            out << ":" << length;
    }
}


void MTree::printTaxa(const char *ofile)
{
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(ofile);
        if (root->isLeaf())
            printTaxa(out, root->neighbors[0]->node);
        else
            printTaxa(out);
        out.close();
        cout << "Taxa list was printed to " << ofile << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, ofile);
    }
}

void MTree::printTaxa(ostream &out, Node *node, Node *dad)
{
    if (!node) node = root;
    if (node->isLeaf())
        out << node->name << endl;
    else
    {
        // internal node
        //for (int i = 0; i < node->neighbors.size(); i++)
        //if (node->neighbors[i]->node != dad)
        FOR_NEIGHBOR_IT(node, dad, it)	{
            printTaxa(out, (*it)->node, node);
        }
    }
}

void MTree::printTaxa(ostream &out, NodeVector &subtree) {
    for (int i = 0; i < leafNum; i++)
        if (subtree[i] != NULL) {
            out << subtree[i]->name << endl;
        }
}

void MTree::readTree(const char *zipFileName, const char *infile, bool &is_rooted) {
    int err;
    zip *archive = zip_open(zipFileName, 0, &err);
    if (!archive) {
        std::cerr << "Can't open zip file: " << zip_strerror(archive) << std::endl;
        return;
    }

    int fileIndex = zip_name_locate(archive, infile, 0);
    if (fileIndex < 0) {
        std::cerr << "Can't read file in zip: " << zip_strerror(archive) << std::endl;
        zip_close(archive);
        return;
    }

    struct zip_stat stat;
    zip_stat_init(&stat);
    zip_stat_index(archive, fileIndex, 0, &stat);

    zip_file *file = zip_fopen_index(archive, fileIndex, 0);

    if (file) {
        std::ofstream outputStream("treeFile.txt", std::ios::binary);
        if (outputStream.is_open()) {
            char buffer[1024];
            zip_int64_t bytesRead;
            while ((bytesRead = zip_fread(file, buffer, sizeof(buffer))) > 0) {
                outputStream.write(buffer, bytesRead);
            }

            outputStream.close();
        } else {
            std::cerr << "Can't open output file" << std::endl;
        }
    } else {
        std::cerr << "Không thể mở file trong zip: " << zip_strerror(archive) << std::endl;
    }

    zip_close(archive);

    readTree("treeFile.txt", is_rooted);
    std::remove("treeFile.txt");
}

void MTree::readTree(const char *infile, bool &is_rooted) {
    ifstream in;
    try {
        in.exceptions(ios::failbit | ios::badbit);
        in.open(infile);
        readTree(in, is_rooted);
        in.close();
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, infile);
    }

    rooted = is_rooted;

    if (verbose_mode >= VB_MED)
        cout << "Tree contains " << leafNum - is_rooted <<
             " taxa and " << nodeNum-1-is_rooted << " branches" << endl;
}


void MTree::readTree(istream &in, bool &is_rooted)
{
    in_line = 1;
    in_column = 1;
    try {
        char ch;
        ch = readNextChar(in);
        if (ch != '(') {
        	cout << in.rdbuf() << endl;
            throw "Tree file not started with an opening-bracket '('";
        }

        leafNum = 0;

        double branch_len;
        Node *node;
        parseFile(in, ch, node, branch_len);
        if (is_rooted || branch_len > 0.0) {
            if (branch_len == -1.0) branch_len = 0.0;
            if (branch_len < 0.0)
                throw ERR_NEG_BRANCH;
            is_rooted = true;
            root = newNode(leafNum, ROOT_NAME);
            root->addNeighbor(node, branch_len);
            node->addNeighbor(root, branch_len);
            leafNum++;
        } else { // assign root to one of the neighbor of node, if any
            FOR_NEIGHBOR_IT(node, NULL, it)
            if ((*it)->node->isLeaf()) {
                root = (*it)->node;
                break;
            }
        }
        // make sure that root is a leaf
        assert(root->isLeaf());

        if (in.eof() || ch != ';')
            throw "Tree file must be ended with a semi-colon ';'";
    } catch (bad_alloc) {
        outError(ERR_NO_MEMORY);
    } catch (const char *str) {
        outError(str, reportInputInfo());
    } catch (string str) {
        outError(str.c_str(), reportInputInfo());
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, reportInputInfo());
    } catch (...) {
        // anything else
        outError(ERR_READ_ANY, reportInputInfo());
    }

    nodeNum = leafNum;
    initializeTree();

    //bool stop = false;
    //checkValidTree(stop);
}

void MTree::initializeTree(Node *node, Node* dad)
{
    if (!node) {
        node = root;
        nodeNum = leafNum;
        branchNum = 0;
    }
    if (!node->isLeaf())
    {
        node->id = nodeNum;
        nodeNum++;
        //node->name = node->id;

    }
    //for (int i = 0; i < node->neighbors.size(); i++)
    //if (node->neighbors[i]->node != dad)
    FOR_NEIGHBOR_IT(node, dad, it) {
        (*it)->id = branchNum;
        (*it)->node->findNeighbor(node)->id = branchNum;
        branchNum++;
        initializeTree((*it)->node, node);
    }
}


void MTree::parseFile(istream &infile, char &ch, Node* &root, double &branch_len)
{
    Node *node;
    int maxlen = 10000;
    char seqname[10000];
    int seqlen;
    double brlen;
    branch_len = -1.0;

    root = newNode();

    if (ch == '(') {
        // internal node
        ch = readNextChar(infile);
        while (ch != ')' && !infile.eof())
        {
            node = NULL;
            parseFile(infile, ch, node, brlen);
            //if (brlen == -1.0)
            //throw "Found branch with no length.";
            //if (brlen < 0.0)
            //throw ERR_NEG_BRANCH;
            root->addNeighbor(node, brlen);
            node->addNeighbor(root, brlen);
            if (infile.eof())
                throw "Expecting ')', but end of file instead";
            if (ch == ',')
                ch = readNextChar(infile);
            else if (ch != ')') {
                string err = "Expecting ')', but found '";
                err += ch;
                err += "' instead";
                throw err;
            }
        }
        if (!infile.eof()) ch = readNextChar(infile);
    }
    // now read the node name
    seqlen = 0;
    char end_ch = 0;
    if (ch == '\'' || ch == '"') end_ch = ch;

    while (!infile.eof() && seqlen < maxlen)
    {
        if (end_ch == 0) {
            if (is_newick_token(ch) || controlchar(ch)) break;
        }
        seqname[seqlen++] = ch;
        ch = infile.get();
        in_column++;
        if (end_ch != 0 && ch == end_ch) {
            seqname[seqlen++] = ch;
            break;
        }
    }
    if ((controlchar(ch) || ch == '[' || ch == end_ch) && !infile.eof())
        ch = readNextChar(infile, ch);
    if (seqlen == maxlen)
        throw "Too long name ( > 100)";
    seqname[seqlen] = 0;
    if (seqlen == 0 && root->isLeaf())
        throw "A taxon has no name.";
    if (seqlen > 0)
        root->name.append(seqname);
    if (root->isLeaf()) {
        // is a leaf, assign its ID
        root->id = leafNum;
        if (leafNum == 0)
            MTree::root = root;
        leafNum++;
    }

    if (ch == ';' || infile.eof())
        return;
    if (ch == ':')
    {
        ch = readNextChar(infile);
        seqlen = 0;
        while (!is_newick_token(ch) && !controlchar(ch) && !infile.eof() && seqlen < maxlen)
        {
            seqname[seqlen] = ch;
            seqlen++;
            ch = infile.get();
            in_column++;
        }
        if ((controlchar(ch) || ch == '[') && !infile.eof())
            ch = readNextChar(infile, ch);
        if (seqlen == maxlen || infile.eof())
            throw "branch length format error.";
        seqname[seqlen] = 0;
        branch_len = convert_double(seqname);
    }
}

/**
	check tree is bifurcating tree (every leaf with level 1 or 3)
*/
void MTree::checkValidTree(bool &stop, Node *node, Node *dad)
{
    if (!node) node = root;
    if (node->degree() != 1 && node->degree() != 3) {
        cout << "Tree is not bifurcating." << endl;
        stop = true;
        return;
    }
    //for (int i = 0; i < node->neighbors.size(); i++)
    //if (node->neighbors[i]->node != dad) {
    FOR_NEIGHBOR_IT(node, dad, it) {
        checkValidTree(stop, (*it)->node, node);
        if (stop)
            return;
    }
}

double MTree::treeLength(Node *node, Node *dad)
{
    if (!node) node = root;
    double sum = 0;
    FOR_NEIGHBOR_IT(node, dad, it) {
        sum += (*it)->length + treeLength((*it)->node, node);
    }
    return sum;
}

double MTree::treeLengthInternal( double epsilon, Node *node, Node *dad)
{
    if (!node) node = root;
    double sum = 0;
    FOR_NEIGHBOR_IT(node, dad, it) {
    	if (!(*it)->node->isLeaf() && !node->isLeaf())
    	{
    		if (treeLength((*it)->node, node) > epsilon) {
    			sum += (*it)->length + treeLengthInternal(epsilon, (*it)->node, node);
    		}
    	}
    	else {
    		if (treeLength((*it)->node, node) > epsilon) {
    			sum += treeLengthInternal(epsilon, (*it)->node, node);
    		}
    	}
    }
    return sum;
}

double MTree::treeDepth(Node *node, Node *dad)
{
    if (!node) node = root;
    double maxsum = 0.0;
    FOR_NEIGHBOR_IT(node, dad, it) {
        double len = (*it)->length;
        if (len < 0.0) len = 0.0;
        double sum = len + treeDepth((*it)->node, node);
        if (sum > maxsum) maxsum = sum;
    }
    return maxsum;
}

void MTree::getNonCherryLeaves(NodeVector &noncherry, NodeVector &cherry, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
    	if (node->isInCherry()) {
    		cherry.push_back(node);
    	} else {
            noncherry.push_back(node);
    	}
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
    	getNonCherryLeaves(noncherry, cherry, (*it)->node, node);
    }
}

void MTree::getTaxa(NodeVector &taxa, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
        taxa.push_back(node);
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
        getTaxa(taxa, (*it)->node, node);
    }
}

void MTree::getAllNodesInSubtree(Node *node, Node *dad, NodeVector &nodeList) {
    assert(node && dad);
    nodeList.push_back(node);
    if (node->isLeaf()) {
        return;
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
        getAllNodesInSubtree((*it)->node, node, nodeList);
    }
}

int MTree::getNumTaxa(Node *node, Node *dad) {
	// This function is WRONG because it will always return 1 if calling from the ROOT
    if (!node) node = root;
    if (node->isLeaf()) {
        return 1;
    }
    int numLeaf = 0;
    FOR_NEIGHBOR_IT(node, dad, it) {
        numLeaf += getNumTaxa((*it)->node, node);
    }
    return numLeaf;
}

void MTree::getInternalNodes(NodeVector &nodes, Node *node, Node *dad) {
    if (!node) node = root;
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it)
    if (!(*it)->node->isLeaf()) {
        getInternalNodes(nodes, (*it)->node, node);
        nodes.push_back((*it)->node);
    }
}

void MTree::getInternalBranches(NodeVector &nodes, NodeVector &nodes2, Node *node, Node *dad) {
    if (!node) node = root;
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it)
    if (!(*it)->node->isLeaf()) {
        getInternalBranches(nodes, nodes2, (*it)->node, node);
        if (!node->isLeaf()) {
            if (node->id < (*it)->node->id) {
                nodes.push_back(node);
                nodes2.push_back((*it)->node);
            } else {
                nodes.push_back((*it)->node);
                nodes2.push_back(node);
            }
        }
    }
}

void MTree::getInBranches(map<string, Branch> &brans, int depth, Node *node, Node *dad) {
    if (depth == 0)
      return;
    assert(isInBran(node, dad));
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (!(*it)->node->isLeaf()) {
            Branch bran(node, (*it)->node);
            brans.insert(pair<string, Branch>(bran.getKey(), bran));
            getInBranches(brans, depth-1, (*it)->node, node);
        }
    }
}

bool MTree::isInBran(Node* node1, Node* node2) {
    return (!node1->isLeaf() && !node2->isLeaf());
}



void MTree::getBranches(NodeVector &nodes, NodeVector &nodes2, Node *node, Node *dad) {
    if (!node) node = root;
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (node->id < (*it)->node->id) {
            nodes.push_back(node);
            nodes2.push_back((*it)->node);
        } else {
            nodes.push_back((*it)->node);
            nodes2.push_back(node);
        }
        getBranches(nodes, nodes2, (*it)->node, node);
    }
}

void MTree::getBranchLengths(DoubleVector &len, Node *node, Node *dad) {
    if (!node) {
        node = root;
        assert(len.size() == branchNum);
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
        len[(*it)->id] = (*it)->length;
        getBranchLengths(len, (*it)->node, node);
    }
}

void MTree::setBranchLengths(DoubleVector &len, Node *node, Node *dad) {
    if (!node) {
        node = root;
        assert(len.size() == branchNum);
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
        (*it)->length = (*it)->node->findNeighbor(node)->length = len[(*it)->id];
        setBranchLengths(len, (*it)->node, node);
    }
}

void MTree::getOrderedTaxa(NodeVector &taxa, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
        if (taxa.empty()) taxa.resize(leafNum);
        taxa[node->id] = node;
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
        getOrderedTaxa(taxa, (*it)->node, node);
    }
}

void MTree::getTaxaName(vector<string> &taxname, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
        if (taxname.empty()) taxname.resize(leafNum);
        taxname[node->id] = node->name;
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
        getTaxaName(taxname, (*it)->node, node);
    }
}


void MTree::getTaxaID(vector<int> &taxa, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
        taxa.push_back(node->id);
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
        getTaxaID(taxa, (*it)->node, node);
    }
}


void MTree::convertSplits(SplitGraph &sg, Split *resp, NodeVector *nodes, Node *node, Node *dad) {
    if (!node) node = root;
    assert(resp->getNTaxa() == leafNum);
    bool has_child = false;
    FOR_NEIGHBOR_IT(node, dad, it) {
        //vector<int> taxa;
        //getTaxaID((*it)->node, node, taxa);

        Split *sp = new Split(leafNum, (*it)->length);
        convertSplits(sg, sp, nodes, (*it)->node, node);
        *resp += *sp;
        if (sp->shouldInvert())
            sp->invert();
		 /* ignore nodes with degree of 2 because such split will be added before */
        if (node->degree() != 2) {
		  sg.push_back(sp);
          if (nodes) nodes->push_back((*it)->node);
        }
        has_child = true;
    }
    if (!has_child)
        resp->addTaxon(node->id);
}

void MTree::convertSplits(vector<string> &taxname, SplitGraph &sg, NodeVector *nodes, Node *node, Node *dad) {
    if (!sg.taxa) {
        sg.taxa = new NxsTaxaBlock();
        for (vector<string>::iterator it = taxname.begin(); it != taxname.end(); it++)
            sg.taxa->AddTaxonLabel(NxsString(it->c_str()));
    }
    if (!sg.splits)
        sg.splits = new MSplitsBlock(&sg);
    if (!sg.pda)
        sg.pda = new MPdaBlock(&sg);

    // make the cycle
    getTaxaID(sg.splits->cycle);
    // make the splits
    Split sp(leafNum);
    convertSplits(sg, &sp, nodes, node, dad);
}

void MTree::convertSplits(SplitGraph &sg, NodeVector *nodes, Node *node, Node *dad) {

    // make the taxa name
    vector<string> taxname;
    taxname.resize(leafNum);
    getTaxaName(taxname);

    convertSplits(taxname, sg, nodes, node, dad);
}

inline int splitnumtaxacmp(const Split* a, const Split* b)
{
    return (a->countTaxa() < b->countTaxa());
}

void MTree::convertToTree(SplitGraph &sg) {
    SplitGraph::iterator it;
    int taxid;
    int count;
	BoolVector has_tax;
	has_tax.resize(sg.getNTaxa(), false);
	// first add trivial splits if not existed
	for (it = sg.begin(); it != sg.end(); it++) {
		taxid = (*it)->trivial();
		if (taxid >= 0) has_tax[taxid] = true;
	}
	for (count = 0; count < has_tax.size(); count++)
		if (!has_tax[count]) {
			Split *sp = new Split(sg.getNTaxa());
			sp->addTaxon(count);
			sg.push_back(sp);
		}
    // sort splits by the number of taxa they contain
    sort(sg.begin(), sg.end(), splitnumtaxacmp);

    // initialize the tree
    rooted = false;
    leafNum = sg.getNTaxa();
    nodeNum = leafNum;

    // create the ground nodes, first as the leaves
    NodeVector leaves;
    vector<Split*> cladetaxa;
    leaves.resize(leafNum, NULL);
    cladetaxa.resize(leafNum, NULL);
    // first add all trivial splits into tree
    for (it = sg.begin(), count = 0; it != sg.end(); it++, count++) {
        //(*it)->report(cout);
        taxid = (*it)->trivial();
        if (taxid < 0) break;
        assert(leaves[taxid] == NULL);
        leaves[taxid] = newNode(taxid, sg.getTaxa()->GetTaxonLabel(taxid).c_str());
        leaves[taxid]->addNeighbor(NULL, (*it)->getWeight());
        cladetaxa[taxid] = (*it);
    }
    // now fill in all missing taxa with zero terminal branch
    for (taxid = 0; taxid < leafNum; taxid++)
        assert(leaves[taxid]);

    // now add non-trivial splits, cotinue with the interrupted iterator
    for (/*it = sg.begin()*/; it != sg.end(); it++) {
        //(*it)->report(cout);
        Split *mysp = *it;
        Node *newnode = newNode(nodeNum);
        int count = 0;

        for (taxid = 0; taxid < leaves.size(); )
            if (cladetaxa[taxid]->subsetOf(*mysp)) // clade is a subset of current split
            {
                count += cladetaxa[taxid]->countTaxa();
                double len = leaves[taxid]->updateNeighbor(NULL, newnode);
                newnode->addNeighbor(leaves[taxid], len);
                leaves[taxid] = leaves.back();
                leaves.pop_back();
                cladetaxa[taxid] = cladetaxa.back();
                cladetaxa.pop_back();
            } else taxid++;
        assert(count == mysp->countTaxa());
        cladetaxa.push_back(mysp);
        leaves.push_back(newnode);

        newnode->addNeighbor(NULL, mysp->getWeight());
        nodeNum++;
    }
    assert(leaves.size() >= 3);
    Node *newnode = newNode(nodeNum);
    for (taxid = 0; taxid < leaves.size(); taxid++) {
        double len = leaves[taxid]->updateNeighbor(NULL, newnode);
        newnode->addNeighbor(leaves[taxid], len);
    }
    root = newnode;
    nodeNum++;
    cladetaxa.clear();
}

Node *MTree::findNodeName(string &name, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->name == name) return node;
    FOR_NEIGHBOR_IT(node, dad, it) {
        Node *res = findNodeName(name, (*it)->node, node);
        if (res) return res;
    }
    return NULL;
}

Node *MTree::findLeafName(string &name, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf() && node->name == name) return node;
    FOR_NEIGHBOR_IT(node, dad, it) {
        Node *res = findLeafName(name, (*it)->node, node);
        if (res) return res;
    }
    return NULL;
}

Node *MTree::findNodeID(int id, Node *node, Node* dad) {
    if (!node) node = root;
    if (node->id == id) return node;
    FOR_NEIGHBOR_IT(node, dad, it) {
        Node *res = findNodeID(id, (*it)->node, node);
        if (res) return res;
    }
    return NULL;
}


void MTree::scaleLength(double norm, bool make_int, Node *node, Node *dad) {
    if (!node) node = root;
    FOR_NEIGHBOR_DECLARE(node, NULL, it) {
        (*it)->length *= norm;
        if (make_int)
            (*it)->length = round((*it)->length);
    }

    FOR_NEIGHBOR(node, dad, it) {
        scaleLength(norm, make_int, (*it)->node, node);
    }
}

void MTree::transformBranchLenRAX(double factor, Node *node, Node *dad) {
    if (!node) node = root;
    FOR_NEIGHBOR_DECLARE(node, NULL, it) {
        (*it)->length /= factor;
        (*it)->length = exp(-(*it)->length);
    }

    FOR_NEIGHBOR(node, dad, it) {
    	transformBranchLenRAX(factor, (*it)->node, node);
    }
}

void MTree::scaleCladeSupport(double norm, bool make_int, Node *node, Node *dad) {
    if (!node) node = root;
    if (!node->isLeaf() && !node->name.empty()) {
        double supp = 0.0;
        try {
            supp = convert_double(node->name.c_str());
        } catch (string str) {
            outError(str);
        }
        supp *= norm;
        if (make_int)
            supp = round(supp);
        node->name = "";
        node->name += supp;
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        scaleCladeSupport(norm, make_int, (*it)->node, node);
    }
}


MTree::~MTree()
{
    if (root != NULL)
        freeNode();
    root = NULL;
}

int MTree::freeNode(Node *node, Node *dad)
{
	if ( root == NULL )
		return 0;
    if (!node) node = root;
    NeighborVec::reverse_iterator it;
    int num_nodes = 1;
    for (it = node->neighbors.rbegin(); it != node->neighbors.rend(); it++)
        if ((*it)->node != dad) {
            num_nodes += freeNode((*it)->node, node);
        }
    delete node;
    return num_nodes;
}

char MTree::readNextChar(istream &in, char current_ch) {
    char ch;
    if (current_ch == '[')
        ch = current_ch;
    else {
        in.get(ch);
        in_column++;
        if (ch == 10) {
            in_line++;
            in_column = 1;
        }
    }
    while (controlchar(ch) && !in.eof()) {
        in.get(ch);
        in_column++;
        if (ch == 10) {
            in_line++;
            in_column = 1;
        }
    }
    // ignore comment
    while (ch=='[' && !in.eof()) {
        while (ch!=']' && !in.eof()) {
            in.get(ch);
            in_column++;
            if (ch == 10) {
                in_line++;
                in_column = 1;
            }
        }
        if (ch != ']') throw "Comments not ended with ]";
        in_column++;
        in.get(ch);
        if (ch == 10) {
            in_line++;
            in_column = 1;
        }
        while (controlchar(ch) && !in.eof()) {
            in_column++;
            in.get(ch);
            if (ch == 10) {
                in_line++;
                in_column = 1;
            }
        }
    }
    return ch;
}

string MTree::reportInputInfo() {
    string str = " (line ";
    str += convertIntToString(in_line) + " column " + convertIntToString(in_column-1) + ")";
    return str;
}


typedef map<int, Neighbor*> IntNeighborMap;

int MTree::sortTaxa(Node *node, Node *dad) {
    if (!node) {
        node = root;
        if (node->isLeaf()) node = node->neighbors[0]->node;
    }
    if (node->isLeaf())
        return node->id;
    IntNeighborMap taxid_nei_map;
    FOR_NEIGHBOR_IT(node, dad, it) {
        int taxid = sortTaxa((*it)->node, node);
        taxid_nei_map.insert(IntNeighborMap::value_type(taxid, (*it)));
    }
    ;
    int i = 0;
    for (IntNeighborMap::iterator it = taxid_nei_map.begin(); it != taxid_nei_map.end(); it++, i++) {
        if (node->neighbors[i]->node == dad) i++;
        node->neighbors[i] = it->second;
    }

    return taxid_nei_map.begin()->first;
}

void MTree::setExtendedFigChar() {
	//fig_char[0] = 179;
	//fig_char[1] = 196;
	fig_char[2] = '/';
	//fig_char[3] = 195;
	fig_char[4] = '\\';
}

void MTree::drawTree(ostream &out, int brtype, double zero_epsilon) {
    IntVector sub_tree_br;
    if (verbose_mode >= VB_DEBUG) {
        printTree(cout);
        cout << endl;
    }
    Node *node = root;
    if (node->isLeaf()) node = node->neighbors[0]->node;
    double scale = 60.0/treeDepth(node);
    //if (verbose_mode >= VB_DEBUG)
    //cout << "Tree depth: " << scale<< endl;
    drawTree2(out, brtype, scale, sub_tree_br, zero_epsilon);
    /*
    if (brtype & WT_INT_NODE)
        drawTree2(out, brtype, scale, sub_tree_br, zero_epsilon);
    else
        drawTree(out, brtype, scale, sub_tree_br, zero_epsilon);
    */
    out << endl;
}

/*
void MTree::drawTree(ostream &out, int brtype, double brscale, IntVector &subtree_br, double zero_epsilon, Node *node, Node *dad) {
    int i, br_len = 3;
    if (!node) {
        node = root;
        if (node->isLeaf()) node = node->neighbors[0]->node;
    } else {

        if (brtype & WT_BR_SCALE) {
            br_len = floor(node->findNeighbor(dad)->length * brscale)-1;
            if (br_len < 3) br_len = 3;
            //if (!node->isLeaf() && br_len < 4) br_len = 4;
        }
        out << '+';
        if ((brtype & WT_INT_NODE) && !node->isLeaf()) {
            string str = convertIntToString(node->id);
            for (i = 0; i < br_len-str.length(); i++) out << '-';
            out << node->id;
        } else
            for (i = 0; i < br_len; i++) out << '-';
    }
    if (node->isLeaf()) {
        out << node->name;
        if (brtype & WT_TAXON_ID)
            out << " (" << node->id << ")";
        out << endl;
        return;
    }
    int descendant_cnt = node->degree();
    if (dad) descendant_cnt--;
    int cnt = 0;
    subtree_br.push_back(br_len);
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (cnt == descendant_cnt-1)
            subtree_br.back() = -subtree_br.back();

        drawTree(out, brtype, brscale, subtree_br, zero_epsilon, (*it)->node, node);
        cnt++;
        if (cnt == descendant_cnt) break;
        for (IntVector::iterator it = subtree_br.begin()+1; it != subtree_br.end(); it++)
        {
            if ((*(it-1)) > 0) out << '|';
            else out << ' ';
            for (i = 0; i < abs(*it); i++) out << ' ';
        }
    }
    subtree_br.pop_back();
}
*/

void MTree::drawTree2(ostream &out, int brtype, double brscale, IntVector &subtree_br, double zero_epsilon, Node *node, Node *dad) {
    int i, br_len = 3;
    IntVector::iterator ii;
    bool zero_length = false;

//	char zero_char = '*'; // for ML
	char zero_char = '-'; // Diep: for MP

    //cout << "DrawTree2!" << endl;
    if (!node) {
        node = root;
        if (node->isLeaf()) node = node->neighbors[0]->node;
    } else {
        if (brtype & WT_BR_SCALE) {
            br_len = floor(node->findNeighbor(dad)->length * brscale)-1;
            if (br_len < 2) br_len = 2;
        }
        if (node->findNeighbor(dad)->length <= zero_epsilon) zero_length = true;
    }
    if (node->isLeaf()) {
        for (ii = subtree_br.begin()+1; ii != subtree_br.end(); ii++) {
            if (abs(*(ii-1)) > 1000) out << ' ';
            else out << fig_char[0];
            int num = abs(*ii);
            if (num > 1000) num -= 1000;
            for (i = 0; i < num; i++) out << ' ';
        }
        out << ((node==dad->neighbors.front()->node) ? fig_char[2] : ((node==dad->neighbors.back()->node) ? fig_char[4] : fig_char[3]));
        for (i = 0; i < br_len; i++)
            out << ((zero_length) ? zero_char : fig_char[1]);
        out << node->name;
        if (brtype & WT_TAXON_ID)
            out << " (" << node->id << ")";
        if (brtype & WT_BR_ID)
            out << " [" << node->neighbors[0]->id << "]";
        if (brtype & WT_BR_LEN)
            out << " " << node->neighbors[0]->length;
        //out << " ";
        //copy (subtree_br.begin(), subtree_br.end(), ostream_iterator<int> (out, " "));
        out << endl;
        return;
    }
    int descendant_cnt = node->degree();
    if (dad) descendant_cnt--;
    int cnt = 0;
    bool first = true;

    br_len = br_len+1000;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (cnt == descendant_cnt-1)
            br_len = -br_len;
        subtree_br.push_back(br_len);

        drawTree2(out, brtype, brscale, subtree_br, zero_epsilon, (*it)->node, node);
        subtree_br.pop_back();
        if (br_len > 1000) br_len -= 1000;
        cnt++;
        if (cnt == descendant_cnt) break;
        if (subtree_br.size() > 1)
            for (ii = subtree_br.begin()+1; ii != subtree_br.end(); ii++) {
                if (abs(*(ii-1)) > 1000) out << ' ';
                else out << fig_char[0];
                if (ii == subtree_br.begin()) continue;
                int num = abs(*ii);
                if (num > 1000) num -= 1000;
                for (i = 0; i < num; i++) out << ' ';
            }
        if (first) {
            if (dad) {
				out << ((node==dad->neighbors.front()->node) ? fig_char[2] : ((node==dad->neighbors.back()->node) ? fig_char[4] : fig_char[3]));
                for (i = 0; i < abs(br_len); i++)
                    out << ((zero_length) ? zero_char : fig_char[1]);
            }
            if (brtype & WT_INT_NODE)
            	out << node->id;
            else
            	out << fig_char[0];
            if (!node->name.empty())
                out << " (" << node->name << ")";
            if (brtype & WT_BR_LEN && dad)
                out << " " << node->findNeighbor(dad)->length;
            if (brtype & WT_BR_ID && dad)
                out << " [" << node->findNeighbor(dad)->id << "]";
            if (!subtree_br.empty()) {
                if (subtree_br.back() >1000)
                    subtree_br.back() -= 1000;
                else if (subtree_br.back() < 0)
                    subtree_br.back() -= 1000;
            }
        } else {
            if (dad) {
                if (abs(subtree_br.back()) > 1000) out << ' ';
                else out << fig_char[0];
                for (i = 0; i < abs(br_len); i++)
                    out << ' ';
            }
            out << fig_char[0];
        }
        //out << " ";
        //copy (subtree_br.begin(), subtree_br.end(), ostream_iterator<int> (out, " "));
        out << endl;
        first = false;
    }
}

bool MTree::equalTopology(MTree *tree) {
	assert(root->isLeaf());
	Node *root2 = tree->findLeafName(root->name);
	if (!root2) return false;
	ostringstream ostr, ostr2;
	printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
	tree->printTree(ostr2, WT_TAXON_ID | WT_SORT_TAXA, root2);
	return ostr.str() == ostr2.str();
}

void MTree::calcDist(char *filename) {
    vector<string> taxname;
    int i, j;

    // allocate memory
    taxname.resize(leafNum);
    double *dist = new double [leafNum * leafNum];
    // calculate the distances
    calcDist(dist);
    // get the taxa name
    getTaxaName(taxname);

    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);

        // now write the distances in phylip .dist format
        out << leafNum << endl;

        for (i = 0; i < leafNum; i++) {
            out << taxname[i] << "   ";
            for (j = 0; j < leafNum; j++) {
                out << dist[i*leafNum + j] << "  ";
            }
            out << endl;
        }
        out.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
    delete [] dist;
}

void MTree::calcDist(double* &dist, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
        calcDist(node, 0.0, dist, node, NULL);
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
        calcDist(dist, (*it)->node, node);
    }
}

void MTree::calcDist(Node *aroot, double cur_len, double* &dist, Node *node, Node *dad) {
    double branch_length;
	if (!node) node = root;
    if (node->isLeaf()) {
        dist[aroot->id * leafNum + node->id] = cur_len;
        dist[node->id * leafNum + aroot->id] = cur_len;
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
    	branch_length = (*it)->length;
        calcDist(aroot, cur_len + branch_length, dist, (*it)->node, node);
    }

}


/*********************************************
	class PDTaxaSet
*********************************************/

void PDTaxaSet::setSubTree(MTree &tree, NodeVector &subtree) {
    stringstream ostr;
    tree.printSubTree(ostr, subtree);
    tree_str = ostr.str();
}

void PDTaxaSet::setTree(MTree &tree) {
    // assign the taxa set
    tree.getTaxa(*this);
    // assign the score
    score = tree.treeLength();

    // assign tree_str
    stringstream ostr;
    tree.printTree(ostr);
    tree_str = ostr.str();
}


void PDTaxaSet::printTaxa(ostream &out) {
    for (iterator it = begin(); it != end(); it++)
        if ((*it)->name != ROOT_NAME)
            out << (*it)->name << endl;
}

void PDTaxaSet::printTaxa(char *filename) {
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);
        printTaxa(out);
        out.close();
        cout << "Taxa list was printed to " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
}

void PDTaxaSet::printTree(ostream &out) {
    if (!tree_str.empty())
        out << tree_str << endl;
}

void PDTaxaSet::printTree(char *filename) {
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);
        printTree(out);
        out.close();
        cout << "Tree was printed to " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
}


void PDTaxaSet::makeIDSet(int ntaxa, Split &id_set) {
    id_set.setNTaxa(ntaxa);
    id_set.setWeight(score);
    for (iterator it = begin(); it != end(); it++)
        id_set.addTaxon((*it)->id);
}

void MTree::writeInternalNodeNames(string &out_file) {
    try {
        ofstream out(out_file.c_str());
        NodeVector nodes;
        getInternalNodes(nodes);
        for (NodeVector::iterator nit = nodes.begin(); nit != nodes.end(); nit++) {
            out  << " " << (*nit)->name;
        }
        out << endl;
        out.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, out_file);
    }
}

void MTree::assignLeafID(Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
        node->id = atoi(node->name.c_str());
        assert(node->id >= 0 && node->id < leafNum);
    }
    FOR_NEIGHBOR_IT(node, dad, it)
    assignLeafID((*it)->node, node);
}

void MTree::getTaxa(Split &taxa, Node *node, Node *dad) {
	if (!node) node = root;
	if (node->isLeaf()) {
		taxa.addTaxon(node->id);
	}
	FOR_NEIGHBOR_IT(node, dad, it)
		getTaxa(taxa, (*it)->node, node);
}


void MTree::extractQuadSubtrees(vector<Split*> &subtrees, Node *node, Node *dad) {
	if (!node) node = root;
	FOR_NEIGHBOR_IT(node, dad, it) {
		extractQuadSubtrees(subtrees, (*it)->node, node);
		if ((*it)->node->isLeaf()) continue;
		// internal branch
		assert(node->degree() == 3 && (*it)->node->degree() == 3);
		int cnt = 0;
		Node *child = (*it)->node;
		FOR_NEIGHBOR_DECLARE(child, node, it2) {
			Split *sp = new Split(leafNum);
			getTaxa(*sp, (*it2)->node, child);
			subtrees.push_back(sp);
			cnt += sp->countTaxa();
		}
		FOR_NEIGHBOR(node, child, it2) {
			Split *sp = new Split(leafNum);
			getTaxa(*sp, (*it2)->node, node);
			subtrees.push_back(sp);
			cnt += sp->countTaxa();
		}
		assert(cnt == leafNum);
	}
}


void MTree::assignBranchSupport(const char *trees_file) {
	cout << "Reading input trees file " << trees_file << endl;
	try {
		ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(trees_file);
        assignBranchSupport(in);
		in.close();
	} catch (ios::failure) {
		outError(ERR_READ_INPUT, trees_file);
	}
}

void MTree::assignBranchSupport(istream &in) {
	SplitGraph mysg;
	NodeVector mynodes;
	convertSplits(mysg, &mynodes, root->neighbors[0]->node);
	vector<Split*> subtrees;
	extractQuadSubtrees(subtrees, root->neighbors[0]->node);
	IntVector decisive_counts;
	decisive_counts.resize(mynodes.size(), 0);
	StrVector occurence_trees; // list of tree IDs where each split occurs
	if (verbose_mode >= VB_MED)
		occurence_trees.resize(mynodes.size());
	SplitGraph::iterator sit;
	for (sit = mysg.begin(); sit != mysg.end(); sit++)
		(*sit)->setWeight(0.0);
	int ntrees, taxid;
	for (ntrees = 1; !in.eof(); ntrees++) {
		MTree tree;
		bool is_rooted = false;

		// read in the tree and convert into split system for indexing
		tree.readTree(in, is_rooted);
		if (verbose_mode >= VB_DEBUG)
			cout << ntrees << " " << endl;
		StrVector taxname;
		tree.getTaxaName(taxname);
		// create the map from taxa between 2 trees
		Split taxa_mask(leafNum);
		for (StrVector::iterator it = taxname.begin(); it != taxname.end(); it++) {
			taxid = mysg.findLeafName(*it);
			if (taxid < 0)
				outError("Taxon not found in full tree: ", *it);
			taxa_mask.addTaxon(taxid);
		}
		// make the taxa ordering right before converting to split system
		taxname.clear();
		int smallid;
		for (taxid = 0, smallid = 0; taxid < leafNum; taxid++)
			if (taxa_mask.containTaxon(taxid)) {
				taxname.push_back(mysg.getTaxa()->GetTaxonLabel(taxid));
				string name = (string)mysg.getTaxa()->GetTaxonLabel(taxid);
				tree.findLeafName(name)->id = smallid++;
			}
		assert(taxname.size() == tree.leafNum);

		SplitGraph sg;
		//NodeVector nodes;
		tree.convertSplits(sg);
		SplitIntMap hash_ss;
		for (sit = sg.begin(); sit != sg.end(); sit++)
			hash_ss.insertSplit((*sit), 1);

		// now scan through all splits in current tree
		int id, qid;
		for (sit = mysg.begin(), id = 0, qid = 0; sit != mysg.end(); sit++, id++)
		if ((*sit)->trivial() < 0) // it is an internal split
		{

			bool decisive = true;
			for (int i = 0; i < 4; i++) {
				if (!taxa_mask.overlap(*subtrees[qid+i])) {
					decisive = false;
					break;
				}
			}
			qid += 4;
			if (!decisive) continue;

			decisive_counts[id]++;
			Split *subsp = (*sit)->extractSubSplit(taxa_mask);
			if (subsp->shouldInvert())
				subsp->invert();
			Split *sp = hash_ss.findSplit(subsp);
			if (sp && sp->trivial() < 0) {
				(*sit)->setWeight((*sit)->getWeight()+1.0);
				if (verbose_mode >= VB_MED)
					occurence_trees[id] += convertIntToString(ntrees) + " ";
				if (verbose_mode >= VB_MAX) {
					for (taxid = 0; taxid < (*sit)->getNTaxa(); taxid++)
						if ((*sit)->containTaxon(taxid))
							cout << " " << mysg.getTaxa()->GetTaxonLabel(taxid);
					cout << " --> ";
					for (taxid = 0; taxid < sp->getNTaxa(); taxid++)
						if (sp->containTaxon(taxid))
							cout << " " << taxname[taxid];
					cout << endl;
				}
			}
			delete subsp;
		}

		char ch;
		in.exceptions(ios::goodbit);
		(in) >> ch;
		if (in.eof()) break;
		in.unget();
		in.exceptions(ios::failbit | ios::badbit);

	}

	cout << ntrees << " trees read" << endl;

	for (int i = 0; i < mysg.size(); i++)
	if (!mynodes[i]->isLeaf())
	{
		stringstream tmp;
		if (!mynodes[i]->name.empty())
			tmp << "/";
		if (mysg[i]->getWeight() == 0.0)
			tmp << "0";
		else
			tmp << round((mysg[i]->getWeight()/decisive_counts[i])*1000)/10;
		if (verbose_mode >= VB_MED)
			tmp << "%" << decisive_counts[i];
		if (!mynodes[i]->isLeaf()) mynodes[i]->name.append(tmp.str());
		if (verbose_mode >= VB_MED) {
			cout << mynodes[i]->name << " " << occurence_trees[i] << endl;
		}
	}
	for (vector<Split*>::reverse_iterator it = subtrees.rbegin(); it != subtrees.rend(); it++)
		delete (*it);
}

void MTree::computeRFDist(const char *trees_file, IntVector &dist) {
	cout << "Reading input trees file " << trees_file << endl;
	try {
		ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(trees_file);
        computeRFDist(in, dist);
		in.close();
	} catch (ios::failure) {
		outError(ERR_READ_INPUT, trees_file);
	}
}

void MTree::computeRFDist(istream &in, IntVector &dist) {
	SplitGraph mysg;
	convertSplits(mysg, NULL, root->neighbors[0]->node);
	SplitGraph::iterator sit;
	for (sit = mysg.begin(); sit != mysg.end(); sit++)
		(*sit)->setWeight(0.0);
	int ntrees, taxid;
	for (ntrees = 1; !in.eof(); ntrees++) {
		MTree tree;
		bool is_rooted = false;

		// read in the tree and convert into split system for indexing
		tree.readTree(in, is_rooted);
		if (verbose_mode >= VB_DEBUG)
			cout << ntrees << " " << endl;
		StrVector taxname;
		tree.getTaxaName(taxname);
		// create the map from taxa between 2 trees
		Split taxa_mask(leafNum);
		for (StrVector::iterator it = taxname.begin(); it != taxname.end(); it++) {
			taxid = mysg.findLeafName(*it);
			if (taxid < 0)
				outError("Taxon not found in full tree: ", *it);
			taxa_mask.addTaxon(taxid);
		}
		// make the taxa ordering right before converting to split system
		taxname.clear();
		int smallid;
		for (taxid = 0, smallid = 0; taxid < leafNum; taxid++)
			if (taxa_mask.containTaxon(taxid)) {
				taxname.push_back(mysg.getTaxa()->GetTaxonLabel(taxid));
				string name = (string)mysg.getTaxa()->GetTaxonLabel(taxid);
				tree.findLeafName(name)->id = smallid++;
			}
		assert(taxname.size() == tree.leafNum);

		SplitGraph sg;
		//NodeVector nodes;
		tree.convertSplits(sg);
		SplitIntMap hash_ss;
		for (sit = sg.begin(); sit != sg.end(); sit++)
			hash_ss.insertSplit((*sit), 1);

		// now scan through all splits in current tree
		int common_splits = 0;
		for (sit = mysg.begin(); sit != mysg.end(); sit++)
		if ((*sit)->trivial() < 0) // it is an internal split
		{

			Split *subsp = (*sit)->extractSubSplit(taxa_mask);
			if (subsp->shouldInvert())
				subsp->invert();
			Split *sp = hash_ss.findSplit(subsp);
			if (sp) {
				common_splits++;
				//(*sit)->setWeight((*sit)->getWeight()+1.0);
				if (verbose_mode >= VB_MAX) {
					for (taxid = 0; taxid < (*sit)->getNTaxa(); taxid++)
						if ((*sit)->containTaxon(taxid))
							cout << " " << mysg.getTaxa()->GetTaxonLabel(taxid);
					cout << " --> ";
					for (taxid = 0; taxid < sp->getNTaxa(); taxid++)
						if (sp->containTaxon(taxid))
							cout << " " << taxname[taxid];
					cout << endl;
				}
			}
			delete subsp;
		}

		//cout << "common_splits = " << common_splits << endl;
		dist.push_back(common_splits);
		char ch;
		in.exceptions(ios::goodbit);
		(in) >> ch;
		if (in.eof()) break;
		in.unget();
		in.exceptions(ios::failbit | ios::badbit);

	}

	cout << ntrees << " trees read" << endl;


}

void MTree::insertTaxa(StrVector &new_taxa, StrVector &existing_taxa) {
	if (new_taxa.empty()) return;
	IntVector id;
	int i;
	id.resize(new_taxa.size());
	for (i = 0; i < id.size(); i++)
		id[i] = i;
	// randomize order before reinsert back into tree
	my_random_shuffle(id.begin(), id.end());

	for (int i = 0; i < new_taxa.size(); i++) {
		Node *old_taxon = findLeafName(existing_taxa[id[i]]);
		assert(old_taxon);
		double len = old_taxon->neighbors[0]->length;
		Node *old_node = old_taxon->neighbors[0]->node;
		Node *new_taxon = newNode(leafNum+i, new_taxa[id[i]].c_str());
		Node *new_node = newNode();
		// link new_taxon - new_node
		new_taxon->addNeighbor(new_node, 0.0);
		new_node->addNeighbor(new_taxon, 0.0);
		// link old_taxon - new_node
		new_node->addNeighbor(old_taxon, 0.0);
		old_taxon->updateNeighbor(old_node, new_node, 0.0);
		// link old_node - new_node
		new_node->addNeighbor(old_node, len);
		old_node->updateNeighbor(old_taxon, new_node, len);
	}

    leafNum = leafNum + new_taxa.size();
    initializeTree();
}
