#include <iostream>
#include <fstream>
#include <getopt.h>
#include <map>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <sstream>
#include <list>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <boost/graph/edmunds_karp_max_flow.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/vector_as_graph.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include "Contig.h"

int DEBUG = 0;
double COOPERATIVE_BURY_THRESHOLD = 80; //used for burying clone into multiple clones..
double BURY_THRESHOLD = 80;
//double EXTRA_FRAGMENT_PERCENTAGE = 50; //if extra fragments of a node (i.e. match nowhere in clone - fragment graph) is more than this percentage then they are added to clone-fragment graph...
double MTP_MST_VERY_HIGH_THRESHOLD = 1; // 3 for barley, 1 for rice.. if a component contains a pair that has an overlap_score_exponent higher than this value then that component is splitted..


void usage() {
	cout << "Compute_MTP: Generates a linear programming model to select MTP clones of contigs in a FPC map..\n";
	cout << "-c ctg_BAC_clone file\n";
	cout << "-f Fingerprinting Method.. 1: HICF, 2: Agarose\n";	
	cout << "-s size file\n";
	cout << "-b bury threshold (default 80)\n";
	cout << "-B cooperative bury threshold (default 80)\n";
	cout << "-t sulston_score_threshold.. Fragments of clones will be matched if their sulston score is lower than this threshold..(Default values agarose:1e-2, HICF:1e-5)\n";
	cout << "-T Tolerance (rice:7, barley:3, cowpea:5)..\n";
	cout << "-G Gellength (rice: 3300, barley:18000, cowpea: 36,000)..\n";
	cout << "-o output_file base file (without extension).. Linear Programming Model file (for GLPK, in MathProg language)\n";
	cout << "-O output ctg 2 ordered_clones (OPTIONAL)\n";
	cout << "-e putative ovErlapping pairs File.. when MST is built.. it's check which clones are overlapping and which are not.. putative overlapping pairs are stored here..(OPTIONAL)\n";
	cout << "-n putative non-overlapping  pairs File.. when MST is built, it's check which clones are overlappig and which are not.. putative non-overlapping pairs are stored here... (OPTIONAL)\n";
	cout << "-v verbose level min 0 max 4\n";
	exit(0);
}


//TODO: Following functions are added for testing purposess only...
//debug-begin
//struct cloneCmp {
//	bool operator()(const Clone *c1, const Clone *c2) const {
//		return c1->getName().compare(c2->getName()) < 0;
//	}
//};
//
//
//typedef pair<int, pair<int, int> > coordinate;
//
//int is_overlapping(coordinate & clone1, coordinate & clone2){
//	int result = 0;
//	if (clone1.first == clone2.first){
//		//cout << clone1.second.first << "\t" << clone1.second.second << endl;
//		//cout << clone2.second.first << "\t" << clone2.second.second << endl;
//		if ((clone1.second.first <= clone2.second.first) && (clone1.second.second >= clone2.second.first)){
//			//cout << " clone 1 overlaps clone 2.. " << endl;
//			result = clone1.second.second - clone2.second.first + 1;
//		}else if ((clone2.second.first <= clone1.second.first) && (clone2.second.second >= clone1.second.first)){
//			//cout << "clone 2 overlaps clone 1 " << endl;
//			result = -(clone2.second.second - clone1.second.first);
//		}else{
//			//cout << "no overlap.." << endl;
//			result = 0;
//		}
//	}else{
//		//no overlap..
//		//cout << "no overlap .." << endl;
//		result = 0;
//	}
//	return result;
//}
//
//
//int is_buried(coordinate & clone1, coordinate & clone2){
//	int result = 0;
//	if (clone1.first == clone2.first){
//		if ((clone1.second.first >= clone2.second.first) && (clone1.second.second <= clone2.second.second)){
//			//clone 1 is buried in clone 2..
//			result = 1;
//		}else if ((clone2.second.first >= clone1.second.first) && (clone2.second.second <= clone1.second.second)){
//			//clone 2 is buried into clone 1..
//			result = -1;
//		}else{
//			//no bury..
//			result = 0;
//		}
//	}else{
//		result = 0;
//	}
//	return result;
//}
//debug-end

struct graph_writer {
	void operator()(std::ostream& out) const {
		out << "graph [size=\"6,6\"]" << std::endl;

		//out << "graph [rankdir=\"TB\"]" << std::endl;
	}
};

struct strCmp {
	bool operator()(const string s1, const string s2) const {
		return s1.compare(s2) < 0;
	}
};

struct pairCmp {
	bool operator()(const pair<string, string> p1, const pair<string, string> p2) const {
		if (p1.first.compare(p2.first) == 0) {
			return p1.second.compare(p2.second) < 0;
		} else {
			return p1.first.compare(p2.first) < 0;
		}
	}
};

struct pairCmp2 {
	bool operator()(const pair<string, int> p1, const pair<string, int> p2) const {
		if (p1.first.compare(p2.first) == 0) {
			return p1.second < p2.second;
		} else {
			return p1.first.compare(p2.first) < 0;
		}
	}
};

struct pairCmp3 {
	bool operator()(const pair<int, int> p1, const pair<int, int> p2) const {
		if (p1.first == p2.first) {
			return p1.second < p2.second;
		} else {
			return p1.first < p2.first;
		}
	}
};

bool isNumeric(const char * s) {
	char * p;
	strtol(s, &p, 10);
	return !*p;
}

namespace boost {
struct edge_component_t {
	enum
	{	num = 555};
	typedef edge_property_tag kind;
} edge_component;
}

int get_copy_id(string cloneName, int size);

typedef map<string, Vertex, strCmp > NameVertexMap;

map<pair<string, int>, int, pairCmp2> fragment_2_copy_id;

typedef adjacency_list < vecS, vecS, undirectedS, property < vertex_name_t, std::string>, property < edge_component_t, std::size_t , property <edge_weight_t, double> > >
		U_Graph;
typedef graph_traits<U_Graph>::vertices_size_type size_type;

typedef property_map < U_Graph, vertex_name_t >::type U_VertexName;
typedef property_map < U_Graph, edge_weight_t >::type U_EdgeWeight;
typedef U_Graph::vertex_descriptor U_Vertex;
typedef U_Graph::edge_descriptor U_Edge;
typedef map<string, U_Vertex, strCmp > NameUVertexMap;

typedef property_map<Graph, Band * VertexProperties::*>::type BandPtr;

map<string, string, strCmp> name_2_color;//map vertex names to colors (each clone will be colored different)...

template<class Name> class custom_label_writer {
public:
	custom_label_writer(Name _name) :
		name(_name) {
	}
	template<class U_Vertex> void operator()(std::ostream& out, const U_Vertex& v) const {
		string clone_name = name[v].substr(0, 8);
		string color = name_2_color[clone_name];
		out << "[label=\"" << name[v] << "\",color=\"" << color
				<< "\",style=filled]";
	}
private:
	Name name;
};

template<class Name> class custom_label_writer_2 {
public:
	custom_label_writer_2(Name _name) :
		name(_name) {
	}
	template<class BandPtr> void operator()(std::ostream& out, const BandPtr& b) const {
		string fragment_name = name[b]->getName();
		//string color = name_2_color[clone_name];
		out << "[label=\"" << fragment_name << "\"]";
	}
private:
	Name name;
};

string min_str(const string & s_1, const string & s_2) {
	if (s_1.compare(s_2) < 0) {
		return s_1;
	} else {
		return s_2;
	}
}

string max_str(const string &s_1, const string & s_2) {
	if (s_1.compare(s_2) > 0) {
		return s_1;
	} else {
		return s_2;
	}
}

/*
 * This class is used to store in a priority queue.. frequency here is # times that a clone occurs in connected components of the MFG..
 */
class CloneFrequency{
	string clone_name;
	int frequency;
public:
	CloneFrequency(string c, int f): clone_name(c), frequency(f) {}
	friend bool operator<(const CloneFrequency& c1, const CloneFrequency &c2) {
		return c1.frequency < c2.frequency;
	}
	string getCloneName() const {return clone_name;}
};

//function prototypes...
void compute_redundant_components(const vector<int> &, const U_VertexName &,
		set<int> &);
string get_clone_name(string);
void compute_clone_component_size(const vector<int>&, const U_VertexName &,
		map<int, int, less<int> > &);
int get_fragment_size(string);
bool is_special_component(const set<int> &, const U_VertexName &);
bool is_fragment_values_not_within_tolerance(const set<int> &,
		const U_VertexName &, const int);
bool
		is_same_clone_together_in_component(const set<int> &,
				const U_VertexName &);
//bool is_loosely_connected_component(const set<int> &, const U_Graph &,
//		const U_VertexName &);
void get_vertices_per_component(const vector<int>&,
		map<int, set<int>, less<int> > &);
double split_component(const set<int> &, const U_Graph &, const U_EdgeWeight &,
		map<int, int, less<int> > &);
int get_max_component_size(map<int, int, less<int> > &, const vector<int> &,
		const U_VertexName &);
void compute_connected_clone_list(map<string, int, strCmp> &,
		map<string, set<string>, strCmp> &);
void tokenize(const string&, vector<string>&, const string& delimiters = "-");
void concatenate_clones(const set<string> &, string &,
		const string& delimiters = "-");
void get_clone_components(const vector<int> &, const U_VertexName &,
		map<string, int, strCmp> &);
void sort_clones_based_on_overlap_score(const string &source_clone,
		const set<string> &connected_clones, vector<string> &sorted_clones,
		map<pair<string, string>, double, pairCmp > & overlap_score_exponent);
void
		get_adjacency_map(
				map<string, set<string>, strCmp >& connected_clones_map,
				map<pair<string, string>, double, pairCmp >& overlap_score_exponent,
				map<pair<string, string>, int, pairCmp> &n_shared_fragments_between_clones,
				map<string, set<string>, strCmp> & adjacency_map,
				const string contig_id);
//void prune_connected_clones_map(
//		map<string, set<string>, strCmp > & connected_clones_map,
//		map<pair<string, string>, double, pairCmp > & overlap_score_exponent);
void traverse_clones(map<string, set<string>, strCmp > & connected_clones_map,
		map<pair<string, string>, double, pairCmp > & overlap_score_exponent,
		vector< vector<string> > & ordered_clones);
void
		get_n_shared_components_between_clones(
				const vector<int> & vertex_ID_2_component_ID,
				const U_VertexName & v_name,
				map<pair<string, string>, int, pairCmp> &n_shared_components_between_clones);
void assign_to_component_map(string clone,
		map<string, set<string>, strCmp > & adjacency_map,
		map<int, set<string>, less<int> > & component_2_clone_map, int index);

bool is_together(map<int, set<int>, less<int> > & component_ID_to_vertices, const U_VertexName & v_name, const set<string> &clones);
bool exists_weak_edge(const set<int> &vertices_in_component, const U_VertexName & v_name, map<pair<string, string>, double, pairCmp > & overlap_score_exponent);
void compute_clones_buried_into_multiple_clones(const vector<int> &vertex_ID_2_component_ID, const U_VertexName &v_name, map<string, Clone *, strCmp> &cloneObjectByName, vector<string> &buried_clones);
void count_number_of_matching_fragments(map<int, set<int>,less<int> > &component_ID_to_vertices, const U_VertexName &v_name, map<string, int, strCmp> &clone_2_n_matching_fragments);
void remove_clone_from_MFG(U_Graph &g, string &clone_name, U_VertexName &v_name);
void generateCandidateSet(const vector<int>& vertex_ID_2_component_ID, const U_VertexName &v_name, const vector<string> &candidate_clones_for_bury, vector<set<string> > & candidate_clone_set, set<string> &clones_in_candidate_clone_set);
void compute_minimum_hitting_set(vector<set<string> > &candidate_clone_set, set<string> &minimum_hitting_set);
bool is_candidate_clone_for_bury(const string &clone_name, const vector<string> candidate_clones_for_bury);


int main(int argc, char**argv) {
	map<string, Clone *, strCmp> cloneObjectByName; // stores the pointer to the Clone object of each clone (by clone name)
	map<string, Contig *, strCmp> contigObjectByName; // stores the pointer to the Contig object of each contig (by contig name)
	map<string, Band *, strCmp> fragmentObjectByName; //stores the pointer to the Fragment object of each fragment (by fragment name)..

	char ctgBACFile[FILENAME_MAX];
	char sizeFile[FILENAME_MAX];
	char outputFile[FILENAME_MAX];
	char outputOrderedClones[FILENAME_MAX];
	char outputPutativeOverlappingClones[FILENAME_MAX];
	char outputPutativeNonOverlappingClones[FILENAME_MAX];

	int ctgBACFilePresent = 0, sizeFilePresent = 0,
			outputOrderedClonesFilePresent = 0,
			output_overlapping_clones_file_present = 0,
	output_non_overlapping_clones_file_present = 0;
	double sulston_score_threshold = -1.0;
	int fingerprinting_method = -1;
	int tolerance = -1, gellen = -1;
	
//	//TODO: added for testing purposes only...
//	char clone_coordinates_file[FILENAME_MAX];

	
	map<pair<string, string>, double, pairCmp > overlap_score_exponent;//stores -ln(sulston(A,B)) between all pairs of clones (except for buried clones, which are totally ignored in the analysis)...
	map<pair<string, string>, int, pairCmp> E_f; //stores the edges between fragments.. (if an fragment is covered by selecting its clone, all fragments adjacent to this fragment in E_f are also covered..
	map<pair<string, string>, int, pairCmp> E; //stores the edges between clones and fragments..
	//map<string, int, strCmp> fragmentDegree; //stores the degree of a fragment.. if degree = 1.. it's connected to a clone and a fragment, if degree = 0, only connected to a clone..
	while (1) {
		char c = getopt(argc, argv, "hc:f:s:b:B:t:T:G:o:O:v:e:n:");

		if (c == -1)
			break;

		switch (c) {
		case 'h':
			usage();
			break;

		case 'c':
			strcpy(ctgBACFile, optarg);
			printf("Using ctgBACFile in %s\n", ctgBACFile);
			ctgBACFilePresent = 1;
			break;
		case 'f':
			fingerprinting_method = atoi(optarg);
			printf("Fingerprinting Method Code: %d\n", fingerprinting_method);
			break;
		case 's':
			strcpy(sizeFile, optarg);
			printf("Using size file %s\n", sizeFile);
			sizeFilePresent = 1;
			break;
		case 'b':
			BURY_THRESHOLD = atof(optarg);
			printf("BURY_THRESHOLD : %f\n", BURY_THRESHOLD);
			break;
		case 'B':
			COOPERATIVE_BURY_THRESHOLD = atof(optarg);
			printf("COOPERATIVE_BURY_THRESHOLD: %f\n", COOPERATIVE_BURY_THRESHOLD);
			break;
		case 't':
			sulston_score_threshold = atof(optarg);
			printf("Sulston score threshold: %e\n", sulston_score_threshold);
			break;
		case 'T':
			tolerance = atoi(optarg);
			printf("Using tolerance value %d\n", tolerance);
			break;
		case 'G':
			gellen = atol(optarg);
			printf("Using gel length value %d\n", gellen);
			break;
		case 'o':
			strcpy(outputFile, optarg);
			printf("Will write output data to '%s'\n", outputFile);
			break;
		case 'O':
			strcpy(outputOrderedClones, optarg);
			printf("Will write ordered clones to '%s'\n", outputOrderedClones);
			outputOrderedClonesFilePresent = 1;
			break;
		case 'v':
			DEBUG = atoi(optarg);
			printf("Verbose level: %d\n", DEBUG);
			break;
//		case 'r':
//			criteria = atoi(optarg);
//			printf("Using criteria: %d\n", criteria);
//			break;
		case 'e':
			strcpy(outputPutativeOverlappingClones, optarg);
			printf("Using putative overlapping clones file: %s\n",
					outputPutativeOverlappingClones);
			output_overlapping_clones_file_present = 1;
			break;
		case 'n':
			strcpy(outputPutativeNonOverlappingClones, optarg);
			printf("Using putative non-overlapping clones file: %s\n",
					outputPutativeNonOverlappingClones);
			output_non_overlapping_clones_file_present = 1;
			break;
//		case 'x'://TODO: added for testing purposes only...
//			strcpy(clone_coordinates_file, optarg);
//			printf ("using coordinate file %s", clone_coordinates_file);
//			break;
		default:
			printf(
					"getopt() returned character code 0%o - not handled in code. Useage instructions follow.\n",
					c);
			usage();
		}
	}

	if (!sizeFilePresent || !ctgBACFilePresent) {
		cerr << "input file missing\n";
		usage();
	}

	//	if (strlen(outputFile) == 0 || strlen(outputCtgMTPFile) == 0) {
	//		cerr << "output file missing\n";
	//		usage();
	//	}

	if (strlen(outputFile) == 0) {
		cerr << "output file missing\n";
		usage();
	}

	
	if (fingerprinting_method == HICF){
		MTP_MST_VERY_HIGH_THRESHOLD = MTP_MST_VERY_HIGH_THRESHOLD_HICF;
		//update connected_component_size_threshold if it's not set by user..
//		if (p_value == -1){
//			p_value = MTP_MST_HICF_Q;
//		}
		//update sulston_score_threshold if not set by the user..
		if (sulston_score_threshold == -1){
			sulston_score_threshold = MTP_MST_HICF_C;
		}
	}else if (fingerprinting_method == AGAROSE){
		MTP_MST_VERY_HIGH_THRESHOLD = MTP_MST_VERY_HIGH_THRESHOLD_AGAROSE;
		//update connected_component_size_threshold if it's not set by user..
//		if (p_value == -1){
//			p_value = MTP_MST_HICF_Q;
//		}
		//update sulston_score_threshold if not set by the user..
		if (sulston_score_threshold == -1){
			sulston_score_threshold = MTP_MST_AGAROSE_C;
		}
	}else{
		cerr << "Invalid fingerprinting code\n";
		usage();
	}
	
	if (sulston_score_threshold< 0 || sulston_score_threshold> 1){
		cerr << "invalid sulston score threshold.. must be between 0 and 1\n";
		usage();
	}
	
	if (tolerance <= 0 || gellen <= 0){
		cerr << "Invalid tolerance/gellen...\n";
		usage();
	}

	if (BURY_THRESHOLD < 0 || BURY_THRESHOLD > 100){
		cerr << "Invalid range for BURY_THRESHOLD.. Should be between [0,100]";
		usage();
	}
	
	
	if (output_overlapping_clones_file_present*output_non_overlapping_clones_file_present < 0){
		cerr << "you need both files at the same time..\n";
		usage();
	}
	
	//read the input file...
	ifstream inCtgBACFile(ctgBACFile, ios::in);

	if (!inCtgBACFile) {
		cerr << "cannot open input file\n";
		exit(1);
	}

	//read the ctg-BAC file and generate clone & contig objects
	char line[90000];
	char cloneName[20];
	char contigName[10];
	int nClones;

	while (inCtgBACFile.getline(line, 90000)) {
		char * tok;
		tok = strtok(line, "\t");
		strcpy(contigName, tok);
		Contig * ctgPtr = new Contig(contigName);
		contigObjectByName[contigName] = ctgPtr;

		tok = strtok(NULL, "\t");
		nClones = atoi(tok);
		//cout << "nClones: " << nClones << endl;
		tok = strtok(NULL, "\t");

		char cloneNames[90000] = "";
		strcpy(cloneNames, tok);
		//cout << cloneNames << endl;
		char * tok2;
		tok2 = strtok(cloneNames, " ");

		while (tok2 != NULL) {
			strcpy(cloneName, tok2);
			//cout << "clonename: " << cloneName << endl;
			if (strlen(cloneName) == 0) {
				cout << "clonename length is zero" << endl;
				exit(1);
			}
			Clone * clonePtr = new Clone(cloneName);
			cloneObjectByName[cloneName] = clonePtr;

			ctgPtr->addToClones(clonePtr);

			tok2 = strtok(NULL, " ");

		}

	}

	cout << "ctg-BAC file is read\n";

	//read the size file and store Fragments...
	ifstream inSizeFile(sizeFile, ios::in);
	if (!inSizeFile) {
		cerr << "Cannot read the size file, " << sizeFile << endl;
		exit(1);
	}

	vector<Band *> tmpFragments;
	bool skip = false;//used for ignoring clones that are not in cloneObjectByName map....
	while (inSizeFile.getline(line, 500)) {
		char * tok;
		tok = strtok(line, "\t");
		//not need the rest of the tokens, if there is any..

		//if the value is a band size of an skipped clone then ignore..
		if (isNumeric(tok) && !skip) {
			int n = atoi(tok);
			if (n == -1) {
				//store this clone..
				Clone * clonePtr = cloneObjectByName[cloneName];

				for (vector<Band *>::iterator iter = tmpFragments.begin(); iter
						!= tmpFragments.end(); iter++) {
					clonePtr->addToFragments(*iter);
				}
				tmpFragments.clear();
			} else {
				//nTotalBands++;
				//discretize(n);
				int copy_id = get_copy_id(cloneName, n);
				Band *b = new Band(cloneName, n, copy_id);
				fragmentObjectByName[b->getName()] = b;
				tmpFragments.push_back(b);
			}
		} else {
			strcpy(cloneName, line);
			/*ignore this clone if it is not in cloneObjectByName map... 
			 (because size file may contain clones whose coordinates are not available)..
			 */
			if (cloneObjectByName.find(cloneName) == cloneObjectByName.end()) {
				skip = true;//all band sizes of this clone will be skipped..
			} else {
				skip = false;
			}
		}
	}

	cout << "size file is read\n";

	if (DEBUG > 3) {
		cout << "Priniting ctg-BAC table\n";
		for (map<string, Contig *, strCmp>::iterator iter =
				contigObjectByName.begin(); iter != contigObjectByName.end(); iter++) {
			iter->second->print();
		}
	}
	
	
//	//TODO following section is added for debugging purposess only...
//	//DEBUG--BEGIN
//	
//	map<Clone *, coordinate, cloneCmp> clone_2_coordinates; //stores the coordinates of clones (chr_no, start, end)..
//	ifstream inCloneCoordinatesFile(clone_coordinates_file, ios::in);
//	
//	if (!inCloneCoordinatesFile) {
//		cerr << "cannot open input file...\n";
//		exit(1);
//	}
//	
//	cout << "reading coordinates...\n";
//	int chr_no;
//	long start, end;
//	while (inCloneCoordinatesFile >> chr_no >> cloneName >> start >> end) {
//		//cout << "cloneName : " << cloneName << endl;
//		if (cloneObjectByName.find(cloneName) != cloneObjectByName.end()){
//			coordinate c;
//			c.first = chr_no;
//			c.second.first = start;
//			c.second.second = end;
//			//cout << cloneName << "\t" << chr_no << "\t" << start << "\t" << end << endl;
//			
//			clone_2_coordinates[cloneObjectByName[cloneName]] = c;
//		}
//	}
//
//	
//	//DEBUG--END

	
	//remove buried clones...
	for (map<string, Contig *, strCmp>::iterator con_iter =
			contigObjectByName.begin(); con_iter != contigObjectByName.end(); con_iter++) {
		set<Clone *, Contig::cloneCmp> clones = con_iter->second->getClones();
		for (set<Clone *, Contig::cloneCmp>::iterator cl_iter_1 =
				clones.begin(); cl_iter_1 != clones.end(); cl_iter_1++) {
			multiset<int> f_clone_1 = (*cl_iter_1)->getFragmentSizes();
			for (set<Clone *, Contig::cloneCmp>::iterator cl_iter_2 =
					clones.upper_bound(*cl_iter_1); cl_iter_2 != clones.end(); cl_iter_2++) {
				multiset<int> f_clone_2 = (*cl_iter_2)->getFragmentSizes();

				int total_n_shared = numberOfSharedBands(tolerance, f_clone_1, f_clone_2);
				double per;

				if (DEBUG > 3) {
					cout << (*cl_iter_1)->getName() << ": " << f_clone_1.size()
							<< "\t" << (*cl_iter_2)->getName() << ": " << f_clone_2.size()
							<< "\t" << total_n_shared << endl;
				}
				if (f_clone_1.size() > f_clone_2.size()) { //try to bury clone_2
					per = (double)total_n_shared/(double)f_clone_2.size()*100;
					if (per >= BURY_THRESHOLD) {
						if (DEBUG > 3) {
							cout << (*cl_iter_2)->getName() << " buried into "
									<< (*cl_iter_1)->getName() << endl;
						}
						(*cl_iter_2)->setParent(*cl_iter_1);
					}
				} else {//try to bury clone_1
					per = (double)total_n_shared/(double)f_clone_1.size()*100;
					if (per >= BURY_THRESHOLD) {
						if (DEBUG > 3) {
							cout << (*cl_iter_1)->getName() << " buried into "
									<< (*cl_iter_2)->getName() << endl;
						}
						(*cl_iter_1)->setParent(*cl_iter_2);
					}
				}
			}
		}
	}

	//ofstream ctgMTPFile(outputCtgMTPFile, ios::out);
	//if (!ctgMTPFile) {
	//	cerr << "Cannot open ctg_2_mtp file\n";
	//	exit(1);
	//}

	ofstream outFile(outputFile, ios::out);
	if (!outFile) {
		cerr << "cannot open output file \n";
		exit(1);
	}

	ofstream outputOrderedClonesFile;
	if (outputOrderedClonesFilePresent){
		outputOrderedClonesFile.open(outputOrderedClones, ios::out);
		if (!outputOrderedClonesFile) {
			cerr << "cannot open output file" << outputOrderedClones << "\n";
			exit(1);
		}
	}
	ofstream outputOverlappingClonesFile, outputNonOverlappingClonesFile;
	
	if (output_overlapping_clones_file_present){
		outputOverlappingClonesFile.open(outputPutativeOverlappingClones, ios::out);
		if (!outputOverlappingClonesFile){
			cerr << "cannot open output file " << outputOverlappingClonesFile << endl;
			exit(1);
		}
		
		outputNonOverlappingClonesFile.open(outputPutativeNonOverlappingClones, ios::out);
		if (!outputNonOverlappingClonesFile){
			cerr << "cannot open output file " << outputNonOverlappingClonesFile << endl;
			exit(1);
		}
	}
	
	//Go over all contigs and find MTP for each of them....
	for (map<string, Contig *, strCmp>::iterator con_iter =
			contigObjectByName.begin(); con_iter != contigObjectByName.end(); con_iter++) {
		cout << "contig id: " << con_iter->first << endl;
		set<Clone *, Contig::cloneCmp> clones = con_iter->second->getClones();
		vector<string> mtp_clones;//stores the mtp clones for the current contig...
		
		int contig_size = 0;
		//map<string, bool, strCmp> clone_2_ismatched;//stores if a clone matches to any other clone in the contig...
		for (set<Clone *, Contig::cloneCmp>::iterator cl_iter_1 =
				clones.begin(); cl_iter_1 != clones.end(); cl_iter_1++) {
			//cout << "clonename : " << (*cl_iter_1)->getName() << endl;
			if ((*cl_iter_1)->isBuried()) {
				continue;
			}
			contig_size++;
			//clone_2_ismatched[(*cl_iter_1)->getName()] = false;//initialization...
			//if ((*cl_iter_1)->isBuried()){cout << "buried " << (*cl_iter_1)->getName() << endl;;continue;}

			//adding edges from clone to its fragments...
			vector <Band *> fragments = (*cl_iter_1)->getFragments();
			for (unsigned int i = 0; i < fragments.size(); i++) {
				pair<string, string> p((*cl_iter_1)->getName(), fragments[i]->getName());
				E[p] = 1;
			}
		}
		
		if (contig_size == 2){
			outFile << con_iter->first << "\t" << contig_size << "\t";
			for (set<Clone *, Contig::cloneCmp>::iterator cl_iter_1 = clones.begin(); cl_iter_1 != clones.end(); cl_iter_1++) {	
				if ((*cl_iter_1)->isBuried()) {continue;}
				outFile << (*cl_iter_1)->getName() << " ";
			}
			outFile << "\n";
			continue;
		}
		

		//compare canonical clones (all to all)

		for (set<Clone *, Contig::cloneCmp>::iterator cl_iter_1 =
				clones.begin(); cl_iter_1 != clones.end(); cl_iter_1++) {
			if ((*cl_iter_1)->isBuried()) {
				continue;
			}
			//bool ismatched = false;
			vector<Band *> bands_1 = (*cl_iter_1)->getFragments();

			//set a color for this clone (used in g_f)....
			int r = 0, g = 0, b = 0;
			r = rand() % 255;
			g = rand() % 255;
			b = rand() % 255;
			stringstream color;
			color << "#";
			//cout.width(2);
			//cout << "r.............";
			color << hex << setw(2) << setfill('0') << r;
			color << hex << setw(2) << setfill('0') << g;
			color << hex << setw(2) << setfill('0') << b;

			//cout << "color: " << color.str() << endl;
			name_2_color[(*cl_iter_1)->getName()] = color.str();
			for (set<Clone *, Contig::cloneCmp>::iterator cl_iter_2 =
					clones.upper_bound(*cl_iter_1); cl_iter_2 != clones.end(); cl_iter_2++) {
				if ((*cl_iter_2)->isBuried()) {
					continue;
				}
				//cout << "clonename_2 : " << (*cl_iter_2)->getName() << endl;
				//compute the sulston score ...
				double score;
				int nShared = numberOfSharedBands(tolerance, (*cl_iter_1)->getFragmentSizes(), (*cl_iter_2)->getFragmentSizes());
				if (nShared == 0) {
					score = 1.0;
				} else {
					fastSulston((*cl_iter_1)->getNFragments(), (*cl_iter_2)->getNFragments(), nShared, tolerance, gellen, score);
				}

				//if the sulston score is lower than the threshold then match the fragments of these clones and store each match in E_fragment..

				//store the score in overlap_score_exponent...We may need overlap score exponent of A & B, event they don't seem to overlapping for this threshold, they may be in the same component (i.e. connected, implicitly overlapping)
				pair<string, string> p;
				p.first = min_str((*cl_iter_1)->getName(), (*cl_iter_2)->getName());
				p.second = max_str((*cl_iter_1)->getName(), (*cl_iter_2)->getName());
				overlap_score_exponent[p] = -log10(score);

				if (score <= sulston_score_threshold) {
					//if (DEBUG > 3) {

					//cout << overlap_score_exponent[p] << endl;
					//clone_2_ismatched[(*cl_iter_1)->getName()] = true;
					//clone_2_ismatched[(*cl_iter_2)->getName()] = true;
					//ismatched = true;
					//generate nodes for each fragment of both clones..
					Graph g;
					NameVertexMap L, R;

					BandPtr band_ptr = get(&VertexProperties::band_ptr, g);

					property_map < Graph, edge_capacity_t >::type capacity =
							get(edge_capacity, g);
					property_map < Graph, edge_reverse_t >::type rev = get(
							edge_reverse, g);
					property_map < Graph, edge_residual_capacity_t >::type
							residual_capacity = get(edge_residual_capacity, g);
					//property_map<Graph, int VertexProperties::*>::type size = get(&VertexProperties::size, g);

					vector<Band *> bands_1 = (*cl_iter_1)->getFragments();
					vector<Band *> bands_2 = (*cl_iter_2)->getFragments();
					//cout << "clonename_1 : " << (*cl_iter_1)->getName() << " sizes" << endl;
					for (unsigned int i = 0; i < bands_1.size(); i++) {
						Vertex u;
						NameVertexMap::iterator pos;
						bool inserted;
						tie(pos, inserted) = L.insert(std::make_pair(bands_1[i]->getName(), Vertex()));
						if (inserted) {
							u = add_vertex(g);
							band_ptr[u] = bands_1[i];
							pos->second = u;
						} else {
							u = pos->second;
						}
					}
					//cout << "clonename_2 : " << (*cl_iter_2)->getName() << " sizes" << endl;
					for (unsigned int i = 0; i < bands_2.size(); i++) {
						Vertex u;
						NameVertexMap::iterator pos;
						bool inserted;
						tie(pos, inserted) = R.insert(std::make_pair(bands_2[i]->getName(), Vertex()));
						if (inserted) {
							u = add_vertex(g);
							band_ptr[u] = bands_2[i];
							pos->second = u;
						} else {
							u = pos->second;
						}
					}
					//add the edges...
					unsigned int lstart = 0;
					for (unsigned int i = 0; i < bands_1.size(); i++) {
						bool updated = false;
						for (unsigned int j = lstart; j < bands_2.size(); j++) {
							int diff = bands_1[i]->getSize()
									- bands_2[j]->getSize();
							if (abs(diff) <= tolerance) {
								if (!updated) {
									lstart = j;
									updated = true;
								}
								//add an edge between i & j..
								Edge e1, e2;
								Vertex l, r;
								l = L[bands_1[i]->getName()];
								r = R[bands_2[j]->getName()];
								bool b1, b2;
								tie(e1, b1) = add_edge(l, r, g);
								tie(e2, b2) = add_edge(r, l, g);
								if (!b1 || !b2) {
									cerr << "cannot add edges\n";
									exit(0);
								}
								capacity[e1] = 1;
								capacity[e2] = 0;
								rev[e1] = e2;
								rev[e2] = e1;
							} else if (diff < 0) {
								break;
							}
						}
					}
					if (DEBUG > 3) {
						stringstream ssGraphvizOutputFileName_t;
						ssGraphvizOutputFileName_t << outputFile;

						string graphizOutputFileName_t =
								ssGraphvizOutputFileName_t.str() + "-"
										+ con_iter->first + "-" + (*cl_iter_1)->getName() + "-" + (*cl_iter_2)->getName() + "_mbp.dot";
						ofstream graphvizOutFile_t(
								graphizOutputFileName_t.c_str(), ios::out);
						if (!graphvizOutFile_t) {
							cerr << "cannot open output file \n";
							exit(1);
						}

						write_graphviz(graphvizOutFile_t, g,
								custom_label_writer_2<BandPtr>(band_ptr),
								default_writer(), graph_writer());
						graphvizOutFile_t.close();
					}
					//print out the whole graph...
					//Graph::edge_iterator i, i_end;
					//for (tie(i, i_end) = edges(g); i != i_end; ++i){
					/*Vertex s = boost::source(*i, g);
					 Vertex t = boost::target(*i, g);*/

					//cout << band_ptr[s]->getSize() << " - " << band_ptr[t]->getSize() << endl;	
					//}


					//compute the maximal matching...
					//first add a source and sink node..
					Vertex source = add_vertex(g);
					Vertex target = add_vertex(g);
					band_ptr[source] = NULL;
					band_ptr[target] = NULL;
					//add an edge from source to L vertices..
					NameVertexMap::iterator iter;
					for (iter = L.begin(); iter != L.end(); iter++) {
						Edge e1, e2;
						bool b1, b2;
						tie(e1, b1) = add_edge(source, iter->second, g);
						tie(e2, b2) = add_edge(iter->second, source, g);
						if (!b1 || !b2) {
							cerr << "cannot edge edge from source to L nodes\n";
							exit(0);
						}
						capacity[e1] = 1;
						capacity[e2] = 0;
						rev[e1] = e2;
						rev[e2] = e1;
					}
					for (iter = R.begin(); iter != R.end(); iter++) {
						Edge e1, e2;
						bool b1, b2;
						tie(e1, b1) = add_edge(iter->second, target, g);
						tie(e2, b2) = add_edge(target, iter->second, g);
						if (!b1 || !b2) {
							cerr << "cannot add edge from R to target\n";
							exit(0);
						}
						capacity[e1] = 1;
						capacity[e2] = 0;
						rev[e1] = e2;
						rev[e2] = e1;
					}

					//compute the maximal matching

					int match = (int)edmunds_karp_max_flow(g, source, target);
					if (DEBUG > 3) {
						cout << (*cl_iter_1)->getName() << "-" << (*cl_iter_2)->getName() << "\tmatch : " << match
								<< endl;
					}
					Graph::edge_iterator ei, ei_end;
					for ( tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
						if (band_ptr[boost::source(*ei, g)] != NULL && band_ptr[boost::target(*ei, g)] != NULL) {
							if (residual_capacity[*ei] == 0 && capacity[*ei]
									> 0) {
								Vertex s = boost::source(*ei, g);
								Vertex t = boost::target(*ei, g);

								//add matched fragments to E_fragment..
								pair<string, string> p(band_ptr[s]->getName(),
										band_ptr[t]->getName());
								if (E_f.find(p) != E_f.end()) {
									cerr
											<< "Error: fragment pair is already added..\n";
									exit(0);
								}
								E_f[p] = 1;
							}
						}
					}
				}
			}
		}

		//If E_F is empty then all clones (except for buried ones should be selected as mtp clones... and we don't need any further analysis...

		if (E_f.size() == 0) {
			cout << "E_f is empty...." << endl;
			for (set<Clone *, Contig::cloneCmp>::iterator cl_iter_1 =
					clones.begin(); cl_iter_1 != clones.end(); cl_iter_1++) {
				if ((*cl_iter_1)->isBuried()) {
					continue;
				}
				mtp_clones.push_back((*cl_iter_1)->getName());
			}
			//write them to the output file..
			outFile << con_iter->first << "\t";
			outFile << mtp_clones.size() << "\t";
			for (unsigned int i = 0; i < mtp_clones.size(); i++) {
				outFile << mtp_clones[i] << " ";
			}
			outFile << endl;
			//skip to the next contig...
			continue;
		}

		//generate clone-fragment graph for E_f.. 
		U_Graph g_f;

		NameUVertexMap name_2_vertex;
		U_VertexName v_name = get(vertex_name, g_f);
		U_EdgeWeight e_weight = get(edge_weight, g_f);
		for (map<pair<string, string>, int, pairCmp >::iterator iter =
				E_f.begin(); iter != E_f.end(); iter++) {
			//cout << iter->first.first << "-" << iter->first.second << endl;
			U_Vertex u, v;
			NameUVertexMap::iterator pos;
			bool inserted;
			tie(pos, inserted) = name_2_vertex.insert(std::make_pair(iter->first.first, U_Vertex()));
			if (inserted) {
				u = add_vertex(g_f);
				v_name[u] = iter->first.first;
				pos->second = u;
			} else {
				u = pos->second;
			}
			tie(pos, inserted) = name_2_vertex.insert(std::make_pair(iter->first.second, U_Vertex()));
			if (inserted) {
				v = add_vertex(g_f);
				v_name[v] = iter->first.second;
				pos->second = v;
			} else {
				v = pos->second;
			}
			//add an edge between u & v..
			U_Edge e;
			bool b;
			tie(e, b) = add_edge(u, v, g_f);
			if (!b) {
				cerr << "cannot add edge between " << v_name[u] << " and "
						<< v_name[v] << endl;
				exit(1);
			}
			string clone_1 = get_clone_name(v_name[u]);
			string clone_2 = get_clone_name(v_name[v]);
			pair<string, string> p;
			p.first = min_str(clone_1, clone_2);
			p.second = max_str(clone_1, clone_2);
			e_weight[e] = overlap_score_exponent[p];
		}

		//generate graphviz output...
		stringstream ssGraphvizOutputFileName;
		ssGraphvizOutputFileName << outputFile;
		

		//compute the connected components..
		std::vector<int> vertex_ID_2_component_ID(num_vertices(g_f), -1);//uses the old number of vertices...
		//map<int, int, less<int> > vertex_ID_2_component_ID;
		connected_components(g_f, make_iterator_property_map(
				vertex_ID_2_component_ID.begin(), get(vertex_index, g_f)));

		map<int, int, less<int> > component_ID_2_size;

		for (std::vector<int>::size_type i = 0; i
				!= vertex_ID_2_component_ID.size(); ++i) {

			if (component_ID_2_size.find(vertex_ID_2_component_ID[i])
					== component_ID_2_size.end()) {
				component_ID_2_size[vertex_ID_2_component_ID[i]] = 1;
			} else {
				component_ID_2_size[vertex_ID_2_component_ID[i]]++;
			}
		}

		if (DEBUG > 2) {
			cout << "Computing which components should be splitted in contig "
					<< con_iter->first << endl;
		}

		map<int, set<int>, less<int> > component_ID_2_vertices;
		get_vertices_per_component(vertex_ID_2_component_ID,
				component_ID_2_vertices);

		int maxID = component_ID_2_vertices.size() + 1;

		for (map<int, set<int>, less<int> >::iterator iter =
				component_ID_2_vertices.begin(); iter
				!= component_ID_2_vertices.end(); iter++) {
			map<int, int, less<int> > color;
			if (component_ID_2_size[iter->first] == 1) {
				continue;
			}

//			if ((sulston_score_threshold >= HIGH_THRESHOLD && is_loosely_connected_component(iter->second, g_f, v_name))
			if (exists_weak_edge(iter->second, v_name, overlap_score_exponent) || is_fragment_values_not_within_tolerance(iter->second, v_name, tolerance) 
					|| is_same_clone_together_in_component(iter->second, v_name)) {
				//if (is_fragment_values_not_within_tolerance(iter->second, v_name) || is_same_clone_together_in_component(iter->second, v_name)){
				if (DEBUG > 2) {
					cout << "splitting component.. ID: " << v_name[*(component_ID_2_vertices[iter->first].begin())] << endl;
				}
				double min_flow = split_component(iter->second, g_f, e_weight, color);
				set<int> component_1;
				set<int> component_2;
				
				for (map<int, int, less<int> >::iterator color_iter =
						color.begin(); color_iter != color.end(); color_iter++) {
					if (color_iter->second == 0) {
						component_1.insert(color_iter->first);
					} else if (color_iter->second == 4) {
						component_2.insert(color_iter->first);
					} else {
						cerr << "Unknown color code: " << color_iter->second
								<< endl;
						exit(0);
					}
				}
				//Compute average min_flow value.. First count # edges to be removed..
				int n_removed_edges = 0;
				for (set<int>::iterator set_iter_1 = component_1.begin(); set_iter_1 != component_1.end(); set_iter_1++) {
					for (set<int>::iterator set_iter_2 = component_2.begin(); set_iter_2 != component_2.end(); set_iter_2++) {
						bool b;
						U_Edge e;
						tie(e, b) = edge(*set_iter_1, *set_iter_2, g_f);
						if (b) {
							n_removed_edges++;
						}
					}
				}	
				//			cout << "n_removed_edges: " << n_removed_edges << endl;
				double average_min_flow = min_flow/(double)n_removed_edges;
				double min_flow_threshold = -1.5*log10(sulston_score_threshold);
				
				//			cout << "average_min_flow: " << average_min_flow << endl;
				//			cout << "min_flow_threshold: " << min_flow_threshold << endl;
				
				//check if any of the components consists of a single clone-fragment..
				//if so check its degree..
				int degree = 0;
				if (component_1.size() == 1 || component_2.size() == 1){
					degree = n_removed_edges;
					//cout  << "degree: " << degree << endl;
				}
				//double connectivity = (double)degree/double(component_1.size() + component_2.size() -1);
				//			if (connectivity != 0){
				//				cout <<"connectivity: " << connectivity << endl;
				//			}

				if ( average_min_flow <= min_flow_threshold && degree < 3){
					
					for (set<int>::iterator set_iter_1 = component_1.begin(); set_iter_1
							!= component_1.end(); set_iter_1++) {
						for (set<int>::iterator set_iter_2 = component_2.begin(); set_iter_2
								!= component_2.end(); set_iter_2++) {
							bool b;
							U_Edge e;
							tie(e, b) = edge(*set_iter_1, *set_iter_2, g_f);
							if (b) {
								//cout << "Removing edge between " << v_name[*set_iter_1] << " and " << v_name[*set_iter_2] << endl;
								remove_edge(*set_iter_1, *set_iter_2, g_f);
							}
						}
					}
				}
				
				if (DEBUG > 2) {
					cout << "Printing new components..." << endl;
					cout << "Component 1: " << endl;
					for (set<int>::iterator t_iter = component_1.begin(); t_iter
							!= component_1.end(); t_iter++) {
						cout << v_name[*t_iter] << endl;
					}
					cout << "Component 2: " << endl;
					for (set<int>::iterator t_iter = component_2.begin(); t_iter
							!= component_2.end(); t_iter++) {
						cout << v_name[*t_iter] << endl;
					}
				}

				//add the new components in case they need to be splitted again....
				component_ID_2_vertices[maxID] = component_1;
				maxID++;
				component_ID_2_vertices[maxID] = component_2;
				maxID++;
			}
		}

		if (DEBUG > 2) {
			cout
					<< "Computing which components should be splitted in contig...Done "
					<< con_iter->first << endl;
			cout << "Computing connected components again.." << endl;
		}

		//compute the connected components again...
		//initiliaze the data structure for connected components...
		vertex_ID_2_component_ID.clear();
		vertex_ID_2_component_ID.resize(num_vertices(g_f), -1);
		component_ID_2_size.clear();
		//compute connected components..
		connected_components(g_f, make_iterator_property_map(
				vertex_ID_2_component_ID.begin(), get(vertex_index, g_f)));
		//generate component_ID_2_size map...
		for (std::vector<int>::size_type i = 0; i
				!= vertex_ID_2_component_ID.size(); ++i) {
			if (component_ID_2_size.find(vertex_ID_2_component_ID[i])
					== component_ID_2_size.end()) {
				component_ID_2_size[vertex_ID_2_component_ID[i]] = 1;
			} else {
				component_ID_2_size[vertex_ID_2_component_ID[i]]++;
			}
		}

		//cout << "Connected components have been computed.." << endl;
		//cout << "Clone component size is being computed\n" << endl;
		map<int, int, less<int> > component_ID_2_clone_component_frequecy;
		compute_clone_component_size(vertex_ID_2_component_ID, v_name,
				component_ID_2_clone_component_frequecy);
		//print for DEBUGGING
		/*cout << "Component_ID_2_clone_component_frequency...\n";
		 for (map<int, int,less<int> >::iterator t_iter = component_ID_2_clone_component_frequecy.begin(); t_iter != component_ID_2_clone_component_frequecy.end(); t_iter++){
		 cout << t_iter->first << "\t" << t_iter->second << endl;
		 }

		 */

			
		
		//*********************************Get rid of spurious components*******************************...
		//TODO: Removal of spurious component is disabled for testing purposes..
		//BEGIN
		
		 
		graph_traits< U_Graph >::vertex_iterator vi, vi_end;
/*
		for (tie(vi, vi_end) = vertices(g_f); vi != vi_end; vi++) {
			int s = component_ID_2_size[vertex_ID_2_component_ID[*vi]];
			int nf = component_ID_2_clone_component_frequecy[vertex_ID_2_component_ID[*vi]];
			char chr_key[20];
			sprintf(chr_key, "%d-%d-%d", contig_size, nf, s);
			string key = chr_key;
		//	if (spurious_component_probability_table.find(key) != spurious_component_probability_table.end()){
//				if (spurious_component_probability_table[key] > p_value){
//					if ((s-1)*nf >= 3){
//						cout << "--------removing spurious component which would not removed by previous approach...!!: " << key << endl;
//						cout << "previous_value: " << (s - 1)*nf << endl;
//					}
//					clear_vertex(*vi, g_f);
//				}else{
//					if ((s-1)*nf < 3){
//						cout << "++++++++++++++keeping component which would be removed by previous approeach.." << endl;
//						cout << "previous_value: " << (s - 1)*nf << endl;
//					}
//				}
//			}else{
				if ((s-1)*nf < 3){
					clear_vertex(*vi, g_f);
//					cout << "C-NF-S: " << key << " does not exist in probability table...\n";
//					cout << "++++++++++++++keeping component which would be removed by previous approeach.." << endl;
//					cout << "previous_value: " << (s - 1)*nf << endl;
				}
//			}
			//		if ((component_ID_2_size[vertex_ID_2_component_ID[*vi]] - 1)*component_ID_2_clone_component_frequecy[vertex_ID_2_component_ID[*vi]] < p_value) {
			//			//cout << "vertex-cleared .." << v_name[*vi];
			//			clear_vertex(*vi, g_f);
			//		}
		}
		

		if (DEBUG > 2) {
			cout << "Spurious components eliminated for contig "
					<< con_iter->first << endl;
		}
		 //FOR DEBUGGING PURPOSES ONLY.. -- BEGIN
		 if (DEBUG > 3){
			 stringstream ssGraphvizOutputFileName_temp1;
			 ssGraphvizOutputFileName_temp1 << outputFile;
			 
			 string graphizOutputFileName_temp1 = ssGraphvizOutputFileName_temp1.str() + "-" + con_iter->first + "_spur_comp_rmvd.dot";
			 ofstream graphvizOutFile_temp1(graphizOutputFileName_temp1.c_str(), ios::out);
			 if (!graphvizOutFile_temp1) {
				 cerr << "cannot open output file \n";
				 exit(1);
			 }
		 
			 write_graphviz(graphvizOutFile_temp1, g_f, custom_label_writer<U_VertexName>(v_name), make_label_writer<U_EdgeWeight>(e_weight), graph_writer());
			 graphvizOutFile_temp1.close();
		 }
		 //FOR DEBUGGING PURPOSES ONLY.. -- END
		//*/
		//END: Disabling removal of spurious components...
		/**********************************************************************/

		
		
		
		
		//*******************************Computing Clones That are Buried into Multiple Clones*******************************...
		 
		//compute the connected components again...
		//initiliaze the data structure for connected components...
		vertex_ID_2_component_ID.clear();
		vertex_ID_2_component_ID.resize(num_vertices(g_f), -1);
		component_ID_2_size.clear();
		//compute connected components..
		connected_components(g_f, make_iterator_property_map(vertex_ID_2_component_ID.begin(), get(vertex_index, g_f)));
		//generate component_ID_2_size map...
		for (std::vector<int>::size_type i = 0; i != vertex_ID_2_component_ID.size(); ++i) {
			if (component_ID_2_size.find(vertex_ID_2_component_ID[i]) == component_ID_2_size.end()) {
				component_ID_2_size[vertex_ID_2_component_ID[i]] = 1;
			} else {
				component_ID_2_size[vertex_ID_2_component_ID[i]]++;
			}
		}
		
		vector<string> candidate_buried_clones;
		compute_clones_buried_into_multiple_clones(vertex_ID_2_component_ID, v_name, cloneObjectByName, candidate_buried_clones);
		
	//	if (candidate_buried_clones.size() > 0){
	//		cout << "Candidate buried clones: " << endl;
	//		for (unsigned int i = 0; i < candidate_buried_clones.size(); i++){
	//			cout << candidate_buried_clones[i] << endl;
	//		}
	//		cout << "----------" << endl;
	//	}

		//compute connected component ID set for candidate buried clones..
		map<string, set<int>, strCmp> cloneName_2_connected_component_IDs;
		for(unsigned int i = 0; i < vertex_ID_2_component_ID.size(); i++){
			string cloneName = get_clone_name(v_name[i]);
			if (cloneName_2_connected_component_IDs.find(cloneName) == cloneName_2_connected_component_IDs.end()){
				set<int> s;
				s.insert(vertex_ID_2_component_ID[i]);
				cloneName_2_connected_component_IDs[cloneName] = s;
			}else{
				cloneName_2_connected_component_IDs[cloneName].insert(vertex_ID_2_component_ID[i]);
			}
		}
		
		
		map<int, set<int>, less<int> > component_ID_to_vertices;
		get_vertices_per_component(vertex_ID_2_component_ID, component_ID_to_vertices);
		map<string, int, strCmp> remove_from_candidate_list;
		
		for (unsigned int i = 0; i < candidate_buried_clones.size(); i++){
			//set<string> canonical_clones;
			set<int> component_IDs = cloneName_2_connected_component_IDs[candidate_buried_clones[i]];
			//cout << "component_IDs.size() " << component_IDs.size() << endl;
			//int n_clusters = 0;
			map<string, int, strCmp> cloneName_2_clusterID;
			map<int, set<string> > clusterID_2_clones;
			int clusterID = 0;
			for (set<int>::iterator s_iter = component_IDs.begin(); s_iter != component_IDs.end(); s_iter++){
				set<int> vertices = component_ID_to_vertices[*s_iter];
				//	cout << "vertices.size(): " << vertices.size() << endl;
				
				int min_clusterID = clusterID;
				for(set<int>::iterator s_iter2 = vertices.begin(); s_iter2 != vertices.end(); s_iter2++){
					if (cloneName_2_clusterID.find(get_clone_name(v_name[*s_iter2])) != cloneName_2_clusterID.end()){
						if (cloneName_2_clusterID[get_clone_name(v_name[*s_iter2])] < min_clusterID){
							min_clusterID = cloneName_2_clusterID[get_clone_name(v_name[*s_iter2])];
						}
					}
				}
				
				for(set<int>::iterator s_iter2 = vertices.begin(); s_iter2 != vertices.end(); s_iter2++){
					if (cloneName_2_clusterID.find(get_clone_name(v_name[*s_iter2])) != cloneName_2_clusterID.end()){
						if (cloneName_2_clusterID[get_clone_name(v_name[*s_iter2])] != min_clusterID){
							//merge...
							int old_clusterID = cloneName_2_clusterID[get_clone_name(v_name[*s_iter2])];
							int new_clusterID = min_clusterID;
							cout << "merging.." << old_clusterID << " to " << new_clusterID << endl;
							
							set<string> clones_in_old_cluster = clusterID_2_clones[old_clusterID];
							for(set<string>::iterator s_iter3 = clones_in_old_cluster.begin(); s_iter3 != clones_in_old_cluster.end(); s_iter3++){
								cloneName_2_clusterID[*s_iter3] = new_clusterID;
							}
							clusterID_2_clones.erase(old_clusterID);
						}
					}else{
						cloneName_2_clusterID[get_clone_name(v_name[*s_iter2])] = min_clusterID;
						if (clusterID_2_clones.find(min_clusterID) != clusterID_2_clones.end()){
							clusterID_2_clones[min_clusterID].insert(get_clone_name(v_name[*s_iter2]));
						}else{
							set<string> s;
							s.insert(get_clone_name(v_name[*s_iter2]));
							clusterID_2_clones[min_clusterID] = s;
						}
					}
				}
				clusterID++;
			}
			
			//cout << "clusterID_2_clones.size() : " << clusterID_2_clones.size() << endl;
			
			if (clusterID_2_clones.size() > 1){
				cout << "remove " << candidate_buried_clones[i] << " from candidate_buried_clones_list..." << endl;
				remove_from_candidate_list[candidate_buried_clones[i]] = 1;
			}
		}
		
		//remove clones from candidate list..
		for (unsigned int i = 0; i < candidate_buried_clones.size(); i++){
			if (remove_from_candidate_list.find(candidate_buried_clones[i]) != remove_from_candidate_list.end()){
				cout << "removing " << candidate_buried_clones.at(i) << endl;
				candidate_buried_clones.erase(candidate_buried_clones.begin() + i);
			}
		}
		
		
		
		
		
		//For each connected component C_i that consists of only candidate buried clones, store the CLONES_C_i.. and the union of all stored CLONES_C_i.
		vector<set<string> > candidate_clone_set;
		set<string> clones_in_candidate_clone_set;//union of sets in the candidate_clone_set vector.. Some clones in candidate_buried_clones may not be in this set.. see below..	
		if (candidate_buried_clones.size() > 0){
			generateCandidateSet(vertex_ID_2_component_ID, v_name, candidate_buried_clones, candidate_clone_set, clones_in_candidate_clone_set);
		}

	//	//printing for debugging...
	//	ostream_iterator<string> output_temp( cout, " ");
	//	if (candidate_clone_set.size() > 0){
	//		cout << "Candidate clone set: " << endl;
	//		for (unsigned int i = 0; i < candidate_clone_set.size(); i++){
	//			copy( candidate_clone_set[i].begin(), candidate_clone_set[i].end(), output_temp);
	//			cout << "\n";
	//		}
	//	}
		
		
		//compute the minimum hitting set for the candidate_clone_set..candidate_clone_set is MODIFIED BY compute_minimum_hitting_set..
		set<string> minimum_hitting_set;
		if (candidate_clone_set.size() > 0){
			compute_minimum_hitting_set(candidate_clone_set, minimum_hitting_set);
		}
		
	//	if (minimum_hitting_set.size() > 0){
	//		copy(minimum_hitting_set.begin(), minimum_hitting_set.end(), output_temp);
	//		cout << "\n";
	//	}
		
		//all elements in candidate_buried_clones that are not in minimum_hitting_set vector can be buried...
		//Two types: a clone occurs in candidate_clone_set, but not in MHT (buried into candidate_buried_clones..)
		//a clone even does not occur in candidate_clone_set (buried into canonical clones)..
		for (unsigned int i = 0; i < candidate_buried_clones.size(); i++){
			//search in minimum_hitting_set_vector.. if not found, set it buried and remove from the MFG...
			if (minimum_hitting_set.find(candidate_buried_clones[i]) == minimum_hitting_set.end()){
	//			cout << "Removing " << candidate_buried_clones[i] << "..." << endl;
				cloneObjectByName[candidate_buried_clones[i]]->setPartiallyBuried(true);
				remove_clone_from_MFG(g_f, candidate_buried_clones[i], v_name);
			}
		}
		
		//*/
		//DEBUG-END
		//**********************************************************/		
		
	
		//check for violation of overlaps... If A-B, A-C are in the graph then B-C cannot be in the graph.. If B-C is also in the graph then one of these three components will be removed..
		//The pair that has the minimum sulston score will be eliminated...
		//FOR DEBUGGING PURPOSES ONLY----BEGIN


		//compute the connected components again...
		//initiliaze the data structure for connected components...
		vertex_ID_2_component_ID.clear();
		vertex_ID_2_component_ID.resize(num_vertices(g_f), -1);
		component_ID_2_size.clear();
		//compute connected components..
		connected_components(g_f, make_iterator_property_map(
				vertex_ID_2_component_ID.begin(), get(vertex_index, g_f)));
		//generate component_ID_2_size map...
		for (std::vector<int>::size_type i = 0; i
				!= vertex_ID_2_component_ID.size(); ++i) {
			if (component_ID_2_size.find(vertex_ID_2_component_ID[i])
					== component_ID_2_size.end()) {
				component_ID_2_size[vertex_ID_2_component_ID[i]] = 1;
			} else {
				component_ID_2_size[vertex_ID_2_component_ID[i]]++;
			}
		}

		//Check for # unmatched fragments for the clones in the graph..
		map<string, int, strCmp> clone_2_total_matched;//clone_name to total number of matched bands...(# times clone occurs in a component of size >= 2)
		for (std::vector<int>::size_type i = 0; i
				< vertex_ID_2_component_ID.size(); i++) {
			string clone_name = get_clone_name(v_name[i]);
			if (component_ID_2_size[vertex_ID_2_component_ID[i]] >= 2) {
				if (clone_2_total_matched.find(clone_name)
						== clone_2_total_matched.end()) {
					clone_2_total_matched[clone_name] = 1;
				} else {
					clone_2_total_matched[clone_name]++;
				}
			}
		}

		for (set<Clone *, Contig::cloneCmp>::iterator iter = clones.begin(); iter
				!= clones.end(); iter++) {
			if ((*iter)->isBuried()) {
				continue;
			}
			vector <Band *> fragments = (*iter)->getFragments();
			//int n_total_bands = fragments.size();
			int total_matched = 0;
			if (clone_2_total_matched.find((*iter)->getName()) != clone_2_total_matched.end()) {
				total_matched = clone_2_total_matched[(*iter)->getName()];
			}
			//int n_unmatched_bands = n_total_bands - total_matched;
			//if (DEBUG > 2) {
			//if ((double)n_unmatched_bands / (double)n_total_bands*100 >= EXTRA_FRAGMENT_PERCENTAGE) {
			if (total_matched == 0) {
				mtp_clones.push_back((*iter)->getName());
			}
		}

		//compute the connected components again...
		//initiliaze the data structure for connected components...
		vertex_ID_2_component_ID.clear();
		vertex_ID_2_component_ID.resize(num_vertices(g_f), -1);
		component_ID_2_size.clear();
		//compute connected components..
		connected_components(g_f, make_iterator_property_map(
				vertex_ID_2_component_ID.begin(), get(vertex_index, g_f)));
		//generate component_ID_2_size map...
		for (std::vector<int>::size_type i = 0; i
				!= vertex_ID_2_component_ID.size(); ++i) {
			if (component_ID_2_size.find(vertex_ID_2_component_ID[i])
					== component_ID_2_size.end()) {
				component_ID_2_size[vertex_ID_2_component_ID[i]] = 1;
			} else {
				component_ID_2_size[vertex_ID_2_component_ID[i]]++;
			}
		}

		//	cout << "Number of vertices in the clone-fragment graph: " << vertex_ID_2_component_ID.size() << endl; 

		if (vertex_ID_2_component_ID.size() == 0) {
			cerr << "Error: No node in clone-fragment graph.. " << endl;
			exit(1);
		}

		if (vertex_ID_2_component_ID.size() == 1) {
			cerr
					<< "Error: Only one vertex in the component.. Should not be..."
					<< endl;
			exit(1);
		}

		//component_ID_to_vertices computed to be used below.....
		component_ID_to_vertices.clear();
		get_vertices_per_component(vertex_ID_2_component_ID, component_ID_to_vertices);
//		cout << "component_ID_to_vertices.size(): " << component_ID_to_vertices.size() << endl;
		
		
		//int max_size = get_max_component_size(component_ID_2_size, vertex_ID_2_component_ID, v_name);
		//max_size = 2;
		//for (int i = 2; i == max_size; i++){
		//compute the adjacency list for each clone...
		map<string, set<string>, strCmp > connected_clones_map;//shows the overlapping clone(s) for each clone based on all components
		map<string, int, strCmp> clone_components;//stores names of underlying clones of each component (A-B occurs once, A-C occur once, etc)..
		get_clone_components(vertex_ID_2_component_ID, v_name, clone_components);
		compute_connected_clone_list(clone_components, connected_clones_map);

		
		
		if (DEBUG) {
			//print out clone_adjacency_map...before elimination...
			std::ostream_iterator<string> output_t(cout, " ");
			for (map<string, set<string>, strCmp >::iterator iter =
					connected_clones_map.begin(); iter
					!= connected_clones_map.end(); iter++) {
				vector<string> sorted_clones;
				sort_clones_based_on_overlap_score(iter->first, iter->second,
						sorted_clones, overlap_score_exponent);
				cout << iter->first << "\t" << iter->second.size() << "\t";
				std::copy(sorted_clones.begin(), sorted_clones.end(), output_t);
				cout << endl;
			}
		}
		//based on G_F compute how many times a pair of clones occur in the same component (i.e. share a fragment)
		map<pair<string, string>, int, pairCmp>
				n_shared_fragments_between_clones;
		get_n_shared_components_between_clones(vertex_ID_2_component_ID,
				v_name, n_shared_fragments_between_clones);
		map<string, set<string>, strCmp> adjacency_map;
		U_Graph g_temp;
		get_adjacency_map(connected_clones_map, overlap_score_exponent,
				n_shared_fragments_between_clones, adjacency_map,
				con_iter->first);

		//		if (DEBUG){
		//			
		//			cout << "Before Pruning" << endl;
		//			std::ostream_iterator<string> output(cout, " ");
		//			for (map<string, set<string>, strCmp >::iterator iter = connected_clones_map.begin(); iter != connected_clones_map.end(); iter++){
		//				vector<string> sorted_clones;
		//				sort_clones_based_on_overlap_score(iter->first, iter->second, sorted_clones, overlap_score_exponent);
		//				cout << iter->first << "\t" << iter->second.size() << "\t";
		//				std::copy(sorted_clones.begin(), sorted_clones.end(), output);
		//				cout << endl;
		//			}		
		//		}
		//select the best two and make sure there are at least 2 boundary clones (i.e. w/ one connection only..)...
		//prune_connected_clones_map(connected_clones_map, overlap_score_exponent);

		//print out clone_adjacency_map after elimination.....


		if (DEBUG) {
			std::ostream_iterator<string> output_after(cout, " ");
			//		bool is_tree = false;

			cout << "After Pruning" << endl;
			for (map<string, set<string>, strCmp >::iterator iter =
					adjacency_map.begin(); iter != adjacency_map.end(); iter++) {
				vector<string> sorted_clones;
				sort_clones_based_on_overlap_score(iter->first, iter->second,
						sorted_clones, overlap_score_exponent);
				cout << iter->first << "\t" << iter->second.size() << "\t";
				//			if (iter->second.size() > 2 && !is_tree){
				////				cout << "TREE" << endl;
				////				is_tree = true;
				////				n_tree++;
				//			}
				std::copy(sorted_clones.begin(), sorted_clones.end(),
						output_after);
				cout << endl;
			}
		}

		vector< vector<string> > ordered_clones;
		traverse_clones(adjacency_map, overlap_score_exponent, ordered_clones);

		if (outputOrderedClonesFilePresent) {
			//			cout << "order_2_clone.size(): " << ordered_clones.size() << endl;
			//			cout << "order: " << endl;
			for (unsigned int i = 0; i < ordered_clones.size(); i++) {
				vector <string> ordered_clones_temp(ordered_clones[i]);
				outputOrderedClonesFile << con_iter->first << "_" << i << "\t"
						<< ordered_clones_temp.size() << "\t";
				for (unsigned int j = 0; j < ordered_clones_temp.size(); j++) {
					outputOrderedClonesFile << ordered_clones_temp[j] << " ";
				}
				outputOrderedClonesFile << endl;
			}
		}

		map<string, string, strCmp > new_connected_clones_map;//for a given clone, this map shows the connected clone w/ max order..
		for (unsigned int index = 0; index < ordered_clones.size(); index++) {
			vector<string> ordered_clones_temp(ordered_clones[index]);

			for (unsigned int d = 1; d < ordered_clones_temp.size(); d++) {
				for (unsigned int i = 0; i < ordered_clones_temp.size() - d; i++) {
					//check if ordered_clones_i & ordered_clones_i+d are overlapping..
					//					if (d == 1 || (new_connected_clones_map.find(ordered_clones_temp[i+1]) != new_connected_clones_map.end() && new_connected_clones_map[ordered_clones_temp[i+1]].compare(ordered_clones_temp[i+d]) == 0 )
					//							&& new_connected_clones_map.find(ordered_clones_temp[i]) != new_connected_clones_map.end() && new_connected_clones_map[ordered_clones_temp[i]].compare(ordered_clones_temp[i+d-1]) == 0){
					bool is_overlapping = false;
					pair<string, string> p;
					p.first = min_str(ordered_clones_temp[i],
							ordered_clones_temp[i+d]);
					p.second = max_str(ordered_clones_temp[i],
							ordered_clones_temp[i+d]);

					bool cr_1 = (d == 1) || ((new_connected_clones_map.find(ordered_clones_temp[i]) != new_connected_clones_map.end() 
								&& new_connected_clones_map[ordered_clones_temp[i]].compare(ordered_clones_temp[i+d-1]) == 0)
							&& (new_connected_clones_map.find(ordered_clones_temp[i+1]) != new_connected_clones_map.end() 
									&& new_connected_clones_map[ordered_clones_temp[i+1]].compare(ordered_clones_temp[i+d]) == 0 ));
//					bool cr_1 =	(d == 1) || ((new_connected_clones_map.find(ordered_clones_temp[i])	!= new_connected_clones_map.end() && new_connected_clones_map[ordered_clones_temp[i]].compare(ordered_clones_temp[i+d-1]) == 0));
					
					set<string> connected_vertices;
					for (unsigned int j = i; j  <= i+d; j++){
						connected_vertices.insert(ordered_clones_temp[j]);
					}
//					cout << "connected_vertices..." << endl;
//					for (set<string>::iterator iter = connected_vertices.begin(); iter != connected_vertices.end(); iter++){
//						cout << *iter << " ";
//					}
//					cout  << endl;
//					cout << "-----" << endl;

					
					bool cr_2 = connected_clones_map[ordered_clones_temp[i]].find(ordered_clones_temp[i+d]) != connected_clones_map[ordered_clones_temp[i]].end();
					
					//bool cr_2 = is_together(component_ID_to_vertices, v_name, connected_vertices);
//					cout << "is_together: " << cr_2 << endl;
					//TODO: 1/2 added for testing purposes..
					bool cr_3 = overlap_score_exponent[p] > -0.4*log10(sulston_score_threshold);
					
					//if clone_i - clone_i+d-1 and clone_i+d-1 - clone_i+d are not overlapping well then clone_i - clone_i+d does not overlap...
					bool cr_4;
					if (d == 1){
						cr_4 = true;
					}else{
						pair<string, string> p2_1;
						p2_1.first = min_str(ordered_clones_temp[i], ordered_clones_temp[i+d-1]);
						p2_1.second = max_str(ordered_clones_temp[i], ordered_clones_temp[i+d-1]);
						
						pair<string, string> p2_2;
						p2_2.first = min_str(ordered_clones_temp[i+d-1], ordered_clones_temp[i+d]);
						p2_2.second = max_str(ordered_clones_temp[i+d-1], ordered_clones_temp[i+d]);
						cr_4 = overlap_score_exponent[p2_1] + overlap_score_exponent[p2_2] > 3*MTP_MST_VERY_HIGH_THRESHOLD;
					}

//					if (!cr_4){
//						cout << cr_1 << "-" << cr_2 << "-" << cr_3 << endl;
//						if (cr_1 && cr_2 && cr_3){
//							cerr << "cr_4 false, others true...!" << endl;
//							cerr << ordered_clones_temp[i] << "-" << ordered_clones_temp[i+d] << endl;
//							//exit(1);
//						}
//					}
					
					bool r = false;
//					if (criteria == 1) {
//					int c1 = cr_1;
//					int c2 = cr_2;
//					int c3 = cr_3;
//					int total = c1 + c2 + c3;
										
					//r = total >= 2;
					//cout << "overlap...\t" << ordered_clones_temp[i] << "\t" << ordered_clones_temp[i+d] << "\t" << cr_1 << "\t" << cr_2 << "\t" << cr_3 << "\t" << cr_4 << endl;
					
				//					} else if (criteria == 2) {
					//r = cr_1 && cr_2;// && cr_3;// && cr_4;
					//TODO r was initially set to cr_1 && cr_2 && cr_3;
					r = cr_1 && cr_2 && cr_3;
					r = cr_1 && cr_2 && cr_4;
					 
//					}else{
//						cerr << "Wrong code.. " << criteria << endl;
//						exit(1);
//					}
					
//					else if (criteria == 3) {
//						r = cr_1 && cr_2;
//					} else if (criteria == 4) {
//						r = cr_2 && cr_3;
//					} else if (criteria == 5) {
//						r = cr_3;
//					} else if (criteria == 6) {
//						r = cr_2;
//					}

					if (r) {
						new_connected_clones_map[ordered_clones_temp[i]]
								= ordered_clones_temp[i+d];
						is_overlapping = true;
					}

					//					if (d == 1 || (new_connected_clones_map.find(ordered_clones_temp[i]) != new_connected_clones_map.end() && new_connected_clones_map[ordered_clones_temp[i]].compare(ordered_clones_temp[i+d-1]) == 0)){
					//					
					//						pair<string, string> p;
					//						p.first = min_str(ordered_clones_temp[i], ordered_clones_temp[i+d]);
					//						p.second = max_str(ordered_clones_temp[i], ordered_clones_temp[i+d]);
					//						if (connected_clones_map[ordered_clones_temp[i]].find(ordered_clones_temp[i+d]) != connected_clones_map[ordered_clones_temp[i]].end() && overlap_score_exponent[p] > -log10(sulston_score_threshold)){
					////						if (connected_clones_map[ordered_clones_temp[i]].find(ordered_clones_temp[i+d]) != connected_clones_map[ordered_clones_temp[i]].end()){
					////						if (overlap_score_exponent[p] > -log10(sulston_score_threshold)){
					//							new_connected_clones_map[ordered_clones_temp[i]] = ordered_clones_temp[i+d];
					//							is_overlapping = true;
					//	//						cout << ordered_clones[i] << "-" << ordered_clones[i+d] << endl;
					//						}
					//					}

					if (output_overlapping_clones_file_present){
						if (is_overlapping) {
							outputOverlappingClonesFile << ordered_clones_temp[i] << "\t"
									<< ordered_clones_temp[i+d] << endl;
						} else {
							//print the non-overlapping_pair...
							outputNonOverlappingClonesFile << ordered_clones_temp[i] << "\t"
									<< ordered_clones_temp[i+d] << endl;
						}
					}
				}
			}
		}
		
		string current_clone = "";
		outFile << con_iter->first << "\t";
		//		/* !!!!!!!!!!!!!!   JUST FOR TESTING  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
		//		if (ordered_clones.size() >= 3 && ordered_clones.size() <= 3) {
		//			mtp_clones.push_back(ordered_clones[1]);
		//		}else{
		for (unsigned int i = 0; i < ordered_clones.size(); i++) {
			vector<string> ordered_clones_temp(ordered_clones[i]);

			int order = 0;
			int limit = ordered_clones_temp.size();
			//		if (ordered_clones.size() == 2) {limit++;}
			//		cout << "limit: " << limit << endl;
			while (order < limit) {
				//			cout << "hello..." << endl;
				//order = iter->first;
				current_clone = ordered_clones_temp[order];

				//			cout << current_clone << "-";
				//outFile << current_clone << " ";
				mtp_clones.push_back(current_clone);
				//go to furthest connected clone to current_clone..
				int new_order = order;
				//string max_clone = "";
				if (new_connected_clones_map.find(current_clone)
						!= new_connected_clones_map.end()) {
					string next_clone = new_connected_clones_map[current_clone];
					for (unsigned int i = 0; i < ordered_clones_temp.size(); i++) {
						if (ordered_clones_temp[i].compare(next_clone) == 0) {
							new_order = i;
							break;
						}
					}
				}

				
				if (new_order == order) {
					order++;
				} else {
					order = new_order;
				}
			}
		}
		outFile << mtp_clones.size() << "\t";
		for (unsigned int i = 0; i < mtp_clones.size(); i++) {
			outFile << mtp_clones[i] << " ";
		}

		outFile << endl;

		
		E.clear();
		E_f.clear();
	}
	
	//cout << "Number of trees: " << n_tree << endl;
	//ctgMTPFile.close();

	return 0;

}

int get_copy_id(string cloneName, int size) {
	pair<string, int> p(cloneName, size);
	if (fragment_2_copy_id.find(p) != fragment_2_copy_id.end()) {
		fragment_2_copy_id[p]++;
	} else {
		fragment_2_copy_id[p] = 1;
	}
	return fragment_2_copy_id[p];
}

string get_clone_name(string fragment_name) {
	char *tok;
	char fragment[20];
	strcpy(fragment, fragment_name.c_str());
	tok = strtok(fragment, "-");
	char clone_name[20];
	strcpy(clone_name, tok);
	string result = clone_name;
	return result;
}

void compute_redundant_components(const vector<int> &vertex_2_component_ID,
		const U_VertexName &v_name, set<int> &redundant_component_ID) {
	map<int, set<string>, less<int> > component_2_vertex_names;
	map<int, set<int>, less<int > > component_2_vertex_IDs;
	get_vertices_per_component(vertex_2_component_ID, component_2_vertex_IDs);

	for (map<int, set<int>, less<int> >::iterator iter =
			component_2_vertex_IDs.begin(); iter
			!= component_2_vertex_IDs.end(); iter++) {
		set <string> t;
		for (set<int>::iterator s_iter = iter->second.begin(); s_iter
				!= iter->second.end(); s_iter++) {
			t.insert(get_clone_name(v_name[*s_iter]));
		}
		component_2_vertex_names[iter->first] = t;
	}

	map<int, set<string>, less<int> >::iterator iter_1, iter_2;
	for (iter_1 = component_2_vertex_names.begin(); iter_1
			!= component_2_vertex_names.end(); iter_1++) {
		if (component_2_vertex_IDs[iter_1->first].size() == 1) {
			continue;
		}//ignore components of size 1..
		for (iter_2 = component_2_vertex_names.upper_bound(iter_1->first); iter_2
				!= component_2_vertex_names.end(); iter_2++) {
			if (component_2_vertex_IDs[iter_2->first].size() == 1) {
				continue;
			}//ignore components of size 1...
			//compute the set_difference betweeen iter_1->second and iter_2->second..
			list<string> diff;
			if (iter_1->second.size() <= iter_2->second.size()) {
				set_difference(iter_1->second.begin(), iter_1->second.end(),
						iter_2->second.begin(), iter_2->second.end(),
						front_inserter(diff));
				if (!diff.size()) {
					redundant_component_ID.insert(iter_2->first);
				}
			} else {
				set_difference(iter_2->second.begin(), iter_2->second.end(),
						iter_1->second.begin(), iter_1->second.end(),
						front_inserter(diff));
				if (!diff.size()) {
					redundant_component_ID.insert(iter_1->first);
				}
			}
		}
	}
}

void get_vertices_per_component(const vector<int>& vertex_ID_2_component_ID,
		map<int, set<int>, less<int> > & component_ID_2_vertices) {
	for (std::vector<int>::size_type i = 0; i
			!= vertex_ID_2_component_ID.size(); ++i) {
		int component_ID = vertex_ID_2_component_ID[i];
		if (component_ID_2_vertices.find(component_ID)
				== component_ID_2_vertices.end()) {
			set<int> s;
			s.insert(i);
			component_ID_2_vertices[component_ID] = s;
		} else {
			component_ID_2_vertices[component_ID].insert(i);
		}
	}
}

void compute_clone_component_size(const vector<int>& vertex_ID_2_component_ID,
		const U_VertexName & v_name,
		map<int, int, less<int> > & component_ID_2_clone_component_frequency) {
	//construct set of vertex ID's for each component...
	map<int, set<int>, less<int> > component_ID_to_vertices;
	get_vertices_per_component(vertex_ID_2_component_ID,
			component_ID_to_vertices);

	//construct componentID to clone-string..
	map<int, string, less<int> > component_ID_to_clone_string;
	for (map<int, set<int>, less<int> >::iterator iter =
			component_ID_to_vertices.begin(); iter
			!= component_ID_to_vertices.end(); iter++) {
		set<string> clone_names;
		for (set<int>::iterator s_iter = iter->second.begin(); s_iter
				!= iter->second.end(); s_iter++) {
			clone_names.insert(get_clone_name(v_name[*s_iter]));
		}
		//concatenate clone names...
		string concatenated_clone_name = "";
		concatenate_clones(clone_names, concatenated_clone_name);
		component_ID_to_clone_string[iter->first] = concatenated_clone_name;
	}

	//construct the clone_string to frequency..
	map<string, int, strCmp> clone_string_to_freq;
	for (map<int, string, less<int> >::iterator iter =
			component_ID_to_clone_string.begin(); iter
			!= component_ID_to_clone_string.end(); iter++) {
		if (clone_string_to_freq.find(iter->second)
				== clone_string_to_freq.end()) {
			clone_string_to_freq[iter->second] = 1;
		} else {
			clone_string_to_freq[iter->second]++;
		}
	}

	//construct component_ID_2_clone_component_frequency..For special components, assign a large number...
	for (map<int, string, less<int> >::iterator iter =
			component_ID_to_clone_string.begin(); iter
			!= component_ID_to_clone_string.end(); iter++) {
		if (is_special_component(component_ID_to_vertices[iter->first], v_name)) {
			component_ID_2_clone_component_frequency[iter->first] = BIG_NUMBER;
		} else {
			component_ID_2_clone_component_frequency[iter->first]
					= clone_string_to_freq[iter->second];
		}
	}
}

//for the clones that we must select, we generate a component that comprises only that clones fragments.. in this component, we don't care if the fragment values are w/in the tolerance, etc.
bool is_special_component(const set<int> & vertices_in_component,
		const U_VertexName & v_name) {
	if (vertices_in_component.size() == 1) {
		return false;
	}

	set<string> clone_names;
	for (set<int>::iterator iter = vertices_in_component.begin(); iter
			!= vertices_in_component.end(); iter++) {
		clone_names.insert(get_clone_name(v_name[*iter]));
	}
	return (clone_names.size() == 1); //if there are multiple vertices from only one clone in the component then it is special..
}

//bool is_loosely_connected_component(const set<int> & vertices_in_component,
//		const U_Graph & g, const U_VertexName & v_name) {
//	int number_of_edges_in_component = 0, ideal_number_of_edges_in_component =
//			0;
//	//first check if this is a special component..
//	if (is_special_component(vertices_in_component, v_name)) {
//		return false;
//	}
//
//	//if not a special component...
//	//compute the ideal number of edges in the component.. it is n.(n-1)/2 where n is the number of vertices in the component..
//	ideal_number_of_edges_in_component = vertices_in_component.size()
//			*(vertices_in_component.size() - 1)/2;
//
//	//compute the number of edges in the component.. it is the total degree of vertices in the component/2...
//	for (set<int>::iterator iter = vertices_in_component.begin(); iter
//			!= vertices_in_component.end(); iter++) {
//		//cout << v_name[*iter] << " degree : " << degree(*iter, g) << endl;
//		number_of_edges_in_component += degree(*iter, g);
//	}
//	number_of_edges_in_component /= 2;
//
//	if ((double)number_of_edges_in_component
//			/(double)ideal_number_of_edges_in_component*100 <= MIN_CONNECTIVITY) {
//		return true;
//	} else {
//		return false;
//	}
//}

//checks the sulston score between all pairs of clones in a given component.. if one of them has a sulston score lower than VERY_HIGH_THRESHOLD, then 
//return true.. If this function returns true for a component then this component has to be splitted...
bool exists_weak_edge(const set<int> &vertices_in_component, const U_VertexName & v_name, map<pair<string, string>, double, pairCmp > & overlap_score_exponent){
	set <string> clones;
	
	
	for (set<int>::iterator s_iter_1 = vertices_in_component.begin(); s_iter_1 != vertices_in_component.end(); s_iter_1++){
		clones.insert(get_clone_name(v_name[*s_iter_1]));
	}
	
	for (set<string>::iterator s_iter_1 = clones.begin(); s_iter_1 != clones.end(); s_iter_1++){
		for (set<string>::iterator s_iter_2 = clones.upper_bound(*s_iter_1); s_iter_2 != clones.end(); s_iter_2++){
			//check the overlap score of *s_iter_1 & *s_iter_2
			pair<string, string> p;
			p.first = min_str(*s_iter_1, *s_iter_2);
			p.second = max_str(*s_iter_1, *s_iter_2);
			
			if (overlap_score_exponent.find(p) == overlap_score_exponent.end()){
				cerr << "Cannot find the overlap score of " << p.first << " - " << p.second << endl;
				exit(1);
			}
			
			//print the overlap score for debuggin purposes
//			cout << "overlap_score_exponent[" << p.first << ", " << p.second << "]: " << overlap_score_exponent[p] << endl; 
			
			if (overlap_score_exponent[p] <= MTP_MST_VERY_HIGH_THRESHOLD){
				return true;
			}
		}
	}
	return false;
}

bool is_same_clone_together_in_component(
		const set<int> & vertices_in_component, const U_VertexName & v_name) {
	//first check if this is a special component..
	if (is_special_component(vertices_in_component, v_name)) {
		return false;
	}

	//if not a special component...
	set<string> clone_names;
	for (set<int>::iterator iter = vertices_in_component.begin(); iter
			!= vertices_in_component.end(); iter++) {
		clone_names.insert(get_clone_name(v_name[*iter]));
	}
	return (clone_names.size() != vertices_in_component.size()); //if these two statements are not equal, then at least one clone occurs more than once in this component  
}

/**
 * computes the difference between the second max and second min.. if it's more than >TOL then return true...)
 */
bool is_fragment_values_not_within_tolerance(
		const set<int> & vertices_in_component, const U_VertexName & v_name, const int tolerance) {
	
	
	//first check if this is a special component..
	if (is_special_component(vertices_in_component, v_name)) {
		return false;
	}

	//if not a special component...
	int min = BIG_NUMBER, max = 0, second_min = BIG_NUMBER, second_max = 0;;
	for (set<int>::iterator iter = vertices_in_component.begin(); iter
			!= vertices_in_component.end(); iter++) {
		//cout << get_fragment_size(v_name[*iter]) << endl;
		if (get_fragment_size(v_name[*iter]) < min) {
			second_min = min;
			min = get_fragment_size(v_name[*iter]);
		}else if (get_fragment_size(v_name[*iter]) < second_min){
			second_min = get_fragment_size(v_name[*iter]);
		}
		
		if (get_fragment_size(v_name[*iter]) > max) {
			second_max = max;
			max = get_fragment_size(v_name[*iter]);
		}else if (get_fragment_size(v_name[*iter]) > second_max) {
			second_max = get_fragment_size(v_name[*iter]);
		}
	}
	
	if (second_min == BIG_NUMBER){
		second_min = min;
	}
	
	if (second_max == 0){
		second_max = max;
	}
	
//	if (((second_max - second_min) > TOL)){
//		cout << "Min: " << min << "\n" << "Second_min: " << second_min << "\n" << "Max: " << max << "\n" << "Second_max: " << second_max << endl;
//	}
	
	return ((second_max - second_min) > tolerance);//if the difference between max and min is greater than TOL then return true...
}

int get_fragment_size(string fragment_name) {
	char f_name[20];
	strcpy(f_name, fragment_name.c_str());
	char* tok = strtok(f_name, "-");
	tok = strtok(NULL, "-");
	return atoi(tok);
}

double split_component(const set<int> & vertices_in_component, const U_Graph & g,
		const U_EdgeWeight & e_weight, map<int, int, less<int> > & min_cut) {
	//build a directed graph G'=(V',E') where V'=vertices_in_component and for each edge (u,v) in G where u,v are in vertices_in_component, add (u,v) and (v,u) with same capacity...
	map<int, int, less<int> > old_to_new_vertex_map;//maps vertices in V to V'...
	map<int, int, less<int> > new_to_old_vertex_map;//maps vertices in V' to V...
	Graph g_new;//G'
	property_map < Graph, edge_capacity_t >::type capacity = get(edge_capacity,
			g_new);
	property_map < Graph, edge_reverse_t >::type rev = get(edge_reverse, g_new);
	property_map < Graph, edge_residual_capacity_t >::type residual_capacity =
			get(edge_residual_capacity, g_new);

	vector<Vertex> g_new_vertices;

	for (set<int>::iterator iter_1 = vertices_in_component.begin(); iter_1
			!= vertices_in_component.end(); iter_1++) {
		//add an vertex to G'
		Vertex u;
		if (old_to_new_vertex_map.find(*iter_1) == old_to_new_vertex_map.end()) {
			old_to_new_vertex_map[*iter_1] = add_vertex(g_new);
			g_new_vertices.push_back(old_to_new_vertex_map[*iter_1]);
		}
		u = old_to_new_vertex_map[*iter_1];
		new_to_old_vertex_map[u] = *iter_1;
		for (set<int>::iterator iter_2 =
				vertices_in_component.upper_bound(*iter_1); iter_2
				!= vertices_in_component.end(); iter_2++) {
			Vertex v;
			if (old_to_new_vertex_map.find(*iter_2)
					== old_to_new_vertex_map.end()) {
				old_to_new_vertex_map[*iter_2] = add_vertex(g_new);
				g_new_vertices.push_back(old_to_new_vertex_map[*iter_2]);
			}
			v = old_to_new_vertex_map[*iter_2];

			//add an edge between u & v if their corresponding vertices in G are adjacent....
			U_Edge e;
			bool is_adjacent;
			tie(e, is_adjacent) = edge(*iter_1, *iter_2, g);
			if (is_adjacent) {
				bool b1, b2;
				Edge e1, e2, e1_rev, e2_rev;
				tie(e1, b1) = add_edge(u, v, g_new);
				tie(e1_rev, b2) = add_edge(v, u, g_new);
				if (!b1 || !b2) {
					cerr << "Cannot add edge in g_new..\n";
					exit(0);
				}

				capacity[e1] = e_weight[e];
				//cout << "capacity[e1]: " << capacity[e1] << endl;
				capacity[e1_rev] = 0;
				rev[e1] = e1_rev;
				rev[e1_rev] = e1;

				tie(e2, b1) = add_edge(v, u, g_new);
				tie(e2_rev, b2) = add_edge(u, v, g_new);
				if (!b1 || !b2) {
					cerr << "Cannot add edge in g_new..\n";
					exit(0);
				}
				capacity[e2] = e_weight[e];
				capacity[e2_rev] = 0;
				rev[e2] = e2_rev;
				rev[e2_rev] = e2;
			}
		}
	}

	//print_graph(g_new);
	//Now G' is built... A source vertex will be selected and edmunds-karp max flow algorithm will be run n-1 times..

	double min_flow = 10000;
	//select a source vertex
	Vertex s = g_new_vertices[0];
	for (unsigned int j = 1; j < g_new_vertices.size(); j++) {
		//run edmunds-karp between s & j.. (s: source, g_new_vertices[j]: sink..)
		vector<default_color_type> color(num_vertices(g_new));
		vector<Edge> pred(num_vertices(g_new));
		//cout << "Computing the max flow.." << endl;
		double flow = edmunds_karp_max_flow(g_new, s, g_new_vertices[j],
				capacity, residual_capacity, rev, &color[0], &pred[0]);
		//cout << "flow is " << flow << endl; 
		if (flow < min_flow) {
			min_flow = flow;
			min_cut.clear();
			for (unsigned int i = 0; i < color.size(); i++) {
				min_cut[new_to_old_vertex_map[i]] = color[i];
			}
		}
	}
	return min_flow;
}

int get_max_component_size(map<int, int, less<int> > & component_ID_2_size,
		const vector<int> &vertex_ID_2_component_ID, const U_VertexName & v_name) {
	int max = 0;
	map<int, set<int>, less<int> > component_ID_to_vertices;
	get_vertices_per_component(vertex_ID_2_component_ID,
			component_ID_to_vertices);
	for (map<int, int, less<int> >::iterator iter = component_ID_2_size.begin(); iter
			!= component_ID_2_size.end(); iter++) {
		if (!is_special_component(component_ID_to_vertices[iter->first], v_name)
				&& iter->second >= max) {
			max = iter->second;
		}
	}
	return max;
}

void get_clone_components(const vector<int> & vertex_ID_2_component_ID,
		const U_VertexName & v_name, map<string, int, strCmp> &clone_components) {
	//construct set of vertex ID's for each component...
	map<int, set<int>, less<int> > component_ID_to_vertices;
	get_vertices_per_component(vertex_ID_2_component_ID,
			component_ID_to_vertices);

	for (map<int, set<int>, less<int> >::iterator iter =
			component_ID_to_vertices.begin(); iter
			!= component_ID_to_vertices.end(); iter++) {
		set<string> clone_names;
		for (set<int>::iterator s_iter = iter->second.begin(); s_iter
				!= iter->second.end(); s_iter++) {
			clone_names.insert(get_clone_name(v_name[*s_iter]));
		}
		//concatenate clone names...
		string concatenated_clone_name = "";
		concatenate_clones(clone_names, concatenated_clone_name);
		clone_components[concatenated_clone_name] = 1;
	}

}

void concatenate_clones(const set<string> & clone_names,
		string & concatenated_clone_name, const string& delimiters) {
	concatenated_clone_name = "";
	for (set<string>::iterator s_iter = clone_names.begin(); s_iter
			!= clone_names.end(); s_iter++) {
		concatenated_clone_name += *s_iter + delimiters;
	}
}

void tokenize(const string& str, vector<string>& tokens,
		const string& delimiters) {
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos) {
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

void compute_connected_clone_list(map<string, int, strCmp> &clone_components,
		map<string, set<string>, strCmp> &connected_clones_map) {

	for (map<string, int, strCmp>::iterator iter = clone_components.begin(); iter
			!= clone_components.end(); iter++) {
		vector<string> clone_names;
		tokenize(iter->first, clone_names);
		//if (clone_names.size() == size){//consider components of a given size only...
		for (unsigned int i = 0; i < clone_names.size(); i++) {
			if (connected_clones_map.find(clone_names[i])
					== connected_clones_map.end()) {
				set<string> s;
				connected_clones_map[clone_names[i]] = s;
			}
			for (unsigned int j = 0; j < clone_names.size(); j++) {
				if (i != j) {
					connected_clones_map[clone_names[i]].insert(clone_names[j]);
				}
			}
		}
		//}
	}
}

/*void check_circular_violation(map<string, set<string>, strCmp > & clone_adjacency_map){
 for (map<string, set<string>, strCmp >::iterator iter = clone_adjacency_map.begin(); iter != clone_adjacency_map.end(); iter++){
 if (iter->second.size() > 2){
 cerr << iter->first << " is adjacent to " << iter->second.size() << " clones..." << endl;
 }else{
 
 }
 
 }
 }*/

void sort_clones_based_on_overlap_score(const string &source_clone,
		const set<string> &connected_clones, vector<string> &sorted_clones,
		map<pair<string, string>, double, pairCmp > & overlap_score_exponent) {
	map<double, list<string>, greater<double> > score_2_clone_name;
	for (set<string>::iterator iter = connected_clones.begin(); iter
			!= connected_clones.end(); iter++) {
		pair<string, string> p;
		p.first = min_str(source_clone, *iter);
		p.second = max_str(source_clone, *iter);
		if (score_2_clone_name.find(overlap_score_exponent[p])
				== score_2_clone_name.end()) {
			list<string> l;
			l.push_back(*iter);
			score_2_clone_name[overlap_score_exponent[p]] = l;
		} else {
			score_2_clone_name[overlap_score_exponent[p]].push_back(*iter);
		}
	}

	//cout << "Printing score_2_clone_name Map...\n";
	for (map<double, list<string>, greater<double> >::iterator iter =
			score_2_clone_name.begin(); iter != score_2_clone_name.end(); iter++) {
		if (DEBUG) {
			cout << iter->first << "\t" << iter->second.size() << endl;
		}
		for (list<string>::iterator l_iter = iter->second.begin(); l_iter
				!= iter->second.end(); l_iter++) {
			sorted_clones.push_back(*l_iter);
		}
	}
}

//based on G_F compute how many times a pair of clones occur in the same component (i.e. share a fragment)
void get_n_shared_components_between_clones(
		const vector<int> & vertex_ID_2_component_ID,
		const U_VertexName & v_name,
		map<pair<string, string>, int, pairCmp> &n_shared_components_between_clones) {
	map<int, set<int>, less<int> > component_ID_2_vertices;
	get_vertices_per_component(vertex_ID_2_component_ID,
			component_ID_2_vertices);
	for (map<int, set<int>, less<int> >::iterator iter =
			component_ID_2_vertices.begin(); iter
			!= component_ID_2_vertices.end(); iter++) {
		for (set<int>::iterator s_iter = iter->second.begin(); s_iter
				!= iter->second.end(); s_iter++) {
			pair<string, string> p;
			string clone_1 = get_clone_name(v_name[iter->first]);
			string clone_2 = get_clone_name(v_name[*s_iter]);
			//int f_size = get_fragment_size(v_name[iter->first]);
			p.first = min_str(clone_1, clone_2);
			p.second = max_str(clone_1, clone_2);
			if (n_shared_components_between_clones.find(p)
					== n_shared_components_between_clones.end()) {
				//n_shared_components_between_clones[p] = f_size;
				n_shared_components_between_clones[p] = 1;
			} else {
				//n_shared_components_between_clones[p] += f_size;
				n_shared_components_between_clones[p] ++;
			}
		}
	}
}

void get_adjacency_map(
		map<string, set<string>, strCmp > & connected_clones_map,
		map<pair<string, string>, double, pairCmp > & overlap_score_exponent,
		map<pair<string, string>, int, pairCmp> &n_shared_components_between_clones,
		map<string, set<string>, strCmp> & adjacency_map, const string contig_id) {

	//build an overlap graph.. 
	U_Graph g;

	NameUVertexMap name_2_vertex;
	U_VertexName v_name = get(vertex_name, g);
	U_EdgeWeight e_weight = get(edge_weight, g);

	//add a vertex for each clone first...
	for (map<string, set<string>, strCmp>::iterator iter =
			connected_clones_map.begin(); iter != connected_clones_map.end(); iter++) {
		U_Vertex u;
		NameUVertexMap::iterator pos;
		bool inserted;
		tie(pos, inserted) = name_2_vertex.insert(std::make_pair(iter->first, U_Vertex()));
		if (inserted) {
			u = add_vertex(g);
			v_name[u] = iter->first;
			pos->second = u;
		} else {
			u = pos->second;
		}
	}
	//add edges between vertices...
	for (map<string, set<string>, strCmp >::iterator iter =
			connected_clones_map.begin(); iter != connected_clones_map.end(); iter++) {
		U_Vertex u = name_2_vertex[iter->first];
		for (set<string>::iterator s_iter = iter->second.begin(); s_iter
				!= iter->second.end(); s_iter++) {
			U_Vertex v = name_2_vertex[*s_iter];
			bool b;
			U_Edge e;
			tie(e, b) = edge(u, v, g);
			//if there is no edge between u  & v then add one...
			if (!b) {
				tie(e, b) = add_edge(u, v, g);
				if (!b) {
					cerr << "cannot add edge between " << v_name[u] << " and "
							<< v_name[v] << endl;
					exit(1);
				}
				pair<string, string> p;
				p.first = min_str(v_name[u], v_name[v]);
				p.second = max_str(v_name[u], v_name[v]);
				//HERE OVERLAP SCORE IS THE LOG OF THE SULSTON SCORE.. 
				//IT'S ALSO POSSIBLE TO USE THE NUMBER (OR TOTAL LENGTH) OF THE SHARED FRAGMENTS AS THE WEIGHT..
				e_weight[e] = -1 * overlap_score_exponent[p];
				//e_weight[e] = -1 * n_shared_components_between_clones[p];
			}
		}
	}

	//Graph g is ready....

	//print the graph g before computing mst...
//	string graphizOutputFileName ="output_small-" + contig_id + ".dot";
//	ofstream graphvizOutFile(graphizOutputFileName.c_str(), ios::out);
//	if (!graphvizOutFile) {
//		cerr << "cannot open output file \n";
//		exit(1);
//	}
//	write_graphviz(graphvizOutFile, g,
//			custom_label_writer<U_VertexName>(v_name), make_label_writer<
//					U_EdgeWeight>(e_weight), graph_writer());
//	graphvizOutFile.close();

	//compute MST of the graph...
	std::vector < U_Edge > spanning_tree;

	kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

	//remove the edges that are not part of the spanning tree...
	//modified July 15, 2008 (serdar bozdag) Begin (Reason: removal of edges that are not in the MST is not needed in the following steps of this function..
	//	graph_traits < U_Graph >::edge_iterator ei, ei_end;
	//	for ( tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
	//		bool found = false;
	//		for (unsigned int i = 0; i < spanning_tree.size(); i++) {
	//			if (source(*ei, g) == source(spanning_tree[i], g) && target(*ei, g)
	//					== target(spanning_tree[i], g)) {
	//				found = true;
	//				break;
	//			}
	//		}
	//		if (!found) {
	//			remove_edge(*ei, g);
	//		}
	//	}
	//modified July 15, 2008 (serdar bozdag) End
	
	//print the MST of g...
//	string graphizOutputFileName_mst ="output_small-" + contig_id + "_mst.dot";
//	ofstream graphvizOutFile_mst(graphizOutputFileName_mst.c_str(), ios::out);
//	if (!graphvizOutFile_mst) {
//		cerr << "cannot open output file \n";
//		exit(1);
//	}
//	write_graphviz(graphvizOutFile_mst, g,
//			custom_label_writer<U_VertexName>(v_name), make_label_writer<
//					U_EdgeWeight>(e_weight), graph_writer());
//	graphvizOutFile_mst.close();

	//this gives the adjacency graph...
	for (std::vector < U_Edge >::iterator ei = spanning_tree.begin(); ei
			!= spanning_tree.end(); ++ei) {
		//add target to source's list..
		if (adjacency_map.find(v_name[source(*ei, g)]) == adjacency_map.end()) {
			set<string> s;
			s.insert(v_name[target(*ei, g)]);
			adjacency_map[v_name[source(*ei, g)]] = s;
		} else {
			adjacency_map[v_name[source(*ei, g)]].insert(v_name[target(*ei, g)]);
		}
		//add source to target's list...
		if (adjacency_map.find(v_name[target(*ei, g)]) == adjacency_map.end()) {
			set<string> s;
			s.insert(v_name[source(*ei, g)]);
			adjacency_map[v_name[target(*ei, g)]] = s;
		} else {
			adjacency_map[v_name[target(*ei, g)]].insert(v_name[source(*ei, g)]);
		}
	}
}


void assign_to_component_map(string clone,
		map<string, set<string>, strCmp > & adjacency_map,
		map<int, set<string>, less<int> > & component_2_clone_map, int index) {
	//assign the clone to the current component..
	if (component_2_clone_map.find(index) != component_2_clone_map.end()) {
		component_2_clone_map[index].insert(clone);
	} else {
		set<string> s;
		s.insert(clone);
		component_2_clone_map[index] = s;
	}

	//retrive adjacent clones of the clone and remove it from the adjaceny map..
	set <string> adjacent_clones = adjacency_map[clone];
	adjacency_map.erase(clone);

	//for each adjacent clone of the "clone", call the this function to add them to the current component... (if they are still in the adjaceny_map..)
	for (set<string>::iterator s_iter = adjacent_clones.begin(); s_iter
			!= adjacent_clones.end(); s_iter++) {
		if (adjacency_map.find(*s_iter) != adjacency_map.end()) {
			assign_to_component_map(*s_iter, adjacency_map,
					component_2_clone_map, index);
		}
	}
}

void traverse_clones(map<string, set<string>, strCmp > & adjacency_map,
		map<pair<string, string>, double, pairCmp > & overlap_score_exponent,
		vector< vector<string> > & ordered_clones) {

	//first generate connected_components...
	//make a copy of the adjacency map..
	map<string, set<string>, strCmp > adjacency_map_copy(adjacency_map.begin(),
			adjacency_map.end());
	//	cout << "adjacency_map_copy.size(): " << adjacency_map_copy.size() << endl;
	//cout << "traverse clones is called..." << endl;
	map<int, set<string>, less<int> > component_2_clone_map;
	int index = 0;
	while (adjacency_map_copy.size() > 0) {
		//		cout << "index: " << index << endl; 
		//retrieve the first element..
		string clone = adjacency_map_copy.begin()->first;
		//assign this clone and all connnected clones to the current component..
		assign_to_component_map(clone, adjacency_map_copy,
				component_2_clone_map, index);
		//increase the component ID
		index++;
	}
	if (DEBUG) {
		//		cout << "component_2_clone_map.size(): " << component_2_clone_map.size() << endl;
		//print the contents of the components for the debugging purposes...
		for (map<int, set<string>, less<int> >::iterator iter =
				component_2_clone_map.begin(); iter
				!= component_2_clone_map.end(); iter++) {
			cout << iter->first << "\t" << iter->second.size() << "\t";
			for (set<string>::iterator s_iter = iter->second.begin(); s_iter
					!= iter->second.end(); s_iter++) {
				cout << *s_iter << " ";
			}
			cout << endl;
		}
	}

	//find the ordering for each component...
	for (map<int, set<string>, less<int> >::iterator c_iter =
			component_2_clone_map.begin(); c_iter
			!= component_2_clone_map.end(); c_iter++) {
		//int order = 0;
		map<string, pair<int, double>, strCmp> distance_to_fork;//stores the distance of the boundary clones to the fork clones... Also total edge weight from boundary clone to the fork.. (fork clone: a clone with a degree > 2.. boundary clone: degree 1..)
		vector<string> ordered_clones_temp;

		//find the distance of each boundary clone to the fork clone (if there is any, if there is no fork clone then distance will be infinity...)
		for (set<string>::iterator s_iter = c_iter->second.begin(); s_iter
				!= c_iter->second.end(); s_iter++) {
			if (adjacency_map[*s_iter].size() == 1) {
				string boundary_clone = *s_iter;
				string current_clone = *s_iter;
				string previous_clone = "";
				//TO DO: DISTANCE TO THE FURTHEST FORK NEED TO BE COMPUTED...
				bool end = false;
				bool fork_found = false;
				int distance = 0;
				double total_edge_weight = 0;
				while (!end) {
					if (adjacency_map[current_clone].size() > 2) {
						//						cout << "Fork found..." << endl;
						fork_found = true;
						break;
					} else {
						end = true;//just for temporarily.. if next clone is found.. it will be false again...
						for (set<string>::iterator s_iter_2 =
								adjacency_map[current_clone].begin(); s_iter_2
								!= adjacency_map[current_clone].end(); s_iter_2++) {
							//							cout << "*s_iter_2 vs. previous clone: " << *s_iter_2 << "\t" << previous_clone << endl;
							if ((*s_iter_2).compare(previous_clone) != 0) {
								previous_clone = current_clone;
								current_clone = *s_iter_2;
								//								cout << "Current clone: " << current_clone << endl;
								distance++;
								pair<string, string> p;
								p.first
										= min_str(current_clone, previous_clone);
								p.second = max_str(current_clone,
										previous_clone);
								//								cout << "Overlap_score_exponent<" << current_clone << ", " << previous_clone << ">: " << overlap_score_exponent[p] << endl;
								total_edge_weight += overlap_score_exponent[p];

								end = false;//if next clone is found, then it is not END... end should be false again..
								break;
							}
						}
					}
				}
				//				cout << "Done traversing a component.." << endl;
				if (fork_found) {
					pair<int, double> p(distance, total_edge_weight);
					distance_to_fork[boundary_clone] = p;;
				} else {
					pair<int, double> p(BIG_NUMBER, total_edge_weight);
					distance_to_fork[boundary_clone] = p;
				}
			}
		}

		//		cout << "Distance_to_fork is computed..." << endl;
		bool end = false;
		//pick a unvisited boundary clone with the max distance..
		int max_first = 0;
		double min_second= BIG_NUMBER;
		string boundary_clone = "";
		for (map<string, pair<int, double>, strCmp>::iterator iter =
				distance_to_fork.begin(); iter != distance_to_fork.end(); iter++) {
			if (iter->second.first > max_first) {
				max_first = iter->second.first;
				min_second = iter->second.second;
				boundary_clone = iter->first;
			} else if (iter->second.first == max_first) {
				if (iter->second.second < min_second) {
					max_first = iter->second.first;
					min_second = iter->second.second;
					boundary_clone = iter->first;
				}
			}
		}
		//		cout << "Selected boundary clone: " << boundary_clone << endl;
		//traverse until a boundary clone is hit...
		string current_clone = boundary_clone;
		string previous_clone = "";

		while (!end) {

			//store it as ordered clone..
			ordered_clones_temp.push_back(current_clone);

			if (adjacency_map[current_clone].size() > 2) {
				//if it's a fork then select the edge w/ min overlap score and make sure that other edges are removed.....
				double min= BIG_NUMBER;
				string min_clone = "";
				for (set<string>::iterator s_iter =
						adjacency_map[current_clone].begin(); s_iter
						!= adjacency_map[current_clone].end(); s_iter++) {
					if ((*s_iter).compare(previous_clone) != 0) {
						pair<string, string> p;
						p.first = min_str(current_clone, *s_iter);
						p.second = max_str(current_clone, *s_iter);
						if (overlap_score_exponent[p] < min) {
							min = overlap_score_exponent[p];
							min_clone = *s_iter;
						}
					}
				}

				previous_clone = current_clone;
				current_clone = min_clone;
				//				cout << "New current_clone after fork: " << current_clone << endl;
			} else {
				end = true;//just for temporarily.. if next clone is found.. it will be false again...
				for (set<string>::iterator s_iter =
						adjacency_map[current_clone].begin(); s_iter
						!= adjacency_map[current_clone].end(); s_iter++) {
					if ((*s_iter).compare(previous_clone) != 0) {
						previous_clone = current_clone;
						current_clone = *s_iter;
						//						cout << "New current clone: " << current_clone << endl;
						end = false;//if next clone is found, then it is not END... end should be false again..
						break;
					}
				}
			}
		}
		ordered_clones.push_back(ordered_clones_temp);
	}
}

//checks if clones are occuring in at least one component together...
bool is_together(map<int, set<int>, less<int> > & component_ID_to_vertices, const U_VertexName & v_name, const set<string> &clones){
	for (map<int, set<int>, less<int> >::iterator iter = component_ID_to_vertices.begin(); iter != component_ID_to_vertices.end(); iter++) {
		set<string> clone_names;
//		cout << " ********Clone Names:\t";
		for (set<int>::iterator s_iter = iter->second.begin(); s_iter != iter->second.end(); s_iter++) {
//			cout << get_clone_name(v_name[*s_iter]) << " ";
			clone_names.insert(get_clone_name(v_name[*s_iter]));
		}
//		cout << endl;
		//compute the difference between clones & iter->second..
		list<string> diff;
		set_difference(clones.begin(), clones.end(), clone_names.begin(), clone_names.end(), front_inserter(diff));
		
		//if all elements in clones are in clone_names (i.e. diff is empty..) then return true...
		if (diff.size() == 0){
			return true;
		}
	}
	
	//if we reach here then it means clones is never completely together in a component in CFG...
	return false;
}
		
/*
 * Computes the clones that are buried into multiple clones.. If X% or more of a clone's fragments are present in the MFG then it is buried into multiple clones.. 
 *
*/
void compute_clones_buried_into_multiple_clones(const vector<int> &vertex_ID_2_component_ID, const U_VertexName &v_name, map<string, Clone *, strCmp> &cloneObjectByName, vector<string> &buried_clones){
	
	//construct set of vertex ID's for each component...
	map<int, set<int>, less<int> > component_ID_to_vertices;
	get_vertices_per_component(vertex_ID_2_component_ID, component_ID_to_vertices);
	
	
	//compute # matching fragments for each clone...
	map<string, int, strCmp> clone_2_n_matching_fragments;
	count_number_of_matching_fragments(component_ID_to_vertices, v_name, clone_2_n_matching_fragments);
	
	//go over each clone in cloneObjectByName map and store the ones that are (1) in the clone_2_n_matching_fragments and (2) percentage of their matching fragments is above a threshold
	for(map<string, Clone *, strCmp>::iterator iter = cloneObjectByName.begin(); iter != cloneObjectByName.end(); iter++){
		if (clone_2_n_matching_fragments.find(iter->first) != clone_2_n_matching_fragments.end()){
			double percentage = (double)clone_2_n_matching_fragments[iter->first]/(double)iter->second->getNFragments()*100;
			//cout << "percentage: " << percentage << endl;
			if (percentage >= COOPERATIVE_BURY_THRESHOLD){
				buried_clones.push_back(iter->first);
			}
		}
	}
}

/*
 * Computes number of matching fragments for each clone in the MFG.. Go over each connected component >= 2 and count # clones in each one..
*/
void count_number_of_matching_fragments(map<int, set<int>,less<int> > &component_ID_to_vertices, const U_VertexName &v_name, map<string, int, strCmp> &clone_2_n_matching_fragments){
	for(map<int, set<int>, less<int> >::iterator iter = component_ID_to_vertices.begin(); iter != component_ID_to_vertices.end(); iter++){
		if (iter->second.size() >= 2){
			for (set<int>::iterator s_iter = iter->second.begin(); s_iter != iter->second.end(); s_iter++){
				string clone_name = get_clone_name(v_name[*s_iter]);
				if (clone_2_n_matching_fragments.find(clone_name) == clone_2_n_matching_fragments.end()){
					clone_2_n_matching_fragments[clone_name] = 1;
				}else{
					clone_2_n_matching_fragments[clone_name]++;
				}
			}
		}
	}
}


/*
 * Removes the all fragments of a clone in the MFG..
 * 
 */

void remove_clone_from_MFG(U_Graph &g, string &clone_name, U_VertexName &v_name){
	graph_traits< U_Graph >::vertex_iterator vi, vi_end;

	for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
		string name = get_clone_name(v_name[*vi]);
		if (clone_name.compare(name) == 0){
			clear_vertex(*vi, g);
		}
	}
}

/*
 * For each connected component C_i in the MFG, if C_i consists of only candidate clones for bury then store the set of clone names in a vector.. otherwise discard C_i.. 
 * 
 * If a candidate clone is always together with non-candidate clone in the MFG, then this clone will not occur in the vector.. Such clones can be buried withouth any worry..
 */
void generateCandidateSet(const vector<int>& vertex_ID_2_component_ID, const U_VertexName &v_name, const vector<string> &candidate_clones_for_bury, vector<set<string> > & candidate_clone_set, set<string> &clones_in_candidate_clone_set){
	
	//construct set of vertex ID's for each component...
	map<int, set<int>, less<int> > component_ID_to_vertices;
	get_vertices_per_component(vertex_ID_2_component_ID, component_ID_to_vertices);
	
	
	//go over each component and generate candidate_clone_set..
	for(map<int, set<int>, less<int> >::iterator iter = component_ID_to_vertices.begin(); iter != component_ID_to_vertices.end(); iter++){
		if (iter->second.size() >=2 ){
			set<string> clone_set;
			bool ignore = false; //do not store the clone_set if true.. If there is a clone in the component that is not a candidate then set to TRUE..
			for (set<int>::iterator s_iter = iter->second.begin(); s_iter != iter->second.end(); s_iter++){
				string clone_name = get_clone_name(v_name[*s_iter]);
				//check if this clone is a candidate clone for bury..
				if (is_candidate_clone_for_bury(clone_name, candidate_clones_for_bury)){
					clone_set.insert(clone_name);
					clones_in_candidate_clone_set.insert(clone_name);
				}else{
					ignore = true;
					break;
				}
			}
			if (!ignore){
				candidate_clone_set.push_back(clone_set);
			}
		}
	}
	
}

/*
 * Searces the first parameter in the vector of the second parameter.. returns T or F accordingly..
 * 
 * */
bool is_candidate_clone_for_bury(const string &clone_name, const vector<string> candidate_clones_for_bury){
	for (unsigned int i = 0; i < candidate_clones_for_bury.size(); i++){
		if (clone_name.compare(candidate_clones_for_bury[i]) == 0){
			return true;
		}
	}
	return false;
}

/*
 * Computes the minimum hitting set of a vector of sets by using a greedy approach..
 * 
 * First the element that occur in maximum number of sets is selected..
 */
void compute_minimum_hitting_set(vector<set<string> > &candidate_clone_set, set<string> &minimum_hitting_set){
	priority_queue< CloneFrequency > clone_frequency_queue;
	bool end = false;
	
	while(!end){//exit the loop when the priority queue is empty...
		//first compute how many times each clone occur in the candidate_clone_set..
		map<string, int, strCmp> clone_frequency; 
		for (unsigned int i = 0; i < candidate_clone_set.size(); i++){
			for (set<string>::iterator s_iter = candidate_clone_set[i].begin(); s_iter != candidate_clone_set[i].end(); s_iter++){
				if (clone_frequency.find(*s_iter) == clone_frequency.end()){
					clone_frequency[*s_iter] = 1;
				}else{
					clone_frequency[*s_iter]++;
				}
			}
		}
		
		
		//store the clone_frequency into the priority queue..
		for (map<string, int, strCmp>::iterator iter = clone_frequency.begin(); iter != clone_frequency.end(); iter++){
			CloneFrequency c_f(iter->first, iter->second);
			clone_frequency_queue.push(c_f);
		}
		
		//exit the while loop iff the priority queue is empty..
		if (clone_frequency_queue.empty()){
			break;
		}
		
		//pop the first element from the priority queue...
		CloneFrequency c_f = clone_frequency_queue.top();
		clone_frequency_queue.pop();
		
		//clear all candidate clone set that contain p.first (i.e. the most frequent clone)
		for (unsigned int i = 0; i < candidate_clone_set.size(); i++){
			if (candidate_clone_set[i].find(c_f.getCloneName()) != candidate_clone_set[i].end()){
				candidate_clone_set[i].clear();//remove the all elements in the set, since it's covered by p.first..
			}
		}
		
		//store p.first in minimum hitting set vector..
		minimum_hitting_set.insert(c_f.getCloneName());
		
		//clear the priority queue for the next iteration..
		while(!clone_frequency_queue.empty()){
			clone_frequency_queue.pop();
		}
	}
}


