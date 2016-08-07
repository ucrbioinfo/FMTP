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
#include <queue>
#include <iomanip>
#include <boost/graph/edmunds_karp_max_flow.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/vector_as_graph.hpp>

#include "Contig.h"

//#define HIGH_THRESHOLD 1e-5 //if cutoff threshold is greater than or equal to this, then it's considered high.. (NOT USED ANYMORE.. EXCEPT FOR ITERATION_ID >= 1)

int DEBUG = 0;
double COOPERATIVE_BURY_THRESHOLD = 80; //used for burying clone into multiple clones..
double BURY_THRESHOLD = 80;
double EXTRA_FRAGMENT_PERCENTAGE = 30; //if extra fragments of a node (i.e. match nowhere in clone - fragment graph) is more than this percentage then they are added to clone-fragment graph...
double MTP_ILP_VERY_HIGH_THRESHOLD = -1; //5 for barley, 2.5 for rice.. if -log(Sulston(x,y)) is lower than this value then x & y are definitely not overlapping...




//#define COOPERATIVE_BURY_THRESHOLD 80


//int WINDOW_SIZE = 0;

//int big_cluster_threshold = 10;//for debugging purposes...
//int MIN_CONNECTIVITY = 75; // number_of_edges/number_of_ideal_edges in a component.. if for a components, this value is lower than this threshold then that component is splitted..

//and we make sure that the clone of these fragment are selected in the solution...

void usage() {
	cout << "Compute_MTP: Generates a linear programming model to select MTP clones of contigs in a FPC map..\n";
	cout << "-f Fingerprinting Method.. 1: HICF, 2: Agarose\n";
	cout << "-c ctg_BAC_clone file\n";
	cout << "-s size file\n";
	cout << "-b bury threshold (default 80)\n";
	cout << "-B cooperative bury threshold (default 80)\n";
	cout << "-t sulston_score_threshold.. Fragments of clones will be matched if their sulston score is lower than this threshold..(Default values agarose:1e-10, HICF:1e-25)\n";
	cout << "-T Tolerance (rice:7, barley:3, cowpea:5)..\n";
	cout << "-G Gellength (rice: 3300, barley:18000, cowpea: 36,000)..\n";
	cout << "-o output_file base file (without extension).. Linear Programming Model file (for GLPK, in MathProg language)\n";
	cout << "-O partial ctg_2_mtp_file\n";
	cout << "-v verbose level min 0 max 4\n";
	exit(0);
}

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

/*struct pairCmp4 {
	bool operator()(const pair<string, int> p1, const pair<string, int> p2) const {
		return p1.second < p2.second;
	}
};*/

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
bool exists_weak_edge(const set<int> &vertices_in_component, const U_VertexName & v_name, map<pair<string, string>, double, pairCmp > & overlap_score_exponent);
void compute_clones_buried_into_multiple_clones(const vector<int> &vertex_ID_2_component_ID, const U_VertexName &v_name, map<string, Clone *, strCmp> &cloneObjectByName, vector<string> &buried_clones);
void count_number_of_matching_fragments(map<int, set<int>,less<int> > &component_ID_to_vertices, const U_VertexName &v_name, map<string, int, strCmp> &clone_2_n_matching_fragments);
void remove_clone_from_MFG(U_Graph &g, string &clone_name, U_VertexName &v_name);
void generateCandidateSet(const vector<int>& vertex_ID_2_component_ID, const U_VertexName &v_name, const vector<string> &candidate_clones_for_bury, vector<set<string> > & candidate_clone_set, set<string> &clones_in_candidate_clone_set);
void compute_minimum_hitting_set(vector<set<string> > &candidate_clone_set, set<string> &minimum_hitting_set);
bool is_candidate_clone_for_bury(const string &clone_name, const vector<string> candidate_clones_for_bury);



//TODO: Following functions are added for testing purposess only...
//debug-begin
//struct cloneCmp {
//	bool operator()(const Clone *c1, const Clone *c2) const {
//		return c1->getName().compare(c2->getName()) < 0;
//	}
//};
//
//void concatenate_clones(const set<string> & clone_names,
//						string & concatenated_clone_name, const string& delimiters) {
//	concatenated_clone_name = "";
//	for (set<string>::iterator s_iter = clone_names.begin(); s_iter
//		 != clone_names.end(); s_iter++) {
//		concatenated_clone_name += *s_iter + delimiters;
//	}
//}
//
//void get_clone_components(const vector<int> & vertex_ID_2_component_ID,
//						  const U_VertexName & v_name, map<string, int, strCmp> &clone_components) {
//	//construct set of vertex ID's for each component...
//	map<int, set<int>, less<int> > component_ID_to_vertices;
//	get_vertices_per_component(vertex_ID_2_component_ID,
//							   component_ID_to_vertices);
//	
//	for (map<int, set<int>, less<int> >::iterator iter =
//		 component_ID_to_vertices.begin(); iter
//		 != component_ID_to_vertices.end(); iter++) {
//		set<string> clone_names;
//		for (set<int>::iterator s_iter = iter->second.begin(); s_iter
//			 != iter->second.end(); s_iter++) {
//			clone_names.insert(get_clone_name(v_name[*s_iter]));
//		}
//		//concatenate clone names...
//		string concatenated_clone_name = "";
//		concatenate_clones(clone_names, concatenated_clone_name,"-");
//		clone_components[concatenated_clone_name] = 1;
//	}
//	
//}
//
//
//void tokenize(const string& str, vector<string>& tokens,
//			  const string& delimiters) {
//	// Skip delimiters at beginning.
//	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
//	// Find first "non-delimiter".
//	string::size_type pos = str.find_first_of(delimiters, lastPos);
//	
//	while (string::npos != pos || string::npos != lastPos) {
//		// Found a token, add it to the vector.
//		tokens.push_back(str.substr(lastPos, pos - lastPos));
//		// Skip delimiters.  Note the "not_of"
//		lastPos = str.find_first_not_of(delimiters, pos);
//		// Find next "non-delimiter"
//		pos = str.find_first_of(delimiters, lastPos);
//	}
//}
//
//
//void compute_connected_clone_list(map<string, int, strCmp> &clone_components,
//								  map<string, set<string>, strCmp> &connected_clones_map) {
//	
//	for (map<string, int, strCmp>::iterator iter = clone_components.begin(); iter
//		 != clone_components.end(); iter++) {
//		vector<string> clone_names;
//		tokenize(iter->first, clone_names,"-");
//		//if (clone_names.size() == size){//consider components of a given size only...
//		for (unsigned int i = 0; i < clone_names.size(); i++) {
//			if (connected_clones_map.find(clone_names[i])
//				== connected_clones_map.end()) {
//				set<string> s;
//				connected_clones_map[clone_names[i]] = s;
//			}
//			for (unsigned int j = 0; j < clone_names.size(); j++) {
//				if (i != j) {
//					connected_clones_map[clone_names[i]].insert(clone_names[j]);
//				}
//			}
//		}
//		//}
//	}
//}
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

int main(int argc, char**argv) {
	map<string, Clone *, strCmp> cloneObjectByName; // stores the pointer to the Clone object of each clone (by clone name)
	map<string, Contig *, strCmp> contigObjectByName; // stores the pointer to the Contig object of each contig (by contig name)
	map<string, Band *, strCmp> fragmentObjectByName; //stores the pointer to the Fragment object of each fragment (by fragment name)..

	char ctgBACFile[FILENAME_MAX];
	char sizeFile[FILENAME_MAX];
	char outputFile[FILENAME_MAX];
	char outputCtgMTPFile[FILENAME_MAX];

	int ctgBACFilePresent = 0, sizeFilePresent = 0;
	double sulston_score_threshold = -1.0;

	int fingerprinting_method = -1;
	int tolerance = -1, gellen = -1;	
	
//	//TODO: added for testing purposes only...
//	char clone_coordinates_file[FILENAME_MAX];

	
	map<pair<string, string>, int, pairCmp> E_f; //stores the edges between fragments.. (if an fragment is covered by selecting its clone, all fragments adjacent to this fragment in E_f are also covered..
	map<pair<string, string>, int, pairCmp> E; //stores the edges between clones and fragments..
	//map<string, int, strCmp> fragmentDegree; //stores the degree of a fragment.. if degree = 1.. it's connected to a clone and a fragment, if degree = 0, only connected to a clone..
	map<pair<string, string>, double, pairCmp > overlap_score_exponent;//stores -ln(sulston(A,B)) between all pairs of "overlapping" clones...	
	while (1) {
		char c = getopt(argc, argv, "hc:s:f:b:B:t:T:G:o:O:v:");

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
		case 'T':
			tolerance = atoi(optarg);
			printf("Using tolerance value %d\n", tolerance);
			break;
		case 'G':
			gellen = atol(optarg);
			printf("Using gel length value %d\n", gellen);
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
			printf("COOPERATIVE_BURY_THRESHOLD : %f\n", COOPERATIVE_BURY_THRESHOLD);
			break;		
		case 't':
			sulston_score_threshold = atof(optarg);
			printf("Sulston score threshold: %e\n", sulston_score_threshold);
			break;
//		case 'i':
//			iteration_ID = atoi(optarg);
//			printf("Iteration ID: %d\n", iteration_ID);
//			break;
		case 'o':
			strcpy(outputFile, optarg);
			printf("Will write output data to '%s'\n", outputFile);
			break;
		case 'O':
			strcpy(outputCtgMTPFile, optarg);
			printf("Will write partial mtp output data to '%s'\n",
					outputCtgMTPFile);
			break;
		case 'v':
			DEBUG = atoi(optarg);
			printf("Verbose level: %d\n", DEBUG);
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

	if (strlen(outputFile) == 0 || strlen(outputCtgMTPFile) == 0) {
		cerr << "output file missing\n";
		usage();
	}

	
	if (fingerprinting_method == HICF){
		MTP_ILP_VERY_HIGH_THRESHOLD = MTP_ILP_VERY_HIGH_THRESHOLD_HICF;
		//update connected_component_size_threshold if it's not set by user..
//		if (p_value == -1){
//			p_value = MTP_ILP_HICF_Q;
//		}
		//update sulston_score_threshold if not set by the user..
		if (sulston_score_threshold == -1){
			sulston_score_threshold = MTP_ILP_HICF_C;
		}
	}else if (fingerprinting_method == AGAROSE){
		MTP_ILP_VERY_HIGH_THRESHOLD = MTP_ILP_VERY_HIGH_THRESHOLD_AGAROSE;
		//update connected_component_size_threshold if it's not set by user..
//		if (p_value == -1){
//			p_value = MTP_ILP_AGAROSE_Q;
//		}
		//update sulston_score_threshold if not set by the user..
		if (sulston_score_threshold == -1){
			sulston_score_threshold = MTP_ILP_AGAROSE_C;
		}
	}else{
		cerr << "Invalid fingerprinting code\n";
		usage();
	}
	
	
	if ((sulston_score_threshold < 0) || (sulston_score_threshold > 1)) {
		cerr << "invalid sulston score threshold.. must be between 0 and 1\n";
		usage();
	}

	if (tolerance < 0 || gellen <= 0){
		cerr << "Invalid tolerance/gellen...\n";
		usage();
	}
	
	if (BURY_THRESHOLD < 0 || BURY_THRESHOLD > 100){
		cerr << "Invalid range for BURY_THRESHOLD.. Should be between [0,100]";
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

		tok2 = strtok (NULL, " ");

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

			for (vector<Band *>::iterator iter = tmpFragments.begin(); iter != tmpFragments.end(); iter++) {
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
//	//DEBUG--END
	
if (DEBUG > 3) {
	cout << "Priniting ctg-BAC table\n";
	for (map<string, Contig *, strCmp>::iterator iter = contigObjectByName.begin(); iter != contigObjectByName.end(); iter++) {
		iter->second->print();
	}
}

	
	
	
//remove buried clones...


	
	//TODO: following for loop i=1 to 30 is added for testing purposes only...
	//TODO: BURY_THRESHOLD_ORG added for testing purposes only..
//	double BURY_THRESHOLD_ORG = BURY_THRESHOLD;
//	for (int i= 0; i < 20; i++){
		//TODO: followint tp, fp.. variables are added for testing purposes only...
	//	int tp = 0, fp = 0, fn = 0, tn = 0;
//		BURY_THRESHOLD = BURY_THRESHOLD_ORG + i;
for (map<string, Contig *, strCmp>::iterator con_iter = contigObjectByName.begin(); con_iter != contigObjectByName.end(); con_iter++) {
	set<Clone *, Contig::cloneCmp> clones = con_iter->second->getClones();
	for (set<Clone *, Contig::cloneCmp>::iterator cl_iter_1 = clones.begin(); cl_iter_1 != clones.end(); cl_iter_1++) {
		multiset<int> f_clone_1 = (*cl_iter_1)->getFragmentSizes();
		//cout << "f_clone_1.size(): " << f_clone_1.size() << endl;
		for (set<Clone *, Contig::cloneCmp>::iterator cl_iter_2 = clones.upper_bound(*cl_iter_1); cl_iter_2 != clones.end(); cl_iter_2++) {
			multiset<int> f_clone_2 = (*cl_iter_2)->getFragmentSizes();
			//cout << "f_clone_2.size(): " << f_clone_2.size() << endl;
			
			int total_n_shared = numberOfSharedBands(tolerance, f_clone_1, f_clone_2);
			if (DEBUG > 3) {
				cout << (*cl_iter_1)->getName() << ": " << f_clone_1.size() << "\t" << (*cl_iter_2)->getName() << ": " << f_clone_2.size() << "\t" << total_n_shared << endl;
			}
			
//			TODO: variable buried_observed and ISBURIED added for testing purposes only..
		//	bool ISBURIED = false;
//			bool buried_observed = false;
			if (f_clone_1.size() > f_clone_2.size()) { //try to bury clone_2
				if ((double)total_n_shared/(double)f_clone_2.size()*100 > BURY_THRESHOLD) {
					if (DEBUG > 3) {
						cout << (*cl_iter_2)->getName() << " buried into " << (*cl_iter_1)->getName() << endl;
					}
					(*cl_iter_2)->setParent(*cl_iter_1);
					//buried_observed = true;
				}
			} else {//try to bury clone_1
				if ((double)total_n_shared/(double)f_clone_1.size()*100 > BURY_THRESHOLD) {
					if (DEBUG > 3) {
						cout << (*cl_iter_1)->getName() << " buried into " << (*cl_iter_2)->getName() << endl;
					}
					(*cl_iter_1)->setParent(*cl_iter_2);
					//buried_observed = true;
				}
			}
			//BEGIN--DEBUG
//			ISBURIED = is_buried(clone_2_coordinates[*cl_iter_1], clone_2_coordinates[*cl_iter_2]);
//			if (ISBURIED){
//				double score;
//				int nShared = numberOfSharedBands(tolerance, (*cl_iter_1)->getFragmentSizes(), (*cl_iter_2)->getFragmentSizes());
//				if (nShared == 0) {
//					score = 1.0;
//				} else {
//					fastSulston((*cl_iter_1)->getNFragments(), (*cl_iter_2)->getNFragments(), nShared, tolerance, gellen, score);
//				}
//				cout << score << endl;
//			}
//			int overlap_size = is_overlapping(clone_2_coordinates[*cl_iter_1], clone_2_coordinates[*cl_iter_2]);
//			if (ISBURIED == false && overlap_size > 60000) continue;
//			if (ISBURIED){
//				if(buried_observed){
//					tp++;
//				}else{
//					fn++;
					//if ((double)total_n_shared/(double)f_clone_2.size()*100 > (double)total_n_shared/(double)f_clone_1.size()*100){
//						cout << "FN____\t" << (double)total_n_shared/(double)f_clone_2.size()*100 << "\t" << f_clone_2.size() <<  endl;
//					}else{
//						cout << "FN____\t" << (double)total_n_shared/(double)f_clone_1.size()*100 << "\t" << f_clone_1.size() << endl;
//					}
//				}
//			}else{
//				if (buried_observed){
//					fp++;
//					if ((double)total_n_shared/(double)f_clone_2.size()*100 > (double)total_n_shared/(double)f_clone_1.size()*100){
//						cout << "FP____\t" << (double)total_n_shared/(double)f_clone_2.size()*100 << "\t" << f_clone_2.size() << endl;
//					}else{
//						cout << "FP____\t" << (double)total_n_shared/(double)f_clone_1.size()*100 << "\t" << f_clone_1.size() << endl;
//					}
//				}else{
//					tn++;
//				}
//			}
		//END--DEBUG
		}
	}
}
	
	//TODO: exit command is added below for testing purpose..
//	double tpr = (double)tp/(double)(tp+fn);
//	double fpr = (double)fp/(double)(fp+tn);
//	cout << BURY_THRESHOLD << "\t" << fpr << "\t" << tpr << endl;
//}
//	exit(-1);

	
ofstream ctgMTPFile(outputCtgMTPFile, ios::out);
if (!ctgMTPFile) {
	cerr << "Cannot open ctg_2_mtp file\n";
	exit(1);
}

//Go over all contigs and find MTP for each of them....
for (map<string, Contig *, strCmp>::iterator con_iter = contigObjectByName.begin(); con_iter != contigObjectByName.end(); con_iter++) {
	cout << "contig id: " << con_iter->first << endl;
	set<Clone *, Contig::cloneCmp> clones = con_iter->second->getClones();
	//map<string, bool, strCmp> clone_2_ismatched;//stores if a clone matches to any other clone in the contig...
	int contig_size = 0;
	for (set<Clone *, Contig::cloneCmp>::iterator cl_iter_1 = clones.begin(); cl_iter_1 != clones.end(); cl_iter_1++) {
		//cout << "clonename : " << (*cl_iter_1)->getName() << endl;
		if ((*cl_iter_1)->isBuried()) {continue;}
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

//if there are only 2 canonical clones then print all of them and ignore this contig in further processing...
	if (contig_size == 2){
		ctgMTPFile << con_iter->first << "\t" << contig_size << "\t";
		for (set<Clone *, Contig::cloneCmp>::iterator cl_iter_1 = clones.begin(); cl_iter_1 != clones.end(); cl_iter_1++) {	
			if ((*cl_iter_1)->isBuried()) {continue;}
			ctgMTPFile << (*cl_iter_1)->getName() << " ";
		}
		ctgMTPFile << "\n";
		continue;
	}
	
	
	//compare canonical clones (all to all)

	for (set<Clone *, Contig::cloneCmp>::iterator cl_iter_1 = clones.begin(); cl_iter_1 != clones.end(); cl_iter_1++) {
		if ((*cl_iter_1)->isBuried()) {continue;}
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
		for (set<Clone *, Contig::cloneCmp>::iterator cl_iter_2 = clones.upper_bound(*cl_iter_1); cl_iter_2 != clones.end(); cl_iter_2++) {
			if ((*cl_iter_2)->isBuried()) {
//				cout << "buried_clone : " << (*cl_iter_2)->getName() << "... buried into " << (*cl_iter_2)->getParent()->getName() << endl;				
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
			//cout << (*cl_iter_1)->getName() << "\t" << (*cl_iter_2)->getName() << "\t" << "score : " << score << endl;
			//if the sulston score is lower than the threshold then match the fragments of these clones and store each match in E_fragment..
			//store the score in overlap_score_exponent...
			pair<string, string> p;
			p.first = min_str((*cl_iter_1)->getName(), (*cl_iter_2)->getName());
			p.second = max_str((*cl_iter_1)->getName(), (*cl_iter_2)->getName());
			overlap_score_exponent[p] = -log10(score);

			if (score <= sulston_score_threshold) {
				if (DEBUG > 3) {
					cout << (*cl_iter_1)->getName() << "\t" << (*cl_iter_2)->getName() << "\t" << "score : " << overlap_score_exponent[p] << endl;
				}

				//cout << overlap_score_exponent[p] << endl;
				//clone_2_ismatched[(*cl_iter_1)->getName()] = true;
				//clone_2_ismatched[(*cl_iter_2)->getName()] = true;
				//ismatched = true;
				//generate nodes for each fragment of both clones..
				Graph g;
				NameVertexMap L, R;

				BandPtr band_ptr = get(&VertexProperties::band_ptr, g);

				property_map < Graph, edge_capacity_t >::type capacity = get(edge_capacity, g);
				property_map < Graph, edge_reverse_t >::type rev = get(edge_reverse, g);
				property_map < Graph, edge_residual_capacity_t >::type residual_capacity = get(edge_residual_capacity, g);
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
						int diff = bands_1[i]->getSize() - bands_2[j]->getSize();
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

					string graphizOutputFileName_t = ssGraphvizOutputFileName_t.str() + "-" + con_iter->first + "-" + (*cl_iter_1)->getName() + "-" + (*cl_iter_2)->getName() + "_mbp.dot";
					ofstream graphvizOutFile_t(graphizOutputFileName_t.c_str(), ios::out);
					if (!graphvizOutFile_t) {
						cerr << "cannot open output file \n";
						exit(1);
					}

					write_graphviz(graphvizOutFile_t, g, custom_label_writer_2<BandPtr>(band_ptr), default_writer(), graph_writer());
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
					tie (e1, b1) = add_edge(source, iter->second, g);
					tie (e2, b2) = add_edge(iter->second, source, g);
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
					cout << (*cl_iter_1)->getName() << "-" << (*cl_iter_2)->getName() << "\tmatch : " << match << endl;
				}
				Graph::edge_iterator ei, ei_end;
				for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
					if (band_ptr[boost::source(*ei, g)] != NULL && band_ptr[boost::target(*ei, g)] != NULL) {
						if (residual_capacity[*ei] == 0 && capacity[*ei] > 0) {
							Vertex s = boost::source(*ei, g);
							Vertex t = boost::target(*ei, g);

							//add matched fragments to E_fragment..
							pair<string, string> p(band_ptr[s]->getName(), band_ptr[t]->getName());
							if (E_f.find(p) != E_f.end()) {
								cerr << "Error: fragment pair is already added..\n";
								exit(0);
							}
							E_f[p] = 1;
						}
					}
				}
			}
		}
	}

	//generate clone-fragment graph for E_f.. 
	U_Graph g_f;

	NameUVertexMap name_2_vertex;
	U_VertexName v_name = get(vertex_name, g_f);
	U_EdgeWeight e_weight = get(edge_weight, g_f);
	for (map<pair<string, string>, int, pairCmp >::iterator iter = E_f.begin(); iter != E_f.end(); iter++) {
		//cout << iter->first.first << "-" << iter->first.second << endl;
		U_Vertex u,v;
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
		tie(e,b) = add_edge(u, v, g_f);
		if (!b) {
			cerr << "cannot add edge between " << v_name[u] << " and " << v_name[v] << endl;
			exit(1);
		}
		string clone_1 = get_clone_name(v_name[u]);
		string clone_2 = get_clone_name(v_name[v]);
		pair<string, string> p;
		p.first = min_str(clone_1, clone_2);
		p.second = max_str(clone_1, clone_2);
		e_weight[e] = overlap_score_exponent[p];
	}

	//generate graphviz output after generating initial MFG...
	stringstream ssGraphvizOutputFileName;
	if (DEBUG > 3){

		ssGraphvizOutputFileName << outputFile;
		string graphizOutputFileName = ssGraphvizOutputFileName.str() + "-" + con_iter->first + "-initialMFG.dot";
		ofstream graphvizOutFile(graphizOutputFileName.c_str(), ios::out);
		if (!graphvizOutFile) {
			cerr << "cannot open output file \n";
			exit(1);
		}

		write_graphviz(graphvizOutFile, g_f, custom_label_writer<U_VertexName>(v_name), make_label_writer<U_EdgeWeight>(e_weight), graph_writer());
		graphvizOutFile.close();
	}
	
	//find biconnected components
	/*property_map < U_Graph, edge_component_t >::type component = get(edge_component, g_f);
	 biconnected_components(g_f, component);
	 //std::cerr << "Found " << num_comps << " biconnected components.\n";

	 std::vector<U_Vertex> art_points;
	 articulation_points(g_f, std::back_inserter(art_points));
	 //std::cerr << "Found " << art_points.size() << " articulation points.\n";
	 
	 */

	/*//output_the_graph...
	 
	 
	 graphvizOutFile_2 << "graph A {\n graph [size=\"6,6\"]\n" << "  node[shape=\"circle\"]\n";
	 
	 graph_traits< U_Graph >::vertex_iterator vi, vi_end;
	 for (tie(vi, vi_end) = vertices(g_f); vi != vi_end; vi++){
	 graphvizOutFile_2 << *vi << "[label=\"" << v_name[*vi] << "\"];\n"; 
	 }
	 
	 for (std::size_t i = 0; i < art_points.size(); ++i) {
	 graphvizOutFile_2 << art_points[i] << " [ style=\"filled\", fillcolor=\"red\" ];"<< std::endl;
	 }
	 
	 graph_traits < U_Graph >::edge_iterator ei, ei_end;
	 for (tie(ei, ei_end) = edges(g_f); ei != ei_end; ++ei){
	 graphvizOutFile_2 << source(*ei, g_f) << " -- " << target(*ei, g_f) << ";\n";
	 }
	 graphvizOutFile_2 << "}\n";
	 */

	//******************************...graph pruning...******************************************
	//Get rid of articulation points (i.e. vertices that is in more than one biconnected components)

	//for (std::size_t i = 0; i < art_points.size(); ++i) {
	//clear_vertex(art_points[i], g_f);
	//remove_vertex(art_points[i], g_f);
	//}


	//cout << "Articulations points have been removed from the graph" << endl;
	//write the pruned graph into a file..
	/*string graphizOutputFileName_3 = ssGraphvizOutputFileName.str() + "-" + con_iter->first + "_bc.dot";
	 ofstream graphvizOutFile_3(graphizOutputFileName_3.c_str(), ios::out);
	 if (!graphvizOutFile_3){
	 cerr << "cannot open output file \n";
	 exit(1);
	 }
	 write_graphviz(graphvizOutFile_3, g_f, custom_label_writer<U_VertexName>(v_name), default_writer(), graph_writer());
	 graphvizOutFile_3.close();
	 */

	//compute the connected components..
	std::vector<int> vertex_ID_2_component_ID(num_vertices(g_f), -1);//uses the old number of vertices...
	//map<int, int, less<int> > vertex_ID_2_component_ID;
	connected_components(g_f, make_iterator_property_map(vertex_ID_2_component_ID.begin(), get(vertex_index, g_f)));

	map<int, int, less<int> > component_ID_2_size;

	for (std::vector<int>::size_type i = 0; i != vertex_ID_2_component_ID.size(); ++i) {

		if (component_ID_2_size.find(vertex_ID_2_component_ID[i]) == component_ID_2_size.end()) {
			component_ID_2_size[vertex_ID_2_component_ID[i]] = 1;
		} else {
			component_ID_2_size[vertex_ID_2_component_ID[i]]++;
		}
	}

	
	//cout << "Computing which components should be splitted in contig " << con_iter->first << endl;
	map<int, set<int>, less<int> > component_ID_2_vertices;
	get_vertices_per_component(vertex_ID_2_component_ID, component_ID_2_vertices);

	int maxID = component_ID_2_vertices.size() + 1;
	int n_cases_weak_overlap = 0, n_cases_extra_fragment = 0, n_cases_unmatched_fragments = 0, n_total_split = 0;
	
	for (map<int, set<int>, less<int> >::iterator iter = component_ID_2_vertices.begin(); iter != component_ID_2_vertices.end(); iter++) {
		map<int, int, less<int> > color;
		if (component_ID_2_size[iter->first] == 1) {continue;}
		bool b1 = exists_weak_edge(iter->second, v_name, overlap_score_exponent);
		bool b2 =  is_fragment_values_not_within_tolerance(iter->second, v_name, tolerance);
		bool b3 = is_same_clone_together_in_component(iter->second, v_name);
		
		//counting number of each cases for MFG pruning..
		if (b1) n_cases_weak_overlap++;
		if (b2) n_cases_unmatched_fragments++;
		if (b3) n_cases_extra_fragment++;
		
		
		if (b1 || b2 || b3) {
			n_total_split++;
			//if (is_fragment_values_not_within_tolerance(iter->second, v_name) || is_same_clone_together_in_component(iter->second, v_name)){
			//if (DEBUG > 2) {
//				cout << "splitting component.. ID: " << v_name[*(component_ID_2_vertices[iter->first].begin())] << endl;
//			cout << "b1: " << b1 << endl;
//			cout << "b2: " << b2 << endl;
//			cout << "b3: " << b3 << endl;
			//}
			double min_flow = split_component(iter->second, g_f, e_weight, color);
			set<int> component_1;
			set<int> component_2;
			
			
			for(map<int, int, less<int> >::iterator color_iter = color.begin(); color_iter != color.end(); color_iter++) {
				if (color_iter->second == 0) {
					component_1.insert(color_iter->first);
				} else if (color_iter->second == 4) {
					component_2.insert(color_iter->first);
				} else {
					cerr << "Unknown color code: " << color_iter->second << endl;
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
			//TODO:degree_th_2 is not used now, but may be used...
			double degree_th_1 = 3;
			double degree_th_2 = 0;
			if (component_1.size() == 1) {
				degree = n_removed_edges;
				degree_th_2 = (double)component_2.size()/2.0;
			}else if (component_2.size() == 1){
				degree = n_removed_edges;
				degree_th_2 = (double)component_1.size()/2.0;
			}
			//double connectivity = (double)degree/double(component_1.size() + component_2.size() -1);
//			if (connectivity != 0){
//				cout <<"connectivity: " << connectivity << endl;
//			}
			
			if ( average_min_flow <= min_flow_threshold && degree < degree_th_1){
//				cout << "splitting..." << endl;
				for (set<int>::iterator set_iter_1 = component_1.begin(); set_iter_1 != component_1.end(); set_iter_1++) {
					for (set<int>::iterator set_iter_2 = component_2.begin(); set_iter_2 != component_2.end(); set_iter_2++) {
						bool b;
						U_Edge e;
						tie(e, b) = edge(*set_iter_1, *set_iter_2, g_f);
						if (b) {
							//cout << "Removing edge between " << v_name[*set_iter_1] << " and " << v_name[*set_iter_2] << endl;
							remove_edge(*set_iter_1, *set_iter_2, g_f);
						}
					}
				}
				if (DEBUG > 2) {
					cout << "Printing new components..." << endl;
					cout << "Component 1: " << endl;
					for (set<int>::iterator t_iter = component_1.begin(); t_iter != component_1.end(); t_iter++) {
						cout << v_name[*t_iter] << endl;
					}
					cout << "Component 2: " << endl;
					for (set<int>::iterator t_iter = component_2.begin(); t_iter != component_2.end(); t_iter++) {
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
	}
//	if (n_total_split > 0){
//		cout << "Printing Number of Cases for MFG Pruning..." << endl;
//		cout << "Total number of splits..." << n_total_split << endl;
//		cout << "Extra Fragment: " << n_cases_extra_fragment << endl;
//		cout << "Unmatched fragments: " << n_cases_unmatched_fragments << endl;
//		cout << "Weak overlap: " << n_cases_weak_overlap << endl;
//	}
	
//print MFG after pruning...
	if (DEBUG > 3){
		stringstream ssGraphvizOutputFileName_temp1;
		ssGraphvizOutputFileName_temp1 << outputFile;
		
		string graphizOutputFileName_temp1 = ssGraphvizOutputFileName_temp1.str() + "-" + con_iter->first + "_prunedMFG.dot";
		ofstream graphvizOutFile_temp1(graphizOutputFileName_temp1.c_str(), ios::out);
		if (!graphvizOutFile_temp1) {
			cerr << "cannot open output file \n";
			exit(1);
		}
		
		write_graphviz(graphvizOutFile_temp1, g_f, custom_label_writer<U_VertexName>(v_name), make_label_writer<U_EdgeWeight>(e_weight), graph_writer());
		graphvizOutFile_temp1.close();	

	}
	
	
	/*****************Get rid of spurious components*****************************************************/
	graph_traits< U_Graph >::vertex_iterator vi, vi_end;
	//TODO: Removal of spurious components is disabled for testing purposes!!!!
	//BEGIN
	/*
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

	//cout << "Connected components have been computed.." << endl;
	//cout << "Clone component size is being computed\n" << endl;
	map<int, int, less<int> > component_ID_2_clone_component_frequecy;
	compute_clone_component_size(vertex_ID_2_component_ID, v_name, component_ID_2_clone_component_frequecy);

	//Get rid of spurious components...


	for (tie(vi, vi_end) = vertices(g_f); vi != vi_end; vi++) {
		int s = component_ID_2_size[vertex_ID_2_component_ID[*vi]];
		int nf = component_ID_2_clone_component_frequecy[vertex_ID_2_component_ID[*vi]];
		char chr_key[20];
		sprintf(chr_key, "%d-%d-%d", contig_size, nf, s);
		string key = chr_key;
		//TODO: using old spurious cc removal procedure for testing...
		//if (spurious_component_probability_table.find(key) != spurious_component_probability_table.end()){
//			if (spurious_component_probability_table[key] > p_value){
				//if ((s-1)*nf >= 3){
////					cout << "--------removing spurious component which would not removed by previous approach...!!: " << key << endl;
////					cout << "previous_value: " << (s - 1)*nf << endl;
////				}
//				clear_vertex(*vi, g_f);
//			}else{
//				//if ((s-1)*nf < 3){
////					cout << "++++++++++++++keeping component which would be removed by previous approeach.." << endl;
////					cout << "previous_value: " << (s - 1)*nf << endl;
////				}
//			}
//		}else{
			if ((s-1)*nf < 3){
//				cout << "C-NF-S: " << key << " does not exist in probability table...\n";
//				cout << "++++++++++++++keeping component which would be removed by previous approeach.." << endl;
//				cout << "previous_value: " << (s - 1)*nf << endl;
				//TODO: this part is added temporarily.. in reality all cases should exist in the probability table..
				clear_vertex(*vi, g_f);
			}
		//}
//		if ((component_ID_2_size[vertex_ID_2_component_ID[*vi]] - 1)*component_ID_2_clone_component_frequecy[vertex_ID_2_component_ID[*vi]] < p_value) {
//			//cout << "vertex-cleared .." << v_name[*vi];
//			clear_vertex(*vi, g_f);
//		}
	}

	
	
	//Printing MFG after removing spurious components...
	if (DEBUG > 3){
		stringstream ssGraphvizOutputFileName_temp2;
		ssGraphvizOutputFileName_temp2 << outputFile;
		
		string graphizOutputFileName_temp2 = ssGraphvizOutputFileName_temp2.str() + "-" + con_iter->first + "_spur_comp_rmvd.dot";
		ofstream graphvizOutFile_temp2(graphizOutputFileName_temp2.c_str(), ios::out);
		if (!graphvizOutFile_temp2) {
			cerr << "cannot open output file \n";
			exit(1);
		}
		
		write_graphviz(graphvizOutFile_temp2, g_f, custom_label_writer<U_VertexName>(v_name), make_label_writer<U_EdgeWeight>(e_weight), graph_writer());
		graphvizOutFile_temp2.close();
	}
	//*/
	//END: Disabling removal of spurious components...

	
	
	
	
	
	/**********************************************************************/
	
	
	//Computing Clones That are Buried into Multiple Clones...
	

	
	
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
	
	if (DEBUG > 3){
		if (candidate_buried_clones.size() > 0){
			cout << "Candidate buried clones: " << endl;
			for (unsigned int i = 0; i < candidate_buried_clones.size(); i++){
				cout << candidate_buried_clones[i] << endl;
			}
			cout << "----------" << endl;
		}
	}

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
	ostream_iterator<string> output_temp( cout, " ");
	if (DEBUG > 4){
		if (candidate_clone_set.size() > 0){
			cout << "Candidate clone set: " << endl;
			for (unsigned int i = 0; i < candidate_clone_set.size(); i++){
				copy( candidate_clone_set[i].begin(), candidate_clone_set[i].end(), output_temp);
				cout << "\n";
			}
		}
	}
	
	//compute the minimum hitting set for the candidate_clone_set..candidate_clone_set is MODIFIED BY compute_minimum_hitting_set..
	set<string> minimum_hitting_set;
	if (candidate_clone_set.size() > 0){
		compute_minimum_hitting_set(candidate_clone_set, minimum_hitting_set);
	}
	
	if (DEBUG > 3){
		if (minimum_hitting_set.size() > 0){
			cout << "Minimum hitting set...\n";
			copy(minimum_hitting_set.begin(), minimum_hitting_set.end(), output_temp);
			cout << "\n";
		}
	}
	
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
	
	
	//FOR DEBUGGING PURPOSES ONLY.. -- BEGIN
	if (DEBUG > 3){
		stringstream ssGraphvizOutputFileName_temp2;
		ssGraphvizOutputFileName_temp2 << outputFile;
	
		string graphizOutputFileName_temp2 = ssGraphvizOutputFileName_temp2.str() + "-" + con_iter->first + "_buried_clones_rmvd.dot";
		ofstream graphvizOutFile_temp2(graphizOutputFileName_temp2.c_str(), ios::out);
		if (!graphvizOutFile_temp2) {
			cerr << "cannot open output file \n";
			exit(1);
		}
	
		write_graphviz(graphvizOutFile_temp2, g_f, custom_label_writer<U_VertexName>(v_name), make_label_writer<U_EdgeWeight>(e_weight), graph_writer());
		graphvizOutFile_temp2.close();
	}
	//FOR DEBUGGING PURPOSES ONLY.. -- END	
	 

	
	//**********************************************************/		
	
	
	
	
	//Adding Mandatory Clones to the MFG...
	
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

	//compute total length of matached fragments..
	map<string, int, strCmp> clone_2_total_matched_length;
	for (std::vector<int>::size_type i = 0; i < vertex_ID_2_component_ID.size(); i++) {
		int fragment_size = get_fragment_size(v_name[i]);
		string clone_name = get_clone_name(v_name[i]);
		if (component_ID_2_size[vertex_ID_2_component_ID[i]] >= 2) {
			if (clone_2_total_matched_length.find(clone_name) == clone_2_total_matched_length.end()) {
				clone_2_total_matched_length[clone_name] = fragment_size;
			} else {
				clone_2_total_matched_length[clone_name] += fragment_size;
			}
		}
	}
	
	//Check for # unmatched fragments for the clones in the graph..
	map<string, int, strCmp> clone_2_total_matched;//clone_name to total number of matched bands...(# times clone occurs in a component of size >= 2)
	for (std::vector<int>::size_type i = 0; i < vertex_ID_2_component_ID.size(); i++) {
		string clone_name = get_clone_name(v_name[i]);
		if (component_ID_2_size[vertex_ID_2_component_ID[i]] >= 2) {
			if (clone_2_total_matched.find(clone_name) == clone_2_total_matched.end()) {
				clone_2_total_matched[clone_name] = 1;
			} else {
				clone_2_total_matched[clone_name]++;
			}
		}
	}

	for (set<Clone *, Contig::cloneCmp>::iterator iter = clones.begin(); iter != clones.end(); iter++) {
		if ((*iter)->isBuried()) {continue;}
		vector <Band *> fragments = (*iter)->getFragments();
		int n_total_bands = fragments.size();
		int total_band_length = 0;
		for (unsigned int idx = 0; idx < fragments.size(); idx++){
			total_band_length += fragments[idx]->getSize();
		}
		int total_matched = 0;
		if (clone_2_total_matched.find((*iter)->getName()) != clone_2_total_matched.end()) {
			total_matched = clone_2_total_matched[(*iter)->getName()];
		}
		
		int total_matched_length = 0;
		if (clone_2_total_matched_length.find((*iter)->getName()) != clone_2_total_matched_length.end()){
			total_matched_length = clone_2_total_matched_length[(*iter)->getName()];
		}
		double unmatched_vs_matched_ratio = (double)(total_band_length - total_matched_length)/(double)total_band_length*100;
		
		int n_unmatched_bands = n_total_bands - total_matched;
		if (DEBUG > 2) {
			cout << (*iter)->getName() << " has " << total_matched << " matching and " << n_unmatched_bands << " unmatched bands.." << endl;
		}
		//if ((double)n_unmatched_bands / (double)n_total_bands*100 >= EXTRA_FRAGMENT_PERCENTAGE) {
		if (unmatched_vs_matched_ratio >= EXTRA_FRAGMENT_PERCENTAGE){
//				cout << "************Not mandatory..." << (*iter)->getName() << endl;
//			}
			int N = 2;
			//go over the fragments of this clone and select two unmatched fragments (it can be singleton in the graph or not in the graph at all...)
			U_Vertex v[N];//stores the unmatched vertices...
			int n = 0;//keep track of the number of unmatched vertices picked.. exit the for loop when n = 2.. (we need two vertices...)
			bool inserted;
			NameUVertexMap::iterator pos;
			for (unsigned int i = 0; i < fragments.size() && n < N; i++) {
				pos = name_2_vertex.find(fragments[i]->getName());//check if there is a vertex for this fragment..
				if (pos == name_2_vertex.end()) {//if this fragment is not represented in the graph, add it to the graph...
					tie(pos, inserted) = name_2_vertex.insert(std::make_pair(fragments[i]->getName(), U_Vertex()));
					if (inserted) {
						v[n] = add_vertex(g_f);
						v_name[v[n]] = fragments[i]->getName();
						pos->second = v[n];
						n++;
					} else {
						cerr << "node with name " << fragments[i]->getName() << "already inserted!! Not possible..." << endl;
						exit(1);
					}
				} else if (degree(pos->second, g_f) == 0) {//if in the graph but a singleton then store this vertex too...
					v[n] = pos->second;
					n++;
				}
			}

			U_Edge e;
			bool b;
			for (int i = 0; i < N-1; i++) {
				tie (e, b) = add_edge(v[i], v[i+1], g_f);
				if (!b) {
					cerr << "Cannot add edge between " << v_name[v[0]] << " and " << v_name[v[1]] << endl;
					exit(0);
				}
			}
		}
		//else if (unmatched_vs_matched_ratio >= EXTRA_FRAGMENT_PERCENTAGE){
//			cout << "!!!!! Mandatory... " << (*iter)->getName() << endl;
//		}
	}

	

	
	
	
	if (DEBUG) {
		string graphizOutputFileName_3 = ssGraphvizOutputFileName.str() + "-" + con_iter->first + "_pruned.dot";
		ofstream graphvizOutFile_3(graphizOutputFileName_3.c_str(), ios::out);
		if (!graphvizOutFile_3) {
			cerr << "cannot open output file \n";
			exit(1);
		}
		write_graphviz(graphvizOutFile_3, g_f, custom_label_writer<U_VertexName>(v_name), make_label_writer<U_EdgeWeight>(e_weight), graph_writer());
		graphvizOutFile_3.close();
	}
	
	


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

	//if this is not the very first iteration and the threshold is high enough then check if max. component size is 2 or not...
//	if (iteration_ID && sulston_score_threshold >= HIGH_THRESHOLD) {
//		int max_component_size = get_max_component_size(component_ID_2_size, vertex_ID_2_component_ID, v_name);
//		if (max_component_size <= 2) {
//			//select all clones as mtp clones....
//			cout << "max component size is 2 for contig " << con_iter->first << endl;
//			ctgMTPFile << con_iter->first << "\t";
//			int nclones = 0;
//			for (set<Clone *, Contig::cloneCmp>::iterator iter = clones.begin(); iter != clones.end(); iter++) {
//				if ((*iter)->isBuried()) {continue;}
//				nclones++;
//			}
//			ctgMTPFile << nclones << "\t";
//			for (set<Clone *, Contig::cloneCmp>::iterator iter = clones.begin(); iter != clones.end(); iter++) {
//				if ((*iter)->isBuried()) {continue;}
//				ctgMTPFile << (*iter)->getName() << " ";
//			}
//			ctgMTPFile << "\n";
//			//now skip this contig...
//			E.clear();
//			E_f.clear();
//			continue;
//		}
//	}

	//remove redundant components...
	set<int> redundant_component_ID;
	compute_redundant_components(vertex_ID_2_component_ID, v_name, redundant_component_ID);

	/*cout << "components to be removed" << endl;
	 for (set<int>::iterator t_iter = redundant_component_ID.begin(); t_iter != redundant_component_ID.end(); t_iter++){
	 cout << *t_iter << endl;
	 }*/

	//remove vertices in the redundant components...
	//graph_traits< U_Graph >::vertex_iterator vi, vi_end;
	for (tie(vi, vi_end) = vertices(g_f); vi != vi_end; vi++) {
		if (redundant_component_ID.find(vertex_ID_2_component_ID[*vi]) != redundant_component_ID.end()) {
			clear_vertex(*vi, g_f);
		}
	}

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

	vector<string> vertices_to_remove;
	//Remove singletons from the graph.. For display purposes only....
	for (std::vector<int>::size_type i = 0; i != vertex_ID_2_component_ID.size(); ++i) {
		if (component_ID_2_size[vertex_ID_2_component_ID[i]] == 1) {//if the vertex in a component of size 1
			vertices_to_remove.push_back(v_name[i]);
		}
	}

	for (unsigned int i = 0; i < vertices_to_remove.size(); i++) {
		remove_vertex(name_2_vertex[vertices_to_remove[i]] - i, g_f);
	}

	if (DEBUG){
		//write the pruned graph into a file..
		string graphizOutputFileName_2 = ssGraphvizOutputFileName.str() + "-" + con_iter->first + "_final.dot";
		ofstream graphvizOutFile_2(graphizOutputFileName_2.c_str(), ios::out);
		if (!graphvizOutFile_2) {
			cerr << "cannot open output file \n";
			exit(1);
		}
		write_graphviz(graphvizOutFile_2, g_f, custom_label_writer<U_VertexName>(v_name), make_label_writer<U_EdgeWeight>(e_weight), graph_writer());
		graphvizOutFile_2.close();
	}
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


	//E_f for this contig is ready.. The edges from clones to fragments are also available..Now we need to generate the linear model...
	//output the current edge...
	stringstream ssOutputFileName;
	ssOutputFileName << outputFile;
	string outputFileName = ssOutputFileName.str() + "-" + con_iter->first + ".txt";
	ofstream outFile(outputFileName.c_str(), ios::out);
	if (!outFile) {
		cerr << "cannot open output file \n";
		exit(1);
	}

	outFile << "/*Sets*/\n";
	outFile << "set C; #clones..\n";
	outFile << "set F; #fragments..\n\n\n";
	outFile << "/*Decision Variables*/\n";
	outFile << "#clone c selected (1) or not (0)..\n";
	outFile << "var x {i in C}, integer, >=0, <= 1;\n\n\n";
	outFile << "/*Parameters*/\n";
	//outFile << "#fragment size\n";
	//outFile << "param s{i in F}, integer, >=0;\n";
	outFile << "#clone-fragment edge\n";
	outFile << "param e{i in C, j in F}, integer, >=0, <=1;\n\n\n";
	outFile << "#fragment-fragment edge\n";
	outFile << "param e_f{i in F, j in F}, integer, >=0, <=1;\n\n\n";

	outFile << "/*Objective function*/\n";
	outFile << "minimize y: sum{c in C} x[c];\n\n\n";
	//This line was added for testing purposes, but now disabled (093010) outFile << "minimize y: sum{c in C} sum{ff in F} x[c]*e[c,ff];\n\n\n";
	
	outFile << "/* Constraints */\n";
	outFile << "s.t. EmptySet: sum{c in C} x[c] >= 1;\n";
	outFile << "s.t. FragSet {f in F}: sum{c in C} x[c]*e[c,f] + sum{c in C} sum{ff in F} x[c]*e[c,ff]*(e_f[ff,f]+e_f[f,ff]) >= 1;\n\n\n";

	outFile << "data;\n\n\n";

	outFile << "set C := ";
	//print out all clones...
	list<string> all_clonename;
	list<string> all_fname;
	for (map<pair<string, string>, int, pairCmp>::iterator iter = E.begin(); iter != E.end(); iter++) {
		all_clonename.push_back(iter->first.first);
	}

	//all fragment nodes remained in g_f will be stores in all_fname vector..
	for (tie(vi, vi_end) = vertices(g_f); vi != vi_end; vi++) {
		if (degree(*vi, g_f) > 0) {
			all_fname.push_back(v_name[*vi]);
		}
	}

	all_clonename.sort();
	all_clonename.unique();
	all_fname.sort();
	all_fname.unique();
	for (list<string>::iterator iter = all_clonename.begin(); iter != all_clonename.end(); iter++) {
		outFile << "\"" << *iter << "\"" << " ";
	}
	outFile << ";\n";

	//cout << "set C: is printed..\n";

	outFile << "set F := ";
	for (list<string>::iterator iter = all_fname.begin(); iter != all_fname.end(); iter++) {
		outFile << *iter << " ";
	}
	outFile << ";\n";

	//cout << "set F is printed\n";
	outFile << "param e: ";
	for (list<string>::iterator iter = all_fname.begin(); iter != all_fname.end(); iter++) {
		outFile << *iter << " ";
	}

	outFile << " :=\n";
	for (list<string>::iterator iter = all_clonename.begin(); iter != all_clonename.end(); iter++) {
		outFile << "           "<< "\"" << *iter << "\"" << " ";//print clonename
		for (list<string>::iterator iter2 = all_fname.begin(); iter2 != all_fname.end(); iter2++) {
			pair<string, string> p(*iter, *iter2);
			if (E.find(p) != E.end()) {
				outFile << E[p] << " ";
			} else {
				outFile << "0 ";
			}
		}
		//cout << nTotal << endl;
		outFile << "\n";
	}
	outFile << ";\n";
	outFile << "param e_f: ";
	for (list<string>::iterator iter = all_fname.begin(); iter != all_fname.end(); iter++) {
		outFile << *iter << " ";
	}
	outFile << " :=\n";

	for (list<string>::iterator iter = all_fname.begin(); iter != all_fname.end(); ++iter) {
		outFile << "           " << *iter << " ";//print frag. name

		//int nTotal = 0;
		for (list<string>::iterator iter2 = all_fname.begin(); iter2 != all_fname.end(); ++iter2) {
			if ((name_2_vertex.find(*iter) == name_2_vertex.end()) || (name_2_vertex.find(*iter2) == name_2_vertex.end())) {
				cerr << "cannot find names " << *iter << " or " << *iter2 << endl;
				exit(1);
			}
			if (vertex_ID_2_component_ID[name_2_vertex[*iter]] == vertex_ID_2_component_ID[name_2_vertex[*iter2]]) {
				outFile << "1 ";
			} else {
				outFile << "0 ";
			}
		}
		outFile << "\n";
	}
	outFile << ";\n";
	outFile << "end;\n";

	outFile.close();

	//clear E & E_f
	E.clear();
	E_f.clear();
}
ctgMTPFile.close();

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
		for (set<string>::iterator s_iter = clone_names.begin(); s_iter
				!= clone_names.end(); s_iter++) {
			concatenated_clone_name += *s_iter + "-";
		}
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
//		cout << "flow is " << flow << endl; 
		if (flow < min_flow){
			min_flow = flow;
			min_cut.clear();
			for (unsigned int i = 0; i < color.size(); i++) {
				min_cut[new_to_old_vertex_map[i]] = color[i];
			}
		}
	}
	//cout << "min_flow: " << min_flow << endl;
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
			//cout << "overlap_score_exponent[" << p.first << ", " << p.second << "]: " << overlap_score_exponent[p] << endl; 
			
			if (overlap_score_exponent[p] <= MTP_ILP_VERY_HIGH_THRESHOLD){
				return true;
			}
		}
	}
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
