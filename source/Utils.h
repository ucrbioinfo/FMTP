#ifndef UTILS_H_
#define UTILS_H_

#include <cmath>
#include <string>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/tuple/tuple.hpp>
#include <set>
#include "Band.h"


using namespace std;
using namespace boost;

#define NBANDS 2000
#define MiN(a,b) (a < b ? a : b)
#define MaX(a,b) (a < b ? b : a)
#define BIG_NUMBER 100000 //a large number...

#define DEBUG_KEYWORD "DEBUG"
#define TOL_KEYWORD "TOLERANCE"
#define GEL_LEN_KEYWORD "GEL_LENGTH"
#define FAST_SUL_ITERATIONS_KEYWORD "FAST_SULSTON_ITERATIONS"
#define COOPERATIVE_BURY_THRESHOLD_KEYWORD "BURY_INTO_MULTIPLE_CLONES_PERCENTAGE"
#define BURY_THRESHOLD_KEYWORD "BURY_INTO_ONE_CLONE_PERCENTAGE"
#define EXTRA_FRAGMENT_PERCENTAGE_KEYWORD "MANDATORY_CLONE_MIN_EXTRA_FRAGMENT_PERCENTAGE"
#define MTP_ILP_VERY_HIGH_THRESHOLD_KEYWORD "MTP_ILP_VERY_HIGH_THRESHOLD"
#define MTP_MST_VERY_HIGH_THRESHOLD_KEYWORD "MTP_MST_VERY_HIGH_THRESHOLD"


//#define TOL 7 //cowpead:5, barley: 3, rice: 7
//#define GEL_LENGTH 3300 //cowpea:36,000, rice: 3300, barley: 18,000 
#define FAST_SUL_ITERATIONS 3

#define HICF 1
#define AGAROSE 2

//default cutoff threshold values for various cases..
#define MTP_ILP_HICF_C 1e-25
#define MTP_ILP_AGAROSE_C 1e-10
#define MTP_MST_HICF_C 1e-5
#define MTP_MST_AGAROSE_C 1e-2

//default MTP_MST_VERY_HIGH_THRESHOLD (C') values for various cases...
#define MTP_MST_VERY_HIGH_THRESHOLD_HICF 3
#define MTP_MST_VERY_HIGH_THRESHOLD_AGAROSE 1

//default MTP_ILP_VERY_HIGH_THRESHOLD (C') values for various cases...
#define MTP_ILP_VERY_HIGH_THRESHOLD_HICF 5
#define MTP_ILP_VERY_HIGH_THRESHOLD_AGAROSE 3

//default Q values for various cases

#define MTP_ILP_HICF_Q 4
#define MTP_ILP_AGAROSE_Q 3
#define MTP_MST_HICF_Q 4
#define MTP_MST_AGAROSE_Q 3



//int TOL = 7; //7 for agarose, 3 for HICF
//int GEL_LENGTH = 3300;//3300 for agarose, 18000 for HICF
//int FAST_SUL_ITERATIONS = 3;
//int DEBUG = 0;
//double COOPERATIVE_BURY_THRESHOLD = 80; //used for burying clone into multiple clones..
//double BURY_THRESHOLD = 80;
//double EXTRA_FRAGMENT_PERCENTAGE = 40; //if extra fragments of a node (i.e. match nowhere in clone - fragment graph) is more than this percentage then they are added to clone-fragment graph...
//double MTP_ILP_VERY_HIGH_THRESHOLD = 2.5; //5 for barley, 2.5 for rice.. if -log(Sulston(x,y)) is lower than this value then x & y are definitely not overlapping...
//double MTP_MST_VERY_HIGH_THRESHOLD = 1; // 3 for barley, 1 for rice.. if a component contains a pair that has an overlap_score_exponent higher than this value then that component is splitted..

/*
#define LEFT -1
#define RIGHT 1

#define NPOOLS 2000

struct EdgeProperties {
  bool matched;
};
*/

/*
typedef adjacency_list_traits < vecS, vecS, directedS > Traits;
	typedef adjacency_list < listS, vecS, directedS,
		property < vertex_name_t, std::string >,
		property < edge_capacity_t, long,
		property < edge_residual_capacity_t, long,
		property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;

	
*/	



struct VertexProperties {
	Band * band_ptr;
};

typedef adjacency_list < vecS, vecS, bidirectionalS, property < vertex_name_t, std::string>, property < edge_name_t, long> > BFS_Graph;
typedef property_map < BFS_Graph, vertex_name_t >::type VertexName;
typedef property_map < BFS_Graph, edge_name_t >::type EdgeName;
typedef adjacency_list < vecS, vecS, bidirectionalS >::vertex_descriptor BFS_Vertex;
typedef adjacency_list < vecS, vecS, bidirectionalS >::edge_descriptor BFS_Edge;

typedef adjacency_list_traits < vecS, vecS, directedS > Traits;
typedef adjacency_list < listS, vecS, directedS,
		VertexProperties, 
		property < edge_capacity_t, double,
		property < edge_residual_capacity_t, double,
		property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;

//*/	
/*typedef adjacency_list < vecS, vecS, directedS,
    no_property, EdgeProperties, VertexProperties> Graph;	
*/

typedef Traits::vertex_descriptor Vertex;
typedef Traits::edge_descriptor Edge;

//typedef property_map < Graph, edge_capacity_t >::type Capacity;
//typedef property_map < Graph, edge_reverse_t >::type Reverse;			

	
	

//function prototypes...
int numberOfSharedBands(int, const multiset<int> &, const multiset<int> &);
int lengthOfSharedBands(int, const multiset<int> &, const multiset<int> &);
void fastSulston (int, int, int, int, long, double&);
void assignIntValue(int *, const int);
void assignDoubleValue(double *, const double);

#endif /*UTILS_H_*/
