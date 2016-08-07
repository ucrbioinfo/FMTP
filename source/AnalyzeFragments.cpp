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

#include "Contig.h"

int WINDOW_SIZE = 0;
int BURY_THRESHOLD = 70;

void usage()
{
	cout << "AnalzeFragments: Generates a linear programming model to select MTP clones of contigs in a FPC map..\n";
	cout << "-c ctg_BAC_clone file\n"; 
	cout << "-s size file\n";
	cout << "-g genome equivalence (default 6)\n";
	cout << "-p phase (1 or 2) (default 1).. \n";
	cout << "-w window size \n";
	cout << "-b bury threshold (default 70)\n";
	cout << "-o output_file base file (without extension).. Linear Programming Model file (for GLPK, in MathProg language)\n";
	exit(0);
}

struct strCmp {
    bool operator()( const string s1, const string s2 ) const {
      return s1.compare(s2) < 0;
    }
};

struct pairCmp {
    bool operator()( const pair<string, string> p1, const pair<string, string> p2 ) const {
      if (p1.first.compare(p2.first) == 0 ){
        return p1.second.compare(p2.second) < 0;
      }else{
        return p1.first.compare(p2.first) < 0;
      }
    }
};

bool isNumeric (const char * s){
    char * p;
    strtol(s, &p, 10);
    return !*p;
}

void discretize(int & n){
	int y = n % WINDOW_SIZE;
	if (y > ceil(WINDOW_SIZE/2)){
		n += WINDOW_SIZE - y;
	}else{
		n -= y;
	}
}


int main(int argc, char**argv){
	
	map<string, Clone *, strCmp> cloneObjectByName; // stores the pointer to the Clone object of each clone (by clone name)
	map<string, Contig *, strCmp> contigObjectByName; // stores the pointer to the Contig object of each contig (by contig name)
	
	char ctgBACFile[FILENAME_MAX];
	char sizeFile[FILENAME_MAX];
	char outputFile[FILENAME_MAX];
	
	int ctgBACFilePresent = 0, sizeFilePresent = 0;
	int genome_equivalence = 6;
	int phase = 1;
	
	map<string, vector<int>, strCmp > fragmentDegree; //key: fragment size-contigID.. value: degree of copies of nodes for this fragment size.. (ex. f_size 700 may have 3 copies.. in this case, the vector will have three elements.. the degree of 700-1 is the first element of the vector, ...
	while (1)    {
      char c = getopt(argc, argv, "vhf:c:s:o:g:p:w:b:");

      if (c == -1)
        break;

      switch (c)
        {
        case 'v':
          break;

        case 'h':
          usage();
          break;

		case 'c':
			strcpy(ctgBACFile, optarg);
			printf("Using ctgBACFile in %s\n", ctgBACFile);
			ctgBACFilePresent = 1;
			break;
        
        case 's':
          strcpy(sizeFile,optarg);
          printf("Using size file %s\n",sizeFile);
          sizeFilePresent = 1;
          break;
        
        case 'g':
        	genome_equivalence = atoi(optarg);
        	printf("Genome equivalence : %d\n", genome_equivalence);
        	break;
        case 'p':
        	phase = atoi(optarg);
        	printf("Phase : %d\n", phase);
        	break;
        case 'w':
        	WINDOW_SIZE = atoi(optarg);
        	printf("WINDOW_SIZE : %d\n", WINDOW_SIZE);
        	break;		
        case 'b':
        	BURY_THRESHOLD = atoi(optarg);
        	printf("BURY_THRESHOLD : %d\n", BURY_THRESHOLD);
        	break;
        case 'o':
          strcpy(outputFile,optarg);
          printf("Will write output data to '%s'\n",outputFile);
          break;
          
        default:
          printf ("getopt() returned character code 0%o - not handled in code. Useage instructions follow.\n", c);
          usage();
        }
    }
    
    if (!sizeFilePresent || !ctgBACFilePresent){
    	cerr << "input file missing\n";
    	usage();
    } 
    
    if (strlen(outputFile) == 0){
    	cerr << "output file missing\n";
    	usage();
    }
    
    if (genome_equivalence <=0){
    	cerr << "invalid genome equivalence\n";
    	usage();
    }
    
    if (WINDOW_SIZE == 0){
    	cerr << "invalid window size..\n";
    	usage();
    }
    
    //read the input file...
    ifstream inCtgBACFile(ctgBACFile, ios::in);
    
    if (!inCtgBACFile){
    	cerr << "cannot open input file\n";
    	exit(1);
    }
    
    
    //read the ctg-BAC file and generate clone & contig objects
    char line[90000];
	char cloneName[20];
	char contigName[10];
	int nClones;
	
	while (inCtgBACFile.getline(line, 90000)){
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
		
		while (tok2 != NULL)  {
			strcpy(cloneName, tok2);
			//cout << "clonename: " << cloneName << endl;
			if (strlen(cloneName) == 0){
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
	if (!inSizeFile){
		cerr << "Cannot read the size file, " << sizeFile << endl;
		exit(1);
	} 
	
	multiset<int> tmpFragments;
	bool skip = false;//used for ignoring clones that are not in cloneObjectByName map....
	while (inSizeFile.getline(line, 500)){
		char * tok;
		tok = strtok(line, "\t");
		//not need the rest of the tokens, if there is any..
		
		//if the value is a band size of an skipped clone then ignore..
		if (isNumeric(tok) && !skip){
			int n = atoi(tok);
			if (n == -1){
				//store this clone..
				Clone * clonePtr = cloneObjectByName[cloneName];
				
				for (multiset<int>::iterator iter = tmpFragments.begin(); iter != tmpFragments.end(); iter++){
					clonePtr->addToFragments(*iter);	
				}
				tmpFragments.clear();
			}else{
				//nTotalBands++;
				discretize(n);
				tmpFragments.insert(n);
			}
		}else{
			strcpy(cloneName, line);
			/*ignore this clone if it is not in cloneObjectByName map... 
			(because size file may contain clones whose coordinates are not available)..
			*/
			if (cloneObjectByName.find(cloneName) == cloneObjectByName.end()){
				skip = true;//all band sizes of this clone will be skipped..
			}else{
				skip = false;
			}
		}
	}
	
	cout << "size file is read\n";
	
	//remove buried clones...
	for (map<string, Contig *, strCmp>::iterator con_iter = contigObjectByName.begin(); con_iter != contigObjectByName.end(); con_iter++){
		set<Clone *, Contig::cloneCmp> clones = con_iter->second->getClones();
		for (set<Clone *, Contig::cloneCmp>::iterator cl_iter_1 = clones.begin(); cl_iter_1 != clones.end(); cl_iter_1++){
			multiset<int> f_clone_1 = (*cl_iter_1)->getFragments();
			for (set<Clone *, Contig::cloneCmp>::iterator cl_iter_2 = clones.upper_bound(*cl_iter_1); cl_iter_2 != clones.end(); cl_iter_2++){
				multiset<int> f_clone_2 = (*cl_iter_2)->getFragments();
				list<int> intersection;
				set_intersection(f_clone_1.begin(), f_clone_1.end(), f_clone_2.begin(), f_clone_2.end(), front_inserter(intersection));
				//cout << (*cl_iter_1)->getName() << "\t" << (*cl_iter_2)->getName() << "\t" << intersection.size() << endl;
				if (f_clone_1.size() > f_clone_2.size()){ //try to bury clone_2
					
					if ((double)intersection.size()/(double)f_clone_2.size()*100 >= BURY_THRESHOLD){
						//cout << (*cl_iter_2)->getName() << "\t" << (*cl_iter_1)->getName() << endl;
						(*cl_iter_2)->setParent(*cl_iter_1);
					}
				}else{//try to bury clone_1
					if ((double)intersection.size()/(double)f_clone_1.size()*100 >= BURY_THRESHOLD){
						//cout << (*cl_iter_1)->getName() << "\t" << (*cl_iter_2)->getName() << endl;
						(*cl_iter_1)->setParent(*cl_iter_2);
					}
				}
			}
		}
	}
	
	//generate fragments for each contig..
	//for each contig, visit each clone in this contig..
	
	for (map<string, Contig *, strCmp>::iterator con_iter = contigObjectByName.begin(); con_iter != contigObjectByName.end(); con_iter++){
		set<Clone *, Contig::cloneCmp> clones = con_iter->second->getClones();
		//map f_contig.. 
		 // key: fragment size f.. value: pair<total m(f) of clones that contain f, max(m(f)) for all clones in this contig>..
		 // check 06-20-07-3.1 for more details of this map.. Basically, the idea is to generate F_contig that contains all unique fragments in clones in this contig with their multiplicity...
		//
		map<unsigned int, pair<unsigned int, unsigned int>, less<int> > f_contig;
		//go over all clones in this contig...
		for (set<Clone *, Contig::cloneCmp>::iterator cl_iter = clones.begin(); cl_iter != clones.end(); cl_iter++){
			//skip if it's a buried clone..
			if ((*cl_iter)->isBuried()){continue;}
			
			//go over all fragments of *s_iter
			multiset<int> f_clone = (*cl_iter)->getFragments();
			list<int> f_clone_unique(f_clone.begin(), f_clone.end());
			f_clone_unique.sort();
			f_clone_unique.unique();
			for (list<int>::iterator f_iter = f_clone_unique.begin(); f_iter != f_clone_unique.end(); f_iter++){
				//cout << *f_iter << endl;
				//insert fragment *f_iter to f_contig.. if it is already inserted update it..
				if (f_contig.find(*f_iter) != f_contig.end()){
					f_contig[*f_iter].first += f_clone.count(*f_iter);
					//if this clone contains more copies of this fragment then update f_contig[frag_size].second...
					if (f_clone.count(*f_iter) > f_contig[*f_iter].second){
						f_contig[*f_iter].second = f_clone.count(*f_iter);
					}
				}else{
					pair<int, int> p;
					p.first = f_clone.count(*f_iter);//nclones contain this fragment
					p.second = f_clone.count(*f_iter); //number of copies of this fragment in this clone...
					f_contig[*f_iter] = p;
				}
			}
		}
		
		//now f_contig is ready...
		//go over each element of f_contig and store fragments of this contig...
		for (map<unsigned int, pair<unsigned int, unsigned int>, less<int> >::iterator iter = f_contig.begin(); iter != f_contig.end(); iter++){
			int size = iter->first;
			//cout << size << "\t" << iter->second.first << "\t" << iter->second.second << endl;
			int m = max(ceil((double)iter->second.first/(double)genome_equivalence), (double)iter->second.second);
			//cout << size << "\t" << m << endl;
			//f.m = max(iter->second.first/GENOME_EQUI, iter->second.second);
			//f.m = max(ceil(2/3), (double)4);
			
			for (int i = 0; i < m; i++){
				con_iter->second->addToFragments(size);
			}
			vector<int> degree(m, 0);//degree of m copies of fragment size...
			
			stringstream s1;
			s1 << size;
						
			string f_name;
			f_name = s1.str() + "-" + con_iter->first;
			
			fragmentDegree[f_name] = degree;
		}
	}
	
	cout << "Fragments for each contig are generated\n";
	
	//print all contigs for debugging...
	/*for (map<string, Contig *, strCmp>::iterator con_iter = contigObjectByName.begin(); con_iter != contigObjectByName.end(); con_iter++){
		//con_iter->second->print();
		//cout << "Number of fragments :" << con_iter->second->getFragments().size() << endl;
	}*/
	
	
	//generate clone-fragment association....
	
	
	for (map<string, Contig *, strCmp>::iterator con_iter = contigObjectByName.begin(); con_iter != contigObjectByName.end(); con_iter++){
		map<pair<string, string>, int, pairCmp > edge;//<clone_name, fragment_name> all clone-fragment pairs in which clone contains the fragment....
		set<Clone *, Contig::cloneCmp> clones = con_iter->second->getClones();
		
		//go over all clones in this contig...
		for (set<Clone *, Contig::cloneCmp>::iterator cl_iter = clones.begin(); cl_iter != clones.end(); cl_iter++){
			//skip if it's a buried clone...
			if ((*cl_iter)->isBuried()){continue;}
			//go over all fragments of *s_iter
			multiset<int> f_clone = (*cl_iter)->getFragments();
			list<int> f_clone_unique(f_clone.begin(), f_clone.end());
			f_clone_unique.unique();
			for (list<int>::iterator f_iter = f_clone_unique.begin(); f_iter != f_clone_unique.end(); f_iter++){
				string f_name;
				std::stringstream s1;
				s1 << *f_iter;
								
				f_name = s1.str() + "-" + con_iter->first; //fragment name: Fragment_size-id-contigID... ex. 750-1-15.. 
				
				int m = f_clone.count(*f_iter);
				
				//for each copy of the fragment selected a unique copy of the fragment...
				vector<int> id_of_previous_selected;
				
				for (int i = 0; i < m; i++){
					//check a copy of f_name starting from id_of_previous_selected + 1 st item in fragmentDegree[f_name] with degree < genome_equivalence..
					bool fragmentSelected = false;
					
					
					int min_index = -1, min_degree = -1;
					for (unsigned int j = 0; j < fragmentDegree[f_name].size(); j++){
						if (min_degree == -1){
							min_index = j;
							min_degree= fragmentDegree[f_name][j];
						}else if(fragmentDegree[f_name][j] < min_degree){
							//ensure that j has not selected before..
							bool selected_before = false;
							for (unsigned int k = 0; k < id_of_previous_selected.size(); k++){
								if (id_of_previous_selected[k] == (int)j){
									selected_before = true;
									break;
								}
							}
							if (!selected_before){
								min_index = j;
								min_degree= fragmentDegree[f_name][j];
							}
						}
					}
					//min_index is found...
					
					//cout << "fragment Degree " << fragmentDegree[f_name][min_index] << " for " << f_name << ". gen_equ: " << genome_equivalence << endl;
					if (fragmentDegree[f_name][min_index] < genome_equivalence){//do not occupy a node completely in the first phase...
						//cout << "entered" << endl;
						//a fragment is selected..
						fragmentSelected = true;
						//update id_of_previous_selected not to select the same fragment node in this clone..
						id_of_previous_selected.push_back(min_index);
						//update the degree of the selected node..
						fragmentDegree[f_name][min_index]++;
						//add an edge between clonename & f_name-j..
						stringstream s;
						s << min_index;
						string final_f_name = f_name + "-" + s.str();
						
						
						pair<string, string> p;
						p.first = (*cl_iter)->getName();
						p.second = final_f_name;
						//cout << p.first << "\t" << p.second << endl;
						edge[p] = 1;
					}
					
					if (!fragmentSelected){
						//cout << "fragment Degree " << fragmentDegree[f_name][min_index] << " for " << f_name << endl;
						cerr << "No node is found for " << *f_iter << " for clonename " << (*cl_iter)->getName() << " in contig " << con_iter->first << endl;
						exit(1);
					}
				}
			}
		}
		//cout << "clone-fragment edges are generated\n";
		//update edge...
		for (map<pair<string, string>, int, pairCmp >::iterator iter = edge.begin(); iter != edge.end(); iter++){
			char *tok;
			string f_name;
			char f_f_name[15];
			//cout << "1" << endl;
			strcpy(f_f_name, iter->first.second.c_str());
			//cout << f_f_name << endl;
			tok = strtok(f_f_name, "-");
			string s1(tok);
			
			tok = strtok(NULL, "-");
			string s2(tok);
			tok = strtok(NULL, "-");
			
			int copy_id = atoi(tok);
			f_name = s1 + "-" + s2;
			//cout << f_name << "\t" << copy_id << "\t" << fragmentDegree[f_name][copy_id] << endl;
			int nclones = contigObjectByName[s2]->getClones().size();
/*			if (phase == 1){
				if (fragmentDegree[f_name][copy_id] <= genome_equivalence/2 && fragmentDegree[f_name][copy_id] < nclones/2){
					iter->second = 0;
				}
			}else if (phase == 2){
				if (fragmentDegree[f_name][copy_id] <= genome_equivalence/5 && fragmentDegree[f_name][copy_id] < nclones/5){
					//cout << iter->first.first << " is removed.." << endl;
					iter->second = 0;
				}
			}
*/
			int threshold_1 = genome_equivalence/phase;
			int threshold_2 = nclones/phase;
			//cout << "th: " << threshold_1 << "\t" << threshold_2 << endl;
			if (fragmentDegree[f_name][copy_id] < threshold_1 && fragmentDegree[f_name][copy_id] < threshold_2){
					//cout << iter->first.first << " is removed.." << endl;
					iter->second = 0;
			}		
		}	
		
		//output the current edge...
		stringstream ssOutputFileName;
		ssOutputFileName << outputFile;
		string outputFileName = ssOutputFileName.str() + "-" + con_iter->first + ".txt";
		ofstream outFile(outputFileName.c_str(), ios::out);
		if (!outFile){
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
		outFile << "#fragment size\n";
		outFile << "param s{i in F}, integer, >=0;\n";
		outFile << "#clone-fragment edge\n";
		outFile << "param e{i in C, j in F}, integer, >=0, <=1;\n\n\n";
		
		outFile << "/*Objective function*/\n";
		outFile << "minimize y: sum{c in C} x[c];\n\n\n";
		
		outFile << "/* Constraints */\n";
		outFile << "s.t. EmptySet: sum{c in C} x[c] >= 1;\n";
		outFile << "s.t. FragSet {f in F}: sum{c in C} x[c]*e[c,f] >= 1;\n\n\n";
		
		outFile << "data;\n\n\n";
	
		outFile << "set C := ";
		//print out all clones...
		list<string> all_clonename;
		list<string> all_fname;
		for (map<pair<string, string>, int, pairCmp>::iterator iter = edge.begin(); iter != edge.end(); iter++){
			all_clonename.push_back(iter->first.first);
			if (iter->second){
				all_fname.push_back(iter->first.second);
			}
		}
		all_clonename.sort();
		all_clonename.unique();
		all_fname.sort();
		all_fname.unique();
		for (list<string>::iterator iter = all_clonename.begin(); iter != all_clonename.end(); iter++){
			outFile << *iter << " ";
		}
		outFile << ";\n";
		
		//cout << "set C: is printed..\n";
		
		outFile << "set F := ";
		for (list<string>::iterator iter = all_fname.begin(); iter != all_fname.end(); iter++){
			outFile << *iter << " ";
		}
		outFile << ";\n";
		
		//cout << "set F is printed\n";
		 
		outFile << "param e: ";
		for (list<string>::iterator iter = all_fname.begin(); iter != all_fname.end(); iter++){
			outFile << *iter << " ";
		}
		
		outFile << " :=\n";
		for (list<string>::iterator iter = all_clonename.begin(); iter != all_clonename.end(); iter++){
			outFile << "           "<< *iter << " ";//print clonename
			//cout << *iter << "\t";
			//int nTotal = 0;
			for (list<string>::iterator iter2 = all_fname.begin(); iter2 != all_fname.end(); iter2++){
				pair<string, string> p(*iter, *iter2);
				if (edge.find(p) != edge.end()){
					//if (edge[p]){nTotal++;}
					outFile << edge[p] << " ";
				}else{
					outFile << "0 ";
				}
			}
			//cout << nTotal << endl;
			outFile << "\n";
		}
		outFile << ";\n";
		outFile << "end;\n";
		
		
		outFile.close();	
		
	}
	
	
	
	
	
	//print fragmentDegree for debugging purposes..
	/*for (map<string, vector<int>, strCmp>::iterator iter = fragmentDegree.begin(); iter != fragmentDegree.end(); iter++){
		for (unsigned int i = 0; i < iter->second.size(); i++){
			cout << iter->first << "\t" << iter->second[i] << endl;
		}
	}*/
	
	//edge vector is ready.. we can output the linear programming model....
	
	cout << "done\n";
	return 0;
		
}
