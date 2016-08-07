#include <float.h>
#include <iostream>
#include <iomanip>
#include <map>
#include <algorithm>
#include "Utils.h"
#include <fstream> 

struct compStr {
    bool operator()( const string s1, const string s2 ) const {
      return s1.compare(s2) < 0;
    }
};

/*int getPoolID(char *_poolName){
	char *tok;
	tok = strtok(_poolName, "-");//ignore the first token..
	if (tok != NULL){
		tok = strtok(NULL, "-");
		return atoi(tok);
	}else{
		cerr << "pool_naming error.. pool-poolID is the correct format.." << endl;
		exit(1);
	}
	
}*/


int numberOfSharedBands(int t, const multiset<int> & bands1, const multiset<int> & bands2){
	int match = 0;
	multiset<int>::iterator lstart = bands2.begin();
	
	for (multiset<int>::iterator iter_1 = bands1.begin(); iter_1 != bands1.end(); iter_1++){
		for (multiset<int>::iterator iter_2 = lstart; iter_2 != bands2.end(); iter_2++){
			//cout << i << "\t" << j << "\t" << bands1[i] << "\t" << bands2[j] << endl;
			int diff = *iter_1 - *iter_2;
			if (abs(diff) <= t){
				match++;
				//lstart = j + 1.. start from the next band in the next iteration..
				lstart = iter_2;
				if (lstart != bands2.end()){	lstart++;}
				break;
			}else if  (diff < 0){
				lstart = iter_2;
				break;
			}
		}
	}
	
	return match;
} 

int lengthOfSharedBands(int t, const multiset<int> & bands1, const multiset<int> & bands2){
	int match = 0;
	int total_length = 0;
	multiset<int>::iterator lstart = bands2.begin();
	
	for (multiset<int>::iterator iter_1 = bands1.begin(); iter_1 != bands1.end(); iter_1++){
		for (multiset<int>::iterator iter_2 = lstart; iter_2 != bands2.end(); iter_2++){
			//cout << i << "\t" << j << "\t" << bands1[i] << "\t" << bands2[j] << endl;
			int diff = *iter_1 - *iter_2;
			if (abs(diff) <= t){
				match++;
				total_length += (*iter_1 + *iter_2)/2;
				//lstart = j + 1.. start from the next band in the next iteration..
				lstart = iter_2;
				if (lstart != bands2.end()){	lstart++;}
				break;
			}else if  (diff < 0){
				lstart = iter_2;
				break;
			}
		}
	}
	
	return total_length;
} 



void fastSulston (int clone1_nbands, int clone2_nbands, int match, int t, long gellen, double& prob){

    register int i, k, n, nbH, nbL;
    double a, c, pp, pa, pb;
    double lim, Logfact[NBANDS+1], psmn=0.0;

	Logfact[0] = 0.0;
	Logfact[1] = 0.0;
	for (i = 2; i <= NBANDS; ++i)
		Logfact[i] = Logfact[i - 1] + log((double) i);
	
	lim = log(DBL_MIN);


    if (t == 0)  psmn = 1 - (double) 1.0 / gellen;
    else  psmn = 1 - (double) (t << 1) / gellen;

    nbL = MiN(clone1_nbands,clone2_nbands);
    nbH = MaX(clone1_nbands,clone2_nbands);

    prob = DBL_MIN;

    pp = (double) pow(psmn, nbH);
    pa = log(pp);
    pb = log(1.0 - pp);

    /* 18june00 starting at nbl 64, nhl 74, if the match is zero,
            this optimization causes the prob to be below cutoff */

    if (match==0) prob =1.0;
    else {
      for (k=0, n = match; n <= nbL && k < FAST_SUL_ITERATIONS; n++, k++)
      {
            c = Logfact[nbL] - Logfact[nbL - n] - Logfact[n];
            a = c + n * pb + (nbL - n) * pa;
            prob +=  (a >= lim) ? exp(a) : exp(lim);
       }
    }

}



void assignIntValue(int *key, const int value){
	if (value != -1){
		*key = value;
	}
}

void assignDoubleValue(double *key, const double value){
	if (value != -1.0){
		*key = value;
	}
}
