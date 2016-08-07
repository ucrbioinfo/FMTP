#ifndef CONTIG_H_
#define CONTIG_H_

#include <iostream>
#include <set>
#include "Clone.h"



class Contig : public Clone
{
public:
	Contig(const char * _s);
	virtual ~Contig();
	struct cloneCmp {
	    bool operator()( const Clone * c1, const Clone * c2 ) const {
	      return c1->getName().compare(c2->getName()) < 0;
	    }
	};
	void addToClones(Clone * _c);
	set<Clone *, cloneCmp> getClones() const {return clones;}
	void print();	
private:
	set<Clone *, cloneCmp> clones;	
};

#endif /*CONTIG_H_*/
