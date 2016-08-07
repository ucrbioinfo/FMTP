#ifndef CLONE_H_
#define CLONE_H_

#include "Utils.h"


#include <vector>
#include <string>
#include <set>

using namespace std;

class Clone{
	
public:
	Clone(const char * = "");
	virtual ~Clone();
	void setName(const char * _name);
	string getName() const {return name;}
	void addToFragments(Band * _f);
	int getNFragments() const {return fragments.size();}
	vector<Band *> getFragments() const {return fragments;}
	multiset<int> getFragmentSizes() const;
	
	void setPartiallyBuried(const bool _buried);
	void setParent(Clone * const _parent);
	Clone * getParent() const {return parent;}
	void unsetParent();
	bool isBuried() const {return buried;}
	
private:
	void setBuried(const bool _buried);

protected:
	string name;
	//set<Fragment, compFragment> fragments;
	bool buried;
	vector<Band *> fragments;
	Clone * parent;
};

#endif /*CLONE_H_*/
