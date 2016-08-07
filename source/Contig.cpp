#include "Contig.h"

Contig::Contig(const char * s){
	setName(s);
}

Contig::~Contig()
{
}

void Contig::addToClones(Clone * _c){
	clones.insert(_c);
}

void Contig::print(){
	//print contig name & clones & fragments..
	cout << "Name: " << name << endl;
	cout << "Clones\t" << clones.size() << "\n";
	for (set<Clone *, cloneCmp>::iterator iter = clones.begin(); iter != clones.end(); iter++){
		cout << (*iter)->getName() << endl;
		vector<Band *> fragments = (*iter)->getFragments();
		cout << "Fragments\t" << fragments.size() << "\n";
		for (vector<Band *>::iterator iter = fragments.begin(); iter != fragments.end(); iter++){
			cout << (*iter)->getSize() << endl;
		}		
	}
	
	
	cout << "--------------------\n\n";
}
