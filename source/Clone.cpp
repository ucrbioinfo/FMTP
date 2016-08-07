#include "Clone.h"

Clone::Clone(const char * s){
	setName(s);
	setBuried(false);
}

Clone::~Clone()
{
}

void Clone::setName(const char * _name){
	name.assign(_name);
}

void Clone::addToFragments(Band * _f){
	fragments.push_back(_f);
}

multiset<int> Clone::getFragmentSizes() const{
	multiset<int> result;
	for(unsigned int i = 0; i < fragments.size(); i++){
		result.insert(fragments[i]->getSize());
	}
	return result;
}

void Clone::setPartiallyBuried(const bool _buried){
	setBuried(_buried);
}

void Clone::setBuried(const bool _buried){
	buried = _buried;
}

void Clone::setParent(Clone * const _parent){
	if (_parent == NULL){
		unsetParent();
	}else{
		parent = _parent;
		setBuried(true);
	}
}

void Clone::unsetParent(){
	parent = NULL;
	setBuried(false);
}

