#include "Band.h"

Band::Band(std::string _clone_name, int _size, int _copy_id)
{
	std::stringstream ss;
	ss << _clone_name;
	ss << "-" << _size << "-" << _copy_id;
	name = ss.str();
	size = _size;
	copy_id = _copy_id;
	clone_name = _clone_name;
}

Band::Band(){
}

Band::~Band()
{
}
