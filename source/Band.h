#ifndef BAND_H_
#define BAND_H_

#include <sstream>
#include <string>

class Band
{
public:
	Band();
	Band(std::string _clone_name, int _size, int _copy_id);
	virtual ~Band();
	std::string getName() const {return name;}
	int getSize() const {return size;}
	int getCopyID() const {return copy_id;}
	std::string getCloneName() const {return clone_name;}
	
protected:
	std::string name;
	std::string clone_name;
	int size;
	int copy_id;
};

#endif /*BAND_H_*/
