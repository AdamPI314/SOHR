#ifndef __MY_DEBUG_H_
#define __MY_DEBUG_H_

#include <map>
#include <string>

class Globaldebug{
	public:
		std::map<std::string, bool> tags_map={{"stage_number", false},{"species_index", false}};

	public:
		void set_value(std::string entry, bool b){
			tags_map.at(entry) = b;
		}

		bool get_value(std::string entry){
			return tags_map.at(entry);
		}

	public:
		Globaldebug(){}
		~Globaldebug(){}

};

//The global variable should be declared extern in a header file included by both source files, and then defined in only one of those source files:
extern Globaldebug gdb;

#endif
