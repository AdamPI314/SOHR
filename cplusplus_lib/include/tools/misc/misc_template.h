/*
 * miscellaneous templates, typedef, etc
 */

#ifndef __MISC_TEMPLATES_H_
#define __MISC_TEMPLATES_H_

#include <vector>
#include <algorithm>

namespace misc_template{
//To store unique elements, derive std::vector<>, insert element without changing their order
template <class T>
class vector_sr:public std::vector<T>{
public:
	bool insert_sr(T t_){
		//not found
		if(std::find(this->begin(), this->end(), t_)==this->end()){
			this->push_back(t_);
			return true;
		}
		//found
		else
			return false;
	}
};//template

/*
 * make vector
 * And use it like this:
 * std::vector<int> v = make_vector<int>() << 1 << 2 << 3;
 * http://stackoverflow.com/questions/8906545/how-to-initialize-a-vector-in-c
 */
template <typename T>
class make_vector {
public:
  typedef make_vector<T> my_type;
  my_type& operator<< (const T& val) {
    data_.push_back(val);
    return *this;
  }
  operator std::vector<T>() const {
    return data_;
  }
private:
  std::vector<T> data_;
};

}//namespace



#endif
