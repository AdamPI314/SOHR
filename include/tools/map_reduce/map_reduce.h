#ifndef __MAP_REDUCE_H_
#define __MAP_REDUCE_H_

#include "../misc/global_macros.h"


#if defined(__USE_MPI_)

#include <boost/mpi.hpp>
#include <boost/serialization/map.hpp>
#include <string>

struct merge_maps {
	std::map<std::string, int> operator() (std::map<std::string, int> &l, const std::map<std::string, int> &r);
};

namespace boost {
	namespace mpi {
		template <>
		struct is_commutative<merge_maps, std::map<std::string, int> > : mpl::true_ {};
	}
}

#endif // __USE_MPI_


#endif
