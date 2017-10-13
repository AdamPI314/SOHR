#ifndef __SUPERSPECIESNETWORK_CPP_
#define __SUPERSPECIESNETWORK_CPP_

#include <queue>
#include <boost/property_tree/json_parser.hpp> //for json_reader
#include <boost/property_tree/xml_parser.hpp> //for write_xml
#include <boost/math/constants/constants.hpp>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/limits.hpp>
#include <boost/graph/grid_graph.hpp>
#include <numeric>

#include <math.h> /* log */

#include "../../../include/tools/misc/misc_template.h"
#include "../../../include/speciesNetwork/superspeciesNetwork/superspeciesNetwork.h"

//infinitesimal dt
#define INFINITESIMAL_DT 1.0E-14

namespace speciesNetwork_sr {



	superSpeciesNetwork::superSpeciesNetwork()
	{
	}

	superSpeciesNetwork::~superSpeciesNetwork()
	{
	}

}/*namespace speciesNetwork_sr*/

#endif
