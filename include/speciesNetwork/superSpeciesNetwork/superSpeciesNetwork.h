#ifndef __SUPERSPECIESNETWORK_H_
#define __SUPERSPECIESNETWORK_H_

#include "../../tools/misc/misc_template.h"
#include "../../random/random.h"

//In the C and C++ programming languages, #pragma once is a non-standard but widely supported
//preprocessor directive designed to cause the current source file to be included only once in a
//single compilation.
//#pragma once

#include <fstream>
#include <utility>                          // for std::pair
#include <algorithm>                        // for std::for_each and std::find_if
#include <iterator>
#include <ctime>
#include <limits> //std::numeric_limits
#include <string> //for std::string
#include <vector> //for std::vector
#include <map> //for std::map
#include <queue> //for std::priority_queue
#include <unordered_set>

#include "../../tools/misc/misc_template.h"
//#include "print_graph.h"
#include "../../tools/misc/graph_bundled.h"
#include "../../search_algorithm/eppstein_algorithm/eppstein_algorithm.h"
#include "../../tools/matrix/matrix_sr.h"

#include "../../pathwayHandler/pathwayHandler.h"

#include <boost/config.hpp>
#include <boost/utility.hpp> // for boost::tie
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp> //for boost::regex
#include <boost/unordered_set.hpp> //for boost::unorder_set
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp> //for boost::random::uniform_real_distribution
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp> //for boost::filesystem
#include <boost/tokenizer.hpp> //for boost::tokenizer
#include <boost/random/random_device.hpp> //for boost::random::random_device
#include <boost/property_tree/ptree.hpp> //for property_tree

#include "../../fileIO/commandLineConfigFileReader/clcf_parser.h"
#include "../../statistics/statistics.h"
//typedef, and relationParser and graph related stuff
#include "../../relationshipParser/relationshipParser.h"
#include "../../propagator/dlsodePropagator/dlsodePropagator.h"
#include "../../propagator/ssaPropagator/ssaPropagator.h"
#include "../../propagator/SOHRPropagator/SOHRPropagator.h"

//Try including nr3.h after all boost and system headers!!!
//#include "../cubicSpline/interp_1d.h"

namespace speciesNetwork_sr {
	namespace mt = misc_template;
	namespace rsp = relationshipParser_sr;
	namespace pgt = propagator_sr;
	namespace eppstein = eppstein_algorithm;


	//when and where type
	typedef std::pair<rsp::my_time_t, std::size_t> when_where_t;
	//species index and flux
	typedef std::pair<std::size_t, double> spe_index_flux_t;

	//more reaction network related thing, like edge pair, edge property, and vertex property etc.
	//vertex type, every vertex in the graph stores the index of species
	typedef rsp::index_int_t vertex_t;
	struct VertexProperties_graph {
		//Look "chem.out" for detailed info, all the species are labeled with integer,
		//here make them start from 0
		vertex_t vertex;
	};

	//Edge correspond to reactions in file "chem.out"
	struct EdgeProperties_graph {
		rsp::index_int_t edge_index; //its edge index in the graph

		EdgeProperties_graph() :edge_index(0) {}
		EdgeProperties_graph(rsp::index_int_t reaction_index_in) :edge_index(0) {}
	};


	class superSpeciesNetwork :public Graph_bundled<VertexProperties_graph, EdgeProperties_graph> {
	protected:
		//current working directory
		std::string cwd;
		////random seed for this core
	protected:
		//configurations, read configuration file named "setting.json"
		boost::property_tree::ptree rnk_pt;

	protected:
		boost::uint32_t random_seed_for_this_core;

		//random number generator
		random_sr::random* rand;

	protected:
		std::vector<rsp::element_info> element_v;
	protected:
		//species and species information
		std::vector<rsp::spe_info_base> species_network_v;
		//species name and species index, stored in a map
		rsp::spe_name_index_map_t spe_name_index_map;
		//dead/end species, species that can not be destroyed, it might have a linkage in the reaction network, but the destruction rate
		//is always zero! set this vector manually. check the last line of file int_drc.dat, if the cumulative destruction rate is zero,
		//the species is a dead species
		std::set<vertex_t> dead_species;
	protected:
		//shared pointer of chattering
		std::shared_ptr<chattering_sr::chattering> sp_chattering_rnk;
	protected:
		/*
		 *	we have two spaces here, one is reaction mechanism space, in which
		 *	we can have duplicated reactions
		 *	we can have same species on both sides of reactions
		 *	we can have two direction reactions,
		 *	the other space, is reaction network space, for convenience,
		 *	we don't have duplicated reaction
		 *	we can not have same species on both sides of reaction
		 *	just have one direction reaction
		*/
		//just need the basic information
		std::vector<rsp::reaction_info_base> reaction_network_v;
		//save a copy of reaction info, full information
		std::vector<edge_iter> edge_index_to_edge_iterator;
	protected:
		//number of edges
		std::size_t num_edges;
		//number of vertices
		std::size_t num_vertices;
	protected:
		//time
		rsp::my_time_t min_time;
		rsp::my_time_t max_time;
		//minimum time below which the program/system will crash, like exponential, or some other operation will be screw up.
		rsp::my_time_t sys_min_time;

	protected:
		//reference time of the whole system, to a combustion system, tau is the ignition delay time
		rsp::my_time_t tau;

	public:
		superSpeciesNetwork();
		~superSpeciesNetwork();

	};//class superSpeciesNetwork




}/*namespace speciesNetwork_sr*/

#endif
