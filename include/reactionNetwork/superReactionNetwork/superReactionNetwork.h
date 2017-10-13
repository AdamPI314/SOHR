#ifndef __SUPERREACTIONNETWORK_H_
#define __SUPERREACTIONNETWORK_H_

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

namespace reactionNetwork_sr {
	namespace mt = misc_template;
	namespace rsp = relationshipParser_sr;
	namespace pgt = propagator_sr;
	namespace eppstein = eppstein_algorithm;


	//when and where type
	typedef std::pair<rsp::my_time_t, std::size_t> when_where_t;
	//species index and flux
	typedef std::pair<std::size_t, double> spe_index_flux_t;

	//atom M-matrix type
	typedef std::map<std::string, matrix_sr::size_t_matrix_t > atom_M_matrix_t;
	//atom R-matrix type
	typedef std::map<std::string, matrix_sr::path_R_matrix_t> atom_R_matrix_t;

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
		rsp::index_int_t reaction_index; //which reaction it is
		//just for search algorithm
		double edge_weight;

		//stoichiometric coefficient, like A+A->B+B, s_coef_reactant= s_coef_product= 2
		rsp::stoichiometric_coef_t s_coef_reactant; // stoichiometric coefficient, default value = 1
		rsp::stoichiometric_coef_t s_coef_product; // stoichiometric coefficient, default value = 1
		EdgeProperties_graph() :edge_index(0), reaction_index(0), edge_weight(0.0), s_coef_reactant(1.0), s_coef_product(1.0) {}
		EdgeProperties_graph(rsp::index_int_t reaction_index_in, rsp::stoichiometric_coef_t s_coef_reactant_in, rsp::stoichiometric_coef_t s_coef_product_in) :
			edge_index(0), reaction_index(reaction_index_in), edge_weight(0.0), s_coef_reactant(s_coef_reactant_in), s_coef_product(s_coef_product_in) {}
	};


	class superReactionNetwork :public Graph_bundled<VertexProperties_graph, EdgeProperties_graph> {
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


	protected:
		//will write a random class

	protected:
		//for search algorithm
		//shift of edge weight
		double delta_edge_weight;
		//self-defined maximum value of edge_weight
		double maximum_edge_weight;
		//self-defined minimum value of edge_weight
		double minimum_edge_weight;
		//shortest distance from dijkstra or bellman_ford algorithm
		double shortest_distance;

	protected:
		// Graph theoretic enumeration of chemical pathways
		// IRPC paper, http://www.tandfonline.com/doi/abs/10.1080/0144235X.2016.1220774
		// M-matrix
		//std::vector<std::vector<std::size_t> > M_matrix;
		atom_M_matrix_t atom_M_matrix;
		// R-matrix
		//matrix_sr::path_R_matrix_t R_matrix;
		atom_R_matrix_t atom_R_matrix;

	public:
		superReactionNetwork(std::vector<double> uncertainties, std::size_t random_seed_for_this_core, std::string cwd_in);
		~superReactionNetwork();

	public:
		//random
		std::size_t return_index_randomly_given_probability_vector(const std::vector<double>& probability) const
		{
			return rand->return_index_randomly_given_probability_vector(probability);
		}

	public:
		////Read initial configuration file named "setting.cfg"
		bool read_init_config();
		//read species and reaction info, parse reaction network, transform all stuff into reaction network format
		//construct the edgeVector and the edgeProperty so that we can build the reaction network
		void read_chem_out_spe_reaction_network(std::vector<VertexPair> &edgeVector, std::vector<EdgeProperties_graph> &edgePro, std::vector<VertexProperties_graph>& vertex_info);
		// update species super atom info
		// super atom is defined as a general atom, which could be any atom, the number equals to the sum of all atoms
		void update_super_atom_info(std::string super_atom = "X");
	public:
		void set_species_initial_concentration();

	public:
		void set_min_time(rsp::my_time_t min_time_in) { this->min_time = min_time_in; }
		void set_max_time(rsp::my_time_t max_time_in) { this->max_time = max_time_in; }
		void set_sys_min_time(rsp::my_time_t sys_min_time_in) { this->sys_min_time = sys_min_time_in; }
		void set_tau(rsp::my_time_t end_time_in) { this->tau = end_time_in; }

		rsp::my_time_t get_max_time() const;

	public:
		//return path end time from file
		rsp::my_time_t return_tau() const;

		//return initial species
		rsp::index_int_t return_initial_spe() const;

	public:
		//86 and 89 are dead species, they transform to each other very fast
		void set_dead_spe();

	public:
		//initiate the graph according to the input Edge vector
		void initGraph(const vector<VertexPair>& edgeVector, const std::vector<EdgeProperties_graph>& edgePro);
		//update the vertex info according the vertex vector, std::vector<vertex_t>.
		//they should be in the same order.
		void update_vertex_info(const std::vector<VertexProperties_graph> & vertex_info);

		/*
		* set this->species.out_reaction_index_s_coef_v
		*/
		void set_spe_out_reaction_info();
		void search_for_out_reaction(vertex_t vertex, mt::vector_sr<rsp::reaction_index_s_coef_t> & reaction_index_s_coef_v);

		/*
		* search edges for out species of a reaction, record the species index and weight
		* the weight is determined by follow str, like atom "H"
		* version 1, weight = stoichiometric coefficient * num of atoms
		* reaction_index_s_coef_t
		*/
		bool search_for_out_spe_v1(rsp::index_int_t reaction_index, std::vector<rsp::spe_index_weight_t > &out_spe_index_weight_v, std::string atom_followed = "H");

		/*
		* search edges for out species of a reaction, record the species index and weight
		* the weight is determined by follow str, like atom "H"
		* version 2, weight = stoichiometric coefficient
		*/
		bool search_for_out_spe_v2(rsp::index_int_t reaction_index, std::vector<rsp::spe_index_weight_t > &out_spe_index_weight_v, std::string atom_followed = "H");

		/*
		* search edges for out species of a reaction, record the species index and weight
		* the weight is determined by follow str, like atom "H"
		* version 2, weight = num of atoms
		*/
		bool search_for_out_spe_v3(rsp::index_int_t reaction_index, std::vector<rsp::spe_index_weight_t > &out_spe_index_weight_v, std::string atom_followed = "H");

		/*
		* set this->reaction.out_spe_index_weight_v
		* the weight is determined by follow str, like atom "H"
		* version 1, weight = stoichiometric coefficient * num of atoms
		*/
		bool set_reaction_out_spe_info(std::string atom_followed);
		void set_reaction_out_spe_info();

		void set_out_spe_index_branching_ratio_map_map(std::string atom_followed);
		void set_out_spe_index_branching_ratio_map_map();


	public:
		/*
		* print reaction network for drawing network use
		* node.csv
		* edge.csv
		*/
		void print_network(std::string filename = "/output/node.csv");
		/*
		* just print, for test use
		*
		*/

		void print();
		/*
		* print initial json file for species labelling
		*/
		void print_initial_spe_label_json(std::string filename = "/output/initial_spe_label.json") const;
		/*
		* print initial json file for reactions labelling
		*/
		void print_initial_reaction_label_json(std::string filename = "/output/initial_reaction_label.json") const;

	public:
		/*set atom followed reaction rates at time t
		1. update all reaction rates at time t
		2. update edges weights
		*/
		/*
		* version 1, 2, 3 are problematic at this moment, 4/19/2016, only the fourth version is working, natural log of branching ratio version
		*/
		//version zero, edge weights equal reaction rates*stoichoimetric coefficient, don't take the negative of edge weight, for test use
		void set_atom_followed_edge_weight_at_time_v0(rsp::my_time_t t, std::string atom, double default_weight);
		//version one, edge weights equal reaction rates*stoichoimetric coefficient
		void set_atom_followed_edge_weight_at_time_v1(rsp::my_time_t t, std::string atom, double default_weight);
		//version two, edge weights equal pseudo-first order rate constant
		void set_atom_followed_edge_weight_at_time_v2(rsp::my_time_t t, std::string atom, double default_weight);
		//version 3, branching ratio
		void set_atom_followed_edge_weight_at_time_v3(rsp::my_time_t t, std::string atom, double default_weight);
		//version 4, natural logarithm of branching ratio, got to be negative of log(0.0)=self_defined magic number
		void set_atom_followed_edge_weight_at_time_v4(rsp::my_time_t t, std::string atom, double default_weight);
		/*
		* @ a time t
		* Dijkstra algorithm, just use the reaction rate as the edge weight, no s_coef_reactant or s_coef_product
		* two versions, first use std::set, second use std::priority_queue
		* it is like a template method, call method overloaded in derived classes
		*/
		/*
		* Modify Dijkstra's Algorithm to get the Shortest Path Between Two Nodes
		* The best way to approach this is to think in terms of paths from each reachable node back to the source,
		* commonly denoted by v. As you update each given node's distance, keep track of its direct parent on that path to v.
		* That node is the given node's parent in the shortest-distance-to-v tree. When you've built that parent map,
		* to find the shortest path from v to w build the path in reverse order: w, parent[w], parent[parent[w]], ...
		*/
		double dijkstra_w_set(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in, std::string atom_followed = "H");
		double dijkstra_w_priority_queue(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in, std::string atom_followed = "H");
		double dijkstra_boost(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in);
		double dijkstra_boost_v2(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in);
		//vectors to store the predecessors(p), vectors to store the distances from the root(d)
		void dijkstra_boost(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in,
			std::vector<Vertex> &predecessors, std::vector<double> &distance);
		void dijkstra_boost_v2(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in,
			std::vector<Vertex> &predecessors, std::vector<Edge> &predecessors_edges, std::vector<double> &distance);
		void bellman_ford_boost(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in,
			std::vector<Vertex> &predecessors, std::vector<double> &distance);
		void bellman_ford_boost_v2(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in,
			std::vector<Vertex> &predecessors, std::vector<Edge> &predecessors_edges, std::vector<double> &distance);
		//dijkstra algorithm on the reverse graph
		void dijkstra_reverse_graph_boost(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in,
			std::vector<reverse_vertex_t>& predecessors, std::vector<double>& distance);
		void dijkstra_reverse_graph_boost_v2(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in,
			std::vector<reverse_vertex_t>& predecessors, std::vector<reverse_edge_t> &predecessors_edges, std::vector<double>& distance);
		//bellman_ford algorithm on the reverse graph
		void bellman_ford_reverse_graph_boost(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in,
			std::vector<reverse_vertex_t>& predecessors, std::vector<double>& distance);
		void bellman_ford_reverse_graph_boost_v2(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in,
			std::vector<reverse_vertex_t>& predecessors, std::vector<reverse_edge_t> &predecessors_edges, std::vector<double>& distance);
		//for print out the path
		void print_dijkstra_algorithm_shortest_path(std::map<vertex_t, vertex_t > vertex_parent_map, const vertex_t s, const vertex_t v) const;

		//k shortest path algorithm, sidetrack_tree algorithm
		/**
		* Computes the K shortest paths (allowing cycles) in a graph from node s to node t in graph G using Eppstein's
		* algorithm. ("Finding the k Shortest Paths", Eppstein)
		* <p/>
		* Some explanatory notes about how Eppstein's algorithm works:
		* - Start with the shortest path in the graph from s to t, which can be calculated using Dijkstra's algorithm.
		* - The second shortest path must diverge away from this shortest path somewhere along the path at node u, with
		* edge (u,v) known as a "sidetrack" edge. It then follows the shortest path from v to t.
		* - In general, let T represent the shortest path tree rooted at target node t, containing for each node v, an edge
		* (v,w) to the parent of node v (w) along its shortest path to node t.
		* - All other edges (u,v) in the graph are "sidetrack" edges. In other words, this includes any edge (u,v) for
		* which node v is not on the shortest path from u to t.
		* - All paths from s to t can be uniquely represented by the sequence of sidetrack edges that appear in the path.
		* - All non-sidetrack edges in the path are represented implicitly in the shortest path tree, T.
		* - All paths from s to t can be represented in a tree, H, where the root node has no sidetrack edges, its children
		* are all of the possible one sidetrack paths, the children of each of these sidetracks are the paths with a
		* second sidetrack edge, and so on.
		* - Each path from s to t corresponds to a path from the root of H to some descendant node in H.
		* - The tree is a heap, as each additional sidetrack creates a new path whose cost is either the same as or greater
		* than the cost of its parent path.
		* - One possible way to find the k shortest paths, then, is to maintain a priority queue of candidate paths where
		* once a path is pulled from the priority queue, its children in the heap are added to the priority queue.
		* - Two observations about this heap:
		* - Each node in H has at most O(E) children, where E is the number of edges in the graph.
		* - If G has looped paths, then this heap is infinite.
		* - Eppstein's algorithm is basically a scheme for re-organizing this heap so that each heap node has at most 4
		* children, which is O(1) instead of O(E) - performing potentially much less work for each path found.
		* - This re-organized heap R has the following form, in terms of the original heap H:
		* - For each heap node in H, select its best child, B.
		* - This child is known as a "cross edge" child.
		* - Child B has N siblings in H.
		* - These N siblings are removed from its parent in heap H and are made instead to be descendants of B in heap R.
		* - The best child, B, has up to three of its siblings placed as direct children (along with its own best "cross
		* edge" child in the original heap H, this yields a total of four possible children in R).
		* - Among its sibling-children, one of its children in R is the root of a new binary heap containing all of the
		* siblings of B in H that are sidetrack edges directed away from the same node in G.
		* - All of the remaining siblings of B in H are placed in a binary heap with B as the root in R.
		* - Because this is a binary heap, B has up to two children of this type.
		* - Thus, B has at most four children in R. Any node in R that is a "cross-edge" child has up to four children,
		* whereas other nodes fall inside binary heaps and are limited to having at most two children.
		*/
		typedef std::deque<vertex_t> path_vertex_t;
		typedef std::deque<Edge> path_edge_t;
		void eppstein_alogrithm(const Vertex s, const Vertex t, double time_in = 0.0, std::string atom = "H");
		//create the sidetrack edge tree, insert the first level nodes
		void initiate_sidetack_tree(eppstein::sidetrack_tree &st_tree, const path_vertex_t &path_v, path_edge_t &path_e, const std::vector<double> &distance_r) const;
		//add node to sidetrack tree recursively
		void add_vertex_to_sidetrack_tree_recursively(eppstein::sidetrack_tree &st_tree, const eppstein::sidetrack_tree::Vertex curr_vertex_index_in_tree, eppstein::vertex_index_t curr_vertex_index_in_original_graph,
			const std::vector<double>& distance, const std::vector<double> &distance_r, std::size_t level_from_current_vertex_in_original_graph = 0, const std::size_t level_from_root = 0) const;
		//add element to priority queue
		void k_shortest_path_with_pq(const eppstein::sidetrack_tree &st_tree, eppstein::len_path_pq_t & lp_pq, std::vector<eppstein::path_info_t > &K_shortest_path_info_vec) const;
		//retrieve path from path_info_vector
		void retrieve_path_info(const std::vector<Vertex> &predecessors_r, const std::vector<reverse_edge_t> &predecessor_edges_r, const std::vector<eppstein::path_info_t > &K_shortest_path_info_vec,
			const Vertex s, const Vertex t, const std::size_t top_k = 5);

	public:
		/*
		* random pick next reaction based on their reaction rates and stoichiometric coefficient
		*/
		rsp::index_int_t random_pick_next_reaction(vertex_t curr_vertex);
		/*
		* known the reaction, randomly pick next species based on their weight
		*/
		vertex_t random_pick_next_spe(rsp::index_int_t reaction_index, std::string atom_followed = "H");

	public:
		// XS1R1S2--> X and S1R1S2
		void split_atom_followed_and_pathway(std::string str_in, std::string &atom_followed, std::string &pathway) const;
		/*
		* parse pathway like S2R2S5R4S1 to two vectors
		* species vector 2->5->1
		* reaction vector 2-4
		*/
		bool parse_pathway_to_vector(std::string pathway_in, std::vector<rsp::index_int_t> &spe_vec, std::vector<rsp::index_int_t> &reaction_vec) const;
		int get_number_of_elements() const;

	public:
		//update out reaction rate of current vertex at a specific time point
		virtual bool update_reaction_rate(double in_time, vertex_t curr_vertex) = 0;
		//update rate at specific time point
		virtual bool update_reaction_rate(double in_time) = 0;
		void set_reaction_rate(vertex_t i, double reaction_rate);
		void set_is_reaction_rate_nonzero_from_setting_file();
		// set is reaction rate nonzero based the real reaction rate
		void set_is_reaction_rate_nonzero_from_previous_iteration();

	public:
		/*
		* reaction time from importance sampling, exact time
		* if reaction_time> tau, let it be, don't cut it off
		*/
		virtual double reaction_time_from_importance_sampling_without_cutoff(rsp::my_time_t curr_time, vertex_t curr_spe, double Y) = 0;
		/*
		* reaction time from importance sampling, exact time
		* if reaction_time> tau, cut it off
		*/
		virtual double reaction_time_from_importance_sampling(rsp::my_time_t curr_time, vertex_t curr_spe, double Y) = 0;

		/*
		* chattering group reaction time from importance sampling, exact time
		* if reaction_time> tau, let it be, don't cut it off
		*/
		virtual double chattering_group_spe_reaction_time_from_importance_sampling_without_cutoff(rsp::my_time_t curr_time, vertex_t index, double Y) = 0;


		//for concreteReactionNetwork and reactionNetworkSolver
	public:
		//set Prob_max(tau^{j}|t+tau^{j-1};S^{j-1}), with pathway_end_time fixed
		virtual bool set_spe_prob_max_at_a_time(double init_time, double end_time, size_t index_t) = 0;
		//set the  prob_min of a species at a given time
		virtual bool set_spe_prob_min_at_a_time(double init_time, double end_time, size_t index_t) = 0;
		virtual double get_spe_prob_max_at_a_time(double init_time, double end_time, size_t index_t) const = 0;
		virtual double get_spe_prob_min_at_a_time(double init_time, double end_time, size_t index_t) const = 0;
	public:
		//defer to subclass
		virtual double evaluate_spe_concentration_at_time(double time, std::size_t index = 0) const = 0;
		virtual double evaluate_chattering_group_spe_ss_prob_at_time(double in_time, size_t index = 0) const = 0;

	public:
		//prob that a spe will react in time range
		double prob_spe_will_react_in_a_time_range(double init_time, double pathway_end_time, size_t curr_spe);
		virtual double prob_chattering_group_spe_will_react_in_a_time_range(double init_time, double pathway_end_time, size_t index) = 0;

	public:
		/*
		* input- reaction_time, next reaction time
		* input- next reaction R
		* input- current species spe
		* input- next species spe
		* output- reaction and spe branching ratio
		* deal with fast reactions, next_spe is not a product of next_reaction, but inter-convert to one product of next_reaction rapidly
		*/
		/* use out_spe_index_branching_ratio_map_map, do not need to search for out species and calculate spe branching ratio each time*/
		double reaction_spe_branching_ratio(double reaction_time, rsp::index_int_t curr_spe, rsp::index_int_t next_reaction, rsp::index_int_t next_spe, std::string atom_followed = "H");

	public:
		/*
		* return when it is, where we are, for MPI
		* move one step
		* if we reach trapped species, randomly select one species from the trapped species pair to proceed
		*/
		when_where_t pathway_move_one_step(double time, vertex_t curr_vertex, std::string &curr_pathway_local, std::string atom_followed = "H");

		/*
		* simulate in a specific time range once, knowing the initial species,  manually set which atom to follow, and return the pathway
		*/
		std::string pathway_sim_once(double init_time, double end_time, vertex_t init_vertex, std::string atom_followed = "H");

	public:
		/*
		* return where we are, when it is, for MPI ,for pathway probability
		* force the next reaction to be  next_reaction, next species to be next_spe, because what we want is just pathway probability
		* directly set the next reaction and next species the ones we want
		*/
		double pathway_prob_sim_move_one_step(double when_time, vertex_t curr_spe, rsp::index_int_t next_reaction, vertex_t next_spe, double &pathway_prob, std::string atom_followed = "H");
		//input a pathway, return its pathway prob
		/*
		//parse pathway to spe and reaction vector
		std::vector<size_t> spe_vec; std::vector<size_t> reaction_vec;
		this->parse_pathway_to_vector(pathway_in, spe_vec, reaction_vec);
		*/
		double pathway_prob_input_pathway_sim_once(const double init_time, const double pathway_end_time, const std::vector<rsp::index_int_t> &spe_vec, const std::vector<rsp::index_int_t> &reaction_vec, std::string atom_followed = "H");


	public:
		// initiate M-Matrix
		void initiate_M_matrix(std::string atom_followed);
		// initiate M-Matrix following all atoms
		void initiate_M_matrix();
		// print out M-matrix
		void print_M_matrix(std::string atom_followed);

		matrix_sr::size_t_matrix_t return_M_matrix(std::string atom_followed);

		// initiate R-Matrix, v1, store reaction index
		void initiate_R_matrix_v1(std::string atom_followed);
		// initiate R-Matrix, v2, store edge index, notice edge index is different than reaction index
		// one reaction might corresponds to multiple edges
		void initiate_R_matrix_v2(std::string atom_followed);
		// initiate R-matrix following all atoms
		void initiate_R_matrix();

		matrix_sr::path_R_matrix_t return_R_matrix(std::string atom_followed);

		// print out R-matrix
		void print_R_matrix(std::string atom_followed);
		// get M_matrix element
		std::size_t get_M_matrix_element(std::string atom_followed, std::size_t i, std::size_t j);
		// get R_matrix element
		matrix_sr::path_R_matrix_element_t get_R_matrix_element(std::string atom_followed, std::size_t i, std::size_t j);

	public:
		// convert R matrix pathway representation to pathway string, like S0R1S1
		std::string R_matrix_path_representation_to_string(matrix_sr::path_t p);
		// tell whether a path contains zero reaction rate reactions or not
		bool contains_zero_reaction_rate_reactions(matrix_sr::path_t p);
		std::vector<std::string> get_path_string_element_i_j(const matrix_sr::path_R_matrix_t &pRm, std::size_t i, std::size_t j);

		// return a vector of pathway string of length n, species i-->j, change matrix element, so that each matrix element can not contain more than topN sub-path
		std::vector<std::string> get_path_string_update_matrix_element_i_j_topN(matrix_sr::path_R_matrix_t &pRm, const std::size_t i, const std::size_t j,
			const std::string atom_followed = "H", const std::size_t topN = 10, const double start_time = 0.0, const double end_time = 1.0);
		// save vector of path string to file
		void path_string_vector_s2f(std::vector<std::string> vs, std::string filename = "./output/pathname.csv");

		// save heuristic path string to file
		void heuristic_path_string_vector_s2f(std::string atom_followed = "H", std::size_t n = 1, std::string filename = "./output/heuristic_pathname.csv");
		// save heuristic path string to memory, n is pathway length
		std::set<std::string> heuristic_path_string_vector_s2m(std::string atom_followed = "H", std::size_t n = 1);
		// sorted by path length
		std::set<std::string> heuristic_path_string_vector_sorted_based_on_path_length(std::string atom_followed = "H", std::size_t n = 1, std::size_t topN = 10);
		// sorted by pathway probability, here topN is a const for all species because of matrix multiplication, matrix power N is the same for all species
		std::set<std::string> heuristic_path_string_vector_sorted_based_on_path_prob(std::string atom_followed = "H", std::size_t n = 1, std::size_t topN = 10, double end_time_ratio = 1.0);


		// path weight for path sorting
		// method 1, by path length
		double calculate_path_weight_path_length(std::string path);
		// method 2, by pathway probability at a time
		double calculate_path_weight_based_on_path_probability(std::string path, std::string atom_followed = "H", double start_time = 0.0, double end_time = 1.0);

		// for a stage number, generate pathway and save to file, this stage number could be time stage number or iteration stage number
		// sort by path length, following all elements
		std::size_t heuristic_path_string_vector_by_stage_number_path_length_all_elements(const std::size_t stage_n, std::string filename = "./output/heuristic_pathname_1.csv", std::size_t = 10);
		// sort by pathway probability at a time, using one trajectory, following all elements
		std::size_t heuristic_path_string_vector_by_stage_number_path_prob_all_elements(const std::size_t stage_n, std::string filename = "./output/heuristic_pathname_1.csv", std::size_t topN = 10, double end_time_ratio = 1.0);
		std::size_t heuristic_path_string_vector_by_stage_number_path_prob_all_elements_s2m(const std::size_t stage_n, std::vector<std::string> &path_all_v, std::size_t topN = 10, double end_time_ratio = 1.0);
		// sort by pathway probability at a time, using one trajectory, following "X", super atom
		std::size_t heuristic_path_string_vector_by_stage_number_path_prob_super_element(const std::size_t stage_n, std::string filename = "./output/heuristic_pathname_1.csv", std::size_t = 10, double end_time_ratio = 1.0);

		// species index with initial concentrations
		std::set<std::size_t> return_species_index_with_initial_concentration() const;
		std::set<std::size_t> return_species_index_without_initial_concentration() const;
		std::set<std::pair<std::size_t, double> > return_species_index_and_initial_concentration() const;

		// generate pathway by running monte carlo trajectories, by following one element
		void generate_path_by_running_monte_carlo_trajectory_s2m(std::vector<statistics > &statistics_v, std::size_t Ntrajectory, std::string atom_followed = "H", double end_time_ratio = 1.0);
		// generate pathway by running monte carlo trajectories, by following all elements
		std::size_t generate_path_by_running_monte_carlo_trajectory_all_elements_s2m(std::vector<statistics > &statistics_v, std::size_t Ntrajectory, double end_time_ratio = 1.0);

		std::vector<rsp::element_info> return_element_vecotr() const;

	};//class reactionNetwork




}/*namespace reactionNetwork_sr*/

#endif