#ifndef __SUPERREACTIONNETWORK_CPP_
#define __SUPERREACTIONNETWORK_CPP_

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
#include "../../../include/reactionNetwork/superReactionNetwork/superReactionNetwork.h"

//infinitesimal dt
#define INFINITESIMAL_DT 1.0E-14

namespace reactionNetwork_sr {

	superReactionNetwork::superReactionNetwork(std::vector<double> uncertainties, std::size_t random_seed_for_this_core, std::string cwd_in)
	{

		//current working directory
		this->cwd = cwd_in;

		//read configuration file "setting.json"
		boost::property_tree::read_json(this->cwd + std::string("/input/setting.json"), rnk_pt, std::locale());
		read_init_config();

		//random seed for this core
		this->random_seed_for_this_core = static_cast<boost::uint32_t>(random_seed_for_this_core);

		rand = new random_sr::random(this->random_seed_for_this_core);

		std::vector<VertexPair> edgeVector; std::vector<EdgeProperties_graph> edgePro; std::vector<VertexProperties_graph> vertex_info;
		this->read_chem_out_spe_reaction_network(edgeVector, edgePro, vertex_info);
		//update super atom info
		this->update_super_atom_info(rnk_pt.get<std::string>("pathway.super_atom"));

		set_species_initial_concentration();

		//initialize graph
		initGraph(edgeVector, edgePro);

		update_vertex_info(vertex_info);
		set_spe_out_reaction_info();
		set_reaction_out_spe_info();

		set_out_spe_index_branching_ratio_map_map();

		//set dead species
		set_dead_spe();

		initiate_M_matrix();
		initiate_R_matrix();

		set_is_reaction_rate_nonzero_from_setting_file();

		//print();
		//print_network();
		//std::cout << "test" << std::endl;
		//print_initial_spe_label_json();

	}

	superReactionNetwork::~superReactionNetwork()
	{
		delete rand;
	}

	////Read initial configuration file named "setting.cfg"
	bool superReactionNetwork::read_init_config()
	{
		set_min_time(this->rnk_pt.get<cf_parser::my_time_t>("time.min_time"));
		set_max_time(this->rnk_pt.get<cf_parser::my_time_t>("time.max_time"));
		set_sys_min_time(this->rnk_pt.get<cf_parser::my_time_t>("time.sys_min_time"));
		set_tau(this->rnk_pt.get<cf_parser::my_time_t>("time.tau"));

		return true;
	}

	void superReactionNetwork::read_chem_out_spe_reaction_network(std::vector<VertexPair> &edgeVector, std::vector<EdgeProperties_graph> &edgePro, std::vector<VertexProperties_graph>& vertex_info)
	{
		/*
		 * read species information
		 */
		std::vector<rsp::spe_info> species_v;
		//read element, species and reaction information
		rsp::relationshipParser::read_chem_out_ele_spe(element_v, species_v, spe_name_index_map, this->cwd + "/input/chem.out");

		for (std::size_t i = 0; i < species_v.size(); ++i) {/*for1*/
			rsp::spe_info_base species_network_temp;
			species_network_temp.prob_max = species_v[i].prob_max;
			species_network_temp.prob_min = species_v[i].prob_min;
			species_network_temp.reaction_k_index_s_coef_v = species_v[i].reaction_k_index_s_coef_v;
			species_network_temp.spe_component = species_v[i].spe_component;
			species_network_temp.spe_conc = species_v[i].spe_conc;
			species_network_temp.spe_index = species_v[i].spe_index;
			species_network_temp.spe_name = species_v[i].spe_name;
			species_network_temp.survival_probability = species_v[i].survival_probability;

			species_network_v.push_back(species_network_temp);
		}/*for1*/

		//update vetex information, here our vertex just store the index of vextex in species space, look the difinition of VertexProperties_graph
		vertex_info.resize(species_network_v.size());
		for (std::size_t i = 0; i < species_network_v.size(); ++i) {
			vertex_info[i].vertex = species_network_v[i].spe_index;
		}


		/*
		 * read reaction information
		 */
		std::vector<rsp::reaction_info> reaction_v;
		//include duplicated reactions
		rsp::relationshipParser::read_chem_out_reaction(species_v, reaction_v, spe_name_index_map, this->cwd + "/input/chem.out");
		rsp::relationshipParser::set_reaction_net_reactant_product(reaction_v);
		//reactionNetwork and chemkin index lookup table, notice chemkin index is Fortran index style
		rsp::reactionNetwork_chemkin_index_map_t reactionNetwork_chemkin_index_map;
		rsp::relationshipParser::read_reactionNetwork_chemkin_index_map(reactionNetwork_chemkin_index_map, this->cwd + "/input/chem.out");


		std::size_t edge_counter = 0;
		//construct A->B transition
		EdgeProperties_graph edgePro_tmp;

		for (rsp::reactionNetwork_chemkin_index_map_t::const_iterator itr = reactionNetwork_chemkin_index_map.begin(); itr != reactionNetwork_chemkin_index_map.end(); ++itr) {
			rsp::reaction_info_base reaction_info_base_temp;
			//just need the first one if have multiple duplicated reactions
			//std::cout << itr->first << "\t-->\t";
			//std::cout << itr->second[0] << "\t";
			//need to convert Fortran style index into C++ style index
			rsp::index_int_t reaction_v_ind = static_cast<rsp::index_int_t>(abs(itr->second[0])) - 1;
			//std::cout << reaction_v[reaction_v_ind].reaction_name << "\t";
			//std::cout << reaction_v[reaction_v_ind].reaction_direction << "\n";

			//forward reaction
			if (itr->second[0] > 0) {
				reaction_info_base_temp.reaction_direction = rsp::forward;

				//edge vector
				for (std::size_t i = 0; i < reaction_v[reaction_v_ind].net_reactant.size(); ++i) {/*for i*/
					for (std::size_t j = 0; j < reaction_v[reaction_v_ind].net_product.size(); ++j) {/*for j*/
						//std::cout << "[" << reaction_v[reaction_v_ind].net_reactant[i].first << "," << reaction_v[reaction_v_ind].net_product[j].first << "]" << "\t";
						edgeVector.push_back(std::make_pair(reaction_v[reaction_v_ind].net_reactant[i].first, reaction_v[reaction_v_ind].net_product[j].first));
						edgePro_tmp.edge_index = edge_counter; ++edge_counter;
						edgePro_tmp.reaction_index = itr->first;
						edgePro_tmp.s_coef_reactant = reaction_v[reaction_v_ind].net_reactant[i].second;
						edgePro_tmp.s_coef_product = reaction_v[reaction_v_ind].net_product[j].second;

						edgePro.push_back(edgePro_tmp);

					}/*for j*/
				}/*for i*/
				//std::cout << std::endl;
			}/*if*/
			//backward reaction
			else if (itr->second[0] < 0) {
				reaction_info_base_temp.reaction_direction = rsp::backward;
				//edge vector
				for (std::size_t i = 0; i < reaction_v[reaction_v_ind].net_product.size(); ++i) {/*for i*/
					for (std::size_t j = 0; j < reaction_v[reaction_v_ind].net_reactant.size(); ++j) {/*for j*/
						//std::cout << "[" << reaction_v[reaction_v_ind].net_product[i].first << "," << reaction_v[reaction_v_ind].net_reactant[j].first << "]" << "\t";
						edgeVector.push_back(std::make_pair(reaction_v[reaction_v_ind].net_product[i].first, reaction_v[reaction_v_ind].net_reactant[j].first));
						edgePro_tmp.edge_index = edge_counter; ++edge_counter;
						edgePro_tmp.reaction_index = itr->first;
						edgePro_tmp.s_coef_reactant = reaction_v[reaction_v_ind].net_product[i].second;
						edgePro_tmp.s_coef_product = reaction_v[reaction_v_ind].net_reactant[j].second;

						edgePro.push_back(edgePro_tmp);

					}/*for j*/

				}/*for i*/


			}/*else if*/
			reaction_info_base_temp.reaction_name = reaction_v[reaction_v_ind].reaction_name;
			this->reaction_network_v.push_back(reaction_info_base_temp);

		}/*for*/
		rsp::relationshipParser::spe_information_s2f(species_v);
		rsp::relationshipParser::reaction_information_s2f(species_v, reaction_v, reactionNetwork_chemkin_index_map);

	}

	void superReactionNetwork::update_super_atom_info(std::string super_atom)
	{
		for (std::size_t i = 0; i < this->species_network_v.size(); ++i) {
			this->species_network_v[i].spe_component[super_atom] = 0;

			for (auto x : this->element_v)
				this->species_network_v[i].spe_component[super_atom] += this->species_network_v[i].spe_component[x.ele_name];
		}

	}

	void superReactionNetwork::set_species_initial_concentration()
	{
		//read with json_parser as property_tree
		for (auto key1 : this->rnk_pt.get_child("chem_init.species_index_concentration")) {
			this->species_network_v[boost::lexical_cast<std::size_t>(key1.first)].spe_conc = key1.second.get_value<double>()*this->rnk_pt.get<double>("SOHR_init.massConservationFactor");
		}

		if (this->rnk_pt.get<std::string>("propagator.normalize_initial_concentration") == "yes") {
			//renormalization
			double total_conc = 0.0;
			for (std::size_t i = 0; i < this->species_network_v.size(); ++i) {
				total_conc += this->species_network_v[i].spe_conc;
			}//for

			for (std::size_t i = 0; i < this->species_network_v.size(); ++i) {
				this->species_network_v[i].spe_conc /= total_conc;
			}//for
		}//if

	}//set_initial_concentration

	rsp::my_time_t superReactionNetwork::get_max_time() const
	{
		return this->max_time;
	}

	rsp::my_time_t superReactionNetwork::return_tau() const
	{
		return this->rnk_pt.get<rsp::my_time_t>("time.tau");
	}

	rsp::index_int_t superReactionNetwork::return_initial_spe() const
	{
		return this->rnk_pt.get<rsp::index_int_t>("pathway.init_spe");
	}

	void superReactionNetwork::set_dead_spe()
	{
		//86 and 89 are dead species, they transform to each other very fast
		std::set<rsp::index_int_t> dead_spe_index;

		for (auto key1 : this->rnk_pt.get_child("pathway.dead_species")) {
			//std::cout<<key1.second.get_value<std::size_t>()<<std::endl;
			dead_spe_index.insert(key1.second.get_value<rsp::index_int_t>());
		}

		this->dead_species = dead_spe_index;
	}


	void superReactionNetwork::initGraph(const vector<VertexPair>& edgeVector, const std::vector<EdgeProperties_graph>& edgePro)
	{
		for (std::size_t i = 0; i < edgeVector.size(); ++i)
		{
			AddEdge(edgeVector[i].first, edgeVector[i].second, edgePro[i]);
		}//for

		//update edge index
		std::size_t count = 0;
		for (edge_range_t er = getEdges(); er.first != er.second; ++er.first)
		{
			edge_index_to_edge_iterator.push_back(er.first);
			properties(*er.first).edge_index = count++;
		}//for
		num_edges = get_num_edges();
		num_vertices = get_num_vertices();
		//edge_index_map=get(edge_index, graph);
	}

	void superReactionNetwork::update_vertex_info(const std::vector<VertexProperties_graph>& vertex_info) {
		vertex_range_t vp;
		int count = 0;
		for (vp = vertices(graph); vp.first != vp.second; ++vp.first, ++count)
		{
			properties(*vp.first).vertex = vertex_info[count].vertex;
			//std::cout << properties(*vp.first).vertex << std::endl;
		}
	}

	void superReactionNetwork::set_spe_out_reaction_info()
	{
		mt::vector_sr<rsp::reaction_index_s_coef_t > reaction_index_s_coef_v;
		for (std::size_t i = 0; i < this->species_network_v.size(); ++i) {
			//dead species
			if (std::find(this->dead_species.begin(), this->dead_species.end(), i) != this->dead_species.end())
				continue;
			reaction_index_s_coef_v.clear();
			search_for_out_reaction(i, reaction_index_s_coef_v);
			this->species_network_v[i].reaction_k_index_s_coef_v = reaction_index_s_coef_v;
		}

	}

	/*
	 * search out edge reactions of a species, record the reactant stoichoimetric coefficient
	 */
	void superReactionNetwork::search_for_out_reaction(vertex_t vertex, mt::vector_sr<rsp::reaction_index_s_coef_t > & reaction_index_s_coef_v) {
		for (out_edge_range_t itr = getOutEdges(vertex); itr.first != itr.second; ++itr.first) {
			reaction_index_s_coef_v.insert_sr(std::make_pair(properties(*itr.first).reaction_index, properties(*itr.first).s_coef_reactant));
		}

	}

	bool superReactionNetwork::search_for_out_spe_v1(rsp::index_int_t reaction_index, std::vector<rsp::spe_index_weight_t>& out_spe_index_weight_v, std::string atom_followed)
	{
		mt::vector_sr<rsp::spe_index_weight_t > out_spe_index_weight_v_tmp;

		boost::property_map<GraphContainer, vertex_index_t>::type vertex_id = get(vertex_index, graph);
		std::pair<vertex_iter, vertex_iter> vp;
		for (edge_range_t er = getEdges(); er.first != er.second; ++er.first) {//for
			//if found
			if (properties(*er.first).reaction_index == reaction_index) {//if
				std::size_t spe_index = get(vertex_id, target(*er.first, graph));
				//ignore duplicate elements
				out_spe_index_weight_v_tmp.insert_sr(
					std::make_pair(spe_index,
						/*spe index*/properties(*er.first).s_coef_product* this->species_network_v[spe_index].spe_component[atom_followed] /*weight*/)
				);

			}//if
		}//for

		out_spe_index_weight_v = out_spe_index_weight_v_tmp;

		return true;
	}

	bool superReactionNetwork::search_for_out_spe_v2(rsp::index_int_t reaction_index, std::vector<rsp::spe_index_weight_t>& out_spe_index_weight_v, std::string atom_followed)
	{
		mt::vector_sr<rsp::spe_index_weight_t > out_spe_index_weight_v_tmp;

		boost::property_map<GraphContainer, vertex_index_t>::type vertex_id = get(vertex_index, graph);
		std::pair<vertex_iter, vertex_iter> vp;
		for (edge_range_t er = getEdges(); er.first != er.second; ++er.first) {//for
			//if found
			if (properties(*er.first).reaction_index == reaction_index) {//if
				std::size_t spe_index = get(vertex_id, target(*er.first, graph));
				//ignore duplicate elements
				out_spe_index_weight_v_tmp.insert_sr(
					std::make_pair(
						spe_index, //spe index
						properties(*er.first).s_coef_product //weight
					)
				);

			}//if
		}//for

		out_spe_index_weight_v = out_spe_index_weight_v_tmp;

		return true;
	}

	bool superReactionNetwork::search_for_out_spe_v3(rsp::index_int_t reaction_index, std::vector<rsp::spe_index_weight_t>& out_spe_index_weight_v, std::string atom_followed)
	{
		mt::vector_sr<rsp::spe_index_weight_t > out_spe_index_weight_v_tmp;

		boost::property_map<GraphContainer, vertex_index_t>::type vertex_id = get(vertex_index, graph);
		std::pair<vertex_iter, vertex_iter> vp;
		for (edge_range_t er = getEdges(); er.first != er.second; ++er.first) {//for
			//if found
			if (properties(*er.first).reaction_index == reaction_index) {//if
				std::size_t spe_index = get(vertex_id, target(*er.first, graph));
				//ignore duplicate elements
				out_spe_index_weight_v_tmp.insert_sr(
					std::make_pair(
						spe_index, //spe index
						this->species_network_v[spe_index].spe_component[atom_followed] //weight
					)
				);

			}//if
		}//for

		out_spe_index_weight_v = out_spe_index_weight_v_tmp;

		return true;
	}

	bool superReactionNetwork::set_reaction_out_spe_info(std::string atom_followed)
	{
		std::vector<rsp::spe_index_weight_t > out_spe_index_weight_v_tmp;
		for (std::size_t i = 0; i < this->reaction_network_v.size(); ++i) {
			out_spe_index_weight_v_tmp.clear();
			//		//version 1
			search_for_out_spe_v1(i, out_spe_index_weight_v_tmp, atom_followed);
			////		//version 2
			//		search_for_out_spe_v2(i, out_spe_index_weight_v_tmp, atom_followed);
			//		//version 3
			////		search_for_out_spe_v3(i, out_spe_index_weight_v_tmp, atom_followed);
			this->reaction_network_v[i].out_spe_index_weight_v_map[atom_followed] = out_spe_index_weight_v_tmp;
		}

		return true;
	}

	void superReactionNetwork::set_reaction_out_spe_info()
	{
		for (auto x : this->element_v)
			this->set_reaction_out_spe_info(x.ele_name);
		//super atom
		this->set_reaction_out_spe_info(rnk_pt.get<std::string>("pathway.super_atom"));
	}

	void superReactionNetwork::set_out_spe_index_branching_ratio_map_map(std::string atom_followed)
	{

		for (std::size_t r_index = 0; r_index < this->reaction_network_v.size(); ++r_index) {
			double prob_total = 0.0;
			//calcualte out spe total weight for a reaction
			for (std::size_t i = 0; i < reaction_network_v[r_index].out_spe_index_weight_v_map[atom_followed].size(); ++i) {
				prob_total += reaction_network_v[r_index].out_spe_index_weight_v_map[atom_followed][i].second;
			}
			//calculate the fraction
			for (std::size_t i = 0; i < reaction_network_v[r_index].out_spe_index_weight_v_map[atom_followed].size(); ++i) {
				reaction_network_v[r_index].out_spe_index_branching_ratio_map_map[atom_followed]
					[reaction_network_v[r_index].out_spe_index_weight_v_map[atom_followed][i].first] =
					reaction_network_v[r_index].out_spe_index_weight_v_map[atom_followed][i].second / prob_total;
			}
		}

	}

	void superReactionNetwork::set_out_spe_index_branching_ratio_map_map()
	{
		for (auto x : this->element_v)
			this->set_out_spe_index_branching_ratio_map_map(x.ele_name);
		//super atom
		this->set_out_spe_index_branching_ratio_map_map(rnk_pt.get<std::string>("pathway.super_atom"));
	}


	void superReactionNetwork::print_network(std::string filename)
	{
		std::ofstream fn((this->cwd + filename).c_str());
		fn << "Id,Label,Size\n";
		for (std::size_t i = 0; i < this->species_network_v.size(); ++i) {
			//		fn<<i<<","<<species_network_v[i].spe_name<<","<<"1.0"<<std::endl;
			fn << i << "," << species_network_v[i].spe_name << std::endl;
		}

		std::ofstream fe((this->cwd + std::string("/output/edge.csv")).c_str());
		fe << "Source,Target,Label,Weight\n";

		//this->update_reaction_rate(target_time);
		std::pair<vertex_iter, vertex_iter> vp;
		boost::property_map<GraphContainer, vertex_index_t>::type vertex_id = get(vertex_index, graph);
		for (edge_range_t er = getEdges(); er.first != er.second; ++er.first) {
			fe << get(vertex_id, source(*er.first, graph)) << "," << get(vertex_id, target(*er.first, graph)) << ",";
			fe << properties(*er.first).edge_index << "," << this->reaction_network_v[properties(*er.first).reaction_index].reaction_rate << std::endl;
		}

		fe.clear(); fe.close();
	}

	void superReactionNetwork::print()
	{

		////test random number generator
		//std::ofstream fout((this->cwd + std::string("/output/random.csv")).c_str());
		//for (size_t i = 0; i < 1000; ++i) {
		//	fout << std::setprecision(10) << this->rand->random01() << std::endl;
		//}
		//fout.clear(); fout.close();

		////	for(size_t i=0; i<10; ++i){
		////		std::cout<<random_min_max(generator, 0.0, 10.0)<<std::endl;
		////	}

		//std::cout << "Current path is : " << cwd << std::endl;
		//for (rsp::spe_name_index_map_t::const_iterator itr = spe_name_index_map.begin(); itr != spe_name_index_map.end(); ++itr) {
		//	std::cout << "spe name:\t" << itr->first << "\tindex:\t" << itr->second << std::endl;
		//}



		////species info
		//for (size_t i = 0; i < species_network_v.size(); ++i) {
		//	std::cout << species_network_v[i].spe_name << "\n";
		//	for (rsp::spe_component_t::const_iterator itr = species_network_v[i].spe_component.begin(); itr != species_network_v[i].spe_component.end(); ++itr)
		//		std::cout << "\t-->" << itr->first << "\t" << itr->second << std::endl;
		//}


		//std::cout << "\nnum_vertices: " << num_vertices << std::endl;
		//std::cout << "num_edges: " << num_edges << std::endl;

		//boost::property_map<GraphContainer, vertex_index_t>::type vertex_id = get(vertex_index, graph);

		//std::cout << "\nvertices(g) = \n";
		//std::pair<vertex_iter, vertex_iter> vp;
		//for (vp = vertices(graph); vp.first != vp.second; ++vp.first) {//for
		//	std::cout << get(vertex_id, *vp.first) << " ";
		//	std::cout << properties(*vp.first).vertex << std::endl;
		//}//for
		//std::cout << std::endl;
		////vertices properties


		//std::cout << "edges(g) = \n";

		//for (edge_range_t er = getEdges(); er.first != er.second; ++er.first) {
		//	std::cout << "(" << get(vertex_id, source(*er.first, graph))
		//		<< "," << get(vertex_id, target(*er.first, graph)) << ") ";
		//	std::cout << properties(*er.first).edge_index << " " << properties(*er.first).reaction_index <<
		//		" " << properties(*er.first).s_coef_reactant << " " << properties(*er.first).s_coef_product << std::endl;
		//}

		//for (edge_range_t er = getEdges(); er.first != er.second; ++er.first) {
		//	std::cout << "(" << species_network_v[get(vertex_id, source(*er.first, graph))].spe_name
		//		<< "," << species_network_v[get(vertex_id, target(*er.first, graph))].spe_name << ") ";
		//	std::cout << properties(*er.first).edge_index << " " << properties(*er.first).reaction_index <<
		//		" " << properties(*er.first).s_coef_reactant << " " << properties(*er.first).s_coef_product << std::endl;
		//}


		/*
		 * generate Mathematica Format graph input file
		 * DirectedEdge["A", "B"]
		 */

		 //	for (edge_range_t er=getEdges(); er.first!=er.second; ++er.first){
		 //		std::cout << "DirectedEdge[\"" <<species_network_v[get(vertex_id, source(*er.first, graph))].spe_name
		 //														<< "\", \"" << species_network_v[get(vertex_id, target(*er.first, graph))].spe_name << "\"]"<<std::endl;
		 //	}


		 //edge_index_map=get(edge_index, graph);
		 //test edge_index and edge number inside graph
		 //edge_iterator is not random access iterator, doesn't support offset dereference operator a[n]
		 //only support a++ or a-- and dereference *a operator
		//std::cout << "MISC edge index:" << std::endl;
		////	edge_range_t er_t= getEdges();
		//size_t edge_index_t = 2;
		//edge_iter iter_e = this->edge_index_to_edge_iterator[edge_index_t];
		//std::cout << "(" << get(vertex_id, source(*iter_e, graph))
		//	<< "," << get(vertex_id, target(*iter_e, graph)) << ") ";
		//std::cout << properties(*iter_e).edge_index << std::endl;



		////species out reaction info
		//for (std::size_t i = 0; i < species_network_v.size(); ++i) {
		//	std::cout << i << "\t" << species_network_v[i].reaction_k_index_s_coef_v.size() << std::endl;
		//}

		////reaction out spe info
		//std::string atom_followed("H");
		//for (std::size_t i = 0; i < reaction_network_v.size(); ++i) {
		//	std::cout << i << "\t" << reaction_network_v[i].out_spe_index_weight_v_map[atom_followed].size() << std::endl;
		//}

		////std::cout << "haha:\t" << std::endl;
		//for (std::size_t i = 0; i < species_network_v[2].reaction_k_index_s_coef_v.size(); ++i) {
		//	std::cout << species_network_v[2].reaction_k_index_s_coef_v[i].first << "\t" << species_network_v[2].reaction_k_index_s_coef_v[i].second << std::endl;
		//}


		////reaction rates test
		//std::cout << "reaction rates:\n";
		//for (auto x : reaction_network_v) {
		//	std::cout << x.reaction_rate << std::endl;
		//}

		//for (size_t i = 0; i < this->species_network_v.size(); i++)
		//{
		//	std::cout << i << "\t" << this->species_network_v[i].spe_name << ", ";
		//	for (auto x : this->species_network_v[i].reaction_k_index_s_coef_v)
		//		std::cout << x.first << "\t" << x.second << ", ";
		//	std::cout << "\n";
		//}

		rsp::reactionNetwork_chemkin_index_map_t reactionNetwork_chemkin_index_map;
		rsp::relationshipParser::read_reactionNetwork_chemkin_index_map(reactionNetwork_chemkin_index_map, this->cwd + "/input/chem.out");

		std::string atom_followed("O");
		std::ofstream fout((this->cwd + std::string("/output/species_reaction_index_helper_") + atom_followed + std::string(".csv")).c_str());

		for (size_t i = 0; i < this->species_network_v.size(); i++)
		{
			fout << this->species_network_v[i].spe_name << "\n";
			for (size_t j = 0; j < this->species_network_v.size(); j++)
			{
				if (i < j) {
					auto p = this->get_R_matrix_element(atom_followed, i, j);
					if (p.size() > 0) {
						fout << "\t-->" << this->species_network_v[j].spe_name << "\n\t\t";
						for (size_t k = 0; k < p.size(); k++)
						{
							//this is edge index, one-step path, a path with only one reaction
							//edge index to reaction index
							auto iter_e = edge_index_to_edge_iterator[p[k][0]];
							auto reaction_index = this->properties(*iter_e).reaction_index;
							//fout << this->reaction_network_v[reaction_index].reaction_name;

							auto reaction_index_in_paper = reactionNetwork_chemkin_index_map[reaction_index].front();

							if (reaction_index_in_paper > 0)
								fout << "R" << abs(reaction_index_in_paper) - 1;
							else
								/*fout << "R" << -1 * (abs(reaction_index_in_paper) - 1);*/
								fout << "R" << abs(reaction_index_in_paper) - 1 << "*";

							if (k != p.size() - 1)
								fout << ",";
							else
								fout << "\n";
						}
					}
				}
			}
		}

		fout.clear(); fout.close();


	}

	void superReactionNetwork::print_initial_spe_label_json(std::string filename) const
	{
		std::ofstream fout((cwd + filename).c_str());
		fout << "{\n";

		for (std::size_t i = 0; i < species_network_v.size(); ++i)
		{
			std::size_t NofAtoms = 0;
			std::size_t TypesofAtoms = 0;
			fout << "\"" << i << "\":\n" << "{";
			fout << "\"name\":" << "\"" << species_network_v[i].spe_name << "\",\n";
			fout << "\"structure\":" << "\"" << species_network_v[i].spe_name << "\",\n";

			for (rsp::spe_component_t::const_iterator itr = species_network_v[i].spe_component.begin(); itr != species_network_v[i].spe_component.end(); ++itr) {
				NofAtoms += itr->second;
				if (itr->second != 0)
					TypesofAtoms += 1;
			}
			fout << "\"TotalNofAtoms\":" << "\"" << NofAtoms << "\",\n";
			fout << "\"TypesofAtoms\":" << "\"" << TypesofAtoms << "\",\n";

			for (rsp::spe_component_t::const_iterator itr = species_network_v[i].spe_component.begin(); itr != species_network_v[i].spe_component.end(); ++itr) {
				if (itr->second != 0)
					fout << "\"" << itr->first << "\"" << ":" << "\"" << itr->second << "\",\n";
			}

			fout << "\"TypesofChemicalMoiety\":" << "\"" << "\",\n";

			for (std::size_t j = 0; j < NofAtoms; ++j) {
				fout << "\"\":\"\"";
				if (j != NofAtoms - 1)
					fout << ",\n";
			}
			fout << std::endl << "}";
			if (i != species_network_v.size() - 1)
				fout << ",";
			fout << "\n";
		}

		fout << "}\n";
		fout.close();

	}

	void superReactionNetwork::superReactionNetwork::print_initial_reaction_label_json(std::string filename) const
	{
	}

	void superReactionNetwork::set_atom_followed_edge_weight_at_time_v0(rsp::my_time_t t, std::string atom, double default_weight)
	{
		//1. update all reaction rates
		update_reaction_rate(t);
		//2. update edges weights
		Edge e;
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			//if the source and target vertex of the edge both contain the followed atom
			if ((species_network_v[source(e, this->graph)].spe_component[atom] != 0) && (species_network_v[target(e, this->graph)].spe_component[atom] != 0))
				properties(e).edge_weight = species_network_v[source(e, this->graph)].spe_component[atom] *
				reaction_network_v[properties(e).reaction_index].reaction_rate;

			else
				properties(e).edge_weight = default_weight;
		}
		this->delta_edge_weight = 0.0;
		this->minimum_edge_weight = 0.0;

	}

	void superReactionNetwork::set_atom_followed_edge_weight_at_time_v1(rsp::my_time_t t, std::string atom, double default_weight)
	{
		//1. update all reaction rates
		update_reaction_rate(t);
		//2. update edges weights
		Edge e;
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			//if the source and target vertex of the edge both contain the followed atom
			if ((species_network_v[source(e, this->graph)].spe_component[atom] != 0) && (species_network_v[target(e, this->graph)].spe_component[atom] != 0))
				properties(e).edge_weight = species_network_v[source(e, this->graph)].spe_component[atom] *
				reaction_network_v[properties(e).reaction_index].reaction_rate *
				-1.0;//take the negative of it

			else
				properties(e).edge_weight = default_weight;
		}

		//find the most negative value
		double min_t = 0.0;
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			if (properties(e).edge_weight < min_t)
				min_t = properties(e).edge_weight;
		}
		//std::cout << min_t << std::endl;
		//set the minimum_edge_weight to be the reverse of the minimum weight
		this->minimum_edge_weight = min_t;
		//don't want it to be zero, add the minimum edge weight
		this->delta_edge_weight = -2.0*this->minimum_edge_weight;
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			if (properties(e).edge_weight != default_weight)
				properties(e).edge_weight += this->delta_edge_weight;
		}

	}

	void superReactionNetwork::set_atom_followed_edge_weight_at_time_v2(rsp::my_time_t t, std::string atom, double default_weight)
	{
		//1. update all reaction rates
		update_reaction_rate(t);
		//2. update edges weights
		Edge e;
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			//if the source and target vertex of the edge both contain the followed atom
			if ((species_network_v[source(e, this->graph)].spe_component[atom] != 0) && (species_network_v[target(e, this->graph)].spe_component[atom] != 0)
				&& evaluate_spe_concentration_at_time(t, source(e, this->graph)) != 0.0)
			{
				//calculate seudo-first order rate constant
				properties(e).edge_weight = species_network_v[source(e, this->graph)].spe_component[atom] *
					reaction_network_v[properties(e).reaction_index].reaction_rate /
					evaluate_spe_concentration_at_time(t, source(e, this->graph)) *
					-1.0;//take the negative of it
			}

			else
				properties(e).edge_weight = default_weight;
		}

		//find the most negative value
		double min_t = 0.0;
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			if (properties(e).edge_weight < min_t)
				min_t = properties(e).edge_weight;
		}

		//set the minimum_edge_weight to be the reverse of the minimum weight
		this->minimum_edge_weight = min_t;
		//don't want it to be zero, add the minimum edge weight
		this->delta_edge_weight = -2.0*this->minimum_edge_weight;

		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			if (properties(e).edge_weight != default_weight)
				properties(e).edge_weight += this->delta_edge_weight;
		}

	}

	void superReactionNetwork::set_atom_followed_edge_weight_at_time_v3(rsp::my_time_t t, std::string atom, double default_weight)
	{
		//1. update all reaction rates
		update_reaction_rate(t);
		//2. update edges weights
		Edge e;
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			//if the source and target vertex of the edge both contain the followed atom
			if ((species_network_v[source(e, this->graph)].spe_component[atom] != 0) && (species_network_v[target(e, this->graph)].spe_component[atom] != 0)
				&& evaluate_spe_concentration_at_time(t, source(e, this->graph)) != 0.0)
			{
				//calculate seudo-first order rate constant
				properties(e).edge_weight = species_network_v[source(e, this->graph)].spe_component[atom] *
					reaction_network_v[properties(e).reaction_index].reaction_rate /
					evaluate_spe_concentration_at_time(t, source(e, this->graph));
			}

			else
				properties(e).edge_weight = 0.0;
		}
		//Vertex v;
		BGL_FORALL_VERTICES(v, this->graph, GraphContainer)
		{
			Edge e;
			double edge_weight_total = 0.0;
			//sum over all out edge weights
			BGL_FORALL_OUTEDGES(v, e, this->graph, GraphContainer)
			{
				edge_weight_total += properties(e).edge_weight;
			}
			if (edge_weight_total != 0) {//if
				//average
				std::cout << "-->\t";
				BGL_FORALL_OUTEDGES(v, e, this->graph, GraphContainer)
				{
					properties(e).edge_weight /= edge_weight_total;
					std::cout << properties(e).edge_weight << "\t";
				}
				std::cout << "<--" << std::endl;
			}//if
			else
			{
				std::cout << "-->\t";
				BGL_FORALL_OUTEDGES(v, e, this->graph, GraphContainer)
				{
					std::cout << properties(e).edge_weight << "\t";
				}
				std::cout << "<--" << std::endl;
			}
		}

		//
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			//edge weight not zero
			if (properties(e).edge_weight != 0)
			{
				properties(e).edge_weight *= -1.0;
			}
			//edge weight zero
			else
			{
				properties(e).edge_weight = default_weight;
			}

		}

		//find the most negative value
		double min_t = 0.0;
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			if (properties(e).edge_weight < min_t)
				min_t = properties(e).edge_weight;
		}

		//set the minimum_edge_weight to be the reverse of the minimum weight
		this->minimum_edge_weight = min_t;
		//don't want it to be zero, add the minimum edge weight
		this->delta_edge_weight = -2.0*this->minimum_edge_weight;

		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			if (properties(e).edge_weight != default_weight)
				properties(e).edge_weight += this->delta_edge_weight;
		}

	}

	void superReactionNetwork::set_atom_followed_edge_weight_at_time_v4(rsp::my_time_t t, std::string atom, double default_weight)
	{
		//1. update all reaction rates
		update_reaction_rate(t);
		//2. update edges weights
		Edge e;
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			//if the source and target vertex of the edge both contain the followed atom
			if ((species_network_v[source(e, this->graph)].spe_component[atom] != 0) && (species_network_v[target(e, this->graph)].spe_component[atom] != 0)
				&& evaluate_spe_concentration_at_time(t, source(e, this->graph)) != 0.0)
			{
				//calculate seudo-first order rate constant
				properties(e).edge_weight = species_network_v[source(e, this->graph)].spe_component[atom] *
					reaction_network_v[properties(e).reaction_index].reaction_rate /
					evaluate_spe_concentration_at_time(t, source(e, this->graph));
			}

			else
				properties(e).edge_weight = 0.0;
		}
		//Vertex v;
		BGL_FORALL_VERTICES(v, this->graph, GraphContainer)
		{
			Edge e;
			double edge_weight_total = 0.0;
			//sum over all out edge weights
			BGL_FORALL_OUTEDGES(v, e, this->graph, GraphContainer)
			{
				edge_weight_total += properties(e).edge_weight;
			}
			if (edge_weight_total != 0) {//if
				//average
				std::cout << "-->\t";
				BGL_FORALL_OUTEDGES(v, e, this->graph, GraphContainer)
				{
					properties(e).edge_weight /= edge_weight_total;
					//take the natural log
					if (properties(e).edge_weight != 0)
						properties(e).edge_weight = log(properties(e).edge_weight);
					else
						//out self defined magic number
						properties(e).edge_weight = boost::math::constants::pi<double>();

					std::cout << properties(e).edge_weight << "\t";
				}
				std::cout << "<--" << std::endl;
			}//if
			else
			{
				std::cout << "-->\t";
				BGL_FORALL_OUTEDGES(v, e, this->graph, GraphContainer)
				{
					//assign a magic number
					properties(e).edge_weight = boost::math::constants::pi<double>();
					std::cout << properties(e).edge_weight << "\t";
				}
				std::cout << "<--" << std::endl;
			}
		}

		//since log(brancing ratio)<=0, only take the negative without shifting will do the job
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			//edge weight not our magic number
			if (properties(e).edge_weight != boost::math::constants::pi<double>())
			{
				properties(e).edge_weight *= -1.0;
			}
			//edge weight is actually our magic number
			else
			{
				properties(e).edge_weight = default_weight;
			}

		}

		this->delta_edge_weight = 0.0;
		this->minimum_edge_weight = 0.0;

	}

	double superReactionNetwork::dijkstra_w_set(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in, std::string atom_followed)
	{
		//a map keep track of the node's direct parent when update nodes distance
		std::map<vertex_t, vertex_t > vertex_parent_map;

		update_reaction_rate(time_in);
		vector<double> min_distance(this->species_network_v.size(), std::numeric_limits<double>::max());
		min_distance[source_v] = 0.0;
		//the second is the vertex, the first is the current distance from the source to the vertex
		std::set<std::pair<double, rsp::index_int_t> > active_vertices;

		active_vertices.insert({ 0.0, source_v });

		while (!active_vertices.empty())
		{
			auto where = active_vertices.begin()->second;
			std::cout << "where:\t" << species_network_v[where].spe_name << "\tdistance:\t" << min_distance[where] << std::endl;

			if (where == target_v) {
				std::cout << "min_distance:\t" << min_distance[where] << std::endl;

				print_dijkstra_algorithm_shortest_path(vertex_parent_map, source_v, target_v);
				std::cout << "\n";

				return min_distance[where];
			}
			active_vertices.erase(active_vertices.begin());

			//assume that each reaction just has one product
			for (auto edge : species_network_v[where].reaction_k_index_s_coef_v)
			{
				auto next_out_vertex = reaction_network_v[edge.first].out_spe_index_weight_v_map[atom_followed][0].first;

				if (min_distance[next_out_vertex] > min_distance[where] + reaction_network_v[edge.first].reaction_rate)
				{
					vertex_parent_map[next_out_vertex] = where;
					active_vertices.erase({ min_distance[next_out_vertex], next_out_vertex });
					min_distance[next_out_vertex] = min_distance[where] + reaction_network_v[edge.first].reaction_rate;
					active_vertices.insert({ min_distance[next_out_vertex], next_out_vertex });
				}

			}

		}

		return 0.0;

	}

	double superReactionNetwork::dijkstra_w_priority_queue(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in, std::string atom_followed)
	{
		//a map keep track of the node's direct parent when update nodes distance
		std::map<vertex_t, vertex_t > vertex_parent_map;

		update_reaction_rate(time_in);
		vector<double> min_distance(this->species_network_v.size(), std::numeric_limits<double>::max());
		//Determines whether the node has been visited or not
		vector<bool> visited(this->species_network_v.size(), false);

		//Custom Comparator for Determining priority for priority queue (shortest edge comes first)
		struct prioritize {
		public: bool operator ()(std::pair<vertex_t, double>&p1, std::pair<vertex_t, double>&p2) const
		{
			return p1.second > p2.second;
		}
		};
		//Priority queue to store vertex,weight pairs
		priority_queue<pair<vertex_t, double>, vector<pair<vertex_t, double> >, prioritize> pq;
		pq.push(make_pair(source_v, min_distance[source_v] = 0.0)); //Pushing the source with distance from itself as 0

		while (!pq.empty())
		{
			auto where = pq.top().first;
			pq.pop();
			if (visited[where] == true)
				continue;
			visited[where] = true;

			std::cout << "where:\t" << species_network_v[where].spe_name << "\tdistance:\t" << min_distance[where] << std::endl;

			if (where == target_v) {
				std::cout << "min_distance:\t" << min_distance[where] << std::endl;

				print_dijkstra_algorithm_shortest_path(vertex_parent_map, source_v, target_v);
				std::cout << "\n";

				return min_distance[where];
			}


			//assume that each reaction just has one product
			for (auto edge : species_network_v[where].reaction_k_index_s_coef_v)
			{
				auto next_out_vertex = reaction_network_v[edge.first].out_spe_index_weight_v_map[atom_followed][0].first;

				if (!visited[next_out_vertex] && min_distance[next_out_vertex] > min_distance[where] + reaction_network_v[edge.first].reaction_rate)
				{
					vertex_parent_map[next_out_vertex] = where;

					min_distance[next_out_vertex] = min_distance[where] + reaction_network_v[edge.first].reaction_rate;
					pq.push(std::make_pair(next_out_vertex, min_distance[next_out_vertex]));
				}

			}

		}

		return 0.0;
	}

	double superReactionNetwork::dijkstra_boost(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in)
	{
		//update_reaction_rate(time_in);
		//create a descriptor for the source node
		Vertex s = vertex(source_v, this->graph);
		Vertex t = vertex(target_v, this->graph);

		//create the property_map from edges to weights
		//doesn't work for bundled properties
		//property_map<GraphContainer, edge_weight_t>::type weightmap = get(edge_weight, this->graph);
		std::map<Edge, double> weightMap;
		boost::associative_property_map< std::map<Edge, double> > property_weight_map(weightMap);
		//auto weight_map = get(&EdgeProperties_graph::edge_index, this->graph);
		Edge e;
		//indexing the edges
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			put(property_weight_map, e, properties(e).edge_weight);
			//put(property_weight_map, e, this->reaction_network_v[properties(e).reaction_index].reaction_rate);
		}

		//create dijkstra visitor
		//predecessor vector and corresponding edges vector
		std::vector<Vertex> predecessors_v; std::vector<Edge> predecessors_edge_v;
		my_visitor vis(graph, predecessors_v, predecessors_edge_v);

		//create vectors to store the predecessors (p)
		std::vector<Vertex> predecessors(this->num_vertices);
		//create vectors to store the distances from the root (d)
		std::vector<double> distance(this->num_vertices);
		//evaluate dijkstra on graph g with source s, predecessor_map p and distance_map d
		//note that predecessor_map(..).distance_map(..) is a bgl_named_params<P, T, R>, so a named parameter
		//http://stackoverflow.com/questions/31145082/bgl-dijkstra-shortest-paths-with-bundled-properties
		//http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/bundles.html
		//you can pass it to dijkstra using direct or named params. Let's do the simplest
		boost::dijkstra_shortest_paths(this->graph, s, boost::no_named_parameters().weight_map(property_weight_map)
			.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
			.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(vertex_index, this->graph)))
			.visitor(vis)
		);

		std::cout << "distances from start vertex:" << std::endl;
		for (auto v = vertices(this->graph).first; v != vertices(this->graph).second; ++v)
		{
			std::cout << "distance(" << properties(*v).vertex << ") = "
				<< distance[*v] << std::endl;

		}
		std::cout << std::endl;

		//Extract the shortest path from s to t.
		typedef std::vector<Edge> path_t;
		path_t path;

		for (auto u = predecessors[t]; u != t; t = u, u = predecessors[t])
		{
			std::pair<Edge, bool> edge_pair = boost::edge(u, t, this->graph);
			path.push_back(edge_pair.first);
		}

		std::cout << std::endl;
		std::cout << "Shortest Path from s to t:" << std::endl;
		for (path_t::reverse_iterator riter = path.rbegin(); riter != path.rend(); ++riter)
		{
			Vertex u_tmp = boost::source(*riter, this->graph);
			Vertex v_tmp = boost::target(*riter, this->graph);
			Edge   e_tmp = boost::edge(u_tmp, v_tmp, this->graph).first;

			std::cout << "  " << species_network_v[properties(u_tmp).vertex].spe_name << " -> "
				<< species_network_v[properties(v_tmp).vertex].spe_name << "    (weight: "
				<< reaction_network_v[properties(e_tmp).reaction_index].reaction_rate << ")" << std::endl;
		}
		t = vertex(target_v, this->graph);
		std::cout << distance[t] << std::endl;
		return distance[t];

	}

	double superReactionNetwork::dijkstra_boost_v2(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in)
	{
		//update_reaction_rate(time_in);
		//create a descriptor for the source node
		Vertex s = vertex(source_v, this->graph);
		Vertex t = vertex(target_v, this->graph);

		//create the property_map from edges to weights
		//doesn't work for bundled properties
		//property_map<GraphContainer, edge_weight_t>::type weightmap = get(edge_weight, this->graph);
		std::map<Edge, double> weightMap;
		boost::associative_property_map< std::map<Edge, double> > property_weight_map(weightMap);
		//auto weight_map = get(&EdgeProperties_graph::edge_index, this->graph);
		Edge e;
		//indexing the edges
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			put(property_weight_map, e, properties(e).edge_weight);
			//put(property_weight_map, e, this->reaction_network_v[properties(e).reaction_index].reaction_rate);
		}

		//create dijkstra visitor
		//create vectors to store the predecessors (p)
		//and corresponding edges vector (e), edges on the shortest path (e)
		std::vector<Vertex> predecessors; std::vector<Edge> predecessors_edges_v;
		my_visitor vis(graph, predecessors, predecessors_edges_v);

		//create vectors to store the distances from the root (d)
		std::vector<double> distance(this->num_vertices);
		//evaluate dijkstra on graph g with source s, predecessor_map p and distance_map d
		//note that predecessor_map(..).distance_map(..) is a bgl_named_params<P, T, R>, so a named parameter
		//http://stackoverflow.com/questions/31145082/bgl-dijkstra-shortest-paths-with-bundled-properties
		//http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/bundles.html
		//you can pass it to dijkstra using direct or named params. Let's do the simplest
		boost::dijkstra_shortest_paths(this->graph, s, boost::no_named_parameters().weight_map(property_weight_map)
			.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
			//.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(vertex_index, this->graph)))
			.visitor(vis)
		);

		std::cout << "distances from start vertex:" << std::endl;
		for (auto v = vertices(this->graph).first; v != vertices(this->graph).second; ++v)
		{
			std::cout << "distance(" << properties(*v).vertex << ") = "
				<< distance[*v] << std::endl;

		}
		std::cout << std::endl;

		//Extract the shortest path from s to t.
		typedef std::vector<Edge> path_t;
		path_t path;

		for (auto u = predecessors[t]; u != t; t = u, u = predecessors[t])
		{
			std::pair<Edge, bool> edge_pair = boost::edge(u, t, this->graph);
			path.push_back(edge_pair.first);
		}

		std::cout << std::endl;
		std::cout << "Shortest Path from s to t:" << std::endl;
		for (path_t::reverse_iterator riter = path.rbegin(); riter != path.rend(); ++riter)
		{
			auto u_tmp = boost::source(*riter, this->graph);
			auto v_tmp = boost::target(*riter, this->graph);
			auto e_tmp = predecessors_edges_v[v_tmp];

			std::cout << "  " << species_network_v[properties(u_tmp).vertex].spe_name << " -> "
				<< species_network_v[properties(v_tmp).vertex].spe_name << "    (weight: "
				<< reaction_network_v[properties(e_tmp).reaction_index].reaction_rate << ")" << std::endl;
		}
		t = vertex(target_v, this->graph);
		std::cout << distance[t] << std::endl;
		return distance[t];
	}

	void superReactionNetwork::dijkstra_boost(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in,
		std::vector<Vertex>& predecessors, std::vector<double>& distance)
	{
		predecessors.resize(this->num_vertices);
		distance.resize(this->num_vertices);

		//update_reaction_rate(time_in);
		//create a descriptor for the source node
		Vertex s = vertex(source_v, this->graph);
		Vertex t = vertex(target_v, this->graph);

		//create the property_map from edges to weights
		//doesn't work for bundled properties
		//property_map<GraphContainer, edge_weight_t>::type weightmap = get(edge_weight, g);
		std::map<Edge, double> weightMap;
		boost::associative_property_map< std::map<Edge, double> > property_weight_map(weightMap);
		Edge e;
		//indexing the edges
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			put(property_weight_map, e, properties(e).edge_weight);
			//put(property_weight_map, e, this->reaction_network_v[properties(e).reaction_index].reaction_rate);
		}

		//evaluate dijkstra on graph g with source s, predecessor_map p and distance_map d
		//note that predecessor_map(..).distance_map(..) is a bgl_named_params<P, T, R>, so a named parameter
		//http://stackoverflow.com/questions/31145082/bgl-dijkstra-shortest-paths-with-bundled-properties
		//http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/bundles.html
		//you can pass it to dijkstra using direct or named params. Let's do the simplest
		//boost::dijkstra_shortest_paths(this->graph, s, boost::no_named_parameters().weight_map(property_weight_map)
		//	.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
		//	.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(vertex_index, this->graph))));
		boost::dijkstra_shortest_paths(this->graph, s, weight_map(property_weight_map)
			.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
			.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(vertex_index, this->graph))));

		//std::cout << "distances from start vertex:" << std::endl;
		//for (auto v = vertices(this->graph).first; v != vertices(this->graph).second; ++v)
		//{
		//	std::cout << "distance(" << properties(*v).vertex << ") = "
		//		<< distance[*v] << std::endl;
		//}
		//std::cout << std::endl;

		////Extract the shortest path from s to t.
		//typedef std::vector<Edge> path_t;
		//path_t path;

		//for (auto u = predecessors[t]; u != t; t = u, u = predecessors[t])
		//{
		//	std::pair<Edge, bool> edge_pair = boost::edge(u, t, this->graph);
		//	path.push_back(edge_pair.first);
		//}

		//std::cout << "Shortest Path from s to t:" << std::endl;
		//for (path_t::reverse_iterator riter = path.rbegin(); riter != path.rend(); ++riter)
		//{
		//	Vertex u_tmp = boost::source(*riter, this->graph);
		//	Vertex v_tmp = boost::target(*riter, this->graph);
		//	Edge   e_tmp = boost::edge(u_tmp, v_tmp, this->graph).first;

		//	std::cout << "  " << species_network_v[properties(u_tmp).vertex].spe_name << " -> "
		//		<< species_network_v[properties(v_tmp).vertex].spe_name << "    (weight: "
		//		<< reaction_network_v[properties(e_tmp).reaction_index].reaction_rate << ")" << std::endl;
		//}
		//t = vertex(target_v, this->graph);
		//std::cout << distance[t] << std::endl;
		//return distance[t];
	}

	void superReactionNetwork::dijkstra_boost_v2(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in, std::vector<Vertex>& predecessors, std::vector<Edge>& predecessors_edges, std::vector<double>& distance)
	{
		//predecessors.resize(this->num_vertices);
		distance.resize(this->num_vertices);

		//update_reaction_rate(time_in);
		//create a descriptor for the source node
		Vertex s = vertex(source_v, this->graph);
		Vertex t = vertex(target_v, this->graph);

		//create the property_map from edges to weights
		//doesn't work for bundled properties
		//property_map<GraphContainer, edge_weight_t>::type weightmap = get(edge_weight, g);
		std::map<Edge, double> weightMap;
		boost::associative_property_map< std::map<Edge, double> > property_weight_map(weightMap);
		Edge e;
		//indexing the edges
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			put(property_weight_map, e, properties(e).edge_weight);
			//put(property_weight_map, e, this->reaction_network_v[properties(e).reaction_index].reaction_rate);
		}

		//create dijkstra visitor
		//create vectors to store the predecessors (p)
		//and corresponding edges vector (e), edges on the shortest path (e)
		my_visitor vis(graph, predecessors, predecessors_edges);

		//evaluate dijkstra on graph g with source s, predecessor_map p and distance_map d
		//note that predecessor_map(..).distance_map(..) is a bgl_named_params<P, T, R>, so a named parameter
		//http://stackoverflow.com/questions/31145082/bgl-dijkstra-shortest-paths-with-bundled-properties
		//http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/bundles.html
		//you can pass it to dijkstra using direct or named params. Let's do the simplest
		//boost::dijkstra_shortest_paths(this->graph, s, boost::no_named_parameters().weight_map(property_weight_map)
		//	.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
		//	.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(vertex_index, this->graph))));
		boost::dijkstra_shortest_paths(this->graph, s, weight_map(property_weight_map)
			.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
			//.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(vertex_index, this->graph)))
			.visitor(vis)
		);
	}

	void superReactionNetwork::bellman_ford_boost(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in, std::vector<Graph_bundled::Vertex>& predecessors, std::vector<double>& distance)
	{
		predecessors.resize(this->num_vertices);
		distance.resize(this->num_vertices);

		//update_reaction_rate(time_in);
		//create a descriptor for the source node
		Vertex s = vertex(source_v, this->graph);
		Vertex t = vertex(target_v, this->graph);

		//create the property_map from edges to weights
		//doesn't work for bundled properties
		//property_map<GraphContainer, edge_weight_t>::type weightmap = get(edge_weight, g);
		std::map<Edge, double> weightMap;
		boost::associative_property_map< std::map<Edge, double> > property_weight_map(weightMap);
		Edge e;
		//indexing the edges
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			put(property_weight_map, e, properties(e).edge_weight);
			//put(property_weight_map, e, this->reaction_network_v[properties(e).reaction_index].reaction_rate);
		}

		//evaluate dijkstra on graph g with source s, predecessor_map p and distance_map d
		//note that predecessor_map(..).distance_map(..) is a bgl_named_params<P, T, R>, so a named parameter
		//http://stackoverflow.com/questions/31145082/bgl-dijkstra-shortest-paths-with-bundled-properties
		//http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/bundles.html
		//you can pass it to dijkstra using direct or named params. Let's do the simplest
		//boost::dijkstra_shortest_paths(this->graph, s, weight_map(property_weight_map)
		//	.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
		//	.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(vertex_index, this->graph))));
		//http://www.boost.org/doc/libs/1_60_0/libs/graph/doc/bellman_ford_shortest.html

		boost::bellman_ford_shortest_paths(this->graph, boost::num_vertices(graph), weight_map(property_weight_map)
			.root_vertex(s)
			.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
			.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(vertex_index, this->graph))));

		//std::cout << "distances from start vertex:" << std::endl;
		//for (auto v = vertices(this->graph).first; v != vertices(this->graph).second; ++v)
		//{
		//	std::cout << "distance(" << properties(*v).vertex << ") = "
		//		<< distance[*v] << std::endl;
		//}
		//std::cout << std::endl;

		////Extract the shortest path from s to t.
		//typedef std::vector<Edge> path_t;
		//path_t path;

		//for (auto u = predecessors[t]; u != t; t = u, u = predecessors[t])
		//{
		//	std::pair<Edge, bool> edge_pair = boost::edge(u, t, this->graph);
		//	path.push_back(edge_pair.first);
		//}

		//std::cout << "Shortest Path from s to t:" << std::endl;
		//for (path_t::reverse_iterator riter = path.rbegin(); riter != path.rend(); ++riter)
		//{
		//	Vertex u_tmp = boost::source(*riter, this->graph);
		//	Vertex v_tmp = boost::target(*riter, this->graph);
		//	Edge   e_tmp = boost::edge(u_tmp, v_tmp, this->graph).first;

		//	std::cout << "  " << species_network_v[properties(u_tmp).vertex].spe_name << " -> "
		//		<< species_network_v[properties(v_tmp).vertex].spe_name << "    (weight: "
		//		<< reaction_network_v[properties(e_tmp).reaction_index].reaction_rate << ")" << std::endl;
		//}
		//t = vertex(target_v, this->graph);
		//std::cout << distance[t] << std::endl;
		//return distance[t];
	}

	void superReactionNetwork::bellman_ford_boost_v2(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in, std::vector<Vertex>& predecessors, std::vector<Edge>& predecessors_edges, std::vector<double>& distance)
	{
		//predecessors.resize(this->num_vertices);
		distance.resize(this->num_vertices);

		//update_reaction_rate(time_in);
		//create a descriptor for the source node
		Vertex s = vertex(source_v, this->graph);
		Vertex t = vertex(target_v, this->graph);

		//create the property_map from edges to weights
		//doesn't work for bundled properties
		//property_map<GraphContainer, edge_weight_t>::type weightmap = get(edge_weight, g);
		std::map<Edge, double> weightMap;
		boost::associative_property_map< std::map<Edge, double> > property_weight_map(weightMap);
		Edge e;
		//indexing the edges
		BGL_FORALL_EDGES(e, this->graph, GraphContainer)
		{
			put(property_weight_map, e, properties(e).edge_weight);
			//put(property_weight_map, e, this->reaction_network_v[properties(e).reaction_index].reaction_rate);
		}

		//create dijkstra visitor
		//create vectors to store the predecessors (p)
		//and corresponding edges vector (e), edges on the shortest path (e)
		my_visitor vis(graph, predecessors, predecessors_edges);

		//evaluate dijkstra on graph g with source s, predecessor_map p and distance_map d
		//note that predecessor_map(..).distance_map(..) is a bgl_named_params<P, T, R>, so a named parameter
		//http://stackoverflow.com/questions/31145082/bgl-dijkstra-shortest-paths-with-bundled-properties
		//http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/bundles.html
		//you can pass it to dijkstra using direct or named params. Let's do the simplest
		//boost::dijkstra_shortest_paths(this->graph, s, weight_map(property_weight_map)
		//	.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
		//	.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(vertex_index, this->graph))));
		//http://www.boost.org/doc/libs/1_60_0/libs/graph/doc/bellman_ford_shortest.html

		boost::bellman_ford_shortest_paths(this->graph, boost::num_vertices(graph), weight_map(property_weight_map)
			.root_vertex(s)
			.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
			//.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(vertex_index, this->graph)))
			.visitor(vis)
		);
	}

	void superReactionNetwork::dijkstra_reverse_graph_boost(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in,
		std::vector<Graph_bundled::reverse_vertex_t>& predecessors, std::vector<double>& distance)
	{

		//dijkstra on reversed graph
		//https://groups.google.com/forum/#!msg/boost-list/j8ufXCgqUdE/2SlbjLD_3pIJ
		predecessors.resize(this->num_vertices);
		distance.resize(this->num_vertices);

		reverse_graph_t reversed_graph(this->graph);

		//update_reaction_rate(time_in);
		//create a descriptor for the source node
		Vertex s = vertex(source_v, this->graph);
		Vertex t = vertex(target_v, this->graph);

		//create the property_map from edges to weights
		//doesn't work for bundled properties
		//property_map<GraphContainer, edge_weight_t>::type weightmap = get(edge_weight, g);
		std::map<reverse_edge_t, double> weightMap;
		boost::associative_property_map< std::map<reverse_edge_t, double> > property_weight_map(weightMap);

		Edge e;
		//indexing the edges
		BGL_FORALL_EDGES(e, reversed_graph, reverse_graph_t)
		{
			//one way to get the edge properties of the reversed graph
			typename property_map<reverse_graph_t, edge_properties_t>::type param = get(edge_properties, reversed_graph);
			put(property_weight_map, e, param[e].edge_weight);
			//put(property_weight_map, e, this->reaction_network_v[param[e].reaction_index].reaction_rate);
		}
		//evaluate dijkstra on graph g with source s, predecessor_map p and distance_map d
		//note that predecessor_map(..).distance_map(..) is a bgl_named_params<P, T, R>, so a named parameter
		//http://stackoverflow.com/questions/31145082/bgl-dijkstra-shortest-paths-with-bundled-properties
		//http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/bundles.html
		//you can pass it to dijkstra using direct or named params. Let's do the simplest
		boost::dijkstra_shortest_paths(reversed_graph, s, weight_map(property_weight_map)
			.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
			.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(vertex_index, this->graph))));

		//std::cout << "distances from start vertex:" << std::endl;
		//for (auto v = vertices(this->graph).first; v != vertices(this->graph).second; ++v)
		//{
		//	std::cout << "distance(" << properties(*v).vertex << ") = "
		//		<< distance[*v] << std::endl;

		//}

		////Extract the shortest path from s to t.
		//typedef std::vector<Edge> path_t;
		//path_t path;

		//for (auto u = predecessors[t]; u != t; t = u, u = predecessors[t])
		//{
		//	std::pair<Edge, bool> edge_pair = boost::edge(u, t, this->graph);
		//	path.push_back(edge_pair.first);
		//}

		//std::cout << std::endl;
		//std::cout << "Shortest Path from t to s:" << std::endl;
		//for (path_t::reverse_iterator riter = path.rbegin(); riter != path.rend(); ++riter)
		//{
		//	Vertex u_tmp = boost::source(*riter, this->graph);
		//	Vertex v_tmp = boost::target(*riter, this->graph);
		//	Edge   e_tmp = boost::edge(v_tmp, u_tmp, this->graph).first;

		//	std::cout << "  " << species_network_v[properties(u_tmp).vertex].spe_name << " -> "
		//		<< species_network_v[properties(v_tmp).vertex].spe_name << "    (weight: "
		//		<< reaction_network_v[properties(e_tmp).reaction_index].reaction_rate << ")" << std::endl;
		//}
		//t = vertex(target_v, this->graph);
		//std::cout << distance[t] << std::endl;
		////return distance[t];
	}

	void superReactionNetwork::dijkstra_reverse_graph_boost_v2(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in, std::vector<reverse_vertex_t>& predecessors_r, std::vector<reverse_edge_t> &predecessors_edges_r, std::vector<double>& distance)
	{
		//dijkstra on reversed graph
		//https://groups.google.com/forum/#!msg/boost-list/j8ufXCgqUdE/2SlbjLD_3pIJ
		//predecessors.resize(this->num_vertices);
		distance.resize(this->num_vertices);

		reverse_graph_t reversed_graph(this->graph);

		//update_reaction_rate(time_in);
		//create a descriptor for the source node
		Vertex s = vertex(source_v, this->graph);
		Vertex t = vertex(target_v, this->graph);

		//create the property_map from edges to weights
		//doesn't work for bundled properties
		//property_map<GraphContainer, edge_weight_t>::type weightmap = get(edge_weight, g);
		std::map<reverse_edge_t, double> weightMap;
		boost::associative_property_map< std::map<reverse_edge_t, double> > property_weight_map(weightMap);

		Edge e;
		//indexing the edges
		BGL_FORALL_EDGES(e, reversed_graph, reverse_graph_t)
		{
			//one way to get the edge properties of the reversed graph
			typename property_map<reverse_graph_t, edge_properties_t>::type param = get(edge_properties, reversed_graph);
			put(property_weight_map, e, param[e].edge_weight);
			//put(property_weight_map, e, this->reaction_network_v[param[e].reaction_index].reaction_rate);
		}

		//create dijkstra visitor
		//create vectors to store the predecessors (p)
		//and corresponding edges vector (e), edges on the shortest path (e)
		my_visitor_r vis(reversed_graph, predecessors_r, predecessors_edges_r);

		//evaluate dijkstra on graph g with source s, predecessor_map p and distance_map d
		//note that predecessor_map(..).distance_map(..) is a bgl_named_params<P, T, R>, so a named parameter
		//http://stackoverflow.com/questions/31145082/bgl-dijkstra-shortest-paths-with-bundled-properties
		//http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/bundles.html
		//you can pass it to dijkstra using direct or named params. Let's do the simplest
		boost::dijkstra_shortest_paths(reversed_graph, s, weight_map(property_weight_map)
			.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
			//.predecessor_map(boost::make_iterator_property_map(predecessors_r.begin(), get(vertex_index, this->graph)))
			.visitor(vis)
		);
	}

	void superReactionNetwork::bellman_ford_reverse_graph_boost(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in, std::vector<Graph_bundled::reverse_vertex_t>& predecessors, std::vector<double>& distance)
	{
		//dijkstra on reversed graph
		//https://groups.google.com/forum/#!msg/boost-list/j8ufXCgqUdE/2SlbjLD_3pIJ
		predecessors.resize(this->num_vertices);
		distance.resize(this->num_vertices);

		reverse_graph_t reversed_graph(this->graph);

		//update_reaction_rate(time_in);
		//create a descriptor for the source node
		Vertex s = vertex(source_v, this->graph);
		Vertex t = vertex(target_v, this->graph);

		//create the property_map from edges to weights
		//doesn't work for bundled properties
		//property_map<GraphContainer, edge_weight_t>::type weightmap = get(edge_weight, g);
		std::map<reverse_edge_t, double> weightMap;
		boost::associative_property_map< std::map<reverse_edge_t, double> > property_weight_map(weightMap);

		Edge e;
		//indexing the edges
		BGL_FORALL_EDGES(e, reversed_graph, reverse_graph_t)
		{
			//one way to get the edge properties of the reversed graph
			typename property_map<reverse_graph_t, edge_properties_t>::type param = get(edge_properties, reversed_graph);
			put(property_weight_map, e, param[e].edge_weight);
			//put(property_weight_map, e, this->reaction_network_v[param[e].reaction_index].reaction_rate);
		}
		//evaluate dijkstra on graph g with source s, predecessor_map p and distance_map d
		//note that predecessor_map(..).distance_map(..) is a bgl_named_params<P, T, R>, so a named parameter
		//http://stackoverflow.com/questions/31145082/bgl-dijkstra-shortest-paths-with-bundled-properties
		//http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/bundles.html
		//you can pass it to dijkstra using direct or named params. Let's do the simplest
		//boost::dijkstra_shortest_paths(reversed_graph, s, weight_map(property_weight_map)
		//	.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
		//	.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(vertex_index, this->graph))));
		//http://www.boost.org/doc/libs/1_60_0/libs/graph/doc/bellman_ford_shortest.html
		boost::bellman_ford_shortest_paths(reversed_graph, boost::num_vertices(graph), weight_map(property_weight_map)
			.root_vertex(s)
			.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
			.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(vertex_index, this->graph))));

		//std::cout << "distances from start vertex:" << std::endl;
		//for (auto v = vertices(this->graph).first; v != vertices(this->graph).second; ++v)
		//{
		//	std::cout << "distance(" << properties(*v).vertex << ") = "
		//		<< distance[*v] << std::endl;

		//}

		////Extract the shortest path from s to t.
		//typedef std::vector<Edge> path_t;
		//path_t path;

		//for (auto u = predecessors[t]; u != t; t = u, u = predecessors[t])
		//{
		//	std::pair<Edge, bool> edge_pair = boost::edge(u, t, this->graph);
		//	path.push_back(edge_pair.first);
		//}

		//std::cout << std::endl;
		//std::cout << "Shortest Path from t to s:" << std::endl;
		//for (path_t::reverse_iterator riter = path.rbegin(); riter != path.rend(); ++riter)
		//{
		//	Vertex u_tmp = boost::source(*riter, this->graph);
		//	Vertex v_tmp = boost::target(*riter, this->graph);
		//	Edge   e_tmp = boost::edge(v_tmp, u_tmp, this->graph).first;

		//	std::cout << "  " << species_network_v[properties(u_tmp).vertex].spe_name << " -> "
		//		<< species_network_v[properties(v_tmp).vertex].spe_name << "    (weight: "
		//		<< reaction_network_v[properties(e_tmp).reaction_index].reaction_rate << ")" << std::endl;
		//}
		//t = vertex(target_v, this->graph);
		//std::cout << distance[t] << std::endl;
		//return distance[t];
	}

	void superReactionNetwork::bellman_ford_reverse_graph_boost_v2(const vertex_t source_v, const vertex_t target_v, const rsp::my_time_t time_in, std::vector<reverse_vertex_t>& predecessors_r, std::vector<reverse_edge_t>& predecessors_edges_r, std::vector<double>& distance)
	{
		//dijkstra on reversed graph
		//https://groups.google.com/forum/#!msg/boost-list/j8ufXCgqUdE/2SlbjLD_3pIJ
		//predecessors.resize(this->num_vertices);
		distance.resize(this->num_vertices);

		reverse_graph_t reversed_graph(this->graph);

		//update_reaction_rate(time_in);
		//create a descriptor for the source node
		Vertex s = vertex(source_v, this->graph);
		Vertex t = vertex(target_v, this->graph);

		//create the property_map from edges to weights
		//doesn't work for bundled properties
		//property_map<GraphContainer, edge_weight_t>::type weightmap = get(edge_weight, g);
		std::map<reverse_edge_t, double> weightMap;
		boost::associative_property_map< std::map<reverse_edge_t, double> > property_weight_map(weightMap);

		Edge e;
		//indexing the edges
		BGL_FORALL_EDGES(e, reversed_graph, reverse_graph_t)
		{
			//one way to get the edge properties of the reversed graph
			typename property_map<reverse_graph_t, edge_properties_t>::type param = get(edge_properties, reversed_graph);
			put(property_weight_map, e, param[e].edge_weight);
			//put(property_weight_map, e, this->reaction_network_v[param[e].reaction_index].reaction_rate);
		}

		//create dijkstra visitor
		//create vectors to store the predecessors (p)
		//and corresponding edges vector (e), edges on the shortest path (e)
		my_visitor_r vis(reversed_graph, predecessors_r, predecessors_edges_r);

		//evaluate dijkstra on graph g with source s, predecessor_map p and distance_map d
		//note that predecessor_map(..).distance_map(..) is a bgl_named_params<P, T, R>, so a named parameter
		//http://stackoverflow.com/questions/31145082/bgl-dijkstra-shortest-paths-with-bundled-properties
		//http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/bundles.html
		//you can pass it to dijkstra using direct or named params. Let's do the simplest
		//boost::dijkstra_shortest_paths(reversed_graph, s, weight_map(property_weight_map)
		//	.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
		//	.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(vertex_index, this->graph))));
		//http://www.boost.org/doc/libs/1_60_0/libs/graph/doc/bellman_ford_shortest.html
		boost::bellman_ford_shortest_paths(reversed_graph, boost::num_vertices(graph), weight_map(property_weight_map)
			.root_vertex(s)
			.distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, this->graph)))
			//.predecessor_map(boost::make_iterator_property_map(predecessors_r.begin(), get(vertex_index, this->graph)))
			.visitor(vis)
		);
	}

	void superReactionNetwork::print_dijkstra_algorithm_shortest_path(std::map<vertex_t, vertex_t > vertex_parent_map, const vertex_t s, const vertex_t v) const
	{
		if (s == v)
			std::cout << "shortest path:\t" << species_network_v[s].spe_name;
		else if (vertex_parent_map.find(v) == vertex_parent_map.end())
			std::cout << "no path from " << species_network_v[s].spe_name << " to " << species_network_v[v].spe_name;
		else
		{
			print_dijkstra_algorithm_shortest_path(vertex_parent_map, s, vertex_parent_map[v]);
			std::cout << species_network_v[v].spe_name;
		}
	}

	void superReactionNetwork::eppstein_alogrithm(const Vertex s, const Vertex t, double time_in, std::string atom)
	{
		//set_atom_followed_edge_weight_at_time_v0(time_in, atom, std::numeric_limits<double>::max());
		//set_atom_followed_edge_weight_at_time_v1(time_in, atom, std::numeric_limits<double>::max());
		//set_atom_followed_edge_weight_at_time_v2(time_in, atom, std::numeric_limits<double>::max());
		//set_atom_followed_edge_weight_at_time_v3(time_in, atom, std::numeric_limits<double>::max());
		set_atom_followed_edge_weight_at_time_v4(time_in, atom, std::numeric_limits<double>::max());

		const auto top_k = rnk_pt.get<std::size_t>("search_algorithm.top_k");

		//shortest path from s to t, original graph
		std::vector<Vertex> predecessors;
		std::vector<double> distance;
		std::vector<Edge> predecessors_edges;

		//bellman_ford_boost(s, t, time, predecessors, distance);
		//bellman_ford_boost_v2(s, t, time_in, predecessors, predecessors_edges, distance);
		//dijkstra_boost(s, t, time, predecessors, distance);
		dijkstra_boost_v2(s, t, time_in, predecessors, predecessors_edges, distance);

		//std::cout << "original graph:\n";
		//print_graph(this->graph, get(vertex_index, this->graph));

		//std::cout << "reversed graph:\n";
		//print_graph(make_reverse_graph(this->graph), get(vertex_index, this->graph));
		//shortest path from t to s, reversed graph
		std::vector<reverse_vertex_t> predecessors_r;
		std::vector<reverse_edge_t> predecessor_edges_r;
		std::vector<double> distance_r;
		//dijkstra_reverse_graph_boost(t, s, time, predecessors_r, distance_r);
		dijkstra_reverse_graph_boost_v2(t, s, time_in, predecessors_r, predecessor_edges_r, distance_r);
		//bellman_ford_reverse_graph_boost(t, s, time, predecessors_r, distance_r);
		//bellman_ford_reverse_graph_boost_v2(t, s, time_in, predecessors_r, predecessor_edges_r, distance_r);

		//next, generate sidetracks, build binary heap
		//Extract the shortest path from s to t.
		//vertex and edge on the shortest path
		path_vertex_t path_v; path_edge_t path_e;
		path_v.push_front(t);
		for (auto n1 = predecessors[t], n2 = predecessors[n1]; n1 != n2; n2 = n1, n1 = predecessors[n2])
		{
			path_v.push_front(n1);
			path_e.push_front(predecessors_edges[path_v[1]]);
		}

		//build sidetrack tree
		eppstein::sidetrack_tree st_tree;
		//initiate the sidetrack tree, only the first level
		initiate_sidetack_tree(st_tree, path_v, path_e, distance_r);

		//the other levels
		//root node vertex is zero
		BGL_FORALL_ADJ(0, v, st_tree.getTree(), eppstein::sidetrack_tree::TreeContainer)
		{
			add_vertex_to_sidetrack_tree_recursively(st_tree, v, st_tree.properties(v).to_vertex, distance, distance_r, 0, 0);
		}

		//build binary heap, just use priority queue
		eppstein::len_path_pq_t len_path_pq;
		//first push the shortest path, k shortest path info
		std::vector<eppstein::path_info_t > k_shortest_path_info_vec;

		eppstein::len_path_t lp_tmp; lp_tmp.cost = 0.0;
		len_path_pq.push(lp_tmp);
		k_shortest_path_with_pq(st_tree, len_path_pq, k_shortest_path_info_vec);

		std::cout << "shortest distance:\t" << distance[t] - this->delta_edge_weight << std::endl;
		//retrieve path info
		retrieve_path_info(predecessors_r, predecessor_edges_r, k_shortest_path_info_vec, s, t, top_k);
	}

	void superReactionNetwork::initiate_sidetack_tree(eppstein::sidetrack_tree& st_tree, const path_vertex_t& path_v, path_edge_t& path_e, const std::vector<double> &distance_r) const
	{
		eppstein::VertexProperties_tree parent_node;
		//add vertex, return its index in the new tree
		auto p_v = st_tree.AddVertex(parent_node);
		st_tree.properties(p_v).vertex_index_in_tree = p_v;

		//iterate over shortest path node(spn), not the last element, viz, the target node
		//for (std::size_t pi = 0; pi < path_v.size() - 1; ++pi)
		//could include the last element
		for (std::size_t pi = 0; pi < path_v.size(); ++pi)
		{
			vertex_t spn = path_v[pi];

			eppstein::VertexProperties_tree child_node;
			//all the out edge of current node
			out_edge_range_t out_e = getOutEdges(spn);

			for (auto e = out_e.first; e != out_e.second; ++e)
			{
				//edge iterator to edge descriptor
				//if this edge is on the shortest path, use arrow guard here
				if (std::find(path_e.begin(), path_e.end(), *e) != path_e.end())
					continue;
				//there exist no path from "to_vertex" to target node (infinity distance)
				if (distance_r[target(*e, this->graph)] == std::numeric_limits<double>::max())
					continue;
				//not on the shortest path, we call it sidetrack edge
				child_node.from_vertex = spn;
				child_node.to_vertex = target(*e, this->graph);
				child_node.reaction_index = properties(*e).reaction_index;
				child_node.cost = properties(*e).edge_weight + distance_r[child_node.to_vertex] - distance_r[child_node.from_vertex];
				//child_node.cost = reaction_network_v[child_node.edge].reaction_rate*properties(*e).s_coef_reactant + distance_r[child_node.to_vertex] - distance_r[child_node.from_vertex];
				child_node.path_length_before_to_vertex = pi;

				//add vertex, return its index in the new tree
				auto c_v = st_tree.AddVertex(child_node);
				st_tree.properties(c_v).vertex_index_in_tree = c_v;

				//std::cout << "add edge here\t" << p_v << "-->" << c_v << "\t" << child_node.cost << "\t" << std::endl;

				st_tree.AddEdge(p_v, c_v, eppstein::EdgeProperties_tree());
			}

		}
	}

	void superReactionNetwork::add_vertex_to_sidetrack_tree_recursively(eppstein::sidetrack_tree& st_tree, const eppstein::sidetrack_tree::Vertex curr_vertex_index_in_tree, const eppstein::vertex_index_t curr_vertex_index_in_original_graph,
		const std::vector<double>& distance, const std::vector<double>& distance_r, std::size_t level_from_current_vertex_in_original_graph, const std::size_t level_from_root) const
	{
		//use arrow guard here
		//max level from the shortest path, cannot be more than this "search_algorithm.max_level"
		if (level_from_root >= rnk_pt.get<std::size_t>("search_algorithm.max_level"))
			return;
		//not the target node itself and exist a path from current node to target node (cost is not infinity)
		if (distance_r[curr_vertex_index_in_original_graph] == 0 || distance_r[curr_vertex_index_in_original_graph] == std::numeric_limits<double>::max())
			return;

		eppstein::VertexProperties_tree child_node;
		//all the out edge of current node
		out_edge_range_t out_e = getOutEdges(curr_vertex_index_in_original_graph);

		for (auto e = out_e.first; e != out_e.second; ++e)
		{
			//there exist no path from "to_vertex" to target node (infinity distance)
			if (distance_r[target(*e, this->graph)] == std::numeric_limits<double>::max())
				continue;
			//new node information
			child_node.from_vertex = curr_vertex_index_in_original_graph;
			child_node.to_vertex = target(*e, this->graph);
			child_node.reaction_index = properties(*e).reaction_index;
			child_node.cost = properties(*e).edge_weight + distance_r[child_node.to_vertex] - distance_r[child_node.from_vertex];
			//if new node already on shortest path of preceding vertex(this must happen)
			//don't add this vertex, visit next level
			if (child_node.cost == 0.0) {
				//if go back to shortest path been visited before, level_from_current_vertex_in_original_graph+1, level_from_root remains the same
				//add_vertex_to_sidetrack_tree_recursively(st_tree, curr_vertex_index_in_tree, child_node.to_vertex, distance, distance_r, level_from_current_vertex_in_original_graph + 1, level_from_root);
				add_vertex_to_sidetrack_tree_recursively(st_tree, curr_vertex_index_in_tree, child_node.to_vertex, distance, distance_r, level_from_current_vertex_in_original_graph + 1, level_from_root + 1);
				//arrow guard
				continue;
			}

			child_node.cost += st_tree.properties(curr_vertex_index_in_tree).cost;
			child_node.path_length_before_to_vertex = st_tree.properties(curr_vertex_index_in_tree).path_length_before_to_vertex + level_from_current_vertex_in_original_graph + 1;

			//add vertex, return its index in the new tree
			auto c_v = st_tree.AddVertex(child_node);
			st_tree.properties(c_v).vertex_index_in_tree = c_v;

			std::cout << "add edge here\t" << "level_from_root:\t" << level_from_root << "\t" << "vertex in original graph:\t" << species_network_v[child_node.from_vertex].spe_name << "-->" << species_network_v[child_node.to_vertex].spe_name << "\t" << curr_vertex_index_in_tree << "-->" << c_v << "\t" << child_node.cost << "\t" << "reaction index:\t" << child_node.reaction_index << std::endl;

			st_tree.AddEdge(curr_vertex_index_in_tree, c_v, eppstein::EdgeProperties_tree());

			//recursive function
			add_vertex_to_sidetrack_tree_recursively(st_tree, c_v, child_node.to_vertex, distance, distance_r, level_from_current_vertex_in_original_graph, level_from_root + 1);
		}

	}

	void superReactionNetwork::k_shortest_path_with_pq(const eppstein::sidetrack_tree& st_tree, eppstein::len_path_pq_t& lp_pq, std::vector<eppstein::path_info_t > &K_shortest_path_info_vec) const
	{
		//if priority queue is empty, return
		if (lp_pq.empty())
			return;

		//pop out current node, push all its children to priority_queue
		auto sidetrack_e = lp_pq.top();	lp_pq.pop();
		auto curr_vertex = sidetrack_e.vertex_index_in_sidetrack_tree;

		eppstein::path_info_t path_info;
		path_info.cost = sidetrack_e.cost;
		path_info.ref_path_index = sidetrack_e.ref_path_index;
		path_info.path_length_before_to_vertex = st_tree.properties(curr_vertex).path_length_before_to_vertex;
		path_info.to_vertex = st_tree.properties(curr_vertex).to_vertex;
		path_info.reaction_index = st_tree.properties(curr_vertex).reaction_index;

		K_shortest_path_info_vec.push_back(path_info);

		//pop out current node, push all its children to priority_queue
		eppstein::len_path_t lp_tmp;
		BGL_FORALL_ADJ(curr_vertex, v, st_tree.getTree(), eppstein::sidetrack_tree::TreeContainer)
		{
			if (static_cast<vertex_index_t>(v) == curr_vertex)
				continue;
			lp_tmp.cost = st_tree.properties(v).cost;
			lp_tmp.vertex_index_in_sidetrack_tree = st_tree.properties(v).vertex_index_in_tree;
			lp_tmp.ref_path_index = K_shortest_path_info_vec.size() - 1;
			lp_pq.push(lp_tmp);
		}
		//
		k_shortest_path_with_pq(st_tree, lp_pq, K_shortest_path_info_vec);


	}

	void superReactionNetwork::retrieve_path_info(const std::vector<Vertex> &predecessors_r, const std::vector<reverse_edge_t> &predecessor_edges_r, const std::vector<eppstein::path_info_t>& K_shortest_path_info_vec,
		const Vertex s, const Vertex t, const std::size_t top_k)
	{
		//vertex
		std::vector<std::vector<vertex_t> > K_shortest_path_v_vec;
		std::vector<std::vector<rsp::index_int_t> > K_shortest_path_e_vec;
		std::vector<vertex_t> path_v; std::vector<rsp::index_int_t> path_e;
		auto n0 = s;
		while (n0 != t)
		{
			path_v.push_back(n0);
			path_e.push_back(properties(predecessor_edges_r[n0]).reaction_index);

			n0 = predecessors_r[n0];
		}
		path_v.push_back(t);
		K_shortest_path_v_vec.push_back(path_v);
		K_shortest_path_e_vec.push_back(path_e);

		//start from the second path info element
		for (std::size_t i = 1; i < K_shortest_path_info_vec.size() && i < top_k; ++i)
		{
			auto v = K_shortest_path_info_vec[i];
			//number of nodes = path length + 1
			std::vector<vertex_t> path_v_t(K_shortest_path_v_vec[v.ref_path_index].begin(),
				K_shortest_path_v_vec[v.ref_path_index].begin() + v.path_length_before_to_vertex + 1);
			std::vector<rsp::index_int_t> path_e_t(K_shortest_path_e_vec[v.ref_path_index].begin(),
				K_shortest_path_e_vec[v.ref_path_index].begin() + v.path_length_before_to_vertex);
			path_e_t.push_back(v.reaction_index);

			auto n = v.to_vertex;
			while (n != static_cast<vertex_index_t>(t))
			{
				path_v_t.push_back(n);
				path_e_t.push_back(properties(predecessor_edges_r[n]).reaction_index);

				n = predecessors_r[n];
			}
			path_v_t.push_back(t);

			K_shortest_path_v_vec.push_back(path_v_t);
			K_shortest_path_e_vec.push_back(path_e_t);

		}
		//print out
		for (std::size_t i = 0; i < K_shortest_path_v_vec.size(); ++i)
		{
			auto v1 = K_shortest_path_v_vec[i];
			auto e1 = K_shortest_path_e_vec[i];
			std::cout << K_shortest_path_info_vec[i].cost - e1.size() << "\t";
			for (std::size_t j = 0; j < v1.size(); ++j)
			{
				std::cout << species_network_v[v1[j]].spe_name;
				if (j + 1 < v1.size()) {
					std::cout << "-->";
					std::cout << e1[j] << "-->";
				}
			}
			std::cout << "\n";
		}
	}

	rsp::index_int_t superReactionNetwork::spe_random_pick_next_reaction(vertex_t curr_spe)
	{
		//probability vector
		std::vector<double> prob(this->species_network_v[curr_spe].reaction_k_index_s_coef_v.size());
		for (std::size_t i = 0; i < prob.size(); ++i) {
			prob[i] = this->species_network_v[curr_spe].reaction_k_index_s_coef_v[i].second* //s_coef_product
				this->reaction_network_v[this->species_network_v[curr_spe].reaction_k_index_s_coef_v[i].first].reaction_rate;
		}

		return this->species_network_v[curr_spe].reaction_k_index_s_coef_v[
			rand->return_index_randomly_given_probability_vector(prob)
		].first;
	}

	vertex_t superReactionNetwork::reaction_random_pick_next_spe(rsp::index_int_t reaction_index, std::string atom_followed)
	{
		//probability vector
		std::vector<double> prob(this->reaction_network_v[reaction_index].out_spe_index_weight_v_map[atom_followed].size());
		for (std::size_t i = 0; i < this->reaction_network_v[reaction_index].out_spe_index_weight_v_map[atom_followed].size(); ++i) {
			prob[i] = this->reaction_network_v[reaction_index].out_spe_index_weight_v_map[atom_followed][i].second;
		}

		return this->reaction_network_v[reaction_index].out_spe_index_weight_v_map[atom_followed][
			rand->return_index_randomly_given_probability_vector(prob)
		].first;
	}

	vertex_t reactionNetwork_sr::superReactionNetwork::spe_random_pick_next_spe(rsp::index_int_t curr_spe, std::string atom_followed)
	{
		auto spe_rxn_c1_c2_map = this->sp_all_species_group_rnk->out_species_rxns.at(curr_spe);

		std::vector<double> prob(spe_rxn_c1_c2_map.size(), 0.0);
		std::vector<std::size_t> spe_index(spe_rxn_c1_c2_map.size(), 0);

		size_t i = 0;
		for (auto s_rxn_c1_c2 : spe_rxn_c1_c2_map)
		{
			auto next_spe = s_rxn_c1_c2.first;
			spe_index[i] = next_spe;
			prob[i] = spe_spe_branching_ratio(s_rxn_c1_c2.second, -1.0, curr_spe, next_spe, atom_followed, false);

			++i;
		}

		return spe_index[rand->return_index_randomly_given_probability_vector(prob)];
	}


	void superReactionNetwork::split_atom_followed_and_pathway(std::string str_in, std::string &atom_followed, std::string &pathway) const
	{
		atom_followed.clear();
		pathway.clear();

		//find first S, record position
		auto found = str_in.find(std::string("S"));
		atom_followed = str_in.substr(0, found);
		//std::cout<< atom_followed << std::endl;
		pathway = str_in.substr(found);
		//std::cout << pathway << std::endl;

		//return true;
	}

	bool superReactionNetwork::parse_pathway_to_vector(std::string pathway_in, std::vector<rsp::index_int_t>& spe_vec, std::vector<rsp::index_int_t>& reaction_vec) const
	{
		spe_vec.resize(0);
		reaction_vec.resize(0);

		//const char* pattern="\\w\\d+";
		const char* pattern = "[-]?\\d+";
		boost::regex re(pattern);

		boost::sregex_iterator it(pathway_in.begin(), pathway_in.end(), re);
		boost::sregex_iterator end;
		std::vector<std::string> reaction_spe;
		for (boost::sregex_iterator it_t = it; it_t != end; ++it_t) {
			//std::cout<<it_t->str()<<std::endl;
			reaction_spe.push_back(it_t->str());
		}
		for (size_t i = 0; i < reaction_spe.size(); ++i) {
			if (i % 2 == 0) {
				spe_vec.push_back(boost::lexical_cast<rsp::index_int_t>(reaction_spe[i]));
			}
			else if (i % 2 == 1) {
				reaction_vec.push_back(boost::lexical_cast<rsp::index_int_t>(reaction_spe[i]));
			}
		}

		return true;
	}

	int reactionNetwork_sr::superReactionNetwork::get_number_of_elements() const
	{
		return this->element_v.size();
	}

	void reactionNetwork_sr::superReactionNetwork::set_reaction_rate(vertex_t i, double reaction_rate)
	{
		this->reaction_network_v[i].reaction_rate = reaction_rate;
	}

	void reactionNetwork_sr::superReactionNetwork::set_is_reaction_rate_nonzero_from_setting_file()
	{
		for (auto key : this->rnk_pt.get_child("pathway.non_zero_rate_reaction")) {
			this->reaction_network_v[key.second.get_value<std::size_t>()].is_reaction_rate_nonzero = true;
		}
	}

	void reactionNetwork_sr::superReactionNetwork::set_is_reaction_rate_nonzero_from_previous_iteration()
	{
		//time, half of pathway end time
		double time = 0.5 * this->rnk_pt.get<double>("time.tau");
		this->update_reaction_rate(time);
		for (std::size_t i = 0; i < this->reaction_network_v.size(); ++i) {
			if (this->reaction_network_v[i].reaction_rate > 0)
				this->reaction_network_v[i].is_reaction_rate_nonzero = true;
			else
				this->reaction_network_v[i].is_reaction_rate_nonzero = false;
		}
	}

	double superReactionNetwork::prob_spe_will_react_in_a_time_range(double init_time, double pathway_end_time, size_t curr_spe)
	{
		////set pathway end time
		//set_pathway_end_time(pathway_end_time);
		set_spe_prob_max_at_a_time(init_time, pathway_end_time, curr_spe);
		return species_network_v[curr_spe].prob_max;

	}

	double reactionNetwork_sr::superReactionNetwork::reaction_spe_branching_ratio(double reaction_time, rsp::index_int_t curr_spe, rsp::index_int_t next_reaction, rsp::index_int_t next_spe, std::string atom_followed, bool update_reaction_rate)
	{
		//update rate in the reaction network
		if (update_reaction_rate == true)
			this->update_reaction_rate(reaction_time, curr_spe);

		//probability
		double prob_total = 0.0, prob_target_reaction = 0.0;
		for (std::size_t i = 0; i < this->species_network_v[curr_spe].reaction_k_index_s_coef_v.size(); ++i) {//for
			//found next reaction
			if (this->species_network_v[curr_spe].reaction_k_index_s_coef_v[i].first == next_reaction) {
				prob_target_reaction = this->species_network_v[curr_spe].reaction_k_index_s_coef_v[i].second* //s_coef_product
					this->reaction_network_v[this->species_network_v[curr_spe].reaction_k_index_s_coef_v[i].first].reaction_rate; //reaction rate
				prob_total += prob_target_reaction;
			}
			//not found next reaction
			else {
				prob_total += this->species_network_v[curr_spe].reaction_k_index_s_coef_v[i].second* //s_coef_product
					this->reaction_network_v[this->species_network_v[curr_spe].reaction_k_index_s_coef_v[i].first].reaction_rate; //reaction rate
			}
		}//for

		double reaction_branching_ratio;
		//it prob_total ==0.0, it must because I set it to be zero artificially
		//it depends
		if (prob_total == 0.0) {
			reaction_branching_ratio = 1.0;
		}
		else {
			reaction_branching_ratio = prob_target_reaction / prob_total;
		}


		double spe_branching_ratio = 0.0;
		//next species found
		if (this->reaction_network_v[next_reaction].out_spe_index_branching_ratio_map_map[atom_followed].count(next_spe) > 0)
			spe_branching_ratio = this->reaction_network_v[next_reaction].out_spe_index_branching_ratio_map_map[atom_followed].at(next_spe);
		////next species not found
		//else {
		//	int chattering_group_id = this->species_network_v[next_spe].chattering_group_id;
		//	if (chattering_group_id != -1) {
		//		//gotta to consider the case the "next_spe" is not found, but species in the same group as "next_spe" is found
		//		for (auto n_s : this->sp_chattering_rnk->species_chattering_group_mat[chattering_group_id]) {
		//			if (this->reaction_network_v[next_reaction].out_spe_index_branching_ratio_map_map[atom_followed].count(n_s) > 0) {
		//				spe_branching_ratio = this->reaction_network_v[next_reaction].out_spe_index_branching_ratio_map_map[atom_followed].at(n_s);
		//				break;
		//			}
		//		}

		//		//treat the special case when the next species is a chattering species
		//		//multiply by the probability of being current species within the chattering group
		//		auto ss_prob_idx = this->sp_chattering_rnk->spe_idx_2_super_group_idx.at(next_spe);
		//		auto chattering_ratio = this->evaluate_chattering_group_ss_prob_at_time(reaction_time, ss_prob_idx);
		//		//check zero case
		//		if (chattering_ratio > 0.0)
		//			spe_branching_ratio *= chattering_ratio;

		//	}//chattering_group_id != -1
		//}

		return reaction_branching_ratio*spe_branching_ratio;
	}

	double reactionNetwork_sr::superReactionNetwork::spe_spe_branching_ratio(const std::vector<species_group_sr::rxn_c1_c2>& rxn_c1_c2_vec, double reaction_time, rsp::index_int_t curr_spe, rsp::index_int_t next_spe, std::string atom_followed, bool update_reaction_rate)
	{
		double ratio_tmp = 0.0;
		for (auto rxn_c1_c2 : rxn_c1_c2_vec) {
			auto reaction_index = rxn_c1_c2.r_idx;
			//here is the time is a dummy variable, since we set not to update reaction rates
			ratio_tmp += reaction_spe_branching_ratio(-1.0, curr_spe, reaction_index, next_spe, atom_followed, update_reaction_rate);
		}
		return ratio_tmp;
	}


	when_where_t superReactionNetwork::pathway_move_one_step(double time, vertex_t curr_spe, std::string & curr_pathway, std::string atom_followed)
	{
		//Monte-Carlo simulation
		//generate the random number u_1 between 0 and 1.0
		double u_1 = 0.0;
		do {
			u_1 = rand->random01();
		} while (u_1 == 1.0);
		when_where_t when_where(0.0, curr_spe);


		int chattering_group_id = this->species_network_v[curr_spe].chattering_group_id;
		//none chattering case
		if (chattering_group_id == -1) {
			time = reaction_time_from_importance_sampling_without_cutoff(time, curr_spe, u_1);
			if (time > this->tau) {
				//if curr_vertex is a dead species, should return here
				when_where.first = time;
				return when_where;
			}

			//update rate in the reaction network
			update_reaction_rate(time, curr_spe);
			rsp::index_int_t next_reaction_index = spe_random_pick_next_reaction(curr_spe);
			//random pick next spe
			vertex_t next_vertex = reaction_random_pick_next_spe(next_reaction_index, atom_followed);

			curr_pathway += "R";
			curr_pathway += boost::lexical_cast<std::string>(next_reaction_index);

			curr_pathway += "S";
			curr_pathway += boost::lexical_cast<std::string>(next_vertex);

			when_where.first = time;
			when_where.second = next_vertex;

			return when_where;
		}
		//chattering case
		//if it is chattering, and it is the first time reach chattering group, "move one step"
		//is actually move two steps, add reaction "G_{group index}"
		else {
			//calculate time from total drc of chattering species
			time = chattering_group_reaction_time_from_importance_sampling_without_cutoff(time, chattering_group_id, u_1);

			//time out of range, stop and return
			if (time > this->tau) {
				when_where.first = time;
				return when_where;
			}

			//actually move two steps, (1) from one chattering species to another chattering species
			//(2) from chattering species to the outside
			/*step */
			//choose chattering species direction randomly based on drc at this time, actually going out from that species
			std::vector<double> drc_prob(this->sp_chattering_rnk->species_chattering_group_mat[chattering_group_id].size(), 0.0);
			for (std::size_t i = 0; i < drc_prob.size(); ++i) {
				drc_prob[i] = this->evaluate_spe_drc_at_time(time,
					this->sp_chattering_rnk->species_chattering_group_mat[chattering_group_id][i]);

				//gonna take steady state concentration, or real concentration of species at this time into consideration
				//can try steady state concentration vs. real equilibrium concentration
				drc_prob[i] *= this->evaluate_spe_concentration_at_time(time,
					this->sp_chattering_rnk->species_chattering_group_mat[chattering_group_id][i]);
			}

			auto next_vertex1 = this->sp_chattering_rnk->species_chattering_group_mat[chattering_group_id][
				rand->return_index_randomly_given_probability_vector(drc_prob)
			];

			curr_pathway += "R";
			//negative reaction index represent chattering group number
			//since there is no -1 * 0, which means, to the first chattering_group_id 0, negative 0 is still 0,
			//negative sign will not show on pathway string, here we make it to be -1*(chattering_group_id+1)
			curr_pathway += boost::lexical_cast<std::string>(-1 * (chattering_group_id + rsp::INDICATOR));

			curr_pathway += "S";
			curr_pathway += boost::lexical_cast<std::string>(next_vertex1);
			/*step 1*/

			/*step 2*/
			//update rate in the reaction network
			update_reaction_rate(time, next_vertex1);
			rsp::index_int_t next_reaction_index2 = spe_random_pick_next_reaction(next_vertex1);
			//random pick next spe
			vertex_t next_vertex2 = reaction_random_pick_next_spe(next_reaction_index2, atom_followed);

			curr_pathway += "R";
			curr_pathway += boost::lexical_cast<std::string>(next_reaction_index2);

			curr_pathway += "S";
			curr_pathway += boost::lexical_cast<std::string>(next_vertex2);

			when_where.first = time;
			when_where.second = next_vertex2;
			/*step 2*/

			return when_where;
		}


	}

	std::string superReactionNetwork::pathway_sim_once(double init_time, double end_time, vertex_t init_spe, std::string atom_followed)
	{
		//set the pathway end time
		set_tau(end_time);
		std::string curr_pathway;
		when_where_t when_where(init_time, init_spe);

		//initial species
		curr_pathway += "S";
		curr_pathway += boost::lexical_cast<std::string>(init_spe);

		while (when_where.first < tau) {
			when_where = pathway_move_one_step(when_where.first, when_where.second, curr_pathway, atom_followed);
		}

		return curr_pathway;
	}

	when_where_t reactionNetwork_sr::superReactionNetwork::species_pathway_move_one_step(double time, vertex_t curr_spe, std::string & curr_pathway, std::string atom_followed)
	{
		//Monte-Carlo simulation
		//generate the random number u_1 between 0 and 1.0
		double u_1 = 0.0;
		do {
			u_1 = rand->random01();
		} while (u_1 == 1.0);
		when_where_t when_where(0.0, curr_spe);


		int chattering_group_id = this->species_network_v[curr_spe].chattering_group_id;
		//none chattering case
		if (chattering_group_id == -1) {
			time = reaction_time_from_importance_sampling_without_cutoff(time, curr_spe, u_1);
			if (time > this->tau) {
				//if curr_vertex is a dead species, should return here
				when_where.first = time;
				return when_where;
			}

			//update rate in the reaction network
			update_reaction_rate(time, curr_spe);
			//random pick next spe
			vertex_t next_vertex = spe_random_pick_next_spe(curr_spe, atom_followed);

			curr_pathway += "S";
			curr_pathway += boost::lexical_cast<std::string>(next_vertex);

			when_where.first = time;
			when_where.second = next_vertex;

			return when_where;
		}
		//chattering case
		//if it is chattering, and it is the first time reach chattering group, "move one step"
		//is actually move two steps, add reaction "G_{group index}"
		else {
			//calculate time from total drc of chattering species
			time = chattering_group_reaction_time_from_importance_sampling_without_cutoff(time, chattering_group_id, u_1);

			//time out of range, stop and return
			if (time > this->tau) {
				when_where.first = time;
				return when_where;
			}

			//actually move two steps, (1) from one chattering species to another chattering species
			//(2) from chattering species to the outside
			/*step */
			//choose chattering species direction randomly based on drc at this time, actually going out from that species
			std::vector<double> drc_prob(this->sp_chattering_rnk->species_chattering_group_mat[chattering_group_id].size(), 0.0);
			for (std::size_t i = 0; i < drc_prob.size(); ++i) {
				drc_prob[i] = this->evaluate_spe_drc_at_time(time,
					this->sp_chattering_rnk->species_chattering_group_mat[chattering_group_id][i]);

				//gonna take steady state concentration, or real concentration of species at this time into consideration
				//can try steady state concentration vs. real equilibrium concentration
				drc_prob[i] *= this->evaluate_spe_concentration_at_time(time,
					this->sp_chattering_rnk->species_chattering_group_mat[chattering_group_id][i]);
			}

			auto next_vertex1 = this->sp_chattering_rnk->species_chattering_group_mat[chattering_group_id][
				rand->return_index_randomly_given_probability_vector(drc_prob)
			];

			curr_pathway += "R";
			//negative reaction index represent chattering group number
			//since there is no -1 * 0, which means, to the first chattering_group_id 0, negative 0 is still 0,
			//negative sign will not show on pathway string, here we make it to be -1*(chattering_group_id+1)
			curr_pathway += boost::lexical_cast<std::string>(-1 * (chattering_group_id + rsp::INDICATOR));

			curr_pathway += "S";
			curr_pathway += boost::lexical_cast<std::string>(next_vertex1);
			/*step 1*/

			/*step 2*/
			//update rate in the reaction network
			update_reaction_rate(time, next_vertex1);

			//random pick next spe
			vertex_t next_vertex2 = spe_random_pick_next_spe(next_vertex1, atom_followed);

			curr_pathway += "S";
			curr_pathway += boost::lexical_cast<std::string>(next_vertex2);

			when_where.first = time;
			when_where.second = next_vertex2;
			/*step 2*/

			return when_where;
		}

	}

	std::string reactionNetwork_sr::superReactionNetwork::species_pathway_sim_once(double init_time, double end_time, vertex_t init_spe, std::string atom_followed)
	{
		//set the pathway end time
		set_tau(end_time);
		std::string curr_pathway;
		when_where_t when_where(init_time, init_spe);

		//initial species
		curr_pathway += "S";
		curr_pathway += boost::lexical_cast<std::string>(init_spe);

		while (when_where.first < tau) {
			when_where = species_pathway_move_one_step(when_where.first, when_where.second, curr_pathway, atom_followed);
		}

		return curr_pathway;
	}

	double superReactionNetwork::pathway_prob_sim_move_one_step(double when_time, vertex_t curr_spe, rsp::index_int_t next_reaction, vertex_t next_spe, double & pathway_prob, std::string atom_followed)
	{
		if (when_time >= (tau - INFINITESIMAL_DT)) {
			return when_time;

		}

		this->set_spe_prob_max_at_a_time(when_time, tau, curr_spe);

		double u_1;
		if (species_network_v[curr_spe].prob_max > 0.0) {
			u_1 = rand->random_min_max(0, species_network_v[curr_spe].prob_max);
		}
		else {
			u_1 = 0.0;
		}

		when_time = reaction_time_from_importance_sampling(when_time, curr_spe, u_1);

		pathway_prob *= reaction_spe_branching_ratio(when_time, curr_spe, next_reaction, next_spe, atom_followed);

		return when_time;
	}

	double superReactionNetwork::pathway_prob_input_pathway_sim_once(double const init_time, const double pathway_end_time, const std::vector<rsp::index_int_t> &spe_vec, const std::vector<rsp::index_int_t> &reaction_vec, std::string atom_followed)
	{
		//set pathway end time
		set_tau(pathway_end_time);

		//basically, we assume there must be a reaction at the beginning, so should multiply be the 1-P_min(tau=0|t;S^{0})
		double pathway_prob = 1.0;
		double when_time = init_time;

		//start from the first reaction
		for (size_t i = 0; i < reaction_vec.size();)
		{
			//none-chattering reaction
			if (reaction_vec[i] >= 0) {
				pathway_prob *= prob_spe_will_react_in_a_time_range(when_time, pathway_end_time, spe_vec[i]);
				when_time = pathway_prob_sim_move_one_step(when_time, spe_vec[i], reaction_vec[i], spe_vec[i + 1], pathway_prob, atom_followed);
				//move one step
				++i;
			}
			//chattering reaction, chattering case
			else {
				int chattering_group_id = this->species_network_v[spe_vec[i]].chattering_group_id;

				//add time delay first, regenerate random number, inverse to get exact time, get steady state time first
				//then calculate steady state ratios
				double chattering_group_prob = prob_chattering_group_will_react_in_a_time_range(when_time, pathway_end_time, chattering_group_id);
				pathway_prob *= chattering_group_prob;

				//avoid problems around boundary
				if (when_time < (tau - INFINITESIMAL_DT)) {
					double u_1 = 1.0;
					if (chattering_group_prob > 0.0) {
						u_1 = rand->random_min_max(0, chattering_group_prob);
					}
					else
						u_1 = 0.0;

					when_time = chattering_group_reaction_time_from_importance_sampling_without_cutoff(when_time, chattering_group_id, u_1);

					/*step 1*/
					//based on drc at this time, calculate probability going out by that direction
					std::vector<double> drc_prob(this->sp_chattering_rnk->species_chattering_group_mat[chattering_group_id].size(), 0.0);
					for (std::size_t i = 0; i < drc_prob.size(); ++i) {
						drc_prob[i] = this->evaluate_spe_drc_at_time(when_time,
							this->sp_chattering_rnk->species_chattering_group_mat[chattering_group_id][i]);

						//gonna take steady state concentration, or real concentration of species at this time into consideration
						//can try steady state concentration vs. real equilibrium concentration
						drc_prob[i] *= this->evaluate_spe_concentration_at_time(when_time,
							this->sp_chattering_rnk->species_chattering_group_mat[chattering_group_id][i]);
					}
					double drc_prob_sum = std::accumulate(drc_prob.begin(), drc_prob.end(), 0.0);
					//make sure there is at least one direction out, there is no, dead end, return 0.0 probability
					if (drc_prob_sum <= 0.0)
						return 0.0;

					//notice out species is spe_vec[i + 1], next_species1
					pathway_prob *= drc_prob[this->sp_chattering_rnk->spe_idx_2_chattering_group_id_idx[spe_vec[i + 1]].second] / drc_prob_sum;
					/*step 1*/

					/*step 2*/
					pathway_prob *= reaction_spe_branching_ratio(when_time, spe_vec[i + 1], reaction_vec[i + 1], spe_vec[i + 2], atom_followed);
					/*step 2*/

				}//boundary time problem

				//move two steps actually
				i += 2;
			}//if chattering case

		}

		//got to multiply by P_min or says (1-P_max)
		set_spe_prob_max_at_a_time(when_time, pathway_end_time, spe_vec.back());

		pathway_prob *= (1 - species_network_v[spe_vec.back()].prob_max);

		return pathway_prob;
	}

	double reactionNetwork_sr::superReactionNetwork::species_pathway_prob_sim_move_one_step(double when_time, vertex_t curr_spe, vertex_t next_spe, double & pathway_prob, std::string atom_followed)
	{
		if (when_time >= (tau - INFINITESIMAL_DT)) {
			return when_time;

		}

		this->set_spe_prob_max_at_a_time(when_time, tau, curr_spe);

		double u_1;
		if (species_network_v[curr_spe].prob_max > 0.0) {
			u_1 = rand->random_min_max(0, species_network_v[curr_spe].prob_max);
		}
		else {
			u_1 = 0.0;
		}

		when_time = reaction_time_from_importance_sampling(when_time, curr_spe, u_1);

		//pathway_prob *= reaction_spe_branching_ratio(when_time, curr_spe, next_reaction, next_spe, atom_followed);

		return when_time;
	}


	void reactionNetwork_sr::superReactionNetwork::initiate_M_matrix(std::string atom_followed)
	{
		//resize and initialization
		this->atom_M_matrix[atom_followed].resize(this->species_network_v.size());
		for (std::size_t i = 0; i < this->atom_M_matrix[atom_followed].size(); ++i)
			this->atom_M_matrix[atom_followed][i].resize(this->species_network_v.size());
		for (std::size_t i = 0; i < this->atom_M_matrix[atom_followed].size(); ++i)
			for (std::size_t j = 0; j < this->atom_M_matrix[atom_followed][i].size(); ++j)
				this->atom_M_matrix[atom_followed][i][j] = 0;

		//actually build M-matrix
		for (std::size_t i = 0; i < this->atom_M_matrix[atom_followed].size(); ++i) {
			for (std::size_t j = 0; j < this->species_network_v[i].reaction_k_index_s_coef_v.size(); ++j) {
				for (std::size_t k = 0; k < this->reaction_network_v[this->species_network_v[i].reaction_k_index_s_coef_v[j].first].out_spe_index_weight_v_map[atom_followed].size(); ++k) {
					//both contains atom followed
					if (this->species_network_v[i].spe_component[atom_followed] != 0 && this->species_network_v[this->reaction_network_v[this->species_network_v[i].reaction_k_index_s_coef_v[j].first].out_spe_index_weight_v_map[atom_followed][k].first].spe_component[atom_followed] != 0) {
						this->atom_M_matrix[atom_followed][i][this->reaction_network_v[this->species_network_v[i].reaction_k_index_s_coef_v[j].first].out_spe_index_weight_v_map[atom_followed][k].first] += 1;
					}
				}
			}
		}

	}

	void reactionNetwork_sr::superReactionNetwork::initiate_M_matrix()
	{
		for (auto x : this->element_v)
			this->initiate_M_matrix(x.ele_name);

		this->initiate_M_matrix(rnk_pt.get<std::string>("pathway.super_atom"));
	}

	void reactionNetwork_sr::superReactionNetwork::print_M_matrix(std::string atom_followed)
	{
		for (std::size_t i = 0; i < this->atom_M_matrix[atom_followed].size(); ++i) {
			std::cout << i << ",\t" << this->species_network_v[i].spe_name << ",\t";
			for (std::size_t j = 0; j < this->atom_M_matrix[atom_followed][i].size(); ++j) {
				std::cout << this->atom_M_matrix[atom_followed][i][j] << "\t";
			}
			std::cout << "\n";
		}

	}

	matrix_sr::size_t_matrix_t reactionNetwork_sr::superReactionNetwork::return_M_matrix(std::string atom_followed)
	{
		return this->atom_M_matrix[atom_followed];
	}

	void reactionNetwork_sr::superReactionNetwork::initiate_R_matrix_v1(std::string atom_followed)
	{
		//resize and initialization
		this->atom_R_matrix[atom_followed].resize(this->species_network_v.size());
		for (std::size_t i = 0; i < this->atom_R_matrix[atom_followed].size(); ++i)
			this->atom_R_matrix[atom_followed][i].resize(this->species_network_v.size());
		for (std::size_t i = 0; i < this->atom_R_matrix[atom_followed].size(); ++i)
			for (std::size_t j = 0; j < this->atom_R_matrix[atom_followed][i].size(); ++j)
				this->atom_R_matrix[atom_followed][i][j] = {};

		//actually build M-matrix
		for (std::size_t i = 0; i < this->atom_R_matrix[atom_followed].size(); ++i) {
			for (std::size_t j = 0; j < this->species_network_v[i].reaction_k_index_s_coef_v.size(); ++j) {
				for (std::size_t k = 0; k < this->reaction_network_v[this->species_network_v[i].reaction_k_index_s_coef_v[j].first].out_spe_index_weight_v_map[atom_followed].size(); ++k) {
					//both contains atom followed
					if (this->species_network_v[i].spe_component[atom_followed] != 0 && this->species_network_v[this->reaction_network_v[this->species_network_v[i].reaction_k_index_s_coef_v[j].first].out_spe_index_weight_v_map[atom_followed][k].first].spe_component[atom_followed] != 0) {
						this->atom_R_matrix[atom_followed][i][this->reaction_network_v[this->species_network_v[i].reaction_k_index_s_coef_v[j].first].out_spe_index_weight_v_map[atom_followed][k].first].push_back({ this->species_network_v[i].reaction_k_index_s_coef_v[j].first });
					}
				}
			}
		}

	}

	void reactionNetwork_sr::superReactionNetwork::initiate_R_matrix_v2(std::string atom_followed)
	{
		//resize and initialization
		this->atom_R_matrix[atom_followed].resize(this->species_network_v.size());
		for (std::size_t i = 0; i < this->atom_R_matrix[atom_followed].size(); ++i)
			this->atom_R_matrix[atom_followed][i].resize(this->species_network_v.size());
		for (std::size_t i = 0; i < this->atom_R_matrix[atom_followed].size(); ++i)
			for (std::size_t j = 0; j < this->atom_R_matrix[atom_followed][i].size(); ++j)
				this->atom_R_matrix[atom_followed][i][j] = {};

		//iterate over all edges
		edge_iter iter_beg, iter_end;
		boost::tie(iter_beg, iter_end) = getEdges();
		//source and target index
		std::size_t s_i, t_i;
		for (; iter_beg != iter_end; ++iter_beg) {
			s_i = boost::source(*iter_beg, this->graph);
			t_i = boost::target(*iter_beg, this->graph);

			if (species_network_v[s_i].spe_component[atom_followed] != 0 && species_network_v[t_i].spe_component[atom_followed] != 0)
				this->atom_R_matrix[atom_followed][s_i][t_i].push_back({ properties(*iter_beg).edge_index });
		}

	}

	void reactionNetwork_sr::superReactionNetwork::initiate_R_matrix()
	{
		for (auto x : this->element_v)
			//this->initiate_R_matrix_v1(x.ele_name);
			this->initiate_R_matrix_v2(x.ele_name);

		this->initiate_R_matrix_v2(this->rnk_pt.get<std::string>("pathway.super_atom"));
	}

	matrix_sr::path_R_matrix_t reactionNetwork_sr::superReactionNetwork::return_R_matrix(std::string atom_followed)
	{
		return this->atom_R_matrix[atom_followed];
	}

	void reactionNetwork_sr::superReactionNetwork::print_R_matrix(std::string atom_followed)
	{
		for (std::size_t i = 0; i < this->atom_R_matrix[atom_followed].size(); ++i) {
			std::cout << i << ",\t" << this->species_network_v[i].spe_name << ",\t";
			for (std::size_t j = 0; j < this->atom_R_matrix[atom_followed][i].size(); ++j) {
				//std::cout << "(" << i << "," << j << ")\t";
				std::cout << "(" << this->species_network_v[i].spe_name << "," << this->species_network_v[j].spe_name << ")\t";
				if (this->atom_R_matrix[atom_followed][i][j].size() == 0) {
					continue;
				}
				for (std::size_t k = 0; k < this->atom_R_matrix[atom_followed][i][j].size(); ++k) {
					for (std::size_t l = 0; l < this->atom_R_matrix[atom_followed][i][j][k].size(); ++l)
						std::cout << this->atom_R_matrix[atom_followed][i][j][k][l] << "\t";
				}
			}
			std::cout << std::endl;
		}
	}

	std::size_t reactionNetwork_sr::superReactionNetwork::get_M_matrix_element(std::string atom_followed, std::size_t i, std::size_t j)
	{
		return this->atom_M_matrix[atom_followed][i][j];
	}

	matrix_sr::path_R_matrix_element_t reactionNetwork_sr::superReactionNetwork::get_R_matrix_element(std::string atom_followed, std::size_t i, std::size_t j)
	{
		return this->atom_R_matrix[atom_followed][i][j];
	}

	std::string reactionNetwork_sr::superReactionNetwork::R_matrix_path_representation_to_string(matrix_sr::path_t p)
	{
		std::string path_t;
		if (p.size() > 0) {
			edge_iter iter_e = edge_index_to_edge_iterator[p[0]];
			path_t += std::string("S") + boost::lexical_cast<std::string>(boost::source(*iter_e, this->graph));

			for (std::size_t i = 0; i < p.size(); ++i) {
				iter_e = edge_index_to_edge_iterator[p[i]];
				path_t += std::string("R") + boost::lexical_cast<std::string>(properties(*iter_e).reaction_index);
				path_t += std::string("S") + boost::lexical_cast<std::string>(boost::target(*iter_e, this->graph));
			}
		}
		return path_t;
	}

	bool reactionNetwork_sr::superReactionNetwork::contains_zero_reaction_rate_reactions(matrix_sr::path_t p)
	{
		//arrow guard
		if (p.size() == 0)
			return false;
		for (auto x : p) {
			edge_iter iter_e = edge_index_to_edge_iterator[x];
			if (reaction_network_v[properties(*iter_e).reaction_index].is_reaction_rate_nonzero == false)
				return false;
		}
		return true;
	}

	std::vector<std::string> reactionNetwork_sr::superReactionNetwork::get_path_string_element_i_j(const matrix_sr::path_R_matrix_t &pRm, std::size_t i, std::size_t j)
	{
		if (pRm.size() == 0)
			return std::vector<std::string>(1, std::string("S") + boost::lexical_cast<std::string>(i));

		matrix_sr::path_R_matrix_element_t p = pRm[i][j];
		std::vector<std::string> vs;

		for (auto x : p) {
			std::string ps = R_matrix_path_representation_to_string(x);
			if (!ps.empty())
				vs.push_back(ps);
		}
		return vs;
	}

	std::vector<std::string> reactionNetwork_sr::superReactionNetwork::get_path_string_update_matrix_element_i_j_topN(matrix_sr::path_R_matrix_t &pRm, const std::size_t i, const std::size_t j,
		const std::string atom_followed, const std::size_t topN, const double start_time, const double end_time)
	{
		if (pRm.size() == 0)
			return std::vector<std::string>(1, std::string("S") + boost::lexical_cast<std::string>(i));

		matrix_sr::path_R_matrix_element_t p = pRm[i][j];
		matrix_sr::path_R_matrix_element_t p_new;
		std::multimap<double, std::pair<std::string, std::size_t>, std::greater<double> > prob_path_map;

		//if there is less than topN path, do nothing
		//if there are more than topN path, delete the unimportant ones
		for (std::size_t k = 0; k < p.size(); ++k) {
			std::string ps = R_matrix_path_representation_to_string(p[k]);
			double prob = calculate_path_weight_based_on_path_probability(ps, atom_followed, start_time, end_time);

			if (prob_path_map.size() < topN) {
				prob_path_map.insert(std::make_pair(prob, std::make_pair(ps, k)));
			}
			else
			{
				if (prob <= prob_path_map.crbegin()->first)
					continue;
				else {
					prob_path_map.erase(std::prev(prob_path_map.end()));
					prob_path_map.insert(std::make_pair(prob, std::make_pair(ps, k)));
				}
			} //if <topN

		}

		for (auto x : prob_path_map) {
			p_new.push_back(p[x.second.second]);
		}

		//update matrix
		if (p_new.size() == 0)
			p_new = {};
		pRm[i][j].clear();
		pRm[i][j] = p_new;

		std::vector<std::string> vs;
		for (auto x : p_new) {
			std::string ps = R_matrix_path_representation_to_string(x);
			if (!ps.empty())
				vs.push_back(ps);
		}

		return vs;
	}

	void reactionNetwork_sr::superReactionNetwork::path_string_vector_s2f(std::vector<std::string> vs, std::string filename)
	{
		std::ofstream fout(filename.c_str());
		for (auto x : vs)
			fout << x << std::endl;

		fout.close(); fout.clear();
	}

	void reactionNetwork_sr::superReactionNetwork::heuristic_path_string_vector_s2f(std::string atom_followed, std::size_t n, std::string filename)
	{
		std::unordered_set<std::string> us;
		for (std::size_t k = 0; k <= n; ++k) {
			auto pRmn = matrix_sr::matrix_power(this->atom_R_matrix[atom_followed], k);
			for (auto key : this->rnk_pt.get_child("chem_init.species_index_concentration")) {
				std::size_t si = boost::lexical_cast<std::size_t>(key.first);
				//doesn't contain atom_followed
				if (species_network_v[si].spe_component.at(atom_followed) == 0)
					continue;
				if (k == 0) {
					us.insert(std::string("S") + boost::lexical_cast<std::string>(si));
					//std::cout << std::string("S") + boost::lexical_cast<std::string>(si) << std::endl;
					continue;
				}
				for (std::size_t sj = 0; sj < this->species_network_v.size(); ++sj) {
					auto vs = this->get_path_string_update_matrix_element_i_j_topN(pRmn, si, sj);
					for (auto s : vs)
						us.insert(s);
					//std::cout << s << std::endl;
				}
			}

		}

		//save to file
		std::ofstream fout(filename.c_str());
		for (auto s : us)
			fout << s << std::endl;
		fout.close(); fout.clear();
	}

	std::set<std::string> reactionNetwork_sr::superReactionNetwork::heuristic_path_string_vector_s2m(std::string atom_followed, std::size_t n)
	{
		std::set<std::string> us;
		for (std::size_t k = 0; k <= n; ++k) {
			auto pRmn = matrix_sr::matrix_power(this->atom_R_matrix[atom_followed], k);
			for (auto key : this->rnk_pt.get_child("chem_init.species_index_concentration")) {
				std::size_t si = boost::lexical_cast<std::size_t>(key.first);
				//doesn't contain atom_followed
				if (species_network_v[si].spe_component.at(atom_followed) == 0)
					continue;
				if (k == 0) {
					us.insert(std::string("S") + boost::lexical_cast<std::string>(si));
					//std::cout << std::string("S") + boost::lexical_cast<std::string>(si) << std::endl;
					continue;
				}
				for (std::size_t sj = 0; sj < this->species_network_v.size(); ++sj) {
					auto vs = this->get_path_string_update_matrix_element_i_j_topN(pRmn, si, sj);
					for (auto s : vs)
						us.insert(s);
					//std::cout << s << std::endl;
				}
			}

		}
		return us;
	}

	std::set<std::string> reactionNetwork_sr::superReactionNetwork::heuristic_path_string_vector_sorted_based_on_path_length(std::string atom_followed, std::size_t n, std::size_t topN)
	{
		std::vector<std::multimap<double, std::string> > p_map_v(this->species_network_v.size());

		matrix_sr::path_R_matrix_t pRmn;
		for (std::size_t k = 0; k <= n; ++k) {
			//auto pRmn = matrix_sr::matrix_power(this->atom_R_matrix[atom_followed], k);
			if (k == 0)
				pRmn = matrix_sr::path_R_matrix_t();
			else if (k == 1)
				pRmn = this->atom_R_matrix[atom_followed];
			else
				pRmn = matrix_sr::matrix_multiplication(pRmn, this->atom_R_matrix[atom_followed]);

			for (auto key : this->rnk_pt.get_child("chem_init.species_index_concentration")) {
				std::size_t si = boost::lexical_cast<std::size_t>(key.first);
				//doesn't contain atom_followed
				if (species_network_v[si].spe_component.at(atom_followed) == 0)
					continue;
				if (k == 0) {
					std::string path_name = std::string("S") + boost::lexical_cast<std::string>(si);
					p_map_v[si].insert(std::make_pair(calculate_path_weight_path_length(path_name), path_name));
					continue;
				}
				for (std::size_t sj = 0; sj < this->species_network_v.size(); ++sj) {
					auto vs = this->get_path_string_update_matrix_element_i_j_topN(pRmn, si, sj);
					for (auto s : vs) {
						double prob = calculate_path_weight_path_length(s);
						if (p_map_v[sj].size() < topN) {
							p_map_v[sj].insert(std::make_pair(prob, s));
						}
						else
						{
							if (prob >= p_map_v[sj].crbegin()->first)
								continue;
							else {
								p_map_v[sj].erase(std::prev(p_map_v[sj].end()));
								p_map_v[sj].insert(std::make_pair(prob, s));
							}
						} //if <topN


					} //auto s
				} //sj
			}

		}
		std::set<std::string> us;
		for (auto pmv : p_map_v)
			for (auto ps : pmv)
				//add atom followed to the beging of path
				us.insert(atom_followed + ps.second);

		return us;
	}

	std::set<std::string> reactionNetwork_sr::superReactionNetwork::heuristic_path_string_vector_sorted_based_on_path_prob(std::string atom_followed, std::size_t n, std::size_t topN, double end_time_ratio)
	{
		std::vector<std::multimap<double, std::string, std::greater<double> > > prob_path_map_v(this->species_network_v.size());
		matrix_sr::path_R_matrix_t pRmn;
		for (std::size_t k = 0; k <= n; ++k) {
			//auto pRmn = matrix_sr::matrix_power(this->atom_R_matrix[atom_followed], k);
			if (k == 0)
				pRmn = matrix_sr::path_R_matrix_t();
			else if (k == 1)
				pRmn = this->atom_R_matrix[atom_followed];
			else
				pRmn = matrix_sr::matrix_multiplication(pRmn, this->atom_R_matrix[atom_followed]);

			auto species_with_initial_concentration = return_species_index_with_initial_concentration();
			auto species_without_initial_concentration = return_species_index_without_initial_concentration();

			//species with initial concentration
			for (auto si : species_with_initial_concentration) {
				//doesn't contain atom_followed
				if (species_network_v[si].spe_component.at(atom_followed) == 0)
					continue;
				if (k == 0) {
					std::string path_name = std::string("S") + boost::lexical_cast<std::string>(si);
					prob_path_map_v[si].insert(std::make_pair(calculate_path_weight_based_on_path_probability(path_name, atom_followed, 0.0, end_time_ratio*this->rnk_pt.get<double>("time.tau")), path_name));
					continue;
				}
				//in the mean time, we should change the matrix element so that it doesn't contain too many elements
				//become too big-->lots of memory
				for (std::size_t sj = 0; sj < this->species_network_v.size(); ++sj) {
					//be a little cautious, a little open, 10*topN
					auto vs = this->get_path_string_update_matrix_element_i_j_topN(pRmn, si, sj, atom_followed, 10 * topN, 0.0, end_time_ratio*this->rnk_pt.get<double>("time.tau"));
					for (auto s : vs) {
						double prob = calculate_path_weight_based_on_path_probability(s, atom_followed, 0.0, end_time_ratio*this->rnk_pt.get<double>("time.tau"));
						if (prob_path_map_v[sj].size() < topN) {
							prob_path_map_v[sj].insert(std::make_pair(prob, s));
						}
						else
						{
							if (prob <= prob_path_map_v[sj].crbegin()->first)
								continue;
							else {
								prob_path_map_v[sj].erase(std::prev(prob_path_map_v[sj].end()));
								prob_path_map_v[sj].insert(std::make_pair(prob, s));
							}
						} //if <topN


					} //auto s
				} //sj
			}

			//species without initial concentration, still need to update matrix elements
			for (auto si : species_without_initial_concentration) {
				for (std::size_t sj = 0; sj < this->species_network_v.size(); ++sj)
					//be a little cautious, a little open, 10*topN
					this->get_path_string_update_matrix_element_i_j_topN(pRmn, si, sj, atom_followed, 10 * topN, 0.0, end_time_ratio*this->rnk_pt.get<double>("time.tau"));
			}

		}
		std::set<std::string> us;
		for (auto pmv : prob_path_map_v)
			for (auto ps : pmv)
				//add atom followed to the beging of path
				us.insert(atom_followed + ps.second);

		return us;
	}

	double reactionNetwork_sr::superReactionNetwork::calculate_path_weight_path_length(std::string path)
	{
		return (double)(std::count(path.begin(), path.end(), 'S'));
	}

	double reactionNetwork_sr::superReactionNetwork::calculate_path_weight_based_on_path_probability(std::string path, std::string atom_followed, double start_time, double end_time)
	{
		std::vector<rsp::index_int_t> spe_vec; std::vector<rsp::index_int_t> reaction_vec;
		double prob = 0.0;
		this->parse_pathway_to_vector(path, spe_vec, reaction_vec);
		prob = pathway_prob_input_pathway_sim_once(start_time, end_time, spe_vec, reaction_vec, atom_followed);

		//take the initial concentration of initial species into account
		if (this->species_network_v[spe_vec[0]].spe_conc != 0)
			prob *= this->species_network_v[spe_vec[0]].spe_conc;

		prob *= this->species_network_v[spe_vec[0]].spe_component[atom_followed];
		prob /= this->species_network_v[spe_vec.back()].spe_component[atom_followed];

		return prob;
	}

	std::size_t reactionNetwork_sr::superReactionNetwork::heuristic_path_string_vector_by_stage_number_path_length_all_elements(const std::size_t stage_n, std::string filename, std::size_t topN)
	{
		assert(stage_n >= 0);
		//fetch path length first
		std::vector<std::size_t> path_n_v;
		for (auto key : this->rnk_pt.get_child("pathway.max_path_length"))
			path_n_v.push_back(key.second.get_value<std::size_t>());

		std::size_t path_n;
		//if iteration n is less than path_n_v lenght, fetch by index, otherwise take the last element
		if (stage_n < path_n_v.size())
			path_n = path_n_v[stage_n];
		else
			path_n = path_n_v.back();

		if (stage_n == 0)
			this->set_is_reaction_rate_nonzero_from_setting_file();
		else
			this->set_is_reaction_rate_nonzero_from_previous_iteration();

		std::vector<set<std::string> > all_path_2;

		for (auto x : this->element_v) {
			all_path_2.push_back(this->heuristic_path_string_vector_sorted_based_on_path_length(x.ele_name, path_n, topN));
		}

		set<std::string> all_path;
		for (auto us : all_path_2)
			for (auto s : us)
				all_path.insert(s);

		//save to file
		std::ofstream fout(filename.c_str());
		for (auto x : all_path)
			fout << x << std::endl;
		fout.close(); fout.clear();

		return this->element_v.size();
	}

	std::size_t reactionNetwork_sr::superReactionNetwork::heuristic_path_string_vector_by_stage_number_path_prob_all_elements(const std::size_t stage_n, std::string filename, std::size_t topN, double end_time_ratio)
	{
		assert(stage_n >= 0);
		//fetch path length first
		std::vector<std::size_t> path_n_v;
		for (auto key : this->rnk_pt.get_child("pathway.max_path_length"))
			path_n_v.push_back(key.second.get_value<std::size_t>());

		std::size_t path_n;
		//if iteration n is less than path_n_v lenght, fetch by index, otherwise take the last element
		if (stage_n < path_n_v.size())
			path_n = path_n_v[stage_n];
		else
			path_n = path_n_v.back();

		if (stage_n == 0)
			this->set_is_reaction_rate_nonzero_from_setting_file();
		else
			this->set_is_reaction_rate_nonzero_from_previous_iteration();

		std::vector<set<std::string> > all_path_2;

		for (auto x : this->element_v) {
			all_path_2.push_back(this->heuristic_path_string_vector_sorted_based_on_path_prob(x.ele_name, path_n, topN, end_time_ratio));
		}

		set<std::string> all_path;
		for (auto us : all_path_2)
			for (auto s : us)
				all_path.insert(s);

		//save to file
		std::ofstream fout(filename.c_str());
		for (auto x : all_path)
			fout << x << std::endl;
		fout.close(); fout.clear();

		return this->element_v.size();
	}

	std::size_t reactionNetwork_sr::superReactionNetwork::heuristic_path_string_vector_by_stage_number_path_prob_all_elements_s2m(const std::size_t stage_n, std::vector<std::string> &path_all_v, std::size_t topN, double end_time_ratio)
	{
		path_all_v.resize(0);

		assert(stage_n >= 0);
		//fetch path length first
		std::vector<std::size_t> path_n_v;
		for (auto key : this->rnk_pt.get_child("pathway.max_path_length"))
			path_n_v.push_back(key.second.get_value<std::size_t>());

		std::size_t path_n;
		//if iteration n is less than path_n_v lenght, fetch by index, otherwise take the last element
		if (stage_n < path_n_v.size())
			path_n = path_n_v[stage_n];
		else
			path_n = path_n_v.back();

		if (stage_n == 0)
			this->set_is_reaction_rate_nonzero_from_setting_file();
		else
			this->set_is_reaction_rate_nonzero_from_previous_iteration();

		std::vector<set<std::string> > all_path_2;

		for (auto x : this->element_v) {
			all_path_2.push_back(this->heuristic_path_string_vector_sorted_based_on_path_prob(x.ele_name, path_n, topN, end_time_ratio));
		}

		set<std::string> all_path;
		for (auto us : all_path_2)
			for (auto s : us)
				all_path.insert(s);
		for (auto p : all_path)
			path_all_v.push_back(p);

		return this->element_v.size();

	}


	std::size_t reactionNetwork_sr::superReactionNetwork::heuristic_path_string_vector_by_stage_number_path_prob_super_element(const std::size_t stage_n, std::string filename, std::size_t topN, double end_time_ratio)
	{
		assert(stage_n >= 0);
		//fetch path length first
		std::vector<std::size_t> path_n_v;
		for (auto key : this->rnk_pt.get_child("pathway.max_path_length"))
			path_n_v.push_back(key.second.get_value<std::size_t>());

		std::size_t path_n;
		//if iteration n is less than path_n_v lenght, fetch by index, otherwise take the last element
		if (stage_n < path_n_v.size())
			path_n = path_n_v[stage_n];
		else
			path_n = path_n_v.back();

		if (stage_n == 0)
			this->set_is_reaction_rate_nonzero_from_setting_file();
		else
			this->set_is_reaction_rate_nonzero_from_previous_iteration();

		auto us = this->heuristic_path_string_vector_sorted_based_on_path_prob(this->rnk_pt.get<std::string>("pathway.super_atom"), path_n, topN, end_time_ratio);

		//save to file
		std::ofstream fout(filename.c_str());
		for (auto x : us)
			fout << x << std::endl;
		fout.close(); fout.clear();

		return 1;
	}

	std::set<std::size_t> reactionNetwork_sr::superReactionNetwork::return_species_index_with_initial_concentration() const
	{
		std::set<std::size_t> species_with_initial_concentration;
		for (auto key : this->rnk_pt.get_child("chem_init.species_index_concentration")) {
			species_with_initial_concentration.insert(boost::lexical_cast<std::size_t>(key.first));
		}
		return species_with_initial_concentration;
	}

	std::set<std::size_t> reactionNetwork_sr::superReactionNetwork::return_species_index_without_initial_concentration() const
	{
		auto species_with_initial_concentration = return_species_index_with_initial_concentration();
		std::set<std::size_t> species_without_initial_concentration;

		for (std::size_t i = 0; i < this->species_network_v.size(); ++i) {
			//not in species_with_initial_concentration
			if (species_with_initial_concentration.count(i) == 0)
				species_without_initial_concentration.insert(i);
		}
		return species_without_initial_concentration;
	}

	std::set<std::pair<std::size_t, double> > reactionNetwork_sr::superReactionNetwork::return_species_index_and_initial_concentration() const
	{
		std::set<std::pair<std::size_t, double> > species_concentration;

		for (auto key : this->rnk_pt.get_child("chem_init.species_index_concentration")) {
			species_concentration.emplace(boost::lexical_cast<std::size_t>(key.first), key.second.get_value<double>());
		}
		return species_concentration;
	}


	void reactionNetwork_sr::superReactionNetwork::generate_path_by_running_monte_carlo_trajectory_s2m(std::vector<statistics > &statistics_v, std::size_t Ntrajectory, std::string atom_followed, double end_time_ratio)
	{
		//std::vector<statistics > statistics_v(this->species_network_v.size());

		auto species_with_initial_concentration = return_species_index_with_initial_concentration();

		//species with initial concentration, initial concentration is not zero
		for (auto si : species_with_initial_concentration) {
			//doesn't contain atom_followed
			if (species_network_v[si].spe_component.at(atom_followed) == 0)
				continue;

			//contain atom_followed
			std::string str_tmp;
			for (std::size_t ti = 0; ti < Ntrajectory; ++ti) {
				str_tmp = this->pathway_sim_once(0.0, end_time_ratio*this->rnk_pt.get<double>("time.tau"), si, atom_followed);
				//put the atom followed in front
				statistics_v[si].insert_pathway_stat(atom_followed + str_tmp);
			} //for
		} //for

	} //generate_path_by_running_monte_carlo_trajectory_s2m

	std::size_t reactionNetwork_sr::superReactionNetwork::generate_path_by_running_monte_carlo_trajectory_all_elements_s2m(std::vector<statistics>& statistics_v, std::size_t Ntrajectory, double end_time_ratio)
	{
		for (auto x : this->element_v) {
			generate_path_by_running_monte_carlo_trajectory_s2m(statistics_v, Ntrajectory, x.ele_name, end_time_ratio);
		}
		return this->element_v.size();
	}

	std::vector<rsp::element_info> reactionNetwork_sr::superReactionNetwork::return_element_vecotr() const
	{
		return this->element_v;
	}



}/*namespace reactionNetwork_sr*/

#endif
