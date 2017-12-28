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
		read_chem_out_spe_for_network_info(edgeVector, edgePro, vertex_info);
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

	void superReactionNetwork::read_chem_out_spe_for_network_info(const std::string & cwd, std::vector<rsp::element_info>& element_v, std::vector<rsp::spe_info_base> &species_network_v, std::vector<rsp::reaction_info_base> &reaction_network_v, rsp::spe_name_index_map_t & spe_name_index_map, std::vector<VertexPair>& edgeVector, std::vector<EdgeProperties_graph>& edgePro, std::vector<VertexProperties_graph>& vertex_info, bool w2f)
	{
		/*
		* read species information
		*/
		std::vector<rsp::spe_info> species_v;
		//read element, species and reaction information
		rsp::relationshipParser::read_chem_out_ele_spe(element_v, species_v, spe_name_index_map, cwd + "/input/chem.out");

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
		rsp::relationshipParser::read_chem_out_reaction(species_v, reaction_v, spe_name_index_map, cwd + "/input/chem.out");
		rsp::relationshipParser::set_reaction_net_reactant_product(reaction_v);
		//reactionNetwork and chemkin index lookup table, notice chemkin index is Fortran index style
		rsp::reactionNetwork_chemkin_index_map_t reactionNetwork_chemkin_index_map;
		rsp::relationshipParser::read_reactionNetwork_chemkin_index_map(reactionNetwork_chemkin_index_map, cwd + "/input/chem.out");


		rsp::index_int_t edge_counter = 0;
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
			reaction_network_v.push_back(reaction_info_base_temp);

		}/*for*/

		if (w2f == true)
		{
			rsp::relationshipParser::spe_information_s2f(species_v);
			rsp::relationshipParser::spe_information_s2json(species_v);

			rsp::relationshipParser::reaction_information_s2f(species_v, reaction_v, reactionNetwork_chemkin_index_map);
			rsp::relationshipParser::reaction_information_s2json(species_v, reaction_v, reactionNetwork_chemkin_index_map);
		}
	}

	void superReactionNetwork::read_chem_out_spe_for_network_info(std::vector<VertexPair> &edgeVector, std::vector<EdgeProperties_graph> &edgePro, std::vector<VertexProperties_graph>& vertex_info)
	{

		read_chem_out_spe_for_network_info(this->cwd, this->element_v, this->species_network_v, this->reaction_network_v, this->spe_name_index_map, edgeVector, edgePro, vertex_info);
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
		rsp::index_int_t count = 0;
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
		for (rsp::index_int_t i = 0; i < (rsp::index_int_t)this->species_network_v.size(); ++i) {
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
		if (vertex >= (vertex_t)getVertexCount())
			return;
		for (out_edge_range_t itr = getOutEdges(vertex); itr.first != itr.second; ++itr.first) {
			reaction_index_s_coef_v.insert_sr(std::make_pair(properties(*itr.first).reaction_index, properties(*itr.first).s_coef_reactant));
		}

	}

	bool superReactionNetwork::search_for_out_spe(rsp::index_int_t reaction_index, std::vector<rsp::spe_index_weight_t>& out_spe_index_weight_v, std::string atom_followed)
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

	bool superReactionNetwork::set_reaction_out_spe_info(std::string atom_followed)
	{
		std::vector<rsp::spe_index_weight_t > out_spe_index_weight_v_tmp;
		for (rsp::index_int_t i = 0; i < (rsp::index_int_t)this->reaction_network_v.size(); ++i) {
			out_spe_index_weight_v_tmp.clear();
			search_for_out_spe(i, out_spe_index_weight_v_tmp, atom_followed);
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
		std::vector<rsp::index_int_t> spe_index(spe_rxn_c1_c2_map.size(), 0);

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

		const char* pattern1 = "(S\\d+(?:R[-]?\\d+)?)";
		boost::regex re1(pattern1);

		boost::sregex_iterator it1(pathway_in.begin(), pathway_in.end(), re1);
		boost::sregex_iterator end1;
		std::vector<std::string> reaction_spe;
		for (; it1 != end1; ++it1) {
			reaction_spe.push_back(it1->str());
		}

		const char* pattern2 = "S(\\d+)(?:R(-?\\d+))?";
		boost::regex re2(pattern2);

		for (size_t i = 0; i < reaction_spe.size(); i++)
		{
			std::vector<std::string> rxn_s_idx_str;

			boost::smatch result2;
			if (boost::regex_search(reaction_spe[i], result2, re2)) {
				for (std::size_t mi = 1; mi < result2.size(); ++mi) {
					string s_rxn_idx(result2[mi].first, result2[mi].second);
					rxn_s_idx_str.push_back(s_rxn_idx);
				}
			}

			spe_vec.push_back(boost::lexical_cast<rsp::index_int_t>(rxn_s_idx_str[0]));

			if (rxn_s_idx_str[1] != std::string(""))
				reaction_vec.push_back(boost::lexical_cast<rsp::index_int_t>(rxn_s_idx_str[1]));
			else {
				//not the last species
				if (i != reaction_spe.size() - 1)
					reaction_vec.push_back(INT_MAX);
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

	double reactionNetwork_sr::superReactionNetwork::spe_spe_branching_ratio(const std::vector<species_group_sr::rxn_c1_c2>& rxn_c1_c2_vec,
		double reaction_time, rsp::index_int_t curr_spe, rsp::index_int_t next_spe, std::string atom_followed, bool update_reaction_rate)
	{
		double ratio_tmp = 0.0;
		for (auto rxn_c1_c2 : rxn_c1_c2_vec) {
			auto reaction_index = rxn_c1_c2.r_idx;
			//whether to update reaction rates, is deferred to sub-routine to decice
			ratio_tmp += reaction_spe_branching_ratio(reaction_time, curr_spe, reaction_index, next_spe, atom_followed, update_reaction_rate);
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

	double superReactionNetwork::pathway_prob_input_pathway_sim_once(double const init_time, const double end_time, const std::vector<rsp::index_int_t> &spe_vec, const std::vector<rsp::index_int_t> &reaction_vec, std::string atom_followed)
	{
		//set pathway end time
		set_tau(end_time);

		//basically, we assume there must be a reaction at the beginning, so should multiply be the 1-P_min(tau=0|t;S^{0})
		double pathway_prob = 1.0;
		double when_time = init_time;

		//start from the first reaction
		for (size_t i = 0; i < reaction_vec.size();)
		{
			//none-chattering reaction
			if (reaction_vec[i] >= 0) {
				pathway_prob *= prob_spe_will_react_in_a_time_range(when_time, end_time, spe_vec[i]);
				when_time = pathway_prob_sim_move_one_step(when_time, spe_vec[i], reaction_vec[i], spe_vec[i + 1], pathway_prob, atom_followed);
				//move one step
				++i;
			}
			//chattering reaction, chattering case
			else {
				int chattering_group_id = this->species_network_v[spe_vec[i]].chattering_group_id;

				//add time delay first, regenerate random number, inverse to get exact time, get steady state time first
				//then calculate steady state ratios
				double chattering_group_prob = prob_chattering_group_will_react_in_a_time_range(when_time, end_time, chattering_group_id);
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
		set_spe_prob_max_at_a_time(when_time, end_time, spe_vec.back());

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

		pathway_prob *= spe_spe_branching_ratio(this->sp_all_species_group_rnk->out_species_rxns.at(curr_spe).at(next_spe),
			when_time, curr_spe, next_spe, atom_followed, true);

		return when_time;
	}

	double reactionNetwork_sr::superReactionNetwork::species_pathway_prob_input_pathway_sim_once(const double init_time, const double end_time, const std::vector<rsp::index_int_t>& spe_vec, const std::vector<rsp::index_int_t>& reaction_vec, std::string atom_followed)
	{
		//set pathway end time
		set_tau(end_time);

		//basically, we assume there must be a reaction at the beginning, so should multiply be the 1-P_min(tau=0|t;S^{0})
		double pathway_prob = 1.0;
		double when_time = init_time;

		//start from the first reaction
		for (size_t i = 0; i < reaction_vec.size();)
		{
			//none-chattering reaction
			if (reaction_vec[i] >= 0) {
				pathway_prob *= prob_spe_will_react_in_a_time_range(when_time, end_time, spe_vec[i]);
				when_time = species_pathway_prob_sim_move_one_step(when_time, spe_vec[i], spe_vec[i + 1], pathway_prob, atom_followed);
				//move one step
				++i;
			}
			//chattering reaction, chattering case
			else {
				int chattering_group_id = this->species_network_v[spe_vec[i]].chattering_group_id;

				//add time delay first, regenerate random number, inverse to get exact time, get steady state time first
				//then calculate steady state ratios
				double chattering_group_prob = prob_chattering_group_will_react_in_a_time_range(when_time, end_time, chattering_group_id);
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
					pathway_prob *= spe_spe_branching_ratio(this->sp_all_species_group_rnk->out_species_rxns.at(spe_vec[i + 1]).at(spe_vec[i + 2]),
						when_time, spe_vec[i + 1], spe_vec[i + 2], atom_followed, true);
					/*step 2*/

				}//boundary time problem

				 //move two steps actually
				i += 2;
			}//if chattering case

		}

		//got to multiply by P_min or says (1-P_max)
		set_spe_prob_max_at_a_time(when_time, end_time, spe_vec.back());

		pathway_prob *= (1 - species_network_v[spe_vec.back()].prob_max);

		return pathway_prob;
	}


	double superReactionNetwork::pathway_AT_sim_move_one_step(double when_time, vertex_t curr_spe)
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

		return when_time;
	}

	double superReactionNetwork::pathway_AT_input_pathway_sim_once(const double init_time, const double end_time, const std::vector<rsp::index_int_t>& spe_vec, const std::vector<rsp::index_int_t>& reaction_vec)
	{
		//set pathway end time
		set_tau(end_time);

		//basically, we assume there must be a reaction at the beginning, so should multiply be the 1-P_min(tau=0|t;S^{0})
		double when_time = init_time;

		//start from the first reaction
		for (size_t i = 0; i < reaction_vec.size();)
		{
			//none-chattering reaction
			if (reaction_vec[i] >= 0) {
				when_time = pathway_AT_sim_move_one_step(when_time, spe_vec[i]);
				//move one step
				++i;
			}
			//chattering reaction, chattering case
			else {
				int chattering_group_id = this->species_network_v[spe_vec[i]].chattering_group_id;

				//add time delay first, regenerate random number, inverse to get exact time, get steady state time first
				//then calculate steady state ratios
				double chattering_group_prob = prob_chattering_group_will_react_in_a_time_range(when_time, end_time, chattering_group_id);

				//avoid problems around boundary
				if (when_time < (tau - INFINITESIMAL_DT)) {
					double u_1 = 1.0;
					if (chattering_group_prob > 0.0) {
						u_1 = rand->random_min_max(0, chattering_group_prob);
					}
					else
						u_1 = 0.0;

					when_time = chattering_group_reaction_time_from_importance_sampling_without_cutoff(when_time, chattering_group_id, u_1);
				}//boundary time problem

				 //move two steps actually
				i += 2;
			}//if chattering case
		}

		return when_time;
	}

	double superReactionNetwork::pathway_AT_no_IT_input_pathway_sim_once(const double init_time, const double end_time, const std::vector<rsp::index_int_t>& spe_vec, const std::vector<rsp::index_int_t>& reaction_vec)
	{
		//set pathway end time
		set_tau(end_time);

		//basically, we assume there must be a reaction at the beginning, so should multiply be the 1-P_min(tau=0|t;S^{0})
		double when_time = init_time;

		bool is_IT = false;
		double IT = 0.0;

		//start from the first reaction
		for (size_t i = 0; i < reaction_vec.size();)
		{
			//none-chattering reaction
			if (reaction_vec[i] >= 0) {
				when_time = pathway_AT_sim_move_one_step(when_time, spe_vec[i]);

				if (is_IT == false) {
					is_IT = true;
					IT = when_time;
				}

				//move one step
				++i;
			}
			//chattering reaction, chattering case
			else {
				int chattering_group_id = this->species_network_v[spe_vec[i]].chattering_group_id;

				//add time delay first, regenerate random number, inverse to get exact time, get steady state time first
				//then calculate steady state ratios
				double chattering_group_prob = prob_chattering_group_will_react_in_a_time_range(when_time, end_time, chattering_group_id);

				//avoid problems around boundary
				if (when_time < (tau - INFINITESIMAL_DT)) {
					double u_1 = 1.0;
					if (chattering_group_prob > 0.0) {
						u_1 = rand->random_min_max(0, chattering_group_prob);
					}
					else
						u_1 = 0.0;

					when_time = chattering_group_reaction_time_from_importance_sampling_without_cutoff(when_time, chattering_group_id, u_1);

					if (is_IT == false) {
						is_IT = true;
						IT = when_time;
					}

				}//boundary time problem

				 //move two steps actually
				i += 2;
			}//if chattering case

		}

		return when_time - IT;
	}

	std::pair<double, double> superReactionNetwork::pathway_AT_with_SP_input_pathway_sim_once(const double init_time, const double end_time, const std::vector<rsp::index_int_t>& spe_vec, const std::vector<rsp::index_int_t>& reaction_vec)
	{
		//set pathway end time
		set_tau(end_time);

		//basically, we assume there must be a reaction at the beginning, so should multiply be the 1-P_min(tau=0|t;S^{0})
		double when_time = init_time;

		//start from the first reaction
		for (size_t i = 0; i < reaction_vec.size();)
		{
			//none-chattering reaction
			if (reaction_vec[i] >= 0) {
				when_time = pathway_AT_sim_move_one_step(when_time, spe_vec[i]);
				//move one step
				++i;
			}
			//chattering reaction, chattering case
			else {
				int chattering_group_id = this->species_network_v[spe_vec[i]].chattering_group_id;

				//add time delay first, regenerate random number, inverse to get exact time, get steady state time first
				//then calculate steady state ratios
				double chattering_group_prob = prob_chattering_group_will_react_in_a_time_range(when_time, end_time, chattering_group_id);

				//avoid problems around boundary
				if (when_time < (tau - INFINITESIMAL_DT)) {
					double u_1 = 1.0;
					if (chattering_group_prob > 0.0) {
						u_1 = rand->random_min_max(0, chattering_group_prob);
					}
					else
						u_1 = 0.0;

					when_time = chattering_group_reaction_time_from_importance_sampling_without_cutoff(when_time, chattering_group_id, u_1);
				}//boundary time problem

				 //move two steps actually
				i += 2;
			}//if chattering case
		}


		//got to multiply by P_min or says (1-P_max)
		set_spe_prob_max_at_a_time(when_time, end_time, spe_vec.back());

		return std::make_pair(when_time, 1 - species_network_v[spe_vec.back()].prob_max);
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

	void superReactionNetwork::heuristic_path_string_vector_si_sj_n_s2f(std::string atom_followed, std::size_t si, std::size_t sj, std::size_t n, std::string filename)
	{
		std::unordered_set<std::string> us;
		for (std::size_t k = 0; k <= n; ++k) {
			auto pRmn = matrix_sr::matrix_power(this->atom_R_matrix[atom_followed], k);

			if (k == 0) {
				us.insert(std::string("S") + boost::lexical_cast<std::string>(si));
				continue;
			}
			auto vs = this->get_path_string_update_matrix_element_i_j_topN(pRmn, si, sj);
			for (auto s : vs)
				us.insert(s);
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
