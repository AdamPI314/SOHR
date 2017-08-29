#ifndef __RELATIONSHIPPARSER_CPP_
#define __RELATIONSHIPPARSER_CPP_

#include "relationshipParser.h"

namespace relationshipParser_sr {
	void relationshipParser::read_chem_out_ele_spe(std::vector<element_info>& element_v, std::vector<spe_info>& species_v, spe_name_index_map_t& spe_name_index_map, std::string file_in)
	{
		std::ifstream fin(file_in.c_str());
		std::string str_t;
		/*
		* search for atom first, get name of all related atoms
		*/
		//std::vector<std::string> atom_vec;
		element_info element_tp;

		//	const char* pattern1= "^\\s*\\d+\\.\\s(\\w+)\\s+\\d+\\.\\d+\\s*$";
		const char* pattern1 = "^\\s*\\d+\\.\\s(\\w+)\\s+(\\d+(?:\\.\\d+)?)\\s*$";
		//	boost::regex regexPattern1(pattern1, boost::regex::extended);
		boost::regex regexPattern1(pattern1);
		boost::smatch what;
		while (!fin.eof()) {
			std::getline(fin, str_t);
			bool isMatchFound = boost::regex_search(str_t, what, regexPattern1);
			if (isMatchFound)
			{
				//			for (size_t i=1; i < what.size(); i++)
				//			{
				//				//std::cout << "WHAT " << i << " " << what[i] << std::endl;
				//				atom_vec.push_back(what[i]);
				//			}
				//what[0] is the whole string that match
				//			std::cout<<what[1]<<"-->"<<boost::lexical_cast<atomic_weight_t>(what[2])<<"\n";
				element_tp.ele_name = what[1]; element_tp.atomic_weight = boost::lexical_cast<atomic_weight_t>(what[2]);
				element_v.push_back(element_tp);
			}

		}
		for (std::size_t i = 0; i < element_v.size(); ++i) {
			element_v[i].ele_index = i;
			//		std::cout<<"ele_index-->"<<element_v[i].ele_index<<"\n";
			//		std::cout<<"ele_name-->"<<element_v[i].ele_name<<"\n";
			//		std::cout<<"atomic weight-->"<<element_v[i].atomic_weight<<"\n";

		}

		//std::copy(atom_vec.begin(), atom_vec.end(), std::ostream_iterator<std::string>(std::cout, "\t")); std::cout<<std::endl;
		/*
		* search for species and get species information, like species component
		*/
		//check this regular expression pattern every time you change the "chem.out"
		std::string repeated_group;
		for (std::size_t i = 0; i < element_v.size(); ++i)
		{
			repeated_group += std::string("(\\d+)\\s*");
		}
		//regular expression with repeated group
		std::string pattern_tmp = std::string("\\s*\\d+\\.\\s([\\w\\(\\)\\-\\_,]+)\\s+([GLS])\\s(\\d+|\\d+\\.\\d+)\\s+(\\d+\\.\\d+)\\s+(\\d+\\.\\d+)\\s(\\d+\\.\\d+)\\s+") + repeated_group;
		const char* pattern2 = pattern_tmp.c_str();
		//boost::regex regexPattern2(pattern2, boost::regex::extended);
		//cannot use extended mode here, use the default mode
		//http://stackoverflow.com/questions/20107388/boostregex-usage-in-c
		boost::regex regexPattern2(pattern2);

		fin.close(); fin.clear();
		fin.open(file_in.c_str());

		while (!fin.eof())
		{//[while
			spe_info species_t;
			std::getline(fin, str_t);
			bool isMatchFound = boost::regex_search(str_t, what, regexPattern2);
			if (isMatchFound)
			{
				//the first matched segment is the whole string
				//std::cout<<what[0]<<"\n";

				for (size_t i = 1; i < what.size(); i++)
				{
					if (i == 1) {
						species_t.spe_name = what[i];
					}
					else if (i == 2) {
						species_t.phase = what[i];
					}
					else if (i == 3) {
						species_t.charge = boost::lexical_cast<charge_t>(what[i]);
					}
					else if (i == 4) {
						species_t.spe_weight = boost::lexical_cast<spe_mass_weight_t>(what[i]);
					}
					else if (i == 5) {
						species_t.temperature_low = boost::lexical_cast<temperature_t>(what[i]);
					}
					else if (i == 6) {
						species_t.temperature_high = boost::lexical_cast<temperature_t>(what[i]);
					}
					//else if(what[i]!=std::string("0")){
					//0,1,2,3,4,5,6 are all taken care of
					else {
						species_t.spe_component[element_v[i - 7].ele_name] = boost::lexical_cast<std::size_t>(what[i]);
					}
				}//for

				 //if founded
				 //insert species
				species_v.push_back(species_t);
			}//if

		}//]while

		for (std::size_t i = 0; i < species_v.size(); ++i)
		{
			//store species name and species index in a map
			spe_name_index_map[species_v[i].spe_name] = i;

			//species index
			species_v[i].spe_index = i;
		}



		fin.close(); fin.clear();

	}/*void relationshipParser::read_chem_out_ele_spe*/


	void relationshipParser::spe_information_s2f(const  std::vector<spe_info>& species_v, std::string file_out)
	{
		//print out the reactions
		std::ofstream fout_t(file_out.c_str());
		for (std::size_t i = 0; i < species_v.size(); ++i)
		{
			fout_t << i << "\t-->\t" << species_v[i].spe_name << "\t";
			fout_t << "phase-->" << species_v[i].phase << "\t";
			fout_t << "charge-->" << species_v[i].charge << "\t";
			fout_t << "molecular weight-->" << species_v[i].spe_weight << "\t";
			fout_t << "temperature low-->" << species_v[i].temperature_low << "\t";
			fout_t << "temperature high-->" << species_v[i].temperature_high << "\t";
			for (spe_component_t::const_iterator itr = species_v[i].spe_component.begin(); itr != species_v[i].spe_component.end(); ++itr) {
				fout_t << itr->first << ":" << itr->second << "\t";
			}
			fout_t << std::endl;
		}

		fout_t.close(); fout_t.clear();
	}


	void relationshipParser::read_chem_out_reaction(const std::vector<spe_info>& species_v, std::vector<reaction_info>& reaction_v, const spe_name_index_map_t& spe_name_index_map, std::string file_in)
	{
		std::ifstream fin(file_in.c_str());
		std::string str_t;
		reaction_info reaction_info_tp;
		/*
		* search reactions and their rate coefficient
		*/
		std::string pattern1_reaction_body = "([\\w+\\(\\)\\-\\_,]+(?:=>|<=|=|<=>)+[\\w+\\(\\)\\-\\_,]+)\\s+([\\d.E+\\-]+)\\s+([\\d.+\\-]+)\\s+([\\d.+\\-]+)";
		std::string pattern1_reaction_Enhanced_by = "([\\w+\\(\\)\\-\\_,]+)\\s+(Enhanced by)\\s+([\\d.E+\\-]+)";
		std::string pattern1_reaction_low_pressure_limit = "(Low pressure limit):\\s+([\\d.E+\\-]+)\\s+([\\d.E+\\-]+)\\s+([\\d.E+\\-]+)";
		std::string pattern1_reaction_TROE_centering = "(TROE centering):\\s+([\\d.E+\\-]+)\\s+([\\d.E+\\-]+)\\s+([\\d.E+\\-]+)";

		std::string pattern1_reaction_info = pattern1_reaction_body + "|" +
			pattern1_reaction_Enhanced_by + "|" +
			pattern1_reaction_low_pressure_limit + "|" +
			pattern1_reaction_TROE_centering;
		boost::regex regexPattern1_reaction_info(pattern1_reaction_info);
		boost::smatch what;

		/*
		* for the match tokens
		* 1,2,3,4 match reaction body		std::cout<<"weight\t"<<reaction_v[reaction_ind].reactant[i].second<<std::endl;
		*
		* 5,6,7 match enhanced by
		* 8,9,10,11 match low pressure limit
		* 12,13,14,15 match TROE centering
		*/
		while (!fin.eof()) {
			std::getline(fin, str_t);
			bool isMatchFound = boost::regex_search(str_t, what, regexPattern1_reaction_info);
			if (isMatchFound) {
				if (std::string(what[0]).find("=") != std::string::npos) {
					//find a new reaction, have to save the old reaction first
					if (reaction_info_tp.reaction_name != "")
						reaction_v.push_back(reaction_info_tp);
					//reinitialize reaction info second
					reaction_info_tp.prefactor_A = boost::lexical_cast<prefactor_A_type_t>(what[2]);
					reaction_info_tp.delta = boost::lexical_cast<delta_t>(what[3]);
					reaction_info_tp.barrier_E = boost::lexical_cast<barrier_E_t>(what[4]);
					reaction_info_tp.reaction_rate = 0.0;
					reaction_info_tp.isPressureDependentReaction = false;
					reaction_info_tp.isThreeBodyReaction = false;
					reaction_info_tp.spe_name_3body_enhancement_coef.clear();

					//				std::cout<<what[0]<<std::endl;
					//				std::cout<<what[1]<<std::endl;
					reaction_info_tp.reaction_name = what[1];
				}

				else if (std::string(what[0]).find("Enhanced by") != std::string::npos) {
					//				std::cout<<what[0]<<std::endl;
					//				std::cout<<what[5]<<std::endl;
					//				std::cout<<what[6]<<std::endl;
					//				std::cout<<what[7]<<std::endl;
					reaction_info_tp.spe_name_3body_enhancement_coef[boost::lexical_cast<spe_name_t>(what[5])] = boost::lexical_cast<stoichiometric_coef_t>(what[7]);
				}
				else if (std::string(what[0]).find("Low pressure limit") != std::string::npos) {
					//				std::cout<<what[0]<<std::endl;
					//				std::cout<<what[8]<<std::endl;
					//				std::cout<<what[9]<<std::endl;
					//				std::cout<<what[10]<<std::endl;
					//				std::cout<<what[11]<<std::endl;
					reaction_info_tp.low_pressure_limit[0] = boost::lexical_cast<double>(what[9]);
					reaction_info_tp.low_pressure_limit[1] = boost::lexical_cast<double>(what[10]);
					reaction_info_tp.low_pressure_limit[2] = boost::lexical_cast<double>(what[11]);
				}
				else if (std::string(what[0]).find("TROE centering") != std::string::npos) {
					//				std::cout<<what[0]<<std::endl;
					//				std::cout<<what[12]<<std::endl;
					//				std::cout<<what[13]<<std::endl;
					//				std::cout<<what[14]<<std::endl;
					//				std::cout<<what[15]<<std::endl;
					reaction_info_tp.TROE_centering[0] = boost::lexical_cast<double>(what[13]);
					reaction_info_tp.TROE_centering[1] = boost::lexical_cast<double>(what[14]);
					reaction_info_tp.TROE_centering[2] = boost::lexical_cast<double>(what[15]);
				}

			}//if

		}//while
		 //add the last reaction
		reaction_v.push_back(reaction_info_tp);

		//parse species to species transition, A->B, record the transition reaction, s_coef_reactant and s_coef_product
		const char* pattern2_reactant = "([\\w\\(\\)\\-\\_,]+)(\\+|[\\w\\(\\)\\-\\_,]+)*(?:<=>|=>|<=|=)"; boost::regex regexPattern2_reactant(pattern2_reactant);
		const char* pattern2_product = "(?:<=>|=|=>|<=)([\\w\\(\\)\\-\\_,]+)(?:\\+([\\w\\(\\)\\-\\_,]+))*"; boost::regex regexPattern2_product(pattern2_product);
		const char* pattern2_spe = "(?!\\bM\\b)((?:[\\w\\-\\_,]|\\(\\w\\))+)"; boost::regex regexPattern2_spe(pattern2_spe);
		//if matched species start from a number, parse the number the the species, like 2H, 2H2O2
		const char* pattern2_num_spe = "^(\\d+)([\\w\\(\\)\\-\\_,]+)"; boost::regex regexPattern2_num_spe(pattern2_num_spe);


		//	std::size_t edge_counter=0;
		for (std::size_t i = 0; i < reaction_v.size(); ++i) {
			reaction_v[i].reaction_index = i;
			/*
			* for every reaction, decompose it to A->B transition, record stoichiometric coefficient
			*/

			//find reactants first
			boost::sregex_token_iterator itr2_reactant(reaction_v[i].reaction_name.begin(), reaction_v[i].reaction_name.end(), regexPattern2_reactant, 0);
			boost::sregex_token_iterator end;
			std::string str_tmp = *(itr2_reactant);
			std::string str_tmp_r = str_tmp;
			//std::cout<<str_tmp<<std::endl;
			// decompose species which are reactants
			spe_name_s_coef_t spe_name_s_coef_reactant;
			boost::sregex_token_iterator itr2_reactant_spe(str_tmp.begin(), str_tmp.end(), regexPattern2_spe, 0);
			boost::smatch what;
			for (; itr2_reactant_spe != end; ++itr2_reactant_spe) {
				//parse species like 2H, 2H2O2
				std::string str_tmp2 = (*itr2_reactant_spe);
				if (boost::regex_search(str_tmp2, what, regexPattern2_num_spe))
					spe_name_s_coef_reactant[what[2]] += boost::lexical_cast<std::size_t>(what[1]);
				else
					spe_name_s_coef_reactant[(*itr2_reactant_spe)] += 1;
			}

			//find products then
			boost::sregex_token_iterator itr2_product(reaction_v[i].reaction_name.begin(), reaction_v[i].reaction_name.end(), regexPattern2_product, 0);
			str_tmp = *(itr2_product);
			//std::cout<<str_tmp<<std::endl;
			// decompose species which are products
			spe_name_s_coef_t spe_name_s_coef_product;
			boost::sregex_token_iterator itr2_product_spe(str_tmp.begin(), str_tmp.end(), regexPattern2_spe, 0);
			for (; itr2_product_spe != end; ++itr2_product_spe) {
				//parse species like 2H, 2H2O2
				std::string str_tmp2 = (*itr2_product_spe);
				if (boost::regex_search(str_tmp2, what, regexPattern2_num_spe)) {
					spe_name_s_coef_product[what[2]] += boost::lexical_cast<std::size_t>(what[1]);
				}
				else
					spe_name_s_coef_product[(*itr2_product_spe)] += 1;
			}

			//we finded the reactants and products, add them to reaction information
			for (spe_name_s_coef_t::iterator itr = spe_name_s_coef_reactant.begin(); itr != spe_name_s_coef_reactant.end(); ++itr) {
				//reaction_v[i].reactant.push_back(spe_index_weight_t(spe_name_index_map[itr->first], itr->second));
				spe_name_index_map_t::const_iterator itr_m_t = spe_name_index_map.find(itr->first);
				if (itr_m_t != spe_name_index_map.end())
					reaction_v[i].reactant.push_back(spe_index_weight_t(itr_m_t->second, itr->second));
			}

			for (spe_name_s_coef_t::iterator itr = spe_name_s_coef_product.begin(); itr != spe_name_s_coef_product.end(); ++itr) {
				//reaction_v[i].product.push_back(spe_index_weight_t(spe_name_index_map[itr->first], itr->second));
				spe_name_index_map_t::const_iterator itr_m_t = spe_name_index_map.find(itr->first);
				if (itr_m_t != spe_name_index_map.end())
					reaction_v[i].product.push_back(spe_index_weight_t(itr_m_t->second, itr->second));
			}

		}//for

		fin.close(); fin.clear();

		//parse reaction direction
		for (std::size_t i = 0; i < reaction_v.size(); ++i) {
			parse_reaction_direction(reaction_v[i]);
		}

	}/*void relationshipParser::read_chem_out_reaction*/

	void relationshipParser::reaction_information_s2f(const std::vector<spe_info> & species_v, const std::vector<reaction_info>& reaction_v, const reactionNetwork_chemkin_index_map_t& reactionNetwork_chemkin_index_map, std::string file_out)
	{
		//print out the reactions
		std::ofstream fout(file_out.c_str());

		for (reactionNetwork_chemkin_index_map_t::const_iterator itr = reactionNetwork_chemkin_index_map.begin(); itr != reactionNetwork_chemkin_index_map.end(); ++itr) {
			fout << itr->first << "\t";
			for (index_int_t i = 0; i < itr->second.size(); ++i) {
				//need to convert Fortran style index into C++ style index
				index_int_t reaction_v_ind = static_cast<index_int_t>(abs(itr->second[i])) - 1;
				fout << itr->second[i] << "\t";
				fout << reaction_v[reaction_v_ind].reaction_name << "\t";
			}
			fout << std::endl;

			//just print the first reaction
			index_int_t reaction_v_ind = static_cast<index_int_t>(abs(itr->second[0])) - 1;
			fout << "reactants\t";
			for (index_int_t i = 0; i < reaction_v[reaction_v_ind].reactant.size(); ++i) {
				fout << species_v[reaction_v[reaction_v_ind].reactant[i].first].spe_name << "\t";
			}
			fout << "products\t";
			for (index_int_t i = 0; i < reaction_v[reaction_v_ind].product.size(); ++i) {
				fout << species_v[reaction_v[reaction_v_ind].product[i].first].spe_name << "\t";
			}
			fout << std::endl;

			fout << "net_reactants\t";
			for (index_int_t i = 0; i < reaction_v[reaction_v_ind].net_reactant.size(); ++i) {
				fout << species_v[reaction_v[reaction_v_ind].net_reactant[i].first].spe_name << "\t";
			}
			fout << "net_products\t";
			for (index_int_t i = 0; i < reaction_v[reaction_v_ind].net_product.size(); ++i) {
				fout << species_v[reaction_v[reaction_v_ind].net_product[i].first].spe_name << "\t";
			}
			fout << std::endl;
		}

		fout.close(); fout.clear();

	}

	//regard reactions only as forward reactions
	void relationshipParser::set_spe_sk_info(std::vector<spe_info>& species_v, const std::vector<reaction_info>& reaction_v)
	{
		//	std::cout<<"HAH"<<std::endl;
		for (std::size_t i = 0; i < species_v.size(); ++i) {
			//		std::cout<<i<<"\t"<<species_v[i].spe_name<<"\t"<<species_v[i].spe_index<<"\n";
			for (std::size_t j = 0; j < reaction_v.size(); ++j) {

				for (std::size_t k = 0; k < reaction_v[j].reactant.size(); ++k) {
					//find species in reaction as reactant
					if (species_v[i].spe_index == reaction_v[j].reactant[k].first) {
						species_v[i].reaction_sk_index_s_coef_v.push_back(reaction_sk_index_s_coef_t(-1,
							reaction_index_s_coef_t(reaction_v[j].reaction_index, reaction_v[j].reactant[k].second)));
					}//if
				}//for k

				for (std::size_t k2 = 0; k2 < reaction_v[j].product.size(); ++k2) {
					//find species in reaction as product
					if (species_v[i].spe_index == reaction_v[j].product[k2].first) {
						species_v[i].reaction_sk_index_s_coef_v.push_back(reaction_sk_index_s_coef_t(+1,
							reaction_index_s_coef_t(reaction_v[j].reaction_index, reaction_v[j].product[k2].second)));
					}//if
				}//for k2

			}//for j
		}//for i

		 //net source term
		for (std::size_t i = 0; i < species_v.size(); ++i) {
			for (std::size_t j = 0; j < species_v[i].reaction_sk_index_s_coef_v.size(); ++j) {
				//source term
				if (species_v[i].reaction_sk_index_s_coef_v[j].first == +1) {
					reaction_index_s_coef_t reaction_index_s_coef_tp = species_v[i].reaction_sk_index_s_coef_v[j].second;
					//in case there is a sink term
					for (std::size_t k = 0; k < species_v[i].reaction_sk_index_s_coef_v.size(); ++k) {//k
						if (species_v[i].reaction_sk_index_s_coef_v[k].first == -1 &&
							species_v[i].reaction_sk_index_s_coef_v[k].second.first == reaction_index_s_coef_tp.first) {
							reaction_index_s_coef_tp.second -= species_v[i].reaction_sk_index_s_coef_v[k].second.second;
						}//if
					}//k
					if (reaction_index_s_coef_tp.second > 0) {//if
						species_v[i].reaction_s_index_s_coef_v.push_back(reaction_index_s_coef_tp);
					}//if

				}//if
			}//j
		}//i

		 //net sink term
		for (std::size_t i = 0; i < species_v.size(); ++i) {
			for (std::size_t j = 0; j < species_v[i].reaction_sk_index_s_coef_v.size(); ++j) {
				//sink term
				if (species_v[i].reaction_sk_index_s_coef_v[j].first == -1) {
					reaction_index_s_coef_t reaction_index_s_coef_tp = species_v[i].reaction_sk_index_s_coef_v[j].second;
					//in case there is a source term
					for (std::size_t k = 0; k < species_v[i].reaction_sk_index_s_coef_v.size(); ++k) {//k
						if (species_v[i].reaction_sk_index_s_coef_v[k].first == +1 &&
							species_v[i].reaction_sk_index_s_coef_v[k].second.first == reaction_index_s_coef_tp.first) {
							reaction_index_s_coef_tp.second -= species_v[i].reaction_sk_index_s_coef_v[k].second.second;
						}//if
					}//k
					if (reaction_index_s_coef_tp.second > 0) {//if
						species_v[i].reaction_k_index_s_coef_v.push_back(reaction_index_s_coef_tp);
					}//if

				}//ifspecies_v[spe_ind].reaction_k_index_s_coef_v
			}//j
		}//i

	}
	/*void relationshipParser::set_spe_sk_info*/

	void relationshipParser::set_reaction_net_reactant_product(reaction_info & reaction)
	{
		//net reactant
		for (std::size_t i = 0; i < reaction.reactant.size(); ++i) {
			bool in_product = false;
			//search to see whether it is in product
			for (std::size_t j = 0; j < reaction.product.size(); ++j) {
				//find it in the product
				if (reaction.reactant[i].first == reaction.product[j].first) {
					in_product = true;
					//if it is a net reactant
					if (reaction.reactant[i].second > reaction.product[j].second)
						reaction.net_reactant.push_back(std::make_pair(reaction.reactant[i].first, reaction.reactant[i].second - reaction.product[j].second));

					//can at most found once
					break;
				}
			}
			//not in product
			if (in_product == false)
				reaction.net_reactant.push_back(reaction.reactant[i]);
		}

		//net product
		for (std::size_t i = 0; i < reaction.product.size(); ++i) {
			bool in_reactant = false;
			//search to see whether it is in reactant
			for (std::size_t j = 0; j < reaction.reactant.size(); ++j) {
				//find in reactant
				if (reaction.product[i].first == reaction.reactant[j].first) {
					in_reactant = true;
					//if it is a net product
					if (reaction.product[i].second > reaction.reactant[j].second)
						reaction.net_product.push_back(std::make_pair(reaction.product[i].first, reaction.product[i].second - reaction.reactant[j].second));

					//can at most found once
					break;
				}
			}
			//not found in reactant
			if (in_reactant == false)
				reaction.net_product.push_back(reaction.product[i]);
		}

	}


	void relationshipParser::set_reaction_net_reactant_product(std::vector<reaction_info>& reaction_v)
	{
		for (std::size_t i = 0; i < reaction_v.size(); ++i) {
			relationshipParser::set_reaction_net_reactant_product(reaction_v[i]);
		}

	}

	void relationshipParser::read_reactionNetwork_chemkin_index_map(reactionNetwork_chemkin_index_map_t & reactionNetwork_chemkin_index_map, std::string infile)
	{
		std::ifstream fin(infile.c_str());
		std::string str_t;

		//reaction pattern without index
		const char* pattern1 = "[\\w+\\(\\)\\-\\_,]+(?:=>|<=|=|<=>)+[\\w+\\(\\)\\-\\_,]+";
		boost::regex regexPattern1(pattern1);
		const char* pattern1_index = "\\d+\\.\\s+[\\w+\\(\\)\\-\\_,]+(?:=>|<=|=|<=>)+[\\w+\\(\\)\\-\\_,]+";
		boost::regex regexPattern1_index(pattern1_index);
		//reaction pattern with index

		//all reactions, no duplicated reactions
		mt::vector_sr<std::string> vec_sr;
		//all reactions, with duplicated reactions, vector of vector of std::string
		std::vector<std::vector<std::string> > vec_vec_str;

		while (!fin.eof())
		{
			std::getline(fin, str_t);
			boost::sregex_iterator itr(str_t.begin(), str_t.end(), regexPattern1);
			boost::sregex_iterator itr_index(str_t.begin(), str_t.end(), regexPattern1_index);
			boost::sregex_iterator end;
			for (; (itr != end) && (itr_index != end); ++itr, ++itr_index)
			{
				//get the reaction like:HO2+H=H2+O2
				//the element already exist.
				//vec_sr.insert_sr(itr->str());
				//not found
				if (std::find(vec_sr.begin(), vec_sr.end(), itr->str()) == vec_sr.end()) {
					vec_sr.push_back(itr->str());
					std::vector<std::string> v_str = mt::make_vector<std::string>() << itr_index->str();
					vec_vec_str.push_back(v_str);
				}
				//found
				else {
					//				std::cout<<"Found:\t"<<itr_index->str()<<"\n";
					vec_vec_str.back().push_back(itr_index->str());
				}
			}
		}//while

		std::size_t edge_counter = 0;
		for (std::size_t i = 0; i < vec_sr.size(); ++i) {
			//std::cout<<i<<":\t"<<vec_sr[i]<<std::endl;

			//parse reaction ChemKin style reaction index
			std::vector<int> vec_index;
			const char* pattern2_index = "(\\d+)\\..+"; boost::regex regexPattern2_index(pattern2_index);
			boost::smatch what;
			for (size_t j = 0; j < vec_vec_str[i].size(); ++j) {
				if (boost::regex_search(vec_vec_str[i][j], what, regexPattern2_index)) {
					vec_index.push_back(boost::lexical_cast<int>(what[1]));
				}
			}


			//find reaction arrow, <=> or + or => or <=
			const char* pattern2_arrow = "(<=>|=>|<=|=)"; boost::regex regexPattern2_arrow(pattern2_arrow);
			boost::sregex_token_iterator itr2_arrow(vec_sr[i].begin(), vec_sr[i].end(), regexPattern2_arrow, 0);

			if ((*itr2_arrow) == std::string("<=>") || (*itr2_arrow) == std::string("=")) {
				//std::cout<<*itr2_arrow<<"\t"<<edge_counter<<"\t"<<edge_counter+1<<std::endl;
				reactionNetwork_chemkin_index_map[edge_counter] = vec_index; //forward reaction index
				for (size_t k = 0; k < vec_index.size(); ++k) {
					vec_index[k] *= (-1);
				}
				reactionNetwork_chemkin_index_map[edge_counter + 1] = vec_index; //backward reaction index

				edge_counter += 2; //two reactions, forward and backward
			}//if

			else if ((*itr2_arrow) == std::string("=>")) {
				//std::cout<<*itr2_arrow<<"\t"<<edge_counter<<std::endl;
				reactionNetwork_chemkin_index_map[edge_counter] = vec_index; //forward reaction index
				++edge_counter; //only forward reaction
			}//else if

			else if ((*itr2_arrow) == std::string("<=")) {
				//std::cout<<*itr2_arrow<<"\t"<<edge_counter<<std::endl;
				for (size_t k = 0; k < vec_index.size(); ++k) {
					vec_index[k] *= -1;
				}
				reactionNetwork_chemkin_index_map[edge_counter + 1] = vec_index; //backward reaction index
				++edge_counter; //only backward reaction
			}//else if


		}//for

		fin.close(); fin.clear();
	}	/*void relationshipParser::set_spe_sk_info*/

	void relationshipParser::read_chem_out_duplicated_reaction(std::vector<int>& duplicate_reaction_ind, std::string str)
	{
		std::ifstream fin(str.c_str());
		std::string str_t;

		const char* pattern = "[\\w|+|(|)]+=[\\w|+|(|)]+";
		boost::regex re(pattern);

		//all reactions
		mt::vector_sr<std::string> vec_sr;
		int count = 0;

		while (!fin.eof())
		{
			std::getline(fin, str_t);
			boost::sregex_iterator it(str_t.begin(), str_t.end(), re);
			boost::sregex_iterator end;
			for (; it != end; ++it)
			{
				//get the reaction like:HO2+H=H2+O2
				//the element already exist.
				if (std::find(vec_sr.begin(), vec_sr.end(), it->str()) != vec_sr.end())
					duplicate_reaction_ind.push_back(count);
				vec_sr.insert_sr(it->str());
				++count;
			}
		}
		fin.clear(); fin.close();

	}

	void relationshipParser::parse_reaction_direction(reaction_info & reaction)
	{

		//find reaction arrow, <=> or = or => or <=
		const char* pattern2_arrow = "(<=>|=>|<=|=)"; boost::regex regexPattern2_arrow(pattern2_arrow);
		boost::sregex_token_iterator itr2_arrow(reaction.reaction_name.begin(), reaction.reaction_name.end(), regexPattern2_arrow, 0);
		if ((*itr2_arrow) == std::string("<=>") || (*itr2_arrow) == std::string("=")) {
			reaction.reaction_direction = both;
		}//if

		else if ((*itr2_arrow) == std::string("=>")) {
			reaction.reaction_direction = forward;
		}//else if

		else if ((*itr2_arrow) == std::string("<=")) {
			reaction.reaction_direction = backward;
		}//else if
	}



}/*namespace relationshipParser_sr*/

#endif
