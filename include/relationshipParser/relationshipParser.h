#ifndef __RELATIONSHIPPARSER_H_
#define __RELATIONSHIPPARSER_H_

#include <map> //for std::map
#include <vector> //for std::vector
#include <string> //for std::string
#include <iostream>
#include <fstream>

#include <boost/regex.hpp> //for boost::regex
#include <boost/lexical_cast.hpp> //for boost::regex
#include <boost/algorithm/string.hpp> //for boost::contains for string, to tell whether a string contains a substring or not



#include "../tools/misc/misc_template.h"
#include "../tools/misc/graph_bundled.h"

namespace relationshipParser_sr {

	namespace mt = misc_template;

	//time type
	typedef double my_time_t;
	//index type or count number/integer
	typedef int index_int_t;

	//a constant, used as indicator
	const index_int_t INDICATOR = 100000;

	/*
	* reaction network reaction index and ChemKin reaction index lookup table, the reason the ChemKin index is a vector
	* is there might exist duplicated reactions ChemKin has its own index of reaction, our reaction network redefine
	* reaction index, here is the corresponding relation between them
	* The former is C/C++ style index, The later is Fortran Style index, actually they are exact index in ChemKin space
	* negative index represents backward reaction
	*/
	typedef std::map<index_int_t, std::vector<int> > reactionNetwork_chemkin_index_map_t;

	//stoichiometric coefficient type
	typedef double stoichiometric_coef_t;
	//concentration type
	typedef double concentration_t;
	//pressure type
	typedef double pressure_t;
	//volume type
	typedef double volume_t;
	//temperature type
	typedef double temperature_t;
	//phase type
	typedef std::string phase_t;
	//charge type
	typedef double charge_t;

	/*
	* k = A T**b exp(-E/RT)
	* NOTE:  A units mole-cm-sec-K, E units cal/mole
	*/
	//rate constant type
	typedef double rate_constant_t;
	//pre-factor A type
	typedef double prefactor_A_type_t;
	//b, or says delta type
	typedef double delta_t;
	//Energy or says barrier type
	typedef double barrier_E_t;
	//reaction rate type
	typedef double reaction_rate_t;
	//typedef reaction name type
	typedef std::string reaction_name_t;

	//elements name
	typedef std::string element_name_t;
	//atomic weight type
	typedef double atomic_weight_t;

	//species name type
	typedef std::string spe_name_t;
	//species component type, like total number of atoms, number of a specific atom
	typedef std::map<std::string, index_int_t> spe_component_t;

	//species name and species stoichoimetric coefficient in reaction/transition.
	typedef std::map<spe_name_t, stoichiometric_coef_t> spe_name_s_coef_t;
	//species name and the third body enhancement coefficient for three-body reaction or pressure dependent reaction
	typedef std::map<spe_name_t, stoichiometric_coef_t> spe_name_3body_enhancement_coef_t;
	/*to pressure dependent reaction
	* there are low pressure limit and
	* TROE centering
	*/
	typedef double low_pressure_limit_t[3];
	typedef double TROE_centering_t[3];
	//sink or source reactions of a species and stoichoimetric coefficient
	typedef std::pair<index_int_t, stoichiometric_coef_t> reaction_index_s_coef_t;
	/*
	* All reactions a species involved in, either sink or source reactions
	* The first number, +1/-1 , to indicate whether it is a source or sink reaction
	* sk-->represents source or sink
	*/
	typedef std::pair<index_int_t, reaction_index_s_coef_t> reaction_sk_index_s_coef_t;
	//species mass weight
	typedef double spe_mass_weight_t;

	//species index and weight type
	typedef std::pair<index_int_t, stoichiometric_coef_t> spe_index_weight_t;

	//element name and element index type map
	typedef std::map<element_name_t, index_int_t> ele_name_index_map_t;
	//species name and species index type map
	typedef std::map<spe_name_t, index_int_t> spe_name_index_map_t;

	//element information struct
	struct element_info {
		index_int_t ele_index;
		element_name_t ele_name;
		atomic_weight_t atomic_weight;

	};


	//probability type
	typedef double probability_t;

	//species information base type, define basic information we are gonna to use in the reaction network
	struct spe_info_base {
		index_int_t spe_index;
		spe_name_t spe_name;
		//species component
		spe_component_t spe_component;
		//species concentration
		concentration_t spe_conc;
		//more reaction network related stuff, like survival probability etc
		probability_t survival_probability;
		//P_min for every species, where P_min is defined so that when P<P_min in a certain time range, no reaction will occur
		//P_min in principle, tell us about the time limitation that reaction has to occur in a certain time range, otherwise, 
		//it is out of intertest
		probability_t prob_min;
		//P_max= 1-P_min
		probability_t prob_max;

		//fast transition group/ chattering group
		int chattering_group_id;

		//net sink terms of this species, cancel the same terms on both sides for auto-catalytic reactions, left sink term
		std::vector<reaction_index_s_coef_t> reaction_k_index_s_coef_v;

		//std::

		spe_info_base() :spe_conc(0.0), survival_probability(0.0), prob_min(0.0), prob_max(0.0), chattering_group_id(-1) {}
	};

	//species information struct
	struct spe_info :public spe_info_base {

		spe_mass_weight_t spe_weight;
		phase_t phase;
		charge_t charge;
		//temperature low
		temperature_t temperature_low;
		//temperature high
		temperature_t temperature_high;

		//reactions it is involved in, source and sink reactions and the stoichiometric coefficient
		/*
		* out reaction and stoichoimetric coefficient of this species
		* like 2A+B->C+2D, stoichoimetric coefficient is 2 for A of this reaction
		* the reason to store this information is it is a fixed, unchangeable relation, no need to search reaction network
		* as soon as it is initialized and need it in pathway generation scheme
		*/
		std::vector<reaction_sk_index_s_coef_t> reaction_sk_index_s_coef_v;
		//net source term, cancel the same terms on both sides for auto-catalytic reactions, left source term
		//got to consider auto-catalytic reactions.
		/*
		* A+X->2X
		* this is not a out reaction for X, a source term instead of a sink term
		*/
		std::vector<reaction_index_s_coef_t> reaction_s_index_s_coef_v;

	};


	//reaction direction
	enum reaction_direction_t { forward = 1, backward = -1, both = 2 };
	//reaction information struct base
	//simplified version, used in reaction network
	struct reaction_info_base {
		reaction_name_t reaction_name;
		reaction_direction_t reaction_direction;
		reaction_rate_t reaction_rate;

		// is reaction zero, used as a flag
		// zero->false, none zero->true
		bool is_reaction_rate_nonzero = false;
		/*
		* out species and corresponding weight
		* like 2A+B->H2O2+H100
		* two out species, weight wrt. H atom is: 2, 100
		* the reason to store this information is it is a fixed, unchangeable relation, no need to search reaction network
		* as soon as it is initialized and need it in pathway generation scheme
		*/
		//might be problematic, got to consider auto-catalytic reactions.
		/*
		* A+X->2X
		* X is a out species for this reaction, and the out weight is just 1 instead of 2
		*/
		std::map<std::string, std::vector<spe_index_weight_t > >  out_spe_index_weight_v_map;

		// out species, branching ratio, by follow a atom
		std::map<std::string, std::map<std::size_t, double> > out_spe_index_branching_ratio_map_map;

		reaction_info_base() :reaction_direction(both), reaction_rate(0.0) {}

	};

	//reaction information struct
	struct reaction_info :public reaction_info_base {
		//reaction_name_t reaction_name;
		index_int_t reaction_index;
		rate_constant_t rate_constant;
		prefactor_A_type_t prefactor_A;
		delta_t delta;
		barrier_E_t barrier_E;
		//reaction_rate_t reaction_rate;

		//three body reaction
		bool isThreeBodyReaction;
		spe_name_3body_enhancement_coef_t spe_name_3body_enhancement_coef;
		//pressure dependent reaction
		bool isPressureDependentReaction;
		low_pressure_limit_t low_pressure_limit;
		TROE_centering_t TROE_centering;

		std::vector<spe_index_weight_t> reactant;
		std::vector<spe_index_weight_t> product;

		//net reactant and net prodect
		std::vector<spe_index_weight_t> net_reactant;
		std::vector<spe_index_weight_t> net_product;


	};


	class relationshipParser {
	public:
		//Read the file named "chem.out", read in the element and species information.
		static void read_chem_out_ele_spe(std::vector<element_info>& element_v, std::vector<spe_info>& species_v, spe_name_index_map_t& spe_name_index_map, std::string file_in = "./input/chem.out");

		//s2f --> save to file
		static void spe_information_s2f(const std::vector<spe_info>& species_v, std::string file_out = "./input/species_labelling.csv");
		//save to json file
		static void spe_information_s2json(const std::vector<spe_info>& species_v, std::string file_out = "./input/species_information.json");

		//Read the file named "chem.out", read in the chemical reactions.
		static void read_chem_out_reaction(const std::vector<spe_info>& species_v, std::vector<reaction_info>& reaction_v, const spe_name_index_map_t& spe_name_index_map, std::string file_in = "./input/chem.out");
		
		static void reaction_information_s2f(const std::vector<spe_info> & species_v, const std::vector<reaction_info> & reaction_v,
			const reactionNetwork_chemkin_index_map_t& reactionNetwork_chemkin_index_map, std::string file_out = "./input/reaction_labelling.csv");
		//save to json file
		static void reaction_information_s2json(const std::vector<spe_info> & species_v, const std::vector<reaction_info> & reaction_v,
			const reactionNetwork_chemkin_index_map_t& reactionNetwork_chemkin_index_map, std::string file_out = "./input/reaction_information.json");

	public:
		/*
		* set species source or sink information, viz. which reaction is it in? participate as a reactant or product or both?
		* s->source
		* k->sink
		*/
		static void set_spe_sk_info(std::vector<spe_info>& species_v, const std::vector<reaction_info>& reaction_v);
		/*
		* set reaction net reactant and net product information
		*/
		static void set_reaction_net_reactant_product(reaction_info& reaction);
		static void set_reaction_net_reactant_product(std::vector<reaction_info>& reaction_v);

		/*
		* read reactionNetwork_chemkin_index_map from "chem.out"
		*/
		static void read_reactionNetwork_chemkin_index_map(reactionNetwork_chemkin_index_map_t& reactionNetwork_chemkin_index_map, std::string infile = "./input/chem.out");

		/*
		* read duplicated reaction from "chem.out"
		*/
		static void read_chem_out_duplicated_reaction(std::vector<int> &duplicate_reaction_ind, std::string str);

		/*
		* parse reaction, type of reaction_info, tell whether it is forward, back or both reaction
		* basically tell the reaction direction
		*/
		static void parse_reaction_direction(reaction_info& reaction);

	};



}/*namespace relationshipParser_sr*/

#endif
