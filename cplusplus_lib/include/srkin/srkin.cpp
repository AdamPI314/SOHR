/*
 * srkin.cpp
 *
 *  Created on: Jul 9, 2014
 *      Author: Shirong
 */

#ifndef __SRKIN_CPP_
#define __SRKIN_CPP_

#include "srkin.h"
#include <boost/property_tree/ptree.hpp> //for property_tree
#include <boost/property_tree/json_parser.hpp> //for json_reader
#include <boost/foreach.hpp> //for BOOST_FOREACH
#include <boost/optional/optional.hpp> //for optional

namespace srkin_sr {

	srkin::srkin(std::string cwd_sr_in, rsp::temperature_t t, rsp::pressure_t p, rsp::volume_t v) :RU_sk(8.3144621E7), RUC_sk(1.9872041), PA_sk(1.01325E6)
	{
		set_cwd(cwd_sr_in);
		set_temperature(t);
		set_pressure(p);
		set_volume(v);


		rsp::relationshipParser::read_chem_out_ele_spe(this->element_v_sk, this->species_v_sk, this->spe_name_index_map_sk, this->cwd_sk + "/input/chem.out");
		rsp::relationshipParser::read_chem_out_reaction(this->species_v_sk, this->reaction_v_sk, this->spe_name_index_map_sk, this->cwd_sk + "/input/chem.out");

		rsp::relationshipParser::set_reaction_net_reactant_product(this->reaction_v_sk);
		rsp::relationshipParser::set_spe_sk_info(this->species_v_sk, this->reaction_v_sk);

		read_initial_spe_concentration(this->species_v_sk.size(), "/input/setting.json");

		//set prefactor_A
		for (std::size_t i = 0; i < this->reaction_v_sk.size(); ++i) {
			this->reaction_v_sk[i].rate_constant = cal_reaction_rate_coef(i, this->temperature_sk, this->pressure_sk);
		}

	}

	srkin::~srkin()
	{
	}

	void srkin::read_initial_spe_concentration(const std::size_t nkk, std::string filename)
	{
		//configurations, read configuration file named "setting.json"
		boost::property_tree::ptree pt;
		//read configuration file "setting.cfg"
		boost::property_tree::read_json(this->cwd_sk + std::string("/input/setting.json"), pt, std::locale());

		this->species_conc_v_sk.resize(nkk);

		//read with json_parser as property_tree
		//nice and easy
		BOOST_FOREACH(boost::property_tree::ptree::value_type const& key1, pt.get_child("chem_init.species_index_concentration")) {
			species_conc_v_sk[boost::lexical_cast<std::size_t>(key1.first)] = key1.second.get_value<double>();
		}

	}

	std::string srkin::get_cwd()
	{
		return this->cwd_sk;
	}

	void srkin::set_cwd(std::string cwd_in)
	{
		this->cwd_sk = cwd_in;
	}

	rsp::temperature_t srkin::get_temperature()
	{
		return this->temperature_sk;
	}

	void srkin::set_temperature(rsp::temperature_t t_in)
	{
		this->temperature_sk = t_in;
	}

	rsp::volume_t srkin::get_volume()
	{
		return this->volume_sk;
	}

	void srkin::set_volume(rsp::volume_t volume_in)
	{
		this->volume_sk = volume_in;
	}

	rsp::pressure_t srkin::get_pressure()
	{
		return this->pressure_sk;
	}

	void srkin::set_pressure(rsp::pressure_t pressure_in)
	{
		this->pressure_sk = pressure_in;
	}

	const std::vector<rsp::concentration_t>& srkin::get_species_conc_v() const
	{
		// TODO: insert return statement here
		return this->species_conc_v_sk;
	}

	void srkin::set_species_conc_v(const std::size_t spe_ind, const double conc)
	{
		this->species_conc_v_sk[spe_ind] = conc;
	}

	void srkin::set_species_conc_v(const double speConcV[])
	{
		for (std::size_t i = 0; i < this->species_conc_v_sk.size(); ++i) {
			this->species_conc_v_sk[i] = speConcV[i];
		}
	}

	double srkin::cal_reaction_rate_coef(std::size_t reaction_ind, rsp::temperature_t T, rsp::pressure_t P)
	{
		if (T == 0)
			return this->reaction_v_sk[reaction_ind].prefactor_A;
		else if (this->reaction_v_sk[reaction_ind].delta == 0)
			return this->reaction_v_sk[reaction_ind].prefactor_A;
		else
			return this->reaction_v_sk[reaction_ind].prefactor_A* pow(T, this->reaction_v_sk[reaction_ind].delta) * exp(-this->reaction_v_sk[reaction_ind].barrier_E / T / this->RUC_sk);
	}

	void srkin::set_reaction_rate_coef(std::size_t reaction_ind, rsp::temperature_t T, rsp::pressure_t P)
	{
		this->reaction_v_sk[reaction_ind].rate_constant = this->cal_reaction_rate_coef(reaction_ind, T, P);
	}

	void srkin::cal_reaction_rate(const std::size_t reaction_ind)
	{

		this->reaction_v_sk[reaction_ind].reaction_rate = this->reaction_v_sk[reaction_ind].rate_constant;
		for (std::size_t i = 0; i < reaction_v_sk[reaction_ind].reactant.size(); ++i) {
			this->reaction_v_sk[reaction_ind].reaction_rate *= pow(this->species_conc_v_sk[reaction_v_sk[reaction_ind].reactant[i].first], reaction_v_sk[reaction_ind].reactant[i].second);
		}
	}

	void srkin::cal_reaction_rate()
	{
		for (std::size_t i = 0; i < this->reaction_v_sk.size(); ++i)
			cal_reaction_rate(i);
	}

	void srkin::cal_reaction_rate(const double * T, const double * C, double * FWDK, double * REVK)
	{
		set_temperature(*T);
		set_species_conc_v(C);
		for (std::size_t i = 0; i < this->reaction_v_sk.size(); ++i) {
			cal_reaction_rate(i);
			FWDK[i] = this->reaction_v_sk[i].reaction_rate;
			REVK[i] = 0.0;
		}
	}

	void srkin_sr::srkin::cal_reaction_rate_discrete(const std::size_t reaction_ind)
	{

		this->reaction_v_sk[reaction_ind].reaction_rate = this->reaction_v_sk[reaction_ind].rate_constant;
		for (std::size_t i = 0; i < reaction_v_sk[reaction_ind].reactant.size(); ++i) {
			//check, number of species must be not less than the stoichoimetric coefficient
			if (this->species_conc_v_sk[reaction_v_sk[reaction_ind].reactant[i].first] < reaction_v_sk[reaction_ind].reactant[i].second) {
				this->reaction_v_sk[reaction_ind].reaction_rate = 0.0;
				break;
			}
			else
				this->reaction_v_sk[reaction_ind].reaction_rate *= pow(this->species_conc_v_sk[reaction_v_sk[reaction_ind].reactant[i].first], reaction_v_sk[reaction_ind].reactant[i].second);
		}
	}

	void srkin_sr::srkin::cal_reaction_rate_discrete()
	{
		for (std::size_t i = 0; i < this->reaction_v_sk.size(); ++i)
			cal_reaction_rate_discrete(i);
	}

	rsp::reaction_rate_t srkin::get_reaction_rate(const std::size_t reaction_ind)
	{
		return this->reaction_v_sk[reaction_ind].reaction_rate;
	}

	rsp::reaction_rate_t srkin::cal_spe_destruction_rate(const std::size_t spe_ind)
	{
		double spe_dr_tp = 0.0;
		for (std::size_t i = 0; i < species_v_sk[spe_ind].reaction_k_index_s_coef_v.size(); ++i) {
			spe_dr_tp += reaction_v_sk[species_v_sk[spe_ind].reaction_k_index_s_coef_v[i].first].reaction_rate *species_v_sk[spe_ind].reaction_k_index_s_coef_v[i].second;
		}
		return spe_dr_tp;
	}

	void srkin::cal_spe_destruction_rate(const double * T, const double * C, double * CDOT, double * DDOT)
	{
		set_temperature(*T);
		set_species_conc_v(C);
		cal_reaction_rate();

		for (std::size_t i = 0; i < this->species_v_sk.size(); ++i) {
			CDOT[i] = 0.0;
			DDOT[i] = cal_spe_destruction_rate(i);
		}
	}

} /*namespace srkin_sr */

#endif /* __SRKIN_CPP_ */



