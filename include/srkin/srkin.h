/*
 * srkin.h
 *
 *  Created on: Jul 9, 2014
 *      Author: Shirong
 */

#ifndef __SRKIN_H_
#define __SRKIN_H_

#include <string>
 //#include <iostream>
 //#include <vector>
 //#include <algorithm>
 //#include "../tools/misc/global_extern_vars.h"
 //#include "../tools/misc/misc_template.h"
#include "../relationshipParser/relationshipParser.h"


//ns represents namespace
namespace srkin_sr {

	namespace rsp = relationshipParser_sr;

	class srkin {
	protected:
		//sr
		std::string cwd_sk;
		/*
		* https://en.wikipedia.org/wiki/Gas_constant
		* universal gas constant in cgs units
		* 8.3144621(75)E7 ergs/(mole*K)
		*/
		const double RU_sk;
		/*
		* universal gas constant used only in conjuction with activation energy
		* preferred units, 1.9872041(18) cal/(mol*K)
		*/
		const double RUC_sk;
		/*
		* pressure of one standard atmosphere
		* cgs units, 1.01325E6 dynes/cm**2
		*/
		const double PA_sk;

		//pressure
		rsp::pressure_t pressure_sk;
		//elements name and index map
		rsp::ele_name_index_map_t element_name_index_map_sk;
		//species name and index map
		rsp::spe_name_index_map_t spe_name_index_map_sk;

		//elements vector
		std::vector<rsp::element_info> element_v_sk;
		//species vector
		std::vector<rsp::spe_info> species_v_sk;
		//reaction vector
		std::vector<rsp::reaction_info> reaction_v_sk;

		//species concentration
		std::vector<double> species_conc_v_sk;
		//volume
		rsp::volume_t volume_sk;
		//temperature
		rsp::temperature_t temperature_sk;

	public:
		srkin(std::string cwd_sr_in, rsp::temperature_t t, rsp::pressure_t p, rsp::volume_t v);
		virtual ~srkin();

	public:
		void read_initial_spe_concentration(const std::size_t nkk, std::string filename = "/input/setting.json");

	public:
		std::string get_cwd();
		void set_cwd(std::string cwd_in);

		rsp::temperature_t get_temperature();
		void set_temperature(rsp::temperature_t t_in);

		rsp::volume_t get_volume();
		void set_volume(rsp::volume_t volume_in);

		rsp::pressure_t get_pressure();
		void set_pressure(rsp::pressure_t pressure_in);

		const std::vector<rsp::concentration_t>& get_species_conc_v() const;
		void set_species_conc_v(const std::size_t spe_ind, const double conc);
		void set_species_conc_v(const double speConcV[]);

	public:
		/*
		* calculate the reaction rate coefficient, or says reaction rate const
		* cakbe in Chemkin
		* notice activation energy is in unit of cal/mol
		*/
		double cal_reaction_rate_coef(std::size_t reaction_ind, rsp::temperature_t T, rsp::pressure_t P);

		/*
		* set the reaction rate coefficient, or says reaction rate const
		*/
		void set_reaction_rate_coef(std::size_t reaction_ind, rsp::temperature_t T, rsp::pressure_t P);

		//default-->continuous version, concentration can be coutinuous
		void cal_reaction_rate(const std::size_t reaction_ind);
		void cal_reaction_rate();
		//calculate the reaction rates given temperature and molar concentration
		void cal_reaction_rate(const double *T, const double *C, double *FWDK, double *REVK);

		//discrete version, plugin species numbers, check species number
		//for example, 2A-->A2, if [A]=1, this reaction should not happen
		void cal_reaction_rate_discrete(const std::size_t reaction_ind);
		void cal_reaction_rate_discrete();


		rsp::reaction_rate_t get_reaction_rate(const std::size_t reaction_ind);
		rsp::reaction_rate_t cal_spe_destruction_rate(const std::size_t spe_ind);

		//Returns the molar creation and destruction rates of the species given temperature(s) and molar concentration
		void cal_spe_destruction_rate(const double *T, const double *C, double *CDOT, double *DDOT);



	};

} /*namespace srkin_sr */



#endif /* __SRKIN_H_ */
