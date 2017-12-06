#ifndef __FILEIO_H_
#define __FILEIO_H_

#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <boost/regex.hpp> //for boost::regex
#include <boost/lexical_cast.hpp>
#include <boost/random/mersenne_twister.hpp> //random seed
#include <boost/random/uniform_real_distribution.hpp> //random number
#include <boost/random/discrete_distribution.hpp> //random number
#include <boost/random/random_device.hpp> //for boost::random::random_device


#include "../../tools/misc/misc_template.h"

namespace fileIO {/*namespace fileIO*/

	namespace mt = misc_template;
	/*
	* reaction network reaction index and ChemKin reaction index lookup table, the reason the ChemKin index is a vector is there might exist duplicated reaction
	* ChemKin has its own index of reaction, our reaction network redefine reaction index, here is the corresponding relation between them
	* The former is C/C++ style index, The later is Fortran Style index, actually they are exact index in ChemKin space
	* negative index represents backward reaction
	*/
	typedef std::map<std::size_t, std::vector<std::size_t> > reaction_network_chemkin_index_t;

	class fileIO {
	public:
		//Record the duplicated reaction index
		static bool read_chem_out_duplicated_reaction(std::vector<std::size_t> &duplicate_reaction_ind, std::string str);
		//Read the file named "chem.out", read in the corresponding reaction index
		static bool read_chem_out_index(reaction_network_chemkin_index_t &reaction_index_network, std::string str = "./input/chem.out");

	public:
		//Read in the uncertainty, use the nominal value of rate const, set uncertainty 1.0
		static bool read_generate_uncertainties_w2f_nominal(std::vector<double> &uncertainties, std::string str_in, std::string str_out = "./output/uncertainties_random.csv", std::string chem_out = "./input/chem.out");
		//Read in the uncertainty, randomly generate rate const according to uncertainty
		static bool read_generate_uncertainties_w2f_random(std::vector<double> &uncertainties, std::string str_in, std::string str_out = "./output/uncertainties_random.csv", std::string chem_out = "./input/chem.out", boost::uint32_t random_seed_for_this_core=0);
		//Read in the uncertainty from file
		static bool read_uncertainties_from_file(std::vector<double> &uncertainties, std::string str_in, std::string str_out = "./output/uncertainties_random.csv", std::string chem_out = "./input/chem.out");

	public:
		//read top N line of csv matrix
		static std::vector<std::vector<double> > read_topN_line_csv_matrix(std::string filename= "./output/pathway_time_candidate.csv", std::size_t topN=1);

	};



}/*namespace fileIO*/


#endif
