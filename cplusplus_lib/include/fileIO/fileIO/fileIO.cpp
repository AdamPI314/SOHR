#ifndef __FILEIO_CPP_
#define __FILEIO_CPP_

#include "fileIO.h"
#include "../CSVFileReader/CSVReader.h"

bool fileIO::fileIO::read_chem_out_duplicated_reaction(std::vector<std::size_t>& duplicate_reaction_ind, std::string str)
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
	return true;
}


//Read the file named "chem.out", read in the corresponding reaction index
bool fileIO::fileIO::read_chem_out_index(reaction_network_chemkin_index_t & reaction_index_network, std::string str)
{
	std::ifstream fin(str.c_str());
	std::string str_t;

	//reaction pattern without index
	const char* pattern1 = "[\\w+\\(\\)\\-,]+(?:=>|<=|=|<=>)+[\\w+\\(\\)\\-,]+";
	boost::regex regexPattern1(pattern1);
	const char* pattern1_index = "\\d+\\.\\s+[\\w+\\(\\)\\-,]+(?:=>|<=|=|<=>)+[\\w+\\(\\)\\-,]+";
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
		std::vector<std::size_t> vec_index;
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
			reaction_index_network[edge_counter] = vec_index; //forward reaction index
			for (size_t k = 0; k < vec_index.size(); ++k) {
				vec_index[k] *= (-1);
			}
			reaction_index_network[edge_counter + 1] = vec_index; //backward reaction index

			edge_counter += 2; //two reactions, forward and backward
		}//if

		else if ((*itr2_arrow) == std::string("=>")) {
			//std::cout<<*itr2_arrow<<"\t"<<edge_counter<<std::endl;
			reaction_index_network[edge_counter] = vec_index; //forward reaction index
			++edge_counter; //only forward reaction
		}//else if

		else if ((*itr2_arrow) == std::string("<=")) {
			//std::cout<<*itr2_arrow<<"\t"<<edge_counter<<std::endl;
			for (size_t k = 0; k < vec_index.size(); ++k) {
				vec_index[k] *= -1;
			}
			reaction_index_network[edge_counter + 1] = vec_index; //backward reaction index
			++edge_counter; //only backward reaction
		}//else if


	}//for

	fin.close(); fin.clear();
	return true;
}


//Read in the uncertainty, use the nominal value of rate const, set uncertainty 1.0
bool fileIO::fileIO::read_generate_uncertainties_w2f_nominal(std::vector<double>& uncertainties, std::string str_in, std::string str_out, std::string chem_out)
{
	//read uncertainty
	std::ifstream fin(str_in.c_str());
	int index; double value;
	while ((!fin.eof()) && (fin >> index)) {
		fin >> value;
		uncertainties.push_back(value);
		//cout<<index<<" "<<value<<endl;
	}
	fin.clear(); fin.close();

	//generate uncertainty
	//boost::mt19937 generator(static_cast<boost::uint32_t>(std::time(0)));
	std::vector<double> rate_Constant_variance(uncertainties.size());

	for (size_t i = 0; i < uncertainties.size(); ++i) {
		//boost::random::uniform_real_distribution<> dist(1.0/uncertainties[i], uncertainties[i]);
		rate_Constant_variance[i] = 1.0;
	}
	std::copy(rate_Constant_variance.begin(), rate_Constant_variance.end(), uncertainties.begin());

	//read in duplicated reactions
	std::vector<std::size_t> duplicate_reaction_ind;
	//	chemkin_cpp::chemkin::read_chem_out_duplicated_rxn(duplicate_reaction_ind, "./input/chem.out");
	fileIO::read_chem_out_duplicated_reaction(duplicate_reaction_ind, chem_out);
	//for duplicated rxn, set their variance the same value
	for (size_t i = 0; i < uncertainties.size(); ++i) {
		if (std::find(duplicate_reaction_ind.begin(), duplicate_reaction_ind.end(), i) != duplicate_reaction_ind.end())
			uncertainties[i] = uncertainties[i - 1];
	}

	std::ofstream fout_uncertainties(str_out.c_str());
	for (size_t i = 0; i < uncertainties.size(); ++i) {
		fout_uncertainties << i + 1 << "," << uncertainties[i] << std::endl;
	}
	fout_uncertainties.clear(); fout_uncertainties.close();

	return true;
}

//Read in the uncertainty, randomly generate rate const according to uncertainty
bool fileIO::fileIO::read_generate_uncertainties_w2f_random(std::vector<double>& uncertainties, std::string str_in, std::string str_out, std::string chem_out, boost::uint32_t random_seed_for_this_core)
{
	uncertainties.resize(0);

	//read uncertainty
	std::ifstream fin(str_in.c_str());
	int index; double value;
	while ((!fin.eof()) && (fin >> index)) {
		fin >> value;
		uncertainties.push_back(value);
		//cout<<index<<" "<<value<<endl;
	}
	fin.clear(); fin.close();

	//generate uncertainty
	//boost::mt19937 generator(static_cast<boost::uint32_t>(std::time(0)));
	boost::random_device rd;
	boost::mt19937 generator(rd() + random_seed_for_this_core);
	//boost::mt19937 generator(boost::random::random_device()() + random_seed_for_this_core);
	std::vector<double> rate_Constant_variance(uncertainties.size());

	for (size_t i = 0; i < uncertainties.size(); ++i) {
		boost::random::uniform_real_distribution<> dist(1.0 / uncertainties[i], uncertainties[i]);
		rate_Constant_variance[i] = dist(generator);
	}
	std::copy(rate_Constant_variance.begin(), rate_Constant_variance.end(), uncertainties.begin());

	//read in duplicated reactions
	std::vector<std::size_t> duplicate_reaction_ind;
	fileIO::read_chem_out_duplicated_reaction(duplicate_reaction_ind, chem_out);
	//for duplicated rxn, set their variance the same value
	for (size_t i = 0; i < uncertainties.size(); ++i) {
		if (std::find(duplicate_reaction_ind.begin(), duplicate_reaction_ind.end(), i) != duplicate_reaction_ind.end())
			uncertainties[i] = uncertainties[i - 1];
	}

	std::ofstream fout_uncertainties(str_out.c_str());
	for (size_t i = 0; i < uncertainties.size(); ++i) {
		fout_uncertainties << i + 1 << "," << uncertainties[i] << std::endl;
	}
	fout_uncertainties.clear(); fout_uncertainties.close();

	return true;
}


//Read in the uncertainty from file
bool fileIO::fileIO::read_uncertainties_from_file(std::vector<double>& uncertainties, std::string str_in, std::string str_out, std::string chem_out)
{
	//read uncertainty
	std::ifstream fin(str_in.c_str());
	int index; double value;
	while ((!fin.eof()) && (fin >> index)) {
		fin >> value;
		uncertainties.push_back(value);
		//cout<<index<<" "<<value<<endl;
	}
	fin.clear(); fin.close();

	//read in duplicated reactions
	std::vector<std::size_t> duplicate_reaction_ind;
	fileIO::read_chem_out_duplicated_reaction(duplicate_reaction_ind, chem_out);
	//for duplicated rxn, set their variance the same value
	for (size_t i = 0; i < uncertainties.size(); ++i) {
		if (std::find(duplicate_reaction_ind.begin(), duplicate_reaction_ind.end(), i) != duplicate_reaction_ind.end())
			uncertainties[i] = uncertainties[i - 1];
	}

	std::ofstream fout_uncertainties(str_out.c_str());
	for (size_t i = 0; i < uncertainties.size(); ++i) {
		fout_uncertainties << i + 1 << "," << uncertainties[i] << std::endl;
	}
	fout_uncertainties.clear(); fout_uncertainties.close();

	return true;
}

std::vector<std::vector<double>> fileIO::fileIO::read_topN_line_csv_matrix(std::string filename, std::size_t topN)
{
	std::vector<std::vector<double> > mat;

	//initialize the time points we are going to calculate
	//read file pathway_time.csv
	std::ifstream fin(filename.c_str());
	size_t icounter = 0;
	std::vector<double> vec;
	for (CSVIterator file_itr(fin); (icounter++ < topN) && (file_itr != CSVIterator()); ++file_itr) {
		vec.clear();
		for (size_t i = 0; i < (*file_itr).size(); ++i) {
			vec.push_back(boost::lexical_cast<double>((*file_itr)[i]));
		}
		mat.push_back(vec);
	}
	return mat;
}

#endif

