#ifndef __CLCF_PARSER_CPP_
#define __CLCF_PARSER_CPP_

#include "clcf_parser.h"


//clcf represents command line and config file

clcf_parser::clcf_parser() :generic("Generic options"), config("Configuration"), hidden("Hidden options")
{
	// Declare a group of options that will be 
	// allowed only on command line
	generic.add_options()
		("version,v", "print version string")
		("help", "produce help message")
		("config,c", po::value<std::string>()->default_value("config.cfg"),
			"name of a file of a configuration.")
		;

	// Declare a group of options that will be 
	// allowed both on command line and in
	// config file
	config.add_options()
		("optimization", po::value<int>()->default_value(10),
			"optimization level")
		("cwd",
			po::value<std::string>(),
			"current working directory")
		("include-path,I",
			po::value< std::vector<std::string> >()->composing(),
			"include path")
		("NumProcessor,P", po::value<int>()->default_value(1),
			"Num of Processor")
		("NumTraj,N", po::value<int>()->default_value(1),
			"Num of trajectory")
		("topN,T", po::value<int>()->default_value(1),
			"topN pathways")
		("pathwayEndWith,E",
			po::value<std::string>()->default_value(std::string("ALL")),
			"pathway ending with a species, specify which species it is")
		;

	// Hidden options, will be allowed both on command line and
	// in config file, but will not be shown to the user.
	hidden.add_options()
		("input-file", po::value< std::vector<std::string> >(), "input file")
		;
}

clcf_parser::~clcf_parser() { ; }

int clcf_parser::parse(int ac, char** av)
{
	try {//try[
		//combinations
		po::options_description cmdline_options;
		cmdline_options.add(generic).add(config).add(hidden);

		po::options_description config_file_options;
		config_file_options.add(config).add(hidden);

		po::options_description visible("Allowed options");
		visible.add(generic).add(config);

		po::positional_options_description p;
		p.add("input-file", -1);

		//store in vm
		store(po::command_line_parser(ac, av).
			options(cmdline_options).positional(p).run(), vm);
		notify(vm);

		std::ifstream ifs(vm["config"].as<std::string>().c_str());
		if (!ifs)
		{
			std::cout << "can not open config file: " << vm["config"].as<std::string>() << "\n";
			return 0;
		}
		else
		{
			store(parse_config_file(ifs, config_file_options), vm);
			notify(vm);
		}

		//print out	
		if (vm.count("help")) {
			std::cout << visible << "\n";
			return 0;
		}

	}//]try
	catch (std::exception& e)
	{
		std::cout << e.what() << "\n";
		return 1;
	}

	return EXIT_SUCCESS;
}

void clcf_parser::print() {
	//print out	
	if (vm.count("version")) {
		std::cout << "Multiple sources example, version 1.0\n";
		return;
	}

	if (vm.count("include-path"))
	{
		std::cout << "Include paths are: "
			<< vm["include-path"].as< std::vector<std::string> >() << "\n";
	}

	if (vm.count("input-file"))
	{
		std::cout << "Input files are: "
			<< vm["input-file"].as< std::vector<std::string> >() << "\n";
	}

	if (vm.count("optimization"))
	{
		std::cout << "Optimization level is "
			<< vm["optimization"].as< int >() << "\n";
	}

}


//cf represents and config file
cf_parser::cf_parser() :config("Configuration")
{
	// Declare a group of options that will be allowed in config file
	config.add_options()
		//time section
		("time.min_time", po::value<my_time_t>()->default_value(0.0),
			"minimum time")
		("time.max_time", po::value<my_time_t>()->default_value(1.0),
			"maximum time")
		("time.critical_time", po::value<my_time_t>()->default_value(1.0),
			"critical time")
		("time.sys_min_time", po::value<my_time_t>()->default_value(1.0e-10),
			"system minimum time")
		("time.tau", po::value<my_time_t>()->default_value(1.0),
			"system reference time tau")

		//chem_init section
		("chem_init.pressure_atm", po::value<double>()->default_value(5.54651),
			"initial pressure in unit of atm")
		("chem_init.init_temperature", po::value<double>()->default_value(1000.0),
			"initial temperature in unit of Kelvin")
		// we can use std::vector<std::string> instead of just std::string, if use std::string, use boost::tokenizer to parse data easily
		("chem_init.species_index_concentration", po::value<std::string>(),
			"species index and corresponding concentration")

		//lsode_init section
		("lsode_init.dt", po::value<double>()->default_value(5.54651),
			"lsode time increment in unit of second")
		("lsode_init.atol", po::value<double>()->default_value(5.54651),
			"absolute tolerance parameter")
		("lsode_init.rtol", po::value<double>()->default_value(5.54651),
			"relative tolerance parameter")
		("lsode_init.mf", po::value<int>()->default_value(1),
			"method parameter")
		("lsode_init.jt", po::value<int>()->default_value(1),
			"method parameter")
		("lsode_init.itask", po::value<int>()->default_value(1),
			"itask parameter")
		("lsode_init.iopt", po::value<int>()->default_value(0),
			"0 to indicate no optional inputs used")
		("lsode_init.itol", po::value<int>()->default_value(1),
			"1 or 2 according as atol (below) is a scalar or array")
		("lsode_init.istate", po::value<int>()->default_value(1),
			"an index used for input and output to specify the state of the calculation.")

		("lsode_init.deltaN1", po::value<int>()->default_value(1),
			"print out data every deltaNN1 points")
		("lsode_init.deltaN2", po::value<int>()->default_value(1),
			"print out data every deltaNN2 points")

		//temperature section
		("T.critical_temperature", po::value<double>()->default_value(888),
			"critical temperature of dlsode temperature propagator")
		("T.end_temperature", po::value<double>()->default_value(900),
			"end temperature of dlsode temperature propagator")
		("T.target_temperature", po::value<double>()->default_value(900),
			"target temperature of pathway generation scheme")

		//pathway
		("pathway.init_spe", po::value<std::size_t>()->default_value(2),
			"initial species index of pathway generation scheme")
		("pathway.atom_followed", po::value<std::string>()->default_value("H"),
			"atom which will be followed")
		("pathway.dead_species", po::value<std::string>(),
			"dead species, species with no out edge or out rate is zero")
		("pathway.trapped_species", po::value<std::string>(),
			"trapped species, species are trapped in some local reaction with fast conversion rate")
		("pathway.fast_reaction", po::value<std::string>(),
			"fast reaction, local reaction with fast inter-conversion rate")

		//sohr section
		//number of time point
		("SOHR_init.timeN1", po::value<int>()->default_value(10),
			"Number of time points before critical time")
		("SOHR_init.timeN2", po::value<int>()->default_value(10),
			"Number of time points after critical time")
		("SOHR_init.N2C", po::value<double>()->default_value(10.0),
			"number of species to concentration conversion factor")
		//number of iterations
		("SOHR_init.iterationN1", po::value<int>()->default_value(10),
			"Number of inner iterations")
		("SOHR_init.iterationN2", po::value<int>()->default_value(10),
			"Number of outer iterations")
		//others
		//dir section
		("dir.include-path,I",
			po::value< std::vector<std::string> >()->composing(),
			"include path")
		;
}

int cf_parser::parse(std::string filename)
{
	try {//try[
		po::options_description config_file_options;
		config_file_options.add(config);

		std::ifstream ifs(filename.c_str());
		if (!ifs)
		{
			std::cout << "can not open config file: " << filename << "\n";
			return 0;
		}
		else
		{
			store(parse_config_file(ifs, config_file_options), vm);
			notify(vm);
		}

	}//]try
	catch (std::exception& e)
	{
		std::cout << e.what() << "\n";
		return 1;
	}

	return EXIT_SUCCESS;
}

void cf_parser::print() {
	//print out	
	if (vm.count("dir.include-path"))
	{
		std::cout << "Include paths are: "
			<< vm["dir.include-path"].as< std::vector<std::string> >() << "\n";
	}

	if (vm.count("time.min_time"))
	{
		std::cout << "The minimum time is "
			<< vm["time.min_time"].as< my_time_t >() << "\n";
	}

}


#endif
