#ifndef __RANDOM_CPP_
#define __RANDOM_CPP_

#include "random.h"
namespace random_sr {
	std::size_t random::return_index_randomly_given_probability_vector(const std::vector<double>& probability)
	{
		//initialize the discrete distribution generator for random integers
		boost::random::discrete_distribution<> dist(probability);
		//pick a vertex randomly, return the iterator

		//re-seed the generator
		this->generator = boost::mt19937(boost::random::random_device()() + this->random_seed_for_this_core);

		return dist(generator);
	}

	std::size_t random::return_0_or_1_evenly_randomly()
	{
		//initialize the discrete distribution generator for random integers
		boost::random::discrete_distribution<> dist(this->prob);
		//pick a vertex randomly, return the iterator
		//re-seed the generator
		this->generator = boost::mt19937(boost::random::random_device()() + this->random_seed_for_this_core);
		return dist(generator);
	}
}


#endif
