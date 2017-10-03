#ifndef __RANDOM_H_
#define __RANDOM_H_

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp> //for boost::random::uniform_real_distribution

#include <boost/random/random_device.hpp> //for boost::random::random_device

namespace random_sr {

	class random {
	private:
		//pseudo-random number generators
		//http://www.boost.org/doc/libs/1_55_0/doc/html/boost_random/reference.html
		boost::mt19937 generator;
	private:
		//random seed for this core
		boost::uint32_t random_seed_for_this_core;

	private:
		//probability, for evenly distribution
		std::vector<double> prob;
	public:
		random(boost::uint32_t random_seed_in) :random_seed_for_this_core(random_seed_in) {
			this->generator = boost::mt19937(boost::random::random_device()() + this->random_seed_for_this_core);
			prob = std::vector<double>(2, 0.5);
		}
	public:
		/*
		* http://stackoverflow.com/questions/25474118/getting-seeds-value-from-boost-random-mt19937-generator
		* Rather than trying to extract what you think of as the seed from mt19937, it's easier set the seed explicitly in both runs for reproducibility.
		* See boost's random_demo.cpp, about 20 lines into the main, for an example of setting the seed.
		* The comment points out that using std::time(0) can inadvertently lead to correlated results from two
		* generators if they are seeded in rapid succession based on time.
		*/
		/*
		* Change seed to something else.
		* http://www.boost.org/doc/libs/1_56_0/libs/random/example/random_demo.cpp
		*
		* Caveat: std::time(0) is not a very good truly-random seed.  When
		* called in rapid succession, it could return the same values, and
		* thus the same random number sequences could ensue.  If not the same
		* values are returned, the values differ only slightly in the
		* lowest bits.  A linear congruential generator with a small factor
		* wrapped in a uniform_smallint (see experiment) will produce the same
		* values for the first few iterations.   This is because uniform_smallint
		* takes only the highest bits of the generator, and the generator itself
		* needs a few iterations to spread the initial entropy from the lowest bits
		* to the whole state.
		*/
		//When calling other functions which take a generator or distribution
		//as a parameter, make sure to always call by reference (or pointer).
		//Calling by value invokes the copy constructor, which means that the
		//sequence of random numbers at the caller is disconnected from the
		//sequence at the callee.

		double random01()
		{
			/*
			* re-seed the generator
			*std::time(0) doesn't work here
			*/
			//generator.seed(static_cast<unsigned int>(std::time(0)));
			this->generator = boost::mt19937(boost::random::random_device()() + this->random_seed_for_this_core);
			boost::uniform_01<boost::mt19937> dist(generator);
			return dist();
		}

		/*
		* 	I.e., if all of the values are the same. This is simply due to repeated re-initialization of the random number generator from the system clock.
		* 	http://www.bnikolic.co.uk/blog/cpp-boost-rand-normal.html
		* 	every time the function random_min_max is called it makes a copy of the generator object
		* 	and advances the state of this copy only, not of the parent object.
		* 	Solution by passing references
		* 	The simplest solution in this case it of course to ensure a single generator object exists and only pass references to it.
		* 	To be sure that you are getting good quality random numbers it is best to advance a single generator for the duration of the numerical computation task.
		* 	But boost random number generators do not hold a global state and it is up to the user to ensure generator states are preserved.
		* 	The easiest and usually sufficient way of doing this is to pass references to a single, top-level, generator object.
		* 	Don't use static inside class with random number generator
		*/

		double random_min_max(double min_t, double max_t)
		{
			//initialize the distribution
			boost::random::uniform_real_distribution<double> dist(min_t, max_t);

			//re-seed the generator
			this->generator = boost::mt19937(boost::random::random_device()() + this->random_seed_for_this_core);

			return dist(generator);
		}
		/*
		* return the index given a vector of probability
		*/
		std::size_t return_index_randomly_given_probability_vector(const std::vector<double>& prob);
		std::size_t return_0_or_1_evenly_randomly();


	};


}



#endif
