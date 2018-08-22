/* mlev-hiv-model - A program for running multi-level HIV-1 simulations
 * Copyright (C) 2017 Christiaan H. van Dorp <chvandorp@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "cohort.hpp"

Cohort::Cohort(std::string identifier) : identifier(identifier) {
	std::cout << "created Cohort (identifier: " << identifier << ")" << std::endl;
	// TODO: other options?
}

Cohort::~Cohort() {
	workerpool.waitForWorkerPoolToExit();
	for ( auto hit = population.begin(); hit != population.end(); ++hit ) {
		Host* host = (*hit);
		delete host;
	}
	delete wildtype;
}

void Cohort::init(unsigned long seed, int cohort_size) {
	/* TODO: options for single virus,
   * virus with a number of random mutations,
   * randomly scrambled + evolved virus
   */
	// seed rng:
	rng.seed(seed);
	// init fitnessfunction:
	ffun.init(rng);
	// make a virus:
	wildtype = new Virus(&ffun, rng, Virus::EVOLVE);
	// make immune responses
	immuneResponses = createImmuneResponses(rng);
	population.resize(cohort_size, nullptr);
	for ( auto hit = population.begin(); hit != population.end(); ++hit ) {
		Virus* mutant = nullptr;
		// make sure that an integer is passed to the constructor
		if ( MULTIPLE_PLEV_FOUNDERS ) {
			Virus virus(&ffun, rng, Virus::EVOLVE);
			mutant = new Virus(virus, rng, int(NUMBER_INITIAL_MUTATIONS));
		} else {
			mutant = new Virus(*wildtype, rng, int(NUMBER_INITIAL_MUTATIONS));
		}
		IndPars ipars(rng);
		double R0 = iLevBasicReproductionNumber(*mutant, ipars);
		if ( R0 > 1.0 ) {
			*hit = new Host(rng, immuneResponses, *mutant, ipars, 0.0);
		}
		else {
			std::cout << "WARNING: failed to infect host: R0 = " << R0 << std::endl;
		}
		delete mutant;
	}
	workerpool.initWorkerPool(NUMBER_OF_THREADS, rng);
}

void Cohort::run() {
	std::cout << "running i-level simulations..." << std::endl;
	for ( auto hit = population.begin(); hit != population.end(); ++hit ) {
		workerpool.addNewJob(*hit);
	}
	workerpool.syncWorkerThreads(true); // boolean indicates that ETA is printed
	std::cout << " done!" << std::endl;
}

void Cohort::print(std::ostream & os) const {
	os << "<cohort identifier='" << identifier << "' >" << std::endl;
	os << xmlStringParameters() << std::endl;
	os << xmlStringMhcLoci(immuneResponses, identifier) << std::endl;
	os << wildtype << std::endl;
	for ( auto hit = population.begin(); hit != population.end(); ++hit ) {
		Host* host = *hit;
		if ( host != nullptr ) {
			os << *host << std::endl;
		}
	}
	os << "</cohort>" << std::endl;
}

std::ostream & operator<<(std::ostream & os, const Cohort & cohort) {
	cohort.print(os);
	return os;
}

// other functions

void highlyDetailedInfectionHistory(unsigned long seed, std::string identifier) {
	Rng rng(seed);
	FitnessFunction ffun; ffun.init(rng);
  std::ofstream offun;
  std::string funfilename = DATA_FOLDER + "fitnessfunction-file-" + identifier + ".xml";
  offun.open(funfilename.c_str());
  if ( offun ) {
    ffun.print(offun);
  } else {
		throw std::logic_error("could not write to file " + funfilename + RIGHT_HERE);
  }

  Virus wildtype(&ffun, rng, Virus::EVOLVE);
	Virus virus(wildtype, rng, int(NUMBER_INITIAL_MUTATIONS));

	std::vector<MhcLocus> irs = createImmuneResponses(rng);
	//IndPars ipars(rng);
	IndPars ipars; // TESTING: use default parameters
  std::cout << "Malthusian fitness = " << iLevBasicMalthusianFitness(virus, ipars) << std::endl;

  if ( iLevBasicReproductionNumber(virus, ipars) <= 1.0 ) {
		// TODO: keep trying, say 10 times, only then throw exception
		throw std::logic_error("cannot infect host: i-level fitness is too small" + RIGHT_HERE);
	}

	// create a Host object
	Host host(rng, irs, virus, ipars, 0.0);

	StopWatch stopwatch; // how fast?
	host.run(rng, true);
	std::cout << "time needed: " << stopwatch << std::endl;

	std::ofstream of;
	std::string filename = DATA_FOLDER + "hd-infection-history-" + identifier + ".xml";
	of.open(filename.c_str());
	if ( of ) {
		of << "<hd_infection_history identifier='" << identifier << "' >" << std::endl;
		of << xmlStringParameters() << std::endl;
		// of << xmlStringMhcLoci(irs, version) << std::endl; // don't need all responses...
		of << host << std::endl;
		of << "</hd_infection_history>";
	} else {
		throw std::logic_error("could not write to file " + filename + RIGHT_HERE);
	}
}
