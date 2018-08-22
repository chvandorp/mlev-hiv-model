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

#ifndef WORLD_HPP_
#define WORLD_HPP_

#include <vector>
#include <list> // for percentiles, ...
#include <iterator> // std::advance (for sampling from a std::list)
#include <cmath> // for computing statistics
#include <iostream>
#include <iomanip> // for printing real-time progress
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm> // std::find function
#include <stdexcept>

#include "macros.hpp"
#include "rng.hpp"
#include "network.hpp"
#include "host.hpp"
#include "virus.hpp"
#include "immune_response.hpp"
#include "calc_stats.hpp"
#include "logistic_model.hpp"
#include "parallelism.hpp"
#include "hashed_definitions.hpp"
#include "parameters.hpp"
#include "eta.hpp"

class World : public Network {
public:
	World(std::string , unsigned long , int , int , double );
	World(const World & , unsigned long , double ); // copy some properties of the old world
	~World();
  void run(double , bool ); // tmax, keep_alive
		// statistics based on a HostTrait...
	void printStatistics(); // not const, write to streams, use rng for bootstrapping
	template<class Numeric>
	double getMeanHostTrait(HostTrait<Numeric> ) const;
	template<class Numeric>
	double getSdHostTrait(HostTrait<Numeric> ) const;
	template<class NumericL, class NumericR>
	double getCorrelationHostTrait(HostTrait<NumericL> , HostTrait<NumericR> ) const;
	template<class Numeric>
	double getHeritabilityHostTrait(HostTrait<Numeric> ) const;
	/* TODO/FIXME: remove the xmlNode functions: simply return the BasicStats object
	 * (that can be printed as an xml node)
	 */
	template<class Numeric>
	std::string xmlNodeHostTrait(HostTrait<Numeric> , std::string ) const;
	// FIXME: ArgType1 and ArgType2 should be identical! (but e.g. T and const T & are in conflict)
	template<class Numeric, class ArgType1, class ArgType2>
	std::string xmlNodeHostArgTrait(HostArgTrait<Numeric, ArgType1> , ArgType2 , std::string ) const;
	template<class Numeric, class Error>
	std::string xmlNodeHostMaybeTrait(HostMaybeTrait<Numeric, Error> , std::string ) const;
	// FIXME: ArgType1 and ArgType2 should be identical! (but e.g. T and const T & are in conflict)
	template<class Numeric, class Error, class ArgType1, class ArgType2>
	std::string xmlNodeHostMaybeArgTrait(HostMaybeArgTrait<Numeric, Error, ArgType1> , ArgType2 , std::string ) const;
	// xml nodes correlation and heritability
	template<class NumericL, class NumericR>
	std::string xmlNodeCorrelationHostTrait(HostTrait<NumericL> , HostTrait<NumericR> , std::string ) const;
	template<class NumericL, class NumericR, class ArgTypeR1, class ArgTypeR2>
	// correlation between left trait and right trait, where the right trait takes an argument
	std::string xmlNodeCorrelationHostRArgTrait(HostTrait<NumericL> ,
			HostArgTrait<NumericR, ArgTypeR1> , ArgTypeR2 , std::string ) const;
	template<class Numeric>
	std::string xmlNodeHeritabilityHostTrait(HostTrait<Numeric> , std::string ) const;
	// other statistics
	int getAncestorTreeSize() const;
	void getPlevEpitopesForAllele(const MhcAllele & , BasicStats & , BasicStats & ) const;
	std::pair<double, bool> getGeneticDiversity() const; // mean Hamming distance (scaled)
	std::pair<double, bool> getMeanDistanceFromWt() const; // mean HD from wildtype (scaled)
	// auxiliary function for genetic diversity per locus
	std::vector<double> getAlleleFreqPerLocus(bool ) const; // pass allele
	std::vector<double> getGeneticDiversityPerLocus() const; // entropy
	/* getConsensus and sampleVirus return wildtype when there are no living hosts,
	 * and false as second return value
	 */
	std::pair<Virus, bool> getConsensus() const; // create a consensus sequence
	std::pair<Virus, bool> sampleVirus(Rng & ) const; // returns a virus sampled from the population (R indiv is used for weighing)
	// printing statistics
	std::string xmlNodeFitnessCosts() const;
	std::string xmlNodeGeneticDiversity() const;
	std::string xmlNodeAidsHazard();
	void printXmlGeneralStats(std::ostream & ) const; // most general statistics
	void printXmlHeritDataNode(std::ostream & ) const;
protected:
	// internal methods
	void updateState(double ); // dt
	/* after a burn-in phase to get enough contacts,
	 * an infection can be introduced
	 */
	void infectRandomAgents(int );
	// a worker pool for executing Jobs in parallel
	WorkerPool workerpool;
	// living hosts (for faster statistics) TODO: add agent (map)?
	std::unordered_set<Host*> livingHosts;
	/* there can be multiple wildtypes, calculate distance and
	 * pre-adaptation w.r.t correct lineage. Virus-es inherit an identifier
	 * from their parent
	 */
	std::map<Id, Virus> wildtypes;
	long int cumulativeIncidence;
	long int failedInfections;
	std::list<double> recentInfectionTimes; // used for incidence
	void updateRecentInfectionTimes(int );
	// keep track of all significant MHC associations
	std::map<AlleleLocusPair, std::list<Association>> mhcPolyAssociationsFET;
	std::map<AlleleLocusPair, Association> findAssociationsMhcPolyFET();
	std::map<AlleleLocusPair, std::list<Association>> mhcPolyAssociationsLRT;
	std::map<AlleleLocusPair, std::list<Association>> mhcPolyAssociationsLRTcorrected;
	std::map<AlleleLocusPair, Association> findAssociationsMhcPolyLRT(bool corrected);
	std::map<MhcID, std::list<Association>> mhcVlAssociations;
	std::map<MhcID, Association> findAssociationsMhcVl();
	// time in years
	double t; // FIXME: use a longer symbol to avoid accidental masking
	double dt_cache;
	double birth_rate; // TODO: move to Network.
	FitnessFunction fitnessFunction;
	/* only a part of all loci can be an epitope.
	 * The others are there for 'compensatory' effects.
	 * Hosts choose alleles (sets of immune responses) from the alleles
	 * in allMhcLoci
	 */
	std::vector<MhcLocus> allMhcLoci;
	mutable Rng rng; // FIXME: remove the mutable, pass external RNG for bootstrap statistics
	Id agentIdCounter; // to ID agents // TODO: move to Network
	Id lineageIdCounter; // to ID wildtypes
	// ancestor tracing
	std::list<Host*> ancestors;
	void updateAncestors();
	// filestreams for logging
	void openFileStreams();
	void closeFileStreams();
    // file streams
	std::ofstream ofstat, ofherit, ofancs, oflong, oflonghd, ofcross,
			offun, offundot, ofimmresp, ofseq, ofcost, ofnet, ofassoc, ofaids;
    // 'loads' for printing at different time intervals
	double plstat, plherit, plancs, pllong, pllonghd, plcross, plseq,
			plcost, plnet, plassoc, plaids;
	// constants
	std::string identifier_base;
	std::string identifier; // for filenames
	int rank; // for multiple worlds
	int population_size; // disease-free steady-state population size // TODO: move to Network
	int cohort_size; // for cross-sectional sampling
};

// template methods

/* supposed to be a general mean function for a trait x...
 * HostTrait is a pointer-to-method for Host of type const void -> Numeric
 * e.g. call getMeanHostTrait(&Host::getLogSpvl)
 */

template<class Numeric>
double World::getMeanHostTrait(HostTrait<Numeric> X) const {
	double meanx = 0.0;
	int n = livingHosts.size();
	for ( Host* host : livingHosts ) {
		if ( host == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		meanx += (host->*X)(); // convert T to double... TODO use std::is_arithmatic<Numeric> and some assertion
	}
	if ( n > 0 ) {
		meanx /= n;
	}
	return meanx;
}

template<class Numeric>
double World::getSdHostTrait(HostTrait<Numeric> X) const {
	double meanx = getMeanHostTrait(X);
	double varx = 0.0;
	for ( Host* host : livingHosts ) {
		if ( host == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		double val = (host->*X)();
		varx += (val-meanx)*(val-meanx);
	}
	int m = livingHosts.size();
	double sdx = 0.0;
	if ( m >= 2 ) {
		varx /= (m-1);
		sdx = sqrt(varx);
	}
	return sdx;
}

template<class NumericL, class NumericR> // left and right
double World::getCorrelationHostTrait(HostTrait<NumericL> X, HostTrait<NumericR> Y) const {
	Correlation corr;
	for ( Host* host : livingHosts ) {
		if ( host == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		double x = (host->*X)(); // recall that X and Y are methods of class Host...
		double y = (host->*Y)();
		corr.addPoint(x, y);
	} // for every living host
	corr.computeStats();
	return corr.getCorrcoef();
}

template<class Numeric>
double World::getHeritabilityHostTrait(HostTrait<Numeric> X) const {
	Correlation corr;
	for ( Host* receiver : livingHosts ) {
		if ( receiver == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		Host* transmitter = receiver->getParent();
		if ( transmitter != nullptr ) { // valid transmission couple!
			double x = (transmitter->*X)();
			double y = (receiver->*X)();
			corr.addPoint(x, y);
		}
	}
	corr.computeStats(); // computes more than just the slope...
	return corr.getSlope();
}

template<class Numeric>
std::string World::xmlNodeHostTrait(HostTrait<Numeric> X, std::string name) const {
	BasicStats sts(name);
	for ( Host* host : livingHosts ) {
		if ( host == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		double x = (host->*X)();
		sts.addPoint(x);
	}
	sts.computeStats();
	std::stringstream ss; ss << sts;
	return ss.str();
}

template<class Numeric, class ArgType1, class ArgType2>
// FIXME: ArgType1 and ArgType2 should be identical!
std::string World::xmlNodeHostArgTrait(HostArgTrait<Numeric, ArgType1> X, ArgType2 arg, std::string name) const {
	BasicStats sts(name);
	for ( Host* host : livingHosts ) {
		if ( host == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		double x = (host->*X)(arg);
		sts.addPoint(x);
	}
	sts.computeStats();
	std::stringstream ss; ss << sts;
	return ss.str();
}

template<class Numeric, class Error, class ArgType1, class ArgType2>
// FIXME: ArgType1 and ArgType2 should be identical!
std::string World::xmlNodeHostMaybeArgTrait(HostMaybeArgTrait<Numeric, Error, ArgType1> X, ArgType2 arg, std::string name) const {
	BasicStats sts(name);
	for ( Host* host : livingHosts ) {
		if ( host == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		std::pair<Numeric, Error> x = (host->*X)(arg);
		if ( x.second ) { // Error type should evaluate as true
			sts.addPoint(double(x.first)); // explicit conversion of Numeric to double
		}
	}
	sts.computeStats();
	std::stringstream ss; ss << sts;
	return ss.str();
}


template<class NumericL, class NumericR>
std::string World::xmlNodeCorrelationHostTrait(HostTrait<NumericL> X, HostTrait<NumericR> Y, std::string name) const {
	Correlation corr(name);
	for ( Host* host : livingHosts ) {
		if ( host == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		double x = (host->*X)();
		double y = (host->*Y)();
		corr.addPoint(x, y);
	}
	corr.computeStats();
	corr.bootstrap(rng, NUMBER_OF_BOOTSTRAPS);
	std::stringstream ss; ss << corr;
	return ss.str();
}

template<class NumericL, class NumericR, class ArgTypeR1, class ArgTypeR2>
std::string World::xmlNodeCorrelationHostRArgTrait(HostTrait<NumericL> X,
			HostArgTrait<NumericR, ArgTypeR1> Y, ArgTypeR2 arg, std::string name) const {
	Correlation corr(name);
	for ( Host* host : livingHosts ) {
		if ( host == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		double x = (host->*X)();
		double y = (host->*Y)(arg);
		corr.addPoint(x, y);
	}
	corr.computeStats();
	corr.bootstrap(rng, NUMBER_OF_BOOTSTRAPS);
	std::stringstream ss; ss << corr;
	return ss.str();
}


template<class Numeric>
std::string World::xmlNodeHeritabilityHostTrait(HostTrait<Numeric> X, std::string name) const {
	Correlation corr(name);
	for ( Host* receiver : livingHosts ) {
		if ( receiver == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		} // else...
		Host* transmitter = receiver->getParent();
		if ( transmitter != nullptr ) { // valid transmission couple!
			double x = (transmitter->*X)();
			double y = (receiver->*X)();
			corr.addPoint(x, y);
		}
	}
	corr.computeStats();
	corr.bootstrap(rng, NUMBER_OF_BOOTSTRAPS);
	std::stringstream ss; ss << corr;
	return ss.str();
}

template<class Numeric, class Error>
std::string World::xmlNodeHostMaybeTrait(HostMaybeTrait<Numeric, Error> X, std::string name) const {
	BasicStats sts(name);
	for ( Host* host : livingHosts ) {
		if ( host == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		std::pair<Numeric, Error> x = (host->*X)();
		if ( x.second ) { // an error should always be equivalent to false. TODO Better test?
			sts.addPoint(x.first);
		}
	}
	sts.computeStats();
	std::stringstream ss; ss << sts;
	return ss.str();
}

#endif /* WORLD_HPP_ */
