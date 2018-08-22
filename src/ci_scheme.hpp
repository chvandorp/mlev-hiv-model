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

#ifndef CISCHEME_HPP_
#define CISCHEME_HPP_

#include <list>
#include <map>
#include <unordered_set> // obsolete?
#include <unordered_map>
#include <algorithm> // std::find
#include <iostream>

#include "macros.hpp"
#include "aux.hpp" // for Id, ...
#include "rng.hpp" // for sampling thresholds
#include "immune_response.hpp"
#include "virus.hpp"
#include "integrator.hpp"
#include "vertices.hpp" // user-defined vertices
#include "hashed_definitions.hpp"
#include "ilevel_equations.hpp" // IndPars

enum VirusEvent { VIRUS_TRANSMIT=0, VIRUS_BIRTH, VIRUS_DEATH };

struct VirusRecord {
	/** during the simulation, keep track of escape and compensation events.
	 * What was the escape rate? Did the mutant have an escape epitope?
	 */
	VirusEvent event;
	double t;
	Id id;
	std::list<Id> ids_parent; // for birth events
	double r, R; // growth rate and fitness at event (in current environment)
	double r_relative; // growth rate relative to parent
	double r_mimresp; // relative to parent, mimicing the parents responses TODO
	bool escape; // was the mutant an immune escape?
};

// VirusPtrHash and VirusPtrEq are declared in Virus.hpp
typedef std::unordered_set<Virus*, VirusPtrHash, VirusPtrEq> VirusSet; // obsolete
typedef std::unordered_map<Virus*, Id, VirusPtrHash, VirusPtrEq> VirusIdMap;
typedef std::unordered_map<Virus*, VirusVertex*, VirusPtrHash, VirusPtrEq> VirusMap;
typedef std::unordered_map<Virus*, std::list<VirusRecord>, VirusPtrHash, VirusPtrEq> VirusRecordMap;
// for now, we only need a simple container for responses
typedef std::map<ImmuneResponse*, ResponseVertex*> ResponseMap;


/** CiScheme uses Graph to integrate a system of ODEs
 * with multiple strains and immune responses
 *
 * TODO/issues:
 * 2) Keep track of deterministic virus using a separate list,
 *    and use this list to quickly compute statistics. Normally, most strains
 *    are stochastic, so the complete set of viruses is large
 */

class CiScheme {
public:
	CiScheme(); // need default constructor for Host member
	// CiScheme(const Virus & virusTF, const ImmuneSystem immsys, IndPars & ipars, Rng & rng);
	/* non-default constructor is now replaced by an init member function.
	 * TODO: check out how move constructors work...:
	 * CiScheme(CiScheme&& ); // move-constructor
	 * CiScheme & operator=(CiScheme&& ); // move-assignment operator
	 */
	~CiScheme();
	void init(const Virus & virusTF, const ImmuneSystem immsys,
			const IndPars & ipars, Rng & rng);
	void updateHostState(double dt, Rng & rng);
	double getTime() const; // returns t
	// statistics for the virus
	int getNumberOfClones() const; // i.e. deterministic virus
	int getNumberOfMutants() const; // i.e. stochastic virus (in ODE system)
	int getNumberOfTracedViruses() const; // viruses.size()
	int getNumberOfCachedViruses() const; // deterVirusCache.size()
	// returns the total VL (per ml blood, using vl_multiplier)
	double getVirusLoad() const;
	double getFitnessCostOfEscape() const; // TODO: compute per-escape (return list)
	// Delta r in current environment
	std::vector<double> getFitnessCosts(bool mimic_resps_parent=false) const;
	// Delta r in current environment w.r.t dominant virus. TODO: make a single function
	std::vector<double> getDomFitnessCosts(bool mimic_resps_parent=false) const;
	// returns map Virus* to rel freq (NB not using vl_multiplier)
	std::map<const Virus*, double> getCloneDistribution() const; // NB: pointers owned by CiScheme
	std::map<const ImmuneResponse*, double> getActiveResponses() const;
	double getTargetValue() const; // T
	double getQuiescentValue() const; // Q
	const Virus & getDominantClone() const; // TODO a private version that returns the vertex
	// statistics for the immune responses
	int getNumberOfActiveResponses() const; // TODO: make more intelligent
	int getNumberOfLatentResponses() const; // TODO: make more intelligent
	int getNumberOfMemoryResponses() const;
	double getTotalResponseSize() const; // sum of effector sizes
	/* analyse virus records. escape and compenstation rate.
	 * return a pair with average replacement rate and the number of relacements.
	 * replacements are escapes if the argument is true (default) and reversions
	 * and compensatations when the argument is false.
	 * When the third optinal argument is true, the average is only taken over
	 * internal nodes (i.e. leafs are excluded)
	 */
	std::pair<double, int> getReplacementRate(bool esc=true,
			bool mimic_resps_parent=false, bool internal_nodes=false) const;
	/* in order to save memory, it is possible to clear the VirusMap viruses
	 * from all stochastic strains. Typically done when the i-level simulation
	 * has ended. TODO: make a method to resume the simulation after clearing
	 */
	void clearAllStochasticViruses();
	// printing
	void printState(std::ostream & ) const; // xml
	void printVirusCache(std::ostream & ) const; // xml
	void printVirusRecords(std::ostream & ) const; // xml
	void printAncestorGraph(std::ostream & os) const; // graphviz
protected:
// aux functions
	/* Stoch -> Deter, Deter -> Stoch
	 * NULL -> Stoch (mutants), Stoch -> NULL (dead parents)
	 */
	void updateViruses(Rng & rng);
	void updateResponses(Rng & rng);
	// more aux functions
	void addMutants(DeterVirusVertex* , Rng & rng); // aux function for updateViruses, but also for init
	void removeOrphans(); // delete any stoch virus without a parent
	void informStochVirus(ActiveResponseVertex* ); // add a new active response to Vertex::edgesIn
	void informLatentResponses(DeterVirusVertex* ); // latent responses need to know about a new deter virus
	Id cacheVirus(Virus* );
	std::pair<Id, bool> getVirusId(Virus* ) const; // snd argument is false when the virus has not been cached
	void addVirusRecord(DeterVirusVertex* , Virus*, VirusEvent );
	std::unordered_set<Id> getInternalNodesAncestorGraph() const; // anyone virus that is not a leaf
// data
	VirusMap viruses; // the virus population (map Virus* -> VirusVertex*)
	std::list<DeterVirusVertex*> newparents; // keep track of new parents, untill they reach a threshold
	ResponseMap responses; // all immune responses (map ImmuneResponse -> ResponseVertex*)
	TargetVertex* target; // Vertex in integrationGraph corresponding to Target cells
	IndPars ipars; // the i-level parameters of the Host
	Graph integrationGraph;
	double t_rem; // remainder, i.e. overshoot after last virus/response update
	double h; // time-step
	int numberOfClones; // the number of deterministic viruses
	int numberOfResponses; // the number of active responses
	VirusIdMap deterVirusCache; // all "seen" deterministic viruses, with a unique ID
	Id nextId;
	/* the pointers in VirusRecordMap are 'owned' by deterVirusCache.
	 * no deletion is needed! TODO: use a pointer-based hash, since all viruses will
	 * be unique anyway.
	 */
	VirusRecordMap virusRecords;
private:
	bool is_initialized; // FIXME: flag to avoid double init calls, make better system
};

std::ostream & operator<<(std::ostream & , const CiScheme & ); // uses CiScheme::print


#endif
