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

#ifndef HOST_HPP
#define HOST_HPP

#include <algorithm> // std::find
#include <iostream>
#include <sstream>
#include <list>
#include <vector>
#include <utility> // std::pair
#include <stdexcept>

#include "macros.hpp"
#include "aux.hpp" // use typedef Id
#include "parallelism.hpp" // let Host inherit the class "Job"
#include "virus.hpp"
#include "immune_response.hpp"
#include "ci_scheme.hpp"
#include "frasers_rates.hpp"
#include "hashed_definitions.hpp"
#include "ilevel_equations.hpp" // for IndivParams
#include "calc_stats.hpp" // CD4 down slope

enum DiseasePhaseLabel { ACUTE_PHASE=0, CHRONIC_PHASE, AIDS_PHASE, DEAD_PHASE };

class DiseasePhase {
	/** the DiseasePhase class is used for determining the length of
	 * the infection. It originates from the stage-of-infection model
   * by Hollingsworth et al. Disabled when TRANS_BY_STAGE_OF_INFECT
   * is false, in that case only the chronic phase is used
	 */
public:
	DiseasePhase();
	void init(Rng & );
	void reset(); // reset loads, label, no re-sampling of the thresholds
	bool update(double , double ); // dt, vl. returns true when the status changed
	DiseasePhaseLabel label;
private:
	double acutePhaseThreshold;
	double chronicPhaseThreshold;
	double aidsPhaseThreshold;
	double acutePhaseLoad;
	double chronicPhaseLoad;
	double aidsPhaseLoad;
	double timeSinceInfection; // for "safety"
};

class HostState : public Printable {
public:
	HostState(const DiseasePhase & , const CiScheme & , const IndPars & , double );
	std::map<const Virus*, double> clone_distribution;
	std::map<const ImmuneResponse*, double> active_responses;
	DiseasePhase diseasePhase;
	ImmuneSystem immuneResponses;
	double timeSinceInfection;
	double dt_cache;
	DiseasePhaseLabel diseasePhaseLabel;
	double totalTargetSize; // the number of target cells
	// viral stats
	double w; // average
	double ws[NUMBER_OF_GENES];
  double r0, r0_generic; // average
  double R0, R0_generic; // average
	double virusLoad;
	double infectionRate;
	int numberOfClones; // deterministic virus
	int numberOfMutants; // stochastic virus (in ODE system)
	double shannonIndex;
	double simpsonIndex;
	// immune stats
	int numberOfResponses; // active responses
	int numberOfLatents; // latent responses
	int numberOfMemories; // latents that were once active
	double totalResponseSize; // sum_e E_e
	double fitnessCostEscape; // average...
	// auxiliary methods
	void print(std::ostream & ) const;
	bool isFinalState() const; // returns final_state_flag
private:
	bool final_state_flag; // indicates whether this state is the final state
};

class HostSnapshot : public Printable {
	/* HostSnapshot is a container for periodically, but harder to compute
	 * host data.
	 */
public:
	HostSnapshot(const CiScheme & , double ); // constructor: cis, time since infection
	void print(std::ostream & ) const;
	// data
	std::vector<double> fitnessCosts; // per locus
	std::vector<double> domFitnessCosts; // per locus w.r.t dominant clone
	std::vector<double> fitnessCostsMimicResps; // per locus mimicing parent responses
	double timeSinceInfection;
private:
	/* empty */
};

class Host : public Job, public Printable {
public:
	Host(Rng & , // for sampling from the immune responses and running the simulations
		 const std::vector<MhcLocus> & , // global set of loci
		 const Virus & , // virus has a pointer to the fitness function
 		 const IndPars & , // a set of i-level parameters for additional host-heterogeneity
		 double , // time of infection
	   Id id=ANON_ID); // optionally, inherit the Id of the Agent
	/* a copy constructor for Host copies all host properties,
	 * but resets all states. Also the disease phase thresholds are copied.
	 */
	Host(const Host & host, const Virus & virus, Rng & rng);
	/* destructor: take care of the ancestor trace
	 */
	~Host();
	// i-level
	void run(Rng & , bool make_host_dump=false);
	// p-level
	const Virus & getVirus() const; // most common virus
	const Virus & getVirus(Rng & ) const; // sample a virus (from current state)
	const Virus & getVirusTF() const;
	const Virus & getVirusGsvl() const;
	const Virus & getVirusSpvl(Rng & ) const; // sample a virus from all states
	bool nextState(); // returns false when there are no more states left
	bool isAlive() const; // returns false when the diseasePhaseLabel is DEAD_PHASE
	// statistics (host traits) TODO: sort according to i-level constant or not
	double getInfectionRate() const;
	double getVl() const; // instantaneous virus load
	double getLogVl() const; // instantaneous log10 virus load
	double getPLevMalthusian(double ) const; // beta(V) * S/N - delta(V). Used for Price equation
	double getTimeOfInfection() const; // this is a typical constant i-level trait
	double getTimeSinceInfection() const; // this is a typical non-constant i-level trait
	double getTimeOfDeath() const;
	double getLengthOfInfection() const; // timeOfDeath - timeOfInfection
	double getLogSpvl() const;
	double getLogGeomSpvl() const; // geometric time average
	double getLogGsvl() const;
	bool isGsvlDefined() const;
	double getEarlyLogSpvl() const; // after acute phase (only makes sense for stage-of-infection model)
	double getLogVlMultiplierTF() const;
	double getWTF() const; // fitness parameter of the transmitter/founder virus
	double getWTF(int ) const; // gene-specific fitness parameter of the transmitter/founder virus
	double getLogWTF() const;
	double getW() const; // average fitness parameter of current strains
	double getW(int ) const; // gene-specific (average) fitness parameter of current strains
	double getFitnessCostOfEscape() const;
	// measures of per-locus fitness cost TODO: combine into fewer methods...
	double getRelBasicMalthusian(int ) const; // return average fitness cost (Delta r_0) of mutation at given locus
	double getEarlyRelMalthusian(int ) const; // return average fitness cost (Delta r) of mutation at given locus early in infection
	double getLateRelMalthusian(int ) const; // return average fitness cost (Delta r) of mutation at given locus late in infection
	double getEarlyRelMalthusianMimicResps(int ) const; // return average fitness cost (Delta r) of mutation at given locus early in infection, mimicing responses of the parent
	double getLateRelMalthusianMimicResps(int ) const; // return average fitness cost (Delta r) of mutation at given locus late in infection, mimicing responses of the parent
	double getEarlyDomRelMalthusian(int ) const; // return fitness cost (Delta r) of mutation of dominant clone at given locus early in infection
	double getLateDomRelMalthusian(int ) const; // return fitness cost (Delta r) of mutation of dominant clone at given locus late in infection
	// cumulative infection rate ($\int\beta$)
	double getIndividualReproductionNumber() const;
	int getNumberOfClones() const; // count the number of clones
	int getNumberOfMutants() const; // count the number of VIABLE mutants (i.e. in ODE system)
	double getR0TF() const; // the basic reproduction number of the TF-virus
	double getGenericR0TF() const; // the basic reproduction number of the TF-virus for the average host
  double getMalthusianTF() const; // the Malthusian fitness of the TF-virus
  double getGenericMalthusianTF() const; // the Malthusian fitness of the TF-virus for the average host
	int getNumberOfResponsesTF() const; // the max number of possible responses against the TF virus
	int getNumberOfDomResponsesTF() const; // the max number of possible dominant responses against the TF virus
	double getWeightedNumberOfResponsesTF() const; // the number of possible responses against the TF virus, weighted by dominance
	int getNumberOfResponses(const Virus & ) const; // the max number of possible responses against a given virus
	int getNumberOfDomResponses(const Virus & ) const; // the max number of possible dominant responses against a given virus
	double getWeightedNumberOfResponses(const Virus & ) const; // the number of possible responses against a given virus, weighted by dominance
	double getPreAdaptationTF(const Virus & ) const; // difference between number of responses against TF virus, and argument
	double getPreAdaptationTF_vm(const std::map<Id, Virus> & ) const; // difference between number of responses against TF virus, and argument
	double getWeightedPreAdaptationTF(const Virus & ) const; // immunodominance-weighted difference between number of responses against TF, and argument
	double getWeightedPreAdaptationTF_vm(const std::map<Id, Virus> & ) const; // immunodominance-weighted difference between number of responses against TF, and argument
	std::pair<double, bool> getRelPreAdaptationTF(const Virus & ) const; // relative difference between number of responses against TF, and argument
	std::pair<double, bool> getRelPreAdaptationTF_vm(const std::map<Id, Virus> & ) const; // relative difference between number of responses against TF, and argument
	std::pair<double, bool> getRelWeightedPreAdaptationTF(const Virus & ) const; // immunodominance-weighted relative difference between number of responses against TF, and argument
	std::pair<double, bool> getRelWeightedPreAdaptationTF_vm(const std::map<Id, Virus> & ) const; // immunodominance-weighted relative difference between number of responses against TF, and argument
	std::pair<double, bool> getEscapeRate() const; // growth rate of deterministic escape mutants
	std::pair<double, bool> getEscapeCost() const; // growth rate of deterministic escape mutants in the context of the parents' responses
	std::pair<double, bool> getCompensationRate() const; // growth rate of deterministic compensation/reversion mutants
	std::pair<double, bool> getCompensationRateTrunk() const; // growth rate of deterministic compensation/reversion mutants in the trunk of the ancestor graph
	double getEscapeDensity() const; // number of escapes / length of infection
	double getCompensationDensity() const; // number of compensatory mutations / length of infection
	double getCompensationDensityTrunk() const; // number of compensatory mutations / length of infection
	double getSlopeCD4() const; // an estimate of the CD4 downslope
	double getLowerSlopeCD4() const; // an estimate of the CD4 downslope (2.5 percentile)
	double getUpperSlopeCD4() const; // an estimate of the CD4 downslope (97.5 percentile)
	// auxiliary trait to compute the distance distribution from e.g. the wildtype virus
	double getScaledHammingDistance(const Virus & ) const; // distance between current dominant ctrain and the argument
	double getScaledHammingDistance_vm(const std::map<Id, Virus> & ) const; // distance between current dominant ctrain and the argument
	// properties
	const MhcHaplotype & getMhcHaplo() const;
	bool hasMhcAllele(const MhcID & ) const;
	const IndPars & getIndPars() const;
	// auxiliary
	void print(std::ostream & ) const;
	// methods for ancestor tracing
	void setAddressParent(Host* );
	void addChild(Host* );
	void removeChild(Host* );
	std::pair<Host*, bool> getChild() const;
    // return the first child, second item of pair is false whenever there were no children
	Host* getParent() const; // also used for heritability calculations
	bool deleteMe() const;
	// true if 1) remove_if_ext_lineage == true, 2) children.empty() is true, 3) is_mrca == false
	void setRemoveIfExtLineageFlag(bool );
	bool getRemoveIfExtLineageFlag() const;
	void setMrcaFlag(bool );
	bool getMrcaFlag() const;
	bool coalesces() const; // equivalent to children.size() > 1
	void printTree(std::ostream & ) const; // print this host and all children (recursively)
  int getTreeSize() const; // total number of grand^n children
protected:
	Id id; // TODO
	Id id_parent; // set when setAddressParent is called
	/* Host can be handled by WorkerPool.
	 * implement the virtual method Job::execute here...
	 */
	bool execute(Rng & );
	// the states and a pointer to the current state
	std::list<HostState> states;
	std::list<HostState>::const_iterator cstit; // Const STate ITerator
	// mutexes and conditional variables
	std::mutex nextStateMutex; // has the next state been computed?
	std::condition_variable nextStateConditionalVar; // for waiting until the next state has been computed
	// periodic data that needs a lower frequency than the HostState
	std::list<HostSnapshot> snapshots;
	// further data
	Virus virusTF; // transmitter/founder virus
	Virus virusGsvl; // dominant strain during GSVL measurement
	IndPars ipars; // individual parameter variation
	MhcHaplotype mhcHaplo;
	ImmuneSystem immuneResponses;
	CiScheme cis;
	DiseasePhase diseasePhase;
	double timeOfInfection;
	double timeOfDeath;
	double lengthChronicPhase;
	double cumulativeSpvl;
	double cumulativeLogSpvl;
	double earlySpvl; // right after acute phase (Hollingsworth model)
	double gsvl; // gold standard virus load
	bool gsvlDefined; // when the Host dies before gsvl could ve determined, this boolean is false
	double cumulativeInfectionRate;
  double R0TF; // R0 of the transmitter-founder virus at infection
  double genericR0TF; // R0 of the transmitter-founder virus at infection for the average host
  double malthusianTF; // the Malthusian fitness for the TF virus at infection
  double genericMalthusianTF; // the Malthusian fitness for the TF virus at infection for the average host
  int numberOfResponsesTF; // the number of responses against the TF virus
	int numberOfDomResponsesTF; // number of dominant responses against TF virus
	double weightedNumberOfResponsesTF; // number of responses against the TF virus, weighted by immuno-dominance
	double escapeRate; // average (relative) growth rate of dominant escape mutants
	double escapeCost; // relative growth rate of dominant escape mutants, using parents' responses
	int numEscapes; // number of escapes that occured during infection
	double compensationRate; // average (relative) growth of dominant reversions/compensations
	int numCompensations; // number of compensations/reversions that have occured during infection
	double compensationRateTrunk; // average (relative) growth of dominant reversions/compensations that are internal nodes of the ancestor graph (TODO: trunk)
	int numCompensationsTrunk; // number of compensations/reversions that have occured during infection (TODO: trunk)
	double slopeCD4; // estimation of the CD4+ T-cell decline during chronic infection
	double lowerSlopeCD4; // 2.5 percentile
	double upperSlopeCD4; // 97.5 percentile
	bool slopeCD4Defined; // we only start sampling CD4 counts after an initial phase, hence the slope might not be defined

	// full i-level simulation logging
	std::stringstream hostdump;

	// pointers for ancestor tracing
	std::list<Host*> children;
	Host* parent;
	bool remove_if_ext_lineage; // when the lineage has ended, it is ok to delete the host
	bool is_mrca; // typically, we don't want to delete such hosts recursively
private:
	void init(Rng & rng); // auxiliary init function used by constructors
};

// auxiliary functions
std::string xmlStringMhcHaplo(const MhcHaplotype & );

/* get the p-level infection rate for a given virus load and disease phase.
 * if TRANS_BY_STAGE_OF_INFECT is true, Hollingsworth and Frasers beta1 and
 * beta3 are used for the acute and AIDS phase (resp.). Otherwise beta2(vl)
 * is used for the entire infectious period.
 */
double pLevInfectionRate(double , DiseasePhaseLabel );

/* World has members that compute statistics about the population.
 * When e.g. the mean of a host trait x is computed, the Host member
 * x needs to be passed to the stat function. Therefore, we define
 * the type HostTrait as a (const) Host::<method> to double
 */
template <typename Numeric>
using HostTrait = Numeric(Host::*)() const; // C++11 equivalent of typedef (allowing templates)

// some traits take additional arguments
template <typename Numeric, typename ArgType>
using HostArgTrait = Numeric(Host::*)(ArgType) const;

// some traits are not always defined
template <typename Numeric, typename Error>
using HostMaybeTrait = std::pair<Numeric, Error>(Host::*)() const;

// trait that can fail and take arguments
template <typename Numeric, typename Error, typename ArgType>
using HostMaybeArgTrait = std::pair<Numeric, Error>(Host::*)(ArgType) const;

#endif /* HOSTCLASS_HPP_ */
