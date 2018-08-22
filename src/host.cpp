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

#include "host.hpp"

// methods for DiseasePhase

DiseasePhase::DiseasePhase() { // set the default parameters
	chronicPhaseThreshold = 1.0;
	acutePhaseThreshold = FRASERS_D1;
	aidsPhaseThreshold = FRASERS_D3;
	chronicPhaseLoad = 0.0;
	acutePhaseLoad = 0.0;
	aidsPhaseLoad = 0.0;
	timeSinceInfection = 0.0;
	label = ACUTE_PHASE;
}

void DiseasePhase::init(Rng & rng) { // randomly generated thresholds
	chronicPhaseThreshold = rng.Gamma(1.0/FRASERS_RHO, FRASERS_RHO);
	if ( TRANS_BY_STAGE_OF_INFECT ) {
		acutePhaseThreshold = FRASERS_D1;
		aidsPhaseThreshold = FRASERS_D3;
	}	else { // just make acute and AIDS phase very short
		acutePhaseThreshold = 0.0;
		aidsPhaseThreshold = 0.0;
	}
	label = ACUTE_PHASE;
}

void DiseasePhase::reset() {
	chronicPhaseLoad = 0.0;
	acutePhaseLoad = 0.0;
	aidsPhaseLoad = 0.0;
	timeSinceInfection = 0.0;
	label = ACUTE_PHASE;
}

bool DiseasePhase::update(double dt, double vl) {
	bool new_disease_phase = false;
	switch ( label ) {
		case ACUTE_PHASE: {
			acutePhaseLoad += dt;
			if ( acutePhaseThreshold <= acutePhaseLoad ) {
				label = CHRONIC_PHASE;
				new_disease_phase = true;
			}
			break;
		}
		case CHRONIC_PHASE: {
			double dr = frasersDelta(vl);
			chronicPhaseLoad += dr * dt;
			if ( chronicPhaseThreshold <= chronicPhaseLoad ) {
				label = AIDS_PHASE;
				new_disease_phase = true;
			}
			break;
		}
		case AIDS_PHASE: {
			aidsPhaseLoad += dt;
			if ( aidsPhaseThreshold <= aidsPhaseLoad ) {
				label = DEAD_PHASE;
				new_disease_phase = true;
			}
			break;
		}
		case DEAD_PHASE: {
			std::cerr << "WARNING: dead host.\n";
			break;
		}
		default: {
			throw std::logic_error("unknown disease phase" + RIGHT_HERE);
			break; // redundant
		}
	}
	// for safety, stop the i-level sim after a maximum duration
	timeSinceInfection += dt;
	if ( timeSinceInfection >= MAX_TIME_ILEV ) {
		label = DEAD_PHASE;
		std::cerr << "WARNING: ilev model cut short as it was taking too long"
		          << RIGHT_HERE << std::endl;
		new_disease_phase = true;
	}
	return new_disease_phase;
}

// methods for HostState

HostState::HostState(const DiseasePhase & diseasePhase, const CiScheme & cis,
				const IndPars & ipars, double timeSinceInfection) {
	this->timeSinceInfection = timeSinceInfection;
	virusLoad = cis.getVirusLoad();
	clone_distribution = cis.getCloneDistribution();
	active_responses = cis.getActiveResponses();
	infectionRate = pLevInfectionRate(virusLoad, diseasePhase.label);

 	diseasePhaseLabel = diseasePhase.label;
	numberOfClones = cis.getNumberOfClones();
	numberOfMutants = cis.getNumberOfMutants();

	// compute diversity indices TODO: make specialized function (in CiScheme?)
	shannonIndex = 0.0;
	simpsonIndex = 0.0;
	for ( auto mit = clone_distribution.begin(); mit != clone_distribution.end(); ++mit ) {
		double p = mit->second;
		shannonIndex += -p*log(p);
		simpsonIndex += p*p;
	}

	// get statistics about immune responses
	numberOfResponses = cis.getNumberOfActiveResponses();
	numberOfLatents = cis.getNumberOfLatentResponses();
	numberOfMemories = cis.getNumberOfMemoryResponses();
	totalResponseSize = cis.getTotalResponseSize();
	fitnessCostEscape = cis.getFitnessCostOfEscape();

	// the number of target cells
	totalTargetSize = cis.getTargetValue();

	w = 0.0;
	for ( int i = 0; i < NUMBER_OF_GENES; ++i ) {
		ws[i] = 0.0;
	}
  r0 = 0.0; r0_generic = 0.0;
  R0 = 0.0; R0_generic = 0.0;
	// compute the average w, r0 and R0
	for ( auto mit = clone_distribution.begin(); mit != clone_distribution.end(); ++mit ) {
		const Virus* vir = mit->first;
		double p = mit->second;
        // w-parameters
		w += p * vir->getW();
		for ( int i = 0; i < NUMBER_OF_GENES; ++i ) {
			ws[i] += p * vir->getWGene(i);
		}
		// (malthusian) fitness parameters
    r0 += p * iLevBasicMalthusianFitness(*vir, ipars); // TODO: handled by CiScheme?
    r0_generic += p * iLevBasicMalthusianFitness(*vir, IndPars());
    R0 += p * iLevBasicReproductionNumber(*vir, ipars); // TODO: handled by CiScheme?
    R0_generic += p * iLevBasicReproductionNumber(*vir, IndPars());
	}
	// set the final_state_flag.
	final_state_flag = ( diseasePhaseLabel == DEAD_PHASE );
}

void HostState::print(std::ostream & os) const {
	os << "<host_state "
		<< "t='" << timeSinceInfection << "' "
		<< "Ttot='" << totalTargetSize << "' "
		<< "VL='" << virusLoad << "' "
		<< "w='" << w << "' ";
	for ( int i = 0; i < NUMBER_OF_GENES; ++i ) {
		os << geneNameMap.at(i) << "='" << ws[i] << "' ";
	}
	os << "r0='" << r0 << "' "
    << "r0_generic='" << r0_generic << "' "
    << "R0='" << R0 << "' "
    << "R0_generic='" << R0_generic << "' "
    << "richness='" << numberOfClones << "' "
		<< "shannon='" << shannonIndex << "' "
		<< "simpson='" << simpsonIndex << "' "
		<< "responses='" << numberOfResponses << "' "
		<< "latents='" << numberOfLatents << "' "
		<< "memories='" << numberOfMemories << "' "
		<< "Etot='" << totalResponseSize << "' "
		<< "cost_escape='" << fitnessCostEscape << "' "
		<< ">" << std::endl;
	// todo: more content?
	os  << "</host_state>";
}

bool HostState::isFinalState() const { return final_state_flag; }

// methods for HostSnapshot

HostSnapshot::HostSnapshot(const CiScheme & cis, double timeSinceInfection) {
	this->timeSinceInfection = timeSinceInfection;
	fitnessCosts = cis.getFitnessCosts();
	domFitnessCosts = cis.getDomFitnessCosts();
	fitnessCostsMimicResps = cis.getFitnessCosts(true); // true means "mimic resps parents"
}
void HostSnapshot::print(std::ostream & os) const {
	bool fst; // for correctly printing spaces
	// opening tag
	os << "<host_snapshot "
	   << "t='" << timeSinceInfection << "' "
		 << ">" << std::endl;

  // fitness costs
	os << "<mean_delta_r>" << std::endl;
	fst = true;
  for ( auto cost : fitnessCosts ) {
		os << (fst ? "" : " ") << cost;
		fst = false;
	}
	os << std::endl << "</mean_delta_r>" << std::endl;

	// fitness costs of dominant virus
	os << "<dom_delta_r>" << std::endl;
	fst = true;
	for ( auto cost : domFitnessCosts ) {
		os << (fst ? "" : " ") << cost;
		fst = false;
	}
	os << std::endl << "</dom_delta_r>" << std::endl;

	// fitness costs mimicing parents responses
	os << "<mean_delta_r_mimresp>" << std::endl;
	fst = true;
	for ( auto cost : fitnessCostsMimicResps ) {
		os << (fst ? "" : " ") << cost;
		fst = false;
	}
	os << std::endl << "</mean_delta_r_mimresp>" << std::endl;

	// closing tag
	os << "</host_snapshot>";
}


// methods for Host:...

Host::Host(Rng & rng, const std::vector<MhcLocus> & allMhcLoci,
		const Virus & virus, const IndPars & ipars, double timeOfInfection, Id id) :
			id(id), id_parent(ANON_ID), virusTF(virus), virusGsvl(virus),
			ipars(ipars), timeOfInfection(timeOfInfection) {
	// neutral SPVL mutation-at-transmission...
	if ( SHIRREFFS_NEUTRAL_MODEL ) {
		virusTF.mutateVlMultiplier(rng);
	}
	// Sample immune responses. TODO: let Agent sample immune responses?
	sampleImmuneResponses(immuneResponses, mhcHaplo, allMhcLoci, rng);
	diseasePhase.init(rng);
	// use the init() method to initialize the remaining members
	init(rng);
}

Host::Host(const Host & host, const Virus & virus, Rng & rng) :
		virusTF(virus), virusGsvl(virus) {
	// copy data from host
	id = host.id;
	id_parent = host.id_parent;
	immuneResponses = host.immuneResponses;
	ipars = host.ipars;
	timeOfInfection = host.timeOfInfection;
	diseasePhase = host.diseasePhase;
	diseasePhase.reset(); // keeps sampled threshods
	// use the init() method to initialize the remaining members
	init(rng);
}

void Host::init(Rng & rng) {
	// initialize the CiScheme
	cis.init(virusTF, immuneResponses, ipars, rng);
	// make a first state to get a valid state iterator
	HostState first_state(diseasePhase, cis, ipars, 0.0);
	states.push_back(first_state);
	cstit = states.begin();
	// parameters that are known right now
	R0TF = iLevBasicReproductionNumber(virusTF, ipars);
	genericR0TF = iLevBasicReproductionNumber(virusTF, IndPars()); // pass the default parameters
  malthusianTF = iLevBasicMalthusianFitness(virusTF, ipars);
  genericMalthusianTF = iLevBasicMalthusianFitness(virusTF, IndPars()); // pass the default parameters
	numberOfResponsesTF = getNumberOfResponses(virusTF);
	numberOfDomResponsesTF = getNumberOfDomResponses(virusTF);
	weightedNumberOfResponsesTF = getWeightedNumberOfResponses(virusTF);
	// parameters known only after within host model is completed
	cumulativeInfectionRate = 0.0;
	cumulativeSpvl = 0.0;
	cumulativeLogSpvl = 0.0;
	earlySpvl = 0.0;
	lengthChronicPhase = 0.0;
	gsvl = 0.0;
	gsvlDefined = false;
	escapeRate = 0.0;
	escapeCost = 0.0;
	numEscapes = 0;
	compensationRate = 0.0;
	numCompensations = 0;
	compensationRateTrunk = 0.0;
	numCompensationsTrunk = 0;
	slopeCD4 = 0.0; // TS median
	lowerSlopeCD4 = 0.0; // 2.5 percentile
	upperSlopeCD4 = 0.0; // 97.5 percentile
	slopeCD4Defined = false;
	/* initiate parameters used for ancestor tracing.
	 * make the parent NULL.
	 * if is_mrca is true, then a pointer to the host is stored by World
	 * In this case: no recursive deletion!
	 */
	parent = nullptr;
	remove_if_ext_lineage = false; // make true whenever host/agent dies
	is_mrca = false;
}


Host::~Host() {
	waitForJobToFinish(); // destroying a host with a running simulation might be a bad idea
	// when tracing ancestors, we may have to delete the parent
	if ( parent != nullptr ) {
		parent->removeChild(this);
		if ( parent->deleteMe() ) {
			delete parent;
		}
	}
	/* there may be children that need to be notified of the pending
	 * deletion of their parent
	 */
	std::list<Host*>::iterator hit = children.begin();
	for ( ; hit != children.end(); ++hit ) {
		Host* child = (*hit);
		if ( child != nullptr ) {
			child->parent = nullptr;
		} else {
			std::cerr << "ERROR: child is NULL" << RIGHT_HERE << std::endl;
			std::terminate();
		}
	}
}


// important members for the p-level simulations

const Virus & Host::getVirusGsvl() const {
	waitForJobToFinish(); // TODO GSVL virus is known earlier
	return virusGsvl;
}

const Virus & Host::getVirusTF() const {
	return virusTF;
}

const Virus & Host::getVirus() const {
	// find the most abundant strain
	const Virus* max_vir = nullptr;
	double max_fr = 0.0;
	for ( auto it = cstit->clone_distribution.begin();
				it != cstit->clone_distribution.end(); ++it ) {
		if ( it->second > max_fr ) {
			max_vir = it->first;
			max_fr = it->second;
		}
	}
	if ( max_vir == nullptr ) {
		throw std::logic_error("unable to find dominant virus, or NULL" + RIGHT_HERE);
	}
	return *max_vir;
}

const Virus & Host::getVirus(Rng & rng) const {
	// choose a number between 0 and 1
	double r = rng.Uniform();
	// choose a number from the strain distribution
	for ( auto it = cstit->clone_distribution.begin();
			it != cstit->clone_distribution.end(); ++it ) {
		r -= it->second;
		if ( r <= 0 ) {
			const Virus* samVirus = it->first;
			if ( samVirus == nullptr ) {
				throw std::logic_error("Virus in clone_distribution is NULL" + RIGHT_HERE);
			}
			return *samVirus;
		}
	}
	throw std::logic_error("unable to sample virus" + RIGHT_HERE);
}

const Virus & Host::getVirusSpvl(Rng & rng) const { // sample a virus from all state
	waitForJobToFinish();
	double r1 = rng.Uniform(0.0, cumulativeInfectionRate);
	for ( auto & state : states ) {
		r1 -= state.infectionRate;
		if ( r1 <= 0.0 ) {
			double r2 = rng.Uniform();
			for ( auto & pair : state.clone_distribution ) {
				r2 -= pair.second;
				if ( r2 <= 0.0 ) {
					const Virus* samVirus = pair.first;
					if ( samVirus == nullptr ) {
						throw std::logic_error("Virus in clone_distribution is NULL" + RIGHT_HERE);
					}
					return *samVirus;
				}
			} // for loop over clone distribution
		}
	} // for loop over states
	throw std::logic_error("unable to sample virus" + RIGHT_HERE);
}

bool Host::nextState() {
	/** Increase the state-iterator cstit to point to the next state.
	 * Only do this when the current state is not the final state
	 * (as can be verified with the HostState::isFinalState() method).
	 * Otherwise, the next state may not have been computed yet.
	 * In this case, wait until the i-level simulation to compute the
	 * next state.
	 */
	if ( cstit->isFinalState() ) {
		return false; // don't increment if the current state is the final state
	}	else { // more states are comming...
		/* claim the nextState mutex to see if there is a valid next state
		 * in the states list. If not, wait until these is one added
		 */
		std::unique_lock<std::mutex> lock(nextStateMutex);
		auto pred = [this]{
			auto cstit_copy = cstit; cstit_copy++;
			return cstit_copy != states.end();
		};
		nextStateConditionalVar.wait(lock, pred);
 		// there is a next state, so increase the state iterator...
		cstit++;
		return true;
	}
}

bool Host::isAlive() const {
	// the iterator cstit is not going to be increased passed the final state
	return ( cstit->diseasePhaseLabel != DEAD_PHASE );
}

// method for running the within-host simulation
void Host::run(Rng & rng, bool make_host_dump) {
	/** Simple loop to run the i-level model. Summaries of the states
	 * are stored in std::list<HostState> states
	 * NB: CiScheme cis owns the virus pointers...
	 */
	double dtL = PLEV_TIMESTEP;
	double dtS = PLEV_TIMESTEP / ILEV_PER_PLEV_STEPS;
	double t = 0.0;

	// parameter used for periodically taking computationally intensive snapshots
	double plsnapshot = 0.0;

	// open an xml element
	if ( make_host_dump ) {
		hostdump << "<host_dump >" << std::endl;
		hostdump << cis << std::endl; // the first state
	}

	// for estimating the CD4 slope, we need a Correlation object
	Correlation corr_cd4;

	// run until DEAD
	while ( diseasePhase.label != DEAD_PHASE ) {
		// update the CI scheme
		if ( make_host_dump ) { // take smaller time steps and record the state
			for ( int i = 0; i < ILEV_PER_PLEV_STEPS; ++i ) {
				cis.updateHostState(dtS, rng);
				hostdump << cis << std::endl;
			}
		}	else { // otherwise take only one step
			cis.updateHostState(dtL, rng);
		}
		// update time
		t += dtL;

		// virus load needed for disease phase and cumulative infection rate, etc.
		double vl = cis.getVirusLoad(); // per ml blood

		// update disease phase
		bool new_disease_phase = diseasePhase.update(dtL, vl);

		// update cumulatives
		if ( diseasePhase.label == CHRONIC_PHASE ) {
			cumulativeSpvl += vl * dtL;
			if ( vl > 0.0 ) {
				cumulativeLogSpvl += log(vl) * M_LOG10E * dtL;
			} else {
				std::cerr << "WARNING: virus load is non-positive, can't compute log" << RIGHT_HERE << std::endl;
			}
			lengthChronicPhase += dtL;
		}
		if ( new_disease_phase && diseasePhase.label == CHRONIC_PHASE ) {
			earlySpvl = vl;
		}
		cumulativeInfectionRate += pLevInfectionRate(vl, diseasePhase.label) * dtL;
		if ( !gsvlDefined && (t >= GSVL_SAMPLE_TIME || diseasePhase.label == DEAD_PHASE) ) {
			/* a host could die before it is time to sample GSVL.
			 * in that case, define the GSLV as the last known VL.
			 * However, we don't set the gsvlDefined flag to true.
			 */
			if ( diseasePhase.label != DEAD_PHASE ) {
				gsvlDefined = true;
			}
			gsvl = vl;
			virusGsvl = cis.getDominantClone();
		}
		// create a new HostState
		HostState state(diseasePhase, cis, ipars, t);
		/* claim the nexState mutex before adding the latest state to states,
		 * then broadcast the nextState conditional variable and finally
		 * unlock the mutex
		 */
		nextStateMutex.lock();
		states.push_back(state);
		nextStateConditionalVar.notify_all();
		nextStateMutex.unlock();

		// every HOST_SNAPSHOT_INTERVAL years, store a HostSnapshot
		plsnapshot += dtL;
		if ( plsnapshot >= HOST_SNAPSHOT_INTERVAL ) {
			plsnapshot -= HOST_SNAPSHOT_INTERVAL; // correct for overshoot
			HostSnapshot snaps(cis, t);
			snapshots.push_back(snaps);
		}

		// add a observation for the CD4 slope (time unit per day, count per uL)
		if ( t > GSVL_SAMPLE_TIME ) { // skip acute phase (TODO: choose good start)
			corr_cd4.addPoint(t * DAYS_IN_YEAR, cis.getTargetValue() / UL_IN_BLOOD);
		}
	} // while alive
	// post-i-level simulation statistics...
	timeOfDeath = timeOfInfection + t;
	corr_cd4.computeStats();
	slopeCD4 = corr_cd4.getTheilSenSlope(); // TS works better when there are outliers
	lowerSlopeCD4 = corr_cd4.getLowerTheilSenSlope();
	upperSlopeCD4 = corr_cd4.getUpperTheilSenSlope();
	slopeCD4Defined = corr_cd4.getSize() >= 2; // we need at least 2 observations for a slope

	// compute statistics using the virus records

	// 1st argument of getReplacementRate: use escape mutants
	auto esc_rate_pair = cis.getReplacementRate(true);
	escapeRate = esc_rate_pair.first;
	numEscapes = esc_rate_pair.second;
	// 2nd argument of getReplacementRate: mimic responses parents
	auto esc_cost_pair = cis.getReplacementRate(true, true);
	escapeCost = esc_cost_pair.first;
	if ( numEscapes != esc_cost_pair.second ) {
		throw std::logic_error("CiScheme::getReplacementRate is not consistant" + RIGHT_HERE);
	}
	auto comp_rate_pair = cis.getReplacementRate(false);
	compensationRate = comp_rate_pair.first;
	numCompensations = comp_rate_pair.second;
	// 3rd argument of getReplacementRate: use only strains in the trunk of the ancestor graph (TODO: as it is: just the internal nodes)
	auto comp_rate_trunk_pair = cis.getReplacementRate(false, false, true);
	compensationRateTrunk = comp_rate_trunk_pair.first;
	numCompensationsTrunk = comp_rate_trunk_pair.second;


	if ( make_host_dump ) {
		// dictionary of virus IDs and the genomes
		cis.printVirusCache(hostdump);
		hostdump << std::endl;
		// graph of i-level ancestors
		hostdump << "<ancestor_graph >" << std::endl;
		cis.printAncestorGraph(hostdump);
		hostdump << std::endl << "</ancestor_graph>" << std::endl;
		// close the xml element
		hostdump << "</host_dump>";
	}

	// cleanup to save memory
	cis.clearAllStochasticViruses();
}

// statistics (host traits)

double Host::getTimeOfInfection() const { return timeOfInfection; }
double Host::getTimeSinceInfection() const { return cstit->timeSinceInfection; }

double Host::getTimeOfDeath() const {
	waitForJobToFinish();
	return timeOfDeath;
}

double Host::getLengthOfInfection() const {
	waitForJobToFinish();
	return timeOfDeath - timeOfInfection;
}

double Host::getLogSpvl() const {
	waitForJobToFinish();
	double average_vl = cumulativeSpvl / lengthChronicPhase;
	if ( average_vl <= 0.0 ) {
		throw std::logic_error("can't compute logarithm" + RIGHT_HERE);
	}
	return log(average_vl) * M_LOG10E;
}

double Host::getLogGeomSpvl() const {
	waitForJobToFinish();
	return cumulativeLogSpvl / lengthChronicPhase;
}


double Host::getEarlyLogSpvl() const {
	waitForJobToFinish(); // early spvl is actually know earlier...
	if ( earlySpvl <= 0.0 ) {
		throw std::logic_error("can't compute logarithm" + RIGHT_HERE);
	}
	return log(earlySpvl) * M_LOG10E;
}

double Host::getLogGsvl() const {
	waitForJobToFinish(); // gsvl is actually know earlier...
	if ( gsvl <= 0.0 ) {
		throw std::logic_error("can't compute logarithm" + RIGHT_HERE);
	}
	return log(gsvl) * M_LOG10E;
}

bool Host::isGsvlDefined() const {
	waitForJobToFinish(); // gsvl is actually know earlier...
	return gsvlDefined;
}

double Host::getLogVlMultiplierTF() const {
	double vlm = virusTF.getVlMultiplier();
	if ( vlm <= 0.0 ) {
		throw std::logic_error("can't compute logarithm" + RIGHT_HERE);
	}
	return  log(vlm) * M_LOG10E;
}

double Host::getWTF() const { return virusTF.getW(); }
double Host::getWTF(int i) const { return virusTF.getWGene(i); }
double Host::getR0TF() const { return R0TF; }
double Host::getGenericR0TF() const { return genericR0TF; }
double Host::getMalthusianTF() const { return malthusianTF; }
double Host::getGenericMalthusianTF() const { return genericMalthusianTF; }
int Host::getNumberOfResponsesTF() const { return numberOfResponsesTF; }
int Host::getNumberOfDomResponsesTF() const { return numberOfDomResponsesTF; }
double Host::getWeightedNumberOfResponsesTF() const {	return weightedNumberOfResponsesTF; }

int Host::getNumberOfResponses(const Virus & virus) const {
	auto respCountFun = [&](const ImmuneResponse & ir) { return ir.hasEpitope(virus); };
	return std::count_if(immuneResponses.begin(), immuneResponses.end(), respCountFun);
}

int Host::getNumberOfDomResponses(const Virus & virus) const {
	auto domRespCountFun = [&](const ImmuneResponse & ir) { return ir.hasEpitope(virus) && ir.isDominant(); };
	return std::count_if(immuneResponses.begin(), immuneResponses.end(), domRespCountFun);
}

double Host::getWeightedNumberOfResponses(const Virus & virus) const {
	double n = 0.0;
	for ( auto & ir : immuneResponses ) {
		if ( ir.hasEpitope(virus) ) n += ir.getDominance();
	}
	return n;
}

double Host::getPreAdaptationTF(const Virus & virus) const {
	return numberOfResponsesTF - getNumberOfResponses(virus);
}

double Host::getPreAdaptationTF_vm(const std::map<Id, Virus> & viruses) const {
	try {
		return getPreAdaptationTF(viruses.at(virusTF.getLineageId()));
	} catch ( const std::out_of_range & oor ) {
		throw std::logic_error("out-of-range exception while finding virus in viruses" + RIGHT_HERE);
	}
}

double Host::getWeightedPreAdaptationTF(const Virus & virus) const {
	return weightedNumberOfResponsesTF - getWeightedNumberOfResponses(virus);
}

double Host::getWeightedPreAdaptationTF_vm(const std::map<Id, Virus> & viruses) const {
	try {
		return getWeightedPreAdaptationTF(viruses.at(virusTF.getLineageId()));
	} catch ( const std::out_of_range & oor ) {
		throw std::logic_error("out-of-range exception while finding virus in viruses" + RIGHT_HERE);
	}
}

std::pair<double, bool> Host::getRelPreAdaptationTF(const Virus & virus) const {
	int n = getNumberOfResponses(virus);
	if ( n == 0 ) { // not well defined
		return std::make_pair(0.0, false);
	}
	return std::make_pair(double(numberOfResponsesTF)/n - 1.0, true);
}

std::pair<double, bool> Host::getRelPreAdaptationTF_vm(const std::map<Id, Virus> & viruses) const {
	try {
		return getRelPreAdaptationTF(viruses.at(virusTF.getLineageId()));
	} catch ( const std::out_of_range & oor ) {
		throw std::logic_error("out-of-range exception while finding virus in viruses" + RIGHT_HERE);
	}
}

std::pair<double, bool> Host::getRelWeightedPreAdaptationTF(const Virus & virus) const {
	int n = getNumberOfResponses(virus); // used to test if well defined
	if ( n == 0 ) { // not well defined
		return std::make_pair(0.0, false);
	}
	double wn = getWeightedNumberOfResponses(virus);
	return std::make_pair(weightedNumberOfResponsesTF/wn - 1.0, true);
}

std::pair<double, bool> Host::getRelWeightedPreAdaptationTF_vm(const std::map<Id, Virus> & viruses) const {
	try {
		return getRelWeightedPreAdaptationTF(viruses.at(virusTF.getLineageId()));
	} catch ( const std::out_of_range & oor ) {
		throw std::logic_error("out-of-range exception while finding virus in viruses" + RIGHT_HERE);
	}
}

double Host::getLogWTF() const { return log(virusTF.getW()); }
double Host::getW() const { return cstit->w; }
double Host::getW(int i) const { return cstit->ws[i]; }
double Host::getVl() const { return cstit->virusLoad; }

double Host::getLogVl() const {
	double vl = getVl();
	if ( vl <= 0.0 ) {
		throw std::logic_error("can't compute logarithm" + RIGHT_HERE);
	}
	return log(vl) * M_LOG10E;
}

double Host::getPLevMalthusian(double s) const {
	/* s = S/N, the fraction of susceptible individuals
	 * Returns: beta(V) * S/N - delta(V).
	 * Trait can be used for the Price equation
	 */
	double vl = getVl();
	return frasersBeta(vl) * s - frasersDelta(vl); // Malthusian fitness
}


double Host::getInfectionRate() const {
	return pLevInfectionRate(cstit->virusLoad, cstit->diseasePhase.label);
}

double Host::getIndividualReproductionNumber() const {
	waitForJobToFinish();
	return cumulativeInfectionRate;
}

std::pair<double, bool> Host::getEscapeRate() const {
	waitForJobToFinish();
	bool ok = (numEscapes > 0);
	return std::make_pair(escapeRate, ok);
}

std::pair<double, bool> Host::getEscapeCost() const {
	waitForJobToFinish();
	bool ok = (numEscapes > 0);
	return std::make_pair(escapeCost, ok);
}

std::pair<double, bool> Host::getCompensationRate() const {
	waitForJobToFinish();
	bool ok = (numCompensations > 0);
	return std::make_pair(compensationRate, ok);
}

std::pair<double, bool> Host::getCompensationRateTrunk() const {
	waitForJobToFinish();
	bool ok = (numCompensationsTrunk > 0);
	return std::make_pair(compensationRateTrunk, ok);
}

double Host::getEscapeDensity() const {
	waitForJobToFinish();
	double duration = timeOfDeath - timeOfInfection;
	return numEscapes / duration;
}

double Host::getCompensationDensity() const {
	waitForJobToFinish();
	double duration = timeOfDeath - timeOfInfection;
	return numCompensations / duration;
}

double Host::getCompensationDensityTrunk() const {
	waitForJobToFinish();
	double duration = timeOfDeath - timeOfInfection;
	return numCompensationsTrunk / duration;
}

double Host::getSlopeCD4() const {
	waitForJobToFinish();
	return slopeCD4;
}

double Host::getLowerSlopeCD4() const {
	waitForJobToFinish();
	return lowerSlopeCD4;
}

double Host::getUpperSlopeCD4() const {
	waitForJobToFinish();
	return upperSlopeCD4;
}

double Host::getScaledHammingDistance(const Virus & virus) const {
	double dH = hammingDistance(getVirus(), virus); // getVirus waits for job to finish
	return dH / GENOME_SIZE;
}

double Host::getScaledHammingDistance_vm(const std::map<Id, Virus> & viruses) const {
	try {
		return getScaledHammingDistance(viruses.at(virusTF.getLineageId()));
	} catch ( const std::out_of_range & oor ) {
		throw std::logic_error("out-of-range exception while finding virus in viruses" + RIGHT_HERE);
	}
}

double Host::getFitnessCostOfEscape() const {
	return cstit->fitnessCostEscape;
}

double Host::getRelBasicMalthusian(int loc) const {
	/* return the average fitness const of a mutation at the given locus loc.
	 * TODO: this is a bit inefficient, since this has all been computed
	 * computed before by CiScheme.
	 */
	double dr0 = 0.0;
	for ( auto & pair : cstit->clone_distribution ) {
		/* pair.first is a pointer to the virus
		 * pair.second is the current frequency of that virus
		 */
		Virus mut(*pair.first, loc);
		double r0_wt = iLevBasicMalthusianFitness(*pair.first, ipars);
		double r0_mut = iLevBasicMalthusianFitness(mut, ipars);
		dr0 += pair.second * (r0_mut - r0_wt); // costs are negative...
	}
	return dr0;
}

double Host::getEarlyRelMalthusian(int loc) const {
	// TODO: return additional boolean to indicate success
	waitForJobToFinish();
	if ( snapshots.empty() ) {
		return 0.0;
	}
	return snapshots.front().fitnessCosts[loc];
}

double Host::getLateRelMalthusian(int loc) const {
	// TODO: return additional boolean to indicate success
	waitForJobToFinish();
	if ( snapshots.empty() ) {
		return 0.0;
	}
	return snapshots.back().fitnessCosts[loc];
}

double Host::getEarlyRelMalthusianMimicResps(int loc) const {
	// TODO: return additional boolean to indicate success
	waitForJobToFinish();
	if ( snapshots.empty() ) {
		return 0.0;
	}
	return snapshots.front().fitnessCostsMimicResps[loc];
}

double Host::getLateRelMalthusianMimicResps(int loc) const {
	// TODO: return additional boolean to indicate success
	waitForJobToFinish();
	if ( snapshots.empty() ) {
		return 0.0;
	}
	return snapshots.back().fitnessCostsMimicResps[loc];
}

double Host::getEarlyDomRelMalthusian(int loc) const {
	// TODO: return additional boolean to indicate success
	waitForJobToFinish();
	if ( snapshots.empty() ) {
		return 0.0;
	}
	return snapshots.front().domFitnessCosts[loc];
}

double Host::getLateDomRelMalthusian(int loc) const {
	// TODO: return additional boolean to indicate success
	waitForJobToFinish();
	if ( snapshots.empty() ) {
		return 0.0;
	}
	return snapshots.back().domFitnessCosts[loc];
}


int Host::getNumberOfClones() const {
	return cstit->numberOfClones;
}

int Host::getNumberOfMutants() const {
	return cstit->numberOfMutants;
}

const MhcHaplotype & Host::getMhcHaplo() const { return mhcHaplo; }

bool Host::hasMhcAllele(const MhcID & mhcID) const {
	return std::find(mhcHaplo.begin(), mhcHaplo.end(), mhcID) != mhcHaplo.end();
}

const IndPars & Host::getIndPars() const { return ipars; }

// auxiliary methods for host

void Host::print(std::ostream & os) const {
	waitForJobToFinish(); // only print when the i-level simulation is over...
	// opening tag
	os << "<host "
	   << "time_of_infection='" << timeOfInfection << "' "
	   << "time_of_death='" << timeOfDeath << "' "
		 << "id='" << id << "' "
		 << "id_parent='" << id_parent << "' "
     << "R0_TF='" << R0TF << "' "
     << "R0_generic_TF='" << genericR0TF << "' "
     << "r0_TF='" << malthusianTF << "' "
     << "r0_generic_TF='" << genericMalthusianTF << "' "
     << "responses_TF='" << numberOfResponsesTF << "' "
     << "dom_responses_TF='" << numberOfDomResponsesTF << "' "
	   << "R_indiv='" << getIndividualReproductionNumber() << "' "
	   << "gsvl='" << getLogGsvl() << "' "
	   << "gsvl_defined='" << std::boolalpha << gsvlDefined << std::noboolalpha << "' "
	   << "spvl='" << getLogSpvl() << "' "
		 << "geom_spvl='" << getLogGeomSpvl() << "' "
		 << "escape_rate='" << escapeRate << "' "
 		 << "escape_cost='" << escapeCost << "' "
		 << "num_escapes='" << numEscapes << "' "
		 << "compensation_rate='" << compensationRate << "' "
		 << "num_compensations='" << numCompensations << "' "
		 << "compensation_rate_trunk='" << compensationRateTrunk << "' "
		 << "num_compensations_trunk='" << numCompensationsTrunk << "' "
		 << "slope_CD4='" << slopeCD4 << "' "
		 << "lower_slope_CD4='" << lowerSlopeCD4 << "' "
 		 << "upper_slope_CD4='" << upperSlopeCD4 << "' "
		 << "slope_CD4_defined='" << std::boolalpha << slopeCD4Defined << std::noboolalpha << "' "
		 << ">" << std::endl;
  // Individual parameter variations
  os << ipars << std::endl;
	// MHC haplotype
	os << xmlStringMhcHaplo(mhcHaplo) << std::endl;
	// Immune responses
	os << immuneResponses << std::endl;
	// TF virus
	os << "<virus_tf >" << "\n"
	   << virusTF << "\n"
	   << "</virus_tf>" << std::endl;
	// GSVL virus (equals virusTF if host died before GSLV was sampled)
	os << "<virus_gsvl >" << "\n"
	   << virusGsvl << "\n"
	   << "</virus_gsvl>" << std::endl;
	// states
	for ( auto sit = states.begin(); sit != states.end(); ++sit ) {
		os << (*sit) << std::endl;
	}
	// snapshots
	for ( auto sit = snapshots.begin(); sit != snapshots.end(); ++sit ) {
		os << (*sit) << std::endl;
	}
	// print the "host dump"
	std::string dumpstring = hostdump.str();
	if ( !dumpstring.empty() ) {
		os << dumpstring << std::endl;
	}
	os << "</host>";
}

void Host::printTree(std::ostream & os) const {
	for ( auto cit = children.begin(); cit != children.end(); ++cit ) {
		Host* child = (*cit);
		if ( child == nullptr ) {
			throw std::logic_error("child is NULL" + RIGHT_HERE);
		} // else...
		child->printTree(os);
	}
	print(os);
}

int Host::getTreeSize() const {
  int result = 1; // namely, this Host
  for ( auto cit = children.begin(); cit != children.end(); ++cit ) {
    Host* child = (*cit);
    if ( child == nullptr ) {
			throw std::logic_error("child is NULL" + RIGHT_HERE);
		} // else...
		result += child->getTreeSize();
  }
  return result;
}

// functions for ancestor tracing

void Host::setAddressParent(Host* parent) {
	if ( parent == nullptr ) {
		throw std::invalid_argument("trying to set address of parent with NULL" + RIGHT_HERE);
	} // else...
	this->parent = parent;
	id_parent = parent->id; // used for plotting the ancestor tree
}

void Host::addChild(Host* child) {
	if ( child == nullptr ) {
		throw std::invalid_argument("passed NULL as a child" + RIGHT_HERE);
	}
  children.push_back(child);
	child->setAddressParent(this);
}

void Host::removeChild(Host* child) {
	children.remove(child);
}

std::pair<Host*, bool> Host::getChild() const {
	if ( children.empty() ) return std::make_pair(nullptr, false);
	else return std::make_pair(children.front(), true);
}

Host* Host::getParent() const {
	return parent; // NB could be nullptr!
}

void Host::setRemoveIfExtLineageFlag(bool x) {
	remove_if_ext_lineage = x;
}

bool Host::getRemoveIfExtLineageFlag() const {
	return remove_if_ext_lineage;
}

void Host::setMrcaFlag(bool x) {
	is_mrca = x;
}

bool Host::getMrcaFlag() const {
	return is_mrca;
}

bool Host::deleteMe() const {
	/* allowed to be deleted recursively: Host is dead,
	 * no offspring lineages, not under control of ancestor tracing
	 */
	return ( !is_mrca && remove_if_ext_lineage && children.empty() );
}

bool Host::coalesces() const {
	// test that this host is a common ancestor of two or more hosts.
	return ( children.size() > 1 );
}

// method inherited from Job...
bool Host::execute(Rng & rng) {
	run(rng);
	return true;
}

// auxiliary functions for class Host

std::string xmlStringMhcHaplo(const MhcHaplotype & mhcHaplo) {
	std::stringstream ss;
	ss << "<mhc_haplotype >" << std::endl;
	for ( auto it = mhcHaplo.begin(); it != mhcHaplo.end(); ++it ) {
		int locus = it->first;
		int allele = it->second;
		ss	<< "<mhc_id "
		    << "locus='" << locus << "' "
				<< "allele='" << allele << "' "
				<< "/>" << std::endl;
	}
	ss << "</mhc_haplotype>";
	return ss.str();
}

double pLevInfectionRate(double vl, DiseasePhaseLabel diseasePhaseLabel) {
	if ( TRANS_BY_STAGE_OF_INFECT ) { // use Hollingsworths transmission model
		switch ( diseasePhaseLabel ) {
			case ACUTE_PHASE: return FRASERS_BETA1;
			case CHRONIC_PHASE: return frasersBeta(vl);
			case AIDS_PHASE: return FRASERS_BETA3;
			case DEAD_PHASE: return 0.0;
			default: {
				throw std::logic_error("unknown disease phase" + RIGHT_HERE);
				return 0.0; // redundant
			}
		}
	}	else { // use Frasers beta_2(V) for the entire infectious period
		return frasersBeta(vl);
	}
}
