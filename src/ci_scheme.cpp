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

#include "ci_scheme.hpp"

/********************** Methods for CiScheme *******************/

CiScheme::CiScheme() { // need default constructor for Host member
	target = nullptr;
	t_rem = 0.0;
	h = INIT_STEPSIZE_CISCHEME;
	numberOfClones = 0;
	numberOfResponses = 0;
	nextId = 0;
	is_initialized = false;
}

void CiScheme::init(const Virus & virusTF, const ImmuneSystem imm, const IndPars & ip, Rng & rng) {
	/* TODO: clear any old content? Yes, but for now,
	 * check that init has not been called before
	 */
	if ( !is_initialized ) {
		is_initialized = true;
	} else {
		throw std::logic_error("called init method a second time" + RIGHT_HERE);
	}
	// initialize members...
	ipars = ip;
	t_rem = 0.0;
	h = INIT_STEPSIZE_CISCHEME;
	nextId = 0; // the ID of the TF virus

	Virus* vir = new Virus(virusTF); // copy of the TF virus
	numberOfClones = 1;
	DeterVirusVertex* dvv = new DeterVirusVertex(vir, &ipars);

	Vec initial = iLevInitialCondition(*vir, ipars).first;
	// second returned value is the growth rate
	dvv->assignValue(initial.first, 0); // TODO: "inoculumsize"
	dvv->assignValue(initial.second, 1);
	dvv->assignThreshold(THRESHOLD_ILEV_ENDANGERED); // redundant
  dvv->endangered = false; // dont consider deleting the clone just yet... TODO set/get functions
  // add the virus to the virus - vertex map
 	viruses[vir] = dvv;
	// when the TF reaches a threshold, start keeping track of mutants
	newparents.push_back(dvv);

	TargetVertex* tvert = new TargetVertex(&ipars);
	tvert->assignValue(ipars.Q0, 0); // quiescent cells
	tvert->assignValue(ipars.T0, 1); // activated cells
	target = tvert;
	integrationGraph.addUndirectedEdge(dvv, tvert); // double edge

	// NB assignInOdeSys only works after adding the vertices to Graph
	integrationGraph.assignInOdeSys(dvv, true);
	integrationGraph.assignInOdeSys(tvert, true);

	// add immune responses to the graph
	for ( auto iit = imm.begin(); iit != imm.end(); ++iit ) {
		ImmuneResponse* resp = new ImmuneResponse(*iit); // copy of the response
		// passing false to the constructor of LatentResponseVertex implies naive (not memory)
		LatentResponseVertex* lrv = new LatentResponseVertex(resp, &ipars, false);
		// add the (response, lrv) pair to the ResponseMap
		responses[resp] = lrv;
		lrv->assignThreshold(rng.Exponential(1.0));
		if ( iit->hasEpitope(*vir) ) {
			integrationGraph.addEdge(dvv, lrv);
		}	else {
			integrationGraph.addVertex(lrv);
		}
		// assignInOdeSys only works after adding the vertices to Graph
		integrationGraph.assignInOdeSys(lrv, false);
	}

	// caching and recording...
	Id vir_id = cacheVirus(vir);
	addVirusRecord(dvv, vir, VIRUS_TRANSMIT);
}

CiScheme::~CiScheme() {
	// delete responses
	for ( auto it = responses.begin(); it != responses.end(); ++it ) {
		ImmuneResponse* resp = it->first;
		delete resp;
		ResponseVertex* rvert = it->second;
		delete rvert;
	}
	// delete viruses
	for ( auto it = viruses.begin(); it != viruses.end(); ++it ) {
		Virus* vir = it->first;
		delete vir;
		VirusVertex* vvert = it->second;
		delete vvert;
	}
	// clear the deterVirusCache
	for ( auto it = deterVirusCache.begin(); it != deterVirusCache.end(); ++it ) {
		Virus* vir = it->first;
		delete vir;
	}
	// delete target
	delete target;
}

void CiScheme::updateHostState(double dt_yr, Rng & rng) {
	// convert units of time
	double dt = dt_yr * DAYS_IN_YEAR;
	// advance the within-host state
	if ( dt < t_rem ) { // we don't have to update viruses or responses
		integrationGraph.integrate(dt);
		t_rem -= dt;
	}
	else {
		if ( t_rem > 0.0 ) { // left over from previous update
			// integrate the first bit
			integrationGraph.integrate(t_rem); // FIRST INTEGRATE, THEN UPDATE.
			updateViruses(rng);
			updateResponses(rng);
			dt -= t_rem;
		}
		while ( dt >= h ) {
			// update viruses and responses
			integrationGraph.integrate(h); // FIRST INTEGRATE, THEN UPDATE.
			updateViruses(rng);
			updateResponses(rng);
			dt -= h;
		}
		// integrate the last bit
		if ( dt > 0.0 ) {
			integrationGraph.integrate(dt);
			t_rem = h - dt;
		}
		else {
			t_rem = 0.0;
		}
	}
}

int CiScheme::getNumberOfClones() const { // i.e. deterministic virus
	return numberOfClones;
}
int CiScheme::getNumberOfMutants() const { // i.e. stochastic virus (in ODE system)
	auto pred = [](std::pair<Virus*, VirusVertex*> pair){
		return pair.second->type == Vertex::STOCH_VIRUS && pair.second->isInOdeSys();
	};
	return std::count_if(viruses.begin(), viruses.end(), pred);
}
int CiScheme::getNumberOfTracedViruses() const {
	return viruses.size();
}
int CiScheme::getNumberOfCachedViruses() const {
	return deterVirusCache.size();
}
double CiScheme::getTime() const { return integrationGraph.getTime(); }

double CiScheme::getVirusLoad() const {
	// we first need the number of target cells.
	double T = getTargetValue();
	// then compute the VL
	double VL = 0.0;
	for ( auto it = viruses.begin(); it != viruses.end(); ++it ) {
		VirusVertex* vvert = it->second;
		// filter deterministic virus
		if ( vvert->getType() == Vertex::DETER_VIRUS ) {
			DeterVirusVertex* dvv = (DeterVirusVertex*) vvert;
			/* the dvv compute the VL. dvv knows its parameters
			 * pass T so that dvv does not have to look for TargetVertex
			 */
			double V = dvv->getQssVL(T);
			VL += V / ML_IN_BLOOD * it->first->getVlMultiplier();
		}
	}
	return VL;
}

double CiScheme::getFitnessCostOfEscape() const {
	/** loop over all stochastic virus
	 * check if they are escape w.r.t a parent.
	 * compute several (TODO) measures of fitness cost, w.r.t parent(s).
	 * take average. (TODO: other stats)
	 * TODO: use (and create) VirusVertex methods to make this easier
	 */
	double expectedCost = 0.0;
	int escapes = 0; // number of escape mutants
	for ( auto it = viruses.begin(); it != viruses.end(); ++it ) {
		Virus* virus = it->first;
		VirusVertex* vvert = it->second;
		if ( vvert->getType() == Vertex::STOCH_VIRUS ) {
			StochVirusVertex* svv = (StochVirusVertex*) vvert;
			// get fitness... for average parameters (TODO ???)
			double fitnessMut = iLevBasicMalthusianFitness(*virus, ipars);
			// count responses
			auto resps = svv->getNumberOfResponses(); // returns pair<int, double>
			int n = resps.first;
			double E = resps.second;
			// get the virus' parents and count their respones
			auto parents = svv->getParents();
			double fitnessWt = 0.0;
			int escapeFromParent = 0;
			for ( auto jt = parents.begin(); jt != parents.end(); ++jt ) {
				auto resps_prime = (*jt)->getNumberOfResponses();
				int n_prime = resps_prime.first;
				double E_prime = resps_prime.second;
				// FIXME: what are the possiblities and rules for escape?
				if ( n < n_prime ) {
					fitnessWt += iLevBasicMalthusianFitness(*((*jt)->vir), ipars);
					escapeFromParent++;
				}
			}
			// FIXME: is it an escape if one of the parents also has the mutation?
			if ( escapeFromParent > 0 ) {
				fitnessWt /= escapeFromParent;
				expectedCost += (fitnessMut - fitnessWt); // costs are negative...
				escapes++;
			}
		}
	}
	if ( escapes > 0 ) {
		return expectedCost / escapes;
	}	else {
		return 0.0;
	}
}

std::vector<double> CiScheme::getFitnessCosts(bool mimic_resps_parent) const {
	/* What is the fitness cost of a mutation at each locus?
	 * Loop through stochastic viruses, find their parents,
	 * find at what locus they differ from their parent,
	 * and get the fitness cost (Delta r) from the Graph.
	 * NB: deterministic strains that have not yet reached the
	 * parent threshold are not counted. See the comments
	 * at CiScheme::addMutants for more details...
	 */
	double T = getTargetValue();
	double Vtot = 0.0;
	/* for the VL, we first need the number of target cells. We have to
	 * normalize with the total VL afterwards. The total VL that we need
	 * here can differ from the "real" VL, because we should not count the
	 * deterministic strains that have not reached the parent threshold.
	 * we therefore count the VL of the parents of every child, and correct
	 * for the length of the genome. TODO: prove that this is correct
	 * TODO: when growthrate is extremely slow, there is no StochVirus
	 * at the first sampling moment. We then return 0 (we check for this condition)
	 * Return a Boolean to signal that everything went well...
	 */
	std::vector<double> costs(GENOME_SIZE, 0.0);
	int numStochVir = 0; // PATCH: when there are no stochastic viruses, return 0
	for ( auto it = viruses.begin(); it != viruses.end(); ++it ) {
		Virus* virus = it->first;
		VirusVertex* vvert = it->second;
		if ( vvert->getType() == Vertex::STOCH_VIRUS ) {
			numStochVir++;
			// find parents
			auto parents = vvert->getParents(); // list of DeterVirusVertex*
			for ( auto dvv : parents ) {
				Virus* parent = dvv->getVirus();
				/* find where the virus differs from the parent.
				 * NB: findDifference returns the first locus that is differs between
				 * parent and child. It does not check for other differences, which
				 * would be an error in the current algorithm. This might change when
				 * recombination is added...
				 */
				auto pr = parent->findDifference(*virus);
				if ( !pr.second ) { // test that parent and child are not identical
					throw std::logic_error("ERROR: parent and child virus are identical" + RIGHT_HERE);
				}
				int loc = pr.first;
				double r_parent = dvv->getR(true); // argument "true" means Malthusian
				double r_child = 0.0;
				// we can compute fitness assuming responses of the parent...
				if ( mimic_resps_parent ) {
					r_child = vvert->getR(true, dvv);
				} else { // or just use the proper responses
					r_child = vvert->getR(true); // second argument is nullptr by default
				}
				// costs are negative (and from the prespective of the parent)
				double delta_r = r_child - r_parent;
				// weigh the fitness const with the parent frequency...
				double V = dvv->getQssVL(T);
				costs[loc] += V * delta_r;
				Vtot += V; // we have to normalize below.
			} // for loop over parents
		} // select stochastic virus...
	} // ...from for loop over all viruses
	if ( numStochVir == 0 ) {
		return costs; // [0, 0, ..., 0]
	} // else: normalize with Vtot (if Vtot > 0)
	if ( Vtot <= 0 ) {
		throw std::logic_error("about to devide by zero" + RIGHT_HERE);
	}
	for ( double & cost : costs ) {
		cost *= (GENOME_SIZE / Vtot); // correct for counting parents 'double'
	}
	return costs;
}

std::vector<double> CiScheme::getDomFitnessCosts(bool mimic_resps_parent) const {
	/* see also comments at CiScheme::getFitnessCosts.
	 * find the domiant strain, loop through its mutants (children)
	 * determine the locus where they differ and compute a fitness difference
	 */
	std::vector<double> costs(GENOME_SIZE, 0.0);
	// find dominant clone (and vertex)
	double T = getTargetValue();
	// find clone with the largest V
	double Vmax = 0.0;
	VirusVertex* parent = nullptr;
	for ( auto it = viruses.begin(); it != viruses.end(); ++it ) {
		VirusVertex* vvert = it->second;
		// filter deterministic virus
		if ( vvert->getType() == Vertex::DETER_VIRUS ) {
			// get the QSS VL
			DeterVirusVertex* dvv = (DeterVirusVertex*) vvert;
			// let dvv compute VL, pass T so that dvv doesn't have to look for it
			double V = dvv->getQssVL(T);
			if ( V > Vmax ) {
				parent = vvert;
				Vmax = V;
			}
		}
	}
	if ( parent == nullptr ) {
		std::cerr << "WARNING: unable to find dominant clone" << RIGHT_HERE << std::endl;
		return costs;
	} // else continue computation...
	// find the children of the dominant clone
	auto children = parent->getAllChildren(); // returns Stoch and Deter Virus
	// NB: there is also a getChildren method that returns only StochVirusVertex
	Virus* virus = parent->getVirus();
	for ( auto child : children ) {
		Virus* mutant = child->getVirus();
		auto pr = virus->findDifference(*mutant);
		if ( !pr.second ) {
			throw std::logic_error("parent and child virus are identical" + RIGHT_HERE);
		}
		int loc = pr.first;
		double r_parent = parent->getR(true); // true = Malthusian
		double r_child = 0.0;
		// we can compute fitness assuming responses of the parent...
		if ( mimic_resps_parent ) {
			r_child = child->getR(true, parent);
		} else { // or just use the proper responses
			r_child = child->getR(true); // second argument is nullptr by default
		}
		double delta_r = r_child - r_parent; // costs are negative and w.r.t parent
		costs[loc] = delta_r;
	}
	return costs;
}


std::map<const Virus*, double> CiScheme::getCloneDistribution() const {
	// for the VL, we first need the number of target cells.
	double T = getTargetValue();
	// now search for deterministic viruses
	std::map<const Virus*, double> clone_distribution;
	double Vtot = 0.0; // normalizing constant
	// loop over the virus population
	for ( auto it = viruses.begin(); it != viruses.end(); ++it ) {
		Virus* virus = it->first;
		VirusVertex* vvert = it->second;
		// filter deterministic virus
		if ( vvert->getType() == Vertex::DETER_VIRUS ) {
			// get the QSS VL
			DeterVirusVertex* dvv = (DeterVirusVertex*) vvert;
			// let dvv compute VL, pass T so that dvv doesn't have to look for it
			double V = dvv->getQssVL(T);
			// find the corresponding pointer in the deterVirusCache
			auto cit = deterVirusCache.find(virus);
			clone_distribution[cit->first] = V;
			Vtot += V;
		}
	}
	// normalize using Vtot
	if ( Vtot == 0.0 && clone_distribution.size() > 0 ) {
		throw std::logic_error("about to devide by zero" + RIGHT_HERE);
	} // else: all is fine...
	for ( auto it = clone_distribution.begin(); it != clone_distribution.end(); ++it ) {
		it->second /= Vtot;
	}
	return clone_distribution; // the compiler should be smart enough not to copy anything here...
}


std::map<const ImmuneResponse*, double> CiScheme::getActiveResponses() const {
	/* loop through all responses, check if they have a active response vertex
	 * and return a map with active responses and their clone sizes
	 */
	std::map<const ImmuneResponse*, double> active_responses;
	for ( auto it = responses.begin(); it != responses.end(); ++it ) {
		const ImmuneResponse* resp = it->first;
		const ResponseVertex* rvert = it->second;
		if ( rvert->getType() == Vertex::ACTIVE_RESPONSE ) {
			double E = rvert->getValue();
			active_responses[resp] = E;
		}
	}
	return active_responses;
}

double CiScheme::getTargetValue() const {
	return target->getValue(1);
}

double CiScheme::getQuiescentValue() const {
	return target->getValue(0);
}

const Virus & CiScheme::getDominantClone() const {
	auto clone_distribution = getCloneDistribution();
	double fr_max = 0.0;
	const Virus* vir_max = nullptr;
	for ( auto it = clone_distribution.begin(); it != clone_distribution.end(); ++it ) {
		double fr = it->second;
		const Virus* vir = it->first;
		if ( fr > fr_max ) {
			fr_max = fr;
			vir_max = vir;
		}
	}
	if ( vir_max == nullptr ) {
		throw std::logic_error("could not find a dominant clone" + RIGHT_HERE);
	}
	return *vir_max;
}

int CiScheme::getNumberOfActiveResponses() const {
	// TODO: make more intelligent: count responses against deterministic virus
	return numberOfResponses;
}

int CiScheme::getNumberOfLatentResponses() const {
	// TODO: make more intelligent: count responses against deterministic virus
	return int(responses.size()) - numberOfResponses;
}

int CiScheme::getNumberOfMemoryResponses() const {
	auto pred = [](std::pair<ImmuneResponse*, ResponseVertex*> pair) {
		if ( pair.second->getType() == Vertex::LATENT_RESPONSE ) {
			return ((LatentResponseVertex*) pair.second)->isMemory();
		}
		else return false;
	};
	return std::count_if(responses.begin(), responses.end(), pred);
}

double CiScheme::getTotalResponseSize() const {
	// compute the sum of all immune responses
	double Etot = 0.0;
	for ( auto it = responses.begin(); it != responses.end(); ++it ) {
		Vertex* vert = it->second;
		if ( vert->getType() == Vertex::ACTIVE_RESPONSE ) {
			Etot += vert->getValue();
		}
	}
	return Etot;
}

/*** methods for updating CiScheme ***/

void CiScheme::updateViruses(Rng & rng) {
	/* check if viruses have reached thresholds:
	 *  1a) stochastic strains may be introduced and become deterministic
	 *  1b) stochastic strains may be viable and added to the system of ODEs
	 *  2) deterministic strains may go extinct
   */
  bool removed_deter_virus = false; // have there been extinctions?
	for ( auto it = viruses.begin(); it != viruses.end(); ++it ) { // don't modify the VirusSet!
		Virus* vir = it->first;
		VirusVertex* vvert = it->second;
		switch ( vvert->getType() ) {
			case Vertex::STOCH_VIRUS: {
				StochVirusVertex* svv = (StochVirusVertex*) vvert; // TODO: typecast needed?
				double R = svv->getR();
				if ( !svv->isInOdeSys() ) {
					// should the stochastic mutant be included to the ODEs?
					if ( R > 1.0 ) { // yes!
						integrationGraph.assignInOdeSys(svv, true);
					}
				} else { // re-evaluate membership of ODE system
					if ( R <= 1.0 ) { // set the stochastic strain "on hold"
						integrationGraph.assignInOdeSys(svv, false);
					}	else { // R > 1.0, now check if threshold has been reached...
						// make the stochastic strain deterministic?
						if ( svv->getValue() > svv->getThreshold() ) {
							/* introduce a new deterministic strain by...
						 	 * 1) let the StochVirusVertex make a deterministic copy
						 	 * 2) remove the stochastic vertex
						 	 * 3) do some book keeping
						 	 */
							// creating a deterministic "copy" of svv
							DeterVirusVertex* dvv = svv->mkDeterVirusVertex(integrationGraph);
							double r = dvv->assignInitialCondition(1.0);
							dvv->assignThreshold(THRESHOLD_ILEV_ENDANGERED); // redundant
              dvv->endangered = false;
							// removing the stochastic node from the graph
							integrationGraph.removeVertex(svv);
							// replacing the pointer associated by vir in viruses
							it->second = dvv;
							// cleanup of stochastic vertex (nobody would own it otherwise...)
							delete svv;
							// add vir to the list of new potentential "parents"
							newparents.push_back(dvv); // TODO: restrict the set of parents?
							// add the new deterministic virus to the cache
							Id vir_id = cacheVirus(vir); //TODO use vir_id?
							// make a VirusRecord with some basic info
							addVirusRecord(dvv, vir, VIRUS_BIRTH);
							// keep the deterministic clone counter up-to-date
							numberOfClones++;
							// inform latent responses about the new deterministic virus
							informLatentResponses(dvv);
						} // if threshold reached
					}
				} // if/else not ssv in system of ODEs
				break;
			} // case STOCH_VIRUS
			case Vertex::DETER_VIRUS: {
				DeterVirusVertex* dvv = (DeterVirusVertex*) vvert;
        // do we need to give the clone an "endangered" status?
        if ( !dvv->endangered ) {
          if ( dvv->getValueSum() < THRESHOLD_ILEV_ENDANGERED ) {
						// FIXME: get the right proportions I1, I2
            double R = dvv->getR();
            if ( R < 1.0 ) {
              double u = rng.Uniform();
              double tau = extinction_threshold(u, dvv->getValueSum(), R);
              dvv->assignThreshold(tau);
              dvv->endangered = true;
            }
          } // else: all is fine...
        } else { // the clone is marked "endangered"... delete it?
          // don't delete the last strain (optional)
          if ( dvv->getValueSum() < dvv->getThreshold() ) {
            if ( numberOfClones > 1 || ALLOW_ILEV_EXTINCTION ) {
							// before removing the DeterVirusVertex, make a record
							addVirusRecord(dvv, vir, VIRUS_DEATH);
            	/* strain is extinct, add a stochastic copy (back mutations)
							 * NB: this is strictly not nescesary when dvv does not have
							 * any parents, but then we would have to remove the iterator
							 * from viruses here.
							 */
              StochVirusVertex* svv = dvv->mkStochVirusVertex(integrationGraph);
              svv->assignThreshold(rng.Exponential(1.0));
              // remove from Graph
              integrationGraph.removeVertex(dvv);
              // map vir to the new vertex
              it->second = svv;
              // the vertex of the virus could still be in the list newparents
              newparents.remove(dvv);
              // cleanup of deterministic vertex (nobody would own it otherwise...)
              delete dvv;
              /* signify that "orphans" should be removed after the loop
               * NB: these orphans could include vir.
               */
              removed_deter_virus = true;
              // keep the deterministic clone counter up-to-date
              numberOfClones--;
            } // if ok to delete (possibly last) clone
          } else if ( dvv->getValueSum() >= THRESHOLD_ILEV_ENDANGERED ) {
            // the clone is "rescued" from endangered status
            dvv->endangered = false;
            dvv->assignThreshold(THRESHOLD_ILEV_ENDANGERED); // redundant
        	}
        }
				break;
			} // case DETER_VIRUS
			default: {
				throw std::logic_error("ERROR: Vertex not of type DETER_VIRUS or STOCH_VIRUS" + RIGHT_HERE);
				break; // redundant
			}
		} // switch
	} // while
	// clear any stoch virus without parents.
	if ( removed_deter_virus ) {
		removeOrphans();
	}
	/* add the mutants of new parents.
	 * Only add mutants of parents that  have reached a certain threshold
	 * population size. NB: when a deter virus is extinct, it must be
	 * removed (std::list::remove) from newparents (see above)
	 */
	auto pit = newparents.begin();
	while ( pit != newparents.end() ) {
		DeterVirusVertex* dvv = (*pit);
		if ( dvv->getValueSum() > THRESHOLD_ILEV_PARENT ) {
			// FIXME: perhaps choose I2?, or use an incidence-based threshold
			addMutants(dvv, rng); // try to add ALL point-mutants
			pit = newparents.erase(pit);
		}	else { // keep in list, just increase the iterator
			++pit;
		}
	} // for loop over new parents (newparents)
} // updateViruses

void CiScheme::addMutants(DeterVirusVertex* pvert, Rng & rng) {
	// get the parent from the vertex
	Virus* parent = pvert->vir;
	for ( int i = 0; i < GENOME_SIZE; ++i ) {
		Virus* mutant = new Virus(*parent, i); // use the mutant constructor!
		// try to add mutant to the virus set
		auto inserted = viruses.insert(std::make_pair(mutant, nullptr));
		// use nullptr as a placehoder...
		if ( !inserted.second ) {
			// the mutant already existed, and should not be added
			delete mutant;
			/* however, the existing mutant has a new parent.
			 * When this mutant is still stochastic,
			 * the new parent contributes to the 'mutational load' of the mutant.
			 * This is handled by adding a new edge from parent->vertex
			 * to mutant->vertex
			 * When the mutant is deterministic, an undirected edge
			 * must be added to handle future extictions.
			 */
			mutant = inserted.first->first; // inserted.first points to a (Virus*, VirusVertex*) pair
			VirusVertex* mvert = inserted.first->second;
			switch ( mvert->getType() ) {
				case Vertex::STOCH_VIRUS: {
					// Graph ignores double edges.
					integrationGraph.addEdge(pvert, mvert);
					break;
				}
				case Vertex::DETER_VIRUS: {
					// Graph ignores double edges.
					integrationGraph.addUndirectedEdge(pvert, mvert);
					break;
				}
				default: {
					// ignore other neighbors
					break; // redundant...
				}
			}
		} else { // if not inserted, add a new vertex (with edges) to the graph
			/* notice that 1) we don't have to check for other parents:
			 * Had there been an alternative parent, then mutant was
			 * already added to the virus set. (NB: this does not hold for
			 * deterministic strains that have not reached the parent threshold)
			 * 2) we don't test here if mutant should be added to the
			 * system of ODEs (and by default, it is not). This is done
			 * in the next update step.
			 */
			StochVirusVertex* svv = new StochVirusVertex(mutant, &ipars);
			// in viruses, the mutant is now mapped to nullptr. Correct this:
			inserted.first->second = svv;
			svv->assignThreshold(rng.Exponential(1.0));
			// add the edge parent to mutant
			integrationGraph.addEdge(pvert, svv);
			// add the edge target to mutant
			integrationGraph.addEdge(target, svv);
			// find responses against mutant, and add edges to Graph
			for ( auto rit = responses.begin(); rit != responses.end(); ++rit ) {
				bool has_ep = rit->first->hasEpitope(*mutant);
				ResponseVertex* rvert = rit->second;
				if ( rvert->type == Vertex::ACTIVE_RESPONSE && has_ep ) {
					// the response is active and relevant for the virus
					integrationGraph.addEdge(rvert, svv);
				}
			} // for loop over responses
		} // else: new mutant added to viruses
	} // for loop over all point-mutants of parent
} // addMutants

void CiScheme::removeOrphans() {
	/* delete any stoch virus without a parent
	 * TODO: this function might benefit from information about a recently
	 * removed virus.
	 */
	auto it = viruses.begin();
	while ( it != viruses.end() ) {
		Virus* virus = it->first;
		VirusVertex* vvert = it->second;
		// only check stochastic viruses
		if ( vvert->getType() == Vertex::STOCH_VIRUS ) {
			// cast to StochVirusVertex to get access to getNumberOfParents (TODO: needed?)
			StochVirusVertex* svv = (StochVirusVertex*) vvert;
			if ( svv->getNumberOfParents() == 0 ) {
				integrationGraph.removeVertex(vvert);
				delete vvert;
				delete virus;
				it = viruses.erase(it);
			}	else {
				/* iterator was increased within if-statement
			   * by erasing the current virus
			   */
				++it;
			} // if/else no parents
		} else { // if not stochastic virus...
			/* iterator was increased within if-statement
			 * by erasing the current virus, or by operator++
			 */
			++it;
		}
	} // while-loop over viruses
} // removeOrphans



void CiScheme::updateResponses(Rng & rng) {
	/* check if responses have reached thresholds:
	 * 1) (memory) latent responses may have become active
	 * 2) active responses may become escaped memory or out-competed
	 */
	for ( auto rit = responses.begin(); rit != responses.end(); ++rit ) {
		ResponseVertex* rvert = rit->second;
		switch ( rvert->getType() ) {
			case Vertex::ACTIVE_RESPONSE: {
				ActiveResponseVertex* arv = (ActiveResponseVertex*) rvert;
          if ( !arv->endangered ) {
            if ( arv->getValue() < THRESHOLD_ILEV_ENDANGERED ) {
              double R = arv->getR();
              if ( R < 1.0 ) {
                double u = rng.Uniform();
                double tau = extinction_threshold(u, arv->getValueSum(), R);
                arv->assignThreshold(tau);
                arv->endangered = true;
              } // the reporduction number is smaller than 1: flag as endangered
          	} // the value is below the danger-threshold
            // else: all is fine...
          } else { // the response is "endangered. delete?"
          if ( arv->getValueSum() < arv->getThreshold() ) {
            /* the response is out-competed, or escaped,
             * replace the active response by a (memory)
             * latent response.
             */
          	LatentResponseVertex* lrv = arv->mkLatentResponseVertex(integrationGraph);
          	integrationGraph.removeVertex(arv);
            delete arv;
            // associate the response with the new vertex
            rit->second = lrv;
            lrv->assignThreshold(rng.Exponential(1.0));
            // keep track of statistics // TODO: soon obsolete?
            --numberOfResponses;
          } else if ( arv->getValueSum() >= THRESHOLD_ILEV_ENDANGERED ) {
            // the response is "rescued" from endangered status
            arv->endangered = false; // TODO: set/get functions
            arv->assignThreshold(THRESHOLD_ILEV_ENDANGERED); // redundant
          }
        }
				break;
			} // case active response
			case Vertex::LATENT_RESPONSE: {
				LatentResponseVertex* lrv = (LatentResponseVertex*) rvert;
				if ( !lrv->isInOdeSys() ) {
					if ( lrv->getR() > 1.0 ) { // should the vertex be integrated?
						integrationGraph.assignInOdeSys(lrv, true);
					}
				}	else { // the vertex is included in the system of ODEs
					if ( lrv->getR() <= 1.0 ) { // re-evaluate membership of ODE system
						integrationGraph.assignInOdeSys(lrv, false);
					}	else {
						// make the respose active?
						if ( lrv->getValueSum() > lrv->getThreshold() ) {
							/* the response will be activated. Any stochastic
							 * virus expressing the epitope need to be informed!
							 */
							ActiveResponseVertex* arv = lrv->mkActiveResponseVertex(integrationGraph);
							integrationGraph.removeVertex(lrv);
							delete lrv;
							arv->assignValue(1.0);
							arv->assignThreshold(THRESHOLD_ILEV_ENDANGERED); // redundant
              arv->endangered = false;
							// associate the response with the new vertex
							rit->second = arv;
							// inform stochastic mutants
							informStochVirus(arv);
							// keep track of statistics
							++numberOfResponses; // TOOD soon obsolete?
						} // if threshold reached
					} // if/else: reproduction number <= 1
				} // if/else not in ODE system
				break;
			} // case latent response
			default: {
				throw std::logic_error("ImmuneResponse::vertex not of type ACTIVE_RESPONSE or LATENT_RESPONSE" + RIGHT_HERE);
				break; // redundant
			}
		} // switch on the Vertex::type
	} // for-loop over immune responses
} // method CiScheme::updateResponses

void CiScheme::informStochVirus(ActiveResponseVertex* arv) {
	/* look for stochastic virus in viruses,
	 * and let arv point to it if the virus expresses
	 * the epitope of 'arv's ImmuneResponse 'resp'.
	 */
	ImmuneResponse* resp = arv->resp; // TODO: get function
	for ( auto it = viruses.begin(); it != viruses.end(); ++it ) {
		Virus* vir = it->first;
		VirusVertex* vvert = it->second;
		if ( vvert->getType() == Vertex::STOCH_VIRUS && resp->hasEpitope(*vir) ) {
			// deterministic virus handled elsewhere
			integrationGraph.addEdge(arv, vvert);
		} // virus must be stochastic, does the virus expres the epitope?
	} // for loop over viruses
}

void CiScheme::informLatentResponses(DeterVirusVertex* dvv) {
	Virus* vir = dvv->vir; // TODO: get function
	for ( auto rit = responses.begin(); rit != responses.end(); ++rit ) {
		ImmuneResponse* resp = rit->first;
		ResponseVertex* rvert = rit->second;
		if ( rvert->getType() == Vertex::LATENT_RESPONSE && resp->hasEpitope(*vir) ) {
			integrationGraph.addEdge(dvv, rvert);
		} // the virus expresses the right epitope
	} // loop over immune responses
}


void CiScheme::printState(std::ostream & os) const {
	os << "<odestate " // tag
	   << "t='" << getTime() << "' " // parameters
	   << ">" << std::endl;

	os << "<target "
	   << "x='" << getTargetValue() << "' " // TODO: a single xml node?
	   << "/>" << std::endl;

  os << "<quiescent "
		 << "x='" << getQuiescentValue() << "' "
	 	 << "/>" << std::endl;

	for ( auto vit = viruses.begin(); vit != viruses.end(); ++vit ) {
		Virus* vir = vit->first;
		VirusVertex* vvert = vit->second;
		if ( vvert->getType() == Vertex::DETER_VIRUS ) {
			auto it = deterVirusCache.find(vir);
			if ( it != deterVirusCache.end() ) {
				Id id = it->second;
				os << "<virus "
			     << "x1='" << vvert->getValue(0) << "' " // TODO: better attrib names? print the vertex?
					 << "x2='" << vvert->getValue(1) << "' "
			     << "id='" << id << "' "
			     << "/>" << std::endl;
			}	else {
				throw std::logic_error("requesting Id of non-cached virus" + RIGHT_HERE);
			}
		}
	}
	for ( auto rit = responses.begin(); rit != responses.end(); ++rit ) {
		ImmuneResponse* resp = rit->first;
		Vertex* rvert = rit->second;
		if ( rvert->getType() == Vertex::ACTIVE_RESPONSE ) {
			Id id = resp->getId();
			os << "<response "
			   << "x='" << rvert->getValue() << "' "
		     << "id='" << id << "' "
		     << "/>" << std::endl;
		}
	}
	// print the graph of interactions (as text)
	os << "<dot>" << std::endl;
	os << integrationGraph << std::endl;
	os << "</dot>" << std::endl;
	os << "</odestate>"; // closing tag
}


void CiScheme::printVirusCache(std::ostream & os) const {
	os << "<virus_cache >" << std::endl;
	for ( auto it = deterVirusCache.begin(); it != deterVirusCache.end(); ++it ) {
		Id id = it->second;
		Virus* vir = it->first;
		os << "<unique_virus id='" << id << "' >" << std::endl;
		os << *vir << std::endl;
		os << "</unique_virus>" << std::endl;
	}
	os << "</virus_cache>";
}


std::ostream & operator<<(std::ostream & os, const CiScheme & oscheme) {
	oscheme.printState(os); return os;
}

Id CiScheme::cacheVirus(Virus* vir) {
	// try to find vir in cachedDeterVirus
	auto it = deterVirusCache.find(vir);
	if ( it == deterVirusCache.end() ) {
		// make a copy of vir
		Virus* vir_copy = new Virus(*vir);
		deterVirusCache.emplace(vir_copy, nextId);
		return nextId++; // increase the next ID, postfit ++ returns the original value
	}	else { // the "clonotype" already exists, return the ID
		return it->second;
	}
}

std::pair<Id, bool> CiScheme::getVirusId(Virus* vir) const {
	auto it = deterVirusCache.find(vir);
	if ( it == deterVirusCache.end() ) {
		return std::make_pair(ANON_ID, false);
	}	else { // the "clonotype" already exists, return the ID
		return std::make_pair(it->second, true);
	}
}

void CiScheme::addVirusRecord(DeterVirusVertex* dvv, Virus* vir, VirusEvent event) {
	auto it = deterVirusCache.find(vir); // TODO: beter system
	if ( it == deterVirusCache.end() ) {
		throw std::logic_error("ERROR: Virus has not been cached" + RIGHT_HERE);
	}
	Virus* cached_vir = it->first;

	VirusRecord record; // TODO: constructor for VirusRecord
	record.event = event;
	record.t = getTime();
	record.r = dvv->getR(true); // boolean indicates Malthusian fitness
	record.R = dvv->getR(); // TODO: redundant?
	record.escape = dvv->isEscape();
	record.id = it->second;
	record.r_relative = 0.0;
	record.r_mimresp = 0.0;
	auto parents = dvv->getParents();
	for ( DeterVirusVertex* pvert : parents ) {
		double r_parent = pvert->getR(true); // Malthusian
		record.r_relative += record.r - r_parent;
		double rprime = dvv->getR(true, pvert); // malthusian, parents responses
		record.r_mimresp += rprime - r_parent;

		// find the parents Id
		auto id_err = getVirusId(pvert->vir); // TODO: make sure that the VirusVertes has the Id
		if ( !id_err.second ) {
			throw std::logic_error("ERROR: Virus' parent has not been cached" + RIGHT_HERE);
		}
		record.ids_parent.push_back(id_err.first);
	}
	if ( !parents.empty() ) {
		record.r_relative /= parents.size();
		record.r_mimresp /= parents.size();
	}

	std::list<VirusRecord> record_list(1, record);
	auto result = virusRecords.insert(std::make_pair(cached_vir, record_list));
	if ( !result.second ) {
		result.first->second.splice(result.first->second.end(), record_list);
	}
}

std::pair<double, int> CiScheme::getReplacementRate(bool esc,
			bool mimic_resps_parent, bool internal_nodes) const {
	// find internal nodes
	std::unordered_set<Id> ids_internal; // empty when internal_nodes is false
	if ( internal_nodes ) {
		ids_internal = getInternalNodesAncestorGraph();
	}
	// average replacement rate
	double replacement_rate = 0.0;
	int number_of_replacements = 0;
	// compute statistics
	for ( auto rit = virusRecords.begin(); rit != virusRecords.end(); ++rit ) {
		const std::list<VirusRecord> & record_list = rit->second;
    double t_birth = 0.0;
    double t_death = getTime(); // current time
    double r_relative = 0.0;
		double r_mimresp = 0.0;
    bool escape = false;
    bool right_censored = true; // becomes false whenever there is a death event
		bool internal = false;
    // there might be multiple births and deaths
		int num_birth_events = 0;
		int num_death_events = 0;
		for ( auto & record : record_list ) {
      switch ( record.event ) {
				case VIRUS_TRANSMIT: // NB: deliberate fallthrough to VIRUS_BIRTH
        case VIRUS_BIRTH: {
					num_birth_events++;
					t_birth = record.t;
					r_relative = record.r_relative;
					r_mimresp = record.r_mimresp;
					escape = record.escape;
					internal = ids_internal.find(record.id) != ids_internal.end();
					break;
        }
        case VIRUS_DEATH: {
					num_death_events++;
          t_death = record.t;
          right_censored = false;
          break;
        }
      } // switch record.event
			// TODO: there can be multiple birth and death events
			if ( num_death_events == 1 ) break; // for now, only use the first birth and death
		} // for record in record_list
    double t_present = t_death - t_birth;
		if ( t_present > VIRUS_RECORD_THRESHOLD || right_censored ) {
			if ( !internal_nodes || internal || right_censored ) { // if internal_nodes is true, then internal must be true
				if ( (esc && escape) || (!esc && !escape)) { // average over escapes or reversions
					if ( mimic_resps_parent ) {
						replacement_rate += r_mimresp;
					} else {
						replacement_rate += r_relative;
					}
					number_of_replacements++;
				}
			}
		}
	}
	if ( number_of_replacements > 0 ) {
		replacement_rate /= number_of_replacements;
	}
	return std::make_pair(replacement_rate, number_of_replacements);
}

void CiScheme::printVirusRecords(std::ostream & os) const {
	// TODO
}

void CiScheme::printAncestorGraph(std::ostream & os) const {
	os << "digraph {\n" // start graph
		 << "rankdir=LR;\n"; // landscape orientation
	for ( auto it = virusRecords.begin(); it != virusRecords.end(); ++it ) {
		auto & records = it->second; // alias
		for ( auto & rec : records ) {
			if ( rec.event == VIRUS_BIRTH || rec.event == VIRUS_TRANSMIT ) {
				os << rec.id << " -> { ";
				for ( auto & pid : rec.ids_parent ) {
					os << pid << " ";
				}
				os << "}\n";
			}
		}
	}
	os << "}";
}

void CiScheme::clearAllStochasticViruses() {
	auto it = viruses.begin();
	while ( it != viruses.end() ) {
		Virus* vir = it->first;
		VirusVertex* vvert = it->second;
		if ( vvert->getType() == Vertex::STOCH_VIRUS ) {
			integrationGraph.removeVertex(vvert);
			delete vvert;
			delete vir;
			it = viruses.erase(it);
		}	else {
			it++;
		}
	}
}

std::unordered_set<Id> CiScheme::getInternalNodesAncestorGraph() const {
	std::unordered_set<Id> idSet;
	for ( auto & pair : virusRecords ) { // Virus* -> std::list<VirusRecord>
		auto & records = pair.second; // get an alias
		for ( auto & record : records ) {
			idSet.insert(record.ids_parent.begin(), record.ids_parent.end());
		}
	}
	return idSet;
	// FIXME: avoid the situation A -> B -> A where B should count as a leaf
}
