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

#include "world.hpp"

World::World(std::string identifier_base, unsigned long seed,
			int population_size, int cohort_size, double t0) {
	rng.seed(seed);
	workerpool.initWorkerPool(NUMBER_OF_THREADS, rng);

	cumulativeIncidence = 0;
	failedInfections = 0;
	t = t0; // time
	dt_cache = 0.0; // length of the last time interval

	agentIdCounter = 0;
	lineageIdCounter = 0;
	this->identifier_base = identifier_base;
	rank = 1;
	identifier = identifier_base; // for filenames...
	this->population_size = population_size;
	this->cohort_size = cohort_size;
	birth_rate = population_size / expectedCarnesSurvival();

	openFileStreams(); // sets printing loads

	fitnessFunction.init(rng);
	Virus wildtype = Virus(&fitnessFunction, rng, Virus::EVOLVE);
	wildtype.setLineageId(lineageIdCounter++); // postfix ++
	wildtypes.emplace(wildtype.getLineageId(), wildtype);
	allMhcLoci = createImmuneResponses(rng);

	// fill population with Agents
	for ( int i = 0; i < population_size; ++i ) {
  	// sample agents according to a disease-free steady state age distribution
    Agent* agent = new Agent(t, rng, agentIdCounter++, Agent::RANDAGE);
    // postfix++ returns the old value
		population.push_back(agent);
	}
	std::cout << "created World (id: " << identifier << ")" << std::endl;
}

World::World(const World & old_world, unsigned long seed, double t0) {
	rng.seed(seed);
	workerpool.initWorkerPool(NUMBER_OF_THREADS, rng);

	cumulativeIncidence = 0;
	failedInfections = 0;
	t = t0; // time
	dt_cache = 0.0; // length of the last time interval

	agentIdCounter = 0;
	lineageIdCounter = 0;
	identifier_base = old_world.identifier_base;
	rank = old_world.rank + 1;
	identifier = identifier_base + "_" + std::to_string(rank); // for filenames...
	population_size = old_world.population_size;
	cohort_size = old_world.cohort_size;
	birth_rate = old_world.birth_rate;

	openFileStreams(); // sets 'printing loads'

	fitnessFunction = old_world.fitnessFunction;
	auto pair = old_world.sampleVirus(rng);
	if ( !pair.second ) {
		throw std::logic_error("unable so sample virus from old world" + RIGHT_HERE);
	}
	// NB: pass LOCAL fitnessFunction to wt (is a copy of the old ffun)
	Virus wildtype(&fitnessFunction, pair.first);
	wildtype.setLineageId(lineageIdCounter++); // postfix ++
	wildtypes.emplace(wildtype.getLineageId(), wildtype);

	// the new World has completely different immune responses
	allMhcLoci = createImmuneResponses(rng);

	// fill population with Agents
	for ( int i = 0; i < population_size; ++i ) {
  	// sample agents according to a disease-free steady state age distribution
    Agent* agent = new Agent(t, rng, agentIdCounter++, Agent::RANDAGE);
    // postfix++ returns the old value
		population.push_back(agent);
	}
	std::cout << "created World (id: " << identifier << ")" << std::endl;
}


World::~World() {
	workerpool.waitForWorkerPoolToExit(); // shuts down worker threads
	closeFileStreams(); // finishes off some logging
	// delete living population
	clearPopulation();
	/* by deleting the living population,
	 * also the ancestoral lineage will be deleted,
	 * except for the MRCAs.
	 */
	for ( auto it = ancestors.begin(); it != ancestors.end(); ++it ) {
		Host* ancestor = (*it);
		delete ancestor;
	}
}

void World::run(double tmax, bool keep_alive) {
	int steps = (tmax-t)/PLEV_TIMESTEP;
	EtaEstimator eta(steps, 0.99);
	/* eta must be initialized with the total number of steps,
	 * and a 'relative weight' of historical steps for computing ETA
	 */
	static bool printed_header(false);
	if ( !printed_header ) {
		std::cout << std::setw(10) << "ETA" << " "
		          << std::setw(8) << "t (y)" << " "
							<< std::setw(8) << "N" << " "
							<< std::setw(8) << "I" << " "
							<< std::setw(8) << "dI (1/y)" << " "
							<< std::endl;
		printed_header = true;
	}
	// start loop
  while ( t < tmax ) {
		if ( keep_alive && livingHosts.empty() ) {
			std::cout << "\t[(re)seeding epidemic...]" << std::endl;
			infectRandomAgents(INITIAL_PREVALENCE);
		}
    updateState(PLEV_TIMESTEP);
    updateAncestors();
    printStatistics();
    eta.update();
    // print some info
    std::cout << "\r" << TERM_CLEAR_LINE
							<< std::setw(10) << eta << " "
              << std::setw(8) << t << " "
              << std::setw(8) << getPopulationSize() << " "
              << std::setw(8) << livingHosts.size() << " " // prevalence
              << std::setw(8) << recentInfectionTimes.size() << " " // incidence
							<< std::flush;
  } // while t < tmax
	workerpool.syncWorkerThreads();
}

void World::updateState(double dt) {
  // randomize the order of the population
  population.sort(compareByRandomIndex);
  // sample a number of newborns
	int birth_incidence = rng.Poisson(birth_rate * dt);
  // number of Agents infected during this timestep
  int incidence = 0;
  // walk through population, and update hosts' internal state
  auto it = population.begin();
  while ( it != population.end() ) {
    Agent* agent = (*it);
		if ( agent == nullptr ) {
			throw std::logic_error("agent in population is NULL" + RIGHT_HERE);
		}
		// update internal state. NB: important to do this first
		agent->updateState(dt, rng);
		if ( !agent->isAlive() ) { // dead agent!
			if ( agent->isInfected() ) { // remove the Host pointer from the hosts set
				livingHosts.erase(agent->getHostPtr());
			}
			// delete the agent (removes any contacts with the agent)
			delete agent;
			// remove iterator from population vector, and let it point to the next agent
			it = population.erase(it);
		}	else { // the Agent is alive after updating the internal state
			/* if the Agent is susceptible, and has accumulated enough
			 * infection "hazard", sample a TF virus from its partners.
			 * NB: at this point, all (infected) partners should still
			 * be "around" to sample from.
       * Notice that we first update the state of the agent,
       * hence dead neighbors will already be dead before updating
       * (and hence before the infectionLoad reaches the threshold)
			 */
			if ( !agent->isInfected() && agent->requiresInfection() ) {
				bool ok = agent->infectByRandomPartner(rng, allMhcLoci);
				// ok might be false when e.g. i-level R0 <= 1.0
				if ( ok ) {
					auto pr = livingHosts.insert(agent->getHostPtr());
					if ( !pr.second ) {
						throw std::logic_error("adding duplicate Host pointer to livingHostSet" + RIGHT_HERE);
					}
					workerpool.addNewJob(agent->getHostPtr());
					incidence++;
					cumulativeIncidence++;
				}	else { // Transmission failed, re-sample threshold
					failedInfections++;
				}
				agent->assignInfectionThreshold(rng);
				/* re-assigning the infection threshold is only needed
				 * in cases when the infection was not succesful.
				 * In a future version, we could allow for co-infection!
				 */
			}
			it++; // let iterator it point to the next agent...
		} // else: Agent is alive...
  } // population for loop
	// add newborns
	for ( int i = 0; i < birth_incidence; ++i ) {
		Agent* agent = new Agent(t, rng, agentIdCounter++, Agent::NEWBORN);
		agent->updateState(dt, rng);
		population.push_back(agent);
	}
	/* by updating the Agent's states, some will have reached
	 * the contact formation threshold.
	 * These Agents will try to form new contacts.
	 */
  makeContacts(t, dt, rng);
	// for some of the contacts, set the transmission_route flag true
	makeTransmissionRoutes(rng);
  // increment time
  t += dt;
  // save length last time step for logging
  dt_cache = dt;
	// update incidence list to get windowed incidence
	updateRecentInfectionTimes(incidence);
} // update state


void World::infectRandomAgents(int k) {
	shufflePopulation(rng); // make sure the population is in random order
	// reverse-sort by TMR degree
	auto pred = [](Agent* la, Agent* ra){return !compareAgentByTmrDegree(la, ra);};
	/* infect most connected agents to jump-start the epidemic
	 * std::list::sort preserves the relative order of identical objects in the list
	 */
	population.sort(pred);
	int n = 0; // will be the number of succesful infections
	for ( auto it = population.begin(); it != population.end() && n < k; ++it ) {
		Agent* agent = (*it);
		if ( !agent->isInfected() ) {
			bool ok = false;
			if ( MULTIPLE_PLEV_FOUNDERS || wildtypes.empty() ) {
				Virus wildtype(&fitnessFunction, rng, Virus::EVOLVE);
				wildtype.setLineageId(lineageIdCounter++); // postfix ++
				// make sure the new wildtype is stored in wildtypes
				wildtypes.emplace(wildtype.getLineageId(), wildtype);
				Virus mutant(wildtype, rng, int(NUMBER_INITIAL_MUTATIONS)); // make sure to pass an int
				// use fVir to infect a new host
				ok = agent->infect(rng, allMhcLoci, mutant);
			} else { // use the same founder virus each time
				// wildtypes.begin()->second should be the only wildtype virus in the population
				Virus mutant(wildtypes.begin()->second, rng, int(NUMBER_INITIAL_MUTATIONS));
				ok = agent->infect(rng, allMhcLoci, mutant);
			}
			// update the ancester trace, add to Job list, add Host to livingHosts
			if ( ok ) {
				Host* host = agent->getHostPtr();
				auto pr = livingHosts.insert(host);
				if ( !pr.second ) {
					throw std::logic_error("adding duplicate Host pointer to livingHostSet" + RIGHT_HERE);
				}
				host->setMrcaFlag(true); // make the host an ancestor
				ancestors.push_back(host);
				workerpool.addNewJob(host); // start the i-level simulation
				n++; // increase local infection counter
				cumulativeIncidence++;
			} else {
				failedInfections++;
			} // if, else ok
		} // if agent is not infected
	} // for loop over the population
}


// World members to handle statistics and logging

void World::openFileStreams() {
	// variables used for logging at different intervals
	plstat = 0.0;
	plherit = 0.0;
	pllong = 0.0;
	pllonghd = 0.0;
	plcross = 0.0;
	plseq = 0.0;
	plcost = 0.0;
	plnet = 0.0;

	// open stats file
	ofstat.open((DATA_FOLDER + "stats-file-" + identifier + ".xml").c_str());
	if ( ofstat ) {
		ofstat << "<epidemic identifier='" << identifier << "' >" << std::endl;
		ofstat << xmlStringParameters() << std::endl;
	}	else {
		throw std::runtime_error("cannot write to stats-file" + RIGHT_HERE);
	}

	// open sem data file
	ofherit.open((DATA_FOLDER + "herit-file-" + identifier + ".xml").c_str());
	if ( ofherit ) {
		ofherit << "<herit_datas identifier='" << identifier << "' >" << std::endl;
		ofherit << xmlStringParameters() << std::endl;
	}	else {
		throw std::runtime_error("cannot write to herit-file" + RIGHT_HERE);
	}

	// open the ancestor trace
	ofancs.open((DATA_FOLDER + "ancs-file-" + identifier + ".xml").c_str());
	if ( ofancs ) {
		ofancs << "<ancestor_trace identifier='" << identifier << "' >" << std::endl; // TODO: make some sort of tag-function?
		ofancs << xmlStringParameters() << std::endl;
	}	else {
		throw std::runtime_error("cannot write to ancestor-file" + RIGHT_HERE);
	}

	// open file for longitudinal samples
	oflong.open((DATA_FOLDER + "longitudinal-sample-file-" + identifier + ".xml").c_str());
	if ( oflong ) {
		oflong << "<longitudinal_sample identifier='" << identifier << "' >" << std::endl;
		oflong << xmlStringParameters() << std::endl;
	}	else {
		throw std::runtime_error("cannot write to longitudinal-sample-file" + RIGHT_HERE);
	}

	// open file for longitudinal HD samples
	oflonghd.open((DATA_FOLDER + "longitudinal-hd-sample-file-" + identifier + ".xml").c_str());
	if ( oflonghd ) {
		oflonghd << "<longitudinal_hd_sample identifier='" << identifier << "' >" << std::endl;
		oflonghd << xmlStringParameters() << std::endl;
	}	else {
		throw std::runtime_error("cannot write to longitudinal-hd-sample-file" + RIGHT_HERE);
	}

	// open file for cross-sectional samples
	ofcross.open((DATA_FOLDER + "crosssectional-sample-file-" + identifier + ".xml").c_str());
	if ( ofcross ) {
		ofcross << "<crosssectional_samples identifier='" << identifier << "' >" << std::endl;
		ofcross << xmlStringParameters() << std::endl;
	}	else {
		throw std::runtime_error("cannot write to crosssectional-sample-file" + RIGHT_HERE);
	}

  // open a file for the fitness function (.xml)
  offun.open((DATA_FOLDER + "fitnessfunction-file-" + identifier + ".xml").c_str());
  if ( offun ) {
		// OK! (write stuff to file when closing...)
  } else {
		throw std::runtime_error("cannot write to fitnessfunction-file (.xml)" + RIGHT_HERE);
  }

  // open a file for the fitness function (.dot)
  offundot.open((DATA_FOLDER + "fitnessfunction-file-" + identifier + ".gv").c_str());
  if ( offundot ) {
		// OK! (write stuff to file when closing...)
  } else {
		throw std::runtime_error("cannot write to fitnessfunction-file (.gv)" + RIGHT_HERE);
  }

  // open a file for p-level immune responses
  ofimmresp.open((DATA_FOLDER + "immresp-file-" + identifier + ".xml").c_str());;
  if ( ofimmresp ) {
      // OK! (write stuff to file when closing...)
  } else {
		throw std::runtime_error("cannot write to immuneresponse-file" + RIGHT_HERE);
  }

  ofseq.open((DATA_FOLDER + "sequence-file-" + identifier + ".xml").c_str());
  if ( ofseq ) {
    ofseq << "<sequence_samples identifier='" << identifier << "' >" << std::endl;
  } else {
		throw std::runtime_error("cannot write to sequence-file" + RIGHT_HERE);
  }

	ofcost.open((DATA_FOLDER + "fitnesscost-file-" + identifier + ".xml").c_str());
	if ( ofcost ) {
		ofcost << "<fitnesscost_samples identifier='" << identifier << "' >" << std::endl;
	} else {
		throw std::runtime_error("cannot write to fitnesscost-file" + RIGHT_HERE);
  }

	ofnet.open((DATA_FOLDER + "contact-network-file-" + identifier + ".xml").c_str());
	if ( ofnet ) {
		ofnet << "<contact_networks identifier='" << identifier << "' >" << std::endl;
	} else {
		throw std::runtime_error("cannot write to contact-network-file" + RIGHT_HERE);
	}

	ofassoc.open((DATA_FOLDER + "associations-file-" + identifier + ".xml").c_str());
	if ( ofassoc ) {
		ofassoc << "<associations identifier='" << identifier << "' >" << std::endl;
	} else {
		throw std::runtime_error("cannot write to associations-file" + RIGHT_HERE);
	}

	ofaids.open((DATA_FOLDER + "aids-hazard-file-" + identifier + ".xml").c_str());
	if ( ofnet ) {
		ofaids << "<aids_hazard identifier='" << identifier << "' >" << std::endl;
	} else {
		throw std::runtime_error("cannot write to aids-hazard-file" + RIGHT_HERE);
	}
} // openFileStreams


void World::printStatistics() {
	// Agent-like data
	plnet += dt_cache;
	if ( ofnet ) {
		if ( plnet >= CONTACT_NETWORK_SAMPLE_INTERVAL ) {
			plnet -= CONTACT_NETWORK_SAMPLE_INTERVAL;
			ofnet << "<contact_network t='" << t << "' >" << std::endl;
			printGraph(ofnet); // TODO: DOT is not xml-escaped
			ofnet << std::endl << "</contact_network>" << std::endl;
		}
	}
	// Host-like data
	if ( livingHosts.empty() ) return;
	// increment printing loads
	plstat += dt_cache;
	plherit += dt_cache;
	pllong += dt_cache;
	pllonghd += dt_cache;
	plcross += dt_cache;
	plseq += dt_cache;
	plcost += dt_cache;
	plassoc += dt_cache;
	plaids += dt_cache;

	if ( ofstat ) {
		if ( plstat >= STATS_SAMPLE_INTERVAL ) {
			plstat -= STATS_SAMPLE_INTERVAL; // reset the plstat load
			printXmlGeneralStats(ofstat);
			ofstat << std::endl;
		}
	} else { // if ofstat
		throw std::runtime_error("cannot write to stats-file" + RIGHT_HERE);
	}

	if ( ofherit ) {
		if ( livingHosts.size() > INITIAL_PREVALENCE && plherit >= HERIT_SAMPLE_INTERVAL ) {
			plherit -= HERIT_SAMPLE_INTERVAL;
			printXmlHeritDataNode(ofherit);
			ofherit << std::endl;
		}
	} else {
		throw std::runtime_error("cannot write to herit-data-file" + RIGHT_HERE);
	}

	/* longitudinal samples. The sampling rate does depend on the
	 * population size. TODO: keep this rate constant?
	 */
	if ( oflong ) {
		if ( pllong >= LONGITUDINAL_SAMPLE_INTERVAL ) {
			pllong -= LONGITUDINAL_SAMPLE_INTERVAL;
			// take a random sample
			int j = rng.Integer(population.size());
			auto it = population.begin();
			std::advance(it, j);
			Agent* agent = (*it);
			Host* host = agent->getHostPtr();
			if ( host != nullptr ) {
				oflong << (*host) << std::endl;
			}
		}
	} else { // if oflong
		throw std::runtime_error("cannot write to longitudinal-sample-file" + RIGHT_HERE);
	}

   	/* HD longidudinal samples. The sampling rate does depend on the
	 * population size. TODO: keep this rate constant?
	 */
	if ( oflonghd ) {
		if ( pllonghd >= LONGITUDINAL_HD_SAMPLE_INTERVAL ) {
			pllonghd -= LONGITUDINAL_HD_SAMPLE_INTERVAL;
			// Sample a random Host from livingHosts
			int idx = rng.Integer(livingHosts.size());
			auto it = livingHosts.begin();
			std::advance(it, idx);
			Host* host = *it;
			if ( host == nullptr ) {
				throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
			}
			std::cout << "\t[running HD host simulation...]" << std::flush;
			// get the current virus
			Virus virus = host->getVirus(rng);
			IndPars ipars = host->getIndPars();
			// create a new host infected with the sampled virus TODO: same immune responses?
			Host host_hd(rng, allMhcLoci, virus, ipars, t);
			host_hd.run(rng, true); // pass true for detailed logging...
			oflonghd << host_hd << std::endl;
		} // time to sample...
	} else { // if oflonghd
		throw std::runtime_error("cannot write to longitudinal-hd-sample-file" + RIGHT_HERE);
	}

	if ( ofcross ) {
		if ( plcross >= CROSSSECTIONAL_SAMPLE_INTERVAL ) {
			plcross -= CROSSSECTIONAL_SAMPLE_INTERVAL;
			shufflePopulation(rng); // make sure that we get a random sample
			ofcross << "<crosssectional_sample t='" << t << "' >" << std::endl; // opening tag
			int j = 0; // number of samples taken
			for ( auto it = population.begin(); it != population.end(); ++it ) {
				Agent* agent = (*it);
				Host* host = agent->getHostPtr();
				if ( host != nullptr ) {
					ofcross << (*host) << std::endl;
					// make sure to have enough samples by only counting infected Agents
					j++;
					if ( j >= cohort_size ) break;
				}
			}
			ofcross << "</crosssectional_sample>" << std::endl; // closing tag
		} // if it is time for another sample
	} else {
		throw std::runtime_error("cannot write to crosssectional-sample-file" + RIGHT_HERE);
	} // if, else ofcross

  if ( ofseq ) {
    if ( plseq >= SEQUENCE_SAMPLE_INTERVAL ) {
      plseq -= SEQUENCE_SAMPLE_INTERVAL;
			// TODO: make a method
			if ( !livingHosts.empty() ) {
				// take a random sample
				int i = rng.Integer(livingHosts.size());
				auto it = livingHosts.begin();
				std::advance(it, i);
				Host* host = *it;
				if ( host == nullptr ) {
					throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
				} // else
				Virus virus = host->getVirus(rng);
				ofseq << "<random_virus t='" << t << "' >" << std::endl;
				ofseq << virus << std::endl;
				ofseq << "</random_virus>" << std::endl;
				// compute the consensus
				auto pair = getConsensus();
				ofseq << "<consensus_virus "
				      << "t='" << t << "' "
							<< "consensus_defined='" << std::boolalpha << pair.second << std::noboolalpha << "' "
							<< ">" << std::endl;
				ofseq << pair.first << std::endl;
				ofseq << "</consensus_virus>" << std::endl;
			} // if there are living hosts...
    } // if it is time for another sample
  } else {
		throw std::runtime_error("cannot write to sequence-file" + RIGHT_HERE);
  }

	if ( ofcost ) {
		if ( plcost >= FITNESSCOST_SAMPLE_INTERVAL ) {
			plcost -= FITNESSCOST_SAMPLE_INTERVAL;
			ofcost << xmlNodeFitnessCosts() << std::endl;
		}
	} else {
		throw std::runtime_error("cannot write to fitnesscost-file" + RIGHT_HERE);
	}

	if ( ofassoc ) {
		if ( plassoc >= ASSOCIATION_SAMPLE_INTERVAL ) {
			// TODO: make auxiliary function for printing...
			plassoc -= ASSOCIATION_SAMPLE_INTERVAL;
			// MHC-polymorphism associations using FET
			ofassoc << "<mhc_poly_associations_fet t='" << t << "' >" << std::endl;
			auto new_mpas_fet = findAssociationsMhcPolyFET();
			for ( auto key_val : new_mpas_fet ) {
				auto key = key_val.first;
				Association mpa = key_val.second;
				int mhc_loc = key.first.first;
				int mhc_allele = key.first.second;
				int virus_loc = key.second;
				ofassoc << "<association "
				        << "locus='" << virus_loc << "' >" << std::endl;
				ofassoc << "<mhc_id "
							  << "locus='" << mhc_loc << "' "
								<< "allele='" << mhc_allele << "' "
								<< "/>" << std::endl; // TODO: there should be a function to print the MhcID xml tag
				ofassoc << mpa << std::endl;
				ofassoc << "</association>" << std::endl;
			}
			ofassoc << "</mhc_poly_associations_fet>" << std::endl;
			// un-corrected MHC-polymorphism associations using LRT
			ofassoc << "<mhc_poly_associations_lrt t='" << t << "' >" << std::endl;
			auto new_mpas_lrt = findAssociationsMhcPolyLRT(false); // false means: DO NOT correct for TF virus
			for ( auto key_val : new_mpas_lrt ) {
				auto key = key_val.first;
				Association mpa = key_val.second;
				int mhc_loc = key.first.first;
				int mhc_allele = key.first.second;
				int virus_loc = key.second;
				ofassoc << "<association "
								<< "locus='" << virus_loc << "' >" << std::endl;
				ofassoc << "<mhc_id "
								<< "locus='" << mhc_loc << "' "
								<< "allele='" << mhc_allele << "' "
								<< "/>" << std::endl; // TODO: there should be a function to print the MhcID xml tag
				ofassoc << mpa << std::endl;
				ofassoc << "</association>" << std::endl;
			}
			ofassoc << "</mhc_poly_associations_lrt>" << std::endl;
			// corrected MHC-polymorphism associations
			ofassoc << "<mhc_poly_associations_lrt_corrected t='" << t << "' >" << std::endl;
			auto new_mpas_clrt = findAssociationsMhcPolyLRT(true); // true means: correct for TF virus
			for ( auto key_val : new_mpas_clrt ) {
				auto key = key_val.first;
				Association mpa = key_val.second;
				int mhc_loc = key.first.first;
				int mhc_allele = key.first.second;
				int virus_loc = key.second;
				ofassoc << "<association "
								<< "locus='" << virus_loc << "' >" << std::endl;
				ofassoc << "<mhc_id "
								<< "locus='" << mhc_loc << "' "
								<< "allele='" << mhc_allele << "' "
								<< "/>" << std::endl; // TODO: there should be a function to print the MhcID xml tag
				ofassoc << mpa << std::endl;
				ofassoc << "</association>" << std::endl;
			}
			ofassoc << "</mhc_poly_associations_lrt_corrected>" << std::endl;
			// MHC-VL associations
			ofassoc << "<mhc_vl_associations t='" << t << "' >" << std::endl;
			auto new_mvas = findAssociationsMhcVl();
			for ( auto key_val : new_mvas ) {
				auto key = key_val.first;
				Association mva = key_val.second;
				int mhc_loc = key.first;
				int mhc_allele = key.second;
				ofassoc << "<association >" << std::endl;
				ofassoc << "<mhc_id "
								<< "locus='" << mhc_loc << "' "
								<< "allele='" << mhc_allele << "' "
								<< "/>" << std::endl; // TODO: there should be a function to print the MhcID xml tag
				ofassoc << mva << std::endl;
				ofassoc << "</association>" << std::endl;
			}
			ofassoc << "</mhc_vl_associations>" << std::endl;
		}
	} else {
		throw std::runtime_error("cannot write to associations-file" + RIGHT_HERE);
	}

	if ( ofaids ) {
		if ( plaids >= AIDS_HAZARD_SAMPLE_INTERVAL ) {
			plaids -= AIDS_HAZARD_SAMPLE_INTERVAL;
			ofaids << xmlNodeAidsHazard() << std::endl;
		}
	} else {
		throw std::runtime_error("cannot write to aids-hazard-file" + RIGHT_HERE);
	}
}



void World::closeFileStreams() {
	// write closing tags and close fstreams
	if ( ofstat ) {
		ofstat << "</epidemic>";
		ofstat.close();
	}	else {
		throw std::runtime_error("cannot write to stats-file" + RIGHT_HERE);
	}

	if ( ofherit ) {
		ofherit << "</herit_datas>";
		ofherit.close();
	}	else {
		throw std::runtime_error("cannot write to herit-data-file" + RIGHT_HERE);
	}

	if ( oflong ) {
		oflong << "</longitudinal_sample>";
		oflong.close();
	}	else {
		throw std::runtime_error("cannot write to longitudinal-sample-file" + RIGHT_HERE);
	}

	if ( oflonghd ) {
		oflonghd << "</longitudinal_hd_sample>";
		oflonghd.close();
	}	else {
		throw std::runtime_error("cannot write to longitudinal-hd-sample-file" + RIGHT_HERE);
	}


	if ( ofcross ) {
		ofcross << "</crosssectional_samples>";
		ofcross.close();
	}	else {
		throw std::runtime_error("cannot write to crosssectional-sample-file" + RIGHT_HERE);
	}

	// finish ancestor trace
	if ( ofancs ) {
		for ( auto it = ancestors.begin(); it != ancestors.end(); ++it ) {
			Host* ancestor = (*it);
			ancestor->printTree(ofancs);
		}
		ofancs << "</ancestor_trace>";
		ofancs.close();
	}	else {
		throw std::runtime_error("cannot write to ancestor-trace-file" + RIGHT_HERE);
	}

	// write fitness function to file...
	if ( offun ) {
		fitnessFunction.print(offun);
    offun.close();
  } else {
		throw std::runtime_error("cannot write to fitnessfunction-file (.xml)" + RIGHT_HERE);
	}

	// write fitness function to file...
	if ( offundot ) {
		fitnessFunction.printGraph(offundot);
    offun.close();
  } else {
		throw std::runtime_error("cannot write to fitnessfunction-file (.dot)" + RIGHT_HERE);
	}

	// write immuneresponses to file...
	if ( ofimmresp ) {
		ofimmresp << xmlStringMhcLoci(allMhcLoci, identifier);
    ofimmresp.close();
  } else {
		throw std::runtime_error("cannot write to immuneresponses-file" + RIGHT_HERE);
	}

  if ( ofseq ) {
    ofseq << "</sequence_samples>";
    ofseq.close();
  } else {
		throw std::runtime_error("cannot write to sequence-file" + RIGHT_HERE);
  }

	if ( ofcost ) {
		ofcost << "</fitnesscost_samples>";
		ofcost.close();
	} else {
		throw std::runtime_error("cannot write to fitnesscost-file" + RIGHT_HERE);
	}

	if ( ofnet ) {
		ofnet << "</contact_networks>";
		ofnet.close();
	} else {
		throw std::runtime_error("cannot write to contact-network-file" + RIGHT_HERE);
	}

	if ( ofassoc ) {
		ofassoc << "</associations>";
		ofassoc.close();
	} else {
		throw std::runtime_error("cannot write to associations-file" + RIGHT_HERE);
	}

	if ( ofaids ) {
		ofaids << "</aids_hazard>";
		ofaids.close();
	} else {
		throw std::runtime_error("cannot write to aids-hazard-file" + RIGHT_HERE);
	}
}

int World::getAncestorTreeSize() const {
  int result = 0;
  for ( auto it = ancestors.begin(); it != ancestors.end(); ++it ) {
    Host* anc = (*it);
    if ( anc == nullptr ) {
			throw std::logic_error("one of the ancestors is NULL" + RIGHT_HERE);
		}
    result += anc->getTreeSize();
  }
  return result;
}

void World::getPlevEpitopesForAllele(const MhcAllele & allele,
			BasicStats & nStats, BasicStats & dwnStats) const {
	for ( auto host : livingHosts ) {
		if ( host == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		int n = 0; double dwn = 0.0;
    allele.countEpitopes(host->getVirus(), n, dwn);
    nStats.addPoint(n);
		dwnStats.addPoint(dwn);
  }
	nStats.computeStats();
	dwnStats.computeStats();
}


std::pair<double, bool> World::getGeneticDiversity() const {
	WARN_DEPRECATED_FUN
	/* compare each genome with every other genome
	 * TODO: use a better measure: transform distance matrix into a covariance
	 * matrix, and take the square-root of the trace (total variance)
	 * use from GSL: int gsl_linalg_cholesky_decomp1(gsl_matrix* )
	 */
	long int H = 0;
	// long int to prevent overflow (TODO: double check!)
	long int n = livingHosts.size();
	long int pairs = n * (n-1) / 2;
	if ( pairs == 0 ) {
		return std::make_pair(0.0, false);
	}
	for ( auto it1 = livingHosts.begin(); it1 != livingHosts.end(); ++it1 ) {
		Host* host1 = (*it1);
		if ( host1 == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		for ( auto it2 = livingHosts.begin(); it2 != it1; ++it2 ) {
			Host* host2 = (*it2);
			if ( host2 == nullptr ) {
				throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
			}
			// get (dominant) viruses from both hosts and compute hamming distance
			H += hammingDistance(host1->getVirus(), host2->getVirus());
		}
	}
	double meanHD = double(H) / (pairs * GENOME_SIZE); // beware of overflow!
	return std::make_pair(meanHD, true);
}

std::vector<double> World::getAlleleFreqPerLocus(bool allele) const {
	// compute diversity at the loci, using the current dominant strains
	std::vector<double> allele_freqs(GENOME_SIZE, 0.0);
	int n = livingHosts.size();
	if ( n == 0 ) {
		return allele_freqs;
	} // else...
	for ( Host* host : livingHosts ) {
		if ( host == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		Virus virus = host->getVirus(); // dominant strain
		for ( int i = 0; i < GENOME_SIZE; ++i ) {
			allele_freqs[i] += ( virus[i] == allele ? 1 : 0 );
		}
	}
	for ( int i = 0; i < GENOME_SIZE; ++i ) {
		allele_freqs[i] /= n;
	}
	return allele_freqs;
}

std::vector<double> World::getGeneticDiversityPerLocus() const {
	// compute diversity at the loci, using the current dominant strains
	auto allele_freqs = getAlleleFreqPerLocus(true); // count ones
	std::vector<double> diversity(GENOME_SIZE, 0.0);
	for ( int i = 0; i < GENOME_SIZE; ++i ) {
		double x = allele_freqs[i];
		if ( 0.0 < x && x < 1.0 ) {
			diversity[i] = -(x*log(x) + (1-x)*log(1-x)) / M_LN2;
		}
	}
	return diversity;
}

std::map<AlleleLocusPair, Association> World::findAssociationsMhcPolyFET() {
	// make a list of contingency tables
	std::map<AlleleLocusPair, ContingencyTable> cTables;
	for ( auto & mhcLocus : allMhcLoci ) {
		for ( auto mhcAllele : mhcLocus ) {
			MhcID mhcId = mhcAllele.getMhcID();
			for ( int ell = 0; ell < GENOME_SIZE; ++ell ) {
				AlleleLocusPair key = std::make_pair(mhcId, ell);
				ContingencyTable ctab;
				for ( Host* host : livingHosts ) {
					if ( host == nullptr ) {
						throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
					}
					/* check that the allele is different from the wildtype allele
					 * find the right wildtype in the map wildtypes
					 */
					auto wtit = wildtypes.find(host->getVirus().getLineageId());
					if ( wtit == wildtypes.end() ) {
						throw std::logic_error("linageId not a key in wildtypes" + RIGHT_HERE);
					} // else...
					bool virAllele = host->getVirus()[ell] != wtit->second[ell];
					// getVirus() gets the most common strain, returns a const reference
					bool hasMhc = host->hasMhcAllele(mhcId);
					// add values to the contingency table
					ctab.a += ( hasMhc && virAllele ? 1 : 0 );
					ctab.b += ( hasMhc && !virAllele ? 1 : 0 );
					ctab.c += ( !hasMhc && virAllele ? 1 : 0 );
					ctab.d += ( !hasMhc && !virAllele ? 1 : 0 );
				} // loop over population
				cTables[key] = ctab;
			} // loop over virus genome
		} // loop over MHC alleles
	} // loop over MHC loci
	// compute p-values
	std::map<AlleleLocusPair, double> pValues;
	for ( auto & pair : cTables ) {
		auto & key = pair.first;
		auto & ctab = pair.second;
		// test for polymorphism
		auto pval_ok = ctab.pFisherExactTest(BOTH_TAILS);
		if ( pval_ok.second && pval_ok.first < 1.0 ) { // FIXME: use a better test!
			pValues[key] = pval_ok.first;
		}
	}
	// compute q-values
	auto qValues = computeQValues(pValues);
	// filter out significant values, make timeseries
	std::map<AlleleLocusPair, Association> new_associations; // return value
	for ( auto & pair : qValues ) {
		AlleleLocusPair key = pair.first;
		double q = pair.second;
		double p = pValues.at(key);
		ContingencyTable ctab = cTables.at(key);
		auto phi_ok = ctab.phiCoeff();
		if ( q < HIGH_FDR_THRESHOLD && p < HIGH_SIGNIF_THRESHOLD && phi_ok.second ) {
			Association mpa(t, p, q, phi_ok.first);
			auto it_ok = mhcPolyAssociationsFET.emplace(key, std::list<Association>());
			it_ok.first->second.push_back(mpa);
			new_associations.emplace(key, mpa); // add to return value
		}
	}
	return new_associations;
}

std::map<AlleleLocusPair, Association> World::findAssociationsMhcPolyLRT(bool corrected) {
	std::map<AlleleLocusPair, double> pValues;
	std::map<AlleleLocusPair, double> logOddss;
	for ( auto & mhcLocus : allMhcLoci ) {
		for ( auto mhcAllele : mhcLocus ) {
			MhcID mhcId = mhcAllele.getMhcID();
			for ( int ell = 0; ell < GENOME_SIZE; ++ell ) {
				AlleleLocusPair key = std::make_pair(mhcId, ell);
				DataListAndModel dm(1); // 1 means single allele
				for ( Host* host : livingHosts ) {
					if ( host == nullptr ) {
						throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
					}
					// check that the allele is different from the wildtype allele
					DataPoint dp(1); // 1 means single allele
					/* check that the allele is different from the wildtype allele
					 * find the right wildtype in the map wildtypes
					 */
					auto wtit = wildtypes.find(host->getVirus().getLineageId());
					if ( wtit == wildtypes.end() ) {
						throw std::logic_error("linageId not a key in wildtypes" + RIGHT_HERE);
					} // else...
					dp.observed_allele = host->getVirus()[ell] != wtit->second[ell];
					if ( corrected ) {
						dp.transmitted_allele = host->getVirusTF()[ell] != wtit->second[ell];
					} else { // assume that the host was infected with the wildtype allele
						dp.transmitted_allele = false; // false means not different from wildtype
					}
					// getVirus() gets the most common strain, returns a const reference
					dp.mhc_alleles[0] = host->hasMhcAllele(mhcId);
					// add values to the data list
					dm.data.push_back(dp);
				} // loop over population
				// select the full model
				if ( sufficientObservarions(dm) ) { // ignore really rare cases (of the full model)
					// compute log-odds using maximum likelihood
					FitResult alt_result = fitLogisticModel(dm);
					// now do NOT use MHC allele as a predictor (null model)
					dm.model[0] = false;
					FitResult null_result = fitLogisticModel(dm);
					// null and alternative differ by 1 parameter
					double pval = pLikelihoodRatioTest(null_result.loglike, alt_result.loglike, 1);
					double lo = alt_result.betas[0];
					pValues.emplace(key, pval);
					logOddss.emplace(key, lo);
				} // if sufficient observations
			} // loop over virus genome
		} // loop over MHC alleles
	} // loop over MHC loci
	// compute q-values
	auto qValues = computeQValues(pValues);
	// filter out significant values, make timeseries
	std::map<AlleleLocusPair, Association> new_associations; // return value
	for ( auto & pair : qValues ) {
		AlleleLocusPair key = pair.first;
		double q = pair.second;
		double p = pValues.at(key);
		double lo = logOddss.at(key);
		if ( q < HIGH_FDR_THRESHOLD && p < HIGH_SIGNIF_THRESHOLD ) {
			Association mpa(t, p, q, lo);
			if ( corrected ) {
				auto it_ok = mhcPolyAssociationsLRTcorrected.emplace(key, std::list<Association>());
				it_ok.first->second.push_back(mpa);
			} else {
				auto it_ok = mhcPolyAssociationsLRT.emplace(key, std::list<Association>());
				it_ok.first->second.push_back(mpa);
			}
			new_associations.emplace(key, mpa); // add to return value
		}
	}
	return new_associations;
}

std::map<MhcID, Association> World::findAssociationsMhcVl() {
	std::map<MhcID, double> pValues;
	std::map<MhcID, double> Ustats; // U-statistics (probabilities)
	for ( auto & mhcLocus : allMhcLoci ) {
		for ( auto mhcAllele : mhcLocus ) {
			MhcID mhcID = mhcAllele.getMhcID();
			std::list<double> positives; // has the allele
			std::list<double> negatives; // does not have the allele
			for ( Host* host : livingHosts ) {
				if ( host == nullptr ) {
					throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
				}
				// add VL of host to positives or negatives
				if ( host->hasMhcAllele(mhcID) ) {
					positives.push_back(host->getLogSpvl());
				} else {
					negatives.push_back(host->getLogSpvl());
				}
			} // loop over living Hosts
			// compute test statistic and p-value for Pr(VL_positive > VL_negative)
			if ( !positives.empty() && !negatives.empty() ) {
				// test if non-empty, otherwise pMannWhitneyUTest would throw an exception...
				auto p_U = pMannWhitneyUTest(positives, negatives, BOTH_TAILS);
				// collect p-values and U-statistic
				pValues.emplace(mhcID, p_U.first);
				Ustats.emplace(mhcID, p_U.second);
			}
		} // loop over MHC alleles
	} // loop over MHC loci
	// tramsform pValues into qValues
	auto qValues = computeQValues(pValues);
	std::map<MhcID, Association> new_associations; // return value
	for ( auto & id_q : qValues ) {
		MhcID key = id_q.first;
		double q = id_q.second;
		double p = pValues.at(key);
		double U = Ustats.at(key);
		if ( q < HIGH_FDR_THRESHOLD && p < HIGH_SIGNIF_THRESHOLD ) {
			Association mva(t, p, q, U);
			auto it_ok = mhcVlAssociations.emplace(key, std::list<Association>());
			it_ok.first->second.push_back(mva);
			new_associations.emplace(key, mva); // add to return value
		} // if significant
	}
	return new_associations;
}


std::pair<Virus, bool> World::getConsensus() const {
	std::list<Virus> viruses;
	for ( Host* host : livingHosts ) {
		if ( host == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		viruses.push_back(host->getVirus()); // current dominant strain from the host
	}
	if ( viruses.empty() ) {
		// try to make a consensus from the known wildtypes
		for ( auto & id_vir : wildtypes ) {
			viruses.push_back(id_vir.second);
		}
		if ( viruses.empty() ) { // if viruses is still empty, throw an exception
			throw std::logic_error("wildtypes is empty" + RIGHT_HERE);
		} // else
		return std::make_pair(Virus(&fitnessFunction, viruses), false); // false indicates failure
	} // else, true indicates success
	return std::make_pair(Virus(&fitnessFunction, viruses), true); // consensus constructor
}

std::pair<Virus, bool> World::sampleVirus(Rng & rng) const { // rng masks World::rng
	// returns a virus sampled from the population (R indiv is used for weighing)
	double Rtot = 0.0;
	for ( auto host : livingHosts ) {
		Rtot += host->getIndividualReproductionNumber();
	}
	if ( Rtot <= 0.0 ) {
		if ( wildtypes.empty() ) {
			throw std::logic_error("wildtypes is empty" + RIGHT_HERE);
		} // else, sample from wildtypes.
		auto wtit = wildtypes.begin();
		std::advance(wtit, rng.Integer(wildtypes.size()));
		return std::make_pair(wtit->second, false);
	} // else, sample a host
	double Rpart = rng.Uniform(0.0, Rtot);
	const Host* sampled_host = nullptr;
	for ( auto host : livingHosts ) {
		if ( host == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		Rpart -= host->getIndividualReproductionNumber();
		sampled_host = host; // make sure that there is always a non-NULL host sampled
		if ( Rpart <= 0.0 ) {
			break;
		}
	}
	Virus sampled_virus = sampled_host->getVirusSpvl(rng);
	return std::make_pair(sampled_virus, true);
}

std::pair<double, bool> World::getMeanDistanceFromWt() const {
	WARN_DEPRECATED_FUN
	int H = 0;
	int n = livingHosts.size();
	if ( n == 0 ) { // first get the trivial case out of the way
		return std::make_pair(0.0, false);
	} // else: n > 0 (we can divide by it)
	for ( Host* host : livingHosts ) {
		if ( host == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		// find the right wildtype
		auto wtit = wildtypes.find(host->getVirus().getLineageId());
		if ( wtit == wildtypes.end() ) {
			throw std::logic_error("lineage ID not a key of wildtypes" + RIGHT_HERE);
		} // else...
		H += hammingDistance(wtit->second, host->getVirus());
	}
	double meanHD = double(H) / (n * GENOME_SIZE);
	return std::make_pair(meanHD, true); // ok
}

/* for each virus locus, get statistics about the fitness cost of
 * mutating the locus.
 */
std::string World::xmlNodeFitnessCosts() const {
	std::stringstream ss;
	ss << "<fitness_costs t='" << t << "' >" << std::endl;
	std::string name; // name for the statistic can be reused
	for ( int i = 0; i < GENOME_SIZE; ++i ) {
		name = "basic_mean_fitness_cost_" + std::to_string(i);
		ss << xmlNodeHostArgTrait(&Host::getRelBasicMalthusian, i, name) << std::endl;
		name = "early_mean_fitness_cost_" + std::to_string(i);
		ss << xmlNodeHostArgTrait(&Host::getEarlyRelMalthusian, i, name) << std::endl;
		name = "late_mean_fitness_cost_" + std::to_string(i);
		ss << xmlNodeHostArgTrait(&Host::getLateRelMalthusian, i, name) << std::endl;
		name = "early_dom_fitness_cost_" + std::to_string(i);
		ss << xmlNodeHostArgTrait(&Host::getEarlyDomRelMalthusian, i, name) << std::endl;
		name = "late_dom_fitness_cost_" + std::to_string(i);
		ss << xmlNodeHostArgTrait(&Host::getLateDomRelMalthusian, i, name) << std::endl;
		name = "early_mean_fitness_cost_mimresp_" + std::to_string(i);
		ss << xmlNodeHostArgTrait(&Host::getEarlyRelMalthusianMimicResps, i, name) << std::endl;
		name = "late_mean_fitness_cost_mimresp_" + std::to_string(i);
		ss << xmlNodeHostArgTrait(&Host::getLateRelMalthusianMimicResps, i, name) << std::endl;
	}
	ss << "</fitness_costs>";
	return ss.str();
}


std::string World::xmlNodeGeneticDiversity() const {
	auto diversities = getGeneticDiversityPerLocus();
	std::stringstream ss;
	ss << "<gen_div >" << std::endl;
	std::string sep = "";
	for ( auto d : diversities ) {
		ss << sep << d;
		sep = " ";
	}
	ss << "\n" << "</gen_div>";
	return ss.str();
}

std::string World::xmlNodeAidsHazard() {
	/** Measure the impact of pre-adaptation in terms of loss of AIDS-free
	 * years. Make copies of Hosts, and infect them with an ancestoral
	 * virus, then compare the durations of the infection.
	 * Take level of pre-adaptation and MHC type into account. (TODO)
	 */
	std::cout << "\t[infecting cloned hosts...]" << std::flush;
	std::map<Host*, Host*> twins;
	int twin_counter = 0;
	for ( Host* host : livingHosts ) {
		// dont test ALL hosts...
		if ( cohort_size <= twin_counter++ ) {
			break;
		}
		// clone the Host and run a simulation for the wildtype
		if ( host == nullptr ) {
			throw std::logic_error("Host in livingHosts is NULL" + RIGHT_HERE);
		}
		// there may be multiple wildtype viruses, choose the right lineage
		auto wtit = wildtypes.find(host->getVirusTF().getLineageId());
		if ( wtit == wildtypes.end() ) {
			throw std::logic_error("can't find key in wildtypes" + RIGHT_HERE);
		} // else...
		Host* clone = new Host(*host, wtit->second, rng);
		twins.emplace(host, clone);
		workerpool.addNewJob(clone);
	}
	// sync the workerpool
	workerpool.syncWorkerThreads(false); // false means: do not show a timer...
	// fetch results
	std::stringstream ss;
	ss << "<twins t='" << t << "' >" << std::endl;
	for ( auto twin : twins ) {
		Host* original = twin.first;
		Host* clone = twin.second;
		double D_original = original->getLengthOfInfection();
		double D_clone = clone->getLengthOfInfection();
		auto wpa_original_ok = original->getRelWeightedPreAdaptationTF(clone->getVirusTF());
		ss << "<twin>" << std::endl;
		ss << "<original "
		   << "infection_length='" << D_original << "' "
			 << "rel_wpa='" << wpa_original_ok.first << "' "
			 << "rel_wpa_ok='" << std::boolalpha << wpa_original_ok.second << std::noboolalpha << "' "
			 << "/>" << std::endl;
		ss << "<clone "
		   << "infection_length='" << D_clone << "' "
			 << "/>" << std::endl;
		ss << "</twin>" << std::endl;
		// cleanup
		delete clone; // NB: DO NOT DELETE THE ORIGINAL
	}
	ss << "</twins>";
	return ss.str();
}


void World::printXmlHeritDataNode(std::ostream & os) const {
	os << "<herit_data t='" << t << "' >" << std::endl;
	for ( Host* rec : livingHosts ) {
		if ( rec == nullptr ) {
			throw std::logic_error("ERROR: Host in livingHosts is NULL" + RIGHT_HERE);
		}
		Host* tra = rec->getParent(); // in most cases a valid pointer (unless Patient Zero)
		if ( tra != nullptr ) {
			os << "<herit_data_point " // opening angle
			   << "spvl_rec='" << rec->getLogSpvl() << "' "
			   << "wTF_rec='" << rec->getWTF() << "' "
				 << "responses_TF_rec='" << rec->getNumberOfResponsesTF() << "' ";
			os << "spvl_tra='" << tra->getLogSpvl() << "' "
			   << "wTF_tra='" << tra->getWTF() << "' "
				 << "responses_TF_tra='" << tra->getNumberOfResponsesTF() << "' ";
			MhcHaplotype hapl_tra = tra->getMhcHaplo();
			MhcHaplotype hapl_rec = rec->getMhcHaplo();
			os << "jaccard_mhc='" << jaccardMhcHaplo(hapl_rec, hapl_tra) << "' ";
			os << "/>" << std::endl; // closing angle
              // TODO: add jaccard index based on RESPONSES against the viruses (?)
		} // if tra != NULL (valid transmission couple)
	}
	os << "</herit_data>";
}

void World::printXmlGeneralStats(std::ostream & os) const {
	// simple statistics
	int dI = recentInfectionTimes.size();
	double frFailedInfections = 0.0;
	if ( cumulativeIncidence > 0 ) {
		frFailedInfections = double(failedInfections)/(failedInfections + cumulativeIncidence);
	}
	int I = livingHosts.size();
	int N = getPopulationSize();
	int S = N - I;
	// fraction of susceptibles
	double s = 0.0;
	if ( N > 0 ) s = double(S)/N; // probably a redundant test

	// print some basic stats as attributes
	os << "<state " // opening tag
		 << "t='" << t << "' "
		 << "dI='" << dI << "' "
		 << "I='" << I << "' "
		 << "S='" << S << "' "
		 << "N='" << N << "' "
		 << "fr_fail='" << frFailedInfections << "' "
		 << "cluster_coeff='" << getGlobalClusterCoeff() << "' ";
	os << ">" << std::endl; // close the <state >

	// xml nodes for AgentTraits
	os << xmlNodeAgentTrait(&Agent::getAge, "age") << std::endl;
	os << xmlNodeAgentArgTrait(&Agent::degree, false, "degree") << std::endl; // all contacts
	os << xmlNodeAgentArgTrait(&Agent::degree, true, "degree_tmr") << std::endl; // only transmission routes

	// xml nodes for heritabilities
	os << xmlNodeHeritabilityHostTrait(&Host::getLogSpvl, "h2_spvl") << std::endl;
	os << xmlNodeHeritabilityHostTrait(&Host::getLogGeomSpvl, "h2_geom_spvl") << std::endl;
	os << xmlNodeHeritabilityHostTrait(&Host::getLogGsvl, "h2_gsvl") << std::endl;
	os << xmlNodeHeritabilityHostTrait(&Host::getWTF, "h2_wTF") << std::endl;
	os << xmlNodeHeritabilityHostTrait(&Host::getMalthusianTF, "h2_r0_TF") << std::endl;
	os << xmlNodeHeritabilityHostTrait(&Host::getGenericMalthusianTF, "h2_r0_generic_TF") << std::endl;
	os << xmlNodeHeritabilityHostTrait(&Host::getSlopeCD4, "h2_slope_CD4") << std::endl;
	// xml nodes for correlations
	os << xmlNodeCorrelationHostTrait(&Host::getWTF, &Host::getLogSpvl, "corr_wTF_spvl") << std::endl;
	os << xmlNodeCorrelationHostTrait(&Host::getNumberOfResponsesTF, &Host::getLogSpvl, "corr_responses_TF_spvl") << std::endl;
	os << xmlNodeCorrelationHostTrait(&Host::getMalthusianTF, &Host::getLogSpvl, "corr_r0_TF_spvl") << std::endl;
	os << xmlNodeCorrelationHostTrait(&Host::getGenericMalthusianTF, &Host::getLogSpvl, "corr_r0_generic_TF_spvl") << std::endl;
	os << xmlNodeCorrelationHostTrait(&Host::getLogSpvl, &Host::getSlopeCD4, "corr_spvl_slope_CD4") << std::endl;
	os << xmlNodeCorrelationHostRArgTrait(&Host::getLogVl, &Host::getPLevMalthusian, s, "corr_log_vl_plev_malthusian") << std::endl;
	// xml nodes of HostTraits
	os << xmlNodeHostTrait(&Host::getLogSpvl, "spvl") << std::endl;
	os << xmlNodeHostTrait(&Host::getLogGeomSpvl, "geom_spvl") << std::endl;
	os << xmlNodeHostTrait(&Host::getLogGsvl, "gsvl") << std::endl;
	os << xmlNodeHostTrait(&Host::getEarlyLogSpvl, "early_spvl")  << std::endl;
	os << xmlNodeHostTrait(&Host::getLogVlMultiplierTF, "log_vl_multiplier_TF") << std::endl;
	os << xmlNodeHostTrait(&Host::getIndividualReproductionNumber, "R_indiv")  << std::endl;
	os << xmlNodeHostTrait(&Host::getVl, "vl")  << std::endl;
	os << xmlNodeHostArgTrait(&Host::getPLevMalthusian, s, "plev_malthusian") << std::endl; // s = S/N
	os << xmlNodeHostTrait(&Host::getWTF, "wTF")  << std::endl;
	os << xmlNodeHostTrait(&Host::getW, "w")  << std::endl;
	os << xmlNodeHostTrait(&Host::getR0TF, "R0_TF") << std::endl;
	os << xmlNodeHostTrait(&Host::getMalthusianTF, "malthusian_TF") << std::endl;
	os << xmlNodeHostTrait(&Host::getGenericMalthusianTF, "generic_malthusian_TF") << std::endl;
	os << xmlNodeHostTrait(&Host::getNumberOfResponsesTF, "responses_TF") << std::endl;
	os << xmlNodeHostTrait(&Host::getNumberOfDomResponsesTF, "dom_responses_TF") << std::endl;
	os << xmlNodeHostTrait(&Host::getFitnessCostOfEscape, "cost_esc") << std::endl;
	os << xmlNodeHostMaybeTrait(&Host::getEscapeRate, "escape_rate") << std::endl;
	os << xmlNodeHostMaybeTrait(&Host::getEscapeCost, "escape_cost") << std::endl;
	os << xmlNodeHostTrait(&Host::getEscapeDensity, "escape_density") << std::endl;
	os << xmlNodeHostMaybeTrait(&Host::getCompensationRate, "compensation_rate") << std::endl;
	os << xmlNodeHostTrait(&Host::getCompensationDensity, "compensation_density") << std::endl;
	os << xmlNodeHostMaybeTrait(&Host::getCompensationRateTrunk, "compensation_rate_trunk") << std::endl;
	os << xmlNodeHostTrait(&Host::getCompensationDensityTrunk, "compensation_density_trunk") << std::endl;
	os << xmlNodeHostTrait(&Host::getSlopeCD4, "slope_CD4") << std::endl;
	// make xml nodes for the w-parameters of the genes
	for ( int i = 0; i < NUMBER_OF_GENES; ++i ) {
		os << xmlNodeHostArgTrait(&Host::getWTF, i, geneNameMap.at(i) + "_TF") << std::endl;
	}
	os << xmlNodeGeneticDiversity() << std::endl; // diversity per locus
	// pre-adaptation
	os << xmlNodeHostArgTrait(&Host::getPreAdaptationTF_vm, wildtypes, "wildtype_pre_adaptation_TF") << std::endl;
	os << xmlNodeHostArgTrait(&Host::getWeightedPreAdaptationTF_vm, wildtypes, "wildtype_weighted_pre_adaptation_TF") << std::endl;
	os << xmlNodeHostMaybeArgTrait(&Host::getRelPreAdaptationTF_vm, wildtypes, "wildtype_rel_pre_adaptation_TF") << std::endl;
	os << xmlNodeHostMaybeArgTrait(&Host::getRelWeightedPreAdaptationTF_vm, wildtypes, "wildtype_rel_weighted_pre_adaptation_TF") << std::endl;
	// rate of viral evolution
	os << xmlNodeHostArgTrait(&Host::getScaledHammingDistance_vm, wildtypes, "dist_wildtype") << std::endl;
	// pre-adaptation w.r.t current consensus
	auto cons_ok = getConsensus();
	if ( !cons_ok.second ) {
		std::cerr << "WARNING: unable to compute consensus virus" << RIGHT_HERE << std::endl;
	}
	os << xmlNodeHostArgTrait(&Host::getPreAdaptationTF, cons_ok.first, "consensus_pre_adaptation_TF") << std::endl;
	os << xmlNodeHostArgTrait(&Host::getWeightedPreAdaptationTF, cons_ok.first, "consensus_weighted_pre_adaptation_TF") << std::endl;
	// distance from the consensus
	os << xmlNodeHostArgTrait(&Host::getScaledHammingDistance, cons_ok.first, "dist_consensus") << std::endl; // serves as measure of genetic diversity


	// write MHC data to the stats file
	for ( auto lit = allMhcLoci.begin(); lit != allMhcLoci.end(); ++lit ) { // loop over loci
		for ( auto ait = lit->begin(); ait != lit->end(); ++ait ) { // loop over alleles
			BasicStats nStats("n");
			BasicStats dwnStats("dwn");

			getPlevEpitopesForAllele(*ait, nStats, dwnStats);
			// write to file
			os << "<mhc_allele_stats >" << std::endl
				 << ait->xmlStringMhcId() << std::endl
				 << nStats << std::endl
				 << dwnStats << std::endl
				 << "</mhc_allele_stats>" << std::endl;
		}
	}
	os << "</state>";
}


void World::updateAncestors() {
	auto it = ancestors.begin();
	while ( it != ancestors.end() ) {
		Host* ancestor = (*it);
		if ( !ancestor->coalesces() && ancestor->getRemoveIfExtLineageFlag() ) {
			// update the MRCA of all living hosts
      auto pr = ancestor->getChild();
      if ( pr.second ) { // there was a child
        Host* new_ancestor = pr.first;
      	if ( new_ancestor != nullptr ) {
          new_ancestor->setMrcaFlag(true); // don't recursively remove the new ancestor
        	ancestors.push_back(new_ancestor);
        } else {
          throw std::logic_error("child is NULL" + RIGHT_HERE);
        }
      } else { // if there is no child, we have an extinct lineage
        std::cout << "\t[extinct lineage]" << std::endl;
      }
			// write the ancestor to a file
      if ( ofancs ) {
        ofancs << (*ancestor) << std::endl;
      }	else {
        throw std::runtime_error("cannot write to ancestor-file" + RIGHT_HERE);
      }
			// cleanup and get the next iterator
			delete ancestor;
			it = ancestors.erase(it);
		}	else {
			++it;
		}
	}
}

void World::updateRecentInfectionTimes(int n) {
	// add infection times to list
	std::list<double> ts(n, t);
  recentInfectionTimes.splice(recentInfectionTimes.end(), ts);
  // remove any times more than a year ago.
  auto pred = [&](double s){ return s < t - INCIDENCE_TIME_UNIT; };
  recentInfectionTimes.remove_if(pred);
}
