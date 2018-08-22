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

#include "network.hpp"

/* methods for class Contact */

Contact::Contact(Agent* left, Agent* right, double t, Rng & rng) : left(left), right(right) {
	start = t;
	stop = t + rng.LogNormalMnSd(MEAN_CONTACT_DURATION, SD_CONTACT_FORMATION);
	transmission_route = false;
	// add the new contact to both Agents
	if ( left != nullptr && right != nullptr ) {
		left->addContact(this);
		right->addContact(this);
	}	else {
		throw std::invalid_argument("initiating contact with NULL pointer(s)" + RIGHT_HERE);
	}
	setRandomIndex(rng);
}

Contact::~Contact() {
	if ( left != nullptr ) {
		left->removeContact(this);
	}
	if ( right != nullptr ) {
		right->removeContact(this);
	}
}

bool Contact::isMember(const Agent* a) const {
	return ( left == a || right == a );
}
Agent* Contact::getLeft() const {
	return left;
}
Agent* Contact::getRight() const {
	return right;
}
Agent* Contact::getOther(const Agent* a) const {
	if ( left == a ) return right;
	else if ( right == a ) return left;
	else return nullptr; // warning...
}
bool Contact::isExpired(double t) const {
  return stop <= t;
}
bool Contact::isTransmissionRoute() const {
	return transmission_route;
}
void Contact::setTransmissionRoute(bool x) {
	transmission_route = x;
}

/* methods of class Agent */

Agent::Agent(double t, Rng & rng, Id id, InitType itp) {
  switch ( itp ) {
    case NEWBORN: {
      age = 0.0;
      timeOfBirth = t-age;
      double lifespan = rng.Carnes();
      timeOfNaturalDeath = timeOfBirth + lifespan;
      break;
    }
    case RANDAGE: {
      age = rng.CarnesAge();
      timeOfBirth = t-age; // can be negative, e.g. if t=0...
      double lifespan = rng.Carnes(age);
      timeOfNaturalDeath = timeOfBirth + lifespan;
      break;
    }
    default: {
      throw std::logic_error("initiating Agent with invalid InitType" + RIGHT_HERE);
      break; // redundant
    }
  } // switch on "init type"
  // choose inividual-level parameters
  ipars = IndPars(rng);
	// choose a random Gender (50-50).
	gender = (rng.Bernoulli(0.5) ? MALE : FEMALE);
  // heterogeneity in contact formation rates
  cfRate = rng.LogNormalMnSd(MEAN_CONTACT_FORMATION_RATE, SD_CONTACT_FORMATION_RATE);
	// set the maximum number of contacts (TODO: change name: transmission routes)
	double promiscuity = cfRate / MEAN_CONTACT_FORMATION_RATE;
	maxNumContacts = rng.PosPoisson(MEAN_NUMBER_TRANSMISSION_STUBS * promiscuity);
  // cfThreshold will be re-sampled after the first contact is made
  assignCfThreshold(rng); // sets cfLoad and cfThreshold
  // determine a infection threshold and set the load to 0
  assignInfectionThreshold(rng);
  // the Agent is not infected yet, hence host is nullptr
  host = nullptr;
  // set the member id for identification during analysis
  this->id = id;
  setRandomIndex(rng);
}

Agent::~Agent() {
	if ( host != nullptr ) {
		/* whenever the host Agent is deleted (because of death or the
		 * end of the simulation), the remove_if_ext_lineage flag of the
		 * Host is set to true, so that is is allowed to be recursively
		 * deleted.
		 */
		host->setRemoveIfExtLineageFlag(true);
		if ( host->deleteMe() ) {
			delete host;
		}
	}
	// remove any contacts with this Agent
	removeAllContacts();
}

void Agent::updateState(double dt, Rng & rng) {
  // update age
  age += dt;
  // remove any expired contacts
  removeExpiredContacts(age+timeOfBirth);
  // update contact formation load
	cfLoad += cfRate * dt;
  // update infection load and threshold.
	infectionLoad += dt * getForceOfInfection();
	if ( isInfected() ) {
		host->nextState();
	}
	// sample a new random_index for shuffling the population list
	setRandomIndex(rng);
}

Id Agent::getId() const { return id; }
Host* Agent::getHostPtr() const { return host; }
double Agent::getAge() const { return age; }
Gender Agent::getGender() const { return gender; }

double Agent::getTimeOfDeath() const {
	if ( !isInfected() ) {
		return getTimeOfNaturalDeath();
	}	else {
		return std::min(host->getTimeOfDeath(), getTimeOfNaturalDeath());
	}
}
double Agent::getTimeOfNaturalDeath() const { return timeOfNaturalDeath; }
double Agent::getTimeOfBirth() const { return timeOfBirth; }

bool Agent::isAlive() const {
	if ( age + timeOfBirth < getTimeOfNaturalDeath() ) {
		if ( isInfected() ) {
			return host->isAlive();
		}	else {
			return true;
		}
	}	else {
		return false;
	}
}

bool Agent::isInfected() const { return host != nullptr; }

double Agent::getTransmissionRate() const {
	return ( isInfected() ? host->getInfectionRate() : 0.0 );
}

double Agent::getForceOfInfection() const {
  double foi = 0.0; // as perceived by this agent
  for ( Contact* c : contacts ) {
		if ( c->isTransmissionRoute() ) {
			Agent* partner = c->getOther(this);
			foi += partner->getTransmissionRate(); // 0.0 if the partner is not infected
		}
	}
	return foi;
}

std::unordered_map<Agent*, bool> Agent::getPartners() const {
	std::unordered_map<Agent*, bool> neighbors;
	for ( auto it = contacts.begin(); it != contacts.end(); ++it ) {
		Contact* contact = *it;
		Agent* nb = contact->getOther(this);
		if ( nb == nullptr ) {
			throw std::logic_error("Contact::getOther returned NULL" + RIGHT_HERE);
		}
		neighbors.emplace(nb, contact->isTransmissionRoute());
	}
	return neighbors;
}

bool Agent::hasContactWith(Agent* other) const {
	for ( auto it = contacts.begin(); it != contacts.end(); ++it ) {
		Contact* c = *it;
		if ( c->isMember(other) ) return true;
	}
	return false;
}

bool Agent::hasContactInCommon(Agent* other) const {
	for ( auto it = contacts.begin(); it != contacts.end(); ++it ) {
		Contact* contact = *it;
		Agent* neighbor = contact->getOther(this);
		if ( neighbor == nullptr ) {
			throw std::logic_error("Contact::getOther returned NULL" + RIGHT_HERE);
		}
    if ( other->hasContactWith(neighbor) ) {
      return true;
    }
	}
  return false;
}

int Agent::getNumSharedPartners(Agent* other) const { // get the exact number of common friends
  int result = 0;
	for ( auto it = contacts.begin(); it != contacts.end(); ++it ) {
		Contact* contact = *it;
		Agent* neighbor = contact->getOther(this);
		if ( neighbor == nullptr ) {
			throw std::logic_error("Contact::getOther returned NULL" + RIGHT_HERE);
		}
  	if ( other->hasContactWith(neighbor) ) {
      result++;
    }
	}
  return result;
}

int Agent::degree(bool transmission_routes) const {
	if ( transmission_routes ) {
		auto pred = [](Contact* c){ return c->isTransmissionRoute(); };
		return std::count_if(contacts.begin(), contacts.end(), pred);
	} else {
		return contacts.size();
	}
}

double Agent::clustercoeff() const {
	int k = degree();
	if ( k <= 2 ) {
		return 1.0;
	} // else, k >= 3
	int t = 0; // count the number of "triangles"
	for ( auto it1 = contacts.begin(); it1 != contacts.end(); ++it1 ) {
		Contact* c1 = *it1;
		Agent* nb1 = c1->getOther(this);
		if ( nb1 == nullptr ) {
			throw std::logic_error("Contact::getOther returned NULL" + RIGHT_HERE);
		}
		for ( auto it2 =contacts.begin(); it2 != it1; ++it2 ) {
			Contact* c2 = *it2;
			Agent* nb2 = c2->getOther(this);
			if ( nb2 == nullptr ) {
				throw std::logic_error("Contact::getOther returned NULL" + RIGHT_HERE);
			}
			// two neighbors are friends: increment the number of "triangles"
			t += ( nb1->hasContactWith(nb2) ? 1 : 0 );
		}
	}
	return double(t) / ((k*(k-1))/2);
}

int Agent::getNumFreeStubs() const {
	return maxNumContacts - degree(true); // passing optional true gives only transmission routes
}

int Agent::getMaxNumContacts() const {
    return maxNumContacts;
}


void Agent::addContact(Contact* c) {
	contacts.push_back(c);
}
void Agent::removeContact(Contact* c) {
	contacts.remove(c);
}

void Agent::removeAllContacts() {
	while ( !contacts.empty() ) {
		Contact* c = contacts.front();
		/* the destructor of Contact calls Agent::removeContact for this
		 * and the other Agent in the contact.
		 */
		delete c;
	}
}

void Agent::removeExpiredContacts(double t) {
  // find expired contacts
  std::list<Contact*> expired;
  for ( auto it = contacts.begin(); it != contacts.end(); ++it ) {
    Contact* c = *it;
    if ( c->isExpired(t) ) expired.push_back(c);
  }
  // delete them
  for ( auto it = expired.begin(); it != expired.end(); ++it ) {
	Contact* c = *it;
	/* the destructor of Contact calls Agent::removeContact for this
	 * and the other Agent in the contact.
	 */
	delete c;
  }
}

int Agent::enableTransmissionRoutes(Rng & rng) {
	/* make a transmission route from a contact,
	 * if the contact satisfies certain conditions...
	 */
	if ( age < MIN_CONTACT_AGE ) {
		return 0;
	} // else, try to fill a free stub
	int enabled_tmrs = 0; // return value
	int free_stubs = getNumFreeStubs();
	shuffleContacts(rng); // walk through the contacts in a random order
	for ( Contact* contact : contacts ) {
		// test if there are any free stubs left, stop otherwise
		if ( enabled_tmrs >= free_stubs ) {
			break;
		}
		Agent* nb = contact->getOther(this);
		if ( nb == nullptr ) {
			throw std::logic_error("Contact::getOther returned NULL" + RIGHT_HERE);
		}
		// test if contact is allowed
		if ( !contact->isTransmissionRoute()
				 && nb->getGender() != gender
				 && nb->getAge() >= MIN_CONTACT_AGE
				 && nb->getNumFreeStubs() > 0 ) {
			contact->setTransmissionRoute(true);
			// and increase the counter
			enabled_tmrs++;
		}
	}
	// END TESTING
	return enabled_tmrs;
}

bool Agent::requiresNewContact() const {
  return cfLoad > cfThreshold;
}

double Agent::getCfProbability(double dt) const {
	return 1.0 - exp(-dt*cfRate);
}

void Agent::assignCfThreshold(Rng & rng) {
  cfLoad = 0.0;
  cfThreshold = rng.Exponential(1.0);
}

bool Agent::requiresInfection() const {
	return infectionLoad > infectionThreshold;
}

void Agent::assignInfectionThreshold(Rng & rng) {
	infectionLoad = 0.0;
	infectionThreshold = rng.Exponential(1.0);
}

bool Agent::infectByRandomPartner(Rng & rng, const std::vector<MhcLocus> & allMhcLoci) {
	bool ok = false;
	double foi = getForceOfInfection();
	// sample a number between 0 and the force of infection to determine the infecting host
	double x = rng.Uniform(0.0, foi);
	for ( auto it = contacts.begin(); it != contacts.end(); ++it ) {
		Contact* contact = *it;
		Agent* partner = contact->getOther(this);
		x -= partner->getTransmissionRate();
		if ( x < 0 ) {
			Host* partners_host_obj = partner->getHostPtr();
			if ( partners_host_obj == nullptr ) {
				throw std::logic_error("expected non-NULL host pointer" + RIGHT_HERE);
			}
			// sample a virus from partner, taking relative clone size into account
			Virus virus = partners_host_obj->getVirus(rng);
			ok = infect(rng, allMhcLoci, virus);
			if ( ok ) {
				host->setAddressParent(partners_host_obj);
				partners_host_obj->addChild(host);
			}
			break;
		}
	}
	return ok;
}

bool Agent::infect(Rng & rng, const std::vector<MhcLocus> & allMhcLoci, const Virus & virus) {
	if ( !isInfected() ) {
		// determine if the virus is capable of infecting this Agent
		double R0 = iLevBasicReproductionNumber(virus, ipars);
		bool infection_successful = false;
		if ( INFECTION_ILEV_R0_DEPENDENT ) {
			/* use the R0-dependent major outbreak probability to determine if
			 * the Host can be infected
			 */
			WARN_UNTESTED_FUN
			double pr = major_outbreak_prob(R0);
			infection_successful = rng.Bernoulli(pr);
		} else if ( R0 > 1.0 ) {
			// use a strict threshold R0 <= 1 or R0 > 1
			infection_successful = true;
		}
		// now create a Host object if infection is possible.
		if ( infection_successful ) {
			host = new Host(rng, allMhcLoci, virus, ipars, timeOfBirth + age, id);
			return true;
		}	else {
			return false;
		}
	}	else {
		return false;
	}
}

void Agent::shuffleContacts(Rng & rng) {
	std::for_each(contacts.begin(), contacts.end(), [&](Contact* c){c->setRandomIndex(rng);});
	contacts.sort(compareByRandomIndex);
}

bool compareAgentByAge(Agent* left_agent, Agent* right_agent) {
	if ( left_agent != nullptr && right_agent != nullptr ) {
		return left_agent->getAge() < right_agent->getAge();
	}	else {
		if ( right_agent == nullptr ) return true;
		else return false;
	}
}

bool compareAgentByDegree(Agent* left_agent, Agent* right_agent, bool transmission_routes) {
	if ( left_agent != nullptr && right_agent != nullptr ) {
		return left_agent->degree(transmission_routes) < right_agent->degree(transmission_routes);
	}	else {
		if ( right_agent == nullptr ) return true;
		else return false;
	}
}

bool compareAgentByDegree(Agent* left_agent, Agent* right_agent) {
	return compareAgentByDegree(left_agent, right_agent, false);
}

bool compareAgentByTmrDegree(Agent* left_agent, Agent* right_agent) {
	return compareAgentByDegree(left_agent, right_agent, true);
}


/** methods of the Network class **/

int Network::getPopulationSize() const { return population.size(); }

void Network::makeContacts(double t, double dt, Rng & rng) {
	// pop_vec is a population vector for random access
	std::vector<Agent*> pop_vec;
	pop_vec.reserve(population.size()); // reserve space
	pop_vec.assign(population.begin(), population.end()); // fill vector
	// walk through the population
  for ( Agent* agent_left : population ) {
		if ( !agent_left->requiresNewContact() ) {
			continue; // skip agents that don't require a new contact.
		} // else...
		// first calculate the probability of forming any contact
		double pr_no_contact = 1.0;
		// and also get thus sum of the contact probabilities for later use
		double sum_pr = 0.0;
  	std::unordered_map<Agent*, double> potential_contacts;
		// an auxiliary closure to update probabilities
		auto update_probs = [&](Agent* agent_left, Agent* agent_right) {
			if ( agent_left != agent_right && !agent_left->hasContactWith(agent_right) ) {
				double pr_age = age_assort_prob(agent_left->getAge(), agent_right->getAge());
				double pr_con = concurrency_assort_prob(agent_left->getMaxNumContacts(),
						agent_right->getMaxNumContacts());
				double pr = pr_age * pr_con;
				auto ok = potential_contacts.emplace(agent_right, pr);
				// emplace retuns true if the agent_right was inserted
				if ( ok.second ) {
					pr_no_contact *= (1.0 - pr);
					sum_pr += pr;
				} // else, the Agent was already in the map
			} // not yet a contact
		}; // auxiliary closure
		// first, look for contacts at distance 2: i.e. friends of friends
		if ( rng.Bernoulli(CONTACT_CLUSTERING_BIAS) ) {
			for ( auto pair1 : agent_left->getPartners() ) {
				Agent* partner = pair1.first;
				if ( partner == nullptr ) {
					throw std::logic_error("Agent has NULL partner" + RIGHT_HERE);
				}
				for ( auto pair2 : partner->getPartners() ) {
					Agent* agent_right = pair2.first;
					if ( agent_right == nullptr ) {
						throw std::logic_error("Agent has NULL partner" + RIGHT_HERE);
					}
					update_probs(agent_left, agent_right);
				} // loop over partners of partners
			} // loop over partners
		} // with probability CONTACT_CLUSTERING_BIAS, first look in social nbh
		/* now choose a contact, based on the age difference, condictional on pr_any_contact
		 * suppose that all ages (concurrencies, etc.) have equal contact formation probability p,
		 * then 1 - (1 - p)^n is about n*p for small enough p... TODO
		 */
		if ( rng.Bernoulli(1.0 - pr_no_contact) ) {
			// sample an Agent from potential_contacts
			double u = rng.Uniform(0.0, sum_pr);
			for ( auto pair : potential_contacts ) {
				Agent* agent_right = pair.first;
				if ( agent_right == nullptr ) {
					throw std::logic_error("Agent in potential_contacts is NULL" + RIGHT_HERE);
				}
				double pr = pair.second;
				u -= pr;
				if ( u < 0.0 ) {
					new Contact(agent_left, agent_right, t, rng); // anonymous pointer :-)
					break;
				} // if u drops below 0.0
			} // for loop over remaining contacts
		} else { // if: should there be a contact in local nbh? else: try in entire population
			/* this for loop should be a loop over the population,
			 * but we have no good way of shuffling it.
			 * Therefore use a vector copy that has fast random access (operator[])
			 */
			int max_tries = population.size(); // TODO What is a good cut-off?
			for ( int i = 0; i < max_tries; ++i ) {
				int idx = rng.Integer(pop_vec.size());
				Agent* agent_right = pop_vec[idx];
				if ( agent_right == nullptr ) {
					throw std::logic_error("Agent in copy of population is NULL" + RIGHT_HERE);
				}
				if ( agent_left != agent_right && !agent_left->hasContactWith(agent_right) ) {
					double pr_age = age_assort_prob(agent_left->getAge(), agent_right->getAge());
					double pr_con = concurrency_assort_prob(agent_left->getMaxNumContacts(),
							agent_right->getMaxNumContacts());
					double pr = pr_age * pr_con;
					if ( rng.Bernoulli(pr) ) {
						new Contact(agent_left, agent_right, t, rng); // anonymous pointer :-)
						break; // break the for loop
					} // if true, make contact
				} // if contact allowed (no self of double links)
			} // for loop until max iter reached or contact was formed
		}
		// reset contact formation parameters (also if no contact was made)
		agent_left->assignCfThreshold(rng);
	} // for loop over available agents
}

int Network::makeTransmissionRoutes(Rng & rng) {
	/** walk through the population in a random order,
	 * then ask each Agent to maximize its transmission routes.
	 * function returns the total number of new transmission routes.
	 */
	shufflePopulation(rng);
	int enabled_tmrs = 0;
	for ( Agent* agent : population ) {
		if ( agent == nullptr ) {
			throw std::logic_error("Agent in population is NULL" + RIGHT_HERE);
		}
		enabled_tmrs += agent->enableTransmissionRoutes(rng);
	}
	return enabled_tmrs;
}

void Network::clearPopulation() {
	for ( Agent* agent : population ) {
		delete agent;
	}
	population.clear();
}

void Network::shufflePopulation(Rng & rng) {
	// assign random indices to the agents
	std::for_each(population.begin(), population.end(), [&](Agent* a){a->setRandomIndex(rng);});
	// sort on random index
	population.sort(compareByRandomIndex);
}

void Network::printGraph(std::ostream & os) const {
	/** Draw a graph using the dot language (GraphViz)
	 * Color is used to differentiate between age classes,
	 * the size of the nodes indicates the promiscuity,
	 * and the shape indicates the gender of the individual
	 * light edges indicate auxiliary contacts and heavy edges
	 * indicate potential transmission routes. Infected host nodes
	 * have a red, thick outline.
	 */
	std::map<Gender, std::string> genderShapeMap = {
		std::make_pair(MALE, "circle"),
		std::make_pair(FEMALE, "diamond")
	};
	double node_size = 0.2;
	std::list<Id> already_drawn;
	// to prevent double edges, keep track of already drawn nodes
 	os << "graph {\n"; // start graph
	for ( auto it = population.begin(); it != population.end(); ++it ) {
		Agent* agent = *it;
		if ( agent == nullptr ) {
			throw std::logic_error("Agent in population is NULL" + RIGHT_HERE);
		}
		if ( PRINT_ISOLATED_VERTICES || agent->degree() > 0 ) {
			double age = agent->getAge();
			std::string color = ageToColor(age);
			std::string shape = genderShapeMap.at(agent->getGender());
			os << "\t" << agent->getId() << " [ "; // opening bracket
			// node properties
			os << "style=\"filled\" "
			   << "shape=\"" << shape << "\" "
				 << "width=\"" << node_size * sqrt(agent->getMaxNumContacts()) << "\" "
				 << "regular=\"true\" "
				 << "fixedsize=\"true\" "
			   << "fillcolor=\"" << color << "\" "
			   << "color=\"" << (agent->isInfected() ? "black" : "gray") << "\" "
				 << "penwidth=\"" << (agent->isInfected() ? "3" : "1") << "\" ";
			// closing bracket
			os << "];" << std::endl;
			// edges: make stringstreams for two types of edges
			std::stringstream edge_ss;
			std::stringstream tmr_edge_ss; // transmission routes
			edge_ss << "\t" << agent->getId() << " -- { ";
			tmr_edge_ss << "\t" << agent->getId() << " -- { ";
			std::unordered_map<Agent*, bool> neighbors = agent->getPartners();
			for ( auto pair : neighbors ) {
				Agent* partner = pair.first;
				bool transmission_route = pair.second;
				if ( std::find(already_drawn.begin(), already_drawn.end(), partner->getId()) == already_drawn.end() ) {
					if ( transmission_route ) {
						tmr_edge_ss << partner->getId() << " ";
					} else {
						edge_ss << partner->getId() << " ";
					}
				}
			}
			edge_ss << "} [color=gray];\n"; // light color for background contacts
			tmr_edge_ss << "} [penwidth=3];\n"; // make transmission routes thicker
			// now write both stringstreams to the ostream
			os << edge_ss.str() << tmr_edge_ss.str();
			already_drawn.push_back(agent->getId());
		}
	}
	os << "}"; // close graph
}

std::map<int, int> Network::getDegreeDistribution() const {
	std::map<int, int> ddist;
	for ( auto it = population.begin(); it != population.end(); ++it ) {
		Agent* agent = *it;
		if ( agent == nullptr ) {
			throw std::logic_error("Agent in population is NULL" + RIGHT_HERE);
		}
		int d = agent->degree();
		auto mit = ddist.find(d);
		if ( mit != ddist.end() ) {
			mit->second++;
		}	else {
			ddist[d] = 1;
		}
	}
	return ddist;
}

double Network::getGlobalClusterCoeff() const {
	double triangles = 0.0;
	int triples = 0;
	for ( Agent* agent : population ) {
		if ( agent == nullptr ) {
			throw std::logic_error("Agent in population is NULL" + RIGHT_HERE);
		}
		int k = agent->degree();
		int trpl = (k*(k-1))/2; // integer division, but k*(k-1) is always even
		triangles += trpl * agent->clustercoeff();
		triples += trpl;
	}
	return triangles / triples;
}
