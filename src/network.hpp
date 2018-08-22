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

#ifndef NETWORK_HPP_
#define NETWORK_HPP_

#include <stdexcept>
#include <list>
#include <unordered_map>

#include "macros.hpp"
#include "aux.hpp" // import Id
#include "rng.hpp"
#include "host.hpp"
#include "hashed_definitions.hpp"

/* Some auxiliary functions for network formation.
 */

inline double age_assort_prob(double left_age, double right_age) {
	/* gives the probability that a contact is formed,
	 * given the ages of the two Agents, and the fact that they both
	 * intend to form a contact.
	 * CHAR_AGE_DIFF is the age difference at which this probability
	 * is 1/2.
	 */
	double u = (left_age - right_age) / CHAR_AGE_DIFF;
	return exp(-(u * u) * M_LN2);
}

inline double concurrency_assort_prob(double left_conc, double right_conc) {
  /* cf. age_assort_prob, but then for number of contacts
   */
  double u = (left_conc - right_conc) / CHAR_CONCURRENCY_DIFF;
  return exp(-(u * u) * M_LN2);
}

/* an Agent can have a Gender, determining the possible sexual contacts
 * Possibly, this could also be used to modulate virus load and
 * transmission rate.
 */

enum Gender {
	MALE=0,
	FEMALE
};

/* the classes Agent and Contact are used to handle the contact network.
 * When infected the host pointer is not nullptr.
 */

class Agent; // forward declaration for Contact definition

class Contact : public RandomIndexable {
public:
	Contact(Agent* , Agent* , double , Rng & ); // inform Agents about construction
	~Contact(); // inform Agents about destuction
// aux functions
	bool isMember(const Agent* ) const;
	Agent* getLeft() const;
	Agent* getRight() const;
	Agent* getOther(const Agent* ) const; // returns NULL if argument is not a member!
  bool isExpired(double ) const; // argument: time t
	bool isTransmissionRoute() const;
	void setTransmissionRoute(bool );
protected:
	Agent* left;
	Agent* right;
	double start;
	double stop;
	bool transmission_route;
};



class Agent : public RandomIndexable {
public:
  enum InitType { NEWBORN=0, RANDAGE }; // determines how the lifespan is sampled
	Agent(double , Rng & , Id , InitType=NEWBORN);
	~Agent();
  // update internal state
  void updateState(double , Rng & ); // dt, RNG (for sampling random_index)
  // get stuff...
  Id getId() const;
	Host* getHostPtr() const;
  double getAge() const;
	Gender getGender() const;
  double getTimeOfDeath() const;
	double getTimeOfNaturalDeath() const;
	double getTimeOfBirth() const; // t
  bool isAlive() const; // based on age and getTimeOfDeath()
  bool isInfected() const; // returns true iff host is not NULL
  double getTransmissionRate() const; // is 0.0 if the Agent is not infected
  double getForceOfInfection() const; // sum of infection rates of partners
  std::unordered_map<Agent*, bool> getPartners() const;
	// test for contact
	bool hasContactWith(Agent* ) const;
  bool hasContactInCommon(Agent* ) const; // do they share a friend?
  int getNumSharedPartners(Agent* ) const; // get the exact number of common friends
	int degree(bool transmission_routes=false) const;
	double clustercoeff() const; // how many of my friends are themselves friendly?
	int getNumFreeStubs() const; // maxNumContacts - degree()
  int getMaxNumContacts() const;
	// add and remove contacts
	void addContact(Contact* );
	void removeContact(Contact* );
	void removeAllContacts();
  void removeExpiredContacts(double ); // t
	int enableTransmissionRoutes(Rng & ); // returns the number of enabled transmission routes
  // new contacts
  bool requiresNewContact() const; // true iff cfLoad > cfThreshold
  double getCfProbability(double ) const; // dt
  void assignCfThreshold(Rng & ); // also sets cfLoad to zero
  // infection of the Agent...
  bool requiresInfection() const; // true if infectionLoad > infectionThreshold
  /* assignInfectionhreshold should only be called once
   * (in the future: co/super infection)
   */
  void assignInfectionThreshold(Rng & );
  /* function that takes care of choosing a infecting partner, when
   * the infection threshold has been reached, a partner must be chosen
   * based on chance and infectionrate
   * The function might fail (return false) when no partner turns
   * out to be infected.
   */
  bool infectByRandomPartner(Rng & , const std::vector<MhcLocus> & );
  bool infect(Rng & , const std::vector<MhcLocus> & , const Virus & );
  // if infect is successful, infectedFlag is made true
protected:
	Id id; // constant ID
	/* an Agent already has some data that will be passed to the Host
	 * object at infection. However, in some cases data is needed to
	 * determine if infection is possible.
	 */
	IndPars ipars; // used to determine if infection with a specific virus is possible
	Host* host; // pointer the the Host object. NULL if not infected
  double age; // age of agent (gets updated during simulation)
	double timeOfNaturalDeath;
	double timeOfBirth;
	// the gender of the agent can be used to model different populations
	Gender gender; // MALE or FEMALE
	// max number of contacts
	int maxNumContacts;
  // contact formation rate, load and threshold
  double cfRate, cfLoad, cfThreshold;
  // infection load and threshold (rate depends on neighbors)
  double infectionLoad, infectionThreshold;
	// list of contacts
	std::list<Contact*> contacts;
	// auxiliary function for shuffling the order of the contacts
	void shuffleContacts(Rng & );
};

// function for sorting Agents by age and other traits.
bool compareAgentByAge(Agent* , Agent* );
bool compareAgentByDegree(Agent* , Agent* );
bool compareAgentByTmrDegree(Agent* , Agent* );

/* typedef for Agent traits. Used for computing Agent statistics
 * cf. HostTrait. TODO: perhaps have one template trait for both Host and
 * Agent?
 */
template <typename T>
using AgentTrait = T(Agent::*)() const;

// some traits take additional arguments
template <typename T, typename R>
using AgentArgTrait = T(Agent::*)(R ) const;

// some traits are not always defined
template <typename Numeric, typename Error>
using AgentMaybeTrait = std::pair<Numeric, Error>(Agent::*)() const;



class Network {
	/** a class that can be inherited by World,
	 * and deals with contact formation.
	 */
public:
	Network() { /* empty */ }
	virtual ~Network() { /* empty */ } // TODO: perhaps warn when population.size() > 0
	double getMeanDegree(bool transmission_routes=false) const; // TODO: general function
	double getSdDegree(bool transmission_routes=false) const; // TODO: general function
	// statistics based on AgentTrait...
	template<class Numeric>
	double getMeanAgentTrait(AgentTrait<Numeric> ) const;
	// FIXME: ArgType1 and ArgType2 should be identical! (but e.g. T and const T & are in conflict)
	template<class Numeric, class ArgType1, class ArgType2>
	double getMeanAgentArgTrait(AgentArgTrait<Numeric, ArgType1> , ArgType2 ) const;
	template<class Numeric>
	std::string xmlNodeAgentTrait(AgentTrait<Numeric> , std::string ) const;
	// FIXME: ArgType1 and ArgType2 should be identical! (but e.g. T and const T & are in conflict)
	template<class Numeric, class ArgType1, class ArgType2>
	std::string xmlNodeAgentArgTrait(AgentArgTrait<Numeric, ArgType1> , ArgType2 , std::string ) const;
	// special Agent stats
	std::map<int, int> getDegreeDistribution() const;
	double getGlobalClusterCoeff() const;
	int getPopulationSize() const;
	// print the current contact structure in the dot format
	void printGraph(std::ostream & ) const;
protected:
	std::list<Agent*> population;
	// re-shuffle the population list (e.g. used when logging)
	void shufflePopulation(Rng & rng);
	/* makeContacts: pass Agents available for contacts,
	 * correction for the length of the time interval (?).
	 */
	void makeContacts(double t, double dt, Rng & rng); // TODO: use or remove dt
	/* some contacts can act as transmission routes (e.g. sexual contacts)
	 * after we make contacts, we have to make some of them transmission routes.
	 * makeTransmissionRoutes returns the number of enabled routes
	 */
	int makeTransmissionRoutes(Rng & rng); // dt
	// method to clear the population (deletes Agents)
	void clearPopulation();
private:
	/* TODO: it might be a good idea to make a private
	 * std::vector<Agent*> pop_vec, that can be used in makeContacts.
	 * Although, pop_vec could become invalid in between calls to
	 * makeContacts, which is unsafe...
	 */
};

/* handle templated Agent traits */

template<class Numeric>
double Network::getMeanAgentTrait(AgentTrait<Numeric> X) const {
	double xbar = 0.0;
	for ( Agent* agent : population ) {
		if ( agent == nullptr ) {
			throw std::logic_error("Agent in population is NULL" + RIGHT_HERE);
		} // else...
		xbar += (agent->*X)();
	}
	if ( population.size() > 0 ) return xbar/population.size();
	else return 0.0;
}

template<class Numeric, class ArgType1, class ArgType2>
double Network::getMeanAgentArgTrait(AgentArgTrait<Numeric, ArgType1> X, ArgType2 arg) const {
	double xbar = 0.0;
	for ( Agent* agent : population ) {
		if ( agent == nullptr ) {
			throw std::logic_error("Agent in population is NULL" + RIGHT_HERE);
		} // else...
		xbar += (agent->*X)(arg);
	}
	if ( population.size() > 0 ) return xbar/population.size();
	else return 0.0;
}

template<class Numeric>
std::string Network::xmlNodeAgentTrait(AgentTrait<Numeric> X, std::string name) const {
	BasicStats sts(name);
	for ( Agent* agent : population ) {
		if ( agent == nullptr ) {
			throw std::logic_error("Agent in population is NULL" + RIGHT_HERE);
		}
		double x = (agent->*X)(); // convert Numeric to double
		sts.addPoint(x);
	}
	sts.computeStats();
	std::stringstream ss; ss << sts;
	return ss.str();
}

template<class Numeric, class ArgType1, class ArgType2>
std::string Network::xmlNodeAgentArgTrait(AgentArgTrait<Numeric, ArgType1> X, ArgType2 arg, std::string name) const {
	BasicStats sts(name);
	for ( Agent* agent : population ) {
		if ( agent == nullptr ) {
			throw std::logic_error("Agent in population is NULL" + RIGHT_HERE);
		}
		double x = (agent->*X)(arg); // convert Numeric to double
		sts.addPoint(x);
	}
	sts.computeStats();
	std::stringstream ss; ss << sts;
	return ss.str();
}

#endif
