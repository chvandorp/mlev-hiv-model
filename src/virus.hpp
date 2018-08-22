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

#ifndef VIRUS_HPP_
#define VIRUS_HPP_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include <cmath> // abs
#include <stdexcept>

#include "macros.hpp"
#include "rng.hpp"
#include "calc_stats.hpp" // for computing a good H0
#include "frasers_rates.hpp"
#include "hashed_definitions.hpp"
#include "parallelism.hpp"
#include "aux.hpp" // logit, logit_inv
#include "ilevel_equations.hpp" // IndPars

class FitnessFunction; // forward declaration for Virus
class Virus; // forward declaration for hashing objects

// start and stop to compute fitness of 'ranges' of the genome
typedef std::pair<int, int> Range;
#define FULL_RANGE Range(0, GENOME_SIZE)

inline bool inRange(int loc, const Range & range) {
	return ( loc >= range.first && loc < range.second );
}
inline int len(const Range & range) {
	return range.second - range.first;
}

// the genes have different roles in the i-level model
/* TODO
 * 1) give genes a random order...
 * 2) select a subset of genes
 */

enum GeneLabel {
	BETA_GENE=0, // scale beta with w_beta
	DELTA_GENE, // death rate of producing cells
	VIRPROD_GENE, // production of virions
	RT_GENE, // vary the mutation rate
	VIRCLEAR_GENE, // clearence of virions
	DOWNREG_GENE, // cf. Nef. scale number of producing cells (I2) with 1/w_downreg
	ECLIPSE_GENE, // how fast do infected cells become virion producing cells?
	ABORT_GENE, // influences the probability of aborted infection
	ACTIV_GENE, // influence target cell activation
	NUMBER_OF_GENES // the last element determines the number of genes
};

// the fitness parameters have names for printing

extern const std::map<int, std::string> geneNameMap;
// can't use map::operator[] on a const map. Instead: use map::at()

inline std::string getGeneColor(GeneLabel gene) {
	return GENE_COLOR_MAP + std::to_string(1 + int(gene));
}

// definition of the Virus class...

/* TODO:
 * enable computation of (Malthusian) fitness,
 * possibly all host rates.
 */


// hasher and equality objects used by CiScheme
class VirusPtrHash {
public:
  size_t operator()(Virus* ) const;
private:
  std::hash<std::vector<bool>> hasher;
};

class VirusPtrEq {
public:
  bool operator()(Virus* , Virus* ) const;
private:
};

// TODO: should I assign a Id to DeterVirusVertex? (or VirusVertex, using ANON_ID by default)

class Virus : public Printable {
public:
	// TODO: add restricted
  enum InitType {
		ZERO=0, // a genome consisting of 0s
		SCRAMBLE, // 0 or 1 alleles with probability 1/2
		EVOLVE, // first scramble, then evolve by choosing the best mutation repeatedly
		SUBOPT, // first scramble, then evolve by choosing random beneficial mutations repeatedly
		METROPOLIS, // first scramble, then evolve using the Metropolis-Hastings algorithm
		ANNEALING // first scramble, then evolve using the Simulated Annealing algorithm
	};
  // constructors
  Virus(const FitnessFunction* , Rng & , InitType=SCRAMBLE);
	Virus(const FitnessFunction* , const Virus & ); // copy virus genome, but use other fitness function
	Virus(const FitnessFunction* , const std::list<Virus> & ); // consensus constructor
  Virus(const Virus & , int); // mutant constructor (pass locus)
  Virus(const Virus & , Rng & , int ); // mutant constructor (number of random mutations)
  Virus(const Virus & , Rng & , double ); // mutant constructor (probability of mutation)
  // TODO: more variations: evolve from given virus etc...
	bool operator[](int ) const; // modulo GENOME_SIZE
	bool operator==(const Virus & ) const; // equality of genomes
	bool operator!=(const Virus & ) const; // inequality of genomes
	int hammingDistance(const Virus & ) const;
	/* find the first locus where this differs from argument,
	 * right return value is false when viruses are identical
	 */
	std::pair<int, bool> findDifference(const Virus & ) const;
	// find all loci where this differs from argument
	std::list<int> findAllDifferences(const Virus & ) const;
	double getW() const; // "global" fitness measure
	double getWMutant(int ) const;
  double getWGene(int ) const;
  double getMalthusian() const; // TODO: implement
  double getMalthusian(const IndPars & ) const;  // TODO: implement
	// artificial traits
	double getVlMultiplier() const;
	void mutateVlMultiplier(Rng & ); // TODO: better method?
  // TODO: getMutant etc. redundant?
	Virus getMutant(int ) const; // return point mutant at <arguments> position
  Virus getMutant(Rng & ) const; // return a random point-mutant
	Virus* getMutantPtr(int ) const; // idem, but then return a pointer to a newly allocated Virus
  Virus* getMutantPtr(Rng & ) const; // return a pointer to a random point-mutant
	// TODO: fitness should not refer to w anymore...
	double getMeanFitnessPointMutants() const; // TODO: segments
	double getSdFitnessPointMutants() const;
	double getSkewFitnessPointMutants() const;
	// get and set Id
	Id getLineageId() const;
	void setLineageId(Id ); // TODO: this function violates the "a Virus is immutable idea"
	// print function
	void print(std::ostream & ) const;
	// friend classes
  friend class VirusPtrHash; // for genome-based hashing
  friend class VirusPtrEq;
protected:
// auxiliary functions (not to be used from the outside)
  void scramble(Rng & );
  void evolve(); // at each st, choose the highest increase
  void evolve_subopt(Rng & ); // at each step, choose a random increase
	// don't allow certain mutations (simulating escapes, TODO)
  void evolve_restricted(Rng & , const std::vector<int> & );
  void evolve_metropolis(Rng & , int ); // Metropolis algorithm
  void evolve_annealing(Rng & , int ); // Simulated Annealing
  void calc_ws();
	// data
	std::vector<bool> genome;
	double w_cache;
	double w_caches[NUMBER_OF_GENES];
	// a pointer to the common fitness function(s)
	const FitnessFunction* ffun; // pointer to const FitnessFunction
	// artificial traits for testing the Shireff model
	double vl_multiplier;
	// a generation counter, used to prevent w-corruption
	long int generation_counter;
	Id lineageId; // inherited from parent
};

typedef double (Virus::*VirusTrait)() const; // e.g. getFitness

int hammingDistance(const Virus & , const Virus & );
double correlationLength(const Virus & , Rng & );

class VirusJob : public Job {
	// wrapper for Virus to evolve in parallel
public:
	VirusJob(FitnessFunction* ffun) : vir(nullptr), ffun(ffun) {}
	~VirusJob();
	Virus getVirus() const;
protected:
	Virus* vir;
	FitnessFunction* ffun;
	bool execute(Rng & );  // override Job::execute
};

/* definition of fitnessFunction and FitnessContribution */

class FitnessContribution : public Printable {
public:
	FitnessContribution();
	void init(int , Rng & );
	double operator()(const Virus & ) const; // compute fitness
	double operator()(const Virus & , int ) const; // compute fitness, given mutation
	void addEdgeOut(const FitnessContribution & fc);
	void print(std::ostream & ) const;
  std::vector<int> getNeighbors() const;
  std::vector<int> getEdgesOut() const;
  int getPos() const;
	double getWeight() const;
protected:
	std::vector<int> neighbors;
	std::vector<double> values;
	std::vector<int> edgesOut; // loci that depend on this locus/allele
	int pos;
	int localDegree; // the number of neighbors for this particular locus
	double weight;
	// auxilliary functions
	int getCode(const Virus & ) const;
};

class FitnessFunction : public Printable {
public:
	FitnessFunction();
  FitnessFunction(Rng & );
	virtual ~FitnessFunction();
	void init(Rng & );
	// compute fitness over the whole range
	double operator()(const Virus & ) const;
  // compute fitness, given a mutation at a locus
	double operator()(const Virus & , int ) const;
  // segment-specific fitness
  double operator()(int , const Virus & ) const;
  // compute segment-specific fitness, given a mutation at a locus
	double operator()(int , const Virus & , int ) const;
  // "general purpose" function, used by the above
	double operator()(bool , int , const Virus & , bool , int ) const;
	// compute the fitness difference with a mutant, minimizing the number of computations
	double delta(const Virus & , int ) const;
	// compute a gene's fitness difference with a mutant, minimizing the number of computations
	double delta(int , const Virus & , int ) const;
	// "general purpose" delta function, used in the above
	double delta(bool , int , const Virus & , int ) const;
	Range getRange(int ) const;
	// auxiliary function to find in which gene a locus resides
	GeneLabel getGene(int loc) const;
	void print(std::ostream & ) const;
  void printGraph(std::ostream & ) const; // print a graph that can be used by GraphViz dot
protected:
	std::vector<FitnessContribution> fitnessContributions;
  Range geneRanges[NUMBER_OF_GENES];
	double H0Full; // normalization constant
	double H0s[NUMBER_OF_GENES]; // a H0 for each gene
	void determineH0(Rng & rng); // TODO multiple methods
};

// TODO: make these member functions of Virus

double iLevBasicReproductionNumber(const Virus & , const IndPars & );
double iLevBasicMalthusianFitness(const Virus & , const IndPars & );
std::pair<Vec, double> iLevInitialCondition(const Virus & , const IndPars & );

#endif
