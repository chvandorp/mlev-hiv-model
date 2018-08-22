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

#ifndef IMMUNERESPONSE_HPP_
#define IMMUNERESPONSE_HPP_

#include <iostream>
#include <fstream>
#include <algorithm> // std::sort
#include <stdexcept>

#include "macros.hpp"
#include "aux.hpp" // Id
#include "virus.hpp"
#include "rng.hpp"
#include "frasers_rates.hpp"
#include "ilevel_equations.hpp"
#include "hashed_definitions.hpp"
#include "parameters.hpp"

class ImmuneResponse {
public:
	ImmuneResponse(); // default
	ImmuneResponse(Rng & ); // for testing...
	ImmuneResponse(int , bool , int , Rng & ); // explicit locus, allele, tcr_code
	// locus and allele equal (NB: a Mhc can only have one immune response per locus+allele
	bool operator==(const ImmuneResponse & ) const;
	bool operator!=(const ImmuneResponse & ) const;
	void setId(Id );
	Id getId() const;
	int getLocus() const;
	bool getAllele() const;
	int getTcrCode() const;
	double getH() const;
  double getK() const;
	double getDominance() const;
	bool hasEpitope(const Virus & ) const;
	bool hasEpitopeMutant(const Virus & , int ) const;
	bool isDominant() const;
	// aux
	void print(std::ostream & ) const;
	// static parameters
	const static double k0_dom_def, k_dom_def, k0_sub_def, k_sub_def;
	const static double h0_dom_def, h_dom_def, h0_sub_def, h_sub_def;
	const static double prob_dom; // expectation of dominance
protected:
	Id id;
	int locus; // the epitope number
	bool allele; // is 1 or 0 the epitope?
	int tcr_code;
	std::vector<bool> tcr_bits;
	// warning: undefined behaviour when reading beyond vector<bool> boundaries
	double dominance; // in (0,1), determines the quality of the response
	double h, k; // Michaelis-Menten parameter, affinity parameter TODO: rename h -> hE
	// auxiliary method used in constructor
	void sample_ode_parameters(Rng & rng);
};

// auxiliary functions and definitions...

std::vector<bool> mkTcrBits(int );

class CompareImmresp {
public:
	bool operator()(const ImmuneResponse & , const ImmuneResponse & ) const;
};
std::ostream & operator<<(std::ostream & , const ImmuneResponse & );
typedef std::vector<ImmuneResponse> ImmuneSystem;
std::ostream & operator<<(std::ostream & , const ImmuneSystem & );

typedef std::pair<int, int> MhcID;
// the class HlaAllele
class MhcAllele {
public:
	MhcAllele(); // set values to defaults
	std::vector<ImmuneResponse> immuneResponses;
	double frequency;
	int locus; // A, B, C,...
	int allele; // 1, 2, 57, 28,...
	// auxiliary
	void print(std::ostream & ) const;
	MhcID getMhcID() const; // return std::pair<int, int>(locus, allele)
	std::string xmlStringMhcId() const;
	int countEpitopes(const Virus & ) const;
	void countEpitopes(const Virus & , int & , double & ) const; // TODO: redundant? return pair?
};

std::ostream & operator<<(std::ostream & , const MhcAllele & );

typedef std::vector<MhcAllele> MhcLocus; // all possible alleles on a locus

// auxiliary function for creating immune responses
std::vector<int> samplePermittedEpitopeLoci(Rng & );
std::vector<MhcLocus> createImmuneResponses(Rng & );
std::vector<MhcLocus> evolveImmuneResponses(Rng & ); // TODO
// document immune responses...
std::string xmlStringMhcLoci(const std::vector<MhcLocus> & , std::string );
// 2nd argument is file id

// haplotypes used by Host
typedef std::vector<MhcID> MhcHaplotype;
// pair of Mhc allele, and virus locus
typedef std::pair<MhcID, int> AlleleLocusPair;
double jaccardMhcHaplo(const MhcHaplotype & , const MhcHaplotype & );

// sample immuneresponses from the population-level alleles
void sampleImmuneResponses(ImmuneSystem & , MhcHaplotype & ,
                           const std::vector<MhcLocus> & allMhcLoci, Rng & );

#endif
