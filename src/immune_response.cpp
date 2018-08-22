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

#include "immune_response.hpp"

// unscaled killing rate of effector cells
const double ImmuneResponse::k0_dom_def = 200.0; // TODO: good value?
// killing rate per effector per productively infected cell
const double ImmuneResponse::k_dom_def = ImmuneResponse::k0_dom_def / IndPars::T0_def;
// unscaled killing rate of effector cells
const double ImmuneResponse::k0_sub_def = 10.0;
// killing rate per effector per productively infected cell
const double ImmuneResponse::k_sub_def = ImmuneResponse::k0_sub_def / IndPars::T0_def;
// unscaled Michaelis-Menten parameter for effector proliferation
const double ImmuneResponse::h0_dom_def = 1e-5;
// Michaelis-Menten parameter for effector proliferation
const double ImmuneResponse::h_dom_def = ImmuneResponse::h0_dom_def * IndPars::T0_def;
// unscaled Michaelis-Menten parameter for effector proliferation
const double ImmuneResponse::h0_sub_def = 1e-1;
// Michaelis-Menten parameter for effector proliferation
const double ImmuneResponse::h_sub_def = ImmuneResponse::h0_sub_def * IndPars::T0_def;
// fraction of dominant responses
const double ImmuneResponse::prob_dom = 0.25; // TODO good value?


ImmuneResponse::ImmuneResponse() {
	locus = 0;
	allele = false;
	tcr_code = 0;
	tcr_bits = mkTcrBits(tcr_code);
	dominance = prob_dom;
	// interpolate the ODE parameters between max and min
	h = weighted_geometric_average(h_dom_def, h_sub_def, dominance);
	k = weighted_geometric_average(k_dom_def, k_sub_def, dominance);
	id = ANON_ID;
}

ImmuneResponse::ImmuneResponse(Rng & rng) {
	locus = rng.Integer(GENOME_SIZE);
	allele = rng.Bit();
	tcr_code = rng.Integer(1 << NUMBER_TCR_CONTACTS);
	tcr_bits = mkTcrBits(tcr_code);
	id = ANON_ID;
	sample_ode_parameters(rng); // determines h, k, dominant
}

ImmuneResponse::ImmuneResponse(int locus, bool allele, int tcr_code, Rng & rng) {
	this->locus = locus;
	this->allele = allele;
	this->tcr_code = tcr_code;
	tcr_bits = mkTcrBits(tcr_code);
	id = ANON_ID;
	sample_ode_parameters(rng); // determines h, k, dominant
}

bool ImmuneResponse::operator==(const ImmuneResponse & ir) const {
	return ( ir.locus == locus && ir.allele == allele && ir.tcr_code == tcr_code );
}

bool ImmuneResponse::operator!=(const ImmuneResponse & ir) const {
	return !(*this == ir);
}

void ImmuneResponse::setId(Id id) { this->id = id; }
Id ImmuneResponse::getId() const { return id; }
int ImmuneResponse::getLocus() const { return locus; }
bool ImmuneResponse::getAllele() const { return allele; }
int ImmuneResponse::getTcrCode() const { return tcr_code; }
double ImmuneResponse::getH() const { return h; } // the Michaelis-Menten parameter
double ImmuneResponse::getK() const { return k; } // the killing-rate
double ImmuneResponse::getDominance() const { return dominance; }

bool ImmuneResponse::hasEpitope(const Virus & virus) const {
	if ( virus[locus] == allele ) {
		for ( int j = 0; j < NUMBER_TCR_CONTACTS; ++j ) {
			// for epitopes, assume a circular genome
			int k = (locus+j+1) % GENOME_SIZE;
      if ( tcr_bits[j] != virus[k] ) return false;
		}
		return true;
	}
	else return false;
}

bool ImmuneResponse::hasEpitopeMutant(const Virus & virus, int i) const {
  WARN_UNTESTED_FUN
	if ( i == locus ) {
    if ( virus[locus] != allele ) { // the mutant would match...
      // check the TCR bits (same as above)
      for ( int j = 0; j < NUMBER_TCR_CONTACTS; ++j ) {
				// for epitopes, assume a circular genome
				int k = (locus+j+1) % GENOME_SIZE;
      	if ( tcr_bits[j] != virus[k] ) return false;
      }
      return true;
    } else {
			return false;
		}
  } else { // the mutation is not at the locus (but could be at the TCR bits)
    if ( virus[locus] == allele ) {
    	// check the TCR bits
    	for ( int j = 0; j < NUMBER_TCR_CONTACTS; ++j ) {
				// for epitopes, assume a circular genome
      	int k = (locus+j+1) % GENOME_SIZE;
      	if ( i != k ) { // the mutation is not at the TCR bit
        	if ( tcr_bits[j] != virus[k] ) return false;
      	} else { // the mutation IS at the TCR bit...
        	if ( tcr_bits[j] == virus[k] ) return false;
      	}
      	return true;
    	}
    	return true;
  	} else {
			return false; // no MHC binding
		}
  }
}

bool ImmuneResponse::isDominant() const {
	return dominance > 1.0 - prob_dom; // TODO: what is a good threshold?
}

void ImmuneResponse::print(std::ostream & os) const {
	os 	<< "<immune_response "
		<< "id='" << id << "' "
		<< "locus='" << locus << "' "
		<< "allele='" << allele << "' "
		<< "tcr_code='" << tcr_code << "' "
		<< "dominance='" << dominance << "' "
		<< "h='" << h << "' "
    << "k='" << k << "' "
    << ">" << std::endl;
    os << "</immune_response>"; // close tag.
}


void ImmuneResponse::sample_ode_parameters(Rng & rng) {
	dominance = rng.Beta(prob_dom, 1.0-prob_dom); // TODO: what is a good distribution?
	// interpolate the ODE parameters between max and min
	h = weighted_geometric_average(h_dom_def, h_sub_def, dominance);
	k = weighted_geometric_average(k_dom_def, k_sub_def, dominance);
	// add some noise
	double fr = 0.1; // TODO: make parameter
	h = rng.LogNormalMnSd(h, fr*h);
	k = rng.LogNormalMnSd(k, fr*k);
}


/************* auxiliary functions for ImmuneResponse *************/

std::vector<bool> mkTcrBits(int tcr_code) {
	std::vector<bool> tcr_bits;
	for ( int i = 0; i < NUMBER_TCR_CONTACTS; ++i ) {
		if ( tcr_code & 1 ) tcr_bits.push_back(true);
		else tcr_bits.push_back(false);
		tcr_code >>= 1;
  }
  return tcr_bits;
}

bool CompareImmresp::operator()(const ImmuneResponse & im1,
                                const ImmuneResponse & im2) const {
	// a smaller h is better... TODO: what about k...
	return ( im1.getDominance() < im2.getDominance() );
}

std::ostream & operator<<(std::ostream & os, const ImmuneResponse & ir) {
	ir.print(os);
	return os;
}

std::ostream & operator<<(std::ostream & os, const ImmuneSystem & imsys) {
	os << "<immune_system >" << std::endl;
	ImmuneSystem::const_iterator it = imsys.begin();
	for ( ; it != imsys.end(); ++it ) {
		os << *it << std::endl;
	}
	os << "</immune_system>";
	return os;
}

MhcAllele::MhcAllele() {
	locus = 0; allele = 0; frequency = 0.0;
}

int MhcAllele::countEpitopes(const Virus & virus) const {
  auto pred = [&](const ImmuneResponse & ir){ return ir.hasEpitope(virus); };
  return std::count_if(immuneResponses.begin(), immuneResponses.end(), pred);
}

void MhcAllele::countEpitopes(const Virus & virus, int & n, double & dwn) const {
	n = 0; dwn = 0.0;
	for ( auto & ir : immuneResponses ) {
		if ( ir.hasEpitope(virus) ) {
			n += 1; dwn += ir.getDominance();
		}
	}
	// no return: return values passed by ref.
}

void MhcAllele::print(std::ostream & os) const {
	os 	<< "<mhc_allele "
		<< "locus='" << locus << "' "
		<< "allele='" << allele << "' "
		<< "frequency='" << frequency << "' "
		<< ">" << std::endl;
	for ( auto it = immuneResponses.begin(); it != immuneResponses.end(); ++it ) {
		os << *it << std::endl;
	}
	os << "</mhc_allele>";
}
std::ostream & operator<<(std::ostream & os, const MhcAllele & mhc) {
	mhc.print(os); return os;
}

std::string MhcAllele::xmlStringMhcId() const {
	std::stringstream ss;
	ss << "<mhc_id locus='" << locus << "' "
	   << "allele='" << allele << "' "
	   << "/>";
	return ss.str();
}

MhcID MhcAllele::getMhcID() const {
	return std::make_pair(locus, allele);
}

std::vector<int> samplePermittedEpitopeLoci(Rng & rng) {
	// make a vector of (unique) loci that are allowed to contain epitopes
	int maxNumEpitopeLoci = int(GENOME_SIZE * MAX_FRACTION_EPITOPE_LOCI);
	// check that the chosen parameters do not cause errors
	if ( maxNumEpitopeLoci < NUMBER_MHC_RESPONSES ) {
		throw std::runtime_error("not enough loci are allowed to contain epitopes" + RIGHT_HERE);
	}
	if ( maxNumEpitopeLoci > GENOME_SIZE ) {
		throw std::runtime_error("too many loci are allowed to contain epitopes" + RIGHT_HERE);
	}

  // sample maxNumEpitopeLoci from all loci on the genome
  std::vector<int> permittedEpitopeLoci(GENOME_SIZE, 0);
	rng.Shuffle(permittedEpitopeLoci.data(), 0, GENOME_SIZE);
  permittedEpitopeLoci.resize(maxNumEpitopeLoci);
	return permittedEpitopeLoci;
}

std::vector<MhcLocus> createImmuneResponses(Rng & rng) {
	std::vector<MhcLocus> allMhcLoci; // to be returned
	allMhcLoci.resize(NUMBER_MHC_LOCI);

	std::vector<int> permittedEpitopeLoci = samplePermittedEpitopeLoci(rng);

	// loop over loci
	int locus_id = 0;
	Id response_id = 0; // an id to uniquely identify immune responses...
	for ( auto lit = allMhcLoci.begin(); lit != allMhcLoci.end(); ++lit, ++locus_id ) {
		// make a frequency distribution
		std::vector<double> frequencies(NUMBER_MHC_ALLELES, 0.0);
    std::vector<double> alphas(NUMBER_MHC_ALLELES, 1.0); // parameter for Dirichlet distribution
    rng.Dirichlet(frequencies.data(), NUMBER_MHC_ALLELES, alphas.data());
		lit->resize(NUMBER_MHC_ALLELES);
		// loop over alleles
		int allele_number = 0;
		for ( auto ait = lit->begin(); ait != lit->end(); ++ait, ++allele_number ) {
			ait->locus = locus_id;
			ait->allele = allele_number;
			ait->frequency = frequencies[allele_number];
			// now select random loci from permittedEpitopeLoci
			std::vector<int> permuted_permitted_loci_idxs(permittedEpitopeLoci.size(), 0);
			rng.Shuffle(permuted_permitted_loci_idxs.data(), 0, permittedEpitopeLoci.size());
      // determine the number of responses for the allele
      double fr = CV_MHC_RESPONSES;
			double num = NUMBER_MHC_RESPONSES;
			int num_mhc_resps = 0;
			while ( num_mhc_resps <= 0 || num_mhc_resps >= int(permittedEpitopeLoci.size())) {
				// don't allow for MHC without any, or with to many responses
      	num_mhc_resps = rng.SymSkellam(num, num*fr);
			}
			// get random indices from permuted_permitted_loci_idxs
			for ( int i = 0; i < num_mhc_resps; ++i ) {
				int idx = permuted_permitted_loci_idxs[i];
				int eplocus = permittedEpitopeLoci[idx]; // select a locus on the genome
				// NB: per MHC molecule, we only permit one epitope allele (i.e. MHC binding)
				bool epallele = rng.Bit(); // choose a random allele
				// the number of TCR contacts determine the number of epitope variants
				int num_variants = 1 << NUMBER_TCR_CONTACTS;
				for ( int tcr_code = 0; tcr_code < num_variants; ++tcr_code ) {
					// create the immune response
					ImmuneResponse immresp(eplocus, epallele, tcr_code, rng);
					immresp.setId(response_id++);
					// and add it to the MHC allele
					ait->immuneResponses.push_back(immresp);
				}
			}
		}
	}
	return allMhcLoci;
}

std::vector<MhcLocus> evolveImmuneResponses(Rng & rng) { // TODO finish this!
	WARN_UNTESTED_FUN
	std::vector<MhcLocus> allMhcLoci; // to be returned
	allMhcLoci.resize(NUMBER_MHC_LOCI);
	std::vector<int> permittedEpitopeLoci = samplePermittedEpitopeLoci(rng);

	// evolve MHC alleles for each locus separately
	for ( auto lit = allMhcLoci.begin(); lit != allMhcLoci.end(); ++lit ) {
		// create an ancesteral MHC allele
		lit->resize(NUMBER_MHC_ALLELES);
		// TODO
	}

	return allMhcLoci;
}



std::string xmlStringMhcLoci(const std::vector<MhcLocus> & mhcLoci, std::string version) {
	std::stringstream ss;
	// open the xml node
	ss << "<mhc_loci version='" << version << "' >" << std::endl;
	ss << xmlStringParameters() << std::endl;

	for ( auto lit = mhcLoci.begin(); lit != mhcLoci.end(); ++lit ) {
		if ( lit->size() > 0 ) {
			// get the locus id
			int locus_id = lit->front().locus;
			// open a tag for the locus
			ss << "<mhc_locus locus='" << locus_id << "' >" << std::endl;
			for ( auto ait = lit->begin(); ait != lit->end(); ++ait ) {
				// document the allele
				ss << *ait << std::endl;
			}
			// closing tag for locus
			ss << "</mhc_locus>" << std::endl;
		}
		else {
			std::cerr << "# WARNING: MHC locus with zero HLA alleles" << std::endl;
		}
	}
	// closing tag
	ss << "</mhc_loci>" << std::endl;
	return ss.str();
}

double jaccardMhcHaplo(const MhcHaplotype & hapl1, const MhcHaplotype & hapl2) {
	int meet = 0;
	MhcHaplotype::const_iterator it1 = hapl1.begin();
	for ( ; it1 != hapl1.end(); ++it1 ) {
		MhcHaplotype::const_iterator it2 = std::find(hapl2.begin(), hapl2.end(), *it1);
		if ( it2 != hapl2.end() ) {
			++meet;
		}
	}
	int join = hapl1.size() + hapl2.size() - meet;
	if ( join > 0 ) {
		return double(meet) / join;
	}
	else {
		std::cerr << "WARNING: ill defined Jaccard index" << std::endl;
		return 1.0;
	}
}

void sampleImmuneResponses(ImmuneSystem & immuneResponses, MhcHaplotype & mhcHaplo,
                           const std::vector<MhcLocus> & allMhcLoci, Rng & rng) {
	immuneResponses.clear();
	immuneResponses.reserve(PLOIDY * NUMBER_MHC_LOCI * NUMBER_MHC_RESPONSES * (1 << NUMBER_TCR_CONTACTS));
	// may be too much (overlap MHC binding, homozygous individuals)
	mhcHaplo.clear();
	mhcHaplo.reserve(PLOIDY * NUMBER_MHC_LOCI);
	// may be too much.. (for homozygous individuals)

	// sample a PLOIDY amount of MhcAlleles from each locus
	int loc_idx = 0;
	for ( auto locus = allMhcLoci.begin(); locus != allMhcLoci.end(); ++locus, ++loc_idx ) {
		// sample \leq PLOIDY alleles_indices (don't add multiples)
		std::list<int> allele_indices;
		for ( int i = 0; i < PLOIDY; ++i ) {
			double u = rng.Uniform();
			// find the right allele index (cumulative frequency close to u)
			int allele_idx = 0;
			for ( allele_idx = 0; allele_idx < NUMBER_MHC_ALLELES; ++allele_idx ) {
				u -= (*locus)[allele_idx].frequency;
				if ( u < 0.0 ) break;
			}
			if ( allele_idx == NUMBER_MHC_ALLELES ) {
				std::cerr << "# WARNING: something went wrong while sampling from MhcAlleles"
				          << RIGHT_HERE << std::endl;
				allele_idx--; // failsafe...
			}
			allele_indices.push_back(allele_idx);
		}
		allele_indices.unique(); // remove any duplicates
		for ( auto it = allele_indices.begin(); it != allele_indices.end(); ++it ) {
			int allele_idx = *it;
			/* add (*locus)[idx].immuneResponses to immuneResponses.
			 * only add non-duplicate responses, replace duplicate when
			 * higher affinity
			 */

			// loop over new responses
			for (auto iit = (*locus)[allele_idx].immuneResponses.begin();
			     iit != (*locus)[allele_idx].immuneResponses.end(); ++iit ) {
				// find response in current immuneResponses
				auto jit = std::find(immuneResponses.begin(), immuneResponses.end(), *iit);
				if ( jit == immuneResponses.end() ) {
					immuneResponses.push_back(*iit); // add the new response
				}
				else { // the response was already added. replace with higher affinity?
					double h_old = jit->getH();
					double h_new = iit->getH();
					if ( h_new < h_old ) { // smaller h <-> higher affinity: replace!
						*jit = *iit;
					}
				}
			}
			/* for recording stats (and comparing hosts)
			 * mhcHaplo has no duplicate items.
			 */
			mhcHaplo.push_back(std::pair<int, int>(loc_idx, allele_idx));
		}
	}
	// arguments are passed by ref and monified...
}
