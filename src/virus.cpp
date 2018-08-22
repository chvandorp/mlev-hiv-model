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

#include "virus.hpp"

// initilize the (extern?) segment-name map

const std::map<int, std::string> geneNameMap = {
	std::make_pair(BETA_GENE, "w_beta"),
	std::make_pair(DELTA_GENE, "w_delta"),
	std::make_pair(VIRPROD_GENE, "w_virprod"),
	std::make_pair(RT_GENE, "w_rt"),
	std::make_pair(VIRCLEAR_GENE, "w_virclear"),
	std::make_pair(DOWNREG_GENE, "w_downreg"),
	std::make_pair(ECLIPSE_GENE, "w_eclipse"),
	std::make_pair(ABORT_GENE, "w_abort"),
	std::make_pair(ACTIV_GENE, "w_activ")
};

/********** members for FitnessContribution class ***********/

FitnessContribution::FitnessContribution() {
	pos = 0;
	localDegree = 0;
	weight = 1.0;
}

void FitnessContribution::init(int pos, Rng & rng) {
	neighbors.clear(); values.clear();
	this->pos = pos;
	neighbors.push_back(pos); // the first neighbor is self.
	// chose a random degree
	localDegree = rng.Integer(GENOME_DEGREE+1);
	// choose a weight
	weight = rng.Exponential(1.0);
	//weight = rng.ParetoMinusOne(2.0); // TESTING mean 1.0
	// get random neighboring positions on the genome
	double range = localDegree * GENOME_LOCALITY;
	if ( localDegree > 0 ) {
		int i = 0;
		while ( i < localDegree ) {
			int nbpos = rng.SymSkellam(pos, range);
			// make sure that nbpos is non-negative
      nbpos = mod(nbpos, GENOME_SIZE);
      if ( std::find(neighbors.begin(), neighbors.end(), nbpos) == neighbors.end() ) {
				neighbors.push_back(nbpos);
				++i;
			} // else keep trying...
		}
	}
	// assign fitness contributions to all of the possible combinations
	auto value_generator = [&]{return rng.Beta(0.5, 0.5)-0.5;}; // select a distribution...
	// TESTING: Beta(0.5, 0.5) puts more weight on extremes
	std::vector<double> possible_values;
	for ( int i = 0; i < KAUFFMAN_Q; ++i ) {
		double val = value_generator();
		possible_values.push_back(val);
	}
	int numberOfCombinations = 1 << (localDegree+1);
	for ( int i = 0; i < numberOfCombinations; ++i ) {
		double val = 0.0;
		if ( rng.Bernoulli(KAUFFMAN_P) ) { // set a fraction of the values to 0
			values.push_back(val);
		}	else {
			if ( possible_values.empty() ) { // unlimited number of values
				val = value_generator();
			}	else { // sample from a restricted number of values
				int j = rng.Integer(possible_values.size());
				val = possible_values[j];
			}
			values.push_back(weight*val);
		}
	}
}

double FitnessContribution::operator()(const Virus & virus) const {
	int code = 0;
	for ( int i = 0; i <= localDegree; ++i ) {
		code <<= 1;
		code += (virus[neighbors[i]] ? 1 : 0);
	}
	return values[code];
}

double FitnessContribution::operator()(const Virus & virus, int loc) const {
	int code = 0;
	for ( int i = 0; i <= localDegree; ++i ) {
		code <<= 1;
		if ( loc != neighbors[i] ) { // default
			code += (virus[neighbors[i]] ? 1 : 0);
		}	else { // change the code...
			code += (virus[neighbors[i]] ? 0 : 1);
		}
	}
	return values[code];
}

void FitnessContribution::addEdgeOut(const FitnessContribution & fc) {
	/** Add fc.pos to edgesOut iff pos is in fc.neigbors
	 * That is: if this is a neighbor of fc
	 */
	auto it = std::find(fc.neighbors.begin(), fc.neighbors.end(), pos);
	if ( it != fc.neighbors.end() ) {
		edgesOut.push_back(fc.pos);
	}
}

std::vector<int> FitnessContribution::getNeighbors() const {
    return neighbors; // NB: the first neignbor is self
}

std::vector<int> FitnessContribution::getEdgesOut() const {
    return edgesOut;
}


int FitnessContribution::getPos() const { return pos; }
double FitnessContribution::getWeight() const { return weight; }

void FitnessContribution::print(std::ostream & os) const {
	os << "<fcon "
     << "pos='" << pos << "' "
	   << "degree='" << localDegree << "' "
	   << "weight='" << weight << "' "
	   << ">" << std::endl;
	// print a list of neighbor nodes
	os << "<neighbors >";
	for ( auto it = neighbors.begin(); it != neighbors.end(); ++it ) {
		os << *it << " ";
	}
	os << "</neighbors>" << std::endl;
	// print a list of values
	os << "<values >";
	for ( auto it = values.begin(); it != values.end(); ++it ) {
		os << *it << " ";
	}
	os << "</values>" << std::endl;
	os << "</fcon>";
}

int FitnessContribution::getCode(const Virus & virus) const {
	int code = 0;
	for ( int i = 0; i <= localDegree; ++i ) {
		code <<= 1;
		code += (virus[neighbors[i]] ? 1 : 0);
	}
	return code;
}

/*
std::ostream & operator<<(std::ostream & os, const FitnessContribution & fctr) {
	fctr.print(os);
	return os;
}
*/


/********** members for FitnessFunction class ***********/

FitnessFunction::FitnessFunction() : H0Full(0.0) {
	std::for_each(H0s, H0s+NUMBER_OF_GENES, [](double & x){x = 0.0;});
}

FitnessFunction::FitnessFunction(Rng & rng) : H0Full(0.0) {
	std::for_each(H0s, H0s+NUMBER_OF_GENES, [](double & x){x = 0.0;});
	init(rng);
}


FitnessFunction::~FitnessFunction() {
	/* empty */
}

void FitnessFunction::init(Rng & rng) {
	fitnessContributions.clear();
	for ( int i = 0; i < GENOME_SIZE; ++i ) {
		FitnessContribution fctr;
		fctr.init(i, rng);
		fitnessContributions.push_back(fctr);
	}
	// find outward edges
	for ( auto it = fitnessContributions.begin(); it != fitnessContributions.end(); ++it ) {
		for ( auto jt = fitnessContributions.begin(); jt != fitnessContributions.end(); ++jt ) {
			if ( it != jt ) {
				it->addEdgeOut(*jt);
			}
		}
	}
  // pre-define genes
  for ( int gene = 0; gene < NUMBER_OF_GENES; ++gene ) {
    int start = (GENOME_SIZE / NUMBER_OF_GENES) * gene;
    int stop = (GENOME_SIZE / NUMBER_OF_GENES) * (gene + 1);
    geneRanges[gene] = Range(start, stop);
  }
	// find a good H0
	determineH0(rng);
}

double FitnessFunction::operator()(bool restrict_gene, int gene, const Virus & virus, bool mutate_loc, int loc) const {
	double H = 0.0;
	Range range = ( restrict_gene ? geneRanges[gene] : FULL_RANGE );
	double H0 = ( restrict_gene ? H0s[gene] : H0Full );
	auto start = fitnessContributions.begin() + range.first; // Range is a synonym for std::pair<int, int>
	auto stop = fitnessContributions.begin() + range.second;
	for ( auto it = start; it != stop; ++it ) {
		// depending on mutate_loc switch, use the 'real' contribution, or the 'mutated' one
		double Hi = ( mutate_loc ? (*it)(virus, loc) : (*it)(virus) );
		H += Hi;
	}
	H /= len(range); // /= sqrt(len(range)); // TODO normalize the variance?
	return exp(-(H - H0) / BOLTZMANN_TEMPERATURE);
}

double FitnessFunction::operator()(const Virus & virus) const {
	return operator()(false, 0, virus, false, 0); // final '0' can be replaced by anything
}

double FitnessFunction::operator()(const Virus & virus, int loc) const {
	return operator()(false, 0, virus, true, loc);
}

double FitnessFunction::operator()(int gene, const Virus & virus) const {
	return operator()(true, gene, virus, false, 0);
}

double FitnessFunction::operator()(int gene, const Virus & virus, int loc) const {
	return operator()(true, gene, virus, true, loc);
}


double FitnessFunction::delta(bool restrict_gene, int gene, const Virus & virus, int loc) const {
	/* we need to re-compute the fitness contributions of everyone
	 * with connections with loc. There are two cases: (1) loc is in the range
	 * or not (2).
	 * 1) we need the fitnessContribution of the locus loc,
	 * and all loci in the range, that have loc as a neighnor.
	 * 2) we just need the fitnessContribution of the loci in the range,
	 * that are neighbors to locus loc
	 */
	double DeltaH = 0.0;
	Range range = ( restrict_gene ? geneRanges[gene] : FULL_RANGE );
	auto nbs = fitnessContributions[loc].getEdgesOut();
	if ( inRange(loc, range) ) {
		DeltaH += fitnessContributions[loc](virus, loc);
		DeltaH -= fitnessContributions[loc](virus);
	}
	for ( auto it = nbs.begin(); it != nbs.end(); ++it ) {
		int nb_loc = *it;
		if ( inRange(nb_loc, range) ) {
			DeltaH += fitnessContributions[nb_loc](virus, loc);
			DeltaH -= fitnessContributions[nb_loc](virus);
		}
	}
	DeltaH /= len(range); // /= sqrt(len(range)); // TODO normalize the variance
	return exp(-DeltaH / BOLTZMANN_TEMPERATURE);
}

double FitnessFunction::delta(const Virus & virus, int loc) const {
	return delta(false, 0, virus, loc);
}

double FitnessFunction::delta(int gene, const Virus & virus, int loc) const {
	return delta(true, gene, virus, loc);
}

Range FitnessFunction::getRange(int gene) const {
	return geneRanges[gene];
}


void FitnessFunction::print(std::ostream & os) const {
	os << "<ffun H0='" << H0Full << "' ";
	for ( int gene = 0; gene < NUMBER_OF_GENES; ++gene ) {
		os << "H0_" << geneNameMap.at(gene) << "='" << H0s[gene] << "' ";
	}
	os << ">" << std::endl;
	for ( auto it = fitnessContributions.begin(); it != fitnessContributions.end(); ++it ) {
		os << *it << std::endl;
	}
	os << "</ffun>";
}

void FitnessFunction::printGraph(std::ostream & os) const {
	os << "digraph {" << std::endl;
	os << "graph [outputorder=edgesfirst];" << std::endl;
  // nodes
  for ( auto it = fitnessContributions.begin(); it != fitnessContributions.end(); ++it ) {
    int pos = it->getPos();
		GeneLabel gene = getGene(pos);
		std::string gene_color = getGeneColor(gene);
		double weight = it->getWeight();
    os << "\t" << pos << " ["
		   << "shape=\"point\" "
			 << "width=\"" << 0.1*sqrt(weight) << "\""
 			 << "fillcolor=\"" << gene_color << "\""
			 << "];" << std::endl;
  }
  // edges
  for ( auto it = fitnessContributions.begin(); it != fitnessContributions.end(); ++it ) {
    auto nbs = it->getNeighbors();
    int pos = it->getPos();
    os << "\t" << pos << " -> {";
    for ( auto jt = nbs.begin(); jt != nbs.end(); ++jt ) {
      if ( *jt != pos ) os << *jt << " ";
    }
    os << "} [color=gray arrowsize=0.5];" << std::endl;
  }
  os << "}";
  // TODO: use weight for color or size?
}

void FitnessFunction::determineH0(Rng & rng) {
	/* auxiliary method to determine the normalizing constant.
	 * Evolve a number of viruses, compute their w parameter
	 * and choose an H0 based on these ws
	 */
	std::cout << "determining fitness function parameters..." << std::endl;
	// initiate a worker pool to evolve Virus in parallel.
	WorkerPool wp;
	wp.initWorkerPool(NUMBER_OF_THREADS, rng);
	std::list<VirusJob*> joblist;
	for ( int i = 0; i < NUMBER_VIRUS_SAMPLES; ++i ) {
		VirusJob* job = new VirusJob(this);
		wp.addNewJob(job); // let the worker pool consruct the virus
		joblist.push_back(job);
	}
	wp.syncWorkerThreads(true); // verbose
	std::cout << " done!" << std::endl;
	// calculate statistics...
	BasicStats Hbars;
	BasicStats Hbarss[NUMBER_OF_GENES];
	for ( VirusJob* job : joblist ) {
		Virus virus = job->getVirus();
		double w = virus.getW();
		double Hbar = -log(w) * BOLTZMANN_TEMPERATURE; // average Hi
		Hbars.addPoint(Hbar);
		for ( int gene = 0; gene < NUMBER_OF_GENES; ++gene ) {
			double w = virus.getWGene(gene);
			double Hbar = -log(w) * BOLTZMANN_TEMPERATURE; // average Hi
			Hbarss[gene].addPoint(Hbar);
		}
	}
	Hbars.computeStats();
	for ( int gene = 0; gene < NUMBER_OF_GENES; ++gene ) {
		Hbarss[gene].computeStats();
	}
	// compute median Hbar and set H0
	H0Full = Hbars.getMedian();
	for ( int gene = 0; gene < NUMBER_OF_GENES; ++gene ) {
		H0s[gene] = Hbarss[gene].getMedian();
	}
	// cleanup
	for ( VirusJob* job : joblist ) {
		delete job;
	}
	wp.waitForWorkerPoolToExit();
}

GeneLabel FitnessFunction::getGene(int loc) const {
	// auxiliary function to find in which gene a locus resides
	for ( int gene = 0; gene < NUMBER_OF_GENES; ++gene ) {
		Range range = geneRanges[gene];
		if ( inRange(loc, range) ) {
			return GeneLabel(gene);
		}
	}
	throw std::out_of_range("argument 'loc' not in any gene range" + RIGHT_HERE);
}

/******** methods for Virus class *********/

Virus::Virus(const FitnessFunction* ffun, Rng & rng, InitType itp) : ffun(ffun) {
  if ( ffun == nullptr ) {
		throw std::invalid_argument("FitnessFunction is NULL" + RIGHT_HERE);
	}
  genome.resize(GENOME_SIZE, false); // the 0000...0 virus by default
  switch ( itp ) {
    case ZERO: {
      calc_ws();
      break;
    }
    case SCRAMBLE: {
      scramble(rng);
      break;
  	}
    case EVOLVE: {
      scramble(rng);
      evolve();
      break;
    }
    case SUBOPT: {
      scramble(rng);
      evolve_subopt(rng);
      break;
    }
    case METROPOLIS: {
			scramble(rng);
			evolve_metropolis(rng, EVOLVE_METROPOLIS_STEPS);
			break;
		}
    case ANNEALING: {
			scramble(rng);
			evolve_annealing(rng, EVOLVE_METROPOLIS_STEPS);
			break;
		}
    default: {
			throw std::logic_error("unknown Virus initialization type" + RIGHT_HERE);
			break; // redundant
		}
  } // switch InitType
	// by default, the VL multiplier is 1
	vl_multiplier = 1.0;
	// newly generated virus:
	generation_counter = 0;
	lineageId = ANON_ID;
}

Virus::Virus(const FitnessFunction* ffun, const Virus & virus) : ffun(ffun) {
	// copy virus genome, but use other fitness function
	if ( ffun == nullptr ) {
		throw std::invalid_argument("FitnessFunction is NULL" + RIGHT_HERE);
	}
	genome = virus.genome;
	calc_ws();
	vl_multiplier = virus.vl_multiplier;
	// increment the parents generation_counter
	generation_counter = virus.generation_counter + 1;
	// inherit the lineageId from the parent
	lineageId = virus.lineageId;
}

Virus::Virus(const FitnessFunction* ffun , const std::list<Virus> & viruses) {
	if ( ffun == nullptr ) {
		throw std::invalid_argument("FitnessFunction is NULL" + RIGHT_HERE);
	}
	genome.resize(GENOME_SIZE, false); // the 0000...0 virus by default

	this->ffun = ffun;
	std::vector<int> ones(GENOME_SIZE, 0.0);
	int n = viruses.size();
	if ( n > 0 ) {
		// count number of ones per locus
		for ( auto & virus : viruses ) {
			for ( int i = 0; i < GENOME_SIZE; ++i ) {
				ones[i] += ( virus[i] ? 1 : 0 );
			}
		}
		// make a consensus genome
		for ( int i = 0; i < GENOME_SIZE; ++i ) {
			double p = double(ones[i]) / n;
			if ( p > 0.5 ) { // slight bias towards 0000...0
				genome[i] = true;
			}
		}
	}
	calc_ws();
	vl_multiplier = 1.0; // TODO: take average of viruses in list?
	// newly generated virus:
	generation_counter = 0; // TODO: take average of viruses in list?
	lineageId = ANON_ID;
}

Virus::Virus(const Virus & virus, int loc) { // the MUTANT constructor :-)
	ffun = virus.ffun;
	genome = virus.genome; // copy the genome
	genome[loc] = !genome[loc]; // apply the mutation
	// increment the parents generation_counter
	generation_counter = virus.generation_counter + 1;
	if ( generation_counter % PREVENT_W_CORRUPTION_INTERVAL != 0 ) {
		// compute the fitness parameters, using virus
		w_cache = virus.w_cache * ffun->delta(virus, loc);
		for ( int gene = 0; gene < NUMBER_OF_GENES; ++gene ) {
			w_caches[gene] = virus.w_caches[gene] * ffun->delta(gene, virus, loc);
		}
	} else { // otherwise, compute the w parameters from scratch
		calc_ws();
	}
	// artificial traits
	vl_multiplier = virus.vl_multiplier;
	// inherit the lineage ID
	lineageId = virus.lineageId;
}

Virus::Virus(const Virus & virus, Rng & rng, int n) {
	ffun = virus.ffun;
	genome = virus.genome; // copy the genome
	for ( int i = 0; i < n; ++i ) {
		int loc = rng.Integer(genome.size()); // choose a random locus
		genome[loc] = !genome[loc]; // apply the mutation
	}
	// compute the fitness parameters
	calc_ws();
	// artificial traits
	vl_multiplier = virus.vl_multiplier;
	// increment the parents generation_counter
	generation_counter = virus.generation_counter + 1;
	// inherit the lineage ID
	lineageId = virus.lineageId;
}

Virus::Virus(const Virus & virus, Rng & rng, double p) {
	WARN_UNTESTED_FUN
	ffun = virus.ffun;
	genome = virus.genome; // copy the genome
	for ( size_t i = 0; i < genome.size(); ++i ) {
		if ( rng.Bernoulli(p) ) {
			genome[i] = !genome[i]; // apply a mutation with probability p
		}
	}
	// compute the fitness parameters
	calc_ws();
	// artificial traits
	vl_multiplier = virus.vl_multiplier;
	// increment the parents generation_counter
	generation_counter = virus.generation_counter + 1;
	// inherit the lineage ID
	lineageId = virus.lineageId;
}

double Virus::getW() const {
    return w_cache;
}
double Virus::getWMutant(int locus) const {
    return (*ffun)(*this, locus);
}
double Virus::getWGene(int gene) const {
    return w_caches[gene];
}

double Virus::getVlMultiplier() const {
	return vl_multiplier;
}
void Virus::mutateVlMultiplier(Rng & rng) {
	vl_multiplier *= pow(10.0, rng.Normal(0.0, SHIRREFF_SIGMA_M));
}

int Virus::hammingDistance(const Virus & v) const {
  if ( genome.size() != v.genome.size() ) {
    throw std::logic_error("genomes of different size" + RIGHT_HERE);
  }
	int h = 0;
	for ( size_t i = 0; i < genome.size(); ++i ) {
		if ( genome[i] != v.genome[i] ) {
			++h;
		}
	}
	return h;
}

std::pair<int, bool> Virus::findDifference(const Virus & v) const {
	// finds the first locus where this differs from argument
	if ( genome.size() != v.genome.size() ) {
    throw std::logic_error("genomes of different size" + RIGHT_HERE);
  }
	for ( size_t i = 0; i < genome.size(); ++i ) {
		if ( genome[i] != v.genome[i] ) {
			return std::make_pair(int(i), true);
		}
	}
	return std::make_pair(int(0), false);
}

std::list<int> Virus::findAllDifferences(const Virus & v) const {
	WARN_UNTESTED_FUN
	// find all loci where this differs from argument
	if ( genome.size() != v.genome.size() ) {
    throw std::logic_error("genomes of different size" + RIGHT_HERE);
  }
	std::list<int> diffs;
	for ( size_t i = 0; i < genome.size(); ++i ) {
		if ( genome[i] != v.genome[i] ) {
			diffs.push_back(i);
		}
	}
	return diffs;
}

bool Virus::operator[](int pos) const {
	if ( pos < 0 || pos >= int(genome.size()) ) {
		throw std::out_of_range("index out of bounds" + RIGHT_HERE);
	} // else...
	return genome[pos];
}

bool Virus::operator==(const Virus & vir) const {
	return ( genome == vir.genome );
}

bool Virus::operator!=(const Virus & vir) const {
	return !( *this == vir );
}

Virus Virus::getMutant(int loc) const {
	return Virus(*this, loc); // make a mutant copy of this
}

Virus Virus::getMutant(Rng & rng) const {
	int loc = rng.Integer(genome.size()); //
	return Virus(*this, loc); // make a mutant copy of this
}

Virus* Virus::getMutantPtr(int loc) const {
	return new Virus(*this, loc);
}

Virus* Virus::getMutantPtr(Rng & rng) const {
	int loc = rng.Integer(genome.size()); //
	return new Virus(*this, loc);
}

double Virus::getMeanFitnessPointMutants() const {
	double mw = 0.0;
    for ( size_t i = 0; i < genome.size(); ++i ) {
        mw += (*ffun)(*this, i);
    }
	if ( genome.size() > 0 ) mw /= genome.size();
	return mw;
}

double Virus::getSdFitnessPointMutants() const {
	double sdw = 0.0;
	double mw = getMeanFitnessPointMutants();
    for ( size_t i = 0; i < genome.size(); ++i ) {
        double x = (*ffun)(*this, i) - mw;
        sdw += x*x;
    }
	if ( genome.size() > 0 ) sdw /= genome.size();
	return sdw;
}

double Virus::getSkewFitnessPointMutants() const {
	double mu3 = 0.0;
	double mu1 = getMeanFitnessPointMutants();
	double sd = getSdFitnessPointMutants(); // note that this function also computes the mean...
	if ( sd > 0.0 && genome.size() > 0 ) {
    for ( size_t i = 0; i < genome.size(); ++i ) {
    	double cw = (*ffun)(*this, i) - mu1;
			mu3 += cw*cw*cw;
		}
		return mu3 / (sd*sd*sd * genome.size());
	}	else {
		// skewness undefined
		std::cerr << "WARNING: undefined skewness" << RIGHT_HERE << std::endl;
		return 0.0;
	}
}

Id Virus::getLineageId() const {
	return lineageId;
}
void Virus::setLineageId(Id id) {
	lineageId = id;
}


void Virus::evolve() {
	w_cache = (*ffun)(*this);
	while ( true ) {
    double w_max = 0.0; int i_max = 0;
    for ( size_t i = 0; i < genome.size(); ++i ) {
      double w = w_cache * ffun->delta(*this, i);
      if ( w_max < w ) {
        w_max = w;
        i_max = i;
      }
    }
  	if ( w_max > w_cache ) {
			genome[i_max] = !genome[i_max];
			w_cache = w_max;
		}
		else break;
	}
	calc_ws();
}

void Virus::evolve_subopt(Rng & rng) {
	std::vector<int> empty_vec;
	evolve_restricted(rng, empty_vec);
}

void Virus::evolve_metropolis(Rng & rng, int n) {
	w_cache = (*ffun)(*this);
	for ( int i = 0; i < n; ++i ) {
		int pos = rng.Integer(genome.size());
		double dw = ffun->delta(*this, pos);
		if ( dw > rng.Uniform() ) {
			genome[pos] = !genome[pos];
			w_cache *= dw;
		}
	}
	calc_ws();
}

void Virus::evolve_annealing(Rng & rng, int n) {
	w_cache = (*ffun)(*this);
	double theta = 1.0; // temperature
	for ( int i = 0; i < n; ++i ) {
		int pos = rng.Integer(genome.size());
		double dw = ffun->delta(*this, pos);
		double u = rng.Uniform();
		if ( dw > pow(u, theta) ) {
			genome[pos] = !genome[pos];
			w_cache *= dw;
		}
		// adjust temperature
		theta -= 1.0/n;
	}
	calc_ws();
}


void Virus::evolve_restricted(Rng & rng, const std::vector<int> & forbidden) {
	w_cache = (*ffun)(*this);
	while ( true ) {
		std::map<int, double> candidate_pos;
    for ( size_t i = 0; i < genome.size(); ++i ) {
      double w = w_cache * ffun->delta(*this, i);
			if ( w_cache < w && std::find(forbidden.begin(), forbidden.end(), i) == forbidden.end() ) {
      	candidate_pos[i] = w;
    	}
    }
   	if ( !candidate_pos.empty() ) {
			auto item = candidate_pos.begin();
			std::advance(item, rng.Integer(candidate_pos.size()));
			int i = item->first;
			genome[i] = !genome[i];
			w_cache = item->second;
		}	else {
			break;
		}
	}
	calc_ws();
}

void Virus::scramble(Rng & rng) {
	for ( size_t i = 0; i < genome.size(); ++i ) {
		genome[i] = rng.Bit();
	}
  calc_ws();
}

void Virus::calc_ws() {
	w_cache = (*ffun)(*this);
	for ( int gene = 0; gene < NUMBER_OF_GENES; ++gene ) {
		w_caches[gene] = (*ffun)(gene, *this);
	}
}

// auxiliary members:...

void Virus::print(std::ostream & os) const {
	os << "<virus "
  	 << "vl_multiplier='" << vl_multiplier << "' "
		 << "F='" << generation_counter << "' "
	   << "w='" << w_cache << "' ";
  for ( int gene = 0; gene < NUMBER_OF_GENES; ++gene ) {
		os << geneNameMap.at(gene) << "='" << w_caches[gene] << "' ";
	}
	os << ">" << std::endl;
	for ( auto it = genome.begin(); it != genome.end(); ++it ) {
		os << *it;
	}
	os << std::endl;
	os << "</virus>";
}

/*
std::ostream & operator<<(std::ostream & os, const Virus & virus) {
	virus.print(os);
	return os;
}
*/

int hammingDistance(const Virus & v1, const Virus & v2) {
	return v1.hammingDistance(v2);
}

double correlationLength(const Virus & vir, Rng & rng) {
	// do a random walk
	int n = GENOME_SIZE;
	std::vector<double> ws(n, 0.0);
	Virus mut = vir;
	for ( int i = 0; i < n; ++i ) {
  	ws[i] = mut.getW();
    mut = mut.getMutant(rng);
	}
	// compute autocorrelations...
	int s = sqrt(n); // quite arbitrary...
	double varw = 0.0; double meanw = 0.0;
	for ( int i = 0; i < n; ++i ) {
		meanw += ws[i];	varw += ws[i]*ws[i];
	}
	meanw /= n;	varw /= n; varw -= meanw*meanw;
	std::vector<double> corrs(s, 0.0);
	for ( int d = 1; d <= s; ++d ) {
		double m1 = 0.0; double m2 = 0.0; double sp = 0.0;
		for ( int i = 0; i <= n-d; ++i ) {
			m1 += ws[i];
			m2 += ws[i+d];
			sp += ws[i]*ws[i+d];
		}
		m1 /= (n-d+1); m2 /= (n-d+1); sp /= (n-d+1);
		double corr = (sp - m1*m2) / varw;
		corrs[d] = corr;
	}
	// estimate slope of autocorrelation
	double meancorr = 0.0; double sp = 0.0;
	double meand = (s+1)/2.0; double vard = (s*s-1)/3.0;
	for ( int d = 1; d <= s; ++d ) {
		double lcor = log(corrs[d]);
		meancorr += lcor;
		sp += d*lcor;
	}
	meancorr/=s; sp/=s;
	double slope = (sp - meand*meancorr) / vard;
	double corrlen = -1.0/slope; // dimension = [hamming distance]
	return corrlen;
}


double iLevBasicReproductionNumber(const Virus & virus, const IndPars & ipars) {
	double w_beta = virus.getWGene(BETA_GENE);
	double w_delta = virus.getWGene(DELTA_GENE);
	double w_virprod = virus.getWGene(VIRPROD_GENE);
	double w_virclear = virus.getWGene(VIRCLEAR_GENE);
	double w_rt = virus.getWGene(RT_GENE);
	double w_eclipse = virus.getWGene(ECLIPSE_GENE);
	double w_abort = virus.getWGene(ABORT_GENE);
	// we don't need w_downreg or w_activ
	return iLevBasicReproductionNumber(w_beta, w_delta, w_virprod, w_virclear, w_rt, w_eclipse, w_abort, ipars);
}

double iLevBasicMalthusianFitness(const Virus & virus, const IndPars & ipars) {
	double w_beta = virus.getWGene(BETA_GENE);
	double w_delta = virus.getWGene(DELTA_GENE);
	double w_virprod = virus.getWGene(VIRPROD_GENE);
	double w_virclear = virus.getWGene(VIRCLEAR_GENE);
	double w_rt = virus.getWGene(RT_GENE);
	double w_eclipse = virus.getWGene(ECLIPSE_GENE);
	double w_abort = virus.getWGene(ABORT_GENE);
	// we don't need w_downreg or w_activ
	return iLevBasicMalthusianFitness(w_beta, w_delta, w_virprod, w_virclear, w_rt, w_eclipse, w_abort, ipars);
}

std::pair<Vec, double> iLevInitialCondition(const Virus & virus, const IndPars & ipars) {
	double w_beta = virus.getWGene(BETA_GENE);
	double w_delta = virus.getWGene(DELTA_GENE);
	double w_virprod = virus.getWGene(VIRPROD_GENE);
	double w_virclear = virus.getWGene(VIRCLEAR_GENE);
	double w_rt = virus.getWGene(RT_GENE);
	double w_eclipse = virus.getWGene(ECLIPSE_GENE);
	double w_abort = virus.getWGene(ABORT_GENE);
	// we don't need w_downreg or w_activ
	return iLevInitialCondition(w_beta, w_delta, w_virprod, w_virclear, w_rt, w_eclipse, w_abort, ipars);
}


/* Virus hashing objects methods */

size_t VirusPtrHash::operator()(Virus* virus) const {
	if ( virus == nullptr ) {
		throw std::invalid_argument("trying to hash a NULL virus pointer" + RIGHT_HERE);
	}
  return hasher(virus->genome);
}

bool VirusPtrEq::operator()(Virus* virus1, Virus* virus2) const {
	if ( virus1 == nullptr || virus2 == nullptr ) {
		throw std::invalid_argument("trying to equate NULL virus pointers" + RIGHT_HERE);
	}
  return (virus1->genome == virus2->genome);
}

/* methods of VirusJob */

VirusJob::~VirusJob() {
	delete vir;
}

Virus VirusJob::getVirus() const {
	waitForJobToFinish();
	return *vir;
}

bool VirusJob::execute(Rng & rng) {// override Job::execute
	vir = new Virus(ffun, rng, Virus::EVOLVE); // TODO: options for init type
	return true;
}
