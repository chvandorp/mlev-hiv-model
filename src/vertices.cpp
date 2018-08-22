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

#include "vertices.hpp"


/** derivatives of the components of the state.
 * For safety, stick to the following rules:
 * 1) use the getValue(const double y[]) method of Vertex to get the correct
 *    value during integration. If the Vertex is not in the system of ODEs,
 *    the constant member value is returned, otherwise, y[index is returned
 * 2) ...
 */


void TargetVertex::deriv(const double y[], double dy[]) const {
	/** We have two target cell populations: Q and T
	 * dQ/dt = s - (d_Q + omega0) * Q - sum_i omega_i V_i / (h_Q + sum_i V_i)
	 * dT/dt = omega * Q - (d_T + sum_i beta_i * V_i) * T + sum_i omega_i V_i / (h_Q + sum_i V_i)
	 */
  // the current number of quiescent cells
  double Q = getValue(y, 0);
	// the current number of target cells
	double T = getValue(y, 1);
	// total virus load
	double Vtot = 0.0;
	// omega-weighted virus load
	double Vomegatot = 0.0;
	// const per-capita death rate
	double dTdt = omega0 * Q - d_T * T;
	// virus-induced activation and infection determined below
	double dQdt = s - (d_Q + omega0) * Q;
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		switch ( (*it)->getType() ) {
			case Vertex::DETER_VIRUS: {
				// cast the pointer to Vertex to VirusVertex
				DeterVirusVertex* vvert = (DeterVirusVertex*) (*it);
				double I2 = vvert->getValue(y, 1); // the second element is the I2 population
				// compute QSS for the virus load
				double V = getQssVirusLoad(vvert->d_V, vvert->p_V, vvert->beta, I2, T);
				// NB: DO NOT USE DeterVirusVertex::getQssVL here (or make a new overloaded method)
				dTdt -= vvert->beta * T * V; // infection of target cells
				// update variables for target cell activation
				Vtot += V;
				Vomegatot += vvert->omega * V;
				break;
			}
			default: {
				throw std::logic_error("ill-formed interaction graph" + RIGHT_HERE);
				break; // redundant
			}
		} // switch Vertex::getType()
	}
	// now compute the influx of (activated) target cells
	double virInducedActivation = Q * Vomegatot / (h_Q + Vtot);
	dQdt -= virInducedActivation;
	dTdt += virInducedActivation;
	assignDeriv(dy, dQdt, 0);
	assignDeriv(dy, dTdt, 1);
}

void DeterVirusVertex::deriv(const double y[], double dy[]) const {
	/* We have two infected cell populations: eclipse and productive
	 * dI1/dt = f * beta * pr_exact * V * T - (gamma + d_I1) * I1
	 * dI2/dt = gamma * I1 - d_I2 * I2 - m * I2 * sum_e k_e * E_e
	 */
	// current number of infected cells
	double I1 = getValue(y, 0);
	double I2 = getValue(y, 1);
	// per-capita death rate
	double dI1dt = -(d_I1 + gamma) * I1;
	double dI2dt = gamma * I1 - d_I2 * I2;
	double T = 0.0; // to be found later
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		Vertex* vert = (*it);
		switch ( vert->getType() ) {
			case Vertex::TARGET: { // find index and upate T
				/* the first element is the number of Quiescent cells.
				 * usually we only have one target cell population...
				 */
				T += vert->getValue(y, 1);
				break;
			}
			case Vertex::ACTIVE_RESPONSE: {
				double E = vert->getValue(y); // E is a scalar
				// cast pointer to get access to response-parameters
        ActiveResponseVertex* arv = (ActiveResponseVertex*) vert;
				dI2dt -= (arv->k) * m * I2 * E;
				break;
			}
			case Vertex::DETER_VIRUS: {
				// don't use "parents"
				break;
			}
			default: {
				throw std::logic_error("ill-formed interaction graph" + RIGHT_HERE);
				break; // redundant
			}
		} // switch
	}
  // compute QSS virus load
  double V = getQssVirusLoad(d_V, p_V, beta, I2, T);
  dI1dt += f * beta * pr_exact * V * T; // newly (and faithfully) infected cells
	assignDeriv(dy, dI1dt, 0);
	assignDeriv(dy, dI2dt, 1);
}


void StochVirusVertex::deriv(const double y[], double dy[]) const {
	double incidence = 0.0; // to be updated below if there are any "parents"
	double deathrate = d_I2; // to be updated below if there are any responses
	// we first need to find the number of target cells
	double T = 0.0;
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		Vertex* vert = (*it);
		if ( vert->getType() == Vertex::TARGET ) {
			/* the += is redundant: only one target population..
			 * the first component is the number of quiescent cells,
			 * the second the number of susceptible targets
			 */
			T += vert->getValue(y, 1);
			break; // and hence we can stop the for-loop
		}
	}
	// now we can compute incidence
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		Vertex* vert = (*it);
		switch ( vert->getType() ) {
			case Vertex::TARGET: {
				// we already handled the target cells above
				break;
			}
			case Vertex::ACTIVE_RESPONSE: {
				double E = vert->getValue(y);
				// cast pointer to get access to response-parameter
        ActiveResponseVertex* arv = (ActiveResponseVertex*) vert;
				deathrate += (arv->k) * m * E;
				break;
			}
			case Vertex::DETER_VIRUS: { // the "mother(s)"
				// cast pointer to get access to fitness parameter
				DeterVirusVertex* dvv = (DeterVirusVertex*) vert;
				// the second element of DeterVirusVertex is the number of productive cells
				double I2_mother = vert->getValue(y, 1);
				// compute QSS virus load
				double V_mother = getQssVirusLoad(dvv->d_V, dvv->p_V, dvv->beta, I2_mother, T);
				/* and add to incidence of THIS MUTANT!
				 * include factor dvv->pr_exact to approximate the probability
				 * of just one mutation
				 * TODO: pre compute mothers product of constants
				 * TODO: what would be the correct incidence? use eigenvalues?
				 * we now take the incidence of I1 cells:
				 */
				incidence += dvv->f * dvv->beta * dvv->pr_mut * dvv->pr_exact * V_mother * T;
				break;
			}
			default: {
				throw std::logic_error("ill-formed interaction graph" + RIGHT_HERE);
				break; // redundant
			}
		} // switch
	}
	// compute QSS virus load
	double VperI2 = getPcQssVirusLoad(d_V, p_V, beta, T);
	// now calculate the birthrate and reproduction number
	double birthrate = f * c * beta * pr_exact * VperI2 * T;
	double R = birthrate / deathrate;
	double dLdt = incidence * major_outbreak_prob(R);
  assignDeriv(dy, dLdt);
}

void ActiveResponseVertex::deriv(const double y[], double dy[]) const {
	/** dE/dt = p_E * E * sum_i m_i I2_i / (h + E + sum_i m_i I2_i) - d_E * E
	 * or:
	 * dE/dt = p_E * E * sum_i m_i I2_i / (h + E + \sum_i I2_i) - d_E * E
	 */
	double E = getValue(y);
	double dEdt = -d_E * E; // total death rate
	double I2tot = 0.0; // producing cells expressing the epitope
	double I2mtot = 0.0; // producing cells weighted by MHC expression
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		Vertex* vert = (*it);
		if ( vert->getType() == Vertex::DETER_VIRUS ) {
			// cast pointer to get access to virus-parameters
			DeterVirusVertex* dvv = (DeterVirusVertex*) vert;
			// the second element of a DeterVirusVertex is the number of producing cells
			double I2 = dvv->getValue(y, 1);
			// producing cells are weighted by how 'visible' they are
			I2tot += I2;
			I2mtot += I2 * dvv->m;
		}	else {
			throw std::logic_error("ill-formed interaction graph" + RIGHT_HERE);
		}
	}
	//dEdt += p_E * E * I2mtot / (h + E + I2mtot);
	dEdt += p_E * E * I2mtot / (h + E + I2tot); // TESTING
	assignDeriv(dy, dEdt);
}

void LatentResponseVertex::deriv(const double y[], double dy[]) const {
	/** a lentent response computes a cumulative rate (load) that needs
	 * to cross a threshold. The rate is a "activation rate" times
	 * the introduction probability 1 - 1/R.
	 * Hence, we need the p.c. birthrate and the p.c. deathrate of the
	 * corresponding active response.
	 */
	double I2tot = 0.0; // producing cells expressing the epitope
	double I2mtot = 0.0; // producing cells weighted by MHC expression
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		Vertex* vert = (*it);
		if ( vert->getType() == Vertex::DETER_VIRUS ) {
			// cast pointer to get access to virus parameters
			DeterVirusVertex* dvv = (DeterVirusVertex*) vert;
			// the second element of DeterVirusVertex is the number of producing cells
			double I2 = dvv->getValue(y, 1);
			// producing cells are weighted by how 'visible' they are
			I2tot += I2;
			I2mtot += I2 * dvv->m;
		}	else {
			throw std::logic_error("ill-formed interaction graph" + RIGHT_HERE);
		}
	}
	//double birthrate = p_E * I2mtot / (h + I2mtot); // because E = 0, and p.c.
	double birthrate = p_E * I2mtot / (h + I2tot); // because E = 0, and p.c. TESTING
	double R = birthrate / d_E;
	double dLdt = lambda * major_outbreak_prob(R);
	assignDeriv(dy, dLdt);
}

/* constructor for TargetVertex */

TargetVertex::TargetVertex(const IndPars* ip) : Vertex(TARGET, 2) {
	this->ip = ip; // TODO/FIXME: obsolete?
	d_T = ip->dT;
	d_Q = ip->dQ;
	s = ip->s;
	h_Q = ip->hQ;
	omega0 = ip->omega0;
	// auxiliary stuff
	inDiagram = false; // for drawing Graph in the dot language (false: don't include)
}

/* constructor for VirusVertex (protected) */

VirusVertex::VirusVertex(VertexType t, int size, Virus* vir, const IndPars* ip) :
		Vertex(t, size) {
	this->ip = ip; // TODO: obsolete?
  this->vir = vir;
  // define parameters
  beta = iLevBeta(vir->getWGene(BETA_GENE), *ip);
	d_I1 = ip->dI1;
  d_I2 = iLevDelta(vir->getWGene(DELTA_GENE), *ip);
  m = mhcExpression(vir->getWGene(DOWNREG_GENE), *ip);
  f = getFractionNonAbortive(vir->getWGene(ABORT_GENE), *ip);
  pr_exact = getExactReplicationProb(vir->getWGene(RT_GENE), *ip);
  pr_mut = getMutationProb(vir->getWGene(RT_GENE), *ip);
  d_V = virionDeathrate(vir->getWGene(VIRCLEAR_GENE), *ip);
  p_V = virionProduction(vir->getWGene(VIRPROD_GENE), *ip);
	gamma = getEclipseRate(vir->getWGene(ECLIPSE_GENE), *ip);
	omega = getActivationRate(vir->getWGene(ACTIV_GENE), *ip);
  // pre-computed compound parameters
	c = gamma / (d_I1 + gamma);
}

DeterVirusVertex::DeterVirusVertex(Virus* vir, const IndPars* ip) :
		VirusVertex(DETER_VIRUS, 2, vir, ip) {
	endangered = false;
	// auxiliary stuff
	inDiagram = true; // for drawing Graph in the dot language
}

StochVirusVertex::StochVirusVertex(Virus* vir, const IndPars* ip) :
		VirusVertex(STOCH_VIRUS, 1, vir, ip) {
	/* empty */
}


/* auxiliary methods for VirusVertex */

double VirusVertex::getR(bool malthusian, VirusVertex* mimic_resps) const {
	/* the reproduction number (is default) or Malthusian fitness.
	 * Notice that this function is largely identical to Vertex::deriv.
	 * However, it uses Vertex::getValue() and not the state vector.
	 * It is to be used in between ODE updates (driver_apply), to
	 * determine if a stochastic mutant should be added/removed
	 * from the system of ODEs.
	 */
	double deathrate = d_I2; // updated below if there are any responses
	// target cells to be found later...
	double T = 0.0;
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		switch ( (*it)->getType() ) {
			case Vertex::TARGET: { // update T
				/* the first element of TargetVertex is the number of quiescent cells,
				 * the second is the number of activated target cells
				 */
				T = (*it)->getValue(1);
				break;
			}
			case Vertex::ACTIVE_RESPONSE: {
				if ( mimic_resps == nullptr ) { // default: use own responses
        	ActiveResponseVertex* arv = (ActiveResponseVertex*) (*it);
					double E = arv->getValue();
					deathrate += (arv->k) * m * E;
				}
				break;
			}
			case Vertex::DETER_VIRUS: { // the "mother(s)"
				// don't use the parents.
				break;
			}
			default: {
				throw std::logic_error("ill-formed interaction graph" + RIGHT_HERE);
				break; // redundant
			}
		} // switch
	}
	// if a VirusVertex is passed, we have to use its immune responses instead
	if ( mimic_resps != nullptr ) {
		auto responses = mimic_resps->getResponses();
		for ( auto arv : responses ) {
			double E = arv->getValue();
			deathrate += (arv->k) * m * E;
		}
	}
	// calculate pc virus load
	double VperI2 = getPcQssVirusLoad(d_V, p_V, beta, T);
	double birthrate = f * beta * pr_exact * VperI2 * T;
	if ( malthusian ) {
		Mat M(-(gamma + d_I1), birthrate, gamma, -deathrate);
		auto eivals_real = eigen_values(M);
		if ( !eivals_real.second ) {
			throw std::logic_error("eigenvalues not real" + RIGHT_HERE);
		}
		Vec eivals = eivals_real.first;
		return std::max(eivals.first, eivals.second);
	} else {
		return c * birthrate / deathrate;
	}
}

std::list<ActiveResponseVertex*> VirusVertex::getResponses() const {
	std::list<ActiveResponseVertex*> resps;
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		if ( (*it)->getType() == Vertex::ACTIVE_RESPONSE ) {
			ActiveResponseVertex* arv = (ActiveResponseVertex*) (*it);
			resps.push_back(arv);
		}
	}
	return resps;
}

std::pair<int, double> VirusVertex::getNumberOfResponses() const {
	double deathrate = 0.0; // death due to immune responses
	int numberOfResponses = 0;
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		if ( (*it)->getType() == Vertex::ACTIVE_RESPONSE ) {
			ActiveResponseVertex* arv = (ActiveResponseVertex*) (*it);
			double E = arv->getValue();
			deathrate += (arv->k) * m * E;
			numberOfResponses++;
		}
	}
	return std::make_pair(numberOfResponses, deathrate);
}

std::list<DeterVirusVertex*> VirusVertex::getParents() const {
	std::list<DeterVirusVertex*> parents;
	// loop over EDGES IN!
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		Vertex* v = (*it);
		if ( v->getType() == Vertex::DETER_VIRUS ) {
			parents.push_back((DeterVirusVertex*) v);
		}
	}
	return parents;
}

std::list<StochVirusVertex*> VirusVertex::getChildren() const {
	std::list<StochVirusVertex*> children;
	// loop over EDGES OUT!
	for ( auto it = edgesOut.begin(); it != edgesOut.end(); ++it ) {
		Vertex* v = (*it);
		if ( v->getType() == Vertex::STOCH_VIRUS ) {
			children.push_back((StochVirusVertex*) v);
		}
	}
	return children;
}

std::list<VirusVertex*> VirusVertex::getAllChildren() const {
	std::list<VirusVertex*> children;
	// loop over EDGES OUT!
	for ( auto it = edgesOut.begin(); it != edgesOut.end(); ++it ) {
		Vertex* v = (*it);
		if ( v->getType() == Vertex::STOCH_VIRUS || v->getType() == Vertex::DETER_VIRUS ) {
			children.push_back((VirusVertex*) v);
		}
	}
	return children;
}

Virus* VirusVertex::getVirus() const {
	return vir;
}

int VirusVertex::getNumberOfParents() const {
	auto pred = [](Vertex* v){ return v->getType() == Vertex::DETER_VIRUS; };
	return std::count_if(edgesIn.begin(), edgesIn.end(), pred);
}

bool VirusVertex::isEscape() const {
	// find responses
	auto respPair = getNumberOfResponses(); // (n, m * sum_e k_e E_e)
	// find parents...
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		Vertex* v = (*it);
		if ( v->getType() == Vertex::DETER_VIRUS ) {
			DeterVirusVertex* parent = (DeterVirusVertex*) v;
			// count responses
			auto parentRespPair = parent->getNumberOfResponses();
			if ( parentRespPair.first > respPair.first ) {
				// escape w.r.t at least one of the parents
				return true;
			}
		}
	}
	// else, NOT an escape w.r.t any parent
	return false;
}

/* auxiliary methods for StochVirusVertex */

DeterVirusVertex* StochVirusVertex::mkDeterVirusVertex(Graph & G) const {
	/** This function makes a new DeterVirusVertex that gets
	 * the "right" connections. The responses, deter viruses (see below),
	 * and targets are informed about the new deterministic virus.
	 */
	DeterVirusVertex* v = new DeterVirusVertex(vir, ip);
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		Vertex* u = (*it);
		switch ( u->getType() ) {
			case Vertex::TARGET: {
				G.addUndirectedEdge(v, u);
				break;
			}
			case Vertex::ACTIVE_RESPONSE: {
				G.addUndirectedEdge(v, u);
				break;
			}
			case Vertex::DETER_VIRUS: {
				/* the parent edges can be used for:
				 * 1) after extinction, a stochastic copy is made that
				 * has to know its parents, therefore we remember these
				 * parents here.
				 * 2) once a strain has become deterministic, its parents
				 * are possible mutants (when they become extinct,
				 * see mkStochVirusVertex).
				 * Hence, the parents need an edge from this vertex.
				 * 3) TODO: when mutation rate is high, the effective
				 * infection rate decreases, but there is an influx of the
				 * mutants.
				 */
				G.addUndirectedEdge(v, u);
				break;
			}
			default: {
				throw std::logic_error("ill-formed interaction graph" + RIGHT_HERE);
				break; // redundant
			}
		} // switch
	} // for
	G.assignInOdeSys(v, true);
	// return the pointer
	return v;
}

/* auxiliary methods for DeterVirusVertex */

StochVirusVertex* DeterVirusVertex::mkStochVirusVertex(Graph & G) const {
	/** This function makes a new StochVirusVertex that gets
	 * the "right" connections. Typically called after extinction
	 * of the deterministic virus to allow for back-mutations (e.g.
	 * after a the disappearence of an immune response).
	 */
	StochVirusVertex* v = new StochVirusVertex(vir, ip);
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		Vertex* u = (*it);
		switch ( u->getType() ) {
			case Vertex::TARGET: {
				G.addEdge(u, v); // stoch depends on the target (but target not on stoch)
				break;
			}
			case Vertex::ACTIVE_RESPONSE: {
				G.addEdge(u, v); // stoch depends on response (but response not on stoch)
				break;
			}
			case Vertex::DETER_VIRUS: {
				/* the parent edges were copied in mkDeterVirusVertex to
				 * be used here: These will be the parents of the new stochastic
				 * virus.
				 */
				G.addEdge(u, v);
				break;
			}
			default: {
				throw std::logic_error("ill-formed interaction graph" + RIGHT_HERE);
				break; // redundant
			}
		} // switch
	} // for loop over edgesIn
	/* this function is typically called when a deterministic clone dies,
	 * hence inOdeSys = false
	 */
	G.assignInOdeSys(v, false);
	// return the pointer
	return v;
}

double DeterVirusVertex::assignInitialCondition(double inoculum) {
	if ( inoculum < 0 ) {
		throw std::invalid_argument("negative initial value" + RIGHT_HERE);
	}
	double B = 0.0;
	double delta = d_I2;
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		switch ( (*it)->getType() ) {
			case Vertex::TARGET: {
				double T = (*it)->getValue(1); // elements are 0: Q, 1: T
				double VperI2 = getPcQssVirusLoad(d_V, p_V, beta, T);
				B = f * pr_exact * beta * VperI2 * T;
				break;
			}
			case Vertex::ACTIVE_RESPONSE: {
				// cast to ActiveResponseVertex to get access to killing rate
				ActiveResponseVertex* arv = (ActiveResponseVertex*) (*it);
				double k = arv->k;
				double E = arv->getValue();
				delta += m * k * E;
				break;
			}
			case Vertex::DETER_VIRUS: {
				// don't need parents
				break;
			}
			default: {
				throw std::logic_error("ill-formed interaction graph" + RIGHT_HERE);
				break; // redundant
			}
		}
	}
	/* the linearized system is
	 * d/dt (I1, I2)' = (-(d_I1 + gamma), B; gamma, delta) (I1, I2)'
	 *                = M (I1, I2)'
	 */
	Mat M(-(d_I1 + gamma), B, gamma, -delta);
	auto eivals_real = eigen_values(M);
	if ( !eivals_real.second ) {
		throw std::logic_error("eigenvalues are not real" + RIGHT_HERE);
	}
	Vec eivals = eivals_real.first;
	double dom_eival = std::max(eivals.first, eivals.second);
	Vec dom_eivec = eigen_vector_sum1(M, dom_eival);
	/* make sure that the eigen vector is positive
	 * TODO: make a mustPositive function
	 */
	double u1 = dom_eivec.first;
	double u2 = dom_eivec.second;
	if ( u1 < 0 && u2 < 0 ) {
		u1 = -u1; u2 = -u2;
	} else if ( u1 <= 0 || u2 <= 0 ) {
		throw std::logic_error("eigenvector is not positive" + RIGHT_HERE);
	}
	assignValue(inoculum * u1, 0); // I1
	assignValue(inoculum * u2, 1); // I2
	return dom_eival;
}

double DeterVirusVertex::getQssVL(double T) const {
	double I2 = getValue(1); // element 0 would be I1 (eclipse)
	return getQssVirusLoad(d_V, p_V, beta, I2, T);
}

double DeterVirusVertex::getQssVL() const {
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		if ( (*it)->getType() == Vertex::TARGET ) {
			TargetVertex* tvert = (TargetVertex*) (*it);
			double T = tvert->getValue(1); // element 0 would be Q (quiescent)
			double I2 = getValue(1); // element 0 would be I1 (eclipse)
			return getQssVirusLoad(d_V, p_V, beta, I2, T);
		}
	} // for loop over edges in
	throw std::logic_error("ill-formed interaction graph" + RIGHT_HERE);
}

/* version of the method that can be used by the integrator.
 * TODO: is there a way to combine these overloaded functions?
 */
double DeterVirusVertex::getQssVL(const double y[]) const {
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		if ( (*it)->getType() == Vertex::TARGET ) {
			TargetVertex* tvert = (TargetVertex*) (*it);
			double T = tvert->getValue(y, 1); // element 0 would be Q (quiescent)
			double I2 = getValue(y, 1); // element 0 would be I1 (eclipse)
			return getQssVirusLoad(d_V, p_V, beta, I2, T);
		}
	} // for loop over edges in
	throw std::logic_error("ill-formed interaction graph" + RIGHT_HERE);
}


/* auxiliary methods for ResponseVertex */

// constructors

ResponseVertex::ResponseVertex(VertexType t, ImmuneResponse* resp,
			const IndPars* ip) : Vertex(t, 1) {
	this->resp = resp;
	this->ip = ip;
	h = resp->getH();
  k = resp->getK();
	p_E = ip->pE;
  d_E = ip->dE;
}

LatentResponseVertex::LatentResponseVertex(ImmuneResponse* resp,
			const IndPars* ip, bool memo) : ResponseVertex(LATENT_RESPONSE, resp, ip) {
	memory = memo;
	lambda = ( memo ? ip->lambda_memo : ip->lambda_naive );
}

ActiveResponseVertex::ActiveResponseVertex(ImmuneResponse* resp,
			const IndPars* ip) : ResponseVertex(ACTIVE_RESPONSE, resp, ip) {
	endangered = false;
	// auxiliary stuff
	inDiagram = true;
}

// other

double ResponseVertex::getR(bool malthusian) const {
	double I2tot = 0.0; // producing cells expressing the epitope
	double I2mtot = 0.0; // producing cells weighted by MHC expression
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		if ( (*it)->getType() == Vertex::DETER_VIRUS ) {
			DeterVirusVertex* dvv = (DeterVirusVertex*) (*it);
			double I2 = dvv->getValue(1); // element 0 would be I1 (eclipse)
			I2tot += I2;
			I2mtot += I2 * dvv->m; // producing cells are weighted by how 'visible' they are
		}	else {
			throw std::logic_error("ill-formed interaction graph" + RIGHT_HERE);
		}
	}
  double E = ( type == ACTIVE_RESPONSE ? getValue() : 0.0 ); // TODO: not so pretty...
	//double birthrate = p_E * I2mtot / (h + E + I2mtot);
	double birthrate = p_E * I2mtot / (h + E + I2tot); // TESTING
    if ( malthusian ) {
      return birthrate - d_E;
    } else {
      return birthrate / d_E;
    }
}

ImmuneResponse* ResponseVertex::getResponse() const {
	return resp;
}

/* auxiliary functions for ActiveResponseVertex */

LatentResponseVertex* ActiveResponseVertex::mkLatentResponseVertex(Graph & G) const {
	// pass true to the constructor of LatentResponseVertex, indicating memory
	LatentResponseVertex* v = new LatentResponseVertex(resp, ip, true);
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		Vertex* w = (*it);
		if ( w->getType() == Vertex::DETER_VIRUS ) {
			// the latent response keeps track of deterministic virus
			G.addEdge(w, v);
		}	else {
			throw std::logic_error("ill-formed interaction graph" + RIGHT_HERE);
		} // if/else
	} // for loop over in-edges
	// typically this function is called whenever the active response dies
	G.assignInOdeSys(v, false);
	// return the pointer
	return v;
}

/* auxiliary functions for LatentResponseVertex */

ActiveResponseVertex* LatentResponseVertex::mkActiveResponseVertex(Graph & G) const {
	ActiveResponseVertex* v = new ActiveResponseVertex(resp, ip);
	// make undirected edges from in-edges (from deter virus)
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		Vertex* w = (*it);
		if ( w->getType() == Vertex::DETER_VIRUS) {
			// the active response depends on deter virus, and vice-versa
			G.addUndirectedEdge(w, v);
		}	else {
			throw std::logic_error("ill-formed interaction graph" + RIGHT_HERE);
		} // if/else
	} // for
	G.assignInOdeSys(v, true); // the active response should be in the system od ODEs
	// return the pointer
	return v;
}

bool LatentResponseVertex::isMemory() const {
	return memory;
}
