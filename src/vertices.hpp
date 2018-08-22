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

#include "integrator.hpp" // defines the Vertex interface
#include "ilevel_equations.hpp"
#include "virus.hpp"
#include "immune_response.hpp"

/** Here we define classes that inherit the Vertex interface and
 * define the system of ODEs.
 * We have the following classes (with dimension):
 *
 * TargetVertex (2)
 * VirusVertex --- DeterVirusVertex (2)
 *              |- StochVirusVertex (1)
 * ResponseVertex --- ActiveResponseVertex (1)
 *                 |- LatentResponseVertex (1)
 *
 * TODO: make parent class for activated responses and deterministic
 * viruses that have the member endangered, and a method isEndangered()
 * to be used in CiScheme
 */

/* VirusVertex and its derived classes
 * DeterVirusVertex and StochVirusVertex
 */

class VirusVertex : public Vertex {
public:
  virtual ~VirusVertex() {}
	/* deriv remains undefined, shared auxiliary functions for Deter-
	 * and StochVirusVertex:
	 */
	double getR(bool malthusian=false, VirusVertex* mimic_resps=nullptr) const;
	std::pair<int, double> getNumberOfResponses() const;
	// list of responses that recognize an epitpe
	std::list<ActiveResponseVertex*> getResponses() const;
	// any deterministic virus at hamming distance 1
  int getNumberOfParents() const;
	// list of deterministic clones of which this is a mutant (edgesIn)
  std::list<DeterVirusVertex*> getParents() const;
	// list of stochastic clones that are mutants of this (edgesOut)
	std::list<StochVirusVertex*> getChildren() const;
	// list of VirusVertex in edgesOut
	std::list<VirusVertex*> getAllChildren() const;
	Virus* getVirus() const; // return pointer to the corresponding Virus object
  bool isEscape() const; // are there parents with more responses that this guy
	// data (TODO make private, NB: needed by other verices)
	double pr_exact, pr_mut; // probability of exact replication, mutantion
  double beta, d_I1, d_I2, m, d_V, p_V, gamma, f, omega; // parameters
	double c; // compound parameters
  Virus* vir;
  // access to host-specific parameters
  const IndPars* ip; // TODO: obsolete?
protected:
	// defined in cpp, calls Vertex(t)
  VirusVertex(VertexType t, int size, Virus* vir, const IndPars* ip);
};

class DeterVirusVertex : public VirusVertex {
public:
	DeterVirusVertex(Virus* vir, const IndPars* ip); // calls VirusVertex(DETER_VIRUS, 2, vir, ip)
	virtual void deriv(const double y[], double dy[]) const;
	// aux functions
	StochVirusVertex* mkStochVirusVertex(Graph & G) const; // modifies G. TODO: remove this from G?
	double assignInitialCondition(double inoculum);
	double getQssVL(double T) const; // give the number of target cells
	double getQssVL() const; // searches for the TargetVertex amoung its in-edges
	double getQssVL(const double y[]) const; // can be used by the integrator (the others can NOT)
	// data (TODO: make private)
  bool endangered; // becomes true when the value drops below a fixed threshold
};

class StochVirusVertex : public VirusVertex {
public:
	StochVirusVertex(Virus* vir, const IndPars* ip); // calls VirusVertex(STOCH_VIRUS, 1, vir, ip)
	virtual void deriv(const double y[], double dy[]) const;
	// aux functions
	DeterVirusVertex* mkDeterVirusVertex(Graph & G) const; // modifies G. TODO: remove this from G?
};

/* ResponseVertex and its derived classes
 * ActiveResponseVertex and LatentResponseVertex
 */

class ResponseVertex : public Vertex {
public:
  virtual ~ResponseVertex() {}
 	double getR(bool malthusian=false) const; // the reproduction number
	ImmuneResponse* getResponse() const; // return a pointer to the corresponding ImmuneResponse object
  // data (TODO make private)
  double h, k;
  ImmuneResponse* resp;
  // access to host-specific parameters
  const IndPars* ip; // TODO: obsolete?
protected:
  ResponseVertex(VertexType t, ImmuneResponse* resp, const IndPars* ip); // calls Vertex(t, 1)
	// data
	double p_E, d_E; // proliferation and death rate parameters of the response
};

class ActiveResponseVertex : public ResponseVertex {
public:
	ActiveResponseVertex(ImmuneResponse* resp, const IndPars* ip); // calls ResponseVertex(ACTIVE_RESPONSE, resp, ip);
	virtual void deriv(const double y[], double dy[]) const;
	// aux functions
	LatentResponseVertex* mkLatentResponseVertex(Graph & G) const; // modifies G. TODO: remove this from G?
	// data (TODO make private)
  bool endangered; // becomes true when the value drops below a fixed threshold
};

class LatentResponseVertex : public ResponseVertex {
public:
	// calls ResponseVertex(LATENT_RESPONSE, resp, ip)
	LatentResponseVertex(ImmuneResponse* resp, const IndPars* ip, bool memo);
	virtual void deriv(const double y[], double dy[]) const;
	// aux functions
	ActiveResponseVertex* mkActiveResponseVertex(Graph & G) const; // modifies G. TODO: remove this from G?
	bool isMemory() const; // returns memory
protected:
	// data
	double lambda; // activation rate
	bool memory; // if true, then the repsonse has been active at some point in time
};

/* TargetVertex */

class TargetVertex : public Vertex {
public:
	TargetVertex(const IndPars* ip); // calls Vertex(TARGET, 2)
	virtual void deriv(const double y[], double dy[]) const;
  // access to host-specific parameters
  const IndPars* ip; // TODO: obsolete?
protected:
	// data
	double d_Q, d_T, s, h_Q, omega0;
};
