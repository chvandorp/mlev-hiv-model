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

#ifndef INTEGRATOR_HPP_
#define INTEGRATOR_HPP_

#include <cmath>
#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include <numeric> // accumulate
#include <functional> // std::plus
#include <map>
#include <unordered_set>
#include <stdexcept>

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>

#include "macros.hpp"
#include "aux.hpp"
#include "hashed_definitions.hpp"


#define INVALID_INDEX -1
#define INVALID_ID -1

/* TODO/issues:
 * 1) override setThreshold, and take Rng as an argument (?)
 * 2) split the file into Vertex + Graph, and "user" defined
 * classes.
 * 3) make Vertex data private
 * 4) make data of specific vertices protected
 */

class Vertex; // used in VertexSet typedef
class Graph; // used by StochVirusVertex::mkDeterVirusVertex and DeterVirusVertex::mkStochVirusVertex
class StochVirusVertex; // used by DeterVirusVertex::mkStochVirusVertex
class DeterVirusVertex; // used by VirusVertex::getParents
class LatentResponseVertex; // used by ActiveResponseVertex::mkLatentResponseVertex
class ActiveResponseVertex; // do we need this?

typedef std::unordered_set<Vertex*> VertexSet;

class Vertex : public Printable {
public:
	typedef long int Id; // Vertex::Id
	enum VertexType { // FIXME: define these in vertices.hpp. How?
		NONE=0,
		QUIESCENT,
		TARGET,
		STOCH_VIRUS,
		DETER_VIRUS,
		LATENT_RESPONSE,
		ACTIVE_RESPONSE
	} type;
	virtual ~Vertex() {}; // sub-classes need to overwrite the destructor
	// add and remove edges... Only to modify this object...
	bool addEdgeIn(Vertex* v);
	bool addEdgeOut(Vertex* v);
	bool removeEdgeIn(Vertex* v);
	bool removeEdgeOut(Vertex* v);
	// remove all edgesIn/Out... Also modifies neighbors
	bool removeAllEdgesIn();
	bool removeAllEdgesOut();
	// integration-related functions
	virtual void deriv(const double y[], double dy[]) const=0;
	/* one vertex can represent multiple ODEs, which can be convenient
	 * when they are closely related (e.g. Target and Quiescent cells)
	 * dimension returns the number of ODEs this Vertex represents
	 */
	int dimension() const;
	/* assignIndex and getIndex are used when integrating the ODEs.
	 * They allow vetices to link neighbors to elements of the state vector.
	 * When the vertex is not inOdeSys, an exception is thrown.
	 * assignIndex assigns dim indices idx, idx+1, idx+2, ..., idx+dim-1
	 */
	int assignIndex(int ); // argument idx, returns dim
	int getIndex() const; // throws exception if not a scalar (dim != 1)
	int getIndex(int i) const;
	void assignValue(double ); // throws exception if not a scalar (dim != 1)
	void assignValue(double , int i);
	void assignDeriv(double dy[], double x) const; // throws exception if not a scalar (dim != 1)
	void assignDeriv(double dy[], double x, int i) const;
	double getValue() const; // throws exception if not a scalar (dim != 1)
	double getValue(const double y[]) const; // throws exception if not a scalar (dim != 1)
	double getValue(int i) const;
	double getValueSum() const; // the sum of the vector value
	double getValue(const double y[], int i) const; // if inOdeSys: return y[index], else value
	void assignThreshold(double );
	double getThreshold() const;
	VertexType getType() const;
	bool isInOdeSys() const;
	void assignInOdeSys(bool );
	void assignId(Id ); // used by Graph to set the Vertex::Id id
	Id getId() const;
	bool isInDiagram() const;
	// data
	VertexSet edgesIn;
	VertexSet edgesOut;
	// printing
	virtual void print(std::ostream & os) const;
protected:
	// special constructor to initialize a specific VertexType
	Vertex(VertexType t, int dim) : type(t), inOdeSys(false),
			dim(dim), index(dim, INVALID_INDEX),	value(dim, 0.0),
			threshold(1.0), id(INVALID_ID), inDiagram(false) {}
	/* bool inOdeSys: determines if the vertex should be included
	 * in the integration routine. For instance, non-viable mutants should
	 * not be included to the system od ODEs, to reduce the dimension.
	 * Keeping them in the graph, however may prove to be more efficient.
	 */
	bool inOdeSys;
  /* the int dim is the number of ODEs that this Vertex represents.
	 * it is the size of the vectors index and value
	 */
	int dim;
	/* int index: the index in the vector y used by GSL integrator.
	 * The Graph object automatically assigns an index to each added vertex.
	 * By default, the index = -1 (which would be invalid...)
	 */
	std::vector<int> index;
	/* The Host has a list of ImmuneResponses and a list of Viruses.
	 * These objects in turn have a pointer to a Vertex,
	 * such that they can set Vertex::type-specific values,
	 * and read the field Vertex::value that can be used as a population
	 * size, or as a cumulative rate.
	 */
	std::vector<double> value;
	/* Many variable have a (random) threshold... Something happens
	 * when the value reaches the threshold
	 */
	double threshold;
	/* Vertex id used for printing graphs in the dot language.
	 * the id is assigned by Graph, when a Vertex is added to Graph
	 */
	Id id;
	/* draw the Vertex in the dot-graph?
	 */
	bool inDiagram;
private:
	/* TODO: make index (and others) private, so that "user-defined"
	 * specializations don't have access to them
	 */
};

/* the Graph class */

class Graph {
public:
  Graph(); // default values for driver and dim_cache
  ~Graph(); // free driver
  double getTime() const; // returns the time t
	// graph-theory related functions
	int degree() const; // |V|
	int size() const; // |E|
	int dimension() const; // |e in E s.t. e.inOdeSys == true|
	// modify the graph
	/* addVertex(v) returns true if the v was not there yet
	 */
	bool addVertex(Vertex* v);
	/* addEdge adds the vertices if they are not in the graph yet.
	 * edges are directed!! v1 -> v2
	 */
	bool addEdge(Vertex* v1, Vertex* v2);
	/* add both v1 -> v2 and v2->v1 to Graph
	 */
	bool addUndirectedEdge(Vertex* v1, Vertex* v2);
	/* removeVertex(v) returns true if the vertex was in the graph
	 */
	bool removeVertex(Vertex* v);
	bool removeEdge(Vertex* v1, Vertex* v2);
	/* integration-related functions.
	 * advance the state dt units of time.
	 * If continuation is true, it is assumed that the previous state
	 * was left unchanged, and vertex indices remain valid.
	 * In this case no gsl_odeiv2_driver_reset (or re-allocation)
	 * is required.
	 */
	// add adges to the system of ODEs
	bool assignInOdeSys(Vertex* v, bool inOdeSys);
	// integrate the ODEs
  void integrate(double dt);
	// vertices
 	VertexSet vertices; // TODO: make protected
 	void printGraph(std::ostream & os) const; // for "dot (GraphViz)"
  // total_deriv needs access to vertices, dim_cahce and dim_capacity
  friend int total_deriv(double t, const double y[], double dy[], void* pars);
protected:
	// current time
	double t;
	/* give vertices an index used for integration,
	 * and resize the y_vec
	 * Returns a boolean indicating that the dimension has changed.
	 * The dimension is written to dim_cache, the dimension capacity
     * is written to dim_capacity.
	 */
	bool assignIndices();
	/* copy Vertex::value to y[Vertex::index] */
	void readValues();
	/* copy y[Vertex::index] to Vertex::value */
	void writeValues();
	/* y_vec is a vector used by the integrator. That is... the "data"
	 * double* y = y_vec.data() is used. std::vector takes care of all
	 * the memory allocation.
	 */
	std::vector<double> y_vec;
  // objects used by GSL integrator
  gsl_odeiv2_system sys;
  gsl_odeiv2_driver* driver;
  /* the previous dimension (i.e. after the last call of assignIndices).
   * Notice that the actual dimension can change when verrices are added or removed.
   */
  int dim_cache;
  /* the actual dimension passed to the driver object can be a bit larger
   * to allow growth of the system withour re-allocation
   */
  int dim_capacity;
  /* the previous integrator step size. This is passed to the allocator
   * of the GSL driver at re-allocation.
   */
  double h_cache;
  /* has the state been changed in between calls to integrate?
   * Methods that change the system om ODEs must set this boolean to false
   */
  bool continuation;
  /* every Vertex that is added to Graph is assigned a Vertex::Id,
   * that can be used for printing the Graph in the dot-language
   */
  Vertex::Id nextVertexId;
};

std::ostream & operator<<(std::ostream & , const Graph & ); // calls printGraph

#endif
