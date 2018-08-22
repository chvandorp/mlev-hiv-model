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

#include "integrator.hpp"

// a map to associated Vertices with colors
std::map<Vertex::VertexType, std::string> vertexColorMap = {
	std::make_pair(Vertex::NONE,            "white"),
	std::make_pair(Vertex::TARGET,          "black"),
	std::make_pair(Vertex::DETER_VIRUS,     "red"),
	std::make_pair(Vertex::STOCH_VIRUS,     "violet"),
	std::make_pair(Vertex::ACTIVE_RESPONSE, "cyan"),
	std::make_pair(Vertex::LATENT_RESPONSE, "yellow"),
};

// a map to associated Vertices with names
std::map<Vertex::VertexType, std::string> vertexNameMap = {
	std::make_pair(Vertex::NONE,            "none"),
	std::make_pair(Vertex::TARGET,          "target"),
	std::make_pair(Vertex::DETER_VIRUS,     "deter_virus"),
	std::make_pair(Vertex::STOCH_VIRUS,     "stoch_virus"),
	std::make_pair(Vertex::ACTIVE_RESPONSE, "active_response"),
	std::make_pair(Vertex::LATENT_RESPONSE, "latent_response"),
};

/** Methods for Node and Graph
 * For the classes that inherit Node,
 * a component of the vector field needs to be defined
 * by overwriting the virtual function dy.
 * This determines the system of ODEs.
 */

/* Functions used to modify the graph */

bool Vertex::addEdgeIn(Vertex* v) {
	auto pair = edgesIn.insert(v);
    return pair.second;
}
bool Vertex::addEdgeOut(Vertex* v) {
	auto pair = edgesOut.insert(v);
    return pair.second;
}
bool Vertex::removeEdgeIn(Vertex* v) {
	auto it = std::find(edgesIn.begin(), edgesIn.end(), v);
	if ( it != edgesIn.end() ) {
		edgesIn.erase(it);
		return true;
	} // if
	else return false;
}
bool Vertex::removeEdgeOut(Vertex* v) {
	auto it = std::find(edgesOut.begin(), edgesOut.end(), v);
	if ( it != edgesOut.end() ) {
		edgesOut.erase(it);
		return true;
	} // if
	else return false;
}

bool Vertex::removeAllEdgesIn() {
	bool ok = true;
	for ( auto it = edgesIn.begin(); it != edgesIn.end(); ++it ) {
		Vertex* v = (*it);
		if ( v == nullptr ) {
			throw std::logic_error("NULL pointers in Vertex neighbors" + RIGHT_HERE);
		}
		ok &= v->removeEdgeOut(this); // direction: v -> this
	}
	edgesIn.clear();
	return ok;
}

bool Vertex::removeAllEdgesOut() {
	bool ok = true;
	for ( auto it = edgesOut.begin(); it != edgesOut.end(); ++it ) {
		Vertex* v = (*it);
		if ( v == nullptr ) {
			throw std::logic_error("NULL pointers in Vertex neighbors" + RIGHT_HERE);
		}
		ok &= v->removeEdgeIn(this); // direction: this -> v
	}
	edgesOut.clear();
	return ok;
}

int Vertex::dimension() const {
	return dim;
}

int Vertex::assignIndex(int idx) {
	if ( !inOdeSys ) {
		throw std::logic_error("assigning index to vertex that is "
				"not part of the system of ODEs" + RIGHT_HERE);
	}
	for ( int i = 0; i < dim; ++i ) {
		index[i] = idx + i;
	}
	return dim; // can be used for the next index to be used by the integrator
}

int Vertex::getIndex() const {
	if ( dim != 1 ) {
		throw std::out_of_range("Vertex is not a scalar" + RIGHT_HERE);
	}
	return getIndex(0);
}

int Vertex::getIndex(int i) const {
	if ( !inOdeSys ) {
		throw std::logic_error("requesting index of vertex that is "
				"not part of the system of ODEs" + RIGHT_HERE);
	}
	if ( i < 0 || i >= dim ) {
		throw std::out_of_range("Vertex has different size than caller expected" + RIGHT_HERE);
	}
	return index[i];
}

void Vertex::assignValue(double x) {
	if ( dim != 1 ) {
		throw std::out_of_range("Vertex is not a scalar" + RIGHT_HERE);
	}
	assignValue(x, 0);
}

void Vertex::assignValue(double x, int i) {
	if ( i < 0 || i >= dim ) {
		throw std::out_of_range("Vertex has different size than caller expected" + RIGHT_HERE);
	}
	value[i] = x;
}

void Vertex::assignDeriv(double dy[], double x) const {
	if ( dim != 1 ) {
		throw std::out_of_range("Vertex is not a scalar" + RIGHT_HERE);
	}
	assignDeriv(dy, x, 0);
}

void Vertex::assignDeriv(double dy[], double x, int i) const {
	if ( i < 0 || i >= dim ) {
		throw std::out_of_range("Vertex has different size than caller expected" + RIGHT_HERE);
	}
	dy[index[i]] = x;
}

double Vertex::getValue() const {
	if ( dim != 1 ) {
		throw std::out_of_range("Vertex is not a scalar" + RIGHT_HERE);
	}
	return getValue(0);
}

double Vertex::getValue(int i) const {
	if ( i < 0 || i >= dim ) {
		throw std::out_of_range("Vertex has different size than caller expected" + RIGHT_HERE);
	}
	return value[i];
}

double Vertex::getValueSum() const { // the sum of the elements in the vector value
	return std::accumulate(value.begin(), value.end(), 0.0, std::plus<double>());
}

double Vertex::getValue(const double y[]) const {
	if ( dim != 1 ) {
		throw std::out_of_range("Vertex is not a scalar" + RIGHT_HERE);
	}
	return getValue(y, 0);
}

double Vertex::getValue(const double y[], int i) const {
	// if inOdeSys: return y[index], else value
	if ( i < 0 || i >= dim ) {
		throw std::out_of_range("Vertex has different size than caller expected" + RIGHT_HERE);
	}
	if ( inOdeSys ) {
		return y[index[i]];
	}	else {
		return value[i];
	}
}

void Vertex::assignThreshold(double threshold) { this->threshold = threshold; }
double Vertex::getThreshold() const { return threshold; }
Vertex::VertexType Vertex::getType() const { return type; }
bool Vertex::isInOdeSys() const { return inOdeSys; }
void Vertex::assignInOdeSys(bool inOdeSys) { this->inOdeSys = inOdeSys; }
void Vertex::assignId(Id id) { this->id = id; }
Vertex::Id Vertex::getId() const { return id; }
bool Vertex::isInDiagram() const { return inDiagram; }

void Vertex::print(std::ostream & os) const {
	os << "<vertex "
	   << "type='" << vertexNameMap[type] << "' "
		 << "dim='" << dim << "' "
	   << "inodesys='" << std::boolalpha << inOdeSys << std::noboolalpha << "' "
	   << "threshold='" << threshold << "' "
	   << ">" << std::endl;
	for ( int i = 0; i < dim; ++i ) {
		os << "<value "
		   << "i='" << i << "' "
			 << "x='" << value[i] << "' "
			 << "/>" << std::endl;
	}
	os << "</vertex>";
	// TODO: add more data: in/out degree etc.
}

/************************ Graph methods *************************/


Graph::Graph() : t(0.0), driver(nullptr), dim_cache(0), dim_capacity(0),
			h_cache(INIT_STEPSIZE_INTEGRATOR), continuation(false), nextVertexId(0) {
  // do not keep elements un-initialized...
	sys = {nullptr, nullptr, size_t(dim_capacity), this};
}

Graph::~Graph() {
  if ( driver != nullptr ) { // otherwise free gives an error
    gsl_odeiv2_driver_free(driver);
  }
}

double Graph::getTime() const { return t; }

int Graph::degree() const { return vertices.size(); } // the number of verices

int Graph::size() const { // the number of edges
	int e = 0;
	for ( auto vit = vertices.begin(); vit != vertices.end(); ++vit ) {
		e += (*vit)->edgesOut.size();
	}
	return e;
}

int Graph::dimension() const {
    // Note that dim_cache may not always be valid
	auto op = [](int dim, Vertex* v){ return (v->isInOdeSys() ? v->dimension() : 0); };
	return std::accumulate(vertices.begin(), vertices.end(), 0, op);
}

bool Graph::addVertex(Vertex* v) {
	continuation = false; // STATE MODIFIED (TODO: be more relaxed)
	if ( v == nullptr ) {
		throw std::invalid_argument("NULL pointers in Graph vertices" + RIGHT_HERE);
	}
	auto pair = vertices.insert(v);
	if ( pair.second ) { // the vertex wasn't in the Graph yet!
		v->assignId(nextVertexId++); // pass nextVertexId, then increment nextVertexId
		return true;
	}	else {
		return false;
	}
}

bool Graph::addEdge(Vertex* v1, Vertex* v2) { // directed!! v1 -> v2
	continuation = false; // STATE MODIFIED (TODO: be more relaxed)
	addVertex(v1); // TODO do something with the returned bool: error checking
  addVertex(v2);
	bool ok1 = v1->addEdgeOut(v2);
	bool ok2 = v2->addEdgeIn(v1);
	if ( ok1 != ok2 ) {
		throw std::logic_error("found unbalanced edges in graph" + RIGHT_HERE);
	}
	return ok1;
}

bool Graph::addUndirectedEdge(Vertex* v1, Vertex* v2) {
	continuation = false; // STATE MODIFIED (TODO:be more relaxed)
	// undirected!! v1 -> v2 and v2 -> v1
	bool ok1 = addEdge(v1, v2);
	bool ok2 = addEdge(v2, v1);
	return (ok1 && ok2);
}

bool Graph::removeVertex(Vertex* v) { // removes all edges too
	continuation = false; // STATE MODIFIED (TODO:be more relaxed)
	auto it = vertices.find(v);
	if ( it != vertices.end() ) {
		// remove v from the neighbors for each of v's neighbors
		if ( v == nullptr ) {
			throw std::logic_error("NULL pointers in Graph vertices" + RIGHT_HERE);
		}
		bool okIn = v->removeAllEdgesIn();
		bool okOut = v->removeAllEdgesOut();
		if ( !okIn || !okOut ) {
			throw std::logic_error("found unbalanced edges in graph" + RIGHT_HERE);
		}
		// remove the vertex from the list of vertices
		vertices.erase(it);
		return true;
	}
	else return false;
}

bool Graph::removeEdge(Vertex* v1, Vertex* v2) { // directed!! v1 -> v2
	continuation = false; // STATE MODIFIED (TODO:be more relaxed)
	auto it1 = vertices.find(v1);
	auto it2 = vertices.find(v2);
	if ( it1 != vertices.end() && it2 != vertices.end() ) {
		if ( v1 == nullptr || v2 == nullptr ) {
			throw std::logic_error("NULL pointers in Graph vertices" + RIGHT_HERE);
		}
		bool ok1 = v1->removeEdgeOut(v2);
		bool ok2 = v2->removeEdgeIn(v1);
		if ( ok1 != ok2 ) {
			throw std::logic_error("found unbalanced edges in graph" + RIGHT_HERE);
		}
		return ok1;
	}
	else return false;
}

bool Graph::assignInOdeSys(Vertex* v, bool inOdeSys) {
	continuation = false; // STATE MODIFIED (TODO:be more relaxed)
	/* add or remove a vertex from the system of ODEs. returns false if
	 * the vertex was not in the Graph, or if the state was unchanged.
	 * TODO: move the vertex from constant vertices to ODE system vertices
	 */
	auto it = vertices.find(v);
	if ( it != vertices.end() ) {
		if ( v == nullptr ) {
			throw std::logic_error("NULL pointers in Graph vertices" + RIGHT_HERE);
		}
		if ( v->isInOdeSys() != inOdeSys ) {
			v->assignInOdeSys(inOdeSys);
			return true;
		} else {
			return false;
		}
	}	else {
		return false; // the vertex was not in the graph (Graph remains unchanged)
	}
}


/* A function used by GSL to calculate the total derivative */


int total_deriv(double t, const double y[], double dy[], void* pars) { // is a friend of Graph
  int status = GSL_SUCCESS;

  // cast the void to Graph
  if ( pars == nullptr ) {
		throw std::invalid_argument("passed NULL pointer as 'pars' argument" + RIGHT_HERE);
	}
  Graph* G = (Graph*) pars;

  // get vectors from y
	for ( auto it = G->vertices.begin(); it != G->vertices.end(); ++it ) {
		/* test if a vertex is included in the ODE system,
		 * then call the vertices' deriv method
		 * TODO: it may be more efficient to put the ODE vertices in a
		 * separate list (that is, the pointers)
		 */
		if ( (*it)->isInOdeSys() ) {
			(*it)->deriv(y, dy);
		}
	}
  /* make sure that the additional capacity does not change the adaptive step size
   * TODO: find out how GSL initializes dy
   */
  for ( int i = G->dim_cache; i < G->dim_capacity; ++i ) {
    dy[i] = 0.0;
  }
  return status;
}

void Graph::integrate(double dt) {
	/** advance the state dt units of time.
	 * the boolean continuation is set true after a succesful call of integrate
	 * but will be made false when modifying the Graph.
	 */
	if ( !continuation ) { // reset of re-allocate the ODE solver
		// assign indices to vertices and values to y_vec
		bool dim_cap_changed = assignIndices(); // give every Vertex a y-index and resize y_vec
		// check that the dimension is positive
		if ( dim_cache == 0 ) return;
		// else: positive dimension, continue resetting/(re-)allocating
		if ( dim_cap_changed || driver == nullptr ) { // the driver needs to be re-allocated
			if ( driver != nullptr ) { // otherwise free gives an error
				gsl_odeiv2_driver_free(driver);
				driver = nullptr; // redundant
			}
			sys = {total_deriv, nullptr, size_t(dim_capacity), this}; // is just a struct
			// the reals are h_start, eps_abs, eps_rel, a_y, a_dydt. TODO: parameters
			driver = gsl_odeiv2_driver_alloc_standard_new(&sys, gsl_odeiv2_step_rkf45, h_cache, 1e-6, 1e-6, 1.0, 1.0); // TODO: relax?
			gsl_odeiv2_driver_set_hmax(driver, MAX_STEPSIZE_INTEGRATOR);
		}	else { // we only have to reset the driver, because the dimension capacity is unchanged
			gsl_odeiv2_driver_reset(driver);
		}
		// the state has somehow changed, and hence we need to fetch the values form the verites
		readValues(); // copy Vertex::value to y[Vertex::index]
	}
	// now assume that the driver has been allocated properly

	int status = GSL_SUCCESS; // or is it...

	double* y = y_vec.data(); // returns a pointer to the array used by std::vector
	double t1 = t + dt;

	status = gsl_odeiv2_driver_apply(driver, &t, t1, y); // y and t get updated...
	if ( status != GSL_SUCCESS ) {
		throw std::logic_error("GSL returned error code" + RIGHT_HERE);
	}
	// save the step size
	h_cache = driver->h;
	writeValues(); // copy y[Vertex::index] to vertex::value
	continuation = true; // calling integrate again skips resetting/reallocating
}

bool Graph::assignIndices() {
	int index = 0;
	for ( auto it = vertices.begin(); it != vertices.end(); ++it) {
		Vertex* v = (*it);
		if ( v == nullptr ) {
			throw std::logic_error("NULL pointers in Graph vertices" + RIGHT_HERE);
		}
		else if ( v->isInOdeSys() ) {
			// assignIndex returns the dimension of the vertex
			index += v->assignIndex(index);
		}
	}
	// index now equals the "effective degree" (same as Graph::dimension())
  if ( dim_cache != index ) {
    dim_cache = index; // update the dim_cache
    // we may need to resize and report this (re-allocate integrator objects)
    if ( dim_cache > dim_capacity ) { // increase dim_capacity if it is too small
      dim_capacity = (4*dim_cache)/3; // TODO: what is a good rule? TODO: parameters
      y_vec.resize(dim_capacity, 0.0); // std::vector handles all the memory issues...
      return true;
    } else if ( 3*dim_cache < 2*dim_capacity ) { // decrease dim_capacity id if is too large
      dim_capacity = (4*dim_cache)/3; // TODO: what is a good rule? TODO: parameters
      y_vec.resize(dim_capacity, 0.0); // std::vector handles all the memory issues...
      return true;
    } else {
			return false;
		}
  } else {
    // no change: keep y_vec as it was, no need to re-allocate the driver object
    return false;
  }
}

void Graph::readValues() {
	for ( auto it = vertices.begin(); it != vertices.end(); ++it ) {
		Vertex* v = (*it);
		if ( v == nullptr ) {
			throw std::logic_error("NULL pointers in Graph vertices" + RIGHT_HERE);
		}	else if ( v->isInOdeSys() ) {
			for ( int i = 0; i < v->dimension(); ++i ) {
				// TODO: pass y_vec to a method of v to do the assigning?
				y_vec[v->getIndex(i)] = v->getValue(i);
			}
		}
	}
}

void Graph::writeValues() {
	for ( auto it = vertices.begin(); it != vertices.end(); ++it ) {
		Vertex* v = (*it);
		if ( v == nullptr ) {
			throw std::logic_error("NULL pointers in Graph vertices" + RIGHT_HERE);
		}	else if ( v->isInOdeSys() ) {
			for ( int i = 0; i < v->dimension(); ++i ) {
				// TODO: pass y_vec to a method of v to do the assigning?
				v->assignValue(y_vec[v->getIndex(i)], i);
			}
		}
	}
  // to be safe, make the extra capacity zero
  for ( int i = dim_cache; i < dim_capacity; ++i ) {
    y_vec[i] = 0.0;
  }
}

void Graph::printGraph(std::ostream & os) const {
	/* TODO: 1) use subgraph to represent the bipartite nature
	 * of the interactions (i.e. responses and viruses).
	 * 2) add the value as an attribute?
	 */
	os << "digraph {\n" // start graph
	   << "rankdir=LR;\n"; // landscape orientation
	for ( auto it = vertices.begin(); it != vertices.end(); ++it ) {
		Vertex* v = (*it);
		if ( v->isInDiagram() ) {
			Vertex::Id iv = v->getId();
			os << iv << " [style=filled color=" << vertexColorMap[v->getType()] << "];\n";
		}
	}
	for ( auto it = vertices.begin(); it != vertices.end(); ++it ) {
		Vertex* v = (*it);
		if ( v->isInDiagram() && !v->edgesOut.empty() ) {
			// TODO: there might not be any neigbors inOdeSys
			Vertex::Id iv = v->getId();
			os << iv << " -> { ";
			for ( auto jt = v->edgesOut.begin(); jt != v->edgesOut.end(); ++jt ) {
				Vertex* w = (*jt);
				if ( w->isInDiagram() ) {
					Vertex::Id iw = w->getId();
					os << iw << " ";
				}
			}
			os << "};\n";
		} // else, don't print adjacency
	}
	os << "}";
}

std::ostream & operator<<(std::ostream & os, const Graph & G) {
	G.printGraph(os);
	return os;
}
