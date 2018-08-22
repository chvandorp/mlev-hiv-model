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

#ifndef AUX_HPP_
#define AUX_HPP_

#include <cmath> // log, exp
#include <string>
#include <stdexcept>

#include "macros.hpp"
#include "rng.hpp"
#include "hashed_definitions.hpp"

/* typedefs */

typedef long int Id; // used for identification of agents (Agent, Host, Virus, etc.)
#define ANON_ID -1

/* interfaces */

class Printable {
public:
	Printable() { /* empty */ }
	virtual ~Printable() { /* empty */ }
	virtual void print(std::ostream & ) const = 0;
};

std::ostream & operator<<(std::ostream & , const Printable & );

/* a class that allows for random shuffling by sorting a list
 * For instant std::list<Agent*> or std::list<Contact*>
 */

class RandomIndexable {
public:
	RandomIndexable() : random_index(0) {}
	virtual ~RandomIndexable() {}
	int getRandomIndex() const;
	// use an RNG to reset the random_index
	void setRandomIndex(Rng & );
private:
 	unsigned long random_index;
};

bool compareByRandomIndex(RandomIndexable* , RandomIndexable* );


/* inline functions */

inline double logit(double x) { return log(x/(1-x)); }
inline double logit_inv(double x) {
	double y = exp(x);
	return y/(1+y);
}
inline double extinction_threshold(double u, double x0, double R0) {
    /* calculate the extinction threshold for a declining population (R0 < 1).
     * The threshold for switching between ODEs and a stochastic description
     * is at x0, and u must be a Uniform(0,1) deviate
     */
    double squ = pow(u, 1.0/x0);
    return x0*(1-squ)/(1-squ*R0);
}
inline double major_outbreak_prob(double R) {
    /* The probability of a "major outbreak" given a reproduction number.
     * this is for adding mutants and activating effectors.
     * Notice that this function has a cusp at R = 1
     */
    return ( R > 1.0 ? 1-1/R : 0.0 );
}
inline double deriv_major_outbreak_prob(double R) {
    /* the derivative of major_ourbreak_prob
     * Notice that 1-1/R is not differentiable at R = 1
     */
    return ( R > 1.0 ? 1/(R*R) : 0.0 );
}
inline double major_outbreak_prob_C1(double R, double epsilon=0.1) {
    /* a C^1 approximation of major_outbreak_prob given by
     * (R-1)^2/(R*(R-1) + epsilon)
     * Notice that when epsilon = 0, we get (R-1)^2/(R*(R-1)) = 1 - 1/R
     */
    return ( R > 1.0 ? (R-1)*(R-1) / (R*(R-1) + epsilon) : 0.0 );
}
inline double deriv_major_outbreak_prob_C1(double R, double epsilon=0.1) {
    // the derivative (w.r.t. R) of major_outbreak_prob_C1
    double nmr = R*(R-1) + epsilon;
    return (R-1)*(2*epsilon+R-1)/(nmr*nmr);
}
inline int mod(int m, int n) {
    // compute m mod n with a representative 0 <= r <= n
    if ( m < 0 ) m += ((-m)/n + 1)*n;
    return m % n;
}

std::string ageToColor(double age);

double quick_power(double x, int n);

/* return a weighted geometric average of a and b,
 * where w is the weight of a
 */
double weighted_geometric_average(double a, double b, double w);


// alias for std::pair<double, double>
typedef std::pair<double, double> Vec;

// a 2x2 matrix
struct Mat {
	Mat() : a(0), b(0), c(0), d(0) {}
	Mat(double a, double b, double c, double d) : a(a), b(b), c(c), d(d) {}
	double a, b, c, d;
};

/* compute the eigen values of a 2x2 matrix. (a, b; c, d)
 * second argument of the pair is true if the eigenvalues are real,
 * and false if the eigenvalues are complex.
 * if the eigenvalues are real, the first element of the pair
 * consists of the two eigenvalues.
 * if the eigenvalues are complex, the real and imaginary part are returned
 */
std::pair<Vec, bool> eigen_values(double a, double b, double c, double d);
std::pair<Vec, bool> eigen_values(const Mat & M);

/* compute the eigen vector of (a, b; c, d) corresponding to eigenvalue lambda
 * this function assumes that the eigenvalue is correct
 * eigen_vector_norm1 is a vector of Euclidean norm.
 * eigen_vector_sum1 is a vector that sums to 1
 */
Vec eigen_vector_norm1(double a, double b, double c, double d, double lambda);
Vec eigen_vector_norm1(const Mat & M, double lambda);

Vec	eigen_vector_sum1(double a, double b, double c, double d, double lambda);
Vec	eigen_vector_sum1(const Mat & M, double lambda);


#endif
