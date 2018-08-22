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

#ifndef RNGCLASS_HPP
#define RNGCLASS_HPP

#include <stdexcept>
#include <cmath>
#include <climits>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_zeta.h>

#include "macros.hpp"
#include "carnes.hpp" // realistic life spans

/* this is just a c++ class wrapper for GNU's gsl_rng and it's
 * gsl_randist, except for (at the moment) the methods Erlang, PosBinomial,
 * SkewNormal, and the Carnes lifespan distribution (and friends).
 * Also, some parameter conventions might be different
 * such as in BetaMnSd and Exponential.
 * TODO: replace with C++11 alternative? At least, make "compatible" with
 * C++11 <random> header.
 */

class Rng {
public:
	Rng();
	Rng(unsigned long );
	void seed(unsigned long ); // seed with 32 bit integer
	virtual ~Rng(); // destructor
// members to make Rng compatible with the <random> library (and e.g. std::shuffle)
	typedef unsigned long result_type;
  result_type operator()();
  result_type min() const;
  result_type max() const;
// distributions
	unsigned long Integer(); // call rand_int32
	unsigned long Integer(unsigned long ); // modulo max
	bool Bit(); // a random bit
	bool Bernoulli(double ); // Bernoulli(p) with values 0 = false, 1 = true
	double Uniform(); // uniform(0,1)
	double Uniform(double , double ); // Uniform(a,b)
	double Normal(double , double ); // mean and sd
	double LogNormal(double , double ); // location, scale parameterization
	double LogNormalMnSd(double , double ); // Mean, Sd parameterization
	double SkewNormal(double , double , double ); // mu, sigma and alpha
	double Exponential(double ); // pass the rate parameter (1/mean)
	double TruncExponential(double , double); // less than second argument, naive rejection method (todo)
	double Weibull(double , double ); // scale, shape
	double Erlang(double , unsigned ); // simple Erlang(rate,shape) (sum of exponentials)
	double Gamma(double , double ); // Gamma(scale,shape) distribution
	double ParetoMinusOne(double ); // Pareto(a, 1)-1
	double Beta(double , double ); // Beta(alpha,beta) distribution
	double BetaMnSd(double , double ); // Beta, but with mean and sd as parameterization
	double Carnes(double a0=0.0); // modern North-American life spans, conditioned on a minimal age (0.0 by default)
	double CarnesAge(); // sample ages given a population in equilibrium and the Carnes life span distribution
	int Binomial(int , double ); // Binom(n,p)
	int PosBinomial(int , double); // Binomial, conditioned on > 0
	int Poisson(double ); // Poisson(lambda) with mean lambda
 	int PosPoisson(double ); // a.k.a. Zero-truncated Poisson distribution
	int Zeta(double ); // the Zipf or zeta distribution (power law)
	int Hypergeometric(int , int , int ); // Hyper(n, K, N) (N = population size, K = successes, n = draws)
  int Skellam(double , double ); // Skellam(lambda_1, lambda_2) = Poisson(lambda_1) - Poisson(lambda_2)
  int SymSkellam(int, double ); // Symmetric Skellam distribution with location parameter and sd
	void Shuffle(int* , int , int ); // a random sample (unique)
  void Dirichlet(double* , int , const double* ); // the Dirichlet distribution
private:
	gsl_rng* gslRng;
	void init();
	// make copy constructor and assignment operator unavailable, they don't make sense
	Rng(const Rng & ); // copy constructor not defined
	void operator=(const Rng & ); // assignment operator not defined
};

// auxilliary functions

double expectedCarnesSurvival();


#endif /* RNGCLASS_HPP_ */
