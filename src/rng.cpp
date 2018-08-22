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

#include "rng.hpp"



/******* members for Rng ********/

Rng::Rng() {
	init(); // uses standard seed
}
Rng::Rng(unsigned long s) {
	init();	seed(s);
}
Rng::~Rng() {
	gsl_rng_free(gslRng);
}
void Rng::init() {
	gslRng = gsl_rng_alloc(gsl_rng_taus);
}
void Rng::seed(unsigned long s) {
	gsl_rng_set(gslRng, s);
}

// compatibility with <random>
Rng::result_type Rng::operator()() {
    return Integer();
}
Rng::result_type Rng::min() const {
	return 0;
}
Rng::result_type Rng::max() const {
	return gsl_rng_max(gslRng) - gsl_rng_min(gslRng);
}
// distributions

unsigned long Rng::Integer() {
	return gsl_rng_get(gslRng) - gsl_rng_min(gslRng);
}
unsigned long Rng::Integer(unsigned long max) { // modulo max
	return gsl_rng_uniform_int(gslRng, max);
}
bool Rng::Bit() {
	return gsl_rng_uniform_int(gslRng, 2) == 1;
}
bool Rng::Bernoulli(double p) {
	if ( p < 0.0 || p > 1.0 ) {
		throw std::invalid_argument("p must be between 0 and 1" + RIGHT_HERE);
	}
	return gsl_ran_bernoulli(gslRng, p);
}

double Rng::Uniform() { // uniform(0,1), excludes 0 and 1
	 return gsl_ran_flat(gslRng, 0.0, 1.0);
}

double Rng::Uniform(double a, double b) { // Uniform(a,b)
	 if ( a <= b ) {
		 if ( a==b ) {
			 return a;
		 } else {
			 return gsl_ran_flat(gslRng, a, b);
		 }
	 } else {
		 throw std::invalid_argument("a must be smaller than b" + RIGHT_HERE);
	 }
}

double Rng::Normal(double mu, double sigma) { // mean and sd
	if ( sigma >= 0.0 ) {
		if ( sigma == 0.0 ) {
			return mu;
		}	else {
			return mu + gsl_ran_gaussian(gslRng, sigma);
		}
	}	else {
		throw std::invalid_argument("sigma must be non-negative" + RIGHT_HERE);
	}
}

double Rng::LogNormal(double m, double s) { // location, scale parameterization
	if ( s < 0.0 ) {
		throw std::invalid_argument("s must be non-negative" + RIGHT_HERE);
	}
	if ( s == 0.0 ) return exp(m);
	return exp(Normal(m, s));
}

double Rng::LogNormalMnSd(double mu, double sigma) { // Mean, Variance parameterization
	if ( mu <= 0.0 || sigma < 0.0 ) {
		throw std::invalid_argument("non-positive mu or negative sigma" + RIGHT_HERE);
	}
	if ( sigma == 0.0 ) return mu;
	// now handle non-trivial cases...
	double x = 1.0 + (sigma*sigma) / (mu*mu);
	double s = sqrt(log(x));
	double m = log(mu/sqrt(x));
	return LogNormal(m, s);
}

double Rng::SkewNormal(double m, double s, double alpha) {
	if ( s < 0.0 ) {
		throw std::invalid_argument("s must be nonnegative" + RIGHT_HERE);
	}
	if ( s == 0.0 ) return m; // in the limit s --> 0
	// now handle non-trivial cases
	double Z1 = Normal(0.0,1.0);
	if ( alpha == 0.0 ) return s * Z1 + m; // not skewed, just normal
	else {
		// X = 1/sqrt{1+alpha^2} alpha |Z1| + Z2
		double Z2 = Normal(0.0,1.0);
		double absZ1 = (Z1 > 0 ? Z1 : -Z1);
		return m + s * ((alpha*absZ1 + Z2) / sqrt(1.0 + alpha*alpha));
	}
}

double Rng::Exponential(double lambda) {
	if ( lambda > 0 ) {
		return gsl_ran_exponential(gslRng, 1.0/lambda);
	} else {
		throw std::invalid_argument("lambda must be positive" + RIGHT_HERE);
	}
}
double Rng::TruncExponential(double lambda, double max) {
	double x;
	do {
		x = Exponential(lambda);
	} while ( x > max );
	return x;
}

double Rng::Weibull(double a, double b) { // scale, shape
	if ( a >= 0.0 && b >= 0.0 ) {
		if ( a == 0 ) {
			return 0.0;
		}	else {
			return gsl_ran_weibull(gslRng, a, b);
		}
	}	else {
		throw std::invalid_argument("a and b must be non-negative" + RIGHT_HERE);
	}
}

double Rng::Erlang(double lambda, unsigned k) { // simple Erlang(rate,shape) (sum of exponentials)
	if ( lambda > 0 ) {
		double ans = 0.0;
		for ( unsigned i = 0; i < k; ++i ) {
			ans += Exponential(lambda);
		}
		return ans;
	} else {
		throw std::invalid_argument("lambda must be positive" + RIGHT_HERE);
	}
}

double Rng::Gamma(double scale, double shape) { // Gamma(scale,shape) distribution {
	if ( scale < 0.0 || shape <= 0.0 ) {
		throw std::invalid_argument("negative scale or nonnegative shape" + RIGHT_HERE);
	}
	if ( scale == 0.0 ) {
		return 0.0;
	}	else {
		return gsl_ran_gamma(gslRng, shape, scale);
	}
}

double Rng::ParetoMinusOne(double a) {
	if ( a <= 0 ) {
		throw std::invalid_argument("a must be positive" + RIGHT_HERE);
	}
	return gsl_ran_pareto(gslRng, a, 1.0) - 1.0;
}

double Rng::Beta(double alpha, double beta) { // Beta(alpha,beta) distribution
	if ( alpha > 0.0 && beta > 0.0 ) {
		return gsl_ran_beta(gslRng, alpha, beta);
	}	else {
		throw std::invalid_argument("alpha and beta must be positive" + RIGHT_HERE);
	}
}

double Rng::BetaMnSd(double M, double SD) {
	if ( M <= 0.0 || M >= 1.0 || SD <= 0.0 ) {
		throw std::invalid_argument("M must be between 0 and 1 and SD must be positive" + RIGHT_HERE);
	}
	double a, sumab;
	sumab = M*(1.0-M)/(SD*SD) - 1.0;
	a = M*sumab;
	return Beta(a,sumab-a);
}

double Rng::Carnes(double a0) {
	return ran_carnes(gslRng, CARNES_STD_U1, CARNES_STD_V1, CARNES_STD_U2, CARNES_STD_V2, a0);
}

double Rng::CarnesAge() {
	return ran_carnes_age(gslRng, CARNES_STD_U1, CARNES_STD_V1, CARNES_STD_U2, CARNES_STD_V2);
}

int Rng::Binomial(int n, double p) { // Binom(n,p)
	if ( p >= 0.0 && p <= 1.0 && n >= 0 ) {
		return gsl_ran_binomial(gslRng, p, unsigned(n));
	}	else {
		throw std::invalid_argument("p must be between 0 and 1 and n must be non-negative" + RIGHT_HERE);
	}
}

int Rng::PosBinomial(int n, double p) {
	/* gets a sample x from the binomial distribution
	 * conditioned on x > 0.
	 */
	// handle some trivial cases
	int x = 0;
	if ( n <= 0 || p <= 0.0 || p > 1.0 ) {
		throw std::invalid_argument("n must be positive and p must be between 0 and 1" + RIGHT_HERE);
	}
	if ( n == 1 ) return 1;
	// handle non-trivial cases
	double np = n*p;
	if ( np > 1 ) {
		// use a naive method for large np
		do {
			x = Binomial(n,p);
		} while ( x == 0 );
	} else {
		if ( np > 1e-6 ) {
			// chop-down method
			int n1 = n; double qn_aux = 1.0-p; double qn = 1.0;
			while ( n1 ) { // fast method for calculating qn := (1-p)^n
				if (n1 & 1) qn *= qn_aux;
				qn_aux *= qn_aux;  n1 >>= 1;
			}
			double pp = p / (1-p);
			double pdf = n * pp * qn; // probability of x = 1
			int ell = 1; n1 = n;
			double ran = Uniform() * (1.0-qn); // condition on x > 0
			while ( ran > pdf && ell < n ) {
				ran -= pdf;
				ell++; n1--;
				pdf *= pp * n1;
				ran *= ell;
			}
			x = ell;
		}	else {
			/* use Poisson method
			 * First make an exp(1) jump, and scale so that it is between 0 and np. Then start
			 * counting exp(1) jumps below np.
			 */
			double u0 = Uniform(); double uell(0.0);
			double expr_np = exp(-np); // don't take logs all the time...
			double expr_j0 = 1.0-u0*(1.0-expr_np); // take products, calc exp(np) and don't calc logs!!
			if ( expr_j0 <= expr_np ) {
				x = 1;
			}	else {
				int ell = 0;
				double expr_jell = expr_j0;
				while ( expr_jell > expr_np ) {
					ell++; // take another exponential jump.
					uell = Uniform();
					expr_jell *= uell;
				}
				x = ell;
			}
		}
	}
	return x;
}

int Rng::Poisson(double lambda) { // Poisson(lambda)
	if ( lambda >= 0.0 ) {
		if ( lambda == 0.0 ) {
			return 0;
		}	else {
			return gsl_ran_poisson(gslRng, lambda);
		}
	}	else {
		throw std::invalid_argument("lambda must be non-negative" + RIGHT_HERE);
	}
}

int Rng::PosPoisson(double lambda) {
  if ( lambda < 0.0 ) {
		throw std::invalid_argument("lambda must be non-negative" + RIGHT_HERE);
	}
  int x = 0.0; // return value
  if ( lambda == 0.0 ) {
    return x;
  } else {
    if ( true ) { // TODO: case for large lambda
    	// Poisson with rejection
    	do {
    		x = gsl_ran_poisson(gslRng, lambda);
    	} while ( x == 0 );
    	return x;
  	}
		// else TODO: case for small lambda
	}
}

int Rng::Zeta(double s) {
	/* dumb implementation: sample a uniform(0,zeta(s)) number u
	 * and start taking n^{-s} until u < 0
	 * TODO: use pareto and a rejection scheme
	 */
	if ( s <= 1.0 ) {
		throw std::invalid_argument("s must be > 1" + RIGHT_HERE);
	}
	double u = Uniform(0.0, gsl_sf_zeta(s));
	int n = 1;
	while ( true ) {
		u -= pow(n, -s);
		if ( u <= 0.0 ) break;
		else ++n;
	}
	return n;
}



int Rng::Hypergeometric(int n, int K, int N) { // draws, successes, population size
	if ( N >= 0 && K >= 0 && n >= 0 && K <= N && n <= N ) {
		return gsl_ran_hypergeometric(gslRng, K, N-K, n);
	}	else {
		throw std::invalid_argument("invalid parameters" + RIGHT_HERE);
	}
}


int Rng::Skellam(double lambda1, double lambda2) {
    if ( lambda1 < 0.0 || lambda2 < 0.0 ) {
			throw std::invalid_argument("lambda1 or lambda2 is negative" + RIGHT_HERE);
		}
    return Poisson(lambda1) - Poisson(lambda2);
}

int Rng::SymSkellam(int loc, double sd) {
    if ( sd < 0 ) {
			throw std::invalid_argument("sd is negative" + RIGHT_HERE);
		}
    double lambda = 0.5*sd*sd;
    return loc + Skellam(lambda, lambda);
}

void Rng::Shuffle(int* xs, int min, int n) { // a random sample (unique)
	if ( n > 0 ) {
		for ( int i = 0; i < n; ++i ) {
			xs[i] = min + i;
		}
		gsl_ran_shuffle(gslRng, xs, size_t(n), sizeof(int));
	}	else {
		throw std::invalid_argument("n must be positive" + RIGHT_HERE);
	}
}

void Rng::Dirichlet(double* xs, int n, const double* alphas) {
  if ( n > 0 ) {
    // TODO: accept zero alphas: sample is then restricted to subspace.
    gsl_ran_dirichlet(gslRng, n, alphas, xs);
  } else {
    throw std::invalid_argument("n must be positive" + RIGHT_HERE);
  }
}


// auxilliary functions

double expectedCarnesSurvival() {
	auto ecs = expected_carnes_survival(CARNES_STD_U1, CARNES_STD_V1, CARNES_STD_U2, CARNES_STD_V2);
	if ( !ecs.second ) {
		throw std::logic_error("unable to compute the expected lifespan" + RIGHT_HERE);
	}
	return ecs.first;
}
