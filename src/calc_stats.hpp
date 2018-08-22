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

#ifndef CALC_STATS_HPP_
#define CALC_STATS_HPP_

#include <list>
#include <vector>
#include <algorithm> // for std::advance, std::count_if
#include <utility>
#include <iostream>
#include <fstream> // TESTING q-values
#include <string>
#include <map> // for re-sampling data, mapping keys to p/q-values
#include <stdexcept>
#include <gsl/gsl_bspline.h> // for q-values
#include <gsl/gsl_multifit.h> // for q-values
#include <gsl/gsl_fit.h> // for q-values
#include <gsl/gsl_randist.h> // statistical tests
#include <gsl/gsl_cdf.h> // statistical tests

#include "macros.hpp"
#include "aux.hpp"
#include "rng.hpp" // for bootstrapping

class Statistics : public Printable { // base class for all sorts of statistical computations (TODO)
public:
	Statistics() { /* empty */ } // inline constructor
	Statistics(std::string name) : name(name) { /* empty */ }
	virtual ~Statistics() { /* empty */ } // inline destructor
	void setName(const std::string & );
	// auxilliary
	virtual void print(std::ostream & ) const=0;
protected:
	std::string name;
	// parameters used to compute default percentiles
	static const double f_lower;
	static const double f_quartile1;
	static const double f_median;
	static const double f_quartile3;
	static const double f_upper;
private:
	/* empty */
};

// really basic, 1-dimensional sample statistics

class BasicStats : public Statistics {
public:
	BasicStats(); // set values to defaults
	BasicStats(const std::string & );
	void addPoint(double );
	template <class InputIterator>
  void addPoints(InputIterator first, InputIterator last);
	void computeStats();
	double getMean() const;
	double getVar() const;
	double getSd() const;
	double getMedian() const;
	double getMin() const;
	double getLower() const;
	double getQuartile1() const;
	double getQuartile3() const;
	double getUpper() const;
	double getMax() const;
  int getSize() const;
	// auxilliary
	void print(std::ostream & ) const; // xml node with attributes
protected:
	int n; // the sample size (set to sample.size() when computeStats is called)
	double mean;
	double var;
	double sd;
	double median; // quartile2
	double min;
	double lower; // 2.5%ile
	double quartile1; // 25%ile
	double quartile3; // 75%ile
	double upper; // 97.5%ile
	double max;
	std::list<double> sample;
	// auxiliary methods
	double computePercentile(std::list<double>::iterator & it_cache, int & idx_cache, double f) const;
};

// template methods

template <class InputIterator>
void BasicStats::addPoints(InputIterator first, InputIterator last) {
	for ( InputIterator it = first; it != last; ++it ) {
		addPoint(*it);
	}
}


// compute correlations, heritability

class Correlation : public Statistics {
public:
	Correlation();
	Correlation(std::string name);
	void addPoint(double x, double y);
	template <class InputIterator>
	void addPoints(InputIterator first, InputIterator last);
	void computeStats();
	void bootstrap(Rng & rng, int B);
	double getCovariance() const;
	double getSlope() const;
	double getLowerSlope() const;
	double getUpperSlope() const;
	double getTheilSenSlope() const;
	double getLowerTheilSenSlope() const;
	double getUpperTheilSenSlope() const;
	double getCorrcoef() const;
	double getLowerCorrcoef() const;
	double getUpperCorrcoef() const;
	double getR2() const;
	double getIntercept() const;
	double getSize() const;
	void print(std::ostream & os) const; // xml node with attributes
protected:
	int n;
	double slope;
	double intercept;
	double corrcoef;
	double covariance;
	double R2;
	BasicStats theilSlope;
	BasicStats bootSlope;
	BasicStats bootCorrcoef;
	std::list<std::pair<double, double>> sample;
	// auxiliary function for computeStats and bootstrap
};

void computeCorr(const std::list<std::pair<double, double>> & sample,
	int & n, double & covariance, double & slope, double & intercept,
	double & corrcoef, double & R2);

template <class InputIterator>
void Correlation::addPoints(InputIterator first, InputIterator last) {
	for ( auto it = first; it != last; ++it ) {
		addPoint(it->first, it->second);
	}
}

/** computeQValues is a function that transforms p-value to q-values
 * it uses an auxiliary funtion computeHatPi0 that takes care of the GSL stuff
 * in order to estimate pi_0, the fraction of truly null features.
 */
std::pair<double, bool> computeHatPi0(const std::vector<double> & , const std::vector<double> & );

template<class KeyType>
std::map<KeyType, double> computeQValues(const std::map<KeyType, double> & pValues) {
	int m = pValues.size();
	// sort keys by p-values
	std::list<KeyType> keys;
	for ( auto pair : pValues ) {
		keys.push_back(pair.first);
	}
	keys.sort([&](KeyType lk, KeyType rk){return pValues.at(lk) < pValues.at(rk);});
	// compute \hat{\pi}_0(\lambda) for \lambda = 0, 0.01, 0.02,...
	int n = 20; // 0,1,...,19 TODO: choose good value (based on m)
	double mesh = 1.0 / n;
	std::vector<double> lambdas(n, 0.0);
	for ( int i = 0; i < n; ++i ) {
		lambdas[i] = i * mesh;
	}
	std::vector<double> cdf(n, 0.0);
	for ( int i = 0; i < n; ++i ) {
		auto pred = [&](std::pair<KeyType, double> pair){return pair.second < lambdas.at(i);};
		cdf[i] = std::count_if(pValues.begin(), pValues.end(), pred) / double(m);
	}
	auto pair = computeHatPi0(lambdas, cdf);
	if ( !pair.second ) {
		throw std::runtime_error("unable to estimate pi_0" + RIGHT_HERE);
	} // else: estimation of pi_0 went fine.
	double hatPi0 = pair.first;
	// compute the q-values, by starting with the largest p-value
	std::map<KeyType, double> qValues; // to-be-returned
	double q_previous = 1.0;
	keys.reverse();
	int i = m; // index in the for loop over the keys
	for ( auto k : keys ) {
		q_previous = std::min(q_previous, (m * hatPi0 * pValues.at(k)) / (i--));
		// postfix ++/-- returns the old value
		qValues[k] = q_previous;
	}
	return qValues;
}

/** compute odds ratios and significance */

enum TailType {
	RIGHT_TAIL=0, // null is greater or equal than statistic
	LEFT_TAIL, // null is smaller or equal than statistic
	BOTH_TAILS // use both tails of the distribution
};

class ContingencyTable : public Printable {
	/* contingency table with elements
	 * --------------
	 * a=k, b   | n1
	 * c,   d   | n2
	 * --------------
	 * t,   n-t | n
	 * --------------
	 */
public:
	ContingencyTable() : a(0), b(0), c(0), d(0) {}
	ContingencyTable(long int a, long int b, long int c, long int d) : a(a), b(b), c(c), d(d) {}
	long int a, b, c, d; // max at lease 2^31 - 1: we have to take products
	bool isPolymorphic() const; // i.e. no rows or columns add up to 0
	bool isPositive() const; // i.e. none of elements are 0
	bool isNonNegative() const; // i.e. all elements are >= 0
	std::pair<double, bool> oddsRatio() const;
	std::pair<double, bool> phiCoeff() const;
	std::pair<double, bool> pFisherExactTest(TailType ) const;
	void print(std::ostream & ) const; // TODO: make xml node
};

struct Association : public Printable {
	Association() : t(0.0), p(1.0), q(1.0), X(0.0) {}
	Association(double t, double p, double q, double X) : t(t), p(p), q(q), X(X) {}
	double t, p, q, X; // time, p-value, q-value, statistic
	void print(std::ostream & ) const; // xml node
};

template <class Numeric>
std::pair<double, double> pMannWhitneyUTest(const std::list<Numeric> & xs,
			const std::list<Numeric> & ys, TailType tail) {
	/** Mann-Whithney U test using sum of ranks and Normal approximation.
	 * The first argument is the p-value.
	 * The second argument is the U statistic, scaled by the sample sizes
	 * U = Pr(x < y)
	 * TODO: correct for ties, check direction.
	 */
	// handle undefined cases
	if ( xs.empty() || ys.empty() ) {
		throw std::invalid_argument("one or more lists is empty" + RIGHT_HERE);
	} // else: xs and ys are not empty
	// make list of tagged values
	std::list<std::pair<Numeric, bool>> tagged_values;
	for ( Numeric x : xs ) {
		tagged_values.push_back(std::make_pair(x, true)); // give xs tag 'false'
	}
	for ( Numeric y : ys ) {
		tagged_values.push_back(std::make_pair(y, false)); // give ys tag 'true'
	}
	// sort tagged values
	auto sort_fun = [](const std::pair<Numeric, bool> & a, const std::pair<Numeric, bool> & b) {
		return a.first < b.first;
	};
	tagged_values.sort(sort_fun);
	// sum the ranks
	long int R1 = 0;
	int rank = 1;
	for ( auto & pair : tagged_values ) {
		if ( pair.second ) { // 'true', so value comes from xs
			R1 += rank;
		}
		rank++; // TODO: correct for ties
	}
	long int n1 = xs.size();
	long int n2 = ys.size();
	long int U1 = R1 - (n1*(n1+1))/2;
	double meanU = 0.5 * n1 * n2;
	double sdU = sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0);
	// TODO: correct for ties
	double z = (U1 - meanU) / sdU; // sdU is never 0
	// compute a p-value
	double p = 1.0;
	switch ( tail ) {
		case LEFT_TAIL: {
			p = gsl_cdf_ugaussian_P(z);
			break;
		}
		case RIGHT_TAIL: {
			p = gsl_cdf_ugaussian_Q(z);
			break;
		}
		case BOTH_TAILS: {
			// determine direction...
			if ( z < 0.0 ) { // left tail is adjacent
				p = 2.0 * gsl_cdf_ugaussian_P(z);
				// factor 2 accouts for opposite tail
			} else { // right tail is adjacent
				p = 2.0 * gsl_cdf_ugaussian_Q(z);
				// factor 2 accouts for opposite tail
			}
			break;
		}
		default: {
			throw std::invalid_argument("invalid TailType given" + RIGHT_HERE);
			break; // redundant
		}
	} // switch
	// scale U1 so that it can be interpreted as a probability
	double U = double(U1) / (n1*n2);
	return std::make_pair(p, U);
}

double pLikelihoodRatioTest(double ll0, double ll1, int ddof);

#endif /* CALC_STATS_HPP_ */
