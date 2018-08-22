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

#include "calc_stats.hpp"

// static members

const double Statistics::f_lower = 0.025;
const double Statistics::f_quartile1 = 0.25;
const double Statistics::f_median = 0.5;
const double Statistics::f_quartile3 = 0.75;
const double Statistics::f_upper = 0.975;

// methods for Statistics

void Statistics::setName(const std::string & name) { this->name = name; }

// methods of BasicStats

BasicStats::BasicStats() : n(0), mean(0.0), var(0.0), sd(0.0), median(0.0),
	min(0.0), lower(0.0),	quartile1(0.0), quartile3(0.0),
	upper(0.0), max(0.0) { /* empty */ }

BasicStats::BasicStats(const std::string & name) : Statistics(name), n(0), mean(0.0),
	var(0.0), sd(0.0), median(0.0), min(0.0), lower(0.0), quartile1(0.0),
	quartile3(0.0),	upper(0.0), max(0.0) { /* empty */ }

void BasicStats::addPoint(double x) {	sample.push_back(x); }

void BasicStats::computeStats() {
	if ( sample.empty() ) {
		return; // don't do anything...
	}
	n = sample.size();
	sample.sort();
	// min and max ae easy after sort
	min = sample.front();
	max = sample.back();
	// advance through list for percentiles
	std::list<double>::iterator it_cache = sample.begin();
	int idx_cache = 0;
	lower = computePercentile(it_cache, idx_cache, f_lower);
	// it_cache and idx_cache are updated as a side-effect of computePercentile
	quartile1 = computePercentile(it_cache, idx_cache, f_quartile1);
	median = computePercentile(it_cache, idx_cache, f_median);
	quartile3 = computePercentile(it_cache, idx_cache, f_quartile3);
	upper = computePercentile(it_cache, idx_cache, f_upper);

	// compute moment statistics
	mean = 0.0;
	for ( auto it = sample.begin(); it != sample.end(); ++it ) {
		double & x = (*it); // alias
		mean += x;
	}
	mean /= n;
	// compute sample variance and standard deviation (only when n > 1)
	var = 0.0; sd = 0.0;
	if ( n < 2 ) {
		return;
	} // else n >= 2
	for ( auto it = sample.begin(); it != sample.end(); ++it ) {
		double & x = (*it); // alias
		var += (x - mean) * (x - mean);
	}
	var /= (n-1); // sample variance: correct for bias
	sd = sqrt(var);
}

double BasicStats::getMean() const { return mean; }
double BasicStats::getVar() const { return var; }
double BasicStats::getSd() const { return sd; }
double BasicStats::getMedian() const { return median; }
double BasicStats::getMin() const { return min; }
double BasicStats::getLower() const { return lower; }
double BasicStats::getQuartile1() const { return quartile1; }
double BasicStats::getQuartile3() const { return quartile3; }
double BasicStats::getUpper() const { return upper; }
double BasicStats::getMax() const { return max; }
int BasicStats::getSize() const { return n; }

void BasicStats::print(std::ostream & os) const {
	os << "<basic_stats name='" << name << "' "
		 << "mean='" << mean << "' "
		 << "median='" << median << "' "
		 << "min='" << min << "' "
		 << "lower='" << lower << "' "
		 << "quartile1='" << quartile1 << "' "
		 << "quartile3='" << quartile3 << "' "
		 << "upper='" << upper << "' "
		 << "max='" << max << "' "
		 << "size='" << n << "' "
		 << "/>";
}

double BasicStats::computePercentile(std::list<double>::iterator & it_cache, int & idx_cache, double f) const {
	int n = sample.size();
	if ( n == 0 ) {
		throw std::logic_error("ERROR: can't compute percentiles from an empty sample" + RIGHT_HERE);
	}
	if ( n == 1 ) {
		return sample.front();
	}
	int idx = int(floor((n-1)*f));
	double delta = (n-1)*f - idx;
	int jump = idx - idx_cache;
	std::advance(it_cache, jump);
	idx_cache = idx;
	double xl = *it_cache;
	double xr = xl;
	auto it_next = std::next(it_cache);
	if ( it_next != sample.end() ) {
		xr = *it_next;
	}
	return (1.0 - delta)*xl + delta*xr; // interpolation (cf. GSL statistics)
}

Correlation::Correlation() :
	n(0), slope(0.0), intercept(0.0), corrcoef(0.0), R2(0.0) { /* empty */ }

Correlation::Correlation(std::string name) : Statistics(name),
	n(0), slope(0.0), intercept(0.0), corrcoef(0.0), R2(0.0) { /* empty */ }

void Correlation::addPoint(double x, double y) {
	sample.push_back(std::make_pair(x, y));
}

void Correlation::computeStats() {
	// use auxiliary function to compute correlation
	computeCorr(sample, n, covariance, slope, intercept, corrcoef, R2);
	// calculate Theil-Sen slope
	for ( auto it = sample.begin(); it != sample.end(); ++it ) {
		double & x1 = it->first; double & y1 = it->second;
		for ( auto jt = sample.begin(); jt != it; ++jt ) {
			double & x2 = jt->first; double & y2 = jt->second;
			if ( x1 != x2 ) { // don't divide by 0
				theilSlope.addPoint((y2-y1)/(x2-x1));
			}
		}
	}
	theilSlope.computeStats();
}

void Correlation::bootstrap(Rng & rng, int B) {
	if ( sample.empty() ) {
		return;
	} // else: the sample is at least not empty...
	for ( int b = 0; b < B; ++b ) {
		// re-sample the data
		std::map<int, int> multiplicity_map;
		for ( int i = 0; i < int(sample.size()); ++i ) {
			int idx = rng.Integer(sample.size());
			auto pr = multiplicity_map.emplace(idx, 1);
			if ( !pr.second ) { // the key idx already existed
				pr.first->second++; // so increase the value at that key
			}
		}
		std::list<std::pair<double, double>> bsample;
		int idx = 0; // keep track of the position in the list
		for ( auto it = sample.begin(); it != sample.end(); ++it, ++idx ) {
			auto mit = multiplicity_map.find(idx);
			if ( mit != multiplicity_map.end() ) {
				std::list<std::pair<double, double>> sam(mit->second, *it);
				bsample.splice(bsample.end(), sam);
			} // else multiplicity is 0, so no samples are taken
		}
		// compute statistics (passed by ref and modified by computeCorr)
		int bn = 0;
		double bcovariance = 0.0;
		double bslope = 0.0; double bintercept = 0.0;
		double bcorrcoef = 0.0; double bR2 = 0.0;
		computeCorr(bsample, bn, bcovariance, bslope, bintercept, bcorrcoef, bR2);
		bootCorrcoef.addPoint(bcorrcoef);
		bootSlope.addPoint(bslope);
		// TODO: it's possible to compute bootstrap CIs for more statistics here
	}
	bootCorrcoef.computeStats();
	bootSlope.computeStats();
}

double Correlation::getCovariance() const { return covariance; }
double Correlation::getSlope() const { return slope; }
double Correlation::getLowerSlope() const { return bootSlope.getLower(); }
double Correlation::getUpperSlope() const { return bootSlope.getUpper(); }
double Correlation::getTheilSenSlope() const { return theilSlope.getMedian(); }
double Correlation::getLowerTheilSenSlope() const { return theilSlope.getLower(); }
double Correlation::getUpperTheilSenSlope() const { return theilSlope.getUpper(); }
double Correlation::getCorrcoef() const { return corrcoef; }
double Correlation::getLowerCorrcoef() const { return bootCorrcoef.getLower(); }
double Correlation::getUpperCorrcoef() const { return bootCorrcoef.getUpper(); }
double Correlation::getR2() const { return R2; }
double Correlation::getIntercept() const { return intercept; }
double Correlation::getSize() const { return n; }

void Correlation::print(std::ostream & os) const {
	os << "<correlation name='" << name << "' "
		 << "cov='" << covariance << "' "
		 << "slope='" << slope << "' "
		 << "lower_slope='" << bootSlope.getLower() << "' "
		 << "upper_slope='" << bootSlope.getUpper() << "' "
		 << "intercept='" << intercept << "' "
		 << "corrcoef='" << corrcoef << "' "
		 << "lower_corrcoef='" << bootCorrcoef.getLower() << "' "
 		 << "upper_corrcoef='" << bootCorrcoef.getUpper() << "' "
		 << "R2='" << R2 << "' "
		 << "size='" << n << "' "
		 << "theil_sen_slope='" << getTheilSenSlope() << "' "
		 << "lower_theil_sen_slope='" << getLowerTheilSenSlope() << "' "
		 << "upper_theil_sen_slope='" << getUpperTheilSenSlope() << "' "
     << "/>";
}

void computeCorr(const std::list<std::pair<double, double>> & sample, int & n,
		double & covariance, double & slope, double & intercept, double & corrcoef, double & R2) {
	n = sample.size();
	if ( n < 2 ) {
		covariance = 0.0;
		slope = 0.0;
		intercept = 0.0;
		corrcoef = 0.0;
		R2 = 0.0;
		return;
	} // else, compute statistics
	double Ex = 0.0; double Ey = 0.0;
	double Varx = 0.0; double Vary = 0.0;
	double Covxy = 0.0;
	double SSerr = 0.0;

	// compute the means
	for ( auto it = sample.begin(); it != sample.end(); it++ ) {
		const double & x = it->first; const double & y = it->second; // alias
		Ex += x; Ey += y;
	}
	Ex /= n; Ey /= n;

	// compute the (co)variances
	for ( auto it = sample.begin(); it != sample.end(); it++ ) {
		const double & x = it->first; const double & y = it->second; // alias
		Varx += (x - Ex) * (x - Ex);
		Vary += (y - Ey) * (y - Ey);
		Covxy += (x - Ex) * (y - Ey);
	}
	Varx /= (n-1); Vary /= (n-1); Covxy /= (n-1);
	// covariance has now been computed
	covariance = Covxy;
	// the linear model for virus load y = slope * x + intercept + NORMAL
	if ( Varx > 0.0 ) {
		slope = Covxy / Varx;
	} // else keep slope equal to 0.0
	intercept = Ey - slope * Ex;
	if ( Varx > 0.0 && Vary > 0.0 ) {
		corrcoef = Covxy / sqrt(Varx * Vary);
	}
	// compute coefficient of determination
	for ( auto it = sample.begin(); it != sample.end(); it++ ) {
		const double & x = it->first; const double & y = it->second;
		double err = y - (slope * x + intercept);
		SSerr += err * err;
	}
	SSerr /= (n-1); // because we use Vary in the denominator below...
	if ( Vary > 0.0 ) {
		R2 = 1.0 - SSerr / Vary;
	} // else keep R2 equal to 0.0
}


std::pair<double, bool> computeHatPi0(const std::vector<double> & lambdas,
		const std::vector<double> & cdf) {
	size_t n = cdf.size();
	// revert cdf
	std::vector<double> rcdf(cdf.crbegin(), cdf.crend());
	std::for_each(rcdf.begin(), rcdf.end(), [](double & x){x = 1 - x;});
	// revert lambda
	std::vector<double> rlambdas(lambdas.crbegin(), lambdas.crend());
	std::for_each(rlambdas.begin(), rlambdas.end(), [](double & x){x = 1 - x;});
	// make a weight vector
	std::vector<double> weights = rlambdas;
	std::for_each(weights.begin(), weights.end(), [](double & x){x = (1-x)*(1-x);});
	// fit a line through the rotated cdf
	double c1, cov11, sumsq; // passed by ref to GSL
	int err_code_fit = gsl_fit_wmul(rlambdas.data(), 1, weights.data(), 1, rcdf.data(), 1, n,
			&c1, &cov11, &sumsq);
	/*
	// TESTING
	std::ofstream file;
	file.open("test-qvals");
	for ( size_t i = 0; i < n; ++i ) {
		double x = rlambdas[i];
		double y, e; // passed by ref to GSL
	  gsl_fit_mul_est(x, c1, cov11, &y, &e);
		file << x << " " << rcdf[i] << " " << y << std::endl;
	}
	file.close();
	// END TESTING
	*/
	// return slope and error
	if ( err_code_fit == GSL_SUCCESS) {
		return std::make_pair(c1, true);
	}	else {
		return std::make_pair(1.0, false);
	}
}

bool ContingencyTable::isPolymorphic() const {
	return (isNonNegative() && a+b != 0 && a+c != 0 && b+d != 0 && c+d != 0);
}

bool ContingencyTable::isPositive() const {
	return (a > 0 && b > 0 && c > 0 && d > 0);
}

bool ContingencyTable::isNonNegative() const {
	return (a >= 0 && b >= 0 && c >= 0 && d >= 0);
}

void ContingencyTable::print(std::ostream & os) const {
	os << a << "," << b << ";" << c << "," << d;
};


std::pair<double, bool> ContingencyTable::oddsRatio() const {
	/** Odds Ratio (OR) for a contingency table
	 * this function can fail: in which case the second return value is false.
	 */
	if ( !isPositive() ) {
		return std::make_pair(1.0, false);
	} // else...
	double OR = double(a*d) / (b*c);
	return std::make_pair(OR, true);
}

std::pair<double, bool> ContingencyTable::phiCoeff() const {
	if ( !isPolymorphic() ) {
		return std::make_pair(0.0, false);
	} // else...
	double rows = (a+b)*(c+d); // convert to double to prevent overflow in taking product below
	double cols = (a+c)*(b+d);
	double phi = (a*d - b*c) / sqrt(rows*cols);
	return std::make_pair(phi, true);
}

std::pair<double, bool> ContingencyTable::pFisherExactTest(TailType tail) const {
	/** p-value from FET for the table
	 * --------------
	 * a=k, b   | n1
	 * c,   d   | n2
	 * --------------
	 * t,   n-t | n
	 * --------------
   */
	if ( !isPolymorphic() ) {
		return std::make_pair(1.0, false);
	}
	// TODO: make sure that there is no overflow converting long int -> int
	int n1 = a + b;
	int n2 = c + d;
	int t = a + c;
	int k = a;
	double p = 1.0; // return value
	switch ( tail ) {
		case LEFT_TAIL: {
			p = gsl_cdf_hypergeometric_P(k, n1, n2, t); // \sum_{i \leq k} p(i)
			break;
		}
		case RIGHT_TAIL: {
			p = gsl_cdf_hypergeometric_P(t-k, n2, n1, t); // reflect the table in the horizontal axis
			break;
		}
		case BOTH_TAILS: {
			/** 1) find the direction of the adjacent tail of k.
			 * 2) possibly reflect the table such that the adjacent tail is on the right.
			 * 3) take the mass under the right tail.
			 * 4) find the boundary of the opposite tail (on the left)
			 * 5) take the mass under the left tail
			 */
			// define the range of allowed k-s
			int k_min = std::max(0, t-n2);
			int k_max = std::min(t, n1);
			double pk = gsl_ran_hypergeometric_pdf(k, n1, n2, t); // the probability of k
			/* find the direction of the adjacent tail, flip values of n1, n2, and k
			 * possibilities:
			 * - if k = k_max, then adjacent tail on the right (and only contains k)
			 * - if k = k_min, then the adjacent tail is on the left (only containing k)
			 * - if p(k) > p(k+1), then the adjacent tail is on the right
			 * - if p(k) <= p(k+1), then the adjacent tail is on the left
			 */
			if ( k != k_max && (k == k_min || pk <= gsl_ran_hypergeometric_pdf(k+1, n1, n2, t)) ) {
				// the adjacent tail is on the left, so flip the table
				std::swap(n1, n2);
				k = t-k; k_max = t-k_min; k_min = t-k_max;
			}
			/* the adjacent tail is now the right rail, and the opposite tail is the left tail.
			 * now find the boundary of the left tail, i.e. find the largest ell
			 * such that p(ell) <= p(k)
			 */
			int ell_lower = k_min; // specify bounds
			int ell_upper = k;
			double p_ell_upper = pk;
			while ( ell_upper - ell_lower > 1 ) {
				int ell_mid = (ell_upper + ell_lower) / 2; // mid point
				double p_ell_mid = gsl_ran_hypergeometric_pdf(ell_mid, n1, n2, t);
				if ( p_ell_mid < pk ) { // we have found a new lower bound for ell
					ell_lower = ell_mid;
				} else { // we have found a new upper bound
					ell_upper = ell_mid;
					p_ell_upper = p_ell_mid;
				}
			}
			// now define ell and add the weight under the left tail to p
			int ell = ell_lower; // 'generic' case
			if ( p_ell_upper == pk ) { // rare case
				ell = ell_upper;
			}
			if ( k - ell > 1 ) {
				p = gsl_cdf_hypergeometric_Q(k, n1, n2, t) + pk
				  + gsl_cdf_hypergeometric_P(ell, n1, n2, t);
			} else {
				p = 1.0;
			}
			break;
		} // case BOTH_TAILS
		default: {
			throw std::invalid_argument("invalid TailType argument" + RIGHT_HERE);
			break; // redundant
		}
	}
	return std::make_pair(p, true);
}


void Association::print(std::ostream & os) const {
	os << "<association_stats "
	   << "t='" << t << "' "
		 << "p='" << p << "' "
		 << "q='" << q << "' "
		 << "X='" << X << "' "
	   << "/>";
}

double pLikelihoodRatioTest(double ll0, double ll1, int ddof) {
	/** the likelihood ratio test. ddof is the difference in degrees
	 * of freedom (i.e. number of parameters) between model 0 and model 1.
	 * ddof must be positive
	 */
	if ( ddof <= 0 ) {
		std::invalid_argument("difference of degrees of freedom must be positive" + RIGHT_HERE);
	}
	if ( ll0 >= ll1 ) {
		return 1.0;
	} else {
		return gsl_cdf_chisq_Q(-2*(ll0 - ll1), ddof);
	}
}
