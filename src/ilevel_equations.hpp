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

#ifndef ILEVEL_EQUATIONS_HPP_
#define ILEVEL_EQUATIONS_HPP_

#include <cmath>
#include <iostream>

#include "rng.hpp" // for generating IndPars
#include "aux.hpp" // quick_power
#include "hashed_definitions.hpp"

/** Instead of the default parameters, hosts can have slight
 * variation in the i-level parameters.
 * This produces a more realistic variation in VL
 */

class IndPars : public Printable {
public:
	IndPars(); // standard values
	IndPars(Rng & ); // sample around the standard values
	double T0, dT, dI1, dI2, alpha, s, beta, pE, dE, gamma, c,
		f, dV, pV, lambda_naive, lambda_memo, mu, hQ, omega0, omega, Q0, dQ; // 20
	static const double T0_def, dT_def, dI1_def, dI2_def, alpha_def, s_def,
		beta0_def, beta_def, pE_def, dE_def, gamma_def, f_def, c_def,
		dV_def, pV_def, lambda_naive_def, lambda_memo_def,
		mu_def, hQ_def, omega0_def, omega_def, Q0_def, dQ_def; // 25
  void print(std::ostream & ) const;
protected:
};

/** Function to compute the default beta,
 * given a desired malthusian growth rate
 */

double defaultIlevBeta(double r0, double dI1, double dI2, double dV,
		double f, double gamma, double pV, double T0);


/** These functions are used to compute rates of the i-level model.
 * they are relatively simple functions that are called many times,
 * and are therefore mane inline
 */

inline double iLevBeta(double w_beta, const IndPars & ip) {
	return pow(ip.beta, 1.0 - VIRUS_PARAM_WEIGHT) *
		pow(ip.beta_def * w_beta, VIRUS_PARAM_WEIGHT);
}

inline double iLevDelta(double w_delta, const IndPars & ip) {
    /* We have d_I1 = d_T and d_I2 = d_T + "apoptosis"
     * The apoptosis is modulated by the virus and the host
     */
	double alpha = pow(ip.alpha, 1.0 - VIRUS_PARAM_WEIGHT) *
		pow(ip.alpha_def / w_delta, VIRUS_PARAM_WEIGHT);
	return ip.dT + alpha;
}

inline double virionDeathrate(double w_virclear, const IndPars & ip) {
    return pow(ip.dV, 1.0 - VIRUS_PARAM_WEIGHT) *
			pow(ip.dV_def / w_virclear, VIRUS_PARAM_WEIGHT);
}

inline double virionProduction(double w_virprod, const IndPars & ip) {
    return pow(ip.pV, 1.0 - VIRUS_PARAM_WEIGHT) *
			pow(ip.pV_def * w_virprod, VIRUS_PARAM_WEIGHT);
}

inline double getPcQssVirusLoad(double d_V, double p_V, double beta, double T) {
  /* We use a QSS for the virus load. The equation would be
	 * dV/dt = p_V * I2 - d_V * V - beta * V * T
	 * Hence:
	 * V = p_V * I2 / (d_V + beta * T)
	 * and V/I2 = p_V / (d_V + beta * T)
	 */
  return p_V / (d_V + beta * T);
}

inline double getQssVirusLoad(double d_V, double p_V, double beta, double I2, double T) {
	/* We use a QSS for the virus load. The equation would be
	 * dV/dt = p_V * I2 - d_V * V - beta * V * T
	 * Hence:
	 * V = I2 * p_V / (d_V + beta * T)
	 */
	return I2 * getPcQssVirusLoad(d_V, p_V, beta, T);
}

inline double mhcExpression(double w_downreg, const IndPars & ip) {
	// TODO: host effect? Better model using within-cell model
	return 1.0 / w_downreg;
}

inline double getMutationProb(double w_rt, const IndPars & ip) {
	double mu = pow(ip.mu, 1.0 - VIRUS_PARAM_WEIGHT) *
		pow(ip.mu_def / w_rt, VIRUS_PARAM_WEIGHT);
	return std::min(1.0, mu); // this limiting should almost never happen.
}

inline double getExactReplicationProb(double w_rt, const IndPars & ip) {
	/* if p is the probability of mutating a locus,
	 * then (1 - p)^N is the probability of exact replication,
	 * where N is the genome length.
	 */
	//double mu = getMutationProb(w_rt, ip);
	//return quick_power(1.0-mu, GENOME_SIZE);
	/* TESTING an alternative formula: add an extra fitness cost for higher
	 * mutation rates (e.g. by assuming many more deleterious unmodeled loci)
	 * approximate quick_power with the exponential function
	 * NB: this is no longer a probability, as it can be larger than 1. TODO: rename
	 */
	double mu_virus = getMutationProb(w_rt, ip);
	double mu_host = ip.mu;
	double mu_ratio = 1.0;
	if ( mu_host > 0 ) {
		mu_ratio = mu_virus / mu_host;
	}
	// the factor 1/2 makes keeps mu*q away from max
	return exp(-mu_ratio/2) * sqrt(M_E);
}

inline double getEclipseRate(double w_eclipse, const IndPars & ip) {
	return pow(ip.gamma, 1.0 - VIRUS_PARAM_WEIGHT) *
		pow(ip.gamma_def * w_eclipse, VIRUS_PARAM_WEIGHT);;
}

inline double getFractionNonAbortive(double w_abort, const IndPars & ip) {
	// FIXME: make sure that it is smaller than 1
	return pow(ip.f ,1.0 - VIRUS_PARAM_WEIGHT) *
		pow(ip.f_def * w_abort, VIRUS_PARAM_WEIGHT);
}

inline double getActivationRate(double w_activ, const IndPars & ip) {
	return pow(ip.omega, 1.0 - VIRUS_PARAM_WEIGHT) *
		pow(ip.omega_def * w_activ, VIRUS_PARAM_WEIGHT);
}

inline std::pair<Vec, double> iLevInitialCondition(double w_beta, double w_delta,
		double w_virprod, double w_virclear, double w_rt, double w_eclipse,
		double w_abort, const IndPars & ip) {
	/* r0 is based on the following ODEs:
	 * dI1/dt = f * pr_exact * V * T - (gamma + d_I1) * I1
	 * dI2/dt = gamma * I1 - d_I2 * I2
	 * where pr_exact is the probability of exact replication
	 * the growth rate is then the dominant eigenvalue of the jacobian J
	 * of the system, where J = (-(gamma + d_I1), B; gamma, -d_I2)
	 */
	double pr_exact = getExactReplicationProb(w_rt, ip);
	// compute beta
	double beta = iLevBeta(w_beta, ip);
	// compute delta
	double d_I2 = iLevDelta(w_delta, ip);
	// compute pc QSS Virus load
	double p_V = virionProduction(w_virprod, ip);
	double d_V = virionDeathrate(w_virclear, ip);
	double VperI2 = getPcQssVirusLoad(d_V, p_V, beta, ip.T0);
	// conversion from eclipse to producting
	double gamma = getEclipseRate(w_eclipse, ip);
	// fraction target-cell infections that is not aborted
	double f = getFractionNonAbortive(w_abort, ip);

	double B = f * pr_exact * beta * VperI2 * ip.T0;
	Mat M(-(gamma + ip.dI1), B, gamma, -d_I2);
	auto eivals_real = eigen_values(M);
	if ( !eivals_real.second ) {
		// the eigenvalues should ALWAYS be real
		throw std::logic_error("eigenvalues are not real" + RIGHT_HERE);
	}
	// get the dominant eigenvalue
	double dom_eival = std::max(eivals_real.first.first, eivals_real.first.second);
	// get the corresponding eigenvector
	auto dom_eivec = eigen_vector_sum1(M, dom_eival);
	// make sure to return a positive vector // TODO: make a Must-function
	double u1 = dom_eivec.first; double u2 = dom_eivec.second;
	if ( u1 < 0 && u2 < 0 ) {
		dom_eivec = std::make_pair(-u1, -u2);
	} else if ( u1 <= 0 || u2 <= 0 ) {
		throw std::logic_error("eigenvector is not positive" + RIGHT_HERE);
	}
	return std::make_pair(dom_eivec, dom_eival);
}

inline double iLevBasicReproductionNumber(double w_beta, double w_delta,
		double w_virprod, double w_virclear, double w_rt, double w_eclipse,
		double w_abort, const IndPars & ip) {
	/* R0 is based on the following ODEs:
	 * dI1/dt = f * beta * pr_exact * V * T0 - (d_I1 + gamma) * I1
	 * dI2/dt = gamma * I1 - d_I2 * I2
	 * where pr_exact is the probability of exact replication
	 * Hence we get the following Jacobian
	 * J = (-(d_I1 + gamma), B; gamma -d_I2)
	 * where B = f * beta * pr_exact * VperI2 * T0
	 * We write J = A - D with
	 * A = (0, B; 0, 0) and D = (d_I1 + gamma, 0; -gamma, d_I2)
	 * the next-generation matrix is then given by M = A D^{-1}
	 * M = (B *gamma, B * (d_I1+gamma); 0, 0) / ((d_I1 + gamma) * d_I2)
	 * and therefore R_0 = B / d_I2 * gamma / (d_I1 + gamma)
	 * which co-incides with the 1-stage (QSS for I1) model.
	 */
	// compute pr_exact
	double pr_exact = getExactReplicationProb(w_rt, ip);
	// compute beta
	double beta = iLevBeta(w_beta, ip);
	// compute delta
	double d_I2 = iLevDelta(w_delta, ip);
	// compute pc QSS Virus load
	double p_V = virionProduction(w_virprod, ip);
	double d_V = virionDeathrate(w_virclear, ip);
	double VperI2 = getPcQssVirusLoad(d_V, p_V, beta, ip.T0);
	// conversion from eclipse to producting
	double gamma = getEclipseRate(w_eclipse, ip);
	double c = gamma / (ip.dI1 + gamma);
		// fraction target-cell infections that is not aborted
	double f = getFractionNonAbortive(w_abort, ip);
	// return the result...
	return f * c * pr_exact * beta * VperI2 * ip.T0 / d_I2;
}

inline double iLevBasicMalthusianFitness(double w_beta, double w_delta,
		double w_virprod, double w_virclear, double w_rt, double w_eclipse,
		double w_abort, const IndPars & ip) {
	auto initial = iLevInitialCondition(w_beta, w_delta, w_virprod,
			w_virclear, w_rt, w_eclipse, w_abort, ip);
	return initial.second;
}



#endif
