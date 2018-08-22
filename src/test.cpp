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

#include "test.hpp"

void testLogisticModel(unsigned long seed) {
	Rng rng(seed);
	DataListAndModel dm(1);
	dm.model[0] = true; // one allele
	double pr_transmitted_allele = 0.5;
	double pr_mhc_allele = 0.3;
	double pr_escape = 0.01;
	double pr_reversion  = 0.0;
	for ( int i = 0; i < 1000; ++i ) {
		DataPoint data_point(1);
		data_point.transmitted_allele = rng.Bernoulli(pr_transmitted_allele);
		data_point.observed_allele = data_point.transmitted_allele;
		data_point.mhc_alleles[0] = (rng.Bernoulli(pr_mhc_allele) ? 1 : 0);
		if ( data_point.transmitted_allele && data_point.mhc_alleles[0] ) {
			data_point.observed_allele = rng.Bernoulli(1-pr_escape);
		} else if ( !data_point.transmitted_allele && !data_point.mhc_alleles[i] ) {
			data_point.observed_allele = rng.Bernoulli(pr_reversion);
		}
		dm.data.push_back(data_point);
	}
	auto result1 = fitLogisticModel(dm);
	std::cout << "ok:   " << result1.succes << std::endl
						<< "beta: " << result1.betas[0] << std::endl
						<< "c:    " << result1.c << std::endl
						<< "ll:   " << result1.loglike << std::endl;

	// now test the null model
	dm.model[0] = false;
	auto result0 = fitLogisticModel(dm);
	std::cout << "ok:   " << result0.succes << std::endl
						<< "beta: " << result0.betas[0] << std::endl
						<< "c:    " << result0.c << std::endl
						<< "ll:   " << result0.loglike << std::endl;

	// compute a p-value using the LRT
	double p = pLikelihoodRatioTest(result0.loglike, result1.loglike, 1);
	std::cout << "p-value: " << p << std::endl;
}



void testFET(unsigned long seed) {
	ContingencyTable ctab;
	ctab.a = 5; ctab.b = 1;
	ctab.c = 2; ctab.d = 2;
	double OR = ctab.oddsRatio().first; // ignore error
	double phi = ctab.phiCoeff().first; // ignore error
	double pval_rs = ctab.pFisherExactTest(RIGHT_TAIL).first;
	std::cout << "p-value (FET rs): " << pval_rs << std::endl;
	double pval_ls = ctab.pFisherExactTest(LEFT_TAIL).first;
	std::cout << "p-value (FET ls): " << pval_ls << std::endl;
	double pval_ts = ctab.pFisherExactTest(BOTH_TAILS).first;
	std::cout << "p-value (FET ts): " << pval_ts << std::endl;
	std::cout << "table: " << ctab << std::endl;
	std::cout << "odds ratio: " << OR << std::endl;
	std::cout << "phi coeff: " << phi << std::endl;
}


void testNkModel(unsigned long seed) {
	Rng rng(seed);
	FitnessFunction fitnessFunction;
	fitnessFunction.init(rng);

	Virus virus(&fitnessFunction, rng, Virus::SCRAMBLE);
	std::cout << "# " << correlationLength(virus, rng) << std::endl;
}

void testCarnes(unsigned long seed) {
	Rng rng(seed);
	int N = 100;
	for ( int n = 0; n < N; ++n ) {
		double ell = rng.Carnes();
		double a0 = rng.CarnesAge();
		double ell_given_a0 = rng.Carnes(a0);
		// ell and ell_given_a0 should be iid
		std::cout << ell << " " << a0 << " " << ell_given_a0 << std::endl;
	}

	double es = expectedCarnesSurvival();
	std::cout << "E(S) = " << es << " b = " << 1/es << std::endl;
}

void testLogNormal(unsigned long seed) {
	Rng rng(seed);
	int N = 1000;
	double mu = 3;
	double sigma = 0.5;
	for ( int n = 0; n < N; ++n ) {
		double x = rng.LogNormalMnSd(mu, sigma);
		std::cout << x << std::endl;
	}
}

void testZeta(unsigned long seed) {
	Rng rng(seed);
	for ( double s = 1.5; s < 5.0; s += 0.5 ) {
		std::cout << s;
		for ( int i = 0; i < 10000; ++i ) {
			int n = rng.Zeta(s);
			std::cout << "\t" << n;
		}
		std::cout << std::endl;
	}
}

void plotTP(unsigned long seed) {
    std::cout << "logVL\tbetat\tD\tTP" << std::endl;
    for ( double logVL = 0.0; logVL < 8.0; logVL += 0.1 ) {
        double VL = pow(10.0, logVL);
        double beta = frasersBeta(VL);
        double delta = frasersDelta(VL);
        std::cout << logVL << "\t" << beta << "\t" << 1/delta << "\t" << beta/delta << std::endl;
    }
}

void testSkellam(unsigned long seed) {
    Rng rng(seed);
    int loc = 10;
    double sd = 5;
    for ( int i=0; i < 10000; ++i) {
        std::cout << rng.SymSkellam(loc, sd) << std::endl;
    }
}

void testResponseHdist(unsigned long seed){
	Rng rng(seed);
	int n = 1000;
	for ( int i = 0; i < n; ++i ) {
		ImmuneResponse resp(rng);
		std::cout << resp.getH() << std::endl;
	}
}

void testQValues(unsigned long seed) {
	Rng rng(seed);
	int n = 100;
	double pi0 = 0.8;
	std::map<int, double> pValues;
	std::list<int> keys;
	for ( int i = 0; i < n; ++i ) {
		keys.push_back(i);
	}
	for ( int k : keys ) {
		if ( rng.Bernoulli(1 - pi0) ) {
			pValues[k] = rng.Beta(1, 5); // skewed towards 0
		} else {
			pValues[k] = rng.Beta(1, 1); // uniform
		}
	}
	auto qValues = computeQValues(pValues);
	keys.sort([&](int a, int b){return pValues.at(a) < pValues.at(b);});
	for ( int k : keys ) {
		std::cout << k << " " << pValues[k] << " " << qValues[k] << std::endl;
	}
}
