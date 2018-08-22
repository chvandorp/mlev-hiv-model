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

#ifndef LOGISTIC_MODEL_HPP_
#define LOGISTIC_MODEL_HPP_

#include <iostream> // TESTING
#include <list>
#include <vector>
#include <algorithm> // count_if, accumulate
#include <gsl/gsl_multimin.h>

/* to have a valid statistical test, we need for each combination at least
 * a certain number of observations
 */
#define MIN_TEST_OBSERVATIONS 5

// data struct
struct DataPoint {
  DataPoint() : transmitted_allele(false), observed_allele(false) {}
  DataPoint(int M) : transmitted_allele(false), observed_allele(false), mhc_alleles(M, 0) {}
  bool transmitted_allele;
  bool observed_allele;
  std::vector<int> mhc_alleles; // 0, 1, (2?), or maybe an index
};

struct DataListAndModel {
  DataListAndModel() {}
  DataListAndModel(int M) : model(M, true) {} // by default, select the full model
  std::list<DataPoint> data;
  std::vector<bool> model; // select a submodel
};

struct FitResult {
  FitResult() : succes(false), c(0.0), loglike(0.0) {}
  FitResult(int M) : succes(false), c(0.0), loglike(0.0), betas(M, 0.0) {}
  bool succes; // the minimization algorithm worked
  double c; // weight of transmitted_allele
  double loglike; // the log-likelihood
  std::vector<double> betas; // weight of mhc_alleles
};

// auxiliary function log(1 + exp(x))
inline double log_one_plus_exp(double x) {
  if ( x > 0.0 ) { // prevent overflow
    return x + log(1.0 + exp(-x)); // = log(exp(x)(1+exp(-x))) = log(exp(x) + 1)
  } else {
    return log(1.0 + exp(x)); // x <= 0, hence exp(x) is small
  }
}
// and the derivative of log(1+exp(x))
inline double dlog_one_plus_exp(double x) {
  if ( x > 0.0 ) {
    return 1.0 / (1.0 + exp(-x));
  } else { // x <= 0, so exp(x) is small
    double y = exp(x);
    return y / (1.0 + y);
  }
}

// function f
double neg_loglike_val(const gsl_vector* v, void* params);
// The gradient of f, df = (df/dx, df/dy)
void neg_loglike_grad(const gsl_vector* v, void* params, gsl_vector* df);
// Compute both f and df together.
void neg_loglike_val_grad(const gsl_vector* v, void* params, double* f, gsl_vector* df);

FitResult fitLogisticModel(const DataListAndModel & );

bool sufficientObservarions(const DataListAndModel & , int min_obs=MIN_TEST_OBSERVATIONS);

#endif
