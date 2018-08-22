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

#include "logistic_model.hpp"

// function f
double neg_loglike_val(const gsl_vector* v, void* p) {
  auto dm = (DataListAndModel const *) p;
  int M = dm->model.size();
  int m = std::count_if(dm->model.begin(), dm->model.end(), [](bool x){return x;});

  std::vector<double> betas(M, 0.0);
  int i = 0;
  for ( int I = 0; I < M; ++I ) {
    if ( dm->model[I] ) {
      betas[I] = gsl_vector_get(v, i++); // postfix ++
    }
  }
  double c = gsl_vector_get(v, m);

  double ll = 0.0;
  for ( auto & data_point : dm->data ) {
    double logodds = c * (2*data_point.transmitted_allele - 1);
    for ( int I = 0; I < M; ++I ) {
      logodds += betas[I] * data_point.mhc_alleles[I];
    }
    if ( data_point.observed_allele ) {
      ll += -log_one_plus_exp(-logodds); // p
    } else {
      ll += -log_one_plus_exp(logodds); // 1-p
    }
  }
  return -ll; // negative log-likelihood
}

// The gradient of f, df = (df/dx, df/dy)
void neg_loglike_grad(const gsl_vector* v, void* p, gsl_vector* df) {
  double f; // dummy variable passed by ref
  neg_loglike_val_grad(v, p, &f, df);
}

// Compute both f and df together.
void neg_loglike_val_grad(const gsl_vector* v, void* p, double* f, gsl_vector* df) {
  auto dm = (DataListAndModel const *) p;
  int M = dm->model.size();
  int m = std::count_if(dm->model.begin(), dm->model.end(), [](bool x){return x;});

  std::vector<double> betas(M, 0.0);
  int i = 0;
  for ( int I = 0; I < M; ++I ) {
    if ( dm->model[I] ) {
      betas[I] = gsl_vector_get(v, i++); // postfix ++
    }
  }
  double c = gsl_vector_get(v, m);

  std::vector<double> dlldbetas(M, 0.0);
  double dlldc = 0.0;
  double ll = 0.0; // store in f at the end

  for ( auto & data_point : dm->data ) {
    double logodds = c * (2*data_point.transmitted_allele - 1); // updated below
    double dlogoddsdc = (2*data_point.transmitted_allele - 1);
    std::vector<double> dlogoddsdbetas(M, 0.0);
    for ( int I = 0; I < M; ++I ) {
      logodds += betas[I] * data_point.mhc_alleles[I];
      dlogoddsdbetas[I] = data_point.mhc_alleles[I];
    }
    if ( data_point.observed_allele ) {
      double dl = -dlog_one_plus_exp(-logodds);
      for ( int I = 0; I < M; ++I ) {
        dlldbetas[I] += -dl * dlogoddsdbetas[I];
      }
      dlldc += -dl * dlogoddsdc;
      ll += -log_one_plus_exp(-logodds); // p
    } else {
      double dl = -dlog_one_plus_exp(logodds);
      for ( int I = 0; I < M; ++I ) {
        dlldbetas[I] += dl * dlogoddsdbetas[I];
      }
      dlldc += dl * dlogoddsdc;
      ll += -log_one_plus_exp(logodds); // 1-p
    }
  }
  // store the results for the gradient
  i = 0;
  for ( int I = 0; I < M; ++I ) {
    if ( dm->model[I] ) {
      gsl_vector_set(df, i++, -dlldbetas[I]); // postfix ++
    }
  }
  gsl_vector_set(df, m, -dlldc);
  // and store the results for the value
  *f = -ll; // function returns negative log-likelihood
}

FitResult fitLogisticModel(const DataListAndModel & dm) {
  /** fit a logistic model to virus allele data
   * using MHC alleles as predictors, and the allele of the
   * transmitted virus as a "correction" (i.e. another predictor)
   */
  int M = dm.model.size();
  int m = std::count_if(dm.model.begin(), dm.model.end(), [](bool x){return x;});
  int n = m + 1;
  gsl_multimin_function_fdf neg_loglike_func = {
    &neg_loglike_val, // f
    &neg_loglike_grad, // df
    &neg_loglike_val_grad, // fdf
    (size_t) n, // n
    (void*) &dm // params
  };
  // allocate a multimin workspace, and a vector of logodds
  gsl_multimin_fdfminimizer* s =
      gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, n);
  gsl_vector* x = gsl_vector_calloc(n); // makes all elements 0
  // initialize the multimin workspace
  gsl_multimin_fdfminimizer_set(s, &neg_loglike_func, x, 1e-3, 1e-4);
  // prepare for iterations
  int max_iter = 1000;
  int iter = 0;
  int status = GSL_CONTINUE;
  // iterate and find the minimum
  while ( status == GSL_CONTINUE && iter < max_iter ) {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);
    if ( status ) { // error
      break;
    }
    status = gsl_multimin_test_gradient(s->gradient, 1e-3);
  }
  // make a FitResult object
  FitResult result(M); // to-be-returned
  result.succes = (status == GSL_SUCCESS);
  int i = 0;
  for ( int I = 0; I < M; ++I ) {
    if ( dm.model[I] ) { // parameter is included the model...
      result.betas[I] = gsl_vector_get(s->x, i++); // postfix ++
    }
  }
  result.c = gsl_vector_get(s->x, m);
  result.loglike = -(s->f); // f is the negative log-likelihood
  // free GSL data
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);
  return result;
}

bool sufficientObservarions(const DataListAndModel & dm, int min_obs) {
  /** check that the data contains enough observations
   * to estimate all logodds
   */
  int M = dm.model.size();
  int m = std::count_if(dm.model.begin(), dm.model.end(), [](bool x){return x;});
  std::vector<int> combins(1 << (m+1), 0); // 2^(m+1) combinations
  for ( const DataPoint & dp : dm.data ) {
    // make a combinations of variables
    int idx = 0; // construct index from bits
    for ( int I = 0; I < M; ++I ) {
      if ( dm.model[I] ) {
        bool b = (dp.mhc_alleles[I] != 0);
        idx = (idx << 1) + b;
      }
    }
    // finally, add the observed allele bit
    bool b = dp.observed_allele;
    idx = (idx << 1) + b;
    // and increase the count of the combination
    combins[idx] += 1;
  }
  return std::all_of(combins.begin(), combins.end(), [&](int x){return x >= min_obs;});
}
