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

#include "carnes.hpp"

// the Carnes lifespan distribution

double carnes_chd(double a, void* params) {
	// NB: function returns cumulative hazard MINUS a given cumulative hazard 'ch'
	CarnesParams* cp = (CarnesParams*) params; // type cast the pointer-to-void
	double u1, u2, v1, v2, ch;
	u1 = cp->u1; u2 = cp->u2;
	v1 = cp->v1; v2 = cp->v2;
	ch = cp->ch;
	return carnes_chazard(a, u1, v1, u2, v2) - ch; // inline function
}

double carnes_h(double a, void* params) {
	CarnesParams* cp = (CarnesParams*) params; // type cast the pointer-to-void
	double u1, u2, v1, v2;
	u1 = cp->u1; u2 = cp->u2;
	v1 = cp->v1; v2 = cp->v2;
	return carnes_hazard(a, u1, v1, u2, v2); // inline function
}

void carnes_chdh(double a, void* params, double* y, double* dy) {
	CarnesParams* cp = (CarnesParams*) params; // type cast the pointer-to-void
	double u1, u2, v1, v2, ch;
	u1 = cp->u1; u2 = cp->u2;
	v1 = cp->v1; v2 = cp->v2;
	ch = cp->ch;
	*y = carnes_chazard(a, u1, v1, u2, v1) - ch; // choose a s.t. this is zero
	*dy = carnes_hazard(a, u1, v1, u2, v2); // derivative used for root finding
}

double carnes_s(double a, void* params) {
	CarnesParams* cp = (CarnesParams*) params; // type cast the pointer-to-void
	double u1, u2, v1, v2;
	u1 = cp->u1; u2 = cp->u2;
	v1 = cp->v1; v2 = cp->v2;
	return carnes_survival(a, u1, v1, u2, v2);
}

double ran_carnes(gsl_rng* rng, double u1, double v1, double u2, double v2, double a0) {
	double u = gsl_rng_uniform_pos(rng); // 0 < u < 1
	double ch = -log(u); // 0 < cl < \infty: the sampled cumulative hazard
	// if a0 > 0 is given, then return a lifespan conditioned on a given age a0
	if ( a0 > 0.0 ) {
		double ch0 = carnes_chazard(a0, u1, v1, u2, v2);
		ch += ch0; // use ch0 as a baseline cumulative hazard
	}
	// now invert the cumulative hazard to find the corresponding lifespan
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fdfsolver_type* T;
	gsl_root_fdfsolver *s;
	double x0, x = 70.0; // TODO estimate mean...
	gsl_function_fdf FDF;
	CarnesParams params = {u1, v1, u2, v2, ch};

	FDF.f = &carnes_chd;
	FDF.df = &carnes_h;
	FDF.fdf = &carnes_chdh;
	FDF.params = &params;

	T = gsl_root_fdfsolver_newton;
	s = gsl_root_fdfsolver_alloc(T);
	gsl_root_fdfsolver_set(s, &FDF, x);
	do {
		iter++;
		status = gsl_root_fdfsolver_iterate(s);
		x0 = x;
		x = gsl_root_fdfsolver_root(s);
		status = gsl_root_test_delta (x, x0, 0, 1e-3);
	} while ( status == GSL_CONTINUE && iter < max_iter );
	// TODO error handling
	gsl_root_fdfsolver_free(s);
	if ( a0 > 0.0 ) return ( x > a0 ? x : a0 );
	else return ( x > 0.0 ? x : 0.0 );
}

double carnes_compfun(double a, double u1, double v1, double u2, double v2) {
	double astar = log(1.0 + u1*M_LN2*exp(-v1))/u1; // assumes that u1 > 0
	double ch = carnes_hazard(astar, u1, v1, u2, v2); // hazard at astar (slope)
	double cch = carnes_chazard(astar, u1, v1, u2, v2); // cumulative hazard
	double a0 = astar - cch/ch;

	if ( a < a0 ) return 1.0;
	else return exp(-ch*(a-a0));
}

double carnes_ccompfun(double a, double u1, double v1, double u2, double v2) {
	/* this function equals
	 * \int_0^a carnes_compfun(a')da'
	 */
	double astar = log(1.0 + u1*M_LN2*exp(-v1))/u1; // assumes that u1 > 0
	double ch = carnes_hazard(astar, u1, v1, u2, v2); // hazard at astar (slope)
	double cch = carnes_chazard(astar, u1, v1, u2, v2); // cumulative hazard
	double a0 = astar - cch/ch;

	if ( a < a0 ) return a;
	else return a0 + (1.0 - exp(-ch*(a-a0)))/ch;
}

double carnes_ccompfun_inv(double u, double u1, double v1, double u2, double v2) {
	/* this function is the inverse of carnes_ccompfun
	 * u should be between 0 and 1.
	 */
	double astar = log(1.0 + u1*M_LN2*exp(-v1))/u1; // assumes that u1 > 0
	double ch = carnes_hazard(astar, u1, v1, u2, v2); // hazard at astar (slope)
	double cch = carnes_chazard(astar, u1, v1, u2, v2); // cumulative hazard
	double a0 = astar - cch/ch;
	double m = a0 + 1.0/ch; // upper bound

	if ( u >= 0.0 && u < 1.0 ) {
		double um = u*m; // NB: rescale u
		if ( um < a0 ) {
			return um;
		}
		else { // a0 <= um < m
			return a0 - log(1.0 - ch*(um-a0))/ch;
		}
	}
	else {
		// TODO error handling
		return 0.0;
	}
}

double ran_carnes_age(gsl_rng* rng, double u1, double v1, double u2, double v2) {
	double u, a, pdf, dpdf, p;
	while ( true ) {
		u = gsl_rng_uniform_pos(rng);
		a = carnes_ccompfun_inv(u, u1, v1, u2, v2);
		pdf = carnes_survival(a, u1, v1, u2, v2);
		dpdf = carnes_compfun(a, u1, v1, u2, v2);
		p = gsl_rng_uniform_pos(rng) * dpdf;
		if ( p < pdf ) return a;
	}
}

std::pair<double, bool> expected_carnes_survival(double u1, double v1, double u2, double v2) {
	CarnesParams params = {u1, v1, u2, v2, 0.0};
	gsl_function F;
	F.function = &carnes_s;	F.params = &params;
	double a0 = 0; // lower integration bound
	double epsabs = 0.0;
	double epsrel = 1e-7;
	size_t limit = 1000;
	gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(limit);
	double result, abserr; // passed by ref
	int err = gsl_integration_qagiu(&F, a0, epsabs, epsrel, limit, workspace, &result, &abserr);
	gsl_integration_workspace_free(workspace);
	return std::make_pair(result, err==GSL_SUCCESS);
}
