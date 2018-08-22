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

#ifndef CARNES_HPP_
#define CARNES_HPP_

#include <utility> // pair
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>

#include "macros.hpp"

#define CARNES_STD_U1 0.1
#define CARNES_STD_V1 -10.5
#define CARNES_STD_U2 -0.4
#define CARNES_STD_V2 -8

inline double carnes_hazard(double a, double u1, double v1, double u2, double v2) {
	return exp(u1*a+v1) + exp(u2*a+v2);
}
inline double carnes_chazard(double a, double u1, double v1, double u2, double v2) {
	return exp(u1*a+v1)/u1 + exp(u2*a+v2)/u2 - (exp(v1)/u1 + exp(v2)/u2);
}
inline double carnes_survival(double a, double u1, double v1, double u2, double v2) {
	return exp(-carnes_chazard(a, u1, v1, u2, v2));
}

struct CarnesParams { double u1, v1, u2, v2, ch; }; // needed for root finding
double carnes_chd(double a, void* params); // cumulative hazard - given ch
double carnes_h(double a, void* params); // the derivative of the cumulative hazard!
double carnes_s(double a, void* params); // survival; ignores params->ch
void carnes_chdh(double a, void* params, double* y, double* dy); // hazard and cumul hazard diff at once...
double ran_carnes(gsl_rng* rng, double u1, double v1, double u2, double v2, double ca=0.0); // life span

double carnes_compfun(double a, double u1, double v1, double u2, double v2);
double carnes_ccompfun(double a, double u1, double v1, double u2, double v2);
double carnes_ccompfun_inv(double u, double u1, double v1, double u2, double v2);
double ran_carnes_age(gsl_rng* rng, double u1, double v1, double u2, double v2);

/* expected_carnes_survival returns \int_0^{\infty} exp(-\Lambda(a))da,
 * where \Lambda is the cumulative hazard
 */
std::pair<double, bool> expected_carnes_survival(double u1, double v1, double u2, double v2);

#endif
