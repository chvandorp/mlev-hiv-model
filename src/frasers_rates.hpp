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

#ifndef FRASERS_RATES_HPP
#define FRASERS_RATES_HPP

#include <cmath> // pow, isnan
#include <map>

#define FRASERS_BETA1 2.76
#define FRASERS_BETAMAX 0.317
#define FRASERS_BETA50 13938.0
#define FRASERS_BETAK 1.02
#define FRASERS_BETA50_POW_BETAK 16868.81088
#define FRASERS_BETA3 0.76
#define FRASERS_D1 0.24
#define FRASERS_DMAX 25.4
#define FRASERS_DMIN 0.1
#define FRASERS_D50 3058.0
#define FRASERS_DK 0.41
#define FRASERS_D50_POW_DK 26.85526036
#define FRASERS_D3 0.75
#define FRASERS_RHO 3.46

extern double TP_SENS_SHIFT; // multiplicative
extern double TP_SENS_BETA_SCALE;

double frasersBeta(double ); // beta(v)
double frasersDelta(double ); // 1/D(v)

#endif /* FRASERS_RATES_HPP_ */
