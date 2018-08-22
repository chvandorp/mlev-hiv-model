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

#include "frasers_rates.hpp"

// parameters for the TP-sensitivity analysis.

/** the following table gives the (multiplicative) shift on the VL scale
 * and the correction factor for beta.
 *
 * shift	mult	corr
 * ---------------------
 * -1.50	0.032	0.435
 * -1.25	0.056	0.487
 * -1.00	0.100	0.550
 * -0.75	0.178	0.629
 * -0.50	0.316	0.726
 * -0.25	0.562	0.848
 * 0.00	1.000	1.000
 * 0.25	1.778	1.190
 * 0.50	3.162	1.427
 * 0.75	5.623	1.724
 * 1.00	10.000	2.096
 * 1.25	17.783	2.561
 * 1.50	31.623	3.140
 * ---------------------
 */

double TP_SENS_STEPSIZE = 0.25; // size of the steps taken below

std::map<int, double> TPShiftMap = {
  std::make_pair(-6, 0.032), // -1.5
  std::make_pair(-5, 0.056), // -1.25
  std::make_pair(-4, 0.100), // -1.0
  std::make_pair(-3, 0.178), // -0.75
  std::make_pair(-2, 0.316), // -0.5
  std::make_pair(-1, 0.562), // -0.25
  std::make_pair(0, 1.0), // 0.0
  std::make_pair(1, 1.778), // 0.25
  std::make_pair(2, 3.162), // 0.5
  std::make_pair(3, 5.623), // 0.75
  std::make_pair(4, 10.000), // 1.0
  std::make_pair(5, 17.783),  // 1.25
  std::make_pair(6, 31.623), // 1.5
};

std::map<int, double> TPScaleMap = {
  std::make_pair(-6, 0.435), // -1.5
  std::make_pair(-5, 0.487), // -1.25
  std::make_pair(-4, 0.550), // -1.0
  std::make_pair(-3, 0.629), // -0.75
  std::make_pair(-2, 0.726), // -0.5
  std::make_pair(-1, 0.848), // -0.25
  std::make_pair(0, 1.0), // 0.0
  std::make_pair(1, 1.190), // 0.25
  std::make_pair(2, 1.427), // 0.5
  std::make_pair(3, 1.724), // 0.75
  std::make_pair(4, 2.096), // 1.0
  std::make_pair(5, 2.561),  // 1.25
  std::make_pair(6, 3.140), // 1.5
};


double TP_SENS_SHIFT = TPShiftMap.at(0);
double TP_SENS_BETA_SCALE = TPScaleMap.at(0);


double frasersBeta(double V) { // beta(v)
  // apply a TP shift...
  V /= TP_SENS_SHIFT; // \log_{10} V \mapsto \log_{10} V - \Delta SPVL
	double power1 = pow(V, FRASERS_BETAK);
	double beta = 0.0;
	if ( !std::isnan(power1) ) // test if not a NaN
		beta = FRASERS_BETAMAX * power1 / (power1 + FRASERS_BETA50_POW_BETAK);
	else
		beta = FRASERS_BETAMAX;
    // apply a correction factor (corresponding to the shift.
	return beta * TP_SENS_BETA_SCALE;
}

double frasersDelta(double V) { // 1/D(v)
	double power1 = pow(V,FRASERS_DK);
	double delta = 0.0;
	if ( !std::isnan(power1) ) // test if not a NaN
		delta = ( power1 + FRASERS_D50_POW_DK ) /
			( FRASERS_DMAX * FRASERS_D50_POW_DK + FRASERS_DMIN * power1 );
	else
		delta = FRASERS_DMIN;
	return delta;
}
