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

#ifndef TEST_HPP
#define TEST_HPP

#include <algorithm> // max

#include "virus.hpp"
#include "rng.hpp"
#include "host.hpp"
#include "immune_response.hpp"
#include "frasers_rates.hpp"
#include "hashed_definitions.hpp"
#include "parallelism.hpp"
#include "integrator.hpp"
#include "eta.hpp"
#include "calc_stats.hpp"
#include "logistic_model.hpp"

void testLogisticModel(unsigned long );
void testFET(unsigned long);
void testQValues(unsigned long);
void testSkellam(unsigned long);
void plotTP(unsigned long);
void testZeta(unsigned long);
void testIndPars(unsigned long );
void testCarnes(unsigned long );
void testLogNormal(unsigned long );
void testResponseHdist(unsigned long);


#endif /* TEST_HPP_ */
