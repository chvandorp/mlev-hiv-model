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

#ifndef MAIN_HPP
#define MAIN_HPP

#include <unistd.h> // for cl_args...

#include <iostream>
#include <string>
#include <sstream>
#include <map> // for coupling strings to modes
#include <stdexcept>

#include "macros.hpp"
#include "test.hpp"
#include "hashed_definitions.hpp"
#include "world.hpp"
#include "cohort.hpp"
#include "virus.hpp"

enum SimulationMode {
	SIM_HELP=0,
	ILEV_SIM,
	COHORT_SIM,
	PLEV_SIM
};

int main(int argc, char* argv[]);

// pass a seed, population size, max_time and burnin_time, an id, number of repeated epidemics
void runPlevModel(unsigned long , int , int , double , double , std::string , int );
// pass a seed, cohort size, and a version id
void runCohortModel(unsigned long , int , std::string );
void runIlevModel(unsigned long , std::string );
// for help messages (e.g. when the user passes option -h)
void printHelp();
void printLicence();
#endif
