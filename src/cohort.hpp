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

#ifndef COHORT_HPP_
#define COHORT_HPP_

#include <iostream>
#include <vector>
#include <string>

#include "macros.hpp"
#include "hashed_definitions.hpp"
#include "parameters.hpp"
#include "virus.hpp"
#include "host.hpp"
#include "immune_response.hpp"
#include "parallelism.hpp"
#include "rng.hpp"
#include "eta.hpp"

class Cohort {
public:
	Cohort(std::string ); // identifier
	~Cohort();
	void init(unsigned long , int );
	void run();
	// auxiliary
	void print(std::ostream & ) const;
protected:
	WorkerPool workerpool;
	std::vector<Host*> population;
	Virus* wildtype;
	FitnessFunction ffun;
	std::vector<MhcLocus> immuneResponses;
	Rng rng;
	std::string identifier;
};

std::ostream & operator<<(std::ostream & , const Cohort & );

void highlyDetailedInfectionHistory(unsigned long , std::string );

#endif
