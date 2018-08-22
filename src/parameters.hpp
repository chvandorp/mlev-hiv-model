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

#ifndef PARAMETERS_HPP_
#define PARAMETERS_HPP_

/* functions to handle user-defined parameters.
 * Right now, most parameters are #defined, so all
 * this header does is define a functions that creates
 * an xml node for adding it to output files.
 * That is also why we have to include so many headers... TODO
 */

// TODO: add network constants...

#include <string>
#include <iostream>
#include <sstream>

#include "hashed_definitions.hpp"
#include "ilevel_equations.hpp" // static parameters
#include "immune_response.hpp" // static parameters
#include "frasers_rates.hpp" // TODO: print these parameters too
#include "virus.hpp" // number of genome segments

std::string xmlStringParameters();

#endif
