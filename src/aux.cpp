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

#include "aux.hpp"

/* functions defined on interfaces */

/* methods for the RandomIndexable class */

int RandomIndexable::getRandomIndex() const { return random_index; }
void RandomIndexable::setRandomIndex(Rng & rng) { random_index = rng.Integer(); }

bool compareByRandomIndex(RandomIndexable* left, RandomIndexable* right) {
	if ( left != nullptr && right != nullptr ) {
		return left->getRandomIndex() < right->getRandomIndex();
	}	else {
		if ( right == nullptr ) return true;
		else return false;
	}
}

/* methods for the Printable class */

std::ostream & operator<<(std::ostream & os, const Printable & obj) {
	obj.print(os);
	return os;
}


// auxiliary functions...

/* NB: inline functions are defined in the header.
 */

std::string ageToColor(double age) {
	std::string color;
	if ( age < MIN_CONTACT_AGE ) {
		color = "lightgray";
	} else {
		int ageclass = NUMBER_OF_AGECLASSES * (age - MIN_CONTACT_AGE) / (MAX_AGE - MIN_CONTACT_AGE);
		ageclass = std::min(NUMBER_OF_AGECLASSES-1, ageclass);
		ageclass = NUMBER_OF_AGECLASSES - ageclass; // revert colors!
		color = AGE_COLOR_MAP + std::to_string(ageclass);
	}
	return color;
}

double quick_power(double x, int n) {
	/** computes x^n by repeatedly squaring x.
	 * by definition, 0^0 = 1
	 */
	if ( n < 0 ) {
		x = 1.0/x;
		n = -n;
	}
  double result = 1.0;
  while ( n ) {
    if ( n & 1 ) {
			result *= x;
		}
    n >>= 1;
    x *= x;
  }
  return result;
}

double weighted_geometric_average(double a, double b, double w) {
	if ( a < 0.0 || b < 0.0 || w < 0.0 || w > 1.0 ) {
		throw std::logic_error("ERROR: invalid parameters" + RIGHT_HERE);
	}
	return pow(a, w) * pow(b, 1.0-w);
}

std::pair<Vec, bool> eigen_values(double a, double b, double c, double d) {
	double det = a*d - b*c;
	double tr = a+d;
	double D = tr*tr - 4*det;
	if ( D >= 0 ) {
		double sqrtD = sqrt(D);
		auto lambda = std::make_pair(0.5*(tr+sqrtD), 0.5*(tr-sqrtD));
		return std::make_pair(lambda, true); // true mean real
	} else {
		double sqrtAbsD = sqrt(abs(D));
		auto reimag = std::make_pair(0.5*tr, 0.5*sqrtAbsD);
		return std::make_pair(reimag, false); // false means complex (and not real)
	}
}

std::pair<Vec, bool> eigen_values(const Mat & M) {
	return eigen_values(M.a, M.b, M.c, M.d);
}

Vec eigen_vector_norm1(double a, double b, double c, double d, double lambda) {
	WARN_UNTESTED_FUN
	double denom = sqrt((lambda-a)*(lambda-a) + b*b);
	double u1 = b / denom;
	double u2 = (lambda - a) / denom;
	return std::make_pair(u1, u2);
}

Vec eigen_vector_norm1(const Mat & M, double lambda) {
	return eigen_vector_norm1(M.a, M.b, M.c, M.d, lambda);
}

Vec eigen_vector_sum1(double a, double b, double c, double d, double lambda) {
	double u1 = -b / (a - lambda - b);
	double u2 = (a - lambda) / (a - lambda - b);
	return std::make_pair(u1, u2);
}

Vec eigen_vector_sum1(const Mat & M, double lambda) {
	return eigen_vector_sum1(M.a, M.b, M.c, M.d, lambda);
}
