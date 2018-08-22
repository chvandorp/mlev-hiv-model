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

#include "main.hpp"


int main(int argc, char* argv[]) {
	/* put everything in a try block, and catch exceptions below.
	 * TODO: move the cl argument stuff to an auxiliary function
	 */
	try {
		printLicence();
		SimulationMode mode = SIM_HELP;
		std::map<std::string, SimulationMode> modeMap = {
			std::make_pair("plev", PLEV_SIM),
			std::make_pair("ilev", ILEV_SIM),
			std::make_pair("cohort", COHORT_SIM),
			std::make_pair("help", SIM_HELP)
		};

		// make sense of the command-line options...
		int opt;
		unsigned long seed = 1729;
		std::string identifier = "test";
		// burnin_time: for getting a virgin population in steady state
		double burnin_time = DEFAULT_BURNIN_TIME;
		// max_time: length of actual simulation (after burnin)
		double max_time = DEFAULT_MAX_TIME;
		int population_size = DEFAULT_POPULATION_SIZE;
		int cohort_size = DEFAULT_COHORT_SIZE;
		int num_worlds = 1; // by default, only one World...
		opterr = 0;
		std::string mode_string = "help"; // will contain the passed mode
		bool cl_panic = false; // -h overrules other mode selections, also errors


		std::string argflags = "m:s:i:l:b:n:c:w:"; // need to have a colon (:)
		std::string switchflags = "h";

		/* TODO: use getopt_long with e.g.
    	struct option long_options[] = {
        {"mode",       required_argument, 0,  'm' },
        {"popsize",    required_argument, 0,  'n' },
        {"cohortsize", required_argument, 0,  'c' },
        {"length",     required_argument, 0,  'l' },
        {"burnin",     required_argument, 0,  'b' },
        {"id",         required_argument, 0,  'i' },
        {0,            0,                 0,   0  }
    	};
    */

		while ( (opt = getopt(argc, argv, (switchflags + argflags).c_str())) != -1 ) {
			switch ( opt ) {
				case 'm': { // select mode
					mode_string = std::string(optarg);
					// only select mode when the flag -h is not passed
					if ( !cl_panic && modeMap.find(mode_string) != modeMap.end() ) {
						mode = modeMap[mode_string];
					}
					break;
				}
				case 'h': {
					mode = SIM_HELP;
					cl_panic = true;
					break;
				}
				case 'i': {
					identifier = std::string(optarg);
					break;
				}
				case 's': {
					std::stringstream ss(optarg);
					ss >> seed;
					break;
				}
				case 'l': {
					std::stringstream ss(optarg);
					ss >> max_time;
					break;
				}
				case 'b': {
					std::stringstream ss(optarg);
					ss >> burnin_time;
					break;
				}
				case 'n': {
					std::stringstream ss(optarg);
					ss >> population_size;
					break;
				}
				case 'c': {
					std::stringstream ss(optarg);
					ss >> cohort_size;
					break;
				}
				case 'w': {
					std::stringstream ss(optarg);
					ss >> num_worlds;
					break;
				}
				case '?': {
					if ( argflags.find(optopt) != std::string::npos ) {
						std::cerr << "# ERROR: option -" << char(optopt) << " requires an argument" << std::endl;
					}	else {
						if ( isprint(optopt) ) {
							std::cerr << "# ERROR: unknown option `-" << char(optopt) << "'" << std::endl;
						}	else {
							std::cerr << "# ERROR: unknown option character" << std::endl;
						}
					}
					mode = SIM_HELP;
					cl_panic = true; // make sure mode is not adjusted
					break;
				}
				default: {
					std::cerr << "# ERROR: unexpected response from getopt" << std::endl;
					mode = SIM_HELP;
					cl_panic = true; // make sure mode is not adjusted
					break;
				}
      }
    }

		// print the interpretation of the command line options

		if ( mode != SIM_HELP ) {
			std::cout << "about to start simulation..." << std::endl
		          	<< "seed        = " << seed << std::endl
		          	<< "id          = " << identifier << std::endl
		          	<< "mode        = " << mode_string << std::endl
		          	<< "burnin      = " << burnin_time << " years" << std::endl
		          	<< "length      = " << max_time << " years" << std::endl
		          	<< "pop size    = " << population_size << std::endl
		          	<< "cohort size = " << cohort_size << std::endl
								<< "num worlds  = " << num_worlds << std::endl
								<< "TP shift    = " << TP_SENS_SHIFT << std::endl // TODO: cl argument
								<< "TP scale    = " << TP_SENS_BETA_SCALE << std::endl; // TODO: cl argument
			std::cout << "> ok? (press enter)";
			std::string str;
			std::getline(std::cin, str); // accepts characters until return
		}

		// run different routines, depending on the command line imput

 		switch ( mode ) {
			case ILEV_SIM: {
				runIlevModel(seed, identifier);
				break;
			}
			case COHORT_SIM: {
				runCohortModel(seed, cohort_size, identifier); // TODO cohort size argument
				break;
			}
			case PLEV_SIM: {
				runPlevModel(seed, population_size, cohort_size, max_time, burnin_time, identifier, num_worlds);
				break;
			}
			case SIM_HELP: {
				printHelp();
				break;
			}
			default: {
				std::cerr << "ERROR: unknown mode." << std::endl;
				break;
			}
		}
		std::cout << "goodbye!" << std::endl;
	} // main try block
	catch ( std::exception & ex ) {
		std::cerr << "ERROR: " << ex.what() << std::endl;
	}
	return 0;
}

void runPlevModel(unsigned long seed, int population_size, int cohort_size,
		double max_time, double burnin_time, std::string identifier, int num_worlds) {
	Rng rng(seed);
	World* old_world = nullptr;
	for ( int i = 0; i < num_worlds; ++i ) {
		unsigned long new_seed = rng.Integer();
		World* new_world = nullptr;
		if ( old_world == nullptr ) {
			new_world = new World(identifier, new_seed, population_size, cohort_size, -burnin_time);
		} else {
			new_world = new World(*old_world, new_seed, -burnin_time);
			delete old_world;
		}
		new_world->run(0.0, false); // burn-in, (false means: don't keep epidemic alive)
		new_world->run(max_time, true); // (true means: keep epidemic alive)
		std::cout << std::endl;
		old_world = new_world;
	}
	delete old_world;
}

void runCohortModel(unsigned long seed, int cohort_size, std::string identifier) {
	Cohort cohort(identifier);
	cohort.init(seed, cohort_size);
	cohort.run();

	std::string fileName = DATA_FOLDER + "cohort-file-" + identifier + ".xml";
	std::ofstream fileHandle; fileHandle.open(fileName.c_str());
	if ( fileHandle ) {
		fileHandle << cohort;
		fileHandle.close();
	}	else {
		std::cerr << "ERROR: can't open output file for cohort" << std::endl;
	}
}

void runIlevModel(unsigned long seed, std::string identifier) {
	highlyDetailedInfectionHistory(seed, identifier);
}

void printLicence() {
	const std::string licence =
		"\nmlev-hiv-sim  Copyright (C) 2017  Christiaan H. van Dorp\n\n"
		"This is free software, covered by the GNU General Public License,\n"
		"and comes WITHOUT ANY WARRANTY WHATSOEVER.\n";
	std::cout << licence << std::endl;
 }

void printHelp() {
	const std::string helpmsg =
		"usage: $ ./mlev-hiv-sim [-h] [-m {ilev,cohort,plev,help}]\n"
		"                        [-s SEED] [-i ID] [-l LENGTH] [-b BURNIN]\n"
		"                        [-n POPULATION_SIZE] [-c COHORT_SIZE]\n"
		"\n"
		"Run individual-level and population-level HIV-1 simulations,\n"
		"using the NK-model for genotype-to-phenotype mapping,\n"
		"an ODE model with stochastic dimension for immuno- and virus dynamics,\n"
		"and a dynamic network model with assortative mixing at the epidemic level.\n"
		"\n"
		"optional arguments:\n"
		"  -h            Print this help message\n"
		"  -m {ilev,cohort,plev,help}\n"
		"                Select simulation mode... 'ilev': individual-level\n"
		"                simulation with highly detailed output. 'cohort':\n"
		"                a large number of individual-level simulations with\n"
		"                less detailed output. 'plev': simulate an epidemic.\n"
		"                'help': print this help message.\n"
		"  -s SEED\n"
		"                Give a seed for the random number generator.\n"
		"                This must be a non-negative integer.\n"
		"  -i ID\n"
		"                A string used to identify output file names.\n"
		"  -l LENGTH\n"
		"                The length of the epidemic (in years).\n"
		"  -b BURNIN\n"
		"                The length of the pre-epidemic phase, used to establish\n"
		"                a stable contact distribution.\n"
		"  -n POPULATION_SIZE\n"
		"                The steady-state population size of the virgin population.\n"
		"                This must be a positive integer.\n"
		"  -w NUMBER_OF_WORLDS\n"
		"                The number of consecutive worlds. The founder virus of the\n"
		"                next world is sampled from the end of the last simulation.\n"
		"                This must be a positive integer.\n"
		"  -c COHORT_SIZE\n"
		"                The size of the cohort used in mode 'cohort'. The cohort\n"
		"                size is also used in mode 'plev' for cross-sectional samples.\n"
		"                This must be a positive integer.\n";
	std::cout << helpmsg << std::endl;
}
