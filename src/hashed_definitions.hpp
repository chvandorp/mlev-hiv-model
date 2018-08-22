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

#ifndef HASHED_DEFINITIONS_
#define HASHED_DEFINITIONS_

#include <gsl/gsl_math.h> // math constants such as M_LOG2E, M_LN2 and M_LOG10E

/* quick and dirty "parameter file"
 * these are typical parameters that shouldn't change too much
 * TODO: merge this with parameters.hpp
 */

#define NUMBER_OF_BOOTSTRAPS 200 // minimal for computing 95% CI
#define HIGH_FDR_THRESHOLD 0.2 // 20% significant calls are false positives (False Discovery Rate)
#define HIGH_SIGNIF_THRESHOLD 0.05 // significance level of 5% (False Positive Rate)
#define LOW_FDR_THRESHOLD 0.05 // 5% significant calls are false positives
#define LOW_SIGNIF_THRESHOLD 0.001 // significance level of 0.1% (False Positive Rate)


#define DAYS_IN_YEAR 365.24 // this is a sign of O.C.D.

#define PLEV_TIMESTEP (7.0 / DAYS_IN_YEAR) // a week
#define ILEV_PER_PLEV_STEPS 7 // one per day

/* parameters determining at what intervals statistics are sampled
 * They shouldn't interfere with the model
 * (apart from the thread-based non-determinism)
 */
#define LONGITUDINAL_SAMPLE_INTERVAL 5 // year
#define LONGITUDINAL_HD_SAMPLE_INTERVAL 10 // year
#define CROSSSECTIONAL_SAMPLE_INTERVAL 20 // take a sample every x years
#define STATS_SAMPLE_INTERVAL 1 // compute stats once per x years
#define HERIT_SAMPLE_INTERVAL 10 // data for the heritability model
#define SEQUENCE_SAMPLE_INTERVAL 1 // sample a sequence every x years
#define FITNESSCOST_SAMPLE_INTERVAL 5 // print stats on locus-specific fitness costs
#define CONTACT_NETWORK_SAMPLE_INTERVAL 10 // print the contact network in the DOT language
#define ASSOCIATION_SAMPLE_INTERVAL 5 // find associations every 5 years
#define AIDS_HAZARD_SAMPLE_INTERVAL 20 // measure the incresed hazard due to pre-adaptation

#define INCIDENCE_TIME_UNIT 1.0 // 1.0: p-level incidence is measured per year
#define GSVL_SAMPLE_TIME 0.5 // time after infection for sampling GSVL (years)
#define HOST_SNAPSHOT_INTERVAL 0.5 // take host snapshots every X years during i-level sim
#define VIRUS_RECORD_THRESHOLD (4 * PLEV_TIMESTEP) // don't use all mutants for escape rate statistics

#define ILEVEL_NEUTRAL false // TODO
#define PLEVEL_NEUTRAL false // TODO
#define UNIFORM_INDIV_R 1.5 // when PLEVEL_NEUTRAL, everybody has this R_indiv.
/* choose how infection rate depends on infection phase
 * i.e. beta_1 for acute, beta_2(V) for chronic, and beta_3 for AIDS
 * OR beta_2(V) for the entire infectious period, and D1 = D3 = 0.
 */
#define TRANS_BY_STAGE_OF_INFECT false
/* We cannot infect a Host when the i-level R0 <= 1.
 * One might argue that this holds moore general and that the major outbreak
 * probability needs to be taken into account. However, Frasers Frasers beta
 * is estimated under the assumption of homogeneous and virus-independent
 * susceptibility. Hence it is unclear if infection should be i-level R0
 * dependent. To test the two scenarios, set the following macro true or false
 */
#define INFECTION_ILEV_R0_DEPENDENT false
/* in order to test a i-level neutral null model, we use Shirreffs
 * SPVL "mutation-at-transmission" model, that needs a mutation parameter
 * sigma_M from Shirreff et al. PLOS Computational Biology (2011)
 */
#define SHIRREFFS_NEUTRAL_MODEL false
#define SHIRREFF_SIGMA_M 0.12

#define NUMBER_OF_THREADS 22

#define DEFAULT_COHORT_SIZE 500 // also used for cross-sectional samples during epidemic
#define DEFAULT_POPULATION_SIZE 10000 // disease-free steady state
#define DEFAULT_BURNIN_TIME 100.0
#define DEFAULT_MAX_TIME 200.0 // epidemic starts after BURNIN_TIME and runs for MAX_TIME years
#define INITIAL_PREVALENCE 10 // p-level inoculum size
/* use a single founder virus if false,
 * else a new virus for each initial infection
 * FIXME use a better system for multiple founders and initial mutations
 */
#define MULTIPLE_PLEV_FOUNDERS false
#define NUMBER_INITIAL_MUTATIONS 0 // number of initial mutations from the wild-type virus
/* parameter for initializing evolved virus using Metropolis or
 * Simulated Annealing
 */
#define EVOLVE_METROPOLIS_STEPS 10000
/* to prevent full evaluation of the fitness function after each point mutation,
 * we use the w-parameters from the parent virus and the w-difference to
 * compute the w-s of the mutants.
 * To prevent any numerical errors to accumulate, every X generations,
 * we compute the w-s from scratch, where X is:
 */
#define PREVENT_W_CORRUPTION_INTERVAL 100

#define MEAN_CONTACT_FORMATION_RATE 2.0 // TODO what do others choose?
#define SD_CONTACT_FORMATION_RATE 1.0
#define MEAN_CONTACT_DURATION 1.0 // TODO: idem
#define SD_CONTACT_FORMATION 1.0
#define MEAN_NUMBER_TRANSMISSION_STUBS 2.0 // determines promiscuity
#define MIN_CONTACT_AGE 15.0 // TODO: what do others choose (check Herbeck)?
#define CHAR_AGE_DIFF 12.0 // determines the 'rate' of age assortative mixing
/* determines in what sense promiscuous individuals
 * make contacts with other promiscuous people
 */
#define CHAR_CONCURRENCY_DIFF 1.0
/* control the level of clustering in the network:
 * the clustering bias parameter determines how much more likely new contacts
 * are friends of friends (compared to random individuals)
 * must be beteen 0 and 1 where
 * 0 = choose random neighbors from population
 * 1 = choose neighbors from friends of friends (if possible)
 */
#define CONTACT_CLUSTERING_BIAS 0.9

#define GENOME_SIZE 540 // 540 = 60 * 9
#define GENOME_DEGREE 10 // local degree uniform between 0 and GENOME_DEGREE
#define GENOME_LOCALITY 3 // should be positive (multiplied with degree)
#define KAUFFMAN_P 0.0 // NKP model parameter (probability neutral)
#define KAUFFMAN_Q 0 // NKQ model parameter (number of values, 0 means infinity)
#define BOLTZMANN_TEMPERATURE 0.1 // determines fitness costs...
#define NUMBER_VIRUS_SAMPLES 1000 // used in the construction of FitnessFunction to determine H0
#define ILEVEL_PARAM_VAR 0.1 // i-level variation of the parameters, between 0 and 1
#define VIRUS_PARAM_WEIGHT 0.5 // each i-level parameter is partially deterined by virus and host

// TCR contact bits are used to model the non-anchor "residues" of the peptide
#define NUMBER_TCR_CONTACTS 2
// fraction of loci that can contain an epitope (don't make this too small)
#define MAX_FRACTION_EPITOPE_LOCI 1.0
#define PLOIDY 2 // :-)
#define NUMBER_MHC_LOCI 3 // A, B, C
#define NUMBER_MHC_ALLELES 25 // the number of MHC alleles per locus
#define CV_MHC_RESPONSES 0.33 // variation of MHC responses per MHC molecule
#define NUMBER_MHC_RESPONSES 10 // the number of responses per MHC allele

/* a cut-off for i-level sims.
 * stop the simulation after a maximum duration.
 * this prints a warning.
 */
#define MAX_TIME_ILEV 100
/* during the acute phase, the growth rate should be around 1.5/day:
 * but we could anticipate cripling mutations...
 */
#define MALTHUSIAN_FITNESS 1.5 // default value 1.5 gives very large R0
#define ALLOW_ILEV_EXTINCTION false // false = the last deterministic virus can't be removed
// don't track the mutants of determinstic virus with a clonesize below this threshold
#define THRESHOLD_ILEV_PARENT 5e1
/* when the number of cells drops below a certain threshold, we need to
 * consider deleting the clone.
 */
#define THRESHOLD_ILEV_ENDANGERED 1e1
/* the virus load (V) in the ODE model must be translated into a VL per ml,
 * since the ODE model describes the full body
 * For CD4 counts, we also define micro liters (uL)
 */
#define ML_IN_BLOOD 1e0 // if 1, then we have the "1 mL of blood model"
#define UL_IN_BLOOD (ML_IN_BLOOD * 1e3) // 1000 uL in a mL

// regulate plotting etc.
#define PRINT_ISOLATED_VERTICES false
#define NUMBER_OF_AGECLASSES 9
#define MAX_AGE 100 // only used for coloring, there is no max age in the model.
#define AGE_COLOR_MAP std::string("/rdylgn9/") // NB: 9 age classes
#define GENE_COLOR_MAP std::string("/pastel19/") // NB: 9 genes

// regulate the integrator/cisceme
/* the first try of the step size h used by the integrator.
 * further on, h is determined adaptively. The time unit is a day
 * TODO: print with the other parameters
 */
#define INIT_STEPSIZE_INTEGRATOR 1e-1
#define MAX_STEPSIZE_INTEGRATOR 2.0
#define INIT_STEPSIZE_CISCHEME 2.0

// fun with the terminal
#define TERM_CLEAR_LINE std::string("\033[2K")

// where to find files
#define DATA_FOLDER std::string("data/")


#endif
