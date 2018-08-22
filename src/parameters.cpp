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

#include "parameters.hpp"

// TODO: add other parameters, handle cl arguments
std::string xmlStringParameters() {
	std::stringstream xs;

	xs << "<parameters "
		<< "PLEVEL_NEUTRAL='" << std::boolalpha << PLEVEL_NEUTRAL << std::boolalpha << "' "
		<< "ILEVEL_NEUTRAL='" << std::boolalpha << ILEVEL_NEUTRAL << std::boolalpha <<"' "
		<< "TRANS_BY_STAGE_OF_INFECT='" << std::boolalpha << TRANS_BY_STAGE_OF_INFECT << std::noboolalpha << "' "
		<< "TP_SENS_SHIFT='" << TP_SENS_SHIFT << "' "
		<< "TP_SENS_BETA_SCALE='" << TP_SENS_BETA_SCALE << "' "
		<< "PLEV_TIMESTEP='" << PLEV_TIMESTEP << "' "
		<< "ILEV_PER_PLEV_STEPS='" << ILEV_PER_PLEV_STEPS << "' "
		<< "ILEVEL_PARAM_VAR='" << ILEVEL_PARAM_VAR << "' "
		<< "VIRUS_PARAM_WEIGHT='" << VIRUS_PARAM_WEIGHT << "' "
		<< "beta_def='" << IndPars::beta_def << "' "
		<< "alpha_def='" << IndPars::alpha_def << "' "
		<< "dI2_def='" << IndPars::dI2_def << "' "
		<< "pE_def='" << IndPars::pE_def << "' "
		<< "dE_def='" << IndPars::dE_def << "' "
		<< "k_dom_def='" << ImmuneResponse::k_dom_def << "' "
		<< "k_sub_def='" << ImmuneResponse::k_sub_def << "' "
		<< "h_dom_def='" << ImmuneResponse::h_dom_def << "' "
		<< "h_sub_def='" << ImmuneResponse::h_sub_def << "' "
		<< "prob_dom='" << ImmuneResponse::prob_dom << "' "
		<< "dT_def='" << IndPars::dT_def << "' "
		<< "s_def='" << IndPars::s_def << "' "
		<< "gamma_def='" << IndPars::gamma_def << "' "
		<< "c_def='" << IndPars::c_def << "' "
		<< "dI1_def='" << IndPars::dI1_def << "' "
		<< "T0_def='" << IndPars::T0_def << "' "
		<< "dV_def='" << IndPars::dV_def << "' "
		<< "pV_def='" << IndPars::pV_def << "' "
		<< "lambda_naive_def='" << IndPars::lambda_naive_def << "' "
		<< "lambda_memo_def='" << IndPars::lambda_memo_def << "' "
		<< "mu_def='" << IndPars::mu_def << "' "
		<< "omega0_def='" << IndPars::omega0_def << "' "
		<< "omega_def='" << IndPars::omega_def << "' "
		<< "dQ_def='" << IndPars::dQ_def << "' "
		<< "Q0_def='" << IndPars::Q0_def << "' "
		<< "MALTHUSIAN_FITNESS='" << MALTHUSIAN_FITNESS << "' "
		<< "NUMBER_TCR_CONTACTS='" << NUMBER_TCR_CONTACTS << "' "
		<< "MAX_FRACTION_EPITOPE_LOCI='" << MAX_FRACTION_EPITOPE_LOCI << "' "
		<< "PLOIDY='" << PLOIDY << "' "
		<< "NUMBER_MHC_LOCI='" << NUMBER_MHC_LOCI << "' "
		<< "NUMBER_MHC_ALLELES='" << NUMBER_MHC_ALLELES << "' "
		<< "NUMBER_MHC_RESPONSES='" << NUMBER_MHC_RESPONSES << "' "
		<< "GENOME_SIZE='" << GENOME_SIZE << "' "
		<< "GENOME_DEGREE='" << GENOME_DEGREE << "' "
		<< "GENOME_LOCALITY='" << GENOME_LOCALITY << "' "
		<< "NUMBER_OF_GENES='" << NUMBER_OF_GENES << "' "
		<< "KAUFFMAN_P='" << KAUFFMAN_P << "' "
		<< "KAUFFMAN_Q='" << KAUFFMAN_Q << "' "
		<< "BOLTZMANN_TEMPERATURE='" << BOLTZMANN_TEMPERATURE << "' "
		<< "INITIAL_PREVALENCE='" << INITIAL_PREVALENCE << "' "
		<< "LONGITUDINAL_SAMPLE_INTERVAL='" << LONGITUDINAL_SAMPLE_INTERVAL << "' "
		<< "CROSSSECTIONAL_SAMPLE_INTERVAL='" << CROSSSECTIONAL_SAMPLE_INTERVAL << "' "
		<< "STATS_SAMPLE_INTERVAL='" << STATS_SAMPLE_INTERVAL << "' "
		<< "HERIT_SAMPLE_INTERVAL='" << HERIT_SAMPLE_INTERVAL << "' "
		<< "SEQUENCE_SAMPLE_INTERVAL='" << SEQUENCE_SAMPLE_INTERVAL << "' "
		<< "THRESHOLD_ILEV_PARENT='" << THRESHOLD_ILEV_PARENT << "' "
		<< "THRESHOLD_ILEV_ENDANGERED='" << THRESHOLD_ILEV_ENDANGERED << "' "
		<< "DAYS_IN_YEAR='" << DAYS_IN_YEAR << "' "
		<< "NUMBER_INITIAL_MUTATIONS='" << NUMBER_INITIAL_MUTATIONS << "' "
		<< "ML_IN_BLOOD='" << ML_IN_BLOOD << "' "
		<< ">" << std::endl;
	   // maybe add some comments as plain text?
	xs << "</parameters>";

	return xs.str();
}
