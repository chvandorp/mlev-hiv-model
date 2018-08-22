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

#include "ilevel_equations.hpp"

// number of target cells when uninfected
const double IndPars::T0_def = 1e6;
// death rate of uninfected target cells
const double IndPars::dT_def = 0.1;
// death rate of infected cells in eclipse phase
const double IndPars::dI1_def = IndPars::dT_def;
// virus induced apoptosis of infected cells in virion production phase
const double IndPars::alpha_def = 0.9;
// death rate of infected cells in virion production phase
const double IndPars::dI2_def = IndPars::dT_def + IndPars::alpha_def;
// maximum growth rate of effector cells
const double IndPars::pE_def = 1.1;
// death rate of effector cells
const double IndPars::dE_def = 0.1;
// exit rate from eclipse phase
const double IndPars::gamma_def = 1.0;
// fraction of infected cells that survive the eclipse phase
const double IndPars::c_def = IndPars::gamma_def / (IndPars::gamma_def + IndPars::dI1_def);
// fraction of infected cells that survive abortive infection
const double IndPars::f_def = 0.1;
// virion clearence rate
const double IndPars::dV_def = 23.0;
// virion prodiction rate
const double IndPars::pV_def = 5e1 * IndPars::dV_def;
// convert malthusian fitness into beta_def using defaultIlevBeta
const double IndPars::beta_def = defaultIlevBeta(MALTHUSIAN_FITNESS,
		IndPars::dI1_def, IndPars::dI2_def, IndPars::dV_def, IndPars::f_def,
		IndPars::gamma_def, IndPars::pV_def, IndPars::T0_def);
// activation rate of naive T cells
const double IndPars::lambda_naive_def = 0.02;
// activation rate of memory T cells
const double IndPars::lambda_memo_def = (50 * IndPars::lambda_naive_def);
// mutation rate (per locus, per target cell infection)
const double IndPars::mu_def = 1e-6;
/** in the absence of infection, we have the following ODEs
 * dQ/dt = s - dQ * Q - omega0 * Q
 * dT/dt = omega0 * Q - dT * T
 * and hence, the steady state is
 * Q0 = s / (dQ + omega0)
 * T0 = omega0 / dT * Q0
 * Here we fix the value of T0, and hence, we get
 * Q0 = dT / omega0 * T0
 * s = Q0 * (dQ + omega0)
 */
// target cell activation rate
const double IndPars::omega0_def = 0.012; // TODO find good value (and ref)
// virus induced additional activation rate
const double IndPars::omega_def = 0.012; // TODO find good value (and ref)
// michaelis menten parameter for target cell activation
const double IndPars::hQ_def = 3.3e4; // an intermediate virus load value
// number of quescent cells when uninfected
const double IndPars::Q0_def = IndPars::T0_def * IndPars::dT_def / IndPars::omega0_def;
// deathrate of quiescent target cells
const double IndPars::dQ_def = 0.001;
// in-flux of quescent cells
const double IndPars::s_def = IndPars::Q0_def * (IndPars::dQ_def + IndPars::omega0_def);


IndPars::IndPars() {
	/* by default, choose the default parameters */
	T0 = T0_def;
	dT = dT_def;
	dI1 = dI1_def;
	dI2 = dI2_def;
	alpha = alpha_def;
	s = s_def;
	beta = beta_def;
	pE = pE_def;
	dE = dE_def;
	gamma = gamma_def;
	c = c_def;
	f = f_def;
	dV = dV_def;
	pV = pV_def;
	lambda_naive = lambda_naive_def;
	lambda_memo = lambda_memo_def;
	mu = mu_def;
	hQ = hQ_def;
	omega0 = omega0_def;
	omega = omega_def;
	Q0 = Q0_def;
	dQ = dQ_def;
}

IndPars::IndPars(Rng & rng) {
	/* use rng to produce random individual variation
	 * on the standad i-level parameters.
	 */
	double fr = ILEVEL_PARAM_VAR;
	// sample deviations to the standad parameters
	s = rng.LogNormalMnSd(s_def, fr*s_def);
	dT = rng.LogNormalMnSd(dT_def, fr*dT_def);
  beta = rng.LogNormalMnSd(beta_def, fr*beta_def);
  alpha = rng.LogNormalMnSd(alpha_def, fr*alpha_def);
	pE = rng.LogNormalMnSd(pE_def, fr*pE_def);
	dE = rng.LogNormalMnSd(dE_def, fr*dE_def);
	gamma = rng.LogNormalMnSd(gamma_def, fr*gamma_def);
	f = rng.LogNormalMnSd(f_def, fr*f_def); // TODO: should be between 0 and 1
  dV = rng.LogNormalMnSd(dV_def, fr*dV_def);
  pV = rng.LogNormalMnSd(pV_def, fr*pV_def);
  lambda_naive = rng.LogNormalMnSd(lambda_naive_def, fr*lambda_naive_def);
	lambda_memo = rng.LogNormalMnSd(lambda_memo_def, fr*lambda_memo_def);
	// NB: user could set mu_def = 0
	mu = (mu_def > 0 ? rng.LogNormalMnSd(mu_def, fr*mu_def) : 0);
	hQ = rng.LogNormalMnSd(hQ_def, fr*hQ_def);
	omega0 = rng.LogNormalMnSd(omega0_def, fr*omega0_def);
	// NB: user could set omega_def = 0
	omega = (omega_def > 0 ? rng.LogNormalMnSd(omega_def, fr*omega_def) : 0);
	dQ = rng.LogNormalMnSd(dQ_def, fr*dQ_def);
	// initiate compound parameters
	Q0 = s / (dQ + omega0);
	T0 = omega0 * Q0 / dT;
	dI1 = dT;
	dI2 = dT + alpha;
	c = gamma / (gamma + dI1);
}

void IndPars::print(std::ostream & os) const {
  os << "<indpars "
     << "T0='" << T0 << "' "
     << "dT='" << dT << "' "
     << "dI1='" << dI1 << "' "
     << "dI2='"	<< dI2 << "' "
     << "alpha='" << alpha << "' "
     << "pE='" << pE << "' "
     << "dE='" << dE << "' "
     << "gamma='" << gamma << "' "
     << "f='" << f << "' "
     << "s='" << s << "' "
     << "beta='" << beta << "' "
     << "c='" << c << "' "
     << "dV='" << dV << "' "
     << "pV='" << pV << "' "
		 << "lambda_naive='" << lambda_naive << "' "
		 << "lambda_memo='" << lambda_memo << "' "
		 << "mu='" << mu << "' "
		 << "hQ='" << hQ << "' "
 		 << "omega0='" << omega0 << "' "
		 << "omega='" << omega << "' "
		 << "Q0='" << Q0 << "' "
     << "dQ='" << dQ << "' "
     << "/>";
}

double defaultIlevBeta(double r0, double dI1, double dI2, double dV,
		double f, double gamma, double pV, double T0) {
  /* in the model with a QSS for I1, we had the following definition of beta:
	 * beta * T0 = dV * (r0 + dI2) / (f * c * pV - (r0 + dI2));
	 * with c = gamma / (dI1 + gamma)
	 * for the 2-stage model, we make use of the characteristic equation
	 * r0^2 + tr(M) * r0 + det(M) = 0,
	 * where M = (-(gamma + dI1), B; gamma, -dI2)
	 * and B = f * beta * T0 * pV / (dV + beta * T0)
	 */
	double B = (r0*r0 + (dI2 + dI1 + gamma)*r0 + dI2*(gamma + dI1)) / gamma;
	double beta = B * dV / (T0 * f * pV * (1 - B / (f * pV)));
	return beta;
}
