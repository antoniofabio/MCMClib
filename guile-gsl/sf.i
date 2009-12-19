%{
#include <gsl/gsl_sf.h>
%}

/*SF: special functions*/

double gsl_sf_exp(const double x);
double gsl_sf_exp_mult(const double x, const double y);
double gsl_sf_expm1(const double x);
double gsl_sf_exprel(const double x);
double gsl_sf_exprel_2(const double x);
double gsl_sf_exprel_n(const int n, const double x);

double gsl_sf_lngamma(const double x);
double gsl_sf_gamma(const double x);
double gsl_sf_gammastar(const double x);
double gsl_sf_gammainv(const double x);
double gsl_sf_taylorcoeff(const int n, const double x);
double gsl_sf_fact(const unsigned int n);
double gsl_sf_doublefact(const unsigned int n);
double gsl_sf_lnfact(const unsigned int n);
double gsl_sf_lndoublefact(const unsigned int n);
double gsl_sf_lnchoose(unsigned int n, unsigned int m);
double gsl_sf_choose(unsigned int n, unsigned int m);
double gsl_sf_lnpoch(const double a, const double x);
double gsl_sf_poch(const double a, const double x);
double gsl_sf_pochrel(const double a, const double x);
double gsl_sf_gamma_inc_Q(const double a, const double x);
double gsl_sf_gamma_inc_P(const double a, const double x);
double gsl_sf_gamma_inc(const double a, const double x);
double gsl_sf_lnbeta(const double a, const double b);
double gsl_sf_beta(const double a, const double b);
double gsl_sf_beta_inc(const double a, const double b, const double x);
const double GSL_SF_GAMMA_XMAX;
const int GSL_SF_FACT_NMAX;
const int GSL_SF_DOUBLEFACT_NMAX;

double gsl_sf_hyperg_0F1(const double c, const double x);
double gsl_sf_hyperg_1F1_int(const int m, const int n, double x);
double gsl_sf_hyperg_1F1(double a, double b, double x);
double gsl_sf_hyperg_U_int(const int m, const int n, const double x);
double gsl_sf_hyperg_U(const double a, const double b, const double x);
double gsl_sf_hyperg_2F1(double a, double b, double c, double x);
double gsl_sf_hyperg_2F1_conj(double aR, double aI, double c, double x);
double gsl_sf_hyperg_2F1_renorm(double a, double b, double c, double x);
double gsl_sf_hyperg_2F1_conj_renorm(double aR, double aI, double c, double x);
double gsl_sf_hyperg_2F0(const double a, const double b, const double x);

double gsl_sf_log(const double x);
double gsl_sf_log_abs(const double x);
double gsl_sf_log_1plusx(const double x);
double gsl_sf_log_1plusx_mx(const double x);

double gsl_sf_sin(const double x);
double gsl_sf_cos(const double x);
double gsl_sf_hypot(const double x, const double y);
double gsl_sf_sinc(const double x);
double gsl_sf_lnsinh(const double x);
double gsl_sf_lncosh(const double x);
double gsl_sf_angle_restrict_symm(const double theta);
double gsl_sf_angle_restrict_pos(const double theta);

double  gsl_sf_pow_int(const double x, const int n);
