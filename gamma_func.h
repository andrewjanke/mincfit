/* 
 *  gamma_func.h
 *
 *	This file provides the code vector_fit_function and minc_fit_function with all of 
 *	the details about fitting a gamma variate function to a timecourse.
 *
 *	See equation 5 of Harrer, 2004 for details.
 */


/*
 *	NVARY is the number of parameters which are optimised in the function fitting.
 *
 *	NPARM is the number of additional parameters passed in from the command line which
 *	are used in the function fitting (eg. additional filenames or doubles).
 *
 *	PARM_LIST is a string containing the names of these input arguments, this is used 
 *	for reporting errors in parsing the command line.
 *
 *	SPECIAL_REC is the record of additional constants only needed for this particular
 *	function.
 *
 *	For this function the only additional input is the time_step (the time between 
 *	acquisitions).
 */

#ifndef GAMMA_FUNC_H
#define GAMMA_FUNC_H

#define FNAME_SIZE 10

#define GAMMA_NVARY       4
#define GAMMA_NCONST      0
#define GAMMA_COMMENT     "use the gamma function\n                f = Cmax * ((t-t0)^alpha) * exp(-(t-t0)/beta)\n                Input parameters are: <t0> <Cmax> <alpha> <beta>"


void gamma_nice_to_nasty(double *gparm);
void gamma_nasty_to_nice(gsl_vector * v, double *gparm);

double gamma_f(const gsl_vector * vary, void *params);
void gamma_df(const gsl_vector * vary, void *params, gsl_vector * df);
void gamma_fdf(const gsl_vector * vary, void *params, double *f, gsl_vector * df);

#endif /* GAMMA_FUNC_H */
