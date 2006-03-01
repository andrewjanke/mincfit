 /* 
  *  t1_func.h
  *
  *   This file provides the code fitting a t1 signal intensity function to 
  * spoiled gradient echo data collected with different flip angles.
  *
  *   See Equation 10.10, p. 351 of "Quantitative MRI of the Brain", Tofts, for details.
  */

/*
 *	T1_NVARY is the number of parameters which are optimised in the function fitting.
 */

#ifndef T1_FUNC_H
#define T1_FUNC_H

#define T1_NVARY       2
#define T1_NCONST      1
#define T1_COMMENT     "use the T1 signal intensity function\n                f = (S0 * sin(theta) * (1 - exp(-TR/T1))) / (1 - cos(theta) * exp(-TR/T1))\n                the x values (flip angles) are in degrees\n                TR and T1 are in milliseconds (only S0 and T1 are optimised)\n                Input parameters are: <S0> <T1> <TR>"

double   t1_f(const gsl_vector * vary, void *params);
void     t1_df(const gsl_vector * vary, void *params, gsl_vector * df);
void     t1_fdf(const gsl_vector * vary, void *params, double *f, gsl_vector * df);

#endif /* T1_FUNC_H */
