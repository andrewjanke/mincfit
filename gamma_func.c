/* 
 *  gamma_func.c
 *
 *	This file provides the code vector_fit_function and minc_fit_function with all of 
 *	the details about fitting a gamma variate function to a timecourse.
 *
 * C  =  Cmax  x  ((t - t0) ^ alpha)  x  exp(-(t-t0) / beta)
 *
 *	See equation 5 of Harrer, 2004 for details.
 */

#include "common_func.h"
#include "gamma_func.h"
#include "math.h"

/*
 *  gamma_f -
 *
 *	create a vector of the function being fit (v_fit) and
 *	return a measure as to how well this function matches
 *	the original data (v) as a sum of the squares.
 */

#define REALLY_BIG_NUMBER 1.0e16
double gamma_f(const gsl_vector * vary, void *params)
{
   FUNC_REC *func_rec;
   gsl_vector *v, *v_fit;
   double   Cmax, t0, alpha, beta, t, sum;
   int      i;
   double   gparm[GAMMA_NVARY];

   func_rec = (FUNC_REC *) params;
   v = func_rec->v;
   v_fit = func_rec->v_fit;

   /* convert from nasty to nice */
   gamma_nasty_to_nice((gsl_vector *) vary, gparm);

   t0 = gparm[0];
   Cmax = gparm[1];
   alpha = gparm[2];
   beta = gparm[3];

   if(alpha < 0 || beta < 0 || Cmax < 0){
      return REALLY_BIG_NUMBER;
      }

   sum = 0.0;
   for(i = 0; i < v->size; i++){
      t = gsl_vector_get(func_rec->x_values, i) - t0;

      if(t > 0){
         gsl_vector_set(v_fit, i, Cmax * pow(t, alpha) * exp(-t / beta));
         }
      else {
         gsl_vector_set(v_fit, i, 0);
         }

      sum += pow(gsl_vector_get(v, i) - gsl_vector_get(v_fit, i), 2.0);
      }

   return (sum);
   }

/************************************************************/

/*
 *  gamma_df -
 *
 *	return a vector of the derivative of f (the sum of the
 *	squares) with respect to each optimisation variable.
 */
void gamma_df(const gsl_vector * vary, void *params, gsl_vector * df)
{
   FUNC_REC *func_rec;
   double   t0, Cmax, alpha, beta;
   double   dt0, dCmax, dalpha, dbeta;
   double   dt02, dCmax2, dalpha2, dbeta2;

   double   t, fit_i, v_i;
   int      i;

   func_rec = (FUNC_REC *) params;

   double   local_df[GAMMA_NVARY];
   double   gparm[GAMMA_NVARY];

   /* init the input df vector */
   gsl_vector_set_zero(df);

   /* convert from nasty to nice */
   gamma_nasty_to_nice((gsl_vector *) vary, gparm);
   t0 = gparm[0];
   Cmax = gparm[1];
   alpha = gparm[2];
   beta = gparm[3];

   /* init variables */
   dt0 = dCmax = dalpha = dbeta = 0.0;

   for(i = 0; i < func_rec->v->size; i++){
      /* retrieve the value of fit_i (set previously in gamma_f) */
      fit_i = gsl_vector_get(func_rec->v_fit, i);
      v_i = gsl_vector_get(func_rec->v, i);
      t = gsl_vector_get(func_rec->x_values, i) - t0;

      if(t > 0){
         dCmax += 2 * (fit_i - v_i) * fit_i / Cmax;
         if(t != 0.0){
            dt0 += 2 * (fit_i - v_i) * fit_i * ((1 / beta) - (alpha / t));
            dalpha += 2 * (fit_i - v_i) * fit_i * log(t);
            dbeta += 2 * (fit_i - v_i) * fit_i * t / (beta * beta);
            }
         }
      }
   gsl_vector_set(df, 0, dt0 + (dCmax * Cmax / beta) - (dbeta / alpha));

   dalpha2 = 1.0 / (gsl_vector_get(vary, 1) * (1.0 - log(2.0)));
   dbeta2 = -dalpha2 * beta / alpha;
   dCmax2 = pow(alpha * beta, -alpha) * exp(alpha);
   dCmax2 = dCmax2 * (1.0 + gsl_vector_get(vary, 1) *
                      (-log(alpha * beta) * dalpha2 - (alpha / beta) * dbeta2));
   gsl_vector_set(df, 1, dCmax * dCmax2 + dalpha * dalpha2 + dbeta * dbeta2);

   dCmax2 = pow(alpha * beta, -alpha) * exp(alpha);
   dCmax2 = -dCmax2 * gsl_vector_get(vary, 1) / beta;
   gsl_vector_set(df, 2, dCmax * dCmax2 + dbeta / alpha);

   dalpha2 = 1.0 / (gsl_vector_get(vary, 3) * (log(2.0) - 1.0));
   dbeta2 = -dalpha2 * beta / alpha;
   dCmax2 = pow(alpha * beta, -alpha) * exp(alpha);
   dCmax2 = dCmax2 * gsl_vector_get(vary, 1) *
      (-log(alpha * beta) * dalpha2 - (alpha / beta) * dbeta2);
   gsl_vector_set(df, 3, dCmax * dCmax2 + dalpha * dalpha2 + dbeta * dbeta2);
   }

/************************************************************/

void gamma_fdf(const gsl_vector * x, void *params, double *f, gsl_vector * df)
{
   *f = gamma_f(x, params);
   gamma_df(x, params, df);
   }

/* convert the gamma function parameters to a form in which */
/*   the optimisation function space is smoother            */
void gamma_nice_to_nasty(double *gparm)
{
   double   t0, Cmax, alpha, beta;

   /* get the original nice values */
   t0 = gparm[0];
   Cmax = gparm[1];
   alpha = gparm[2];
   beta = gparm[3];

   /* overwrite with the nasty ones */
   if((alpha > 0.0) && (beta > 0.0) && (Cmax > 0.0)){
      gparm[1] = Cmax * pow(alpha * beta, alpha) * exp(-alpha);
      gparm[2] = alpha * beta + t0;
      gparm[3] = exp(alpha * (log(2.0) - 1.0) + log(gparm[1]));
      }
   else {
      gparm[1] = 1.0;
      gparm[2] = 1.0;
      gparm[3] = 1.0;
      }
   }

/* convert the gamma function parameters back */
void gamma_nasty_to_nice(gsl_vector * v, double *gparm)
{
   double   t0, Cmax, alpha, beta;

   t0 = gsl_vector_get(v, 0);

   if((gsl_vector_get(v, 3) > 0.0) && (gsl_vector_get(v, 1) > 0.0)){
      alpha = (log(gsl_vector_get(v, 3)) - log(gsl_vector_get(v, 1)))
         / (log(2.0) - 1.0);

      beta = (gsl_vector_get(v, 2) - t0) / alpha;

      if((alpha > 0.0) && (beta > 0.0)){
         Cmax = gsl_vector_get(v, 1) / (pow(alpha * beta, alpha)
                                        * exp(-alpha));
         }
      else {
         Cmax = 1.0;
         }
      }
   else {
      Cmax = 1.0;
      alpha = 1.0;
      beta = 1.0;
      }

   gparm[0] = t0;
   gparm[1] = Cmax;
   gparm[2] = alpha;
   gparm[3] = beta;
   }
