/* 
 * t1_func.c
 *
 *	This file provides the code fitting a t1 signal intensity function to 
 * spoiled gradient echo data collected with different flip angles.
 *
 *
 *       S0 x sin(theta) x (1 - E1)
 * S  =  --------------------------
 *          1 - cos(theta) * E1
 *
 *
 * where E1 = exp(-TR/T1), 
 * see Equation 10.10, p. 351 of "Quantitative MRI of the Brain", Tofts.
 */

#include "common_func.h"
#include "t1_func.h"
#include <math.h>

/*
 * t1_f -
 *
 *	create a vector of the function being fit (v_fit) and
 *	return a measure as to how well this function matches
 *	the original data (v) as a sum of the squares.
 */
double t1_f(const gsl_vector * vary, void *params)
{
   FUNC_REC *func_rec;
   gsl_vector *v, *v_fit;
   double   S0, T1, TR, theta, sum;
   int      i;

   func_rec = (FUNC_REC *) params;
   v = func_rec->v;
   v_fit = func_rec->v_fit;
   
   S0 = gsl_vector_get(vary, 0);
   T1 = gsl_vector_get(vary, 1);
   TR = func_rec->constants[0];
   
   sum = 0.0;
   for(i = 0; i < v->size; i++){
      theta = gsl_vector_get(func_rec->x_values, i) * M_PI / 180; 
      gsl_vector_set(v_fit, i, (S0 * sin(theta) * (1 - exp(-TR/T1))) / (1 - cos(theta) * exp(-TR/T1)));
      
      sum += pow(gsl_vector_get(v, i) - gsl_vector_get(v_fit, i), 2.0);
      }

   return (sum);
   }

/************************************************************/

/*
 *  t1_df -
 *
 *	return a vector of the derivative of f (the sum of the
 *	squares) with respect to each optimisation variable.
 */
void t1_df(const gsl_vector * vary, void *params, gsl_vector * df)
{
   FUNC_REC *func_rec;
   double   S0, T1, TR, theta;
   double   fac;
   int      i;
   
   func_rec = (FUNC_REC *) params;
   
   double local_df[T1_NVARY];
   
   S0 = gsl_vector_get(vary, 0);
   T1 = gsl_vector_get(vary, 1);
   TR = func_rec->constants[0];
   
   local_df[0] = 0;
   local_df[1] = 0;

   for(i = 0; i < func_rec->v->size; i++){
      theta = gsl_vector_get(func_rec->x_values, i) * M_PI / 180;
   
      gsl_vector_set(func_rec->v_fit, i, (S0 * sin(theta) * (1 - exp(-TR/T1))) / (1 - cos(theta) * exp(-TR/T1)));
      
      fac = 2.0 * (gsl_vector_get(func_rec->v, i) - gsl_vector_get(func_rec->v_fit, i));
      local_df[0] -= fac * (sin(theta) * (1 - exp(-TR/T1))) / (1 - cos(theta) * exp(-TR/T1));
      local_df[1] -= fac * (S0 * sin(theta) * (cos(theta) - 1) * (TR/(T1*T1)) * exp(-TR/T1)) / pow(1 - cos(theta) * exp(-TR/T1), 2);
         
      }

   gsl_vector_set(df, 0, local_df[0]);
   gsl_vector_set(df, 1, local_df[1]);
   }

/************************************************************/

void t1_fdf(const gsl_vector * x, void *params, double *f, gsl_vector * df)
{
   *f = t1_f(x, params);
   t1_df(x, params, df);
   }
