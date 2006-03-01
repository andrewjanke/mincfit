/* mincfit.c                                                                  */
/*                                                                            */
/* do various regressions on MINC files                                       */
/*                                                                            */
/* Andrew Janke - a.janke@gmail.com                                           */
/* Mark Griffin - mgriffin@bic.mni.mcgill.ca                                  */
/*                                                                            */
/* Copyright Andrew Janke and Mark Griffin, McConnell Brain Imaging Centre    */
/* Permission to use, copy, modify, and distribute this software and its      */
/* documentation for any purpose and without fee is hereby granted,           */
/* provided that the above copyright notice appear in all copies.  The        */
/* author and McGill University make no representations about the             */
/* suitability of this software for any purpose.  It is provided "as is"      */
/* without express or implied warranty.                                       */

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics_double.h>

#include <ParseArgv.h>
#include <time_stamp.h>
#include <voxel_loop.h>
#include <volume_io.h>
#include <minc_vector_io.h>

#include "common_func.h"
#include "gamma_func.h"
#include "t1_func.h"

#define  GAMMA  1
#define  T1     2

#ifndef  FALSE
#define  FALSE  0
#endif

#ifndef  TRUE
#define  TRUE   1
#endif

#define MAX_ITER 1000
#define FNAME_SIZE 10

const char gamma_fnames[GAMMA_NVARY][FNAME_SIZE] = {"t0", "Cmax", "alpha", "beta"};

/* for arguments */
double gamma_vary[GAMMA_NVARY + GAMMA_NCONST] = {0,0,0,0};

const char t1_fnames[T1_NVARY][FNAME_SIZE] = {"S0", "T1"};

/* for arguments */
double t1_vary[T1_NVARY+T1_NCONST] = {0,0,0};


/* typedefs */
typedef struct {
   /* input parameters */
   int      masking;
   double   mask_val;
   int      mask_idx;
   
   /* the parameters */
   gsl_vector *vary;
   
   /* minimiser functions */
   gsl_multimin_function_fdf func;
   const gsl_multimin_fdfminimizer_type *T;
   gsl_multimin_fdfminimizer *s;
   
   } Loop_Data;

/* function prototypes */
void     do_loop(void *caller_data, long num_voxels, int input_num_buffers,
                 int input_vector_length, double *input_data[],
                 int output_num_buffers, int output_vector_length,
                 double *output_data[], Loop_Info * loop_info);
void     print_version_info(void);

/* Argument variables and table */
static int verbose = FALSE;
static int clobber = FALSE;
static int max_buffer = 4 * 1024;
static char *mask_fname = NULL;
static char *x_values_fname = NULL;
int fit_type;
static double tol = 1.0;
static Loop_Data md = {
   FALSE, 1.0, 0,
   NULL,
   NULL, NULL, NULL
   };

static ArgvInfo argTable[] = {
   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "General options:"},
   {"-version", ARGV_FUNC, (char *)print_version_info, (char *)NULL,
    "print version info and exit."},
   {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
    "print out extra information."},
   {"-clobber", ARGV_CONSTANT, (char *)TRUE, (char *)&clobber,
    "clobber existing files."},
   {"-max_buffer", ARGV_INT, (char *)1, (char *)&max_buffer,
    "maximum size of buffers (in kbytes)"},
   {"-mask", ARGV_STRING, (char *)1, (char *)&mask_fname,
    "only fit the time course for the voxels within the specified mask"},
   {"-mask_val", ARGV_FLOAT, (char *)1, (char *)&(md.mask_val),
    "mask value to use"},
   {"-tolerance", ARGV_FLOAT, (char *)1, (char *)&tol,
    "fitting tolerance"},

   {"-x_values", ARGV_STRING, (char *)1, (char *)&x_values_fname,
    "the vector in each voxel is a set of points y = f(x) (eg. x is time)"},

   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "\nFitting functions"},
   {"-gamma", ARGV_FLOAT, (char *)GAMMA_NVARY + GAMMA_NCONST, (char *)&gamma_vary,
    GAMMA_COMMENT},
   {"-T1", ARGV_FLOAT, (char *)T1_NVARY + T1_NCONST, (char *)&t1_vary,
    T1_COMMENT},


   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "\nOutput Options:"},

   {NULL, ARGV_HELP, NULL, NULL, ""},
   {NULL, ARGV_END, NULL, NULL, NULL}
   };


/************************************************************/

int main(int argc, char *argv[])
{
   char   **infiles;
   char   **outfiles;
   char    *out_base;
   int      n_infiles, n_outfiles;
   char    *arg_string;
   Loop_Options *loop_opts;
   int      i, j;
   
   int      n_frames;
   double   *vary_ptr;
   char     *out_fnames_ptr;
   
   
   MINC_Vector *tmp;
   
   FUNC_REC *func_rec;
   
   
   /* Save time stamp and args */
   arg_string = time_stamp(argc, argv);

   /* get arguments */
   if(ParseArgv(&argc, argv, argTable, 0) || (argc < 2)){
      fprintf(stderr, "\nUsage: %s [options] <frame1.mnc> <frame2.mnc> ... <out_base>\n",
              argv[0]);
      fprintf(stderr, "       %s -help\n\n", argv[0]);
      exit(EXIT_FAILURE);
      }
   
   /* get infile names */
   n_infiles = argc - 2;
   infiles = (char **)malloc(sizeof(char *) * (n_infiles + 1));   /* + 1 for mask */
   for(i = 0; i < n_infiles; i++){
      infiles[i] = argv[i + 1];
      }
   n_frames = n_infiles;

   /* get output base */
   out_base = argv[argc-1];

   /* set up mask file */
   if(mask_fname != NULL){
      infiles[n_infiles] = mask_fname;
      md.masking = TRUE;
      md.mask_idx = n_infiles;
      n_infiles++;
      }

   /* check for the infile(s) */
   for(i = 0; i < n_infiles; i++){
      if(access(infiles[i], F_OK) != 0){
         fprintf(stderr, "%s: Couldn't find input file: %s\n\n", argv[0], infiles[i]);
         exit(EXIT_FAILURE);
         }
      }
   
   /* check and set up which type of fitting we are doing */
   if(gamma_vary[0] != 0 ||
      gamma_vary[1] != 0 ||
      gamma_vary[2] != 0 ||
      gamma_vary[3] != 0){
      fit_type = GAMMA;
      }
   else if(t1_vary[0] != 0 ||
      t1_vary[1] != 0 ||
      t1_vary[2] != 0){
      fit_type = T1;
      }
   else{
      fprintf(stderr, "\nThis shouldn't happen, find Toto and get out of Kansas\n\n");
      exit(EXIT_FAILURE);
      }
  
   
   switch (fit_type) {
   case GAMMA:
      
      /* set up the optimisation functions and pointers */
      md.func.f = &gamma_f;
      md.func.df = &gamma_df;
      md.func.fdf = &gamma_fdf;
      md.func.n = GAMMA_NVARY;
      vary_ptr = gamma_vary;
      out_fnames_ptr = gamma_fnames;
      
      /* convert the gamma parameters to a more optimisable form */
      gamma_nice_to_nasty(vary_ptr);
      break;
      
     case T1:
      
      md.func.f = &t1_f;
      md.func.df = &t1_df;
      md.func.fdf = &t1_fdf;
      md.func.n = T1_NVARY;
      vary_ptr = t1_vary;
      out_fnames_ptr = t1_fnames;
      break;
      
      }

   /* set up func_rec */
   func_rec = (FUNC_REC*)malloc(sizeof(FUNC_REC));
   func_rec->v = gsl_vector_alloc(n_frames);
   func_rec->v_fit = gsl_vector_alloc(n_frames);
   func_rec->df = (double*)malloc(md.func.n * sizeof(double));
   func_rec->constants = &vary_ptr[md.func.n];
   md.func.params = (void *)func_rec;
   
   md.T = gsl_multimin_fdfminimizer_conjugate_fr;
   md.s = gsl_multimin_fdfminimizer_alloc(md.T, md.func.n);
   
   /* read in the x_values vector */
   tmp = new_MINC_Vector(0);
   if(input_MINC_Vector(x_values_fname, tmp) == OK){

      /* check # of x_values == # input frames */
      if(tmp->size != n_frames){
         fprintf(stderr,
                 "%s: # of input x_values %s (%d) must equal the number of frames (%d) \n\n",
                 argv[0], x_values_fname, tmp->size, n_frames);
         exit(EXIT_FAILURE);
         }
      
      /* convert to a gsl_vector */
      func_rec->x_values = gsl_vector_alloc(tmp->size);
      for(i = 0; i < tmp->size; i++){
         gsl_vector_set(func_rec->x_values, i, tmp->V[i]);
         }
      }
   else {
      fprintf(stderr, "%s: Couldn't read in: %s\n\n", argv[0], x_values_fname);
      exit(EXIT_FAILURE);
      }
      
   // set up initial guesss from C/L args //
   md.vary = gsl_vector_alloc(md.func.n);
   for(i=md.func.n; i--;){
      gsl_vector_set(md.vary, i, vary_ptr[i]);
      }

   
     /* set up and check for output files */
   n_outfiles = md.func.n;
   outfiles = (char **)malloc(sizeof(char *) * n_outfiles);
   for(i = 0; i < n_outfiles; i++){
      outfiles[i] = (char *)malloc((strlen(out_base) + FNAME_SIZE + 5)
                                   * sizeof(char));
      sprintf(outfiles[i], "%s.%s.mnc", out_base, t1_fnames[i]);

      if(access(outfiles[i], F_OK) == 0 && !clobber){
         fprintf(stderr, "%s: %s exists, use -clobber to overwrite\n\n",
                 argv[0], outfiles[i]);
         exit(EXIT_FAILURE);
         }
      }
   
   /* set up voxel loop options */
   loop_opts = create_loop_options();
   set_loop_verbose(loop_opts, FALSE);
   set_loop_clobber(loop_opts, clobber);
   set_loop_buffer_size(loop_opts, (long)1024 * max_buffer);

   /* do the sampling */
   voxel_loop(n_infiles, infiles, n_outfiles, outfiles, arg_string,
              loop_opts, do_loop, (void *)&md);

   /* tidy up */
   free_loop_options(loop_opts);
   
   
   // lots of gsl_free
   gsl_multimin_fdfminimizer_free(md.s);

   gsl_vector_free(md.vary);


//   gsl_vector_free(f->v);
//   gsl_vector_free(f->v_fit);
//   free(f->df);
//   free(f);

   return (EXIT_SUCCESS);
   }



/* get points from file(s), write out to an input FP */
void do_loop(void *caller_data, long num_voxels, int input_num_buffers,
             int input_vector_length, double *input_data[], int output_num_buffers,
             int output_vector_length, double *output_data[], Loop_Info * loop_info)
{
   Loop_Data *md = (Loop_Data *) caller_data;
   FUNC_REC *func_rec = (FUNC_REC *) md->func.params;
   int      i, idim, ivox;
   int      n_infiles;
   int      do_sample;
   double   mask_value;
   int      dim_index;
   long     index[MAX_VAR_DIMS];
   int iter;
   int status;
   
   double tmp[4];

   /* shut the compiler up */
   // (void)output_num_buffers;
   // (void)output_vector_length;

   n_infiles = (md->masking) ? input_num_buffers - 1 : input_num_buffers;

   /* for each voxel */
   for(ivox = 0; ivox < num_voxels * input_vector_length; ivox++){
      /* nasty way that works for masking or not */
      
      mask_value = 0;
      if(!md->masking ||
         (md->masking && fabs(input_data[md->mask_idx][ivox] - md->mask_val) < 0.5)){

        // fprintf(stderr, "%d\n", ivox);
        
         /* set up the input vector */
         fprintf(stdout, "INPUT vector: ");
         for(i = 0; i < n_infiles; i++){
            fprintf(stdout, "%g ", input_data[i][ivox]);
            gsl_vector_set(func_rec->v, i, input_data[i][ivox]);
            }
         fprintf(stdout, "\n");
         
         
         gsl_multimin_fdfminimizer_set(md->s, &(md->func), md->vary, 0.1, 0.01 * tol);
   
         iter = 0;
         status = GSL_CONTINUE;
         while(status == GSL_CONTINUE && iter++ < MAX_ITER){
         
            fprintf(stderr, "ITER: %d\n", iter);
        //    fprintf(stderr, "start of iteration\n");
            
            if(gsl_multimin_fdfminimizer_iterate(md->s)){
               break;
               }
            
            status = gsl_multimin_test_gradient(md->s->gradient, 0.5 * tol);
            
        //    if(status == GSL_SUCCESS){
        //       printf ("iter: %5d  - %.5f %.5f %.5f %.5f - FIT: %10.5f\n", iter,
        //            gsl_vector_get (md->s->x, 0), 
        //            gsl_vector_get (md->s->x, 1), 
        //            gsl_vector_get (md->s->x, 2), 
        //            gsl_vector_get (md->s->x, 3), 
        //            md->s->f);
        //       }
            
            
            
            }
         
           switch (fit_type) {
           case GAMMA:
      
            /* convert the gamma parameters to a more optimisable form */
            gamma_nasty_to_nice(md->s->x, tmp);
         
            /* write out the resulting parameters */
            fprintf(stdout, "OUTPUT parameters: ");
            for(i = 0; i < output_num_buffers; i++){
             //  output_data[i][ivox] = gsl_vector_get(md->s->x, i);
               output_data[i][ivox] = tmp[i];
               fprintf(stdout, "%g ", tmp[i]);
               }
            fprintf(stdout, "\n");
            break;
            
            case T1:
            default:
            fprintf(stdout, "OUTPUT parameters: ");
               for(i = 0; i < output_num_buffers; i++){
                  output_data[i][ivox] = gsl_vector_get(md->s->x, i);
               fprintf(stdout, "%g ", output_data[i][ivox]);
                  }
            fprintf(stdout, "\n");
               break;
            }    

         }
      }

   }

void print_version_info(void)
{
   fprintf(stdout, "\n");
   fprintf(stdout, "%s version %s\n", PACKAGE, VERSION);
   fprintf(stdout, "Comments to %s\n", PACKAGE_BUGREPORT);
   fprintf(stdout, "\n");
   exit(EXIT_SUCCESS);
   }
