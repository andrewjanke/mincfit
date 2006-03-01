/* minc_vector_io.h */

#ifndef MINC_VECTOR_IO_H
#define MINC_VECTOR_IO_H



/* Structure for MINC Vector information */
typedef struct {
   int      size;
   double  *V;
   } MINC_Vector;

/* create a MINC Vector (does allocation) */
MINC_Vector *new_MINC_Vector(int size);

/* free a MINC Vector */
int   free_MINC_Vector(MINC_Vector * mv);

/* input a MINC Vector from a file */
int   input_MINC_Vector(char *input_file, MINC_Vector * mv);

/* output a MINC Vector to a file */
int   output_MINC_Vector(char *output_file, MINC_Vector * mv, char **comments,
                            int clobber);

/* debugging function to print a MINC Vector */
void     print_MINC_Vector(const char *name, MINC_Vector * mv);

#endif
