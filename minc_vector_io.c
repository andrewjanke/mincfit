/* minc_vector_io.c */

#include "minc_vector_io.h"
#include <volume_io.h>

#define MINC_VECTOR_CHUNK_SIZE 20

static const STRING MINC_VECTOR_FILE_HEADER = "MINC Vector File";
static const STRING MINC_VECTOR_TYPE = "MINC_Vector_Type";
static const STRING NORMAL_MINC_VECTOR = "Normal_MINC_Vector";
static const STRING MINC_VECTOR = "MINC_Vector";

/* create a MINC Vector */
MINC_Vector *new_MINC_Vector(int size)
{
   MINC_Vector *tmp;

   /* allocate space for the struct */
   tmp = (MINC_Vector *) malloc(sizeof(MINC_Vector));

   /* allocate sapce for the vector */
   SET_ARRAY_SIZE(tmp->V, 0, size, MINC_VECTOR_CHUNK_SIZE);
   tmp->size = size;

   return tmp;
   }

/* free a MINC Vector */
int free_MINC_Vector(MINC_Vector * mv)
{

   FREE(mv->V);
   FREE(mv);

   return (OK);
   }

/* input a MINC Vector from a file */
int input_MINC_Vector(char *input_file, MINC_Vector * mv)
{
   STRING   line;
   STRING   type_name;
   STRING   str;
   FILE    *file;

   /* parameter checking */
   if(input_file == NULL){
      print_error("input_MINC_Vector(): passed NULL FILE.\n");
      return (ERROR);
      }

   file = fopen(input_file, "r");
   if(file == NULL){
      print_error("input_MINC_Vector(): error opening %s\n", input_file);
      return (ERROR);
      }

   /* read the header */
   if(mni_input_string(file, &line, (char)0, (char)0) != OK){
      delete_string(line);
      print_error("input_MINC_Vector(): could not read header in file.\n");
      return (ERROR);
      }

   if(!equal_strings(line, MINC_VECTOR_FILE_HEADER)){
      delete_string(line);
      print_error("input_MINC_Vector(): invalid header in file.\n");
      return (ERROR);
      }

   /* read the type of the MINC_VECTOR */
   if(mni_input_keyword_and_equal_sign(file, MINC_VECTOR_TYPE, FALSE) != OK)
      return (ERROR);

   if(mni_input_string(file, &type_name, (char)';', (char)0) != OK){
      print_error("input_MINC_Vector(): missing MINC_VECTOR type.\n");
      return (ERROR);
      }

   if(mni_skip_expected_character(file, (char)';') != OK)
      return (ERROR);

   /* read the beginning of the vector */
   if(!equal_strings(type_name, NORMAL_MINC_VECTOR)){
      print_error("input_MINC_Vector(): invalid MINC_VECTOR type.\n");
      delete_string(type_name);
      return (ERROR);
      }
   delete_string(type_name);

   if(mni_input_string(file, &str, (char)'=', (char)0) != OK)
      return (ERROR);

   if(!equal_strings(str, MINC_VECTOR)){
      print_error("Expected %s =\n", MINC_VECTOR);
      delete_string(str);
      return (ERROR);
      }
   delete_string(str);

   if(mni_skip_expected_character(file, (char)'=') != OK)
      return (ERROR);

   /* input the data */
   return (mni_input_reals(file, &mv->size, &mv->V));
   }

/* output a MINC Vector to a file */
int output_MINC_Vector(char *output_file, MINC_Vector * mv, char **comments,
                          int clobber)
{
   int      i;
   FILE    *file;
   Status   status;

   /* parameter checking */
   if(output_file == NULL){
      print_error("output_MINC_Vector(): passed NULL FILE.\n");
      return (ERROR);
      }

   if(file_exists(output_file) && !clobber){
      print_error("output_MINC_Vector(): output file exists and NOCLOBBER.\n");
      return (ERROR);
      }

   file = fopen(output_file, "w");
   if(file == NULL){
      print_error("output_MINC_Vector(): error opening %s\n", output_file);
      return (ERROR);
      }

   /* output the header */
   status = output_string(file, MINC_VECTOR_FILE_HEADER);
   status &= output_newline(file);

   /* write out comments */
   i = 0;
   while(comments[i] != NULL){
      status &= output_string(file, "% ");
      status &= output_string(file, comments[i]);
      status &= output_newline(file);
      i++;
      }

   status &= output_newline(file);
   status &= output_newline(file);

   /* output vector type */
   status &= output_string(file, MINC_VECTOR_TYPE);
   status &= output_string(file, " = ");
   status &= output_string(file, NORMAL_MINC_VECTOR);
   status &= output_string(file, ";");
   status &= output_newline(file);

   /* output the vector itself */
   status &= output_string(file, MINC_VECTOR);
   status &= output_string(file, " =");
   status &= output_newline(file);

   /* now read the lines of the MINC_VECTOR */
   for(i = 0; i < mv->size; i++){

      fprintf(file, "%g", mv->V[i]);
/*      status &= output_real(file, mv->V[i]); */

      if(i == mv->size - 1){
         status &= output_string(file, ";");
         }
      status &= output_newline(file);
      }

   fclose(file);
   return status;
   }

/* debugging function to print a MINC Vector */
void print_MINC_Vector(const char *name, MINC_Vector * mv)
{
   int      i;

   fprintf(stdout, "%s[%d]: ", name, mv->size);
   for(i = 0; i < mv->size; i++){
      fprintf(stdout, "%g ", mv->V[i]);
      }
   fprintf(stdout, "\n");
   }
