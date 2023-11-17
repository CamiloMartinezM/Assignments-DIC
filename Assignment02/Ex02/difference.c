#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                   ABSOLUTE DIFFERENCE OF TWO IMAGES                      */
/*                                                                          */
/*                 (Copyright Joachim Weickert, 11/2020)                    */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
  Reads two pgm (P5) greyscale images, computes the absolute differences
  in each pixel, as well as their mean and their maximum 
*/

/*--------------------------------------------------------------------------*/

void alloc_double_matrix

     (double ***matrix,  /* matrix */
      long   n1,         /* size in direction 1 */
      long   n2)         /* size in direction 2 */

/*
  allocates memory for a double format matrix of size n1 * n2
*/

{
long i;    /* loop variable */

*matrix = (double **) malloc (n1 * sizeof(double *));

if (*matrix == NULL)
   {
   printf("alloc_double_matrix: not enough memory available\n");
   exit(1);
   }

for (i=0; i<n1; i++)
    {
    (*matrix)[i] = (double *) malloc (n2 * sizeof(double));
    if ((*matrix)[i] == NULL)
       {
       printf("alloc_double_matrix: not enough memory available\n");
       exit(1);
       }
    }

return;

}  /* alloc_double_matrix */

/*--------------------------------------------------------------------------*/

void free_double_matrix

     (double  **matrix,   /* matrix */
      long    n1,         /* size in direction 1 */
      long    n2)         /* size in direction 2 */

/*
  frees memory for a double format matrix of size n1 * n2
*/

{
long i;   /* loop variable */

for (i=0; i<n1; i++)
free(matrix[i]);

free(matrix);

return;

}  /* free_double_matrix */

/*--------------------------------------------------------------------------*/

void read_string

     (char *v)         /* string to be read */

/*
  reads a string v
*/

{
if (fgets (v, 80, stdin) == NULL)
{
   printf("could not read string, aborting\n");
   exit(1);
}

if (v[strlen(v)-1] == '\n')
   v[strlen(v)-1] = 0;

return;

}  /* read_string */

/*--------------------------------------------------------------------------*/

void skip_white_space_and_comments

     (FILE *inimage)  /* input file */

/*
  skips over white space and comments while reading the file
*/

{

int   ch = 0;   /* holds a character */
char  row[80];  /* for reading data */

/* skip spaces */
while (((ch = fgetc(inimage)) != EOF) && isspace(ch));

/* skip comments */
if (ch == '#')
   {
   if (fgets(row, sizeof(row), inimage))
      skip_white_space_and_comments (inimage);
   else
      {
      printf("skip_white_space_and_comments: cannot read file\n");
      exit(1);
      }
   }
else
   fseek (inimage, -1, SEEK_CUR);

return;

} /* skip_white_space_and_comments */

/*--------------------------------------------------------------------------*/

void read_pgm_to_double
   
     (const char  *file_name,    /* name of pgm file */
      long        *nx,           /* image size in x direction, output */
      long        *ny,           /* image size in y direction, output */
      double      ***u)          /* image, output */

/* 
  reads a greyscale image that has been encoded in pgm format P5 to
  an image u in double format;
  allocates memory for the image u;
  adds boundary layers of size 1 such that
  - the relevant image pixels in x direction use the indices 1,...,nx
  - the relevant image pixels in y direction use the indices 1,...,ny
*/

{  
char  row[80];      /* for reading data */
long  i, j;         /* image indices */
long  max_value;    /* maximum color value */
FILE  *inimage;     /* input file */

/* open file */
inimage = fopen (file_name, "rb");
if (inimage == NULL)
   {
   printf ("read_pgm_to_double: cannot open file '%s'\n", file_name);
   exit(1);
   }

/* read header */
if (fgets(row, 80, inimage) == NULL)
   {
   printf ("read_pgm_to_double: cannot read file\n");
   exit(1);
   }

/* image type: P5 */
if ((row[0] == 'P') && (row[1] == '5'))
   {
   /* P5: grey scale image */
   }
else
   {
   printf ("read_pgm_to_double: unknown image format\n");
   exit(1);
   }

/* read image size in x direction */
skip_white_space_and_comments (inimage);
if (!fscanf (inimage, "%ld", nx))
   {
   printf ("read_pgm_to_double: cannot read image size nx\n");
   exit(1);
   }

/* read image size in x direction */
skip_white_space_and_comments (inimage);
if (!fscanf (inimage, "%ld", ny))
   {
   printf ("read_pgm_to_double: cannot read image size ny\n");
   exit(1);
   }

/* read maximum grey value */
skip_white_space_and_comments (inimage);
if (!fscanf (inimage, "%ld", &max_value))
   {
   printf ("read_pgm_to_double: cannot read maximal value\n");
   exit(1);
   }
fgetc(inimage);

/* allocate memory */
alloc_double_matrix (u, (*nx)+2, (*ny)+2);

/* read image data row by row */
for (j=1; j<=(*ny); j++)
 for (i=1; i<=(*nx); i++)
     (*u)[i][j] = (double) getc(inimage);

/* close file */
fclose (inimage);

return;

}  /* read_pgm_to_double */

/*--------------------------------------------------------------------------*/

void comment_line

     (char* comment,       /* comment string (output) */
      char* lineformat,    /* format string for comment line */
      ...)                 /* optional arguments */

/*
  Adds a line to the comment string comment. The string line can contain
  plain text and format characters that are compatible with sprintf.
  Example call:
  print_comment_line(comment, "Text %lf %ld", double_var, long_var).
  If no line break is supplied at the end of the input string, it is
  added automatically.
*/

{
char     line[80];
va_list  arguments;

/* get list of optional function arguments */
va_start (arguments, lineformat);

/* convert format string and arguments to plain text line string */
vsprintf (line, lineformat, arguments);

/* add line to total commentary string */
strncat (comment, line, 80);

/* add line break if input string does not end with one */
if (line[strlen(line)-1] != '\n')
   strncat (comment, "\n", 1);

/* close argument list */
va_end (arguments);

return;

}  /* comment_line */

/*--------------------------------------------------------------------------*/

void write_double_to_pgm

     (double  **u,          /* image, unchanged */
      long    nx,           /* image size in x direction */
      long    ny,           /* image size in y direction */
      char    *file_name,   /* name of pgm file */
      char    *comments)    /* comment string (set 0 for no comments) */

/*
  writes a greyscale image in double format into a pgm P5 file
*/

{
FILE           *outimage;  /* output file */
long           i, j;       /* loop variables */
double         aux;        /* auxiliary variable */
unsigned char  byte;       /* for data conversion */

/* open file */
outimage = fopen (file_name, "wb");
if (NULL == outimage)
   {
   printf("could not open file '%s' for writing, aborting\n", file_name);
   exit(1);
   }

/* write header */
fprintf (outimage, "P5\n");                  /* format */
if (comments != 0)
   fputs (comments, outimage);               /* comments */
fprintf (outimage, "%ld %ld\n", nx, ny);     /* image size */
fprintf (outimage, "255\n");                 /* maximal value */

/* write image data */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     {
     aux = u[i][j] + 0.499999;    /* for correct rounding */
     if (aux < 0.0)
        byte = (unsigned char)(0.0);
     else if (aux > 255.0)
        byte = (unsigned char)(255.0);
     else
        byte = (unsigned char)(aux);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }

/* close file */
fclose (outimage);

return;

}  /* write_double_to_pgm */

/*--------------------------------------------------------------------------*/

int main ()

{
char    in1[80];              /* for reading data */
char    in2[80];              /* for reading data */
char    out[80];              /* for reading data */
double  **u;                  /* greyscale image 1 */
double  **v;                  /* greyscale image 2 */
double  **w;                  /* absolute difference image */
long    i, j;                 /* loop variables */
long    nx, ny;               /* image size in x, y direction */
double  max;                  /* largest absolute difference */
double  mean;                 /* average absolute difference */
char    comments[1600];       /* string for comments */


printf ("\n");
printf ("ABSOLUTE DIFFERENCE BETWEEN TWO GREYSCALE (PGM) IMAGES\n\n");
printf ("**************************************************\n\n");
printf ("    Copyright 2020 by Joachim Weickert            \n");
printf ("    Dept. of Mathematics and Computer Science     \n");
printf ("    Saarland University, Saarbruecken, Germany    \n\n");
printf ("    All rights reserved. Unauthorised usage,      \n");
printf ("    copying, hiring, and selling prohibited.      \n\n");
printf ("    Send bug reports to                           \n");
printf ("    weickert@mia.uni-saarland.de                  \n\n");
printf ("**************************************************\n\n");


/* ---- read input images and output image (pgm format P5) ---- */

printf ("input image 1 (pgm):             ");
read_string (in1);
read_pgm_to_double (in1, &nx, &ny, &u);

printf ("input image 2 (pgm):             ");
read_string (in2);
read_pgm_to_double (in2, &nx, &ny, &v);

printf ("output image (pgm):              ");
read_string (out);


/* ---- allocate memory ---- */

alloc_double_matrix (&w, nx+2, ny+2);


/* ---- compute the absolute differences ---- */

/* compute absolute difference */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     w[i][j] = abs (u[i][j] - v[i][j]);

/* check maximum and average absolute difference */
max  = w[1][1];
mean = 0.0;
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     {
     if (w[i][j] > max) 
        max = w[i][j];
     mean = mean + w[i][j];
     }
mean = mean / (nx * ny);
printf ("\n");
printf ("average absolute difference:     %6.4lf \n",   mean);
printf ("maximum absolute difference:     %1.0lf \n\n", max);


/* ---- scale difference image w to [0, 255] ---- */

if (max > 0.0)
   for (j=1; j<=ny; j++)
    for (i=1; i<=nx; i++)
        w[i][j] = w[i][j] * 255.0 / max;


/* ---- write output image (pgm format P5) ---- */

/* write parameter values in comment string */
comments[0] = '\0';
comment_line (comments, "# absolute differences between 2 images\n");
comment_line (comments, "# scaled to [0, 255] for visualisation\n");
comment_line (comments, "# image 1:   %s\n", in1);
comment_line (comments, "# image 2:   %s\n", in2);
comment_line (comments, "# aver. abs. difference:  %6.4lf\n", mean);
comment_line (comments, "# max.  abs. difference:  %1.0lf\n", max);

/* write image data */
write_double_to_pgm (w, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory ---- */

free_double_matrix (u, nx+2, ny+2);
free_double_matrix (v, nx+2, ny+2);
free_double_matrix (w, nx+2, ny+2);

return(0);

}  /* main */
