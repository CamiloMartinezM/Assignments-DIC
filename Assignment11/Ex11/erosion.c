#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                         EROSION WITH FLAT DISK                           */
/*                                                                          */
/*                 (Copyright by Joachim Weickert, 1/2022)                  */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/* 
  Rouy-Tourin scheme (explicit, one-sided differences). 
  Satisfies extremum principle for time steps ht <= 0.707.
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

void read_long

     (long *v)         /* value to be read */

/*
  reads a long value v
*/

{
char   row[80];    /* string for reading data */

if (fgets (row, 80, stdin) == NULL)
{
   printf("could not read string, aborting\n");
   exit(1);
}

if (row[strlen(row)-1] == '\n')
   row[strlen(row)-1] = 0;
sscanf(row, "%ld", &*v);

return;

}  /* read_long */

/*--------------------------------------------------------------------------*/

void read_double

     (double *v)         /* value to be read */

/*
  reads a double value v
*/

{
char   row[80];    /* string for reading data */

if (fgets (row, 80, stdin) == NULL)
{
   printf("could not read string, aborting\n");
   exit(1);
}

if (row[strlen(row)-1] == '\n')
   row[strlen(row)-1] = 0;
sscanf(row, "%lf", &*v);

return;

}  /* read_double */

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

void dummies_double

     (double **u,        /* image */
      long   nx,         /* size in x direction */
      long   ny)         /* size in y direction */

/*
  creates dummy boundaries for a double format image u by mirroring
*/

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
    {
    u[i][0]    = u[i][1];
    u[i][ny+1] = u[i][ny];
    }

for (j=0; j<=ny+1; j++)
    {
    u[0][j]    = u[1][j];
    u[nx+1][j] = u[nx][j];
    }

return;

}  /* dummies_double */

/*--------------------------------------------------------------------------*/

void analyse_grey_double

     (double  **u,         /* image, unchanged */
      long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in y direction */
      double  *min,        /* minimum, output */
      double  *max,        /* maximum, output */
      double  *mean,       /* mean, output */
      double  *std)        /* standard deviation, output */

/*
  computes minimum, maximum, mean, and standard deviation of a greyscale
  image u in double format
*/

{
long    i, j;       /* loop variables */
double  help1;      /* auxiliary variable */
double  help2;      /* auxiliary variable */

/* compute maximum, minimum, and mean */
*min  = u[1][1];
*max  = u[1][1];
help1 = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     if (u[i][j] < *min) *min = u[i][j];
     if (u[i][j] > *max) *max = u[i][j];
     help1 = help1 + u[i][j];
     }
*mean = help1 / (nx * ny);

/* compute standard deviation */
*std = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     help2  = u[i][j] - *mean;
     *std = *std + help2 * help2;
     }
*std = sqrt (*std / (nx * ny));

return;

}  /* analyse_grey_double */

/*--------------------------------------------------------------------------*/

double max3

       (double  a,           /* first argument */
        double  b,           /* second argument */
        double  c)           /* third argument */

/*
  returns the maximum of a, b, and c
*/

{
double result;

result = a;
if (b > result) result = b;
if (c > result) result = c;

return (result);

}  /* max3 */

/*--------------------------------------------------------------------------*/

double sqr 

       (double  a)           /* argument */

/* 
  returns the square of a 
*/

{
double result;

result = a * a;

return (result);

}  /* sqr */

/*--------------------------------------------------------------------------*/

void erosion 

     (double   tau,       /* time step size, 0 < tau <= 1/sqrt(2) */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      double   hx,        /* pixel size in x direction */
      double   hy,        /* pixel size in y direction */
      double   **u)       /* input: original image; output: processed */

/*
  one erosion step with Rouy-Tourin scheme and time step size tau 
*/

{
long    i, j;                  /* loop variables */
double  **f;                   /* u at old time level */
double  fxp, fxm, fyp, fym;    /* one-sided finite differences */


/* ---- allocate memory ---- */

alloc_double_matrix (&f, nx+2, ny+2);


/* ---- copy u into f ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     f[i][j] = u[i][j];


/* ---- create reflecting dummy boundaries for f ---- */

dummies_double (f, nx, ny);


/* ---- loop ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* one-sided spatial derivatives */
     fxp = (f[i+1][j] - f[i][j]) / hx;
     fxm = (f[i][j] - f[i-1][j]) / hx;
     fyp = (f[i][j+1] - f[i][j]) / hy;
     fym = (f[i][j] - f[i][j-1]) / hy;

     /* update step */
/*
     SUPPLEMENT YOUR CODE HERE
*/
     u[i][j] = f[i][j] - tau * sqrt ( sqr (max3 (0.0, fxm, -fxp))
                                    + sqr (max3 (0.0, fym, -fyp)) );

     }


/* ---- free memory ---- */

free_double_matrix (f, nx+2, ny+2);

return;

}  /* erosion */

/*--------------------------------------------------------------------------*/

int main ()

{
char    in[80];               /* for reading data */
char    out[80];              /* for reading data */
double  **u;                  /* image */
long    nx, ny;               /* image size in x, y direction */
long    k;                    /* iteration counter */
long    kmax;                 /* number of iterations */
double  tau;                  /* time step size */
double  max, min;             /* largest, smallest grey value */
double  mean;                 /* average grey value */
double  std;                  /* standard deviation */
char    comments[1600];       /* string for comments */


printf ("\n");
printf ("ROUY-TOURIN SCHEME FOR EROSION WITH FLAT DISK\n\n");
printf ("**********************************************\n\n");
printf ("    Copyright 2022 by Joachim Weickert        \n");
printf ("    Dept. of Mathematics and Computer Science \n");
printf ("    Saarland University, Saarbruecken, Germany\n\n");
printf ("    All rights reserved. Unauthorized usage,  \n");
printf ("    copying, hiring, and selling prohibited.  \n\n");
printf ("    Send bug reports to                       \n");
printf ("    weickert@mia.uni-saarland.de              \n\n");
printf ("**********************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf ("input image (pgm):                ");
read_string (in);
read_pgm_to_double (in, &nx, &ny, &u);


/* ---- read parameters ---- */

printf ("time step size (<=0.707):         ");
read_double (&tau);

printf ("number of iterations:             ");
read_long (&kmax);

printf ("output image (pgm):               ");
read_string (out);

printf ("\n");


/* ---- analyse initial image ---- */

analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);
printf ("initial image\n");
printf ("minimum:           %8.2lf \n", min);
printf ("maximum:           %8.2lf \n", max);
printf ("mean:              %8.2lf \n", mean);
printf ("standard dev.:     %8.2lf \n\n", std);


/* ---- process image ---- */

for (k=1; k<=kmax; k++)
    {
    /* perform one iteration */
    erosion (tau, nx, ny, 1.0, 1.0, u);
    
    /* analyse filtered image */
    analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);
    printf ("iteration number:  %8ld \n", k);
    printf ("minimum:           %8.2f \n", min);
    printf ("maximum:           %8.2f \n", max);
    printf ("mean:              %8.2f \n", mean);
    printf ("standard dev.:     %8.2f \n\n", std);
    }


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
comment_line (comments, "# erosion with flat disk, Rouy-Tourin scheme\n");
comment_line (comments, "# initial image:  %s\n", in);
comment_line (comments, "# tau:            %8.4lf\n", tau);
comment_line (comments, "# iterations:     %8ld\n", kmax);
comment_line (comments, "# min:            %8.2lf\n", min);
comment_line (comments, "# max:            %8.2lf\n", max);
comment_line (comments, "# mean:           %8.2lf\n", mean);
comment_line (comments, "# standard dev.:  %8.2lf\n", std);

/* write image */
write_double_to_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

free_double_matrix (u, nx+2, ny+2);

return(0);
}
