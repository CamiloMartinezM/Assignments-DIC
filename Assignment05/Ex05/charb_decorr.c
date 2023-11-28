#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*  PERONA-MALIK FILTER WITH CHARBONNIER DIFFUSIVITY FOR GREYSCALE IMAGES   */
/*                                                                          */
/*                 (Copyright by Joachim Weickert, 12/2020)                 */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/* 
  explicit finite difference scheme;
  stopped with decorrelation criterion of Mrazek and Navara (IJCV 2003)
*/

/*--------------------------------------------------------------------------*/

void alloc_double_vector

     (double **vector,   /* vector */
      long   n1)         /* size */

/* 
  allocates memory for a double format vector of size n1
*/

{
*vector = (double *) malloc (n1 * sizeof(double));

if (*vector == NULL)
   {
   printf("alloc_double_vector: not enough memory available\n");
   exit(1);
   }

return;

}  /* alloc_double_vector */

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

void free_double_vector

     (double  *vector,    /* vector */
      long    n1)         /* size */

/* 
  frees memory for a double format vector of size n1
*/

{

free(vector);
return;

}  /* free_double_vector */

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
*std = sqrt(*std / (nx * ny));

return;

}  /* analyse_grey_double */

/*--------------------------------------------------------------------------*/

void charbonnier 

     (double   tau,       /* time step size, 0 < tau < 0.25 */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      double   hx,        /* pixel size in x direction */
      double   hy,        /* pixel size in y direction */
      double   lambda,    /* contrast parameter */
      double   **u)       /* input: current image; output: filtered image */

/*
  performs one explicit iteration step of a Perona-Malik filter with 
  Charbonnier diffusivity; 
*/

{
long    i, j;         /* loop variables */
double  **f;          /* u at old time level */
double  **g;          /* diffusivities */
double  fxp, fxm;     /* forward / backward differences in x direction */
double  fyp, fym;     /* forward / backward differences in y direction */
double  inv_hx;       /* 1.0 / hx */
double  inv_hy;       /* 1.0 / hy */
double  aux;          /* time saver */
double  rxx, ryy;     /* time savers */
double  grad_sqr;     /* |grad(f)|^2 */


/* ---- allocate memory ---- */

alloc_double_matrix (&f, nx+2, ny+2);
alloc_double_matrix (&g, nx+2, ny+2);


/* ---- copy u to f ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     f[i][j] = u[i][j];


/* ---- compute diffusivities g ---- */

/* compute time savers */
inv_hx = 1.0 / hx;    
inv_hy = 1.0 / hy;
aux    = 1.0 / (lambda * lambda);

/* create dummy boundaries of f by reflection */
dummies_double (f, nx, ny);

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* forward / backward differences of f in x- and y-direction */
     fxp = (f[i+1][j] - f[i][j]) * inv_hx;
     fxm = (f[i][j] - f[i-1][j]) * inv_hx;
     fyp = (f[i][j+1] - f[i][j]) * inv_hy;
     fym = (f[i][j] - f[i][j-1]) * inv_hy;

     /* compute squared gradient of f */
     grad_sqr = 0.5 * (fxp * fxp + fxm * fxm + fyp * fyp + fym * fym);

     /* compute Chrabonnier diffusivity g */
     g[i][j] = 1.0 / sqrt (1.0 + grad_sqr * aux);
     }


/* ---- perform one diffusion step ---- */

/* compute time savers */
rxx = tau / (2.0 * hx * hx);
ryy = tau / (2.0 * hy * hy);

/* create dummy boundaries of f and g by reflection */
dummies_double (f, nx, ny);
dummies_double (g, nx, ny);

/* diffuse */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     u[i][j] = f[i][j]
        + rxx * ( (g[i+1][j] + g[i][j]) * (f[i+1][j] - f[i][j])
                + (g[i-1][j] + g[i][j]) * (f[i-1][j] - f[i][j]) )
        + ryy * ( (g[i][j+1] + g[i][j]) * (f[i][j+1] - f[i][j])
                + (g[i][j-1] + g[i][j]) * (f[i][j-1] - f[i][j]) );


/* ---- free memory ---- */

free_double_matrix (f, nx+2, ny+2);
free_double_matrix (g, nx+2, ny+2);

return;

}  /* charbonnier */

/*--------------------------------------------------------------------------*/

double corr

      (long     nx,         /* image dimension in x direction */
       long     ny,         /* image dimension in y direction */
       double   **u,        /* input: first image; unchanged */
       double   **v)        /* input: second image; unchanged */

/*
  computes the correlation coefficient of u and v
*/

{
long    i, j;          /* loop variables */
double  mu_u, mu_v;    /* means of u and v */
double  var_u, var_v;  /* variances of u and v */
double  cov;           /* covariance of u and v */
double  aux;           /* denominator */ 


/* ---- compute the means of u and of v ---- */

mu_u = 0.0;
mu_v = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     mu_u = mu_u + u[i][j];
     mu_v = mu_v + v[i][j];
     }
mu_u = mu_u / (nx * ny);
mu_v = mu_v / (nx * ny);


/* ---- compute the variances of u and v ---- */

var_u = 0.0;
var_v = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     var_u = var_u + (u[i][j] - mu_u) * (u[i][j] - mu_u);
     var_v = var_v + (v[i][j] - mu_v) * (v[i][j] - mu_v);
     } 
var_u = var_u / (nx * ny);
var_v = var_v / (nx * ny);


/* ---- compute the covariance of u and v ---- */

/*
SUPPLEMENT CODE
*/

/* ---- compute the correlation coefficient of u and v ---- */

aux = sqrt (var_u * var_v);
if (aux > 0.0)
   return (cov / aux);
else
   /* one image is constant; correlation coefficient undefined */
   return (0.0);

}  /* corr */

/*--------------------------------------------------------------------------*/

int main ()

{
char    in[80];             /* for reading data */
char    out[80];            /* for reading data */
double  **f;                /* original image */
double  **uold;             /* image at old iteration */
double  **unew;             /* image at new iteration */
double  **v;                /* difference image ("method noise") */
long    i, j;               /* loop variables */
long    k;                  /* iteration counter */
long    nx, ny;             /* image size in x, y direction */ 
double  lambda;             /* contrast parameter */ 
double  tau;                /* time step size */ 
double  ccold, ccnew;       /* old and new correlation coefficients */
double  max, min;           /* largest, smallest grey value */
double  mean;               /* average grey value */
double  std;                /* standard deviation */
char    comments[1600];     /* string for comments */


printf ("\n");
printf ("PERONA-MALIK FILTER WITH CHARBONNIER DIFFUSIVITY\n");
printf ("EXPLICIT SCHEME, STOPPED WITH DECORRELATION CRITERION\n\n");
printf ("**************************************************\n\n");
printf ("    Copyright 2020 by Joachim Weickert            \n");
printf ("    Dept. of Mathematics and Computer Science     \n");
printf ("    Saarland University, Saarbruecken, Germany    \n\n");
printf ("    All rights reserved. Unauthorised usage,      \n");
printf ("    copying, hiring, and selling prohibited.      \n\n");
printf ("    Send bug reports to                           \n");
printf ("    weickert@mia.uni-saarland.de                  \n\n");
printf ("**************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf ("input image (pgm):                   ");
read_string (in);
read_pgm_to_double (in, &nx, &ny, &f);


/* ---- read parameters ---- */

printf ("contrast parameter lambda (> 0.0):   ");
read_double (&lambda);

printf ("time step size (< 0.25)):            ");
read_double (&tau);

printf ("output image (pgm):                  ");
read_string (out);
printf ("\n");


/* ---- allocate memory ---- */

alloc_double_matrix (&uold, nx+2, ny+2);
alloc_double_matrix (&unew, nx+2, ny+2);
alloc_double_matrix (&v,    nx+2, ny+2);


/* ---- analyse initial image ---- */

analyse_grey_double (f, nx, ny, &min, &max, &mean, &std);
printf ("initial image\n");
printf ("minimum:           %10.2lf \n", min);
printf ("maximum:           %10.2lf \n", max);
printf ("mean:              %10.2lf \n", mean);
printf ("standard dev.:     %10.2lf \n\n", std);


/* ---- initialisations for the loop ---- */

/* set iteration counter to 0 */
k = 0;

/* dummy setting to allow correlation reduction in first iteration */
ccnew = 1.0; 

/* copy f to unew */ 
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     unew[i][j] = f[i][j];


/* ---- iterate as long as the correlation coefficent decreases ---- */

do    {
      /* save ccnew in ccold, and unew in uold */
      ccold = ccnew;
      for (i=1; i<=nx; i++)
       for (j=1; j<=ny; j++)
           uold[i][j] = unew[i][j];

      /* perform one iteration with unew and increment k */
      charbonnier (tau, nx, ny, 1.0, 1.0, lambda, unew);
      k = k + 1;
 
      /* compute correlation coefficient ccnew between unew and f-unew */
      for (i=1; i<=nx; i++)
       for (j=1; j<=ny; j++)
           v[i][j] = f[i][j] - unew[i][j];
      ccnew = corr (nx, ny, unew, v);
      printf ("iteration: %5ld       correl. coeff.:   %8.6lf \n", 
              k, ccnew);
      }
while (fabs(ccnew) < fabs(ccold));


/* ---- analyse optimal result uold ---- */

analyse_grey_double (uold, nx, ny, &min, &max, &mean, &std);

printf ("\n");
printf ("correlation minimising result\n");
printf ("iteration:         %10ld   \n", k-1);
printf ("minimum:           %10.2lf \n", min);
printf ("maximum:           %10.2lf \n", max);
printf ("mean:              %10.2lf \n", mean);
printf ("standard dev.:     %10.2lf \n", std);
printf ("correl. coeff.:    %10.6lf \n\n", ccold);


/* ---- write uold to the output image (pgm format P5) ---- */

/* write parameter values in comment string */
comments[0] = '\0';
comment_line (comments, "# Perona-Malik filter, Charbonnier diffusivity\n");
comment_line (comments, "# explicit scheme\n");
comment_line (comments, "# stopped with decorrelation criterion\n");
comment_line (comments, "# initial image:  %s\n", in);
comment_line (comments, "# lambda:         %10.2lf\n", lambda);
comment_line (comments, "# tau:            %10.4lf\n", tau);
comment_line (comments, "# iterations:     %10ld\n", k-1);
comment_line (comments, "# min:            %10.2lf\n", min);
comment_line (comments, "# max:            %10.2lf\n", max);
comment_line (comments, "# mean:           %10.2lf\n", mean);
comment_line (comments, "# stand. dev.:    %10.2lf\n", std);
comment_line (comments, "# corr. coeff.:   %10.6lf\n", ccold);

/* write optimally denoised image data uold */
write_double_to_pgm (uold, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

free_double_matrix (f,    nx+2, ny+2);
free_double_matrix (uold, nx+2, ny+2);
free_double_matrix (unew, nx+2, ny+2);
free_double_matrix (v,    nx+2, ny+2);

return(0);

}  /* main */
