#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*      CORNER DETECTION WITH RESCALED AFFINE MORPHOLOGICAL SCALE-SPACE     */
/*                                                                          */
/*     (Copyright by Stephan Didas, 6/2008 and Joachim Weickert, 1/2022)    */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/* 
  Explicit scheme with central spatial differences. 
  Experimentally stable for time steps ht <= 0.1. 
  No numerical guarantee for an extremum principle. 
  Detects the maxima of the modulus of gradient norm weighted curvature.
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

void amss 

     (double   tau,       /* time step size */
      double   alpha,     /* dissipativity parameter in delta stencil */
      double   gamma,     /* nonnegativity parameter in delta stencil */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      double   h,         /* pixel size in x and y direction */
      double   **u)       /* image, changed */

/*
  performs one explicit step for AMSS with delta stencil
*/

{
long     i, j;                  /* loop variables */
double   **f;                   /* image u at old time level */
double   fx, fy;                /* Sobel derivatives */
double   f00, f11, f22, f33;    /* 2nd derivatives in principal directions */
double   inv_4h, inv_8h;        /* time savers */
double   inv_hh, inv_2hh;       /* time savers */
double   aux;                   /* time saver */
double   a, b, c;               /* weights for u_{xx}, u_{xy}, u_{yy} */
double   delta;                 /* free parameter in the delta stencil */


/* ---- allocate memory ---- */

alloc_double_matrix (&f, nx+2, ny+2);


/* ---- copy u to f ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     f[i][j] = u[i][j];


/* ---- reflecting dummy boundaries ---- */

dummies_double (f, nx, ny);


/* ---- compute time savers ---- */

inv_4h  = 1.0 / (4.0 * h);
inv_8h  = 1.0 / (8.0 * h);
inv_hh  = 1.0 / (h * h);
inv_2hh = 1.0 / (2.0 * h * h);
aux     = gamma * (1 - 2.0 * alpha);


/* ---- loop ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* compute Sobel derivatives */
     fx  =   (f[i+1][j+1] - f[i-1][j+1]) * inv_8h
           + (f[i+1][j]   - f[i-1][j]  ) * inv_4h
           + (f[i+1][j-1] - f[i-1][j-1]) * inv_8h;
     fy  =   (f[i+1][j+1] - f[i+1][j-1]) * inv_8h
           + (f[i]  [j+1] - f[i]  [j-1]) * inv_4h
           + (f[i-1][j+1] - f[i-1][j-1]) * inv_8h;

     /* compute important intermediate expressions */ 
     a =   fy * fy;
     b = - fx * fy;
     c =   fx * fx; 
     delta = alpha * (a + c) + aux * fabs(b);

     /* compute second order derivatives in the 4 principal directions */
     f00 = (f[i+1][j]   - 2.0 * f[i][j] + f[i-1][j])   * inv_hh; 
     f22 = (f[i][j+1]   - 2.0 * f[i][j] + f[i][j-1])   * inv_hh;
     f11 = (f[i+1][j+1] - 2.0 * f[i][j] + f[i-1][j-1]) * inv_2hh;   
     f33 = (f[i-1][j+1] - 2.0 * f[i][j] + f[i+1][j-1]) * inv_2hh;   
     
     /* explicit update step for AMSS */
     u[i][j]  = f[i][j] 
              + tau * cbrt ( (a - delta) * f00 + (delta + b) * f11
                           + (c - delta) * f22 + (delta - b) * f33 );
     }


/* ---- free memory ---- */
 
free_double_matrix (f, nx+2, ny+2);

return;

}  /* amss */

/*--------------------------------------------------------------------------*/

void find_max_curvature

     (long    nx,         /* image size in x direction */
      long    ny,         /* image size in y direction */
      long    hx,         /* spatial step size in x direction */
      long    hy,         /* spatial step size in y direction */
      long    *max_x,     /* x coordinate of maximum */
      long    *max_y,     /* y coordinate of maximum */
      double  **f)        /* input image, unchanged */

/* 
  Searches the point in the image with maximal 
  gradient norm weighted curvature 
*/

{
long    i, j;                    /* loop variables */
double  fx, fy, fxx, fxy, fyy;   /* derivatives */
double  two_hx;                  /* 2.0 * hx, time saver */
double  two_hy;                  /* 2.0 * hx, time saver */
double  hx_sqr;                  /* hx * hx, time saver */
double  hy_sqr;                  /* hy * hy, time saver */
double  two_hx_hy;               /* 2.0 * hx * hy, time saver */
double  grad_norm;               /* norm of the gradient */
double  max_curv;                /* maximal curvature modulus */
double  help;                    /* time saver */
long    old_max_x, old_max_y;


/* ---- create reflecting dummy boundaries for f ---- */

dummies_double (f, nx, ny);


/* ---- process ---- */

/* compute some time savers */
two_hx = 2.0 * hx;
two_hy = 2.0 * hy;
hx_sqr = hx * hx;
hy_sqr = hy * hy;
two_hx_hy = 2.0 * hx * hy;

max_curv = 0.0;
old_max_x = *max_x;
old_max_y = *max_y;

/* loop */
for (i=old_max_x-1; i<=old_max_x+1; i++)
 for (j=old_max_y-1; j<=old_max_y+1; j++)
     {
     /* central spatial derivatives */
     fx  = (f[i+1][j] - f[i-1][j]) / two_hx;
     fy  = (f[i][j+1] - f[i][j-1]) / two_hy;
     fxx = (f[i+1][j] - 2.0 * f[i][j] + f[i-1][j]) / hx_sqr;
     fyy = (f[i][j+1] - 2.0 * f[i][j] + f[i][j-1]) / hy_sqr;
     if (fx * fy < 0.0)
        fxy = (   f[i+1][j+1] - f[i][j+1] - f[i+1][j] + f[i][j]
                + f[i-1][j-1] - f[i][j-1] - f[i-1][j] + f[i][j] )
              / two_hx_hy;
     else
        fxy = ( - f[i-1][j+1] + f[i][j+1] + f[i+1][j] - f[i][j]
                - f[i+1][j-1] + f[i][j-1] + f[i-1][j] - f[i][j] )
              / two_hx_hy;

     /* find maximum curvature */
     grad_norm = sqrt (fx*fx + fy*fy);
     help = fabsf (fx * fx * fyy + fy * fy * fxx - 2.0 * fx * fy * fxy)
     / pow (grad_norm, 1.5);

     if (help > max_curv)
        {
        max_curv = help;
        *max_x = i;
        *max_y = j;
        }
     }

return;

}  /* find_max_curvature */

/*--------------------------------------------------------------------------*/

int main ()

{
char    in[80];               /* for reading data */
char    out[80];              /* for reading data */
double  **u;                  /* image */
long    nx, ny;               /* image size in x, y direction */ 
long    k;                    /* loop variable */ 
long    kmax;                 /* largest iteration number */
double  tau;                   /* time step size */
double  alpha;                /* dissipativity parameter in delta stencil */
double  gamma;                /* nonnegativity parameter in delta stencil */
double  max, min;             /* largest, smallest grey value */
double  mean;                 /* average grey value */
double  std;                  /* standard deviation */
char    comments[1600];       /* string for comments */
long    max_x, max_y;
double  time;
double  sumt;
double  sumt2;
long    samples;
double  sumx1;
double  sumx2;
double  sumx1t;
double  sumx2t;
double  d1, d2;
double  x1, x2;
double  velocity;

printf ("\n");
printf ("CORNER DETECTION WITH AFFINE MORPHOLOGICAL SCALE-SPACE\n\n");
printf ("****************************************************\n\n");
printf ("    Copyright 2008-2022 by                          \n");
printf ("    Joachim Weickert and Stephan Didas              \n");
printf ("    Dept. of Mathematics and Computer Science       \n");
printf ("    Saarland University, Saarbruecken, Germany      \n\n");
printf ("    All rights reserved. Unauthorized usage,        \n");
printf ("    copying, hiring, and selling prohibited.        \n\n");
printf ("    Send bug reports to                             \n");
printf ("    weickert@mia.uni-saarland.de                    \n\n");
printf ("****************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf ("input image (pgm):                         ");
read_string (in);
read_pgm_to_double (in, &nx, &ny, &u);


/* ---- read parameters ---- */

printf ("dissipativity parameter alpha (<=0.5):     ");
read_double (&alpha);

printf ("nonnegativity parameter gamma (in [-1,1]:  ");
read_double (&gamma);

printf ("time step size tau (<=0.01):               ");
read_double (&tau);

printf ("number of iterations (>0):                 ");
read_long (&kmax);

printf ("output image (pgm):                        ");
read_string (out);


/* ---- process image ---- */

sumt    = 0.0;
sumt2   = 0.0;
sumx1   = 0.0;
sumx2   = 0.0;
sumx1t  = 0.0;
sumx2t  = 0.0;
samples = 0;

max_x = 128;
max_y = 66; 

for (k=1; k<=kmax; k++)
    {
    /* perform one iteration */
    printf ("iteration number: %5ld \n", k);
    amss (tau, alpha, gamma, nx, ny, 1.0, u);

    time = pow (4.0/3.0 * tau * (double) k, 0.75) ;

    /* find maximum curvature */
    find_max_curvature (nx, ny, 1.0, 1.0, &max_x, &max_y, u);
    printf ("time:          %8.2f \n", time);
    printf ("point:       (%3ld, %3ld)\n\n", max_x, max_y);

    /* update counters and parameters */
    samples = samples + 1;
    sumt    = sumt + time;
    sumt2   = sumt2 + time * time;
    sumx1   = sumx1 + max_x;
    sumx2   = sumx2 + max_y;
    sumx1t  = sumx1t + max_x * time;
    sumx2t  = sumx2t + max_y * time;
    }


/* ---- analyse filtered image ---- */

analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);
printf ("filtered image\n");
printf ("minimum:          %8.2lf \n", min);
printf ("maximum:          %8.2lf \n", max);
printf ("mean:             %8.2lf \n", mean);
printf ("standard dev.:    %8.2lf \n", std);


/* ---- compute least squares fit ---- */

sumt  = sumt  / (double) samples;
sumx1 = sumx1 / (double) samples;
sumx2 = sumx2 / (double) samples;

d1 =   ( sumx1t - sumt * sumx1 * (double) samples )
     / ( sumt2  - sumt * sumt  * (double) samples );
d2 =   ( sumx2t - sumt * sumx2 * (double) samples )
     / ( sumt2  - sumt * sumt  * (double) samples );

x1 = sumx1 - d1 * sumt;
x2 = sumx2 - d2 * sumt;

velocity = sqrt(d1*d1 + d2*d2);

printf ("corner location:    (%4.2lf, %4.2lf)\n", x1, x2);
printf ("evolution vector:   (%4.2lf, %4.2lf)\n", d1, d2);
printf ("velocity:           %6.2lf\n", velocity);
printf ("corner angle:       %6.2lf\n\n", 
         360.0 * atan (1.0 / (velocity*velocity)) / M_PI);


/* ---- write output image (pgm format P5) ---- */

/* write parameter values in comment string */
comments[0] = '\0';
comment_line (comments, "# AMSS, explicit scheme with delta stencil\n");
comment_line (comments, "# initial image:    %s\n", in);
comment_line (comments, "# alpha:            %8.2lf\n", alpha);
comment_line (comments, "# gamma:            %8.2lf\n", gamma);
comment_line (comments, "# tau:              %8.2lf\n", tau);
comment_line (comments, "# kmax:             %8ld\n", kmax);
comment_line (comments, "# min:              %8.2lf\n", min);
comment_line (comments, "# max:              %8.2lf\n", max);
comment_line (comments, "# mean:             %8.2lf\n", mean);
comment_line (comments, "# standard dev.:    %8.2lf\n", std);
comment_line (comments, "# corner location:    (%4.2lf, %4.2lf)\n", x1, x2);
comment_line (comments, "# evolution vector:   (%4.2lf, %4.2lf)\n", d1, d2);
comment_line (comments, "# velocity:           %6.2lf\n", velocity);
comment_line (comments, "# corner angle:       %6.2lf\n\n",
                         360.0 * atan (1.0 / (velocity*velocity)) / M_PI); 

/* write image data */
write_double_to_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

free_double_matrix (u, nx+2, ny+2);

return(0);

}  /* main */
