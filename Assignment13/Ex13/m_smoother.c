#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                      PDE EVOLUTIONS FOR M-SMOOTHERS                      */
/*                                                                          */
/*                   (Copyright Joachim Weickert, 1/2022)                   */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/* 
  computes PDE evolutions of M-smoothers;
  uses Osher-Rudin minmod scheme for stabilised inverse diffusion
  and an upwind scheme for mean curvature motion 
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

free (vector);
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
    free (matrix[i]);

free (matrix);

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

void multiple_layer_dummies_double

     (double  **f,          /* input image, unchanged */
      long    nx,           /* size in x direction */
      long    ny,           /* size in y direction */
      long    m,            /* size of boundary layer */
      double  **u)          /* output image, changed */

/*
  Copies a double format image f with pixel range [1, nx] * [1, ny] 
  to a double format image u with pixel range [m, nx+m-1] * [m, ny+m-1] 
  and adds a boundary layer of size m.
  This requires that the memory for u is allocated in a larger range:
  alloc_double_matrix (&u, nx+2*m, ny+2*m)
*/

{
long  i, j, n;   /* loop variables */


/* ---- copy f to u (with shift) ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     u[i+m-1][j+m-1] = f[i][j];
/* now u is specified for i=m..nx+m-1 and j=m..ny+m-1 */


/* ---- create dummy layer of size m for u ---- */

for (n=1; n<=m; n++)
    /* create an image in the range i=m-n..nx+m+n-1 and j=m-n..ny+m+n-1 */
    {
    /* copy column m+n-1 into m-n, and column nx+m-n into nx+m+n-1 */
    for (j=m+1-n; j<=ny+m+n-2; j++)
        {
        u[m-n][j]      = u[m+n-1][j];
        u[nx+m+n-1][j] = u[nx+m-n][j];
        }

    /* copy row m+n-1 into m-n, and row ny+m-n into ny+m+n-1 */
    for (i=m-n; i<=nx+m+n-1; i++)
        {
        u[i][m-n]      = u[i][m+n-1];
        u[i][ny+m+n-1] = u[i][ny+m-n];
        }
    } /* for n */

return;

}  /* multiple_layer_dummies_double */

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

double sign 

       (double  a)           /* argument of the sign function */

/* 
  sign function
*/

{
double result;

if (a > 0)
   result = 1.0;
else if (a < 0)
   result = -1.0;
else
   result = 0.0;

return (result);

}  /* sign */

/*--------------------------------------------------------------------------*/

double min3 

       (double  a,           /* first argument */
        double  b,           /* second argument */
        double  c)           /* third argument */

/* 
  returns the minimum of a, b, and c 
*/

{
double result;

result = a;
if (b < result) result = b;
if (c < result) result = c;

return (result);

}  /* min3 */

/*--------------------------------------------------------------------------*/

double max2

       (double  a,           /* first argument */
        double  b)           /* second argument */

/*
  returns the maximum of a and b
*/

{
double result;

result = a;
if (b > result) result = b;

return (result);

}  /* max2 */

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
   
       (double  a)

/* 
  returns the square of a 
*/

{
double result;

result = a * a;

return (result);

}  /* sqr */

/*--------------------------------------------------------------------------*/

double minmod

       (double  a,           /* first argument */
        double  b,           /* second argument */
        double  c)           /* third argument */

/*
  minmod function:
  = 0,                       if any of the a, b, c are 0 or of opposite sign
  = sign(a) min(|a|,|b|,|c|) else
*/

{
double result;

if ((sign(a) == sign(b)) && (sign(b) == sign(c)) && (a != 0.0))
   result = sign(a) * min3 (fabs(a), fabs(b), fabs(c));
else
   result = 0.0;

return (result);

}  /* minmod */

/*--------------------------------------------------------------------------*/

void diff_regular 

     (double   tau,       /* time step size */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      double   h,         /* pixel size in x, y direction */
      double   delta,     /* diagonal weight */
      double   **u)       /* input: original image ;  output: smoothed */

/*
  one explicit homogeneous diffusion step on a 9-point-stencil;
  regular discretisation
*/

{
long     i, j;        /* loop variables */
double   aux1, aux2;  /* time savers */
double   fxp, fxm;    /* one-sided differences, direction x=(1,0) */
double   fyp, fym;    /* one-sided differences, direction y=(0,1) */
double   fdp, fdm;    /* one-sided differences, direction d=(1,1) */
double   fep, fem;    /* one-sided differences, direction e=(1,-1) */
double   **f;         /* work copy of u */
 
/* allocate memory */
alloc_double_matrix (&f, nx+2, ny+2);

/* copy u to f */
for (i=1; i<=nx; i++) 
 for (j=1; j<=ny; j++)
     f[i][j] = u[i][j];
 
/* assign dummy boundary conditions */
dummies_double (f, nx, ny);

/* compute time savers */
aux1 = tau * (1.0 - delta) / (h * h);
aux2 = tau * delta / (2.0 * h * h);

/* diffusion step */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* one-sided differences */
     fxp = f[i+1][j]   - f[i]  [j];
     fxm = f[i]  [j]   - f[i-1][j];
     fyp = f[i]  [j+1] - f[i]  [j];
     fym = f[i]  [j]   - f[i]  [j-1];
     fdp = f[i+1][j+1] - f[i]  [j];
     fdm = f[i]  [j]   - f[i-1][j-1];
     fep = f[i+1][j-1] - f[i]  [j];
     fem = f[i]  [j]   - f[i-1][j+1];

     u[i][j] = f[i][j] + aux1 * (fxp - fxm + fyp - fym) 
                       + aux2 * (fdp - fdm + fep - fem);
     }

/* free memory */
free_double_matrix (f, nx+2, ny+2);

return;

}  /* diff_regular */

/*--------------------------------------------------------------------------*/

void diff_minmod

     (double   tau,       /* time step size */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      double   h,         /* pixel size in x, y direction */
      double   delta,     /* diagonal weight */
      double   **u)       /* input: original image; output: smoothed */

/*
  One explicit step for stabilised homogeneous diffusion in a 5 * 5 stencil; 
  Osher-Rudin minmod scheme.
  To perform backward diffusion, it needs a negative time step size.
*/

{
long     i, j;                    /* loop variables */
double   aux1, aux2;              /* time savers */
double   fxpp, fxp, fxm, fxmm;    /* one-sided differences, direction x */
double   fypp, fyp, fym, fymm;    /* one-sided differences, direction y */
double   fdpp, fdp, fdm, fdmm;    /* one-sided differences, direction (1,1) */
double   fepp, fep, fem, femm;    /* one-sided differences, direction (1,-1) */
double   **f;                     /* work copy of u */

/* allocate memory */
alloc_double_matrix (&f, nx+4, ny+4);

/* compute time savers */
aux1 = tau * (1.0 - delta) / (h * h);
aux2 = tau * delta / (2.0 * h * h);

/* copy u to f (shifted) and add a dummy layer of size 2 */
multiple_layer_dummies_double (u, nx, ny, 2, f);
/* f is now active in [2, nx+1] * [2, ny+1] */

/* stabilised minmod diffusion step */
for (i=2; i<=nx+1; i++)
 for (j=2; j<=ny+1; j++)
     {
     /* one-sided differences */
     fxpp = f[i+2][j]   - f[i+1][j];
     fxp  = f[i+1][j]   - f[i]  [j];
     fxm  = f[i]  [j]   - f[i-1][j];
     fxmm = f[i-1][j]   - f[i-2][j];
     fypp = f[i]  [j+2] - f[i]  [j+1];
     fyp  = f[i]  [j+1] - f[i]  [j];
     fym  = f[i]  [j]   - f[i]  [j-1];
     fymm = f[i]  [j-1] - f[i]  [j-2];
     fdpp = f[i+2][j+2] - f[i+1][j+1];
     fdp  = f[i+1][j+1] - f[i]  [j];
     fdm  = f[i]  [j]   - f[i-1][j-1];
     fdmm = f[i-1][j-1] - f[i-2][j-2];
     fepp = f[i+2][j-2] - f[i+1][j-1];
     fep  = f[i+1][j-1] - f[i]  [j];
     fem  = f[i]  [j]   - f[i-1][j+1];
     femm = f[i-1][j+1] - f[i-2][j+2];

     /* diffusion step with back shift */
/*
     SUPPLEMENT CODE HERE:
     u[i-1][j-1] = f[i][j] + aux1 * (...) + aux2 * (...);
*/
     }

/* free memory */
free_double_matrix (f, nx+4, ny+4);

return;

}  /* diff_minmod */

/*--------------------------------------------------------------------------*/

void curvature

     (double   **f,       /* image, unchanged */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      double   h,         /* pixel size in x, y direction */
      double   limit,     /* curvature limiter */
      double   **curv)    /* curvature, output */

/*
  finite difference approximation of the curvature using
  the second order directional derivative along the gradient;
  limits its range to [-limit, limit]
*/

{
long     i, j;                        /* loop variables */
double   grad_sqr;                    /* |grad(f)|^2, time saver */
double   fx, fy, fxx, fxy, fyy;       /* central differences */
double   inv_2_h, inv_h_sqr;          /* time savers */
double   inv_2_h_sqr, inv_4_h_sqr;    /* time savers */

/* create reflecting dummy boundaries */
dummies_double (f, nx, ny);

/* compute time savers */
inv_2_h     = 1.0 / (2.0 * h);
inv_h_sqr   = 1.0 / (h * h); 
inv_2_h_sqr = 1.0 / (2.0 * h * h);
inv_4_h_sqr = 1.0 / (4.0 * h * h);

/* compute curvature */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* central differences */
     fx  = (f[i+1][j] - f[i-1][j]) * inv_2_h;
     fy  = (f[i][j+1] - f[i][j-1]) * inv_2_h;
     fxx = (f[i+1][j] - 2.0 * f[i][j] + f[i-1][j]) * inv_h_sqr;
     fyy = (f[i][j+1] - 2.0 * f[i][j] + f[i][j-1]) * inv_h_sqr;
     if (fx * fy < 0.0)
        fxy = (   f[i+1][j+1] - f[i][j+1] - f[i+1][j] + f[i][j] 
                + f[i-1][j-1] - f[i][j-1] - f[i-1][j] + f[i][j] ) 
              * inv_2_h_sqr;
     else if (fx * fy > 0.0)
        fxy = ( - f[i-1][j+1] + f[i][j+1] + f[i-1][j] - f[i][j] 
                - f[i+1][j-1] + f[i][j-1] + f[i+1][j] - f[i][j] )
              * inv_2_h_sqr;
     else
        fxy = ( f[i+1][j+1] - f[i+1][j-1] - f[i-1][j+1] + f[i-1][j-1] )
              * inv_4_h_sqr;

     /* compute (slightly regularised) curvature */
     grad_sqr   = fx * fx + fy * fy + 1.0e-10;
     curv[i][j] = (fx * fx * fyy + fy * fy * fxx - 2.0 * fx * fy * fxy)
                  / (grad_sqr * sqrt (grad_sqr));
     }


/* ---- limit curvature ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     if (curv[i][j] > limit)
        curv[i][j] = limit;
     else if (curv[i][j] < -limit)
        curv[i][j] = -limit;

return;

}  /* curvature */

/*--------------------------------------------------------------------------*/

void mcm 

     (double   tau,       /* time step size */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      double   h,         /* pixel size in x, y direction */
      double   delta,     /* diagonal weight */
      double   **u)       /* input: original image; output: smoothed */

/*
  one explicit mean curvature motion step;
  upwind scheme with curvature limiter
*/

{
long     i, j;          /* loop variables */
double   fxp, fxm;      /* one-sided differences, direction x=(1,0) */
double   fyp, fym;      /* one-sided differences, direction y=(0,1) */
double   fdp, fdm;      /* one-sided differences, direction d=(1,1) */
double   fep, fem;      /* one-sided differences, direction e=(1,-1) */
double   aux1, aux2;    /* time savers */
double   **curv;        /* isophote curvature */
double   **f;           /* work copy of u */
 
/* allocate memory */
alloc_double_matrix (&f,    nx+2, ny+2);
alloc_double_matrix (&curv, nx+2, ny+2);

/* copy u to f */
for (i=1; i<=nx; i++) 
 for (j=1; j<=ny; j++)
     f[i][j] = u[i][j];
 
/* assign dummy boundary conditions */
dummies_double (f, nx, ny);

/* compute curvature and limit it to [-2/h, 2/h] */
curvature (f, nx, ny, h, 2.0 / h, curv);

/* compute time savers */
aux1 = tau * (1.0 - delta) / h;
aux2 = tau * delta / (h * sqrt (2.0));

/* MCM step */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* one-sided differences */
     fxp = f[i+1][j]   - f[i]  [j];
     fxm = f[i]  [j]   - f[i-1][j];
     fyp = f[i]  [j+1] - f[i]  [j];
     fym = f[i]  [j]   - f[i]  [j-1];
     fdp = f[i+1][j+1] - f[i]  [j];
     fdm = f[i]  [j]   - f[i-1][j-1];
     fep = f[i+1][j-1] - f[i]  [j];
     fem = f[i]  [j]   - f[i-1][j+1];

     /* MCM step */
     if (tau * curv[i][j] >= 0.0)
        u[i][j] = f[i][j] + curv[i][j] *
                  ( aux1 * sqrt ( sqr (max3 (0.0, -fxm, fxp)) 
                                + sqr (max3 (0.0, -fym, fyp)))
                  + aux2 * sqrt ( sqr (max3 (0.0, -fdm, fdp))
                                + sqr (max3 (0.0, -fem, fep))) ); 
     else
        u[i][j] = f[i][j] + curv[i][j] *
                  ( aux1 * sqrt ( sqr (max3 (0.0, fxm, -fxp)) 
                                + sqr (max3 (0.0, fym, -fyp)))
                  + aux2 * sqrt ( sqr (max3 (0.0, fdm, -fdp))
                                + sqr (max3 (0.0, fem, -fep))) );
     }

/* free memory */
free_double_matrix (f,    nx+2, ny+2);
free_double_matrix (curv, nx+2, ny+2);

return;

}  /* mcm */

/*--------------------------------------------------------------------------*/

void m_smoother 

     (double   tau,       /* time step size */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      double   h,         /* pixel size in x and y direction */
      double   wi,        /* isophote weight */
      double   wf,        /* flowline weight */
      double   delta,     /* diagonal weight */
      double   **u)       /* input: original image; output: smoothed */

/*
  computes one time step of the PDE evolution for M-smoothers; 
  uses upwind scheme for mean curvature motion,
  standard five point stencil for homogeneous forward diffusion,
  and Osher-Rudin discretisation of stabilised inverse diffusion;
*/

{

/* mcm step */
mcm ((wi-wf) * tau, nx, ny, h, delta, u);

/* diffusion step */
if (wf > 0.0)
   /* (forward) diffusion step with regular stencil */
   diff_regular (wf * tau, nx, ny, h, delta, u);
else if (wf < 0.0)
   /* (backward) diffusion step with minmod stencil */
   diff_minmod (wf * tau, nx, ny, h, delta, u);

return;

} /* m_smoother */

/*------------------------------------------------------------------------*/

int main ()

{
char    in[80];               /* for reading data */
char    out[80];              /* for reading data */
double  **u;                  /* original image */
double  wi, wf;               /* isophote / flowline weight */
double  tau;                  /* time step size */
double  tau_max;              /* time step size limit */
long    k;                    /* loop variable */
long    kmax;                 /* number of iterations */
long    nx, ny;               /* image size in x, y direction */
double  delta;                /* diagonal weight */
double  max, min;             /* largest, smallest grey value */
double  mean;                 /* average grey value */
double  std;                  /* standard deviation */
double  h;                    /* grid size */
double  aux;                  /* auxiliary variable */
char    comments[1600];       /* string for comments */


printf ("\n");
printf ("PDE EVOLUTION OF M-SMOOTHERS FOR GREYSCALE (PGM) IMAGES\n");
printf ("UPWIND/MINMOD SCHEME WITH MAXIMUM-MINIMUM PRINCIPLE\n\n");
printf ("*********************************************************\n\n");
printf ("    Copyright 2022 by Joachim Weickert                \n");
printf ("    Faculty of Mathematics and Computer Science       \n");
printf ("    Saarland University, Germany                      \n\n");
printf ("    All rights reserved. Unauthorized usage, copying, \n");
printf ("    hiring, and selling prohibited.                   \n\n");
printf ("    Send bug reports to weickert@mia.uni-saarland.de  \n\n");
printf ("*********************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf("input image (pgm):                  ");
read_string (in);
read_pgm_to_double (in, &nx, &ny, &u);


/* ---- read other parameters ---- */

printf ("isophote weight (for u_xi_xi):      ");
read_double (&wi);

printf ("flowline weight (for u_eta_eta):    ");
read_double (&wf);

printf ("diagonal weight delta (in [0, 1]):  ");
read_double (&delta);

/* compute time step size limit */
h = 1.0;
aux = max2 (fabs (wi - wf) * 2.0 * (sqrt(2.0) - (sqrt(2.0) - 1.0) * delta),
            fabs (wf) * (4.0 - 2.0 * delta));
tau_max = h * h / aux;

printf ("time step size (in (0, %6.4lf)):    ", 0.9999 * tau_max);
read_double (&tau);
 
printf ("number of iterations (> 0):         ");
read_long (&kmax);

printf ("output image (pgm):                 ");
read_string (out);

printf ("\n");


/* ---- analyse initial image ---- */

analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);
printf ("initial image\n");
printf ("minimum:          %8.2lf \n", min);
printf ("maximum:          %8.2lf \n", max);
printf ("mean:             %8.2lf \n", mean);
printf ("standard dev.:    %8.2lf \n\n", std);


/* ---- PDE evolution ---- */

for (k=1; k<=kmax; k++)
    {
    m_smoother (tau, nx, ny, h, wi, wf, delta, u);

    /* check minimum, maximum, mean, standard deviation */
    analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);
    printf ("iteration:        %8ld \n", k);
    printf ("minimum:          %8.2lf \n", min);
    printf ("maximum:          %8.2lf \n", max);
    printf ("mean:             %8.2lf \n", mean);
    printf ("standard dev.:    %8.2lf \n\n", std);
    }

/* ---- write output image (pgm format P5) ---- */

/* write parameter values in comment string */
comments[0] = '\0';
comment_line (comments, "# PDE evolution for M-smoothers\n");
comment_line (comments, "# initial image:      %s\n", in);
comment_line (comments, "# isophote weight:   %9.4lg\n", wi);
comment_line (comments, "# flowline weight:   %9.4lg\n", wf);
comment_line (comments, "# delta:             %9.4lg\n", delta);
comment_line (comments, "# tau:               %9.4lg\n", tau);
comment_line (comments, "# kmax:              %9ld\n", kmax);
comment_line (comments, "# min:               %9.2lf\n", min);
comment_line (comments, "# max:               %9.2lf\n", max);
comment_line (comments, "# mean:              %9.2lf\n", mean);
comment_line (comments, "# standard dev.:     %9.2lf\n", std);

/* write image data */
write_double_to_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory ---- */

free_double_matrix (u, nx+2, ny+2);

return(0);
}
