#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*       EDGE AND COHERENCE ENHANCING ANISOTROPIC DIFFUSION FILTERING       */
/*                                                                          */
/*                   Copyright Joachim Weickert, 11/2021                    */
/*                         All rights reserved.                             */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 features:
 - explicit scheme with stability in the Euclidean norm
 - presmoothing at noise scale:  convolution-based, Neumann b.c.
 - presmoothing at integration scale: convolution-based, Dirichlet b.c.
 - space discretisation with the delta stencil 
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

void analyse_grey_double

     (double  **u,         /* image, unchanged */
      long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in y direction */
      double  *min,        /* minimum, output */
      double  *max,        /* maximum, output */
      double  *mean,       /* mean, output */
      double  *std,        /* standard deviation, output */
      double  *norm)       /* Euclidean norm, output */

/*
  computes minimum, maximum, mean, standard deviation, and Euclidean norm
  of a greyscale image u in double format
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

/* compute standard deviation and Euclidean norm */
*std  = 0.0;
*norm = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     *norm = *norm + u[i][j] * u[i][j];
     help2 = u[i][j] - *mean; 
     *std  = *std + help2 * help2;
     } 
*std  = sqrt (*std   / (nx * ny));
*norm = sqrt (*norm  / (nx * ny));

return;

}  /* analyse_grey_double */

/*--------------------------------------------------------------------------*/

void dummies_Dirichlet

     (double **v,        /* image */
      long   nx,         /* size in x direction */
      long   ny)         /* size in y direction */

/* 
  creates Dirichlet dummy boundaries 
*/

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
    {
    v[i][0]    = 0.0;
    v[i][ny+1] = 0.0;
    }

for (j=0; j<=ny+1; j++)
    {
    v[0][j]    = 0.0;
    v[nx+1][j] = 0.0;
    }

return;

}  /* dummies_Dirichlet */

/*--------------------------------------------------------------------------*/

void dummies_Neumann

     (double **v,        /* image */
      long   nx,         /* size in x direction */
      long   ny)         /* size in y direction */

/* 
   creates dummy boundaries by mirroring 
*/

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
    {
    v[i][0]    = v[i][1];
    v[i][ny+1] = v[i][ny];
    }

for (j=0; j<=ny+1; j++)
    {
    v[0][j]    = v[1][j];
    v[nx+1][j] = v[nx][j];
    }

return;

}  /* dummies_Neumann */

/*--------------------------------------------------------------------------*/

void gauss_conv

    (double   sigma,     /* standard deviation of the Gaussian */
     long     btype,     /* type of boundary condition */
     double   prec,      /* cutoff at precision * sigma */
     long     nx,        /* image dimension in x direction */
     long     ny,        /* image dimension in y direction */
     double   hx,        /* pixel size in x direction */
     double   hy,        /* pixel size in y direction */
     double   **u)       /* input: original image;  output: smoothed */


/*
  Gaussian convolution with a truncated and resampled Gaussian
*/


{
long    i, j, k, l, p;        /* loop variables */
long    length;               /* convolution vector: 0..length */
long    pmax;                 /* upper bound for p */
double  aux1, aux2;           /* time savers */
double  sum;                  /* for summing up */
double  *conv;                /* convolution vector */
double  *help;                /* row or column with dummy boundaries */


/* ----------------------- convolution in x direction -------------------- */

/* compute length of convolution vector */
length = (long)(prec * sigma / hx) + 1;

/* allocate memory for convolution vector */
alloc_double_vector (&conv, length+1);

/* compute entries of convolution vector */
aux1 = 1.0 / (sigma * sqrt(2.0 * 3.1415927));
aux2 = (hx * hx) / (2.0 * sigma * sigma);
for (i=0; i<=length; i++)
    conv[i] = aux1 * exp (- i * i * aux2);

/* normalisation */
sum = conv[0];
for (i=1; i<=length; i++)
    sum = sum + 2.0 * conv[i];
for (i=0; i<=length; i++)
    conv[i] = conv[i] / sum;

/* allocate memory for a row */
alloc_double_vector (&help, nx+length+length);

for (j=1; j<=ny; j++)
    {
    /* copy u in row vector */
    for (i=1; i<=nx; i++)
        help[i+length-1] = u[i][j];

    /* extend signal according to the boundary conditions */
    k = length;
    l = length + nx - 1;
    while (k > 0)
          {
          /* pmax = min (k, nx) */
          if (k < nx)
             pmax = k;
          else
             pmax = nx;

          /* extension on both sides */
          if (btype == 0)
             /* reflecting b.c.: symmetric extension */
             for (p=1; p<=pmax; p++)
                 {
                 help[k-p] = help[k+p-1];
                 help[l+p] = help[l-p+1];
                 }
          else
             /* Dirichlet b.c.: antisymmetric extension */
             for (p=1; p<=pmax; p++)
                 {
                 help[k-p] = - help[k+p-1];
                 help[l+p] = - help[l-p+1];
                 }

          /* update k and l */
          k = k - nx;
          l = l + nx;
          }

    /* convolution step */
    for (i=length; i<=nx+length-1; i++)
        {
        /* compute convolution */
        sum = conv[0] * help[i];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[i+p] + help[i-p]);
        /* write back */
        u[i-length+1][j] = sum;
        }
    } /* for j */

/* free memory */
free_double_vector (help, nx+length+length);
free_double_vector (conv, length + 1);


/* ----------------------- convolution in y direction -------------------- */

/* compute length of convolution vector */
length = (long)(prec * sigma / hy) + 1;

/* allocate memory for convolution vector */
alloc_double_vector (&conv, length + 1);

/* compute entries of convolution vector */
aux1 = 1.0 / (sigma * sqrt(2.0 * 3.1415927));
aux2 = (hy * hy) / (2.0 * sigma * sigma);
for (j=0; j<=length; j++)
    conv[j] = aux1 * exp (- j * j * aux2);

/* normalisation */
sum = conv[0];
for (j=1; j<=length; j++)
    sum = sum + 2.0 * conv[j];
for (j=0; j<=length; j++)
    conv[j] = conv[j] / sum;

/* allocate memory for a row */
alloc_double_vector (&help, ny+length+length);

for (i=1; i<=nx; i++)
    {
    /* copy u in column vector */
    for (j=1; j<=ny; j++)
        help[j+length-1] = u[i][j];

    /* extend signal according to the boundary conditions */
    k = length;
    l = length + ny - 1;
    while (k > 0)
          {
          /* pmax = min (k, ny) */
          if (k < ny)
             pmax = k;
          else
             pmax = ny;

          /* extension on both sides */
          if (btype == 0)
             /* reflecting b.c.: symmetric extension */
             for (p=1; p<=pmax; p++)
                 {
                 help[k-p] = help[k+p-1];
                 help[l+p] = help[l-p+1];
                 }
          else
             /* Dirichlet b.c.: antisymmetric extension */
             for (p=1; p<=pmax; p++)
                 {
                 help[k-p] = - help[k+p-1];
                 help[l+p] = - help[l-p+1];
                 }

          /* update k and l */
          k = k - ny;
          l = l + ny;
          }

    /* convolution step */
    for (j=length; j<=ny+length-1; j++)
        {
        /* compute convolution */
        sum = conv[0] * help[j];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[j+p] + help[j-p]);
        /* write back */
        u[i][j-length+1] = sum;
        }
    } /* for i */

/* free memory */
free_double_vector (help, ny+length+length);
free_double_vector (conv, length+1);

return;

} /* gauss_conv */

/* ------------------------------------------------------------------------ */

void struct_tensor 

     (double   **v,       /* image !! gets smoothed on exit !! */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      double   hx,        /* pixel size in x direction */
      double   hy,        /* pixel size in y direction */
      double   sigma,     /* noise scale */
      double   rho,       /* integration scale */
      double   **dxx,     /* element of structure tensor, output */
      double   **dxy,     /* element of structure tensor, output */
      double   **dyy)     /* element of structure tensor, output */

/*
  Computes the structure tensor at the staggered grid points
  [i+1/2][j+1/2] for i=1,...,nx-1 and j=1,...,ny-1.
*/

{
long     i, j;                 /* loop variables */
double   vxo, vxp, vyo, vyp;   /* derivatives of v */
double   w1, w2;               /* time savers */


/* ---- smoothing at noise scale, reflecting b.c. ---- */

if (sigma > 0.0) 
   gauss_conv (sigma, 0, 3.0, nx, ny, hx, hy, v);  


/* ---- building tensor product ---- */

w1 = 1.0 / hx;
w2 = 1.0 / hy;
dummies_Neumann (v, nx, ny);

for (i=1; i<=nx-1; i++)
 for (j=1; j<=ny-1; j++)
     {
     vxo = w1 * (v[i+1][j]   - v[i][j]); 
     vxp = w1 * (v[i+1][j+1] - v[i][j+1]);
     vyo = w2 * (v[i][j+1]   - v[i][j]);
     vyp = w2 * (v[i+1][j+1] - v[i+1][j]);
     dxx[i][j] = 0.5 * (vxo * vxo + vxp * vxp);
     dyy[i][j] = 0.5 * (vyo * vyo + vyp * vyp);  
     dxy[i][j] = 0.25 * (vxo + vxp) * (vyo + vyp);
     }


/* ---- smoothing at integration scale, Dirichlet b.c. ---- */

if (rho > 0.0)
   {
   gauss_conv (rho, 1, 3.0, nx-1, ny-1, hx, hy, dxx);  
   gauss_conv (rho, 1, 3.0, nx-1, ny-1, hx, hy, dxy);  
   gauss_conv (rho, 1, 3.0, nx-1, ny-1, hx, hy, dyy);  
   }

return;

}  /* struct_tensor */

/* ------------------------------------------------------------------------ */

void PA_trans 

     (double  a11,        /* coefficient (1,1) of a (2*2)-matrix */ 
      double  a12,        /* coefficient (1,2) of a (2*2)-matrix */ 
      double  a22,        /* coefficient (2,2) of a (2*2)-matrix */ 
      double  *c,         /* component 1 of first eigenvector, output */ 
      double  *s,         /* component 2 of first eigenvector, output */ 
      double  *lam1,      /* larger  eigenvalue, output */
      double  *lam2)      /* smaller eigenvalue, output */

/*
  Principal axis transformation, checked for correctness. 
*/

{
double  aux, norm;    /* time savers */ 


/* ---- compute eigenvalues and eigenvectors ---- */

aux = sqrt (pow (a22-a11, 2.0) + 4.0 * a12 * a12);

if (aux == 0.0)
   /* isotropic situation, eigenvectors arbitrary */
   {
   *lam1 = *lam2 = a11;
   *c = 1.0;
   *s = 0.0;
   }

else if (a11 > a22)
   { 
   *lam1 = 0.5 * (a11 + a22 + aux); 
   *lam2 = 0.5 * (a11 + a22 - aux); 
   *c = a11 - a22 + aux;
   *s = 2.0 * a12;
   }

else
   {
   *lam1 = 0.5 * (a11 + a22 + aux); 
   *lam2 = 0.5 * (a11 + a22 - aux); 
   *c = 2.0 * a12;
   *s = a22 - a11 + aux;
   }


/* ---- normalise eigenvectors ---- */

norm = sqrt (*c * *c + *s * *s);
if (norm >= 0.000001) 
   {
   *c = *c / norm;   
   *s = *s / norm;   
   }
else  
   {
   *c = 1.0;
   *s = 0.0;
   }

return;

}  /* PA_trans */

/* ----------------------------------------------------------------------- */

void PA_backtrans

     (double  c,        /* component 1 of first eigenvector */
      double  s,        /* component 2 of first eigenvector */
      double  lam1,     /* first eigenvalue */
      double  lam2,     /* second eigenvalue */
      double  *a11,     /* coefficients of (2*2)-matrix, output */
      double  *a12,     /* coefficients of (2*2)-matrix, output */
      double  *a22)     /* coefficients of (2*2)-matrix, output */


/*
  Principal axis backtransformation of a symmetric (2*2)-matrix.
  A = U * diag(lam1, lam2) * U_transpose with U = (v1 | v2).
  v1 = (c, s) is first eigenvector.
*/

{

*a11 = c * c * lam1 + s * s * lam2;
*a22 = lam1 + lam2 - *a11;             /* trace invariance */
*a12 = c * s * (lam1 - lam2);

return;

} /* PA_backtrans */

/*--------------------------------------------------------------------------*/

void diff_tensor 
     
     (double   lambda,   /* contrast parameter */
      long     nx,       /* image dimension in x direction */
      long     ny,       /* image dimension in y direction */
      double   **dxx,    /* in: structure tensor el., out: diff. tensor el. */
      double   **dxy,    /* in: structure tensor el., out: diff. tensor el. */ 
      double   **dyy)    /* in: structure tensor el., out: diff. tensor el. */ 

/*
  Computes the diffusion tensor of EED by means of the structure tensor.
*/

{
long    i, j;          /* loop variables */
double  aux;           /* time saver */
double  c, s;          /* specify first eigenvector */
double  mu1, mu2;      /* eigenvalues of structure tensor */
double  lam1, lam2;    /* eigenvalues of diffusion tensor */


aux = 1.0 / (lambda * lambda);

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* principal axis transformation */
     PA_trans (dxx[i][j], dxy[i][j], dyy[i][j], &c, &s, &mu1, &mu2);

     /* compute eigenvalues with Weickert diffusivity */
     if (mu1 > 0.0) 
        lam1 = 1.0 - exp (-3.31488 / pow (aux * mu1, 4.0));
     else 
        lam1 = 1.0;
     lam2 = 1.0;

     /* principal axis backtransformation */
     PA_backtrans (c, s, lam1, lam2, &dxx[i][j], &dxy[i][j], &dyy[i][j]); 
     }


/* ---- assign dummy boundaries (no flux) ---- */

dummies_Dirichlet (dxx, nx, ny);
dummies_Dirichlet (dxy, nx, ny);
dummies_Dirichlet (dyy, nx, ny);

return;

}  /* diff_tensor */

/*--------------------------------------------------------------------------*/

double sgn 

       (double    x)           /* argument */

/*
  sign function
*/

{
double  sign;                  /* auxiliary variable */

if (x > 0)
   sign = 1.0;
else if (x < 0)
   sign = -1.0;
else
   sign = 0.0;

return (sign);

}  /* sgn */

/*--------------------------------------------------------------------------*/

void weights 

     (double   **dxx,       /* entry [1,1] of structure tensor, unchanged */
      double   **dxy,       /* entry [1,2] of structure tensor, unchanged */
      double   **dyy,       /* entry [2,2] of structure tensor, unchanged */
      long     nx,          /* image dimension in x direction */ 
      long     ny,          /* image dimension in y direction */ 
      double   h,           /* pixel size in x- and y-direction */
      long     dtype,       /* type of discretisation */
      double   alpha,       /* dissipativity parameter */
      double   **woo,       /* weights in [i,j], output */
      double   **wpo,       /* weights in [i+1,j], output */
      double   **wmo,       /* weights in [i-1,j], output */
      double   **wop,       /* weights in [i,j+1], output */
      double   **wom,       /* weights in [i,j-1], output */
      double   **wpp,       /* weights in [i+1,j+1], output */
      double   **wmm,       /* weights in [i-1,j-1], output */
      double   **wpm,       /* weights in [i+1,j-1], output */
      double   **wmp)       /* weights in [i-1,j+1], output */

/*
  Computes delta stencil weights for the discrete divergence expression.
  The diffusion tensor entries are given on the staggered grid [i+1/2,j+1/2],
  and the stencil weights live on the regular grid. 
*/

{
long    i, j;       /* loop variables */
double  **delta;    /* space-variant parameter on staggered grid */
double  aux;        /* time saver */


/* ---- allocate memory ---- */

alloc_double_matrix (&delta, nx+1, ny+1);


/* ---- initialisations ---- */

aux = 1.0 / (2.0 * h * h);

for (i=0; i<=nx+1; i++)
 for (j=0; j<=ny+1; j++)
     {
     woo[i][j] = 0.0;
     wpo[i][j] = 0.0;
     wmo[i][j] = 0.0;
     wop[i][j] = 0.0;
     wom[i][j] = 0.0;
     wpp[i][j] = 0.0;
     wmm[i][j] = 0.0;
     wpm[i][j] = 0.0;
     wmp[i][j] = 0.0;
     }


/* ---- compute space-variant free parameter delta on staggered grid ---- */

if (dtype == 0)
   /* standard discretisation */
/*
   SUPPLEMENT CODE
*/
else
   /* WWW stencil with beta = (1.0 - 2.0 * alpha) * sgn (dxy[i][j]) */
/*
   SUPPLEMENT CODE
*/
 

/* ---- weights ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* weight in [i+1][j] with data from [i+1/2, j+1/2] and [i+1/2, j-1/2] */
     wpo[i][j] = aux * ( dxx[i][j]   - delta[i][j]
                       + dxx[i][j-1] - delta[i][j-1] );

     /* weight in [i-1][j] with data from [i-1/2, j+1/2] and [i-1/2, j-1/2] */
     wmo[i][j] = aux * ( dxx[i-1][j]   - delta[i-1][j]
                       + dxx[i-1][j-1] - delta[i-1][j-1] );

     /* weight in [i][j+1] with data from [i+1/2, j+1/2] and [i-1/2, j+1/2] */
     wop[i][j] = aux * ( dyy[i][j]   - delta[i][j] 
                       + dyy[i-1][j] - delta[i-1][j] );

     /* weight in [i][j-1] with data from [i+1/2, j-1/2] and [i-1/2, j-1/2] */
     wom[i][j] = aux * ( dyy[i][j-1]   - delta[i][j-1] 
                       + dyy[i-1][j-1] - delta[i-1][j-1] );

     /* weight in [i+1][j+1] with data from [i+1/2, j+1/2] */
     wpp[i][j] = aux * (dxy[i][j] + delta[i][j]);

     /* weight in [i-1][j-1] with data from [i-1/2, j-1/2] */
     wmm[i][j] = aux * (dxy[i-1][j-1] + delta[i-1][j-1]);

     /* weight in [i-1][j+1] with data from [i-1/2, j+1/2] */
     wmp[i][j] = aux * (delta[i-1][j] - dxy[i-1][j]);

     /* weight in [i+1][j-1] with data from [i+1/2, j-1/2] */
     wpm[i][j] = aux * (delta[i][j-1] - dxy[i][j-1]);
     }

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     woo[i][j] = - wpo[i][j] - wmo[i][j]
                 - wop[i][j] - wom[i][j]
                 - wpp[i][j] - wmm[i][j]
                 - wmp[i][j] - wpm[i][j];


/* ---- free memory ---- */

free_double_matrix (delta, nx+1, ny+1);

return;

}  /* weights */

/*--------------------------------------------------------------------------*/

void eced 

     (double   tau,         /* time step size */
      long     nx,          /* image dimension in x direction */ 
      long     ny,          /* image dimension in y direction */ 
      double   h,           /* pixel size in x and y direction */
      double   lambda,      /* contrast parameter */
      double   sigma,       /* noise scale */
      double   rho,         /* integration scale */
      long     dtype,       /* type of space discretisation */
      double   alpha,       /* dissipativity parameter */
      double   **u)         /* input: original image;  output: smoothed */


/* 
  Edge and coherence enhancing anisotropic diffusion filtering. 
  Explicit discretisation that is stable in the Euclidean norm. 
*/

{
long    i, j;                 /* loop variables */
double  **f;                  /* work copy of u */
double  **dxx, **dxy, **dyy;  /* entries of structure / diffusion tensor */
double  **woo;                /* weights for [i,j] */
double  **wpo;                /* weights for [i+1,j] */
double  **wmo;                /* weights for [i-1,j] */
double  **wop;                /* weights for [i,j+1] */
double  **wom;                /* weights for [i,j-1] */
double  **wpp;                /* weights for [i+1,j+1] */
double  **wmm;                /* weights for [i-1,j-1] */
double  **wpm;                /* weights for [i+1,j-1] */
double  **wmp;                /* weights for [i-1,j+1] */


/* ---- allocate memory ---- */

alloc_double_matrix (&f,   nx+2, ny+2);
alloc_double_matrix (&dxx, nx+1, ny+1);
alloc_double_matrix (&dxy, nx+1, ny+1);
alloc_double_matrix (&dyy, nx+1, ny+1);
alloc_double_matrix (&woo, nx+2, ny+2);
alloc_double_matrix (&wpo, nx+2, ny+2);
alloc_double_matrix (&wmo, nx+2, ny+2);
alloc_double_matrix (&wop, nx+2, ny+2);
alloc_double_matrix (&wom, nx+2, ny+2);
alloc_double_matrix (&wpp, nx+2, ny+2);
alloc_double_matrix (&wmm, nx+2, ny+2);
alloc_double_matrix (&wpm, nx+2, ny+2);
alloc_double_matrix (&wmp, nx+2, ny+2);


/* ---- copy u into f ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     f[i][j] = u[i][j];


/* ---- compute structure tensor on staggered grid (alters f!!!) ---- */

struct_tensor (f, nx, ny, h, h, sigma, rho, dxx, dxy, dyy);


/* ---- compute diffusion tensor on staggered grid ---- */

diff_tensor (lambda, nx-1, ny-1, dxx, dxy, dyy);


/* ---- compute stencil weights (on original grid) ---- */

weights (dxx, dxy, dyy, nx, ny, h, dtype, alpha, 
         woo, wpo, wmo, wop, wom, wpp, wmm, wpm, wmp);


/* ---- copy u into f and assign dummy boundaries ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     f[i][j] = u[i][j];
dummies_Neumann (f, nx, ny);


/* ---- explicit diffusion ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     u[i][j] = f[i][j] + tau * 
             ( woo[i][j] * f[i][j] 
             + wpo[i][j] * f[i+1][j] 
             + wmo[i][j] * f[i-1][j]  
             + wop[i][j] * f[i][j+1] 
             + wom[i][j] * f[i][j-1]
             + wpp[i][j] * f[i+1][j+1] 
             + wmm[i][j] * f[i-1][j-1] 
             + wpm[i][j] * f[i+1][j-1]
             + wmp[i][j] * f[i-1][j+1] ); 


/* ---- free memory ---- */

free_double_matrix (f,   nx+2, ny+2);
free_double_matrix (dxx, nx+1, ny+1);
free_double_matrix (dxy, nx+1, ny+1);
free_double_matrix (dyy, nx+1, ny+1);
free_double_matrix (woo, nx+2, ny+2);
free_double_matrix (wpo, nx+2, ny+2);
free_double_matrix (wmo, nx+2, ny+2);
free_double_matrix (wop, nx+2, ny+2);
free_double_matrix (wom, nx+2, ny+2);
free_double_matrix (wpp, nx+2, ny+2);
free_double_matrix (wmm, nx+2, ny+2);
free_double_matrix (wpm, nx+2, ny+2);
free_double_matrix (wmp, nx+2, ny+2);

return;

}  /* eced */

/*--------------------------------------------------------------------------*/

int main ()

{
char    in[80];               /* for reading data */
char    out[80];              /* for reading data */
double  **u;                  /* image */
long    k;                    /* loop variable */
long    nx, ny;               /* image size in x, y direction */
double  tau;                  /* time step size */
long    kmax;                 /* largest iteration number */
double  sigma;                /* noise scale */
double  rho;                  /* integration scale */
double  lambda;               /* contrast parameter */
long    dtype;                /* type of space discretisation */
double  alpha;                /* dissipativity parameter */
double  max, min;             /* largest, smallest grey value */
double  mean;                 /* average grey value */
double  std;                  /* standard deviation */
double  norm;                 /* Euclidean norm */
char    comments[1600];       /* string for comments */


printf("\n");
printf("EDGE AND COHERENCE-ENHANCING ANISOTROPIC DIFFUSION (ECED)\n");
printf("EXPLICIT SCHEME WITH DELTA STENCIL\n\n");
printf("***************************************************\n\n");
printf("    Copyright 2021 by Joachim Weickert         \n");
printf("    Faculty of Mathematics and Computer Science\n");
printf("    Saarland University, Germany               \n\n");
printf("    All rights reserved. Unauthorized usage,   \n");
printf("    copying, hiring, and selling prohibited.   \n\n");
printf("    Send bug reports to                        \n");
printf("    weickert@mia.uni-saarland.de               \n\n");
printf("***************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf ("input image (pgm):                       ");
read_string (in);
read_pgm_to_double (in, &nx, &ny, &u);


/* ---- read other parameters ---- */

printf ("contrast parameter lambda (>0):          ");
read_double (&lambda);

printf ("noise scale sigma (>=0):                 ");
read_double (&sigma);

printf ("integration scale rho (>=0):             ");
read_double (&rho);

printf ("type of space discretisation:\n");
printf (" (0) standard discretisation\n");
printf (" (1) WWW stencil with optimal beta\n");
printf ("your choice:                             ");
read_long (&dtype);

alpha = 0;  /* dummy setting to avoid compiler warnings */
if (dtype == 1)
   {
   printf ("dissipativity parameter alpha (<=0.5):   ");
   read_double (&alpha);
   }

printf ("time step size tau (<%5.3lf):             ", 1.0/(4.0*(1.0-alpha)));
read_double (&tau);

printf ("number of iterations:                    ");
read_long (&kmax);

printf ("output image:                            ");
read_string (out);

printf("\n");
printf("***************************************************\n\n");


/* ---- process image ---- */

analyse_grey_double (u, nx, ny, &min, &max, &mean, &std, &norm);
printf ("initial image\n");
printf ("minimum:          %8.2f \n", min);
printf ("maximum:          %8.2f \n", max);
printf ("mean:             %8.2f \n", mean);
printf ("standard dev.:    %8.2f \n", std);
printf ("Euclidean norm:   %8.2f \n\n", norm);

for (k=1; k<=kmax; k++)
    {
    /* perform one iteration */
    eced (tau, nx, ny, 1.0, lambda, sigma, rho, dtype, alpha, u);

    /* analyse result */
    analyse_grey_double (u, nx, ny, &min, &max, &mean, &std, &norm);
    printf ("iteration number: %8ld \n", k);
    printf ("minimum:          %8.2f \n", min);
    printf ("maximum:          %8.2f \n", max);
    printf ("mean:             %8.2f \n", mean);
    printf ("standard dev.:    %8.2f \n", std);
    printf ("Euclidean norm:   %8.2f \n\n", norm);
    } /* for */


/* ---- write output image (pgm format P5) ---- */

/* write parameter values in comment string */
comments[0] = '\0';
comment_line (comments, "# ECED, explicit scheme, stable in Euclidean norm\n");
if (dtype == 0)
   comment_line (comments, "# standard discretisation\n");
else
   comment_line (comments, "# WWW stencil with free alpha and optimal beta\n");
comment_line (comments, "# initial image:  %s\n", in);
comment_line (comments, "# lambda:         %8.2f\n", lambda);
comment_line (comments, "# sigma:          %8.2f\n", sigma);
comment_line (comments, "# rho:            %8.2f\n", rho);
if (dtype == 1)
   comment_line (comments, "# alpha:          %8.2f\n", alpha);
comment_line (comments, "# tau:            %8.2f\n", tau);
comment_line (comments, "# iterations:     %8ld\n",  kmax);
comment_line (comments, "# minimum:        %8.2f\n", min);
comment_line (comments, "# maximum:        %8.2f\n", max);
comment_line (comments, "# mean:           %8.2f\n", mean);
comment_line (comments, "# standard dev.:  %8.2f\n", std);
comment_line (comments, "# Euclidean norm: %8.2f\n", norm);

/* write image data */
write_double_to_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

free_double_matrix (u, nx+2, ny+2);

return(0);

}  /* main */
