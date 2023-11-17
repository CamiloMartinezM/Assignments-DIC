#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*               FILTER OF CATTE ET AL. FOR GREYSCALE IMAGES                */
/*                                                                          */
/*                 (Copyright by Joachim Weickert, 11/2022)                 */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/* 
  Perona-Malik filter with regularisation of Catte et al. (1992);
  explicit finite difference scheme;
  presmoothing by convolution with a truncated Gaussian
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
   //sprintf (comment, "%s\n", comment);

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
   //fprintf (outimage, comments);             /* comments */
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

double MSE

       (double   **u,       /* processed image */
        double   **f,       /* reference image */
        long     nx,        /* pixel number in x direction */
        long     ny)        /* pixel number in y direction */

/*
  mean squared error
*/

{
long    i, j;                 /* loop variables */
double  sum;                  /* mean squared error */

sum = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     sum = sum + (u[i][j] - f[i][j]) * (u[i][j] - f[i][j]);
sum = sum / (nx * ny);

return (sum);

}  /* MSE */

/*--------------------------------------------------------------------------*/

void gauss_conv

    (double   sigma,     /* standard deviation of the Gaussian */
     double   prec,      /* cutoff at precision * sigma */
     long     nx,        /* image dimension in x direction */
     long     ny,        /* image dimension in y direction */
     double   hx,        /* pixel size in x direction */
     double   hy,        /* pixel size in y direction */
     double   **u)       /* input: original image ;  output: smoothed */


/*
  Gaussian convolution with a truncated and resampled Gaussian
  and reflecting boundary conditions
*/


{
long    i, j, k, l, p;   /* loop variables */
long    length;          /* convolution vector: 0..length */
long    pmax;            /* upper bound for p */
double  aux1, aux2;      /* time savers */
double  sum;             /* for summing up */
double  *conv;           /* convolution vector */
double  *help;           /* row or column with dummy boundaries */


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
alloc_double_vector (&help, nx+2*length);

for (j=1; j<=ny; j++)
    {
    /* copy u in row vector */
    for (i=1; i<=nx; i++)
        help[i+length-1] = u[i][j];

    /* extend signal by mirroring */
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
          for (p=1; p<=pmax; p++)
              {
              help[k-p] = help[k+p-1];
              help[l+p] = help[l-p+1];
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
free_double_vector (help, nx+2*length);
free_double_vector (conv, length+1);


/* ----------------------- convolution in y direction -------------------- */

/* compute length of convolution vector */
length = (long)(prec * sigma / hy) + 1;

/* allocate memory for convolution vector */
alloc_double_vector (&conv, length+1);

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
alloc_double_vector (&help, ny+2*length);

for (i=1; i<=nx; i++)
    {
    /* copy u in column vector */
    for (j=1; j<=ny; j++)
        help[j+length-1] = u[i][j];

    /* extend signal by mirroring */
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
          for (p=1; p<=pmax; p++)
              {
              help[k-p] = help[k+p-1];
              help[l+p] = help[l-p+1];
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
free_double_vector (help, ny+2*length);
free_double_vector (conv, length+1);

return;

} /* gauss_conv */

/*--------------------------------------------------------------------------*/

void catte 

     (double   tau,       /* time step size, 0 < tau < 0.25 */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      double   hx,        /* pixel size in x direction */
      double   hy,        /* pixel size in y direction */
      double   lambda,    /* contrast parameter */
      double   sigma,     /* noise scale */
      double   **u)       /* input: current image; output: filtered image */

/*
  performs one explicit iteration step of a Perona-Malik filter with 
  exponential diffusivity and Catte regularisation; 
*/

{
long    i, j;         /* loop variables */
double  **f;          /* u at old time level */
double  **v;          /* Gaussian-smoothed version of u */
double  **g;          /* diffusivities */
double  vx, vy;       /* central differences in x, y direction */
double  inv_two_hx;   /* 1.0 / (2.0 * hx) */
double  inv_two_hy;   /* 1.0 / (2.0 * hy) */
double  aux;          /* time saver */
double  rxx, ryy;     /* time savers */
double  grad_sqr;     /* |grad(f)|^2 */


/* ---- allocate memory ---- */

alloc_double_matrix (&f, nx+2, ny+2);
alloc_double_matrix (&v, nx+2, ny+2);
alloc_double_matrix (&g, nx+2, ny+2);


/* ---- copy u into f and v ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     v[i][j] = f[i][j] = u[i][j];


/* ---- regularise v by Gaussian convolution ---- */

if (sigma > 0.0)
   gauss_conv (sigma, 5.0, nx, ny, hx, hy, v);


/* ---- compute diffusivities g ---- */

/* compute time savers */
inv_two_hx = 1.0 / (2.0 * hx);    
inv_two_hy = 1.0 / (2.0 * hy);
aux        = 1.0 / (2.0 * lambda * lambda);

/* create dummy boundaries of v by reflection */
dummies_double (v, nx, ny);

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* central differences of v in x- and y-direction */
     vx = (v[i+1][j] - v[i-1][j]) * inv_two_hx;
     vy = (v[i][j+1] - v[i][j-1]) * inv_two_hy;

     /* compute squared gradient of v */
     /*
       SUPPLEMENT CODE
     */

     /* compute exponential Perona-Malik diffusivity g */
     /*
       SUPPLEMENT CODE
     */
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
free_double_matrix (v, nx+2, ny+2);
free_double_matrix (g, nx+2, ny+2);

return;

} /* catte */

/*--------------------------------------------------------------------------*/

int main ()

{
char    in1[80], in2[80];   /* for reading data */
char    out[80];            /* for reading data */
double  **u;                /* evolving noisy image */
double  **f;                /* noise-free original image */
long    k;                  /* loop variable */
long    nx, ny;             /* image size in x, y direction */ 
double  lambda;             /* contrast parameter */ 
double  sigma;              /* noise scale */ 
double  tau;                /* time step size */ 
long    kmax;               /* number of iterations */
double  max, min;           /* largest, smallest grey value */
double  mean;               /* average grey value */
double  std;                /* standard deviation */
double  mse;                /* mean squared error */
char    comments[1600];     /* string for comments */


printf ("\n");
printf ("ISOTROPIC NONLINEAR DIFFUSION FILTERING, EXPLICIT SCHEME\n\n");
printf ("**************************************************\n\n");
printf ("    Copyright 2022 by Joachim Weickert            \n");
printf ("    Dept. of Mathematics and Computer Science     \n");
printf ("    Saarland University, Saarbruecken, Germany    \n\n");
printf ("    All rights reserved. Unauthorised usage,      \n");
printf ("    copying, hiring, and selling prohibited.      \n\n");
printf ("    Send bug reports to                           \n");
printf ("    weickert@mia.uni-saarland.de                  \n\n");
printf ("**************************************************\n\n");


/* ---- read input images (pgm format P5) ---- */

printf ("noisy input image (pgm):             ");
read_string (in1);
read_pgm_to_double (in1, &nx, &ny, &u);

printf ("noise-free original image (pgm):     ");
read_string (in2);
read_pgm_to_double (in2, &nx, &ny, &f);


/* ---- read parameters ---- */

printf ("contrast parameter lambda (> 0.0):   ");
read_double (&lambda);

printf ("noise scale sigma (> 0.0):           ");
read_double (&sigma);

printf ("time step size (< 0.25)):            ");
read_double (&tau);

printf ("number of iterations:                ");
read_long (&kmax);

printf ("output image (pgm):                  ");
read_string (out);
printf ("\n");


/* ---- analyse initial image ---- */

analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);
mse = MSE (u, f, nx, ny);
printf ("initial image\n");
printf ("minimum:          %8.2lf \n", min);
printf ("maximum:          %8.2lf \n", max);
printf ("mean:             %8.2lf \n", mean);
printf ("standard dev.:    %8.2lf \n", std);
printf ("MSE:              %8.2lf \n\n", mse);


/* ---- analyse processed image ---- */

for (k=1; k<=kmax; k++)
    {
    /* perform one iteration */
    catte (tau, nx, ny, 1.0, 1.0, lambda, sigma, u);

    /* analyse result */
    analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);
    mse = MSE (u, f, nx, ny);
    printf ("iteration:        %8ld   \n", k);
    printf ("minimum:          %8.2lf \n", min);
    printf ("maximum:          %8.2lf \n", max);
    printf ("mean:             %8.2lf \n", mean);
    printf ("standard dev.:    %8.2lf \n", std);
    printf ("MSE:              %8.2lf \n\n", mse);
    }


/* ---- write output image (pgm format P5) ---- */

/* write parameter values in comment string */
comments[0] = '\0';
comment_line (comments, "# isotropic nonlinear diffusion filter\n");
comment_line (comments, "# exponential diffusvity, explicit scheme\n");
comment_line (comments, "# initial image:    %s\n", in1);
comment_line (comments, "# reference image:  %s\n", in2);
comment_line (comments, "# lambda:           %8.2lf\n", lambda);
comment_line (comments, "# sigma:            %8.2lf\n", sigma);
comment_line (comments, "# tau:              %8.4lf\n", tau);
comment_line (comments, "# kmax:             %8ld\n", kmax);
comment_line (comments, "# min:              %8.2lf\n", min);
comment_line (comments, "# max:              %8.2lf\n", max);
comment_line (comments, "# mean:             %8.2lf\n", mean);
comment_line (comments, "# stand. dev.:      %8.2lf\n", std);
comment_line (comments, "# MSE:              %8.2lf\n", mse);

/* write image data */
write_double_to_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

free_double_matrix (u, nx+2, ny+2);

return(0);

}  /* main */
