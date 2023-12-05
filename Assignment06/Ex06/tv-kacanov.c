#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                 TV RESTORATION WITH THE KACANOV METHOD                   */
/*                                                                          */
/*                (Copyright by Joachim Weickert, 11/2021)                  */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/*
  features:
  - TV diffusivity with epsilon regularisation;
  - uses Gauss-Seidel as linear system solver
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
*std = sqrt(*std / (nx * ny));

return;

}  /* analyse_grey_double */

/*--------------------------------------------------------------------------*/

void diffusivity 

     (long     nx,        /* image dimension in x direction */ 
      long     ny,        /* image dimension in y direction */ 
      double   hx,        /* pixel size in x direction */
      double   hy,        /* pixel size in y direction */
      double   eps,       /* parameter epsilon */
      double   **u,       /* original image, unchanged */
      double   **g)       /* diffusivity, output */


/* 
  computes epsilon-regularised TV diffusivities 
*/

{
long    i, j;                 /* loop variables */
double  ux, uy;               /* derivatives */
double  grad_sqr;             /* |grad(u)|^2 */
double  eps_sqr;              /* epsilon * epsilon */
double  aux1, aux2;           /* time savers */

dummies_double (u, nx, ny);

aux1 = 1.0 / (2.0 * hx);
aux2 = 1.0 / (2.0 * hy);
eps_sqr = eps * eps;

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     ux = aux1 * (u[i+1][j] - u[i-1][j]);
     uy = aux2 * (u[i][j+1] - u[i][j-1]);
     grad_sqr = ux * ux + uy * uy;
     g[i][j] = 1.0 / sqrt (eps_sqr + grad_sqr);
     }

return;

}  /* diffusivity */

/*--------------------------------------------------------------------------*/

void gauss_seidel

    (double  alpha,  /* regularisation parameter */
     long    nx,     /* image dimension in x direction */ 
     long    ny,     /* image dimension in y direction */ 
     double  hx,     /* pixel size in x direction */
     double  hy,     /* pixel size in y direction */
     double  **f,    /* original image, unchanged on exit */
     double  **g,    /* diffusivity, with boundary modifications on exit */
     double  **u)    /* old and new solution */


/*
  performs one Gauss-Seidel iteration
*/

{
long    i, j;                  /* loop variables */
double  wxp, wyp, wxm, wym;    /* weights */
double  r1, r2;                /* time savers */

dummies_double (u, nx, ny);
dummies_double (g, nx, ny);

r1 = alpha / (2.0 * hx);   
r2 = alpha / (2.0 * hy);   

for (i=1; i<=nx; i++) 
 for (j=1; j<=ny; j++) 
     { 
     wxp = g[i+1][j] + g[i][j]; 
     wxm = g[i-1][j] + g[i][j]; 
     wyp = g[i][j+1] + g[i][j]; 
     wym = g[i][j-1] + g[i][j]; 

     u[i][j] = (f[i][j] + r1 * (wxp * u[i+1][j] + wxm * u[i-1][j])
                        + r2 * (wyp * u[i][j+1] + wym * u[i][j-1]))
               / (1.0 + r1 * (wxp + wxm) + r2 * (wyp + wym)); 
     } /* for i,j */ 

return;

}  /* gauss_seidel */

/*--------------------------------------------------------------------------*/

int main ()

{
char    in[80];               /* for reading data */
char    out[80];              /* for writing data */
double  **u;                  /* evolving image */
double  **g;                  /* diffusivity */
double  **f;                  /* original image */
long    i, j, k, m;           /* loop variables */
long    nx, ny;               /* image size in x, y direction */
double  hx, hy;               /* pixel size in x, y direction */
double  alpha;                /* regularisation parameter */
double  eps;                  /* parameter epsilon */
long    maxfp;                /* number of fixed point iterations */
long    maxgs;                /* mumber of Gauss-Seidel iterations */
double  max, min;             /* largest, smallest grey value */
double  mean;                 /* average grey value */
double  std;                  /* standard deviation */
char    comments[1600];       /* string for comments */


printf("\n");
printf("\n");
printf("TV IMAGE RESTORATION WITH THE KACANOV METHOD\n\n");
printf("*********************************************************\n\n");
printf("    Copyright 2021 by Joachim Weickert                   \n");
printf("    Faculty of Mathematics and Computer Science          \n");
printf("    Saarland University, Saarbruecken, Germany           \n\n");
printf("    All rights reserved. Unauthorized usage, copying,    \n");
printf("    hiring, and selling prohibited.                      \n\n");
printf("    Send bug reports to                                  \n");
printf("    weickert@mia.uni-saarland.de                         \n\n");
printf("*********************************************************\n\n");
printf("- outer iterations:  fixed point\n");
printf("- inner iterations:  Gauss-Seidel\n\n");


/* ---- read input image (pgm format P5) ---- */

/* read image name */
printf ("input image (pgm):                       ");
read_string (in);
read_pgm_to_double (in, &nx, &ny, &f);


/* ---- read other parameters ---- */

printf ("regularisation parameter alpha (>=0.0):  ");
read_double (&alpha);

printf ("parameter epsilon (>0.0):                ");
read_double (&eps);

printf ("number of fixed point iterations:        ");
read_long (&maxfp);

printf ("number of Gauss-Seidel iterations:       ");
read_long (&maxgs);

printf ("output image (pgm):                      ");
read_string (out);

printf("\n");


/* ---- allocate memory ---- */

alloc_double_matrix (&u,  nx+2, ny+2);
alloc_double_matrix (&g, nx+2, ny+2);


/* ---- initialisations ---- */

/* unit grid sizes */
hx = 1.0;
hy = 1.0;

/* copy f to u */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     u[i][j] = f[i][j];


/* ---- analyse initial image ---- */

analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);
printf ("initial image \n");
printf ("minimum:                 %8.2lf \n", min);
printf ("maximum:                 %8.2lf \n", max);
printf ("mean:                    %8.2lf \n", mean);
printf ("standard deviation:      %8.2lf \n\n", std);


/* ---- process image ---- */

for (m=1; m<=maxfp; m++)
  
    /* ---- outer iterations: fixed point ---- */ 

    {     
    /* update diffusivity */
    diffusivity (nx, ny, hx, hy, eps, u, g);

    /* inner iterations: Gauss-Seidel */
    for (k=1; k<=maxgs; k++)
        gauss_seidel (alpha, nx, ny, hx, hy, f, g, u);

    /* analyse image */
    analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);
    printf ("fixed point iteration:   %8ld\n", m);
    printf ("minimum:                 %8.2lf \n", min);
    printf ("maximum:                 %8.2lf \n", max);
    printf ("mean:                    %8.2lf \n", mean);
    printf ("standard deviation:      %8.2lf \n\n", std);
    } /* for m */


/* ---- write output image (pgm format P5) ---- */

/* write parameter values in comment string */
comments[0] = '\0';
comment_line (comments, "# TV image restoration with Kacanov method\n");
comment_line (comments, "# outer iterations:  fixed point (FP)\n");
comment_line (comments, "# inner iterations:  Gauss-Seidel (GS)\n");
comment_line (comments, "# initial image:     %s\n", in);
comment_line (comments, "# alpha:             %8.2lf\n", alpha);
comment_line (comments, "# epsilon:           %8.4lf\n", eps);
comment_line (comments, "# FP iterations:     %8ld\n", maxfp);
comment_line (comments, "# GS iterations:     %8ld\n", maxgs);
comment_line (comments, "# minimum:           %8.2lf\n", min);
comment_line (comments, "# maximum:           %8.2lf\n", max);
comment_line (comments, "# mean:              %8.2lf\n", mean);
comment_line (comments, "# standard dev.:     %8.2lf\n", std);

/* write image data */
write_double_to_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory ---- */

free_double_matrix (f,  nx+2, ny+2);
free_double_matrix (u,  nx+2, ny+2);
free_double_matrix (g, nx+2, ny+2);

return(0);

}  /* main */
