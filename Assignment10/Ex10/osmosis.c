#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                            OSMOSIS FILTERING                             */
/*                                                                          */
/*                  (Copyright Joachim Weickert, 1/2021)                    */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/* 
  explicit scheme for greyscale images;
  reconstructs an image from its canonical drift vectors 
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

void canonical_drift_vectors 

     (double   **v,     /* guidance image, unchanged */
      long     nx,      /* image dimension in x direction */ 
      long     ny,      /* image dimension in y direction */ 
      double   hx,      /* pixel size in x direction */
      double   hy,      /* pixel size in y direction */
      double   **d1,    /* drift vector, x component in [i+1/2,j], output */
      double   **d2)    /* drift vector, y component in [i,j+1/2], output */


/*
  computes the canonical drift vector field that allows to reconstruct the 
  guidance image up to a multiplicative constant
*/

{
long    i, j;             /* loop variables */


/* ---- dummy boundaries for v ---- */

dummies_double (v, nx, ny);

printf("The value of myFloatValue is: %ld\n", nx);
printf("The value of myFloatValue is: %ld\n", ny);

/* ---- initialise drift vector field with 0 ---- */

for (i=0; i<=nx+1; i++)
 for (j=0; j<=ny+1; j++)
     d1[i][j] = d2[i][j] = 0.0;

for (i=0; i<=nx+1; i++)
 for (j=0; j<=ny+1; j++)
     {
/* ---- compute x component of canonical drift vector field ---- */

/* index [i,j] refers to intergrid location [i+1/2,j] */
      d1[i][j] = (d1[i][j] + d1[i+1][j])/(2*hx);

/* SUPPLEMENT CODE */


/* ---- compute y component of canonical drift vector field ---- */

/* index [i,j] refers to intergrid location [i,j+1/2] */
       d2[i][j] = (d2[i][j] + d2[i][j+1])/(2*hy);

/* SUPPLEMENT CODE */
     }


/* ---- modification at the shadow boundaries between i=128 and i=129 ---- */

/* SUPPLEMENT CODE */

for (i=128; i<=129; i++)
 for (j=0; j<=ny+1; j++)
 {
    d1[i][j] = 0;
    d2[i][j] = 0;
 }
 
return;

}  /* canonical_drift_vectors */

/*--------------------------------------------------------------------------*/

void osmosis_weights 

     (double   tau,     /* time step size, 0 < tau < 0.125 */
      long     nx,      /* image dimension in x direction */
      long     ny,      /* image dimension in y direction */
      double   hx,      /* pixel size in x direction */
      double   hy,      /* pixel size in y direction */
      double   **d1,    /* drift vector, x component in [i+1/2,j], unchanged */
      double   **d2,    /* drift vector, y component in [i,j+1/2], unchanged */
      double   **woo,   /* osmosis weight for pixel [i][j],   output */
      double   **wpo,   /* osmosis weight for pixel [i+1][j], output */
      double   **wmo,   /* osmosis weight for pixel [i-1][j], output */
      double   **wop,   /* osmosis weight for pixel [i][j+1], output */
      double   **wom)   /* osmosis weight for pixel [i][j-1], output */

/*
  computes the weights for osmosis filtering
*/

{
long    i, j;             /* loop variables */
double  rx, rxx;          /* time savers */
double  ry, ryy;          /* time savers */


/* ---- initialise all osmosis weights ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     woo[i][j] = 1.0;
     wpo[i][j] = wmo[i][j] = wop[i][j] = wom[i][j] = 0.0; 
     }


/* ---- specify them from the drift vector field ---- */

/* compute time savers */
rx  = tau / (2.0 * hx);
ry  = tau / (2.0 * hy);
rxx = tau / (hx * hx);
ryy = tau / (hy * hy);

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* osmosis weight for pixel [i][j] */
     woo[i][j] = 1.0 - 2.0 * (rxx + ryy)
                 + rx * (d1[i-1][j] - d1[i][j])
                 + ry * (d2[i][j-1] - d2[i][j]);

     /* osmosis weight for pixel [i+1][j] */
     wpo[i][j] = rxx - rx * d1[i][j];

     /* osmosis weight for pixel [i-1][j] */
     wmo[i][j] = rxx + rx * d1[i-1][j];

     /* osmosis weight for pixel [i][j+1] */
     wop[i][j] = ryy - ry * d2[i][j];

     /* osmosis weight for pixel [i][j-1] */
     wom[i][j] = ryy + ry * d2[i][j-1];
     }

return;  

}  /* weights */

/*--------------------------------------------------------------------------*/

void osmosis 

     (long     nx,      /* image dimension in x direction */ 
      long     ny,      /* image dimension in y direction */ 
      double   **woo,   /* osmosis weight for pixel [i][j],   output */
      double   **wpo,   /* osmosis weight for pixel [i+1][j], output */
      double   **wmo,   /* osmosis weight for pixel [i-1][j], output */
      double   **wop,   /* osmosis weight for pixel [i][j+1], output */
      double   **wom,   /* osmosis weight for pixel [i][j-1], output */
      double   **u)     /* input: original image;  output: filtered */

/* 
  performs one explicit step of an osmosis filter 
*/

{
long    i, j;             /* loop variables */
double  **f;              /* work copy of u */
      

/* ---- allocate memory ---- */

alloc_double_matrix (&f, nx+2, ny+2);


/* ---- copy u into f ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     f[i][j] = u[i][j];


/* ---- dummy boundaries for f ---- */

dummies_double (f, nx, ny);


/* ---- compute explicit osmosis of u ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     u[i][j] = woo[i][j] * f[i][j]  
             + wpo[i][j] * f[i+1][j] 
             + wmo[i][j] * f[i-1][j] 
             + wop[i][j] * f[i][j+1] 
             + wom[i][j] * f[i][j-1];


/* ---- free memory ---- */

free_double_matrix (f, nx+2, ny+2);

return;

}  /* osmosis */

/*--------------------------------------------------------------------------*/

int main ()

{
char    in1[80], in2[80];     /* for reading data */
char    out[80];              /* for writing data */
double  **u;                  /* evolving image */
double  **v;                  /* guidance image */
double  **d1;                 /* drift vector field, x component */
double  **d2;                 /* drift vector field, y component */
double  **woo;                /* osmosis weight for pixel [i][j] */
double  **wpo;                /* osmosis weight for pixel [i+1][j] */
double  **wmo;                /* osmosis weight for pixel [i-1][j] */
double  **wop;                /* osmosis weight for pixel [i][j+1] */
double  **wom;                /* osmosis weight for pixel [i][j-1] */
long    i, j, k;              /* loop variables */
long    kmax;                 /* largest iteration number */
long    nx, ny;               /* image size in x, y direction */
double  tau;                  /* time step size */
double  offset;               /* greyscale offset */
double  max, min;             /* largest, smallest grey value */
double  mean;                 /* average grey value */
double  std;                  /* standard deviation */
char    comments[1600];       /* string for comments */


printf("\n");
printf("OSMOSIS FILTERING, EXPLICIT SCHEME\n");
printf("RECONSTRUCTION FROM THE CANONICAL DRIFT VECTOR FIELD\n\n");
printf("*************************************************\n\n");
printf("    Copyright 2021 by Joachim Weickert           \n");
printf("    Dept. of Mathematics and Computer Science    \n");
printf("    Saarland University, Germany                 \n\n");
printf("    All rights reserved. Unauthorized usage,     \n");
printf("    copying, hiring, and selling prohibited.     \n\n");
printf("    Send bug reports to                          \n");
printf("    weickert@mia.uni-saarland.de                 \n\n");
printf("*************************************************\n\n");


/* ---- read initial image (pgm format P5) ---- */

printf ("initial image f:                  ");
read_string (in1);
read_pgm_to_double (in1, &nx, &ny, &u);


/* ---- read guidance image (pgm format P5) ---- */

printf ("guidance image v:                 ");
read_string (in2);
read_pgm_to_double (in2, &nx, &ny, &v);


/* ---- read other parameters ---- */

printf ("greyscale offset (>0.0):          ");
read_double (&offset);

printf ("time step size (<0.125):          ");
read_double (&tau);

printf ("number of iterations:             ");
read_long (&kmax);

printf ("output image:                     ");
read_string (out);
printf("\n");


/* ---- allocate memory ---- */

alloc_double_matrix (&d1,  nx+2, ny+2);
alloc_double_matrix (&d2,  nx+2, ny+2);
alloc_double_matrix (&woo, nx+2, ny+2);
alloc_double_matrix (&wpo, nx+2, ny+2);
alloc_double_matrix (&wmo, nx+2, ny+2);
alloc_double_matrix (&wop, nx+2, ny+2);
alloc_double_matrix (&wom, nx+2, ny+2);


/* ---- add offset in order to make data positive ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     u[i][j] = u[i][j] + offset;
     v[i][j] = v[i][j] + offset;
     }


/* ---- analyse initial image ---- */

analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);
printf ("initial image\n");
printf ("minimum with offset:  %8.2lf \n", min);
printf ("maximum with offset:  %8.2lf \n", max);
printf ("mean with offset:     %8.2lf \n", mean);
printf ("standard deviation:   %8.2lf \n\n", std);


/* ---- perform osmosis filtering ---- */ 

/* compute canonical drift vectors of the guidance image */
canonical_drift_vectors (v, nx, ny, 1.0, 1.0, d1, d2);

/* compute resulting osmosis weights */
osmosis_weights (tau, nx, ny, 1.0, 1.0, d1, d2, woo, wpo, wmo, wop, wom);

/* perform kmax osmosis iterations */
for (k=1; k<=kmax; k++)
    {
    /* perform one iteration */
    osmosis (nx, ny, woo, wpo, wmo, wop, wom, u);

    /* check minimum, maximum, mean, standard deviation */
    analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);
    printf ("iteration number:     %8ld \n", k);
    printf ("minimum with offset:  %8.2lf \n", min);
    printf ("maximum with offset:  %8.2lf \n", max);
    printf ("mean with offset:     %8.2lf \n", mean);
    printf ("standard deviation:   %8.2lf \n\n", std);
    } /* for */


/* ---- subtract offset ----*/

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     u[i][j] = u[i][j] - offset;


/* ---- write output image (pgm format P5) ---- */

/* write parameter values in comment string */
comments[0] = '\0';
comment_line (comments, "# linear osmosis, explicit scheme\n"); 
comment_line (comments, "# initial image:  %s\n", in1);
comment_line (comments, "# guidance image: %s\n", in2);
comment_line (comments, "# tau:                 %8.4lf\n", tau);
comment_line (comments, "# iterations:          %8ld\n",   kmax);
comment_line (comments, "# offset:              %8.2lf\n", offset);
comment_line (comments, "# minimum with offset: %8.2lf\n", min);
comment_line (comments, "# maximum with offset: %8.2lf\n", max);
comment_line (comments, "# mean with offset:    %8.2lf\n", mean);
comment_line (comments, "# standard deviation:  %8.2lf\n", std);

/* write image data */
write_double_to_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory ---- */

free_double_matrix (u,   nx+2, ny+2);
free_double_matrix (v,   nx+2, ny+2);
free_double_matrix (d1,  nx+2, ny+2);
free_double_matrix (d2,  nx+2, ny+2);
free_double_matrix (woo, nx+2, ny+2);
free_double_matrix (wpo, nx+2, ny+2);
free_double_matrix (wmo, nx+2, ny+2);
free_double_matrix (wop, nx+2, ny+2);
free_double_matrix (wom, nx+2, ny+2);

return(0);

}  /* main */

