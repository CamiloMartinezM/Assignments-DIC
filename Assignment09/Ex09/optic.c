#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>

#define PI 3.14159265359


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*               ISOTROPIC FLOW-DRIVEN OPTIC FLOW ESTIMATION                */
/*                                                                          */
/*                 (Copyright by Joachim Weickert, 1/2021)                  */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/*
 features:
 - parabolic approach with Charbonnier diffusivity
 - explicit scheme mit implicit stabilisation
*/

/*--------------------------------------------------------------------------*/

void alloc_double_matrix

(double*** matrix,  /* matrix */
 long   n1,         /* size in direction 1 */
 long   n2)         /* size in direction 2 */

/*
  allocates memory for a double format matrix of size n1 * n2
*/

{
    long i;    /* loop variable */

    *matrix = (double**)malloc(n1 * sizeof(double*));

    if (*matrix == NULL)
    {
        printf("alloc_double_matrix: not enough memory available\n");
        exit(1);
    }

    for (i = 0; i < n1; i++)
    {
        (*matrix)[i] = (double*)malloc(n2 * sizeof(double));
        if ((*matrix)[i] == NULL)
        {
            printf("alloc_double_matrix: not enough memory available\n");
            exit(1);
        }
    }

    return;

}  /* alloc_double_matrix */

/*--------------------------------------------------------------------------*/

void alloc_double_cubix

(double**** cubix,  /* cubix */
 long   n1,         /* size in direction 1 */
 long   n2,         /* size in direction 2 */
 long   n3)         /* size in direction 3 */

/*
  allocates memory for a double format cubix of size n1 * n2 * n3
*/

{
    long i, j;  /* loop variables */

    *cubix = (double***)malloc(n1 * sizeof(double**));

    if (*cubix == NULL)
    {
        printf("alloc_double_cubix: not enough memory available\n");
        exit(1);
    }

    for (i = 0; i < n1; i++)
    {
        (*cubix)[i] = (double**)malloc(n2 * sizeof(double*));
        if ((*cubix)[i] == NULL)
        {
            printf("alloc_double_cubix: not enough memory available\n");
            exit(1);
        }
        for (j = 0; j < n2; j++)
        {
            (*cubix)[i][j] = (double*)malloc(n3 * sizeof(double));
            if ((*cubix)[i][j] == NULL)
            {
                printf("alloc_double_cubix: not enough memory available\n");
                exit(1);
            }
        }
    }

    return;

}  /* alloc_double_cubix */

/*--------------------------------------------------------------------------*/

void free_double_matrix

(double** matrix,   /* matrix */
 long    n1,         /* size in direction 1 */
 long    n2)         /* size in direction 2 */

/*
  frees memory for a double format matrix of size n1 * n2
*/

{
    long i;   /* loop variable */

    for (i = 0; i < n1; i++)
        free(matrix[i]);

    free(matrix);

    return;

}  /* free_double_matrix */

/*--------------------------------------------------------------------------*/

void free_double_cubix

(double*** cubix,   /* cubix */
 long   n1,         /* size in direction 1 */
 long   n2,         /* size in direction 2 */
 long   n3)         /* size in direction 3 */

/*
  frees memory for a double format cubix of size n1 * n2 * n3
*/

{
    long i, j;   /* loop variables */

    for (i = 0; i < n1; i++)
        for (j = 0; j < n2; j++)
            free(cubix[i][j]);

    for (i = 0; i < n1; i++)
        free(cubix[i]);

    free(cubix);

    return;

}  /* free_double_cubix */

/*--------------------------------------------------------------------------*/

void read_string

(char* v)         /* string to be read */

/*
  reads a string v
*/

{
    if (fgets(v, 80, stdin) == NULL)
    {
        printf("could not read string, aborting\n");
        exit(1);
    }

    if (v[strlen(v) - 1] == '\n')
        v[strlen(v) - 1] = 0;

    return;

}  /* read_string */

/*--------------------------------------------------------------------------*/

void read_long

(long* v)         /* value to be read */

/*
  reads a long value v
*/

{
    char   row[80];    /* string for reading data */

    if (fgets(row, 80, stdin) == NULL)
    {
        printf("could not read string, aborting\n");
        exit(1);
    }

    if (row[strlen(row) - 1] == '\n')
        row[strlen(row) - 1] = 0;
    sscanf(row, "%ld", &*v);

    return;

}  /* read_long */

/*--------------------------------------------------------------------------*/

void read_double

(double* v)         /* value to be read */

/*
  reads a double value v
*/

{
    char   row[80];    /* string for reading data */

    if (fgets(row, 80, stdin) == NULL)
    {
        printf("could not read string, aborting\n");
        exit(1);
    }

    if (row[strlen(row) - 1] == '\n')
        row[strlen(row) - 1] = 0;
    sscanf(row, "%lf", &*v);

    return;

}  /* read_double */

/*--------------------------------------------------------------------------*/

void skip_white_space_and_comments

(FILE* inimage)  /* input file */

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
            skip_white_space_and_comments(inimage);
        else
        {
            printf("skip_white_space_and_comments: cannot read file\n");
            exit(1);
        }
    }
    else
        fseek(inimage, -1, SEEK_CUR);

    return;

} /* skip_white_space_and_comments */

/*--------------------------------------------------------------------------*/

void read_pgm_to_double

(const char* file_name,    /* name of pgm file */
 long* nx,           /* image size in x direction, output */
 long* ny,           /* image size in y direction, output */
 double*** u)          /* image, output */

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
    FILE* inimage;     /* input file */

    /* open file */
    inimage = fopen(file_name, "rb");
    if (inimage == NULL)
    {
        printf("read_pgm_to_double: cannot open file '%s'\n", file_name);
        exit(1);
    }

 /* read header */
    if (fgets(row, 80, inimage) == NULL)
    {
        printf("read_pgm_to_double: cannot read file\n");
        exit(1);
    }

 /* image type: P5 */
    if ((row[0] == 'P') && (row[1] == '5'))
    {
    /* P5: grey scale image */
    }
    else
    {
        printf("read_pgm_to_double: unknown image format\n");
        exit(1);
    }

 /* read image size in x direction */
    skip_white_space_and_comments(inimage);
    if (!fscanf(inimage, "%ld", nx))
    {
        printf("read_pgm_to_double: cannot read image size nx\n");
        exit(1);
    }

 /* read image size in x direction */
    skip_white_space_and_comments(inimage);
    if (!fscanf(inimage, "%ld", ny))
    {
        printf("read_pgm_to_double: cannot read image size ny\n");
        exit(1);
    }

 /* read maximum grey value */
    skip_white_space_and_comments(inimage);
    if (!fscanf(inimage, "%ld", &max_value))
    {
        printf("read_pgm_to_double: cannot read maximal value\n");
        exit(1);
    }
    fgetc(inimage);

    /* allocate memory */
    alloc_double_matrix(u, (*nx) + 2, (*ny) + 2);

    /* read image data row by row */
    for (j = 1; j <= (*ny); j++)
        for (i = 1; i <= (*nx); i++)
            (*u)[i][j] = (double)getc(inimage);

       /* close file */
    fclose(inimage);

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
    va_start(arguments, lineformat);

    /* convert format string and arguments to plain text line string */
    vsprintf(line, lineformat, arguments);

    /* add line to total commentary string */
    strncat(comment, line, 80);

    /* add line break if input string does not end with one */
    if (line[strlen(line) - 1] != '\n')
        strncat(comment, "\n", 1);

     /* close argument list */
    va_end(arguments);

    return;

}  /* comment_line */

/*--------------------------------------------------------------------------*/

void write_double_to_pgm_or_ppm

(double*** u,         /* colour image, unchanged */
 long    nc,           /* number of channels */
 long    nx,           /* size in x direction */
 long    ny,           /* size in y direction */
 char* file_name,   /* name of ppm file */
 char* comments)    /* comment string (set 0 for no comments) */

/*
  writes a double format image into a pgm P5 (greyscale) or
  ppm P6 (colour) file;
*/

{
    FILE* outimage;  /* output file */
    long           i, j, m;    /* loop variables */
    double         aux;        /* auxiliary variable */
    unsigned char  byte;       /* for data conversion */

    /* open file */
    outimage = fopen(file_name, "wb");
    if (NULL == outimage)
    {
        printf("Could not open file '%s' for writing, aborting\n", file_name);
        exit(1);
    }

 /* write header */
    if (nc == 1)
        fprintf(outimage, "P5\n");               /* greyscale format */
    else if (nc == 3)
        fprintf(outimage, "P6\n");               /* colour format */
    else
    {
        printf("unsupported number of channels\n");
        exit(0);
    }
    if (comments != 0)
        fputs(comments, outimage);               /* comments */
    fprintf(outimage, "%ld %ld\n", nx, ny);     /* image size */
    fprintf(outimage, "255\n");                 /* maximal value */

    /* write image data */
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++)
            for (m = 0; m <= nc - 1; m++)
            {
                aux = u[m][i][j] + 0.499999;    /* for correct rounding */
                if (aux < 0.0)
                    byte = (unsigned char)(0.0);
                else if (aux > 255.0)
                    byte = (unsigned char)(255.0);
                else
                    byte = (unsigned char)(aux);
                fwrite(&byte, sizeof(unsigned char), 1, outimage);
            }

       /* close file */
    fclose(outimage);

    return;

}  /* write_double_to_pgm_or_ppm */

/*--------------------------------------------------------------------------*/

void dummies_double

(double** u,        /* image */
 long   nx,         /* size in x direction */
 long   ny)         /* size in y direction */

/*
  creates dummy boundaries for a double format image u by mirroring
*/

{
    long i, j;  /* loop variables */

    for (i = 1; i <= nx; i++)
    {
        u[i][0] = u[i][1];
        u[i][ny + 1] = u[i][ny];
    }

    for (j = 0; j <= ny + 1; j++)
    {
        u[0][j] = u[1][j];
        u[nx + 1][j] = u[nx][j];
    }

    return;

}  /* dummies_double */

/*--------------------------------------------------------------------------*/

void analyse_grey_double

(double** u,         /* image, unchanged */
 long    nx,          /* pixel number in x direction */
 long    ny,          /* pixel number in y direction */
 double* min,        /* minimum, output */
 double* max,        /* maximum, output */
 double* mean,       /* mean, output */
 double* std)        /* standard deviation, output */

/*
  computes minimum, maximum, mean, and standard deviation of a greyscale
  image u in double format
*/

{
    long    i, j;       /* loop variables */
    double  help1;      /* auxiliary variable */
    double  help2;      /* auxiliary variable */

    /* compute maximum, minimum, and mean */
    *min = u[1][1];
    *max = u[1][1];
    help1 = 0.0;
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
        {
            if (u[i][j] < *min) *min = u[i][j];
            if (u[i][j] > *max) *max = u[i][j];
            help1 = help1 + u[i][j];
        }
    *mean = help1 / (nx * ny);

    /* compute standard deviation */
    *std = 0.0;
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
        {
            help2 = u[i][j] - *mean;
            *std = *std + help2 * help2;
        }
    *std = sqrt(*std / (nx * ny));

    return;

}  /* analyse_grey_double */

/*--------------------------------------------------------------------------*/

void cartesian_to_polar

(double  x,           /* cartesian coordinate */
 double  y,           /* cartesian coordinate */
 double* r,          /* radius, output */
 double* phi)        /* polar angle, in [0, 2*pi), output */

/*
  transforms cartesian coordinates into polar coordinates
*/

{
    *r = sqrt(x * x + y * y);

    if (x == 0.0)
        if (y >= 0.0)
            *phi = 0.5 * PI;
        else
            *phi = 1.5 * PI;
    else if (x > 0.0)
        if (y >= 0.0)
            *phi = atan(y / x);
        else
            *phi = 2.0 * PI + atan(y / x);
    else
        *phi = PI + atan(y / x);

    return;

}  /* cartesian_to_polar */

/*--------------------------------------------------------------------------*/

void polar_to_RGB

(double  radius,      /* radius, in [0,1] */
 double  phi,         /* angle, unchanged */
 double* R,          /* red   component, in [0,255], output */
 double* G,          /* green component, in [0,255], output */
 double* B)          /* blue  component, in [0,255], output */

/*
  Transforms polar coordinates into RGB values.
  The angle determines the colour, and the radius gives the saturation.
*/

{
    double   alpha, beta;    /* weights for linear interpolation */

    phi = phi / 2.0;

    if ((phi >= 0.0) && (phi < 0.25 * PI))
       /* interpolation between red (0) and blue (0.25 * PI) */
    {
        beta = phi / (0.25 * PI);
        alpha = 1.0 - beta;
        *R = radius * (alpha * 255.0 + beta * 0.0);
        *G = radius * (alpha * 0.0 + beta * 0.0);
        *B = radius * (alpha * 0.0 + beta * 255.0);
    }

    if ((phi >= 0.25 * PI) && (phi < 0.5 * PI))
       /* interpolation between blue (0.25 * PI) and green (0.5 * PI) */
    {
        beta = (phi - 0.25 * PI) / (0.25 * PI);
        alpha = 1.0 - beta;
        *R = radius * (alpha * 0.0 + beta * 0.0);
        *G = radius * (alpha * 0.0 + beta * 255.0);
        *B = radius * (alpha * 255.0 + beta * 0.0);
    }

    if ((phi >= 0.5 * PI) && (phi < 0.75 * PI))
       /* interpolation between green (0.5 * PI) and yellow (0.75 * PI) */
    {
        beta = (phi - 0.5 * PI) / (0.25 * PI);
        alpha = 1.0 - beta;
        *R = radius * (alpha * 0.0 + beta * 255.0);
        *G = radius * (alpha * 255.0 + beta * 255.0);
        *B = radius * (alpha * 0.0 + beta * 0.0);
    }

    if ((phi >= 0.75 * PI) && (phi < PI))
       /* interpolation between yellow (0.75 * PI) and red (PI) */
    {
        beta = (phi - 0.75 * PI) / (0.25 * PI);
        alpha = 1.0 - beta;
        *R = radius * (alpha * 255.0 + beta * 255.0);
        *G = radius * (alpha * 255.0 + beta * 0.0);
        *B = radius * (alpha * 0.0 + beta * 0.0);
    }

    return;

}  /* polar_to_RGB */

/*--------------------------------------------------------------------------*/

void flow

(double   tau,         /* time step size, 0 < tau <= 0.25 */
 long     nx,          /* image dimension in x direction */
 long     ny,          /* image dimension in y direction */
 double   hx,          /* pixel size in x direction */
 double   hy,          /* pixel size in y direction */
 double** fx,        /* x derivative of image */
 double** fy,        /* y derivative of image */
 double** fz,        /* time derivative of image */
 double   alpha,       /* smoothness weight */
 double   lambda,      /* contrast parameter */
 double** u,         /* x component of optic flow */
 double** v)         /* y component of optic flow */

/*
  Optic flow iteration by regarding it as a coupled system of
  two nonlinear diffusion-reaction equations.
  Isotropic nonlinear diffusion with modified explicit discretisation.
*/

{
    long    i, j;             /* loop variables */
    double** uold, ** vold;   /* u and v at old iteration level */
    double** g;              /* diffusivity */
    double  ux, vx;           /* central differences of u, v in x direction */
    double  uy, vy;           /* central differences of u, v in y direction */
    double  aux, rxx, ryy;    /* time savers */
    double  inv_two_hx;       /* time saver */
    double  inv_two_hy;       /* time saver */
    double  inv_lambda_sqr;   /* time saver */
    double  grad_sqr;         /* |grad(u)|^2 + |grad(v)|^2 */


    /* ---- allocate memory ---- */

    alloc_double_matrix(&uold, nx + 2, ny + 2);
    alloc_double_matrix(&vold, nx + 2, ny + 2);
    alloc_double_matrix(&g, nx + 2, ny + 2);


    /* ---- copy u, v into uold, vold ---- */

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
        {
            uold[i][j] = u[i][j];
            vold[i][j] = v[i][j];
        }


   /* ---- compute diffusivity ---- */

   /* create reflecting dummy boundaries for uold, vold */
    dummies_double(uold, nx, ny);
    dummies_double(vold, nx, ny);

    /* compute time savers */
    inv_two_hx = 1.0 / (2.0 * hx);
    inv_two_hy = 1.0 / (2.0 * hy);
    inv_lambda_sqr = 1.0 / (lambda * lambda);

    /* compute |grad u|^2 + |grad v|^2 and diffusivity */
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
        {
        /* compute derivatives of u and v */
            ux = (uold[i + 1][j] - uold[i - 1][j]) * inv_two_hx;
            uy = (uold[i][j + 1] - uold[i][j - 1]) * inv_two_hy;
            vx = (vold[i + 1][j] - vold[i - 1][j]) * inv_two_hx;
            vy = (vold[i][j + 1] - vold[i][j - 1]) * inv_two_hy;

            /* compute |grad u|^2 + |grad v|^2 */
            grad_sqr = ux * ux + uy * uy + vx * vx + vy * vy;

            /* compute Charbonnier diffusivity */
            g[i][j] = 1.0 / sqrt(1.0 + grad_sqr * inv_lambda_sqr);
        }


   /* ---- compute modified explicit step ---- */

   /* compute time savers */
    aux = tau / alpha;
    rxx = tau / (2.0 * hx * hx);
    ryy = tau / (2.0 * hy * hy);

    /* create reflecting dummy boundaries for g */
    dummies_double(g, nx, ny);

    /* compute modified explicit step */
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
        {
            {
            /* SUPPLEMENT CODE */

            /* u[i][j] = ... */
                u[i][j] = 1 / (1 + aux * fx[i][j] * fx[i][j]) * (uold[i][j]
                    + rxx * ((g[i + 1][j] + g[i][j]) * (uold[i + 1][j] - uold[i][j])
                        + (g[i - 1][j] + g[i][j]) * (uold[i - 1][j] - uold[i][j]))
                    + ryy * ((g[i][j + 1] + g[i][j]) * (uold[i][j + 1] - uold[i][j])
                        + (g[i][j - 1] + g[i][j]) * (uold[i][j - 1] - uold[i][j]))
                    - aux * fx[i][j] * (fy[i][j] * vold[i][j] + fz[i][j]));

                 /* v[i][j] = ... */
                v[i][j] = 1 / (1 + aux * fy[i][j] * fy[i][j]) * (vold[i][j]
                    + rxx * ((g[i + 1][j] + g[i][j]) * (vold[i + 1][j] - vold[i][j])
                        + (g[i - 1][j] + g[i][j]) * (vold[i - 1][j] - vold[i][j]))
                    + ryy * ((g[i][j + 1] + g[i][j]) * (vold[i][j + 1] - vold[i][j])
                        + (g[i][j - 1] + g[i][j]) * (vold[i][j - 1] - vold[i][j]))
                    - aux * fy[i][j] * (fx[i][j] * uold[i][j] + fz[i][j]));
            }
        }


   /* ---- free memory ---- */

    free_double_matrix(uold, nx + 2, ny + 2);
    free_double_matrix(vold, nx + 2, ny + 2);
    free_double_matrix(g, nx + 2, ny + 2);

    return;

}  /* flow */

/*--------------------------------------------------------------------------*/

int main()

{
    char    in1[80], in2[80];     /* input file names */
    char    out[80];              /* output file name */
    double** f1, ** f2;           /* images */
    double** fx, ** fy, ** fz;     /* image derivatives */
    double** u, ** v;             /* optic flow components */
    double** w;                  /* optic flow magnitude */
    double*** rgb;               /* colour representation of (u,v) */
    long    i, j, k;              /* loop variables */
    long    kmax;                 /* number of iterations */
    long    nx, ny;               /* image size in x, y direction */
    double  hx, hy;               /* pixel sizes */
    double  tau;                  /* time step size */
    double  alpha;                /* smoothness weight */
    double  lambda;               /* contrast parameter */
    double  min, max;             /* smallest, largest value */
    double  mean;                 /* average value */
    double  std;                  /* standard deviation */
    double  radius, phi;          /* polar coordinates */
    char    comments[1600];       /* string for comments */


    printf("\n");
    printf("ISOTROPIC FLOW-DRIVEN OPTIC FLOW COMPUTATION\n");
    printf("SPATIAL MODEL WITH CHARBONNIER REGULARISATION\n");
    printf("PARABOLIC APPROACH WITH MODIFIED EXPLICIT SCHEME\n\n");
    printf("************************************************************\n\n");
    printf("    Copyright 2021 by Joachim Weickert                      \n");
    printf("    Dept. of Mathematics and Computer Science               \n");
    printf("    Saarland University, Saarbruecken, Germany              \n\n");
    printf("    All rights reserved. Unauthorized usage,                \n");
    printf("    copying, hiring, and selling prohibited.                \n\n");
    printf("    Send bug reports to                                     \n");
    printf("    weickert@mia.uni-saarland.de                            \n\n");
    printf("************************************************************\n\n");


    /* ---- read input images f1 and f2 (pgm format P5) ---- */

    printf("input image 1 (pgm):                    ");
    read_string(in1);
    read_pgm_to_double(in1, &nx, &ny, &f1);

    printf("input image 2 (pgm):                    ");
    read_string(in2);
    read_pgm_to_double(in2, &nx, &ny, &f2);


    /* ---- read parameters ---- */

    printf("smoothness weight alpha (>0.0):         ");
    read_double(&alpha);

    printf("contrast parameter lambda (>0.0):       ");
    read_double(&lambda);

    printf("time step size (<=0.25):                ");
    read_double(&tau);

    printf("number of iterations (>0):              ");
    read_long(&kmax);

    printf("colour-coded output flow field (ppm):   ");
    read_string(out);

    printf("\n");


    /* ---- initialisations ---- */

    /* allocate memory */
    alloc_double_matrix(&fx, nx + 2, ny + 2);
    alloc_double_matrix(&fy, nx + 2, ny + 2);
    alloc_double_matrix(&fz, nx + 2, ny + 2);
    alloc_double_matrix(&u, nx + 2, ny + 2);
    alloc_double_matrix(&v, nx + 2, ny + 2);
    alloc_double_matrix(&w, nx + 2, ny + 2);
    alloc_double_cubix(&rgb, 3, nx + 2, ny + 2);

    /* set grid size to 1 */
    hx = 1.0;
    hy = 1.0;

    /* compute image derivatives fx, fy and fz */
    dummies_double(f1, nx, ny);
    dummies_double(f2, nx, ny);
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
        {
            fx[i][j] = (f1[i + 1][j] - f1[i - 1][j] + f2[i + 1][j] - f2[i - 1][j])
                / (4.0 * hx);
            fy[i][j] = (f1[i][j + 1] - f1[i][j - 1] + f2[i][j + 1] - f2[i][j - 1])
                / (4.0 * hy);
            fz[i][j] = f2[i][j] - f1[i][j];   /* frame distance 1 assumed */
        }

   /* initialise (u,v) with 0 */
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
        {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
        }


   /* ---- iterative optic flow computation ---- */

    for (k = 1; k <= kmax; k++)
    {
    /* perform one iteration */
        flow(tau, nx, ny, hx, hy, fx, fy, fz, alpha, lambda, u, v);

        /* compute flow magnitude */
        for (i = 1; i <= nx; i++)
            for (j = 1; j <= ny; j++)
                w[i][j] = sqrt(u[i][j] * u[i][j] + v[i][j] * v[i][j]);

           /* analyse flow magnitude */
        analyse_grey_double(w, nx, ny, &min, &max, &mean, &std);
        printf("\n");
        printf("iteration number:         %8ld \n", k);
        printf("maximal flow magnitude:   %8.3lf\n", max);
        printf("average flow magnitude:   %8.3lf\n", mean);
        printf("standard deviation:       %8.3lf\n", std);
    }


/* ---- convert optic flow field (u,v) to RGB representation ---- */

    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++)
        {
            cartesian_to_polar(u[i][j], v[i][j], &radius, &phi);
            polar_to_RGB(radius / max, phi,
                          &rgb[0][i][j], &rgb[1][i][j], &rgb[2][i][j]);
        }


   /* ---- write optic flow in colour representation (ppm format P6) ---- */

   /* write parameter values in comment string */
    comments[0] = '\0';
    comment_line(comments, "# isotropic flow-driven optic flow computation\n");
    comment_line(comments, "# spatial approach with Charbonnier reguliser\n");
    comment_line(comments, "# modified explicit scheme\n");
    comment_line(comments, "# frame 1:           %s\n", in1);
    comment_line(comments, "# frame 2:           %s\n", in2);
    comment_line(comments, "# alpha:             %8.2lf\n", alpha);
    comment_line(comments, "# lambda:            %8.2lf\n", lambda);
    comment_line(comments, "# tau:               %8.2lf\n", tau);
    comment_line(comments, "# iterations:        %8ld\n", kmax);
    comment_line(comments, "# max. magnitude:    %8.3lf\n", max);
    comment_line(comments, "# aver. magnitude:   %8.3lf\n", mean);
    comment_line(comments, "# stand. deviation:  %8.3lf\n", std);

    /* write flow data */
    write_double_to_pgm_or_ppm(rgb, 3, nx, ny, out, comments);
    printf("\n");
    printf("output flow field %s successfully written\n\n", out);


    /* ---- free memory ---- */

    free_double_matrix(f1, nx + 2, ny + 2);
    free_double_matrix(f2, nx + 2, ny + 2);
    free_double_matrix(fx, nx + 2, ny + 2);
    free_double_matrix(fy, nx + 2, ny + 2);
    free_double_matrix(fz, nx + 2, ny + 2);
    free_double_matrix(u, nx + 2, ny + 2);
    free_double_matrix(v, nx + 2, ny + 2);
    free_double_matrix(w, nx + 2, ny + 2);
    free_double_cubix(rgb, 3, nx + 2, ny + 2);

    return(0);

}  /* main */