#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*     AMBROSIO-TORTORELLI APPROXIMATION OF THE MUMFORD-SHAH FUNCTIONAL     */
/*                                                                          */
/*                  (Copyright Joachim Weickert, 12/2020)                   */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/*
explicit scheme with implicitly stabilised reaction term
*/

/*--------------------------------------------------------------------------*/

void alloc_double_matrix

(double*** matrix, /* matrix */
 long n1,          /* size in direction 1 */
 long n2)          /* size in direction 2 */

/*
allocates memory for a double format matrix of size n1 * n2
*/

{
    long i; /* loop variable */

    *matrix = (double**)malloc(n1 * sizeof(double*));

    if (*matrix == NULL) {
        printf("alloc_double_matrix: not enough memory available\n");
        exit(1);
    }

    for (i = 0; i < n1; i++) {
        (*matrix)[i] = (double*)malloc(n2 * sizeof(double));
        if ((*matrix)[i] == NULL) {
            printf("alloc_double_matrix: not enough memory available\n");
            exit(1);
        }
    }

    return;

} /* alloc_double_matrix */

/*--------------------------------------------------------------------------*/

void free_double_matrix

(double** matrix, /* matrix */
 long n1,         /* size in direction 1 */
 long n2)         /* size in direction 2 */

/*
frees memory for a double format matrix of size n1 * n2
*/

{
    long i; /* loop variable */

    for (i = 0; i < n1; i++)
        free(matrix[i]);

    free(matrix);

    return;

} /* free_double_matrix */

/*--------------------------------------------------------------------------*/

void read_string

(char* v) /* string to be read */

/*
reads a string v
*/

{
    if (fgets(v, 80, stdin) == NULL) {
        printf("could not read string, aborting\n");
        exit(1);
    }

    if (v[strlen(v) - 1] == '\n')
        v[strlen(v) - 1] = 0;

    return;

} /* read_string */

/*--------------------------------------------------------------------------*/

void read_long

(long* v) /* value to be read */

/*
reads a long value v
*/

{
    char row[80]; /* string for reading data */

    if (fgets(row, 80, stdin) == NULL) {
        printf("could not read string, aborting\n");
        exit(1);
    }

    if (row[strlen(row) - 1] == '\n')
        row[strlen(row) - 1] = 0;
    sscanf(row, "%ld", &*v);

    return;

} /* read_long */

/*--------------------------------------------------------------------------*/

void read_double

(double* v) /* value to be read */

/*
reads a double value v
*/

{
    char row[80]; /* string for reading data */

    if (fgets(row, 80, stdin) == NULL) {
        printf("could not read string, aborting\n");
        exit(1);
    }

    if (row[strlen(row) - 1] == '\n')
        row[strlen(row) - 1] = 0;
    sscanf(row, "%lf", &*v);

    return;

} /* read_double */

/*--------------------------------------------------------------------------*/

void skip_white_space_and_comments

(FILE* inimage) /* input file */

/*
skips over white space and comments while reading the file
*/

{

    int ch = 0;   /* holds a character */
    char row[80]; /* for reading data */

    /* skip spaces */
    while (((ch = fgetc(inimage)) != EOF) && isspace(ch))
        ;

      /* skip comments */
    if (ch == '#') {
        if (fgets(row, sizeof(row), inimage))
            skip_white_space_and_comments(inimage);
        else {
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

(const char* file_name, /* name of pgm file */
 long* nx,              /* image size in x direction, output */
 long* ny,              /* image size in y direction, output */
 double*** u)           /* image, output */

/*
reads a greyscale image that has been encoded in pgm format P5 to
an image u in double format;
allocates memory for the image u;
adds boundary layers of size 1 such that
- the relevant image pixels in x direction use the indices 1,...,nx
- the relevant image pixels in y direction use the indices 1,...,ny
*/

{
    char row[80];   /* for reading data */
    long i, j;      /* image indices */
    long max_value; /* maximum color value */
    FILE* inimage;  /* input file */

    /* open file */
    inimage = fopen(file_name, "rb");
    if (inimage == NULL) {
        printf("read_pgm_to_double: cannot open file '%s'\n", file_name);
        exit(1);
    }

    /* read header */
    if (fgets(row, 80, inimage) == NULL) {
        printf("read_pgm_to_double: cannot read file\n");
        exit(1);
    }

    /* image type: P5 */
    if ((row[0] == 'P') && (row[1] == '5')) {
      /* P5: grey scale image */
    }
    else {
        printf("read_pgm_to_double: unknown image format\n");
        exit(1);
    }

    /* read image size in x direction */
    skip_white_space_and_comments(inimage);
    if (!fscanf(inimage, "%ld", nx)) {
        printf("read_pgm_to_double: cannot read image size nx\n");
        exit(1);
    }

    /* read image size in x direction */
    skip_white_space_and_comments(inimage);
    if (!fscanf(inimage, "%ld", ny)) {
        printf("read_pgm_to_double: cannot read image size ny\n");
        exit(1);
    }

    /* read maximum grey value */
    skip_white_space_and_comments(inimage);
    if (!fscanf(inimage, "%ld", &max_value)) {
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

} /* read_pgm_to_double */

/*--------------------------------------------------------------------------*/

void comment_line

(char* comment,    /* comment string (output) */
 char* lineformat, /* format string for comment line */
 ...)              /* optional arguments */

/*
Adds a line to the comment string comment. The string line can contain
plain text and format characters that are compatible with sprintf.
Example call:
print_comment_line(comment, "Text %lf %ld", double_var, long_var).
If no line break is supplied at the end of the input string, it is
added automatically.
*/

{
    char line[80];
    va_list arguments;

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

} /* comment_line */

/*--------------------------------------------------------------------------*/

void write_double_to_pgm

(double** u,      /* image, unchanged */
 long nx,         /* image size in x direction */
 long ny,         /* image size in y direction */
 char* file_name, /* name of pgm file */
 char* comments)  /* comment string (set 0 for no comments) */

/*
writes a greyscale image in double format into a pgm P5 file
*/

{
    FILE* outimage;     /* output file */
    long i, j;          /* loop variables */
    double aux;         /* auxiliary variable */
    unsigned char byte; /* for data conversion */

    /* open file */
    outimage = fopen(file_name, "wb");
    if (NULL == outimage) {
        printf("could not open file '%s' for writing, aborting\n", file_name);
        exit(1);
    }

    /* write header */
    fprintf(outimage, "P5\n"); /* format */
    if (comments != 0)
        fputs(comments, outimage);            /* comments */
    fprintf(outimage, "%ld %ld\n", nx, ny); /* image size */
    fprintf(outimage, "255\n");             /* maximal value */

    /* write image data */
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++) {
            aux = u[i][j] + 0.499999; /* for correct rounding */
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

} /* write_double_to_pgm */

/*--------------------------------------------------------------------------*/

void dummies_double

(double** u, /* image */
 long nx,    /* size in x direction */
 long ny)    /* size in y direction */

/*
creates dummy boundaries for a double format image u by mirroring
*/

{
    long i, j; /* loop variables */

    for (i = 1; i <= nx; i++) {
        u[i][0] = u[i][1];
        u[i][ny + 1] = u[i][ny];
    }

    for (j = 0; j <= ny + 1; j++) {
        u[0][j] = u[1][j];
        u[nx + 1][j] = u[nx][j];
    }

    return;

} /* dummies_double */

/*--------------------------------------------------------------------------*/

void analyse_grey_double

(double** u,   /* image, unchanged */
 long nx,      /* pixel number in x direction */
 long ny,      /* pixel number in y direction */
 double* min,  /* minimum, output */
 double* max,  /* maximum, output */
 double* mean, /* mean, output */
 double* std)  /* standard deviation, output */

/*
computes minimum, maximum, mean, and standard deviation of a greyscale
image u in double format
*/

{
    long i, j;    /* loop variables */
    double help1; /* auxiliary variable */
    double help2; /* auxiliary variable */

    /* compute maximum, minimum, and mean */
    *min = u[1][1];
    *max = u[1][1];
    help1 = 0.0;
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++) {
            if (u[i][j] < *min)
                *min = u[i][j];
            if (u[i][j] > *max)
                *max = u[i][j];
            help1 = help1 + u[i][j];
        }
    *mean = help1 / (nx * ny);

    /* compute standard deviation */
    *std = 0.0;
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++) {
            help2 = u[i][j] - *mean;
            *std = *std + help2 * help2;
        }
    *std = sqrt(*std / (nx * ny));

    return;

} /* analyse_grey_double */

/*--------------------------------------------------------------------------*/

void ambrosio

(double tau,   /* time step size */
 long nx,      /* image dimension in x direction */
 long ny,      /* image dimension in y direction */
 double hx,    /* pixel size in x direction */
 double hy,    /* pixel size in y direction */
 double beta,  /* weight of the data term */
 double alpha, /* weight of the boundary term */
 double c,     /* gamma convergence parameter */
 double** f,   /* input: initial image, unchanged */
 double** u,   /* input: actual image; output: smoothed */
 double** v)   /* input: actual edge map; output: smoothed */

/*
One explicit update step of the Ambrosio-Tortorelli functional.
Implicit stabilisation of the reaction terms.
*/

{
    long i, j;                     /* loop variables */
    double rxx, ryy;               /* time savers */
    double aux1, aux2, aux3;       /* time savers */
    double aux4, aux5, aux6;       /* time savers */
    double** dc;                   /* diffusivity */
    double** uo, ** vo;             /* old values of u, v */
    double duo_dx, duo_dy;         /* derivatives of uo */
    double inv_two_hx, inv_two_hy; /* time savers */
    double** grad_sqr;             /* |grad(uo)|^2 */

    /* ---- allocate memory ---- */

    alloc_double_matrix(&uo, nx + 2, ny + 2);
    alloc_double_matrix(&vo, nx + 2, ny + 2);
    alloc_double_matrix(&dc, nx + 2, ny + 2);
    alloc_double_matrix(&grad_sqr, nx + 2, ny + 2);

    /* ---- create copies of u and v and assign dummy b.c. ---- */

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++) {
            uo[i][j] = u[i][j];
            vo[i][j] = v[i][j];
        }
    dummies_double(uo, nx, ny);
    dummies_double(vo, nx, ny);

    /* ---- compute diffusivity dc and |grad(uo)|^2 ---- */

    inv_two_hx = 1.0 / (2.0 * hx);
    inv_two_hy = 1.0 / (2.0 * hy);

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++) {
            dc[i][j] = vo[i][j] * vo[i][j];
            duo_dx = (uo[i + 1][j] - uo[i - 1][j]) * inv_two_hx;
            duo_dy = (uo[i][j + 1] - uo[i][j - 1]) * inv_two_hy;
            grad_sqr[i][j] = duo_dx * duo_dx + duo_dy * duo_dy;
        }

      /* ---- compute updates of u and v ---- */

      /* dummy boundaries for the diffusivity */
    dummies_double(dc, nx, ny);

    /* compute time savers */
    rxx = tau / (2.0 * hx * hx);
    ryy = tau / (2.0 * hy * hy);
    aux1 = tau * beta;
    aux2 = 1.0 / (1.0 + aux1);
    aux3 = c * tau / (hx * hx);
    aux4 = c * tau / (hy * hy);
    aux5 = tau / (4.0 * c);
    aux6 = tau / alpha;
    double tau_lap;
    /* update u and v */
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++) {
          /*
          SUPPLEMENT CODE
          */
            tau_lap = 2 * rxx * (uo[i + 1][j] - 2 * uo[i][j] + uo[i - 1][j]) +
                2 * ryy * (uo[i][j + 1] - 2 * uo[i][j] + uo[i][j - 1]);

    //   A = rxx * ((vo[i + 1][j] + vo[i - 1][j]) * (uo[i + 1][j] - uo[i][j]) -
    //              (vo[i - 1][j] + vo[i][j]) * (uo[i][j] - uo[i - 1][j])) +
    //       ryy * ((vo[i][j + 1] + vo[i][j]) * (uo[i][j + 1] - uo[i][j]) -
    //              (vo[i][j - 1] + vo[i][j]) * (uo[i][j] - uo[i][j - 1]));

            u[i][j] = aux2 * (2 * vo[i][j] *
                                  (rxx * (vo[i + 1][j] - vo[i - 1][j]) *
                                      (uo[i + 1][j] - uo[i - 1][j]) +
                                      ryy * (vo[i][j + 1] - vo[i][j - 1]) *
                                      (uo[i][j + 1] - uo[i][j - 1])) +
                              dc[i][j] * tau_lap + aux1 * f[i][j] + uo[i][j]);

            v[i][j] = vo[i][j] + aux3 * (vo[i + 1][j] - 2 * vo[i][j] + vo[i - 1][j]) +
                aux4 * (vo[i][j + 1] - 2 * vo[i][j] + vo[i][j - 1]) -
                aux6 * grad_sqr[i][j] * vo[i][j] + aux5 * (1 - vo[i][j]);
        }

      /* ---- free memory ---- */

    free_double_matrix(dc, nx + 2, ny + 2);
    free_double_matrix(uo, nx + 2, ny + 2);
    free_double_matrix(vo, nx + 2, ny + 2);
    free_double_matrix(grad_sqr, nx + 2, ny + 2);

    return;

} /* ambrosio */

/*--------------------------------------------------------------------------*/

int main()

{
    char in[80];         /* for reading data */
    char out1[80];       /* for writing data */
    char out2[80];       /* for reading data */
    double** f;          /* initial image */
    double** u;          /* image */
    double** v;          /* edge map */
    long i, j, k;        /* loop variables */
    long nx, ny;         /* image size in x, y direction */
    double tau;          /* time step size */
    long kmax;           /* largest iteration number */
    double beta;         /* weight of the data term */
    double alpha;        /* weight of the boundary term */
    double c;            /* gamma convergence parameter */
    double max, min;     /* largest, smallest grey value */
    double mean;         /* average grey value */
    double std;          /* standard deviation */
    char comments[1600]; /* string for comments */

    printf("\n");
    printf("AMBROSIO-TORTORELLI FUNCTIONAL FOR GREYSCALE (PGM) IMAGES;\n");
    printf("STABILISED EXPLICIT GRADIENT DESCENT SCHEME\n\n");
    printf("**************************************************\n\n");
    printf("    Copyright 2020 by Joachim Weickert            \n");
    printf("    Dept. of Mathematics and Computer Science     \n");
    printf("    Saarland University, Saarbruecken, Germany    \n\n");
    printf("    All rights reserved. Unauthorised usage,      \n");
    printf("    copying, hiring, and selling prohibited.      \n\n");
    printf("    Send bug reports to                           \n");
    printf("    weickert@mia.uni-saarland.de                  \n\n");
    printf("**************************************************\n\n");

    /* ---- read input image f (pgm format P5) ---- */

    // printf("input image (pgm):                  ");
    // read_string(in);
    strcpy(in, "house.pgm");
    read_pgm_to_double(in, &nx, &ny, &f);

    /* ---- read other parameters ---- */

    // printf("data term weight beta (>0.0):       ");
    // read_double(&beta);

    // printf("boundary term weight alpha (>0.0):  ");
    // read_double(&alpha);

    // printf("parameter c (>0.0, <<1.0):          ");
    // read_double(&c);

    // printf("time step size tau (<=0.25):        ");
    // read_double(&tau);

    // printf("number of iterations:               ");
    // read_long(&kmax);

    // printf("output image u:                     ");
    // read_string(out1);

    // printf("edge map v:                         ");
    // read_string(out2);
    // printf("\n");

    beta = 0.007;
    alpha = 1;
    c = 0.05;
    tau = 0.2;
    kmax = 50000;
    strcpy(out1, "out1u.pgm");
    strcpy(out2, "out1v.pgm");

    /* ---- allocate memory ---- */

    alloc_double_matrix(&u, nx + 2, ny + 2);
    alloc_double_matrix(&v, nx + 2, ny + 2);

    /* ---- initialise u and v ---- */

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++) {
            u[i][j] = f[i][j];
            v[i][j] = 0.0;
        }

      /* ---- analyse initial image ---- */

    analyse_grey_double(f, nx, ny, &min, &max, &mean, &std);
    printf("initial image\n");
    printf("minimum of f:        %8.2lf \n", min);
    printf("maximum of f:        %8.2lf \n", max);
    printf("mean of f:           %8.2lf \n", mean);
    printf("standard dev. of f:  %8.2lf \n\n", std);

    /* ---- process image u and edge map v ---- */

    for (k = 1; k <= kmax; k++) {
      /* perform one iteration */
        ambrosio(tau, nx, ny, 1.0, 1.0, beta, alpha, c, f, u, v);

        /* analyse image u */
        analyse_grey_double(u, nx, ny, &min, &max, &mean, &std);
        printf("iteration:           %8ld   \n", k);
        printf("minimum of u:        %8.2lf \n", min);
        printf("maximum of u:        %8.2lf \n", max);
        printf("mean of u:           %8.2lf \n", mean);
        printf("standard dev. of u:  %8.2lf \n\n", std);
    } /* for */

    /* ---- write regularised image u (pgm format P5) ---- */

    /* write parameter values in comment string */
    comments[0] = '\0';
    comment_line(comments, "# regularised image u of Ambrosio-Tortorelli\n");
    comment_line(comments, "# stablised explicit gradient descent scheme\n");
    comment_line(comments, "# initial image:  %s\n", in);
    comment_line(comments, "# beta:         %8.4lf\n", beta);
    comment_line(comments, "# alpha:        %8.4lf\n", alpha);
    comment_line(comments, "# c:            %8.4lf\n", c);
    comment_line(comments, "# tau:          %8.4lf\n", tau);
    comment_line(comments, "# kmax:         %8ld\n", kmax);
    comment_line(comments, "# min:          %8.2lf\n", min);
    comment_line(comments, "# max:          %8.2lf\n", max);
    comment_line(comments, "# mean:         %8.2lf\n", mean);
    comment_line(comments, "# st. dev.:     %8.2lf\n", std);

    /* write image data */
    write_double_to_pgm(u, nx, ny, out1, comments);
    printf("regularised image %s successfully written\n", out1);

    /* ---- write edge map image v (pgm format P5) ---- */

    /* write parameter values in comment string */
    comments[0] = '\0';
    comment_line(comments, "# edge map v of Ambrosio-Tortorelli\n");
    comment_line(comments, "# stablised explicit gradient descent scheme\n");
    comment_line(comments, "# initial image:  %s\n", in);
    comment_line(comments, "# beta:         %8.4lf\n", beta);
    comment_line(comments, "# alpha:        %8.4lf\n", alpha);
    comment_line(comments, "# c:            %8.4lf\n", c);
    comment_line(comments, "# tau:          %8.4lf\n", tau);
    comment_line(comments, "# kmax:         %8ld\n", kmax);

    /* rescale v from [0, 1] to [0, 255] */
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            v[i][j] = 255.0 * v[i][j];

        /* write edge map data */
    write_double_to_pgm(v, nx, ny, out2, comments);
    printf("edge map image %s successfully written\n\n", out2);

    /* ---- free memory ---- */

    free_double_matrix(f, nx + 2, ny + 2);
    free_double_matrix(u, nx + 2, ny + 2);
    free_double_matrix(v, nx + 2, ny + 2);

    return (0);

} /* main */
