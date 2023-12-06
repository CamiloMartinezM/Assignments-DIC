#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                    FISTA METHOD FOR TV REGULARISATION                    */
/*                                                                          */
/*  (Copyright by Sven Grewenig, Simon Setzer, Joachim Weickert, 11/2021)   */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

void alloc_double_matrix

    (double ***matrix, /* matrix */
     long n1,          /* size in direction 1 */
     long n2)          /* size in direction 2 */

/*
allocates memory for a double format matrix of size n1 * n2
*/

{
  long i; /* loop variable */

  *matrix = (double **)malloc(n1 * sizeof(double *));

  if (*matrix == NULL) {
    printf("alloc_double_matrix: not enough memory available\n");
    exit(1);
  }

  for (i = 0; i < n1; i++) {
    (*matrix)[i] = (double *)malloc(n2 * sizeof(double));
    if ((*matrix)[i] == NULL) {
      printf("alloc_double_matrix: not enough memory available\n");
      exit(1);
    }
  }

  return;

} /* alloc_double_matrix */

/*--------------------------------------------------------------------------*/

void free_double_matrix

    (double **matrix, /* matrix */
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

    (char *v) /* string to be read */

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

    (long *v) /* value to be read */

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

    (double *v) /* value to be read */

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

    (FILE *inimage) /* input file */

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
  } else
    fseek(inimage, -1, SEEK_CUR);

  return;

} /* skip_white_space_and_comments */

/*--------------------------------------------------------------------------*/

void read_pgm_to_double

    (const char *file_name, /* name of pgm file */
     long *nx,              /* image size in x direction, output */
     long *ny,              /* image size in y direction, output */
     double ***u)           /* image, output */

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
  FILE *inimage;  /* input file */

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
  } else {
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

    (char *comment,    /* comment string (output) */
     char *lineformat, /* format string for comment line */
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

    (double **u,      /* image, unchanged */
     long nx,         /* image size in x direction */
     long ny,         /* image size in y direction */
     char *file_name, /* name of pgm file */
     char *comments)  /* comment string (set 0 for no comments) */

/*
writes a greyscale image in double format into a pgm P5 file
*/

{
  FILE *outimage;     /* output file */
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

void dummies_Dirichlet_double

    (double **u, /* image */
     long nx,    /* size in x direction */
     long ny)    /* size in y direction */

/*
creates Dirichlet dummy boundaries for a double format image u
*/

{
  long i, j; /* loop variables */

  for (i = 1; i <= nx; i++) {
    u[i][0] = 0.0;
    u[i][ny + 1] = 0.0;
  }

  for (j = 0; j <= ny + 1; j++) {
    u[0][j] = 0.0;
    u[nx + 1][j] = 0.0;
  }

  return;

} /* dummies_Dirichlet_double */

/*--------------------------------------------------------------------------*/

void dummies_Neumann_double

    (double **u, /* image */
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

} /* dummies_Neumann_double */

/*--------------------------------------------------------------------------*/

void analyse_grey_double

    (double **u,   /* image, unchanged */
     long nx,      /* pixel number in x direction */
     long ny,      /* pixel number in y direction */
     double *min,  /* minimum, output */
     double *max,  /* maximum, output */
     double *mean, /* mean, output */
     double *std)  /* standard deviation, output */

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

void proj

    (double alpha, /* regularisation weight */
     long nx,      /* pixel number in x direction */
     long ny,      /* pixel number in y direction */
     double **dx,  /* image derivatives in x direction */
     double **dy)  /* image derivatives in y direction */

/*
computes the projection into the set C
*/

{
  long i, j;   /* loop variables */
  double grad; /* gradient magnitude */

  for (i = 1; i < nx + 1; i++)
    for (j = 1; j < ny + 1; j++) {
      grad = sqrt(dx[i][j] * dx[i][j] + dy[i][j] * dy[i][j]);

      /* compute projection if gradient magnitude is larger than alpha */
      if (grad > alpha) {
        dx[i][j] = dx[i][j] * alpha / grad;
        dy[i][j] = dy[i][j] * alpha / grad;
      }
    }

  return;

} /* proj */

/*--------------------------------------------------------------------------*/

void mult_D

    (long nx,     /* pixel number in x direction */
     long ny,     /* pixel number in y direction */
     double **u,  /* input vector */
     double **dx, /* matrix-vector product, x direction, output */
     double **dy) /* matrix-vector product, y direction, output */

/*
computes the matrix-vector product D*u,
i.e. computes the forward differences w.r.t. each direction
*/

{
  long i, j; /* loop variables */

  dummies_Neumann_double(u, nx, ny);

  /* forward differences */
  for (i = 1; i <= nx; i++)
    for (j = 1; j <= ny; j++) {
      dx[i][j] = u[i + 1][j] - u[i][j];
      dy[i][j] = u[i][j + 1] - u[i][j];
    }

  return;

} /* mult_D */

/*--------------------------------------------------------------------------*/

void mult_DT

    (long nx,     /* pixel number in x direction */
     long ny,     /* pixel number in y direction */
     double **u,  /* resulting vector, output */
     double **dx, /* derivative vector, x direction */
     double **dy) /* derivative vector, y direction */

/*
computes the matrix-vector product u = D^T*d,
i.e. computes the negative backward differences
*/

{
  long i, j; /* loop variables */

  dummies_Dirichlet_double(dx, nx, ny);
  dummies_Dirichlet_double(dy, nx, ny);

  /* negative backward differences */
  for (i = 1; i <= nx; i++)
    for (j = 1; j <= ny; j++)
      u[i][j] = dx[i - 1][j] - dx[i][j] + dy[i][j - 1] - dy[i][j];

  return;

} /* mult_DT */

/*--------------------------------------------------------------------------*/

void FISTA

    (long kmax,    /* number of iterations */
     double alpha, /* regularisation weight */
     long nx,      /* pixel number in x direction */
     long ny,      /* pixel number in y direction */
     double **u)   /* in: input image  out: result */

/*
fast iterative shrinkage thresholding algorithm (FISTA)
*/

{
  long i, j, k;     /* loop variables */
  double **b_oldx;  /* old iteration b^k (x) */
  double **b_oldy;  /* old iteration b^k (y) */
  double **b_newx;  /* new iteration b^{k+1} (x) */
  double **b_newy;  /* new iteration b^{k+1} (y) */
  double **bp_oldx; /* old iteration b^k (x), proj. */
  double **bp_oldy; /* old iteration b^k (y), proj. */
  double **bp_newx; /* new iteration b^{k+1} (x), proj. */
  double **bp_newy; /* new iteration b^{k+1} (y), proj. */
  double **v;       /* auxiliary variable */
  double tau;       /* step size */
  double theta;     /* step size */
  double tk, tk1;   /* step sizes */

  /* ---- allocate memory ---- */

  alloc_double_matrix(&b_oldx, nx + 2, ny + 2);
  alloc_double_matrix(&b_oldy, nx + 2, ny + 2);
  alloc_double_matrix(&b_newx, nx + 2, ny + 2);
  alloc_double_matrix(&b_newy, nx + 2, ny + 2);
  alloc_double_matrix(&bp_oldx, nx + 2, ny + 2);
  alloc_double_matrix(&bp_oldy, nx + 2, ny + 2);
  alloc_double_matrix(&bp_newx, nx + 2, ny + 2);
  alloc_double_matrix(&bp_newy, nx + 2, ny + 2);
  alloc_double_matrix(&v, nx + 2, ny + 2);

  /* ---- initialisations ---- */

  /* compute b^0 = D*u */
  mult_D(nx, ny, u, b_newx, b_newy);

  /* initialise step size tau */
  tau = 0.125;

  /* t_{-1} = 0, t_0 = 1 */
  tk = 0.0;
  tk1 = 1.0;

  /* b^{k+1} = 0 */
  for (i = 0; i < nx + 2; i++)
    for (j = 0; j < ny + 2; j++)
      bp_newx[i][j] = bp_newy[i][j] = 0.0;

  /* ---- loop ---- */

  for (k = 0; k < kmax; k++) {

    /* update step sizes t and theta */
    tk1 = (1 + sqrt(1 + 4 * tk * tk))/2;
    theta = (tk - 1)/tk1;

    /*
    SUPPLEMENT CODE
    */

    /* create b^{k} and projected b^{k} */
    for (i = 1; i < nx + 1; i++)
      for (j = 1; j < ny + 1; j++) {
        /*
        SUPPLEMENT CODE
        */
        b_oldx[i][j] = b_newx[i][j];
        b_oldy[i][j] = b_newy[i][j];
        bp_oldx[i][j] = bp_newx[i][j];
        bp_oldy[i][j] = bp_newy[i][j];
      }

    /* compute v = D^T*b^k */
    mult_DT(nx, ny, v, b_oldx, b_oldy);

    /*
    SUPPLEMENT CODE
    */

    /* v = D^T*b^k - u */
    for (i = 1; i < nx + 1; i++)
      for (j = 1; j < ny + 1; j++)
        v[i][j] = v[i][j] - u[i][j];

    /* compute D*v^k */
    mult_D (nx, ny, v, b_newx, b_newy);

    /*
    SUPPLEMENT CODE
    */

    /* compute the vector to be projected */
    for (i = 1; i < nx + 1; i++)
      for (j = 1; j < ny + 1; j++) {
        /*
        SUPPLEMENT CODE
        */
        bp_newx[i][j] = b_oldx[i][j] - tau * b_newx[i][j];
        bp_newy[i][j] = b_oldy[i][j] - tau * b_newy[i][j];
      }

    /* get projection b^{k+1} */
    proj(alpha, nx, ny, bp_newx, bp_newy);

    /* relaxation step */
    for (i = 1; i < nx + 1; i++)
      for (j = 1; j < ny + 1; j++) {
        /*
        SUPPLEMENT CODE
        */
        b_newx[i][j] = bp_newx[i][j] + theta * (bp_newx[i][j] - bp_oldx[i][j]);
        b_newy[i][j] = bp_newy[i][j] + theta * (bp_newy[i][j] - bp_oldy[i][j]);
      }
  } /* for k */

  /* ---- compute the solution  u - D^T*b ---- */

  mult_DT(nx, ny, v, b_newx, b_newy);

  for (i = 1; i < nx + 1; i++)
    for (j = 1; j < ny + 1; j++)
      u[i][j] = u[i][j] - v[i][j];

  /* ---- free memory ---- */

  free_double_matrix(b_oldx, nx + 2, ny + 2);
  free_double_matrix(b_oldy, nx + 2, ny + 2);
  free_double_matrix(b_newx, nx + 2, ny + 2);
  free_double_matrix(b_newy, nx + 2, ny + 2);
  free_double_matrix(bp_oldx, nx + 2, ny + 2);
  free_double_matrix(bp_oldy, nx + 2, ny + 2);
  free_double_matrix(bp_newx, nx + 2, ny + 2);
  free_double_matrix(bp_newy, nx + 2, ny + 2);
  free_double_matrix(v, nx + 2, ny + 2);

  return;

} /* FISTA */

/*--------------------------------------------------------------------------*/

int main()

{
  char in[80];         /* for reading data */
  char out[80];        /* for writing data */
  double **u;          /* image */
  long kmax;           /* number of iterations */
  long nx, ny;         /* image size in x, y direction */
  double alpha;        /* regularisation parameter */
  double max, min;     /* largest, smallest grey value */
  double mean;         /* average grey value */
  double std;          /* standard deviation */
  char comments[1600]; /* string for comments */

  printf("\n");
  printf("TV REGULARISATION WITH FISTA\n\n");
  printf("*********************************************************\n\n");
  printf("    Copyright 2013-2021 by                               \n");
  printf("    Sven Grewenig, Simon Setzer, and Joachim Weickert    \n");
  printf("    Faculty of Mathematics and Computer Science          \n");
  printf("    Saarland University, Saarbruecken, Germany           \n\n");
  printf("    All rights reserved. Unauthorized usage, copying,    \n");
  printf("    hiring, and selling prohibited.                      \n\n");
  printf("    Send bug reports to                                  \n");
  printf("    weickert@mia.uni-saarland.de                         \n\n");
  printf("*********************************************************\n\n");

  /* ---- read input image (pgm format P5) ---- */

  printf("input image (pgm):                      ");
  read_string(in);
  read_pgm_to_double(in, &nx, &ny, &u);

  /* ---- read other parameters ---- */

  printf("regularisation parameter alpha (>0.0):  ");
  read_double(&alpha);

  printf("number of iterations:                   ");
  read_long(&kmax);

  printf("output image (pgm):                     ");
  read_string(out);

  printf("\n");

  /* ---- analyse input image ---- */
  analyse_grey_double(u, nx, ny, &min, &max, &mean, &std);
  printf("input image\n");
  printf("minimum:          %8.2lf \n", min);
  printf("maximum:          %8.2lf \n", max);
  printf("mean:             %8.2lf \n", mean);
  printf("standard dev.:    %8.2lf \n\n", std);

  /* ---- process image ---- */

  FISTA(kmax, alpha, ny, ny, u);

  /* ---- analyse processed image ---- */

  analyse_grey_double(u, nx, ny, &min, &max, &mean, &std);
  printf("processed image\n");
  printf("minimum:          %8.2lf \n", min);
  printf("maximum:          %8.2lf \n", max);
  printf("mean:             %8.2lf \n", mean);
  printf("standard dev.:    %8.2lf \n\n", std);

  /* ---- write output image (pgm format P5) ---- */

  /* write parameter values in comment string */
  comments[0] = '\0';
  comment_line(comments, "# TV regularisation with FISTA\n");
  comment_line(comments, "# input image:    %s\n", in);
  comment_line(comments, "# alpha:          %8.2f\n", alpha);
  comment_line(comments, "# iterations:     %8ld\n", kmax);
  comment_line(comments, "# minimum:        %8.2f\n", min);
  comment_line(comments, "# maximum:        %8.2f\n", max);
  comment_line(comments, "# mean:           %8.2f\n", mean);
  comment_line(comments, "# standard dev.:  %8.2lf\n", std);

  /* write image data */
  write_double_to_pgm(u, nx, ny, out, comments);
  printf("output image %s successfully written\n\n", out);

  /* ---- free memory ---- */

  free_double_matrix(u, nx + 2, ny + 2);

  return (0);

} /* main */
