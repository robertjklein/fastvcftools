/*
 * fastvcftools.c
 *
 * Fast implementation of some vcftools routines.
 *
 * Some of the speed comes from making assumptions about the data, and aborting if assumptions are not met.
 *
 * Initial implementation only does haplotype-based r2 calculation with default parameters.  If it works well, I'll grow it.
 * 
 * Robert J. Klein
 * July 24, 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>

#define MAXSNPS 10000
#define DIST 1000000
#define R2 0.1

#define MallocOrDie(x)     sre_malloc(__FILE__, __LINE__, (x))
#define ReallocOrDie(x,y)  sre_realloc(__FILE__, __LINE__, (x), (y))

void
Die(char *format, ...)
{
  va_list  argp;
                                /* format the error mesg */
  fprintf(stderr, "\nFATAL: ");
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);
                                /* exit  */
  exit(1);
}




void *
sre_malloc(char *file, int line, size_t size)
{
  void *ptr;

  if ((ptr = malloc (size)) == NULL)
    Die("malloc of %ld bytes failed: file %s line %d", size, file, line);
  return ptr;
}

void *
sre_realloc(char *file, int line, void *p, size_t size)
{
  void *ptr;

  if ((ptr = realloc(p, size)) == NULL)
    Die("realloc of %ld bytes failed: file %s line %d", size, file, line);
  return ptr;
}




typedef struct _vcf_file_t {
  FILE *f;
  int n;
  char **names;
} vcf_file_t;

vcf_file_t *open_and_initialize (char *filename) {
  vcf_file_t *retval;
  char *cp;
  char *cp2;
  char buf[65535];
  int i;
  int read_header = 0;

  retval = MallocOrDie(sizeof(vcf_file_t));
  
  if (filename[0] == '-' && filename[1] == '\0')
    retval->f = stdin;
  else {
    cp = strstr(filename, ".gz");
    if (cp == NULL || (cp[3] != '\0' && !isspace(cp[3])))
      retval->f = fopen(filename, "r");
    else {
      sprintf (buf, "gzcat %s", filename);
      retval->f = popen(buf, "r");
    }
  }
  
  if (retval->f == NULL) 
    Die("Could not open %s\n", filename);

  while (read_header == 0 && fgets(buf, 65534, retval->f)) {
    if (buf[0] == '#' && buf[1] == '#') continue;
    if (buf[0] == '#') {
      cp = buf;
      i = 0;
      while (i < 9) {
	while (!isspace(*cp)) cp++;
	while (isspace(*cp)) cp++;
	i++;
      }
      i = 0;
      cp2 = cp;
      while (*cp != '\0' && !isspace(*cp)) {
	while (!isspace(*cp)) cp++;
	while (isspace(*cp)) cp++;
	i++;
      }
      retval->n = i;
      retval->names = MallocOrDie(sizeof(char *)*retval->n);
      cp = cp2;
      i = 0;
      while (*cp != '\0' && !isspace(*cp)) {
	cp2 = cp;
	while (!isspace(*cp2)) cp2++;
	retval->names[i] = MallocOrDie(sizeof(char)*(cp2-cp+1));
	strncpy (retval->names[i], cp, cp2-cp);
	retval->names[i][cp2-cp] = '\0';
	while (!isspace(*cp)) cp++;
	while (isspace(*cp) && *cp != '\0') cp++;
	i++;
      }
      if (i != retval->n) { Die("i!=n\n"); }
      read_header = 1;
    }
  }
  return(retval);
}

typedef struct _vcf_line_t {
  char *chr;
  int pos;
  unsigned int *zeros;
  unsigned int *ones;
  unsigned int num_ints;
  struct _vcf_line_t *next;
} vcf_line_t;

void set_bits (unsigned int *zeros, unsigned int *ones, int i, char a, char b) {
  int index1, mask1, index2, mask2;

  index1 = (2*i) / sizeof(unsigned int);
  mask1 = 1 << ((2*i) % sizeof(unsigned int));

  index2 = (2*i + 1) / sizeof(unsigned int);
  mask2 = 1 << ((2*i + 1) % sizeof(unsigned int));
  
  if (a == '0')
    zeros[index1] |= mask1;
  else if (a == '1')
    ones[index1] |= mask1;
  
  if (b == '0') 
    zeros[index2] |= mask2;
  else if (b == '1')
    ones[index2] |= mask2;
}
  

vcf_line_t *read_vcf_line (vcf_file_t *vcf) {
  vcf_line_t *retval;
  char buf[65536];
  char *cp;
  int i;
  int num_ints;
  FILE *f;

  f=vcf->f;

  if (! fgets(buf, 65535, f))
      return(NULL);

  retval = MallocOrDie(sizeof(vcf_line_t));
  retval->next = NULL;
  cp = buf;
  while (!isspace(*cp)) cp++;
  retval->chr = MallocOrDie(sizeof(char)*(cp-buf+1));
  strncpy(retval->chr, buf, cp-buf);
  retval->chr[cp-buf]='\0';
  while (isspace(*cp))cp++;
  retval->pos=atoi(cp);
  for (i=0; i<8; i++) {
    while (!isspace(*cp)) cp++;
    while(isspace(*cp)) cp++;
  }
  
  num_ints = ((2 * vcf->n) / sizeof(unsigned int))+1;
  retval->zeros=MallocOrDie(sizeof(unsigned int)*num_ints);
  retval->ones=MallocOrDie(sizeof(unsigned int)*num_ints);
  retval->num_ints = num_ints;
  for (i=0; i<num_ints; i++) {
    retval->zeros[i] = 0;
    retval->ones[i] = 0;
  }
  
  for (i=0; i<vcf->n; i++) {
    if (cp[1] != '|' || cp[3] != ':') 
      Die("gt error cp = %s\n\tbuf = %s\n", cp, buf);
    set_bits(retval->zeros, retval->ones, i, cp[0], cp[2]);
    while (!isspace(*cp)) cp++;
    while(isspace(*cp)) cp++;
  }
  return(retval);
}

unsigned int count (unsigned int *p, unsigned int n) {
  int i;
  unsigned int sum = 0;
  for (i=0; i<n; i++)
    sum += __builtin_popcount(p[i]);
  return(sum);
}

inline double min (double a, double b) {
  if (a<b) return(a);
  else return(b);
}

void compute_and_print_r2 (vcf_line_t *a, vcf_line_t *b, unsigned int *workspace) {
  int x11, x12, x21, x22, tot;
  double p1, p2, q1, q2, x11r, x12r, x21r, x22r;
  int i;
  double r2, D, Dmax, Dp;


  for (i=a->num_ints-1; i>=0; i--)
    workspace[i] = a->zeros[i] & b->zeros[i];
  x11 = count(workspace, a->num_ints);

  for (i=a->num_ints-1; i>=0; i--)
    workspace[i] = a->zeros[i] & b->ones[i];
  x12 = count(workspace, a->num_ints);

  for (i=a->num_ints-1; i>=0; i--)
    workspace[i] = a->ones[i] & b->zeros[i];
  x21 = count(workspace, a->num_ints);

  for (i=a->num_ints-1; i>=0; i--)
    workspace[i] = a->ones[i] & b->ones[i];
  x22 = count(workspace, a->num_ints);

  tot = x11+x12+x21+x22;

  x11r = 1.0 * x11 / tot;
  x12r = 1.0 * x12 / tot;
  x21r = 1.0 * x21 / tot;
  x22r = 1.0 * x22 / tot;
  
  p1 = x11r + x12r;
  p2 = x21r + x22r;
  q1 = x11r + x21r;
  q2 = x12r + x22r;

  D = x11r - p1 * q1;

  if (D < 0) {
    Dmax = min (p1 * q1, p2 * q2);
  } else {
    Dmax = min (p1 * q2, p2 * q1);
  };
  Dp = D / Dmax;

  r2 = (D * D) / (p1 * p2 * q1 * q2);
  
  if (r2 >= R2)
    printf ("%s\t%d\t%d\t%d\t%f\t%f\t%f\n", a->chr, a->pos, b->pos, tot, r2, D, Dp);
}

void do_r2 (vcf_file_t *vcf_file) {
  vcf_line_t *front;
  vcf_line_t *back;
  vcf_line_t *cur;
  int flag = 0;
  unsigned int *workspace;
  
  front=read_vcf_line(vcf_file);
  back = front;

  workspace = MallocOrDie(sizeof (unsigned int) * front->num_ints);

  while (front != NULL) {
    cur = front->next;
    flag = 0;
    do {
      if (cur == NULL) {
	back->next = read_vcf_line(vcf_file);
	cur = back->next;
	if (cur != NULL) back = back->next;
      }
      if (cur == NULL) {
	flag = 1;
      } else {
	if ((strcmp(front->chr, cur->chr) != 0)
	    || (cur->pos - front->pos > DIST)) {
	  flag = 1;
	} else {
	  compute_and_print_r2 (front, cur, workspace);
	}
      }
      if (cur != NULL) 
	cur = cur->next;
    } while (flag == 0);
    if (front != NULL) {
      cur = front->next;
      free(front->zeros);
      free(front->ones);
      free(front);
      front = cur;
    }
  }
}

int main (int argc, char **argv) {
  vcf_file_t *vcf_file;

  if (argc != 2)
    Die("USAGE: fastvcftools <file>\n");

  vcf_file = open_and_initialize (argv[1]);

  do_r2 (vcf_file);
}
