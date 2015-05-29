/* options.c: Lesen von Optionen aus der Kommandozeile */
/* options.c: Read options from the command line       */
/* Volker Ahlers, November 1998                        */
/* Extended by Fernando Perez-Diaz, May 2015           */

void opt_error(char *message);
void opt_cat(char *opts, va_list arguments, int numOpts);
void opt_rm(int *argc, char **argv, int pos, int no);
int opt_ignored(int *argc, char **argv);
int opt_save(int *argc, char **argv, char *file, char *mode);
int opt_find(int *argc, char **argv, va_list arguments, int numOpts);
int opt_flag(int *argc, char **argv, int numOpts, ...);
int opt_int(int *argc, char **argv, int *val, int numOpts, ...);
int opt_float(int *argc, char **argv, float *val, int numOpts, ...);
int opt_double(int *argc, char **argv, double *val, int numOpts, ...);
int opt_string(int *argc, char **argv, char *val, int numOpts, ...);
int arg_string(int *argc, char **argv, char *nostart, char *val);
