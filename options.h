/* options.h: Lesen von Optionen aus der Kommandozeile */
/* Volker Ahlers, November 1998                        */

void opt_error(char *message);
void opt_rm(int *argc,char **argv,int pos,int no);
int opt_ignored(int *argc,char **argv);
int opt_save(int *argc,char **argv,char *file,char *mode);
int opt_flag(int *argc,char **argv,char *opt);
int opt_int(int *argc,char **argv,char *opt,int *val);
int opt_float(int *argc,char **argv,char *opt,float *val);
int opt_double(int *argc,char **argv,char *opt,double *val);
int opt_string(int *argc,char **argv,char *opt,char *val);
int arg_string(int *argc,char **argv,char *nostart,char *val);
