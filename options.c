/* options.c: Lesen von Optionen aus der Kommandozeile */
/* options.c: Read options from the command line       */
/* Volker Ahlers, November 1998                        */
/* Extended by Fernando Perez-Diaz, May 2015           */ 

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>
#include "options.h"

void opt_error(char *message)
{
    fprintf(stderr, "\nError in library `options':\n");
    fprintf(stderr, "%s\n\n", message);
    exit(1);
}

void opt_cat(char *opts, va_list arguments, int numOpts)
{
    int j;

    for(j = 0; j < numOpts; j++) {
        strcat(opts, va_arg(arguments, char*));
        strcat(opts, " ");
    }
}

void opt_rm(int *argc, char **argv, int pos, int no)
{
    int i;
  
    if(pos<1 || pos>*argc-no) opt_error("0 < pos <= argc-no required in function opt_rm");
    for(i=pos;i<*argc-no;i++) argv[i] = argv[i+no];
    *argc -= no;
}


int opt_ignored(int *argc, char **argv)
{
    int i;
  
    if(*argc > 1) {
        fprintf(stderr, "\nIgnored option(s):");
        for(i = 1; i < *argc; i++) fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n");
    }
    return *argc - 1;
}


int opt_save(int *argc, char **argv, char *file, char *mode)
{
    int i, stat;
    FILE *out;
  
    out = fopen(file, mode);
    if(out == NULL) stat = 1;
    else {
        stat = fprintf(out, "#");
        for(i = 0; i < *argc; i++) stat = fprintf(out, " %s", argv[i]);
        stat = fprintf(out, "\n");
        if(stat > 0) stat = fclose(out);
        else fclose(out);
    }
    return stat;
}

 
int opt_find(int *argc, char **argv, va_list arguments, int numOpts)
{ 
    int i, j, cond = 1;
    char *opt;

    for(j = 0; j < numOpts; j++) {
        i = 0;
        opt = va_arg(arguments, char*);
        while(i < *argc && strcmp(argv[i], opt)) i++;
        if(i < *argc) return i;
    }
    return 0;
}


int opt_flag(int *argc, char **argv, int numOpts, ...)
{
    int i;
    va_list arguments;
 
    va_start(arguments, numOpts);
    i = opt_find(argc, argv, arguments, numOpts);
    va_end(arguments);

    if(i == 0) return 0;
    opt_rm(argc, argv, i, 1);
    return 1;
}


int opt_int(int *argc, char **argv, int *val, int numOpts, ...)
{
    int i;
    char *rest, message[128], opts[] = "";
    va_list arguments;

    va_start(arguments, numOpts);
    i = opt_find(argc, argv, arguments, numOpts);
    va_end(arguments);

    if(i == 0) return 0;
    if(i < *argc - 1) {
        *val = (int)strtol(argv[i+1], &rest, 10);
        if(!strlen(rest) && !errno) {
            opt_rm(argc, argv, i, 2);
            return 1;
        }
    }

    va_start(arguments, numOpts);
    opt_cat(opts, arguments, numOpts);
    va_end(arguments);
    sprintf(message, "option %s needs an int-type argument", opts);
    opt_error(message);
}


int opt_long(int *argc, char **argv, long *val, int numOpts, ...)
{
    int i;
    char *rest, message[128], opts[] = "";
    va_list arguments;

    va_start(arguments, numOpts);
    i = opt_find(argc, argv, arguments, numOpts);
    va_end(arguments);

    if(i == 0) return 0;
    if(i < *argc - 1) {
        *val = strtol(argv[i+1], &rest, 10);
        if(!strlen(rest) && !errno) {
            opt_rm(argc, argv, i, 2);
            return 1;
        }
    }

    va_start(arguments, numOpts);
    opt_cat(opts, arguments, numOpts);
    va_end(arguments);
    sprintf(message, "option %s needs a long-type argument", opts);
    opt_error(message);
}


int opt_float(int *argc, char **argv, float *val, int numOpts, ...)
{
    int i;
    char *rest, message[128], opts[] = "";
    va_list arguments;

    va_start(arguments, numOpts);
    i = opt_find(argc, argv, arguments, numOpts);
    va_end(arguments);

    if(i == 0) return 0;
    if(i < *argc - 1) {
        *val = (float)strtod(argv[i+1], &rest);
        if(!strlen(rest) && !errno) {
            opt_rm(argc, argv, i, 2);
            return 1;
        }
    }

    va_start(arguments, numOpts);
    opt_cat(opts, arguments, numOpts);
    va_end(arguments);
    sprintf(message, "option %s needs a float-type argument", opts);
    opt_error(message);
}


int opt_double(int *argc, char **argv, double *val, int numOpts, ...)
{
    int i;
    char *rest, message[128], opts[] = "";
    va_list arguments;

    va_start(arguments, numOpts);
    i = opt_find(argc, argv, arguments, numOpts);
    va_end(arguments);

    if(i == 0) return 0;
    if(i < *argc - 1) {
        *val = strtod(argv[i+1], &rest);
        if(!strlen(rest) && !errno) {
            opt_rm(argc, argv, i, 2);
            return 1;
        }
    }

    va_start(arguments, numOpts);
    opt_cat(opts, arguments, numOpts);
    va_end(arguments);
    sprintf(message, "option %s needs a double-type argument", opts);
    opt_error(message);
}


int opt_string(int *argc, char **argv, char *val, int numOpts, ...)
{
    int i;
    char *rest, message[128], opts[] = "";
    va_list arguments;

    va_start(arguments, numOpts);
    i = opt_find(argc, argv, arguments, numOpts);
    va_end(arguments);

    if(i == 0) return 0;
    if(i < *argc-1) {
        strcpy(val, argv[i+1]);
        opt_rm(argc, argv, i, 2);
        return 1;
    }

    va_start(arguments, numOpts);
    opt_cat(opts, arguments, numOpts);
    va_end(arguments);
    sprintf(message, "option %s needs a string-type argument", opts);
    opt_error(message);
}


int arg_string(int *argc, char **argv, char *nostart, char *val)
{
    int i = 1;
  
    if(*argc < 2) return 0;
    while(i < *argc && strchr(nostart, argv[i][0])!=NULL) i++;
    if(i == *argc) return 0;
    strcpy(val, argv[i]);
    opt_rm(argc, argv, i, 1);
    return 1;
}
  
