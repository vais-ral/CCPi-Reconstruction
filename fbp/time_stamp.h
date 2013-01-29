/* ----------------------------------------------------------------------
 * timestamp.h
 * ----------------------------------------------------------------------
 * Linux timer based on gettimeofday(). This can be compiled on windows
 * but it does not provide any timer functionality. 
 *
 * TODO: Implement windows functionality.
 * ----------------------------------------------------------------------
 */
#ifndef TIMESTAMP_H
#define TIMESTAMP_H
#include <stdio.h>
#if defined(__linux__)
#include <unistd.h>
#include <sys/signal.h>
#include <sys/time.h>
#else
#include <sys/timeb.h>
#include <winsock.h>
#endif
#include <time.h>
#include <errno.h>
#include "constants.h"


#if defined(_WIN32) || defined(_WIN64)
extern struct _timeb timebuffer;
#endif
extern struct timeval mystime; /* starting time for the whole program */
extern struct timeval now; /* starting time for the whole program */
extern const char *myname;  /* name of the program */
extern int verbose;

extern void init_timestamp();
extern void timestamp(const char *stampmsg, const int vlevel);
#endif //TIMESTAMP_H
