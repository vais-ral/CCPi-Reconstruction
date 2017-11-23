/* ----------------------------------------------------------------------
 * time_stamp.cpp
 * ----------------------------------------------------------------------
 * Provide a timestamp function based on gettimeofday() on linux. The 
 * code will compile on windows so no conditional compilation macros are
 * required in calling code.
 *
 * ----------------------------------------------------------------------
 */


#include "time_stamp.h"

#if defined(_WIN32) || defined(_WIN64)
struct _timeb timebuffer;
#endif
struct timeval mystime, now; /* starting time for the whole program */
const char *myname;  /* name of the program */
int verbose = 4;


void init_timestamp()
{
#if defined(__linux__)
	gettimeofday(&mystime,NULL);
#else
	_ftime(&timebuffer);
	mystime.tv_usec = (long)(timebuffer.millitm);
	mystime.tv_sec = (long)(timebuffer.time);
#endif
}

void timestamp(const char *stampmsg, const int vlevel){
	char message[MAX_MESSAGE];
	time_t etimes;
	double nowtime, etime,inctime;
        static double lasttime=0;
	if(vlevel <= verbose){
#if defined(__linux__)
		suseconds_t etimeu;
		gettimeofday(&now, NULL);	
#else
		long etimeu;
		_ftime(&timebuffer);
		now.tv_usec = (long)timebuffer.millitm;
		now.tv_sec = (long)(timebuffer.time);
#endif
		etimes = now.tv_sec - mystime.tv_sec;
		etimeu = now.tv_usec - mystime.tv_usec;
		etime = (double)(etimes)+((double)(etimeu))/(1e6);
		nowtime = (double)(now.tv_sec) + ((double)(now.tv_usec) / 1.0e6);
                inctime = nowtime - lasttime;
                lasttime=nowtime;
#if defined(__linux__)
		snprintf(message, MAX_MESSAGE, "Timestamp: %f %f %f %s\n", etime,inctime, nowtime, stampmsg);
#else
		_snprintf(message, MAX_MESSAGE, "Timestamp: %f  %f %f %s\n", etime,inctime, nowtime, stampmsg);
#endif
		printf(message);
	}
}
