/* Prototypes for the signal handler functions */
#ifdef __linux__
#ifndef HANDLERS_H
#   define HANDLERS_H

#   ifdef GNU_TRAP_FPE
#      include <fenv.h>
#      include <signal.h>
void enable_fpe_traps ();
#   endif
       /* GNU_TRAP_FPE */

extern void float_error (int signo);    /* signal handler for sigfpe */
extern void float_action (int signo,siginfo_t * myinfo,void * mycontext);
extern void catchit (int signo);        /* signal handler -- jump to finish_bb */
#endif
#endif /*__linux__*/
