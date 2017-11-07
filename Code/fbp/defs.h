#ifndef DEFS_DIAMOND_H
#define DEFS_DIAMOND_H

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

using namespace std;

#include "ipp.h"

#include <cufft.h>
#include <cutil.h>


#if defined(_WIN32) || defined(_WIN64)

#include <windows.h>
#include <direct.h>
#include <sys/timeb.h>

#define STAT_TYPE struct _stat
#define STAT_FUNC _stat
#define MKDIR_FUNC(_fname_) _mkdir(_fname_)
#define FSEEK_FUNC _fseeki64


#elif defined(__linux__)

#include <unistd.h>
#include <sys/signal.h>
#include <sys/time.h>

#define STAT_TYPE struct stat
#define STAT_FUNC stat
#define MKDIR_FUNC(_fname_) mkdir(_fname_, S_IRWXU|S_IRWXG|S_IRWXO )
#define FSEEK_FUNC fseek

#ifndef __min
#define __min(a, b) ((a)<(b)?(a):(b))
#endif

#ifndef __max
#define __max(a, b) ((a)>(b)?(a):(b))
#endif

#else

#error "Unsupported platform"
#endif


typedef float2 Complex; 

#include "constants.h"
#include "classes.h"
#include "declarations.h"


#endif // DEFS_DIAMOND_H
