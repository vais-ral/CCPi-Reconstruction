// AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY.
#ifndef CCPIRECONSTRUCTION_API_EXPORT_MACRO_H
#define CCPIRECONSTRUCTION_API_EXPORT_MACRO_H

#ifdef CCPIRECONSTRUCTION_STATIC
#   define CCPIRECONSTRUCTION_API
#else
#   ifdef _WIN32
#       ifdef CCPIRECONSTRUCTION_EXPORTS
#           define CCPIRECONSTRUCTION_API __declspec(dllexport)
#       else
#           define CCPIRECONSTRUCTION_API __declspec(dllimport)
#       endif
#   else
#       define CCPIRECONSTRUCTION_API __attribute__ ((visibility("default")))
#   endif
#endif

#endif
