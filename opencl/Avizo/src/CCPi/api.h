// AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY.
#ifndef CCPI_API_EXPORT_MACRO_H
#define CCPI_API_EXPORT_MACRO_H

#ifdef CCPI_STATIC
#   define CCPI_API
#else
#   ifdef _WIN32
#       ifdef CCPI_EXPORTS
#           define CCPI_API __declspec(dllexport)
#       else
#           define CCPI_API __declspec(dllimport)
#       endif
#   else
#       define CCPI_API __attribute__ ((visibility("default")))
#   endif
#endif

#endif
