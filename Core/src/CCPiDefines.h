#ifndef CCPIDEFINES_H
#define CCPIDEFINES_H

#if defined(_WIN32) || defined(__WIN32__)
  #if defined(CCPiReconstructionIterative_EXPORTS)  // add by CMake 
    #define  CCPI_EXPORT __declspec(dllexport)
    #define EXPIMP_TEMPLATE
  #else
    #define  CCPI_EXPORT __declspec(dllimport)
    #define EXPIMP_TEMPLATE extern
  #endif /* CCPi_EXPORTS */
#elif defined(linux) || defined(__linux) || defined(__APPLE__)
 #define CCPI_EXPORT
#endif

#endif