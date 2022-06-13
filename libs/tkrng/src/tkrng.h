#pragma once

#ifdef WIN32
  #define tkrng_EXPORT __declspec(dllexport)
#else
  #define tkrng_EXPORT
#endif

tkrng_EXPORT void tkrng();
