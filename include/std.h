#ifndef STD_H

#define STD_H

using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include <vector>
#include <deque>
#include <math.h>
#include <ctype.h>
#include <sys/time.h>
#include <assert.h>
#include <iomanip>
#include <boost/functional/hash/hash.hpp>
#include <boost/algorithm/string.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
 #include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cdf.h>
     
#define BOOST_UBLAS_NDEBUG 1
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
     
#include <getopt.h>
     
#define UMILLION 1000000ULL
#ifdef WINDOWS
	void gettimeofday( struct timeval* p, void* );
#endif

#ifdef _WIN32
#define rdtsc __rdtsc
#else

unsigned long long rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc": "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}
#endif


#endif
