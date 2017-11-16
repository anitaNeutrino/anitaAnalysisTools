#ifndef ACCLAIM_OPENMP_H
#define ACCLAIM_OPENMP_H

#ifdef ACCLAIM_OPENMP
#include <omp.h>
#endif




namespace Acclaim {

  /**
   * @namespace OpenMP
   * @brief Wraps the OpenMP functions with the necessary include guards in one place
   * Perhaps a little overengineered
   */

  namespace OpenMP {


#ifdef ACCLAIM_OPENMP
    const bool isEnabled = true; /// Defined as true if compiled with OpenMP support, false otherwise
#else
    const bool isEnabled = false;/// Defined as true if compiled with OpenMP support, false otherwise
#endif

    /** 
     * Access the currently existing number of threads
     * To know how many threads there will be, or how many
     * you should prepare for, use getMaxThreads()
     * @return the current number of threads
     */
    
    inline int getNumThreads(){
#ifdef ACCLAIM_OPENMP
      return omp_get_num_threads();
#else
      return 1;
#endif
    }



    /** 
     * Gets the maxmimum number of threads, (probably the same as the number of machine cores)
     * This is how many threads you should exect in normal circumstances
     * 
     * @return The maximum number of threads
     */    
    inline int getMaxThreads(){
#ifdef ACCLAIM_OPENMP
      return omp_get_max_threads();
#else
      return 1;
#endif
    }


    /** 
     * Which thread am I in?
     * 
     * @return the thread index, from 0 to getNumThreads() - 1
     */
    inline int thread(){
#ifdef ACCLAIM_OPENMP
      return omp_get_thread_num();
#else
      return 0;
#endif
    }
  }
}

#endif //ACCLAIM_OPENMP_H
