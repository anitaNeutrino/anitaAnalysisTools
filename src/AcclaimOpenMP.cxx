#include "AcclaimOpenMP.h"
// #include  <algorithm>
// #include  <numeric>

#ifdef ACCLAIM_OPENMP
#include <omp.h>
#endif


int Acclaim::OpenMP::getNumThreads(){
#ifdef ACCLAIM_OPENMP
  return omp_get_num_threads();
#else
  return 1;
#endif  
}


int Acclaim::OpenMP::thread(){
#ifdef ACCLAIM_OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif  
}

int Acclaim::OpenMP::getMaxThreads(){
#ifdef ACCLAIM_OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif  
}



// int Acclaim::OpenMP::setNumThreads(int n){
// #ifdef ACCLAIM_OPENMP
//   omp_set_thread_num(n);
// #else
//   (void) n;
// #endif  

// }


// template <class T>
// Acclaim::OpenMP::ThreadVar<T>::ThreadVar(T initialValue){
//   values.resize(getNumThreads(), initialValue);
// }

// template <class T>
// void Acclaim::OpenMP::ThreadVar<T>::operator++(int){
//   values[thread()]++;
// }


// template <class T>
// T Acclaim::OpenMP::ThreadVar<T>::sum(){
//   return std::accumulate(values.begin(), values.end());
// }
