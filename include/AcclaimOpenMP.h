#ifndef ACCLAIM_OPENMP_H
#define ACCLAIM_OPENMP_H

#include <vector>
// #include <algorithm>
// #include <numeric>
#include <iostream>

namespace Acclaim {

  /**
   * @namespace OpenMP
   * @brief Wraps the OpenMP functions with the necessary include guards in one place
   * Perhaps a little overengineered
   */

  namespace OpenMP {
    int getNumThreads();
    int getMaxThreads();
    int thread();

    // This is an incomplete reimplementation of this...
    // https://computing.llnl.gov/tutorials/openMP/#REDUCTION
    enum ReductionMethod {
      kNone,
      kSum,
      kProduct,
      kMax,
      kMin
    };
    
    template <class T>
    class ThreadVar {
    public:
      ThreadVar(T initialValue = 0, ReductionMethod m=kSum)
	: fReductionMethod(m),
	  fValues(getMaxThreads(), initialValue),
	  fTouched(getMaxThreads(), false)
      {}

      inline void operator++(int){
	touched();
	fValues[thread()]++;
      }
      inline void operator--(int){
	touched();
	fValues[thread()]--;
      }
      inline void operator+=(const T& rhs){
	touched();
	fValues[thread()]+=rhs;
      }
      inline void operator-=(const T& rhs){
	touched();	
	fValues[thread()]-=rhs;
      }
      inline void operator*=(const T& rhs){
	touched();	
	fValues[thread()]*=rhs;
      }
      inline void operator/=(const T& rhs){
	touched();	
	fValues[thread()]/=rhs;
      }
      inline void operator=(const T& rhs){
	touched();	
	fValues[thread()]=rhs;
      }

      /** 
       * If you're in the pragma loop, return the thread local value
       * If you're out of the prama loop, return the reduced value
       * Is this a good idea?
       */
      operator T() const {
	if(getNumThreads() > 1){
	  return fValues[thread()];
	}
	else{
	  return reduction();
	}
      }

    private:

      inline void touched(){
	fTouched[thread()] = true;
      }

      /** 
       * Combine the thread local data to a single value
       * Exactly how that is done depends on the fReductionValue enum
       */
      T reduction() const {
	switch (fReductionMethod){
	  case kSum:		return sum();
	  case kNone: 	        return fValues[0];
	  case kProduct:	return product();
	  case kMax:		return max();
	  case kMin:		return min();	  	    	    
    	  default:		return sum();
	}
      }
      

      T sum() const {
	T rv = 0;
	for(unsigned i=0; i < fValues.size(); i++){
	  rv += fValues[i];
	}
	return rv;
      }

      T product() const {
	T rv = 1;
	for(unsigned i=0; i < fValues.size(); i++){
	  rv *= fValues[i];
	}
	return rv;
      }      


      T max() const {
	T mv = fValues[0];
	for(unsigned i=1; i < fValues.size(); i++){
	  if(fTouched[i] && fValues[i] > mv){
	    mv = fValues[i];
	  }
	}
	return mv;
      }      

      T min() const {
	T mv = fValues[0];
	for(unsigned i=1; i < fValues.size(); i++){
	  if(fTouched[i] && fValues[i] < mv){
	    mv = fValues[i];
	  }
	}
	return mv;
      }
      
      ReductionMethod fReductionMethod;
      std::vector<T> fValues;
      std::vector<bool> fTouched;
      
    };
  }

  
}

#endif //ACCLAIM_OPENMP_H
