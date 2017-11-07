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
    // void setNumThreads(int n);

    // const int default_num_threads = 8;
    
    template <class T>
    class ThreadVar {
    public:
      ThreadVar(T initialValue = 0) : values(getMaxThreads(), initialValue) {
	// setNumThreads(default_num_threads);
	std::cout << __PRETTY_FUNCTION__ << values.size() << ", " << initialValue << std::endl;	
      }

      void operator++(int){
	// std::cout << "thread " << thread() << std::endl;
	values[thread()]++;
      }
      void operator--(int){
	values[thread()]--;
      }      
      void operator+=(const T& rhs){
	values[thread()]+=rhs;
      }
      void operator-=(const T& rhs){
	values[thread()]-=rhs;
      }
      void operator*=(const T& rhs){
	values[thread()]*=rhs;
      }
      void operator/=(const T& rhs){
	values[thread()]/=rhs;
      }
      void operator=(const T& rhs){
	values[thread()]=rhs;
      }
      operator T() const {
	return sum();
      }      
      
      T sum() const {
	T sumOfVals = 0;
	for(unsigned i=0; i < values.size(); i++){
	  sumOfVals += values[i];
	}
	return sumOfVals;
      }

      private:
      std::vector<T> values;
    };
  }

  
}
