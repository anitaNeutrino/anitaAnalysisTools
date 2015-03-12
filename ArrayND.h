/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Makes c++11 matrices easier to read!
	     To compile on older machines I have switched to using boost library implementation of array, 
	     now installed on my Mac.
	     You don't actually need this, it's just for fun.
*************************************************************************************************************** */

#ifndef __CINT__
#ifndef ARRAYND_H
#define ARRAYND_H

// #include <array>
#include <boost/array.hpp>

template <class T1, size_t dim11>
// using Array1D = std::array<T1, dim11>;
using Array1D = boost::array<T1, dim11>;
//#define at((a), (b)) at(a).at(b) /* Is this a bad idea? */

template <class T2, size_t dim21, size_t dim22>
// using Array2D = std::array<std::array<T2, dim22>, dim21>;
using Array2D = boost::array<boost::array<T2, dim22>, dim21>;
#define at2(a, b) at(a).at(b) /* Is this a bad idea? */


template <class T3, size_t dim31, size_t dim32, size_t dim33>
// using Array3D = std::array<std::array<std::array<T3, dim33>, dim32>, dim31>;
using Array3D = boost::array<boost::array<boost::array<T3, dim33>, dim32>, dim31>;
#define at3(a,b,c) at(a).at(b).at(c) /* Is this a bad idea? */

#endif //ARRAYND_H
#endif //__CINT__
