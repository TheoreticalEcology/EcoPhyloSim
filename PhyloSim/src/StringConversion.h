/*
 * StringConversion.h
 *
 *  Created on: 30.07.2015
 *      Author: Betim Musa <musab@informatik.uni-freiburg.de>
 */

#include <sstream>
#include <string>

#ifndef STRINGCONVERSION_H_
#define STRINGCONVERSION_H_

template<typename T>
inline std::string to_string(T number) {
    std::stringstream ss;
    ss << number;
    return ss.str();
}

#endif /* STRINGCONVERSION_H_ */
