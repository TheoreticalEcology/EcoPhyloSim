/*
 * RandomGen.h
 *
 *  Created on: 10.07.2014
 *      Author: Paul
 *              Betim Musa <musab@informatik.uni-freiburg.de>
 */

#ifndef RANDOMGEN_H_
#define RANDOMGEN_H_

#include <vector>


class RandomGen {
public:
    RandomGen();

    ~RandomGen();

    int multinomialDraw(const std::vector<double> &cumWeights, int n);

    double randomDouble(double min, double max);

    int randomInt(int min, int max);

    int randomPoisson(double mean);
};

#endif /* RANDOMGEN_H_ */
