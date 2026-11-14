/*
 * RandomGen.cpp
 *
 *  Created on: 08.07.2014
 *      Author: Paul
 */

#include "RandomGen.h"

#include <algorithm>
#include <Rcpp.h>

// Random numbers come from R's generator, so that results for a seed do not depend on the compiler
// (fixes #14); "Writing R Extensions" advises against the C++ <random> library, whose engines and
// distributions are implementation-defined.

RandomGen::RandomGen() {

}

RandomGen::~RandomGen() {

}

double RandomGen::randomDouble(double min, double max) {
    return min + (max - min) * unif_rand();
}


int RandomGen::randomInt(int min, int max) {
    return min + (int) (unif_rand() * (max - min + 1.0)); // bias at most (max - min + 1) / 2^32, faster than R_unif_index()
}

int RandomGen::randomPoisson(double mean) {
    return (int) R::rpois(mean);
}

// Draws an index i in [0, n) with probability w_i / sum(w), where cumWeights[i] = w_0 + ... + w_i
// (e.g. from std::partial_sum) and the sum is cumWeights[n - 1].
int RandomGen::multinomialDraw(const std::vector<double> &cumWeights, int n) {
    double x = randomDouble(0.0, cumWeights[n - 1]);
    int i = std::upper_bound(cumWeights.begin(), cumWeights.begin() + n, x) - cumWeights.begin();
    return std::min(i, n - 1); // stays in range if the sum is 0, Inf or NaN
}
