/*
 * RandomGen.cpp
 *
 *  Created on: 08.07.2014
 *      Author: Paul
 */

#include "RandomGen.h"

#include <algorithm>
#include <Rcpp.h>

// Random numbers come from R's generator (a copy of its Mersenne-Twister, see RandomGen.h), so that results
// for a seed do not depend on the compiler (fixes #14); "Writing R Extensions" advises against the C++
// <random> library, whose default engine and distributions are implementation-defined.

RandomGen::RandomGen() {

}

RandomGen::~RandomGen() {

}

// R's rpois() uses R's own generator: hand the state over to R and take it back
int RandomGen::randomPoisson(double mean) {
    RNGSync::save();
    int n = (int) R::rpois(mean);
    RNGSync::load();
    return n;
}

// The functions toUnif() and nextBlock() are adapted from R's src/main/RNG.c (GPL (>= 2)):
//   Copyright (C) 1997--2026 The R Core Team
//   Copyright (C) 1995, 1996 Robert Gentleman and Ross Ihaka
//   MT19937: Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura; M. Matsumoto and T. Nishimura,
//   "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform Pseudo-Random Number Generator", ACM
//   Transactions on Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, pp 3--30.
namespace mersenneTwister {
    uint32_t state[625];
    double block[624];
    int position = 624;

    // tempering and conversion to (0, 1) as MT_genrand() and fixup() in R's src/main/RNG.c
    static double toUnif(uint32_t y) {
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= (y >> 18);
        double u = (double) y * 2.3283064365386963e-10;
        return u <= 0.0 ? 0.5 * 2.328306437080797e-10 : u; // u < 1 always
    }

    // next block of 624 numbers, as MT_genrand() in R's src/main/RNG.c
    void nextBlock() {
        static const uint32_t mag01[2] = {0x0, 0x9908b0df};
        uint32_t *mt = state + 1, y;
        int k;
        for (k = 0; k < 624 - 397; k++) {
            y = (mt[k] & 0x80000000) | (mt[k + 1] & 0x7fffffff);
            mt[k] = mt[k + 397] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (; k < 623; k++) {
            y = (mt[k] & 0x80000000) | (mt[k + 1] & 0x7fffffff);
            mt[k] = mt[k + (397 - 624)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[623] & 0x80000000) | (mt[0] & 0x7fffffff);
        mt[623] = mt[396] ^ (y >> 1) ^ mag01[y & 0x1];
        for (int i = 0; i < 624; i++) block[i] = toUnif(mt[i]);
        position = 0;
    }
}

void RNGSync::load() {
    PutRNGstate(); // R's state -> .Random.seed
    Rcpp::IntegerVector seed = Rcpp::Environment::global_env()[".Random.seed"];
    if (seed.size() != 626 || seed[0] % 100 != 3) Rcpp::stop("the simulation needs RNGkind(\"Mersenne-Twister\")");
    for (int i = 0; i < 625; i++) mersenneTwister::state[i] = (uint32_t) seed[i + 1];
    if (mersenneTwister::state[0] > 624) Rcpp::stop("invalid .Random.seed"); // also negative values
    mersenneTwister::position = mersenneTwister::state[0];
    for (int i = mersenneTwister::position; i < 624; i++) {
        mersenneTwister::block[i] = mersenneTwister::toUnif(mersenneTwister::state[i + 1]);
    }
}

void RNGSync::save() {
    Rcpp::Environment global = Rcpp::Environment::global_env();
    Rcpp::IntegerVector old = global[".Random.seed"];
    Rcpp::IntegerVector seed(626); // a new vector, so that copies of the old one (e.g. to restore it) stay unchanged
    seed[0] = old[0];
    mersenneTwister::state[0] = mersenneTwister::position;
    for (int i = 0; i < 625; i++) seed[i + 1] = (int) mersenneTwister::state[i];
    global[".Random.seed"] = seed;
    GetRNGstate(); // .Random.seed -> R's state
}

// Draws an index i in [0, n) with probability w_i / sum(w), where cumWeights[i] = w_0 + ... + w_i
// (e.g. from std::partial_sum) and the sum is cumWeights[n - 1].
int RandomGen::multinomialDraw(const std::vector<double> &cumWeights, int n) {
    double x = randomDouble(0.0, cumWeights[n - 1]);
    int i = std::upper_bound(cumWeights.begin(), cumWeights.begin() + n, x) - cumWeights.begin();
    return std::min(i, n - 1); // stays in range if the sum is 0, Inf or NaN
}
