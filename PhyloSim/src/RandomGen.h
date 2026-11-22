/*
 * RandomGen.h
 *
 *  Created on: 10.07.2014
 *      Author: Paul
 *              Betim Musa <musab@informatik.uni-freiburg.de>
 */

#ifndef RANDOMGEN_H_
#define RANDOMGEN_H_

#include <cstdint>
#include <vector>

// Copy of R's Mersenne-Twister (unif_rand() with RNGkind("Mersenne-Twister")): the same numbers, but inline
// (fixes #207: half of the time of global neutral runs went into unif_rand()). The state is shared by all
// RandomGen objects and synchronized with R's .Random.seed by RNGSync; RandomGen draws give R's numbers only
// while an RNGSync object exists (in callModel(), multinomialDrawCounts() and binnedParentCounts()).
namespace mersenneTwister {
    extern uint32_t state[625]; // as R's dummy[]: state[0] = position in the block, state[1..624] = block
    extern double block[624];   // numbers of the current block
    extern int position;        // next number in block
    void nextBlock();

    inline double unif() {
        if (position >= 624) nextBlock();
        return block[position++];
    }
}

// While an RNGSync object exists, random numbers come from the copy; R's state is taken over at its
// construction and written back at its destruction (also when an error unwinds the stack)
class RNGSync {
public:
    RNGSync() { load(); }
    ~RNGSync() { save(); }
    static void load();
    static void save();
};

class RandomGen {
public:
    RandomGen();

    ~RandomGen();

    int multinomialDraw(const std::vector<double> &cumWeights, int n);

    double randomDouble(double min, double max) { return min + (max - min) * mersenneTwister::unif(); }

    // bias at most (max - min + 1) / 2^32, faster than R_unif_index()
    int randomInt(int min, int max) { return min + (int) (mersenneTwister::unif() * (max - min + 1.0)); }

    int randomPoisson(double mean);
};

#endif /* RANDOMGEN_H_ */
