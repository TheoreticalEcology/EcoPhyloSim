/*
 * individual.h
 *
 *  Created on: 20.06.2014
 *      Author: Paul
 *              Betim Musa <musab@informatik.uni-freiburg.de>
 */

#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

#include "./RandomGen.h"
//#include "Utils/RandomGen.h"
#include "Species.h"

#include <vector>

/**
 * Options for density dependence and recruitment.
 * kernel: 1 = absolute difference of the competition trait, 2 = Gaussian kernel of width
 * 'width', 3 = species identity. table: density term as a function of the similarity S to the
 * neighbors on an equally spaced grid over [0, 1]; empty for the linear response.
 */
struct DensityOptions {
    bool recruitment = false;
    int kernel = 1;
    double width = 0.1;
    std::vector<double> table;
    double posStrength = 0; // positive density dependence, with its own radius, kernel and width
    int posCutoff = 1;
    int posKernel = 3;
    double posWidth = 0.1;
};

/**
 * An individual is a realization of a species with individual properties.
 */
class Individual {
private:
    double kernel(double distance);

public:
    Individual();

    Individual(const Individual &ind);

    ~Individual();

    void operator=(const Individual &ind);

    // Compute the euclidian_distance from the origin (0, 0) to the given position (ref_x, ref_y)
    double euclidian_distance(int ref_x, int ref_y);

    // The species to which this individual belongs.
    Species *m_Species;

    // The position of this individual.
    int m_X_coordinate;
    int m_Y_coordinate;

    double m_LocalDensity;
    unsigned long long m_Age;
    int m_incip_Age;

//	double m_DensityStrength;
//	double m_FitnessWeight;
//	int m_Weight;

    double m_nicheWidth;
    double m_Mean;
    double m_CompetitionMarker;
    double m_NeutralMarker;
    double m_dispersalDistance;

    double m_envStrength;
    double m_compStrength;
    const std::vector<double> *m_ddTable; // response table of DensityOptions, NULL for the linear response
    double m_posStrength; // strength of positive density dependence
    double m_LocalPositive; // dissimilarity to the neighbors for positive density dependence

    RandomGen m_RandomGenerator;

    void reportDeath(int generation);

    void reportBirth(); // should make this private later, but at the moment used in init

    void evolve();

    void evolveDuringSpeciation();

    double getSeedsTo(int rel_x, int rel_y, int dispersal_type, double temp, bool env, bool dd, int generation,
                      double redQueenStrength, double redQueen);

    double getFitness(double temp, bool env, bool dd, int generation, double redQueenStrength, double redQueen);

    // fitness terms of the environment and of density dependence, each in [0, 1]
    double envTerm(double temp);

    double ddTerm(double localDensity);

    double posTerm(double localPositive);

    // probability to die when selected, if several mechanisms act on mortality
    double deathProbability(double temp, bool env, bool dd, unsigned int mortalityStrength, int generation,
                            double redQueenStrength, double redQueen, double maxRate);


    double dispersal(int dispersal_type, double distance);  // 1 for kernel, 2 for nearest neighbor, 3 for global

    void printInfo();
};

#endif /* INDIVIDUAL_H_ */
