/*
 * grid.h
 *
 *  Created on: 20.06.2014
 *      Author: Paul
 *              Betim Musa <musab@informatik.uni-freiburg.de>
 */

#ifndef GRID_H_
#define GRID_H_

#include <utility>  // needed for std::pair
#include <vector>

#include "Phylogeny.h"
#include "RandomGen.h"
#include "Individual.h"

/**
 * The landscape is the main component of our world. It is separated into global environment and
 * local environment. When the species are globally independently (their reproduction doesn't depend
 * on the other species) distributed then we choose the global environment.
 * In the other case, when the distribution of the species is dependent to the local environmental
 * conditions, then we choose a local environment.
 */
class Landscape {
public:
    int m_Cutoff;
    int m_DensityCutoff;
    int m_Dispersal_type;
    // The size of the landscape.
    int m_Xdimensions, m_Ydimensions;
    // Counter for all the species that are born.
    unsigned long long m_Global_Species_Counter;
    unsigned int m_mortalityStrength;

    // Indicates whether the species dispersion depends on local conditions or not.
    bool m_Neutral;
    bool m_DD; // negative density dependence
    bool m_pDD; // positive density dependence
    bool m_Env;
    // switches to determine wheather fitness should affect mortality or reproduction
    bool m_mortality;
    bool m_reproduction;
    bool m_recruitment; // fitness acts on the recruit in the target cell
    static const int maxRecruitTries = 1000;

    // The size of the landscape is defined as m_Xdimensions * m_Ydimensions
    unsigned int m_LandscapeSize;

    // local dispersal: offsets of the cells within m_Cutoff (without the target cell), their dispersal weights
    // and cumulative weights, computed once; buffers for weights that depend on fitness
    std::vector<std::pair<int, int> > m_dispersalOffsets;
    std::vector<double> m_dispersalWeights, m_cumDispersalWeights;
    std::vector<double> m_fitnessWeights, m_cumFitnessWeights, m_cumRecruit;

    // offsets of the cells within the cutoff of negative [0] and positive [1] density dependence, without the
    // focal cell (empty if the mechanism is off); callModel() guarantees 2 * cutoff + 1 <= min(x, y)
    // for active mechanisms, so that no offset reaches the focal cell again
    std::vector<std::pair<int, int> > m_neighborOffsets[2];
    // cells whose stored densities may be outdated (index x * m_Ydimensions + y)
    std::vector<char> m_densityOutdated;

    double m_AirTemperature;
    double m_Speciation_Rate;

    double m_maxDeathRate;
    double m_nicheWidth;
    double m_envStrength;
    double m_compStrength;
    int m_fission;
    double m_redQueen;
    double m_redQueenStrength;
    int m_protracted;
    std::vector<double> airmat;
    std::vector<double> soilmat;
    int m_ddKernel;
    double m_ddWidth;
    int m_posCutoff;
    int m_posKernel;
    double m_posWidth;
    std::vector<double> m_ddTable;


    Landscape(int xsize, int ysize, int type, bool neutral,
              bool dd, bool env, bool mort, bool repro, unsigned int simulationEnd, double specRate,
              int dispersalCutoff, int densityCutoff, unsigned int mortalityStrength, double envStrength,
              double compStrength, int fission, double redQueen, double redQueenStrength, int protracted,
              std::vector<double> airmat, std::vector<double> soilmat, double nicheWidth,
              DensityOptions ddOptions);

    virtual ~Landscape();

    // TODO(Betim): Should it really be public?
    Phylogeny m_Phylogeny;

    // A 2-D array of individuals.
    Individual **m_Individuals;
    std::vector <std::pair<double, double> > m_Environment;
    RandomGen m_RandomGenerator;

    unsigned int m_SimulationEnd;

    // Start the reproduction of the species.
    virtual void reproduce(unsigned int generation);

    // Increase the age of each individual by 1.
    void increaseAge(unsigned int generation);

    void speciation(unsigned int generation);

    // marks the stored densities around cell (x, y) as outdated; must be called whenever the occupant, species or
    // competition marker of a cell changes, otherwise stored densities stay outdated without any error
    void densityUpdate(int x, int y);

    // recomputes the stored densities of cell (x, y) if they are outdated; call before reading them
    void ensureDensity(int x, int y);

    // stored local densities are read only when fitness acts on mortality or reproduction
    bool needsDensity() const { return (m_DD || m_pDD) && (m_mortality || m_reproduction); }

    // offsets (dx, dy) of the cells within the cutoff, ordered by dx, then dy
    static std::vector<std::pair<int, int> > circleOffsets(int cutoff, bool withCenter);

    // coordinate on the torus, for coordinates between -m_Xdimensions (-m_Ydimensions) and 2 * m_Xdimensions - 1
    int wrapX(int x) const { return x < 0 ? x + m_Xdimensions : (x >= m_Xdimensions ? x - m_Xdimensions : x); }
    int wrapY(int y) const { return y < 0 ? y + m_Ydimensions : (y >= m_Ydimensions ? y - m_Ydimensions : y); }

    // dissimilarity of ind to the neighbors of cell (x, y), without the cell itself, for negative or
    // positive density dependence
    double localDissimilarity(Individual &ind, int x, int y, bool positive);

    // sum over the neighbors in localDissimilarity(), for kernel 1 (absolute), 2 (gaussian) or 3 (species)
    template <int kernel>
    double dissimilaritySum(const Individual &ind, int x, int y, bool positive, double width) const;

    // fitness of the offspring of the individual in parent cell p in the target cell (x, y)
    double recruitFitness(int p_x, int p_y, int x, int y);

};

class GlobalEnvironment : public Landscape {
public:
    GlobalEnvironment(int xsize, int ysize, int type, bool neutral, bool dd, bool env, bool mort, bool repro,
                      unsigned int runs, double specRate, int dispersalCutoff, int densityCutoff,
                      unsigned int mortalityStrength,
                      double envStrength, double compStrength, int fission, double redQueen, double redQueenStrength,
                      int protracted, std::vector<double> airmat, std::vector<double> soilmat, double nicheWidth,
                      DensityOptions ddOptions);

    virtual ~GlobalEnvironment();

    void reproduce(unsigned int generation);

    // With the environment acting on recruitment, parents are also proposed from bins of their environmental
    // trait (binnedParent()), after uniformTriesBeforeBins rejected uniform proposals; empty if not used.
    static const int uniformTriesBeforeBins = 10;
    std::vector<std::vector<int> > m_traitBins; // cells by bin of the environmental trait of their individual
    std::vector<int> m_binOfCell, m_positionInBin;
    std::vector<double> m_binBound;  // upper bound of recruitFitness() for a distance of d bins to the target
    std::vector<double> m_cumBinWeights;

    int traitBin(double trait) const;
    void binCell(int cell); // moves the cell into the bin of the trait of its individual
    void binAllCells();
    // parent for the target cell with probability ~ recruitFitness(), or -1 after maxRecruitTries rejections
    int binnedParent(int target);
};

class LocalEnvironment : public Landscape {
public:
    LocalEnvironment(int xsize, int ysize, int type, bool neutral, bool dd, bool env, bool mort, bool repro,
                     unsigned int runs, double specRate, int dispersalCutoff, int densityCutoff,
                     unsigned int mortalityStrength,
                     double envStrength, double compStrength, int fission, double redQueen, double redQueenStrength,
                     int protracted, std::vector<double> airmat, std::vector<double> soilmat, double nicheWidth,
                     DensityOptions ddOptions);

    virtual ~LocalEnvironment();

    void reproduce(unsigned int generation);
};

#endif /* GRID_H_ */

