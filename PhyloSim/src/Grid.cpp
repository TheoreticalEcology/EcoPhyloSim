/*
 * gridMod.cpp
 *
 *  Created on: 12.01.2015
 *      Author: Paul
 *              Betim Musa <musab@informatik.uni-freiburg.de>
 */

#include <Rcpp.h>
#include <cfloat>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>
#include <numeric>
#include<algorithm>

#include "Grid.h"
#include "Individual.h"
#include "Species.h"
#include "debug.h"
#include "./RandomGen.h"

const int Landscape::maxRecruitTries; // definition, as it is passed by reference to std::min


Landscape::Landscape(int xsize, int ysize, int type, bool neutral, bool dd, bool env,
                     bool mort, bool repro, unsigned int runs, double specRate, int dispersalCutoff,
                     int DensityCutoff, unsigned int mortalityStrength,
                     double envStrength, double compStrength, int fission, double redQueen, double redQueenStrength,
                     int protracted,
                     std::vector<double> airmat, std::vector<double> soilmat, double nicheWidth,
                     DensityOptions ddOptions) {

    m_Cutoff = dispersalCutoff;
    m_Neutral = neutral;
    m_DD = dd;
    m_Env = env;
    m_mortality = mort;
    m_reproduction = repro;
    m_SimulationEnd = runs;
    m_DensityCutoff = DensityCutoff;
    m_Dispersal_type = type;
    m_Global_Species_Counter = 1;
    m_Speciation_Rate = specRate;
    m_LandscapeSize = xsize * ysize;
    m_Xdimensions = xsize;
    m_Ydimensions = ysize;
    m_mortalityStrength = mortalityStrength;
    m_envStrength = envStrength;
    m_compStrength = compStrength;
    m_fission = fission;
    m_redQueen = redQueen;
    m_redQueenStrength = redQueenStrength;
    m_protracted = protracted;
    m_nicheWidth = nicheWidth;
    m_recruitment = ddOptions.recruitment;
    m_ddKernel = ddOptions.kernel;
    m_ddWidth = ddOptions.width;
    m_ddTable = ddOptions.table;
    m_pDD = ddOptions.posStrength > 0;
    m_posCutoff = ddOptions.posCutoff;
    m_posKernel = ddOptions.posKernel;
    m_posWidth = ddOptions.posWidth;

    // Construct grid of individuals
    this->m_Individuals = new Individual *[m_Xdimensions];
    for (int cols = 0; cols < m_Xdimensions; cols++) {
        m_Individuals[cols] = new Individual[m_Ydimensions];
    }

    // Initialization

    // Create new initial species
    Species *spec = new Species(1, 1, 0, std::make_pair<int, int>(0, 0), m_SimulationEnd);

    // Add individuals of this new species across the grid
    for (int cols = 0; cols < m_Xdimensions; cols++) {
        for (int rows = 0; rows < m_Ydimensions; rows++) {
            this->m_Individuals[cols][rows].m_Species = spec;
            this->m_Individuals[cols][rows].m_X_coordinate = cols;
            this->m_Individuals[cols][rows].m_Y_coordinate = rows;
            this->m_Individuals[cols][rows].reportBirth();
            this->m_Individuals[cols][rows].m_dispersalDistance = m_Cutoff / 2.0;
            this->m_Individuals[cols][rows].m_envStrength = m_envStrength;
            this->m_Individuals[cols][rows].m_compStrength = m_compStrength;
            this->m_Individuals[cols][rows].m_nicheWidth = m_nicheWidth;
            this->m_Individuals[cols][rows].m_ddTable = m_ddTable.empty() ? NULL : &m_ddTable;
            this->m_Individuals[cols][rows].m_posStrength = ddOptions.posStrength;

        }
    }
    // TODO - should the species be updated?

    // Write the last individual in the phylogeny
    m_Phylogeny.updatePhylogeny(m_Individuals[xsize - 1][ysize - 1].m_Species);

    // Set up Environment
    // if no environment (airmat) is passed, generate one
    if ((airmat.size() == 1) && (soilmat.size() == 1)) {
        for (int i = 0; i < m_Xdimensions; i++) {
            // fixes #142 (gradient above 1 for odd x); unchanged for even x
            m_AirTemperature = 1.0 - std::fabs(2.0 * i / m_Xdimensions - 1.0);
            for (int j = 0; j < m_Ydimensions; j++) {
                std::pair<double, double> envi;
                envi.first = m_AirTemperature;
                envi.second = 1.0;
                m_Environment.push_back(envi);
            }
        }
    }

    // if environment (airmat) is passed use it
    if (airmat.size() > 1) {
        for (int i = 0; i < m_Xdimensions; i++) {
            for (int j = 0; j < m_Ydimensions; j++) {
                std::pair<double, double> envi;
                envi.first = airmat[i * m_Ydimensions + j]; // fixes #122 (airmat index for x != y)
                if (soilmat.size() > 1) {
                    envi.second = soilmat[i * m_Ydimensions + j];
                } else envi.second = 1.0;
                m_Environment.push_back(envi);
            }
        }
    }
    // END Set up Environment

    // Grid Geometry calculations
    for (int positive = 0; positive <= 1; positive++) {
        if (positive ? !m_pDD : !m_DD) continue;
        int cutoff = positive ? m_posCutoff : m_DensityCutoff;
        m_neighborOffsets[positive] = circleOffsets(cutoff, false);
    }
    m_densityOutdated.assign(m_LandscapeSize, 0);
    // END Grid Geometry calculations

    // largest product of the death rates of the mechanisms acting on mortality (Individual::deathProbability)
    double m = m_mortalityStrength;
    double minTable = 1.0;
    for (double t : m_ddTable) minTable = std::min(minTable, t);
    m_maxDeathRate = 1.0;
    if (m_Env) m_maxDeathRate *= 1 + (m - 1) * m_envStrength;
    if (m_DD) m_maxDeathRate *= 1 + (m - 1) * (m_ddTable.empty() ? m_compStrength : 1 - minTable);
    if (m_pDD) m_maxDeathRate *= 1 + (m - 1) * ddOptions.posStrength;

}

Landscape::~Landscape() {
    for (int i = 0; i < m_Xdimensions; i++) {
        for (int j = 0; j < m_Ydimensions; j++) {
            // TODO: FH 20.4.17 : I had suspected that the individuals need to be deleted here, but on reflection this seems fine?
            //delete &m_Individuals[i][j]; // delete stored pointer
        }
        delete[] m_Individuals[i];
    }
    delete[] m_Individuals;
}

void Landscape::reproduce(unsigned int /* generation */) {
    // (Betim): This method is needed in order to use polymorphism method calls.
    // (e.g. Model->updateEnvironment, the argument Landscape needs the method "reproduce"
    // to be defined, although it is void.)
}


void Landscape::increaseAge(unsigned int /* generation */) {
    for (int rows = 0; rows < this->m_Xdimensions; rows++) {
        for (int cols = 0; cols < this->m_Ydimensions; cols++) {
            this->m_Individuals[rows][cols].m_Age += 1;
        }
    }
}


GlobalEnvironment::GlobalEnvironment(int xsize, int ysize, int type, bool neutral, bool dd, bool env, bool mort,
                                     bool repro, unsigned int runs, double specRate, int dispersalCutoff,
                                     int densityCutoff, unsigned int mortalityStrength, double envStrength,
                                     double compStrength, int fission, double redQueen, double redQueenStrength,
                                     int protracted, std::vector<double> airmat, std::vector<double> soilmat, double nicheWidth,
                                     DensityOptions ddOptions) :
        Landscape(xsize, ysize, type, neutral, dd, env, mort, repro, runs, specRate, dispersalCutoff, densityCutoff,
                  mortalityStrength, envStrength, compStrength, fission, redQueen, redQueenStrength, protracted, airmat,
                  soilmat, nicheWidth, ddOptions) {

    // fixes #207 (up to N fitness evaluations per death while few individuals fit): parents are proposed from bins
    // of the environmental trait and accepted with fitness / bound, i.e. with the same probability as before;
    // the bound needs fitness terms between 0 and 1 (envStrength, compStrength in [0, 1]; posTerm, ddTable <= 1)
    if (m_recruitment && m_Env && !m_reproduction && envStrength >= 0 && envStrength <= 1 && compStrength >= 0 &&
        compStrength <= 1) {
        int nBins = (int) std::min(512.0, std::max(1.0, std::ceil(4.0 / nicheWidth)));
        double binWidth = 1.0 / nBins;
        for (int d = 0; d < nBins; d++) {
            // smallest distance between a trait in a bin d bins away and the environment of the target cell,
            // minus a margin for rounding in traitBin()
            double distance = std::max(0.0, (d - 1.0 - 1e-6) * binWidth);
            m_binBound.push_back((envStrength * exp(-0.5 * pow(distance / nicheWidth, 2.0)) + 1 - envStrength) *
                                 (1 + 1e-9)); // the other fitness terms are at most 1
        }
        m_traitBins.resize(nBins);
        m_cumBinWeights.resize(nBins);
    }
}

int GlobalEnvironment::traitBin(double trait) const {
    int nBins = m_traitBins.size();
    return std::min(nBins - 1, std::max(0, (int) (trait * nBins)));
}

void GlobalEnvironment::binCell(int cell) {
    if (m_binOfCell[cell] >= 0) { // remove it from its old bin
        std::vector<int> &old = m_traitBins[m_binOfCell[cell]];
        int last = old.back();
        old[m_positionInBin[cell]] = last;
        m_positionInBin[last] = m_positionInBin[cell];
        old.pop_back();
    }
    int bin = traitBin(m_Individuals[cell / m_Ydimensions][cell % m_Ydimensions].m_Mean);
    m_binOfCell[cell] = bin;
    m_positionInBin[cell] = m_traitBins[bin].size();
    m_traitBins[bin].push_back(cell);
}

void GlobalEnvironment::binAllCells() {
    for (auto &bin : m_traitBins) bin.clear();
    m_binOfCell.assign(m_LandscapeSize, -1);
    m_positionInBin.resize(m_LandscapeSize);
    for (unsigned int cell = 0; cell < m_LandscapeSize; cell++) binCell(cell);
}

int GlobalEnvironment::binnedParent(int target) {
    int x = target / m_Ydimensions, y = target % m_Ydimensions;
    int targetBin = traitBin(m_Environment[target].first);
    int nBins = m_traitBins.size();
    double sum = 0.0;
    for (int b = 0; b < nBins; b++) {
        sum += m_traitBins[b].size() * m_binBound[std::abs(b - targetBin)];
        m_cumBinWeights[b] = sum;
    }
    if (!(sum > 0)) return -1;
    for (int tries = 0; tries < maxRecruitTries; tries++) {
        int b = m_RandomGenerator.multinomialDraw(m_cumBinWeights, nBins);
        const std::vector<int> &cells = m_traitBins[b];
        if (cells.empty()) continue; // possible only through rounding of the cumulative weights
        int cell = cells[m_RandomGenerator.randomInt(0, cells.size() - 1)];
        if (cell == target) continue; // fixes #81
        double bound = m_binBound[std::abs(b - targetBin)];
        double fitness = recruitFitness(cell / m_Ydimensions, cell % m_Ydimensions, x, y);
        if (fitness > bound) Rcpp::stop("recruitment fitness above its bound in binnedParent(): fitness terms must "
                                        "be between 0 and 1 (internal error, please report)");
        if (m_RandomGenerator.randomDouble(0.0, 1.0) * bound < fitness) return cell;
    }
    return -1;
}

GlobalEnvironment::~GlobalEnvironment() {

}

void GlobalEnvironment::reproduce(unsigned int generation) {

    /////////////////////////////////////////////////////
    // NEUTRAL CASE
    if (m_Neutral && (m_redQueen == 0) && (m_redQueenStrength == 0)) {

#ifdef DEBUG
        Rcpp::Rcout << "In global neutral\n";
#endif

        int x_coordinate, y_coordinate, x_parent, y_parent;

        for (unsigned int event = 0; event < m_LandscapeSize; event++) {

            // choose random individual to die
            x_coordinate = m_RandomGenerator.randomInt(0, m_Xdimensions - 1);
            y_coordinate = m_RandomGenerator.randomInt(0, m_Ydimensions - 1);

            // choose random parent to "replace" individual
            do { // fixes #81 (dying individual must not be its own parent)
                x_parent = m_RandomGenerator.randomInt(0, m_Xdimensions - 1);
                y_parent = m_RandomGenerator.randomInt(0, m_Ydimensions - 1);
            } while (x_parent == x_coordinate && y_parent == y_coordinate);


            m_Individuals[x_coordinate][y_coordinate].reportDeath(generation);

            m_Individuals[x_coordinate][y_coordinate] = m_Individuals[x_parent][y_parent]; // overloaded operator, copy !
            m_Individuals[x_coordinate][y_coordinate].m_X_coordinate = x_coordinate;
            m_Individuals[x_coordinate][y_coordinate].m_Y_coordinate = y_coordinate;
            // report birth and evolve is NOT automatic
            m_Individuals[x_coordinate][y_coordinate].evolve(m_traitEvolution, m_traitPull);
        }
    }

    /////////////////////////////////////////////////////
    // NON-NEUTRAL CASE
    // Density dependence and / or Environmental dependence
    else
    {

#ifdef DEBUG
        Rcpp::Rcout << "In global non-neutral\n";
#endif

        int x_parent = 0;
        int y_parent = 0;

        unsigned int numberOfRuns = 0;
        int new_parent = 0;
        unsigned int array_length;

        std::vector<double> weights(m_LandscapeSize); // fixes #149 (heap instead of stack)
        std::vector<double> cumWeights(m_LandscapeSize);

        // all possible parent permutations
        std::vector <std::pair<int, int> > parents(m_LandscapeSize);
        int count = 0;
        for (int x = 0; x < m_Xdimensions; x++) {
            for (int y = 0; y < m_Ydimensions; y++) {
                parents[count].first = x;
                parents[count].second = y;
                count++;
            }
        }

        // IF FITNESS ACTS ON REPRODUCTION, calculate average fitness for lottery competition
        if (m_reproduction) {

#ifdef DEBUG
            Rcpp::Rcout << "In global non-neutral, reproduction fitness\n";
#endif

            array_length = 0;
            for (int kernel_x = 0; kernel_x < m_Xdimensions; kernel_x++) {
                for (int kernel_y = 0; kernel_y < m_Ydimensions; kernel_y++) {
                    ensureDensity(kernel_x, kernel_y);
                    weights[array_length] = m_Individuals[kernel_x][kernel_y].getFitness(
                            m_Environment[kernel_x * m_Ydimensions + kernel_y].first, m_Env, m_DD, generation,
                            m_redQueenStrength, m_redQueen);
                    array_length++;
                }
            }

            std::partial_sum(weights.begin(), weights.end(), cumWeights.begin());

        } // END  if (m_reproduction)

        // TODO could make numberofruns a parameter

        numberOfRuns = m_LandscapeSize;
        // births re-bin their own cell; re-binning all cells once per generation would also pick up trait changes
        // elsewhere, e.g. if speciation changes traits (#114)
        if (!m_traitBins.empty()) binAllCells();

        int event = 0;
        unsigned int numberDeath = 0;

        // TODO - seems to me the next for loop is not yet corrected if fitness goes on reproduction,

        while (numberDeath < numberOfRuns) {
            event++;

            int x_coordinate = m_RandomGenerator.randomInt(0, m_Xdimensions - 1);
            int y_coordinate = m_RandomGenerator.randomInt(0, m_Ydimensions - 1);

            // check if individual's fitness is high enough to survive
            // every m_mortalityStrength'th step the randomly chosen individual dies, no matter its fitness
            // with one mechanism, the old algorithm below gives the same ratios of death rates and keeps
            // the results of earlier versions
            if (m_mortality) ensureDensity(x_coordinate, y_coordinate);
            if (m_mortality && m_Env + m_DD + m_pDD > 1) { // proportional hazards, see Individual::deathProbability
                double p = m_Individuals[x_coordinate][y_coordinate].deathProbability(
                        m_Environment[x_coordinate * m_Ydimensions + y_coordinate].first, m_Env, m_DD,
                        m_mortalityStrength, generation, m_redQueenStrength, m_redQueen, m_maxDeathRate);
                if (m_RandomGenerator.randomDouble(0.0, 1.0) >= p) continue;
            } else if (m_mortality && event % m_mortalityStrength != 0) {
                double weight = m_Individuals[x_coordinate][y_coordinate].getFitness(
                        m_Environment[x_coordinate * m_Ydimensions + y_coordinate].first, m_Env, m_DD, generation,
                        m_redQueenStrength, m_redQueen);
                // important!! the frequency in relation to the base mortality controls the intensity of the mechanisms
                double chanceOfDeath = m_RandomGenerator.randomDouble(0.0, 1.0);
                if (weight > chanceOfDeath) {
                    continue;
                }
            }

            // continue with a real death here
            numberDeath++;

            m_Individuals[x_coordinate][y_coordinate].reportDeath(generation);

            int target = x_coordinate * m_Ydimensions + y_coordinate;
            auto drawParent = [&]() {
                int p;
                do { // fixes #81
                    if (m_reproduction) {
                        p = m_RandomGenerator.multinomialDraw(cumWeights, m_LandscapeSize);
                    } else {
                        p = m_RandomGenerator.randomInt(0, m_LandscapeSize - 1);
                    }
                } while (p == target);
                return p;
            };
            new_parent = drawParent();

            // rejection sampling: parent with probability ~ proposal * fitness of its recruit here;
            // after at most maxRecruitTries rejections, the same distribution is drawn from all weights
            for (int tries = 1; m_recruitment && m_RandomGenerator.randomDouble(0.0, 1.0) >=
                    recruitFitness(parents[new_parent].first, parents[new_parent].second, x_coordinate, y_coordinate);
                 tries++) {
                if (tries == uniformTriesBeforeBins && !m_traitBins.empty()) {
                    int p = binnedParent(target);
                    if (p >= 0) {
                        new_parent = p;
                        break;
                    }
                }
                if (tries < std::min(maxRecruitTries, (int) m_LandscapeSize - 1)) {
                    new_parent = drawParent();
                    continue;
                }
                std::vector<double> cumRecruit(m_LandscapeSize);
                double sum = 0.0;
                for (unsigned int i = 0; i < m_LandscapeSize; i++) {
                    if ((int) i != target) {
                        sum += (m_reproduction ? weights[i] : 1.0) *
                               (recruitFitness(parents[i].first, parents[i].second, x_coordinate, y_coordinate) +
                                DBL_MIN * 100.0);
                    }
                    cumRecruit[i] = sum;
                }
                // if all products underflow, keep the proposal (as in the local fallback)
                if (sum > 0) {
                    do {
                        new_parent = m_RandomGenerator.multinomialDraw(cumRecruit, m_LandscapeSize);
                    } while (new_parent == target);
                }
                break;
            }

            x_parent = parents[new_parent].first;
            y_parent = parents[new_parent].second;

            m_Individuals[x_coordinate][y_coordinate] = m_Individuals[x_parent][y_parent];  // overloaded operator, deep copy ! 


            // TODO: old comment was: 'report birth and evolve is automatic' -> this is wrong, neither evolve nor report birth are automatic. Evolve is below, what about report birth
            m_Individuals[x_coordinate][y_coordinate].m_X_coordinate = x_coordinate;
            m_Individuals[x_coordinate][y_coordinate].m_Y_coordinate = y_coordinate;

            m_Individuals[x_coordinate][y_coordinate].evolve(m_traitEvolution, m_traitPull);
            if (!m_traitBins.empty()) binCell(target);

            parents[x_coordinate * m_Ydimensions + y_coordinate].first = x_coordinate;
            parents[x_coordinate * m_Ydimensions + y_coordinate].second = y_coordinate;


            if (!m_reproduction && needsDensity()) densityUpdate(x_coordinate, y_coordinate);

            // only if env and dens affect reproduction
            if (m_reproduction) {

              Rcpp::Rcout << "REPRODUCTION AT THE MOMENT NOT IMPLEMENTED\n";
            } // END if (m_reproduction)
        } // END while (numberDeath < numberOfRuns)
    }
}


LocalEnvironment::LocalEnvironment(int xsize, int ysize, int type, bool neutral, bool dd, bool env, bool mort,
                                   bool repro, unsigned int runs, double specRate, int dispersalCutoff,
                                   int densityCutoff, unsigned int mortalityStrength, double envStrength,
                                   double compStrength, int fission, double redQueen, double redQueenStrength,
                                   int protracted, std::vector<double> airmat, std::vector<double> soilmat, double nicheWidth,
                                   DensityOptions ddOptions) :
        Landscape(xsize, ysize, type, neutral, dd, env, mort, repro, runs, specRate, dispersalCutoff, densityCutoff,
                  mortalityStrength, envStrength, compStrength, fission, redQueen, redQueenStrength, protracted, airmat,
                  soilmat, nicheWidth, ddOptions) {

    // dispersal weights depend only on the offset; computed once instead of for every death
    m_dispersalOffsets = circleOffsets(m_Cutoff, false);
    for (const auto &o : m_dispersalOffsets) {
        m_dispersalWeights.push_back(m_Individuals[0][0].dispersal(m_Dispersal_type,
                                                                   m_Individuals[0][0].euclidian_distance(o.first, o.second)));
    }
    int nKernel = m_dispersalOffsets.size();
    m_cumDispersalWeights.resize(nKernel);
    std::partial_sum(m_dispersalWeights.begin(), m_dispersalWeights.end(), m_cumDispersalWeights.begin());
    m_fitnessWeights.resize(nKernel);
    m_cumFitnessWeights.resize(nKernel);
    m_cumRecruit.resize(nKernel);
}

LocalEnvironment::~LocalEnvironment() {

}

void LocalEnvironment::reproduce(unsigned int generation) {

#ifdef DEBUG
    Rcpp::Rcout << "In local\n";
#endif

    int x_coordinate;
    int y_coordinate;
    int x_parent = 0;
    int y_parent = 0;
    int new_parent;
    unsigned int numberOfRuns = 0;
    int nKernel = m_dispersalOffsets.size();

    numberOfRuns = m_LandscapeSize;

    // TODO
    // What todo?

    int event = 0;
    unsigned int numberDeath = 0;

    while (numberDeath < numberOfRuns) {
        event++;

        // Take a new coordinate

        x_coordinate = m_RandomGenerator.randomInt(0, m_Xdimensions - 1);
        y_coordinate = m_RandomGenerator.randomInt(0, m_Ydimensions - 1);

        // If fitness acts on mortality, break here with a chance ~ fitness

        // important!! the frequency in relation to the base mortality controls the intensity of the mechanisms
        if (m_mortality) ensureDensity(x_coordinate, y_coordinate);
        if (m_mortality && m_Env + m_DD + m_pDD > 1) { // proportional hazards, see Individual::deathProbability
            double p = m_Individuals[x_coordinate][y_coordinate].deathProbability(
                    m_Environment[x_coordinate * m_Ydimensions + y_coordinate].first, m_Env, m_DD,
                    m_mortalityStrength, generation, m_redQueenStrength, m_redQueen, m_maxDeathRate);
            if (m_RandomGenerator.randomDouble(0.0, 1.0) >= p) continue;
        } else if (m_mortality && event % m_mortalityStrength != 0) {
            double weight = m_Individuals[x_coordinate][y_coordinate].getFitness(
                    m_Environment[x_coordinate * m_Ydimensions + y_coordinate].first, m_Env, m_DD, generation,
                    m_redQueenStrength, m_redQueen);
            double chanceOfDeath = m_RandomGenerator.randomDouble(0.0, 1.0);
            if (weight > chanceOfDeath) continue;
        }

        // If we continue, count up number of deaths and perform replacement

        numberDeath++;

        m_Individuals[x_coordinate][y_coordinate].reportDeath(generation);

        ////////////////////////////////////////////
        // DISPERSAL

        // parent cells; callModel() and checkPar() guarantee 2 * m_Cutoff + 1 <= min(x, y), so that no offset reaches
        // the target cell again
        auto parentX = [&](int i) { return wrapX(x_coordinate + m_dispersalOffsets[i].first); };
        auto parentY = [&](int i) { return wrapY(y_coordinate + m_dispersalOffsets[i].second); };

        if (m_reproduction) {
            for (int i = 0; i < nKernel; i++) {
                int p_x = parentX(i), p_y = parentY(i);
                ensureDensity(p_x, p_y);
                m_fitnessWeights[i] = m_Individuals[p_x][p_y].getSeedsTo(m_dispersalOffsets[i].first,
                                                                         m_dispersalOffsets[i].second, m_Dispersal_type,
                                                                         m_Environment[p_x * m_Ydimensions + p_y].first,
                                                                         m_Env, m_DD, generation, m_redQueenStrength,
                                                                         m_redQueen);
            }
            std::partial_sum(m_fitnessWeights.begin(), m_fitnessWeights.end(), m_cumFitnessWeights.begin());
        }
        const std::vector<double> &weights = m_reproduction ? m_fitnessWeights : m_dispersalWeights;
        const std::vector<double> &cumWeights = m_reproduction ? m_cumFitnessWeights : m_cumDispersalWeights;
        new_parent = m_RandomGenerator.multinomialDraw(cumWeights, nKernel);

        // fixes #165 (rejection sampling: parent with probability ~ weight * fitness of its recruit here);
        // after at most maxRecruitTries rejections, the same distribution is drawn from all weights
        for (int tries = 1; m_recruitment && m_RandomGenerator.randomDouble(0.0, 1.0) >=
                recruitFitness(parentX(new_parent), parentY(new_parent), x_coordinate, y_coordinate);
             tries++) {
            if (tries < std::min(maxRecruitTries, nKernel)) {
                new_parent = m_RandomGenerator.multinomialDraw(cumWeights, nKernel);
                continue;
            }
            double sum = 0.0;
            for (int i = 0; i < nKernel; i++) {
                sum += weights[i] * (recruitFitness(parentX(i), parentY(i), x_coordinate, y_coordinate) +
                                     DBL_MIN * 100.0);
                m_cumRecruit[i] = sum;
            }
            // if all products underflow (fitness about 0 for reproduction and recruitment), keep the proposal
            if (sum > 0) new_parent = m_RandomGenerator.multinomialDraw(m_cumRecruit, nKernel);
            break;
        }
        x_parent = parentX(new_parent);
        y_parent = parentY(new_parent);

        m_Individuals[x_coordinate][y_coordinate] = m_Individuals[x_parent][y_parent];  // overloaded operator, deep copy ! // report birth and evolve is automatic
        m_Individuals[x_coordinate][y_coordinate].m_X_coordinate = x_coordinate;
        m_Individuals[x_coordinate][y_coordinate].m_Y_coordinate = y_coordinate;

        m_Individuals[x_coordinate][y_coordinate].evolve(m_traitEvolution, m_traitPull);

        // UPDATE RELATEDNESS

        if (needsDensity()) densityUpdate(x_coordinate, y_coordinate);
    } // END while (numberDeath < numberOfRuns)
}

void Landscape::densityUpdate(int x, int y) {
    m_densityOutdated[x * m_Ydimensions + y] = 1;
    for (int positive = 0; positive <= 1; positive++) {
        for (const auto &o : m_neighborOffsets[positive]) {
            m_densityOutdated[wrapX(x + o.first) * m_Ydimensions + wrapY(y + o.second)] = 1;
        }
    }
}

// fixes #207 (densities recomputed only when read, not after every death): the stored value depends only on
// the current neighborhood, and every change in it marks the cell, so the value read is the same as before
void Landscape::ensureDensity(int x, int y) {
    char &outdated = m_densityOutdated[x * m_Ydimensions + y];
    Individual &ind = m_Individuals[x][y];
#ifdef ECOPHYLOSIM_CHECK_DENSITY
    // debug build (PKG_CPPFLAGS = -DECOPHYLOSIM_CHECK_DENSITY): a stored value that is not marked as outdated must
    // equal a fresh computation, otherwise a change of the neighborhood was not passed to densityUpdate()
    if (!outdated && ((m_DD && ind.m_LocalDensity != localDissimilarity(ind, x, y, false)) ||
                      (m_pDD && ind.m_LocalPositive != localDissimilarity(ind, x, y, true))))
        Rcpp::stop("stored density at cell (%d, %d) is outdated: densityUpdate() missing", x, y);
#endif
    if (!outdated) return;
    if (m_DD) ind.m_LocalDensity = localDissimilarity(ind, x, y, false);
    if (m_pDD) ind.m_LocalPositive = localDissimilarity(ind, x, y, true);
    outdated = 0;
}

std::vector<std::pair<int, int> > Landscape::circleOffsets(int cutoff, bool withCenter) {
    std::vector<std::pair<int, int> > offsets;
    for (int i = -cutoff; i <= cutoff; i++) {
        int yLim = floor(sqrt(cutoff * cutoff - i * i)); // avoid diagonal bias
        for (int j = -yLim; j <= yLim; j++) {
            if (withCenter || i != 0 || j != 0) offsets.push_back(std::make_pair(i, j));
        }
    }
    return offsets;
}

template <int kernel>
double Landscape::dissimilaritySum(const Individual &ind, int x, int y, bool positive, double width) const {
    double sum = 0.0;
    for (const auto &o : m_neighborOffsets[positive]) {
        const Individual &neighbor = m_Individuals[wrapX(x + o.first)][wrapY(y + o.second)];
        double diff = std::abs(ind.m_CompetitionMarker - neighbor.m_CompetitionMarker);
        if (kernel == 1) sum += diff;
        else if (kernel == 2) sum += exp(-0.5 * pow(diff / width, 2.0));
        else sum += (ind.m_Species == neighbor.m_Species); // each species ID has one Species object
    }
    return sum;
}

double Landscape::localDissimilarity(Individual &ind, int x, int y, bool positive) {
    int kernel = positive ? m_posKernel : m_ddKernel;
    double width = positive ? m_posWidth : m_ddWidth;
    // kernel as template parameter, so that the loop has no branch per neighbor
    double sum = kernel == 1 ? dissimilaritySum<1>(ind, x, y, positive, width)
               : kernel == 2 ? dissimilaritySum<2>(ind, x, y, positive, width)
               : dissimilaritySum<3>(ind, x, y, positive, width);
    double mean = sum / m_neighborOffsets[positive].size(); // mean over the neighbors, without the focal cell
    return kernel == 1 ? mean : 1.0 - mean;
}

double Landscape::recruitFitness(int p_x, int p_y, int x, int y) {
    Individual &parent = m_Individuals[p_x][p_y];
    double out = 1.0;
    if (m_Env) out *= parent.envTerm(m_Environment[x * m_Ydimensions + y].first);
    if (m_DD) out *= parent.ddTerm(localDissimilarity(parent, x, y, false));
    if (m_pDD) out *= parent.posTerm(localDissimilarity(parent, x, y, true));
    return out;
}


void Landscape::speciation(unsigned int generation) {

    // In this function is the main part of the implementation of the speciation mechanisms.

    int specRate = m_RandomGenerator.randomPoisson(m_Speciation_Rate);

    // trait jump of a new species, shared by all its individuals (#114); drawn only for traits with a
    // jump, so that results without a jump stay unchanged
    double jump[3];
    auto drawJump = [&]() {
        for (int t = 0; t < 3; t++)
            jump[t] = m_specJump[t] > 0 ? m_RandomGenerator.randomDouble(-m_specJump[t], m_specJump[t]) : 0.0;
    };

    // Each speciation event starts an incipient lineage with its own tag. Offspring inherit the tag of their
    // parent, so that the whole lineage completes speciation together, `protracted` generations after the
    // event (#143); an event on an incipient individual starts a lineage nested in its own.
    auto newTag = [&](unsigned long long parent, unsigned int start) {
        m_incipient[++m_incipEvents] = {parent, start};
        return m_incipEvents;
    };
    auto outermost = [&](unsigned long long t) {
        while (m_incipient.at(t).parent != 0) t = m_incipient.at(t).parent;
        return t;
    };
    auto carries = [&](unsigned long long t, unsigned long long lineage) {
        while (t != 0 && t != lineage) t = m_incipient.at(t).parent;
        return t != 0;
    };

    // the incipient lineage `lineage` becomes a new species, with the individual at (k, j) as its founder
    auto complete = [&](int k, int j, unsigned long long lineage) {
        Individual &founder = m_Individuals[k][j];
        m_Global_Species_Counter += 1;
        founder.m_Species->m_Children.push_back(m_Global_Species_Counter);
        founder.reportDeath(generation);
        founder.m_Species = new Species(m_Global_Species_Counter, founder.m_Species->get_species_ID(),
                                        generation, std::make_pair(k, j), m_SimulationEnd);
        founder.m_Species->m_Date_of_Split = m_incipient.at(lineage).start;

        drawJump();
        founder.evolveDuringSpeciation(jump);
        m_Phylogeny.updatePhylogeny(founder.m_Species);
        if (needsDensity()) densityUpdate(k, j);
        if (founder.m_incipTag == lineage) founder.m_incipTag = 0;

        // the other members of the lineage, including those of lineages nested in it, which keep their own tag;
        // without protracted speciation and fission, the founder is the only member. The founder is handled
        // first, so that species trait values and results with protracted = 0 stay the same
        if (m_protracted > 0 || m_fission != 0) {
            for (int z = 0; z < m_Xdimensions; z++) {
                for (int q = 0; q < m_Ydimensions; q++) {
                    Individual &ind = m_Individuals[z][q];
                    if (&ind == &founder || !carries(ind.m_incipTag, lineage)) continue;
                    ind.reportDeath(generation);
                    ind.m_Species = founder.m_Species;
                    ind.evolveDuringSpeciation(jump);
                    if (needsDensity()) densityUpdate(z, q);
                    if (ind.m_incipTag == lineage) ind.m_incipTag = 0;
                }
            }
        }
        for (auto &l : m_incipient) if (l.second.parent == lineage) l.second.parent = 0;
        m_incipient.erase(lineage);
    };

    for (int i = 0; i < specRate; i++) {

        int x = m_RandomGenerator.randomInt(0, m_Xdimensions - 1); // rand() % xdimensions;
        int y = m_RandomGenerator.randomInt(0, m_Ydimensions - 1); // rand() % ydimensions;

#ifdef DEBUG
        Rcpp::Rcout << "new species at:\nx: " << x << "; y: " << y << "\n";
#endif

        // without fission, the hit individual starts an incipient lineage; a second event on it in the same
        // generation still gives one species (#208)
        if (m_fission == 0) {
            Individual &ind = m_Individuals[x][y];
            if (ind.m_incipTag == 0 || m_incipient.at(ind.m_incipTag).start != generation)
                ind.m_incipTag = newTag(ind.m_incipTag, generation);
        }

        // fixes #208 (lost events): an event splits only the group of the hit individual, i.e. the
        // individuals of its species in the same incipient lineage, or those not incipient
        if (m_fission != 0) {
            unsigned long long hitGroup = m_Individuals[x][y].m_incipTag;
            std::vector<std::pair<int, int> > cells; // in grid order
            for (int k = 0; k < m_Xdimensions; k++) {
                for (int j = 0; j < m_Ydimensions; j++) {
                    if (m_Individuals[k][j].m_Species == m_Individuals[x][y].m_Species &&
                        m_Individuals[k][j].m_incipTag == hitGroup)
                        cells.push_back(std::make_pair(k, j));
                }
            }

            // fission 2: the half with the smaller coordinates along the axis with the larger extent
            if (m_fission == 2 && cells.size() > 1) {
                int key[2][2]; // per axis: span of the group and first coordinate after its largest gap
                for (int axis = 0; axis < 2; axis++) {
                    int dim = axis ? m_Ydimensions : m_Xdimensions;
                    std::vector<char> occupied(dim, 0);
                    for (const auto &c : cells) occupied[axis ? c.second : c.first] = 1;
                    // largest run of empty coordinates on the torus; ties keep the gap across the edge
                    int first = 0, last = dim - 1;
                    while (!occupied[first]) first++;
                    while (!occupied[last]) last--;
                    int gap = first + dim - 1 - last, start = first;
                    for (int v = first, prev = first; v <= last; v++) {
                        if (!occupied[v]) continue;
                        if (v - prev - 1 > gap) {
                            gap = v - prev - 1;
                            start = v;
                        }
                        prev = v;
                    }
                    key[axis][0] = dim - gap;
                    key[axis][1] = start;
                }
                int axis = key[1][0] > key[0][0];
                int dim = axis ? m_Ydimensions : m_Xdimensions, start = key[axis][1];
                auto pos = [&](const std::pair<int, int> &c) { return ((axis ? c.second : c.first) - start + dim) % dim; };
                std::stable_sort(cells.begin(), cells.end(),
                                 [&](const std::pair<int, int> &a, const std::pair<int, int> &b) { return pos(a) < pos(b); });
            }

            // the first, third, ... individuals (fission 1) or the first half (fission 2) form the new species;
            // a group of one individual is replaced as a whole; a split incipient lineage keeps its start
            unsigned long long tag = hitGroup == 0 ? newTag(0, generation)
                                                   : newTag(0, m_incipient.at(hitGroup).start);
            for (unsigned int z = 0; z < cells.size(); z++) {
                if (m_fission == 1 ? z % 2 : 2 * z >= cells.size()) continue;
                m_Individuals[cells[z].first][cells[z].second].m_incipTag = tag;
            }
            // without protracted speciation, the next event already sees the new species
            if (m_protracted == 0) complete(cells[0].first, cells[0].second, tag);
        }
    } // END  for (int i = 0; i < specRate; i++)

    // Lineages whose speciation started `protracted` generations ago become new species, in grid order of
    // their first member; nested lineages after the lineage they are nested in
    for (int k = 0; k < m_Xdimensions; k++) {
        for (int j = 0; j < m_Ydimensions; j++) {
            while (m_Individuals[k][j].m_incipTag != 0) {
                unsigned long long lineage = outermost(m_Individuals[k][j].m_incipTag);
                if (generation - m_incipient.at(lineage).start < (unsigned int) m_protracted) break;
                complete(k, j, lineage);
            }
        }
    }
    // lineages that died out before completing speciation (after a fission split, tags are not ordered by start)
    for (auto it = m_incipient.begin(); it != m_incipient.end();) {
        if (generation - it->second.start >= (unsigned int) m_protracted) it = m_incipient.erase(it);
        else ++it;
    }
}
        



