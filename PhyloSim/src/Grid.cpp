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
        m_densityOffsets[positive] = circleOffsets(cutoff, true);
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
            this->m_Individuals[rows][cols].m_incip_Age += 1;

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
        if (m_RandomGenerator.randomDouble(0.0, 1.0) * m_binBound[std::abs(b - targetBin)] <
            recruitFitness(cell / m_Ydimensions, cell % m_Ydimensions, x, y)) return cell;
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
            m_Individuals[x_coordinate][y_coordinate].evolve();
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

            m_Individuals[x_coordinate][y_coordinate].evolve();
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

    // fixes #207 (dispersal weights depend only on the offset; computed once instead of for every death)
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

        m_Individuals[x_coordinate][y_coordinate].evolve();

        // UPDATE RELATEDNESS

        if (needsDensity()) densityUpdate(x_coordinate, y_coordinate);
    } // END while (numberDeath < numberOfRuns)
}

void Landscape::densityUpdate(int x, int y) {
    for (int positive = 0; positive <= 1; positive++) {
        for (const auto &o : m_densityOffsets[positive]) {
            m_densityOutdated[wrapX(x + o.first) * m_Ydimensions + wrapY(y + o.second)] = 1;
        }
    }
}

// fixes #207 (densities recomputed only when read, not after every death): the stored value depends only on
// the current neighborhood, and every change in it marks the cell, so the value read is the same as before
void Landscape::ensureDensity(int x, int y) {
    char &outdated = m_densityOutdated[x * m_Ydimensions + y];
    if (!outdated) return;
    Individual &ind = m_Individuals[x][y];
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
    // fixes #207 (kernel as template parameter, so that the loop has no branch per neighbor)
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


    for (int i = 0; i < specRate; i++) {

        int x = m_RandomGenerator.randomInt(0, m_Xdimensions - 1); // rand() % xdimensions;
        int y = m_RandomGenerator.randomInt(0, m_Ydimensions - 1); // rand() % ydimensions;

#ifdef DEBUG
        Rcpp::Rcout << "new species at:\nx: " << x << "; y: " << y << "\n";
#endif

        // In every case the counter for the incipient ages set to zero.
        // This value is updated every time step. Which individuals are chosen is dependent on the
        // other speciation mechanisms (fission or no fission).

        if (m_fission == 0) {
            m_Individuals[x][y].m_incip_Age = 0;
        }

        if (m_fission != 0) {
            std::vector<int> xvec;
            std::vector<int> yvec;
            for (int k = 0; k < m_Xdimensions; k++) {
                for (int j = 0; j < m_Ydimensions; j++) {
                    if (m_Individuals[k][j].m_Species == m_Individuals[x][y].m_Species) {
                        xvec.push_back(k);
                        yvec.push_back(j);
                    }
                }
            }


            if (m_fission == 1) {
                for (unsigned int z = 0; z < xvec.size(); z += 2) {
                    m_Individuals[xvec[z]][yvec[z]].m_incip_Age = 0;
                }
            }

            if (m_fission == 2) {
                int sum = std::accumulate(xvec.begin(), xvec.end(), 0);
                int mean = sum / xvec.size();
                if (sum == 0) mean = -1;

                for (unsigned int z = 0; z < xvec.size(); z++) {
                    if (xvec[z] < mean) {
                        m_Individuals[xvec[z]][yvec[z]].m_incip_Age = 0;
                    }
                }
            }
        }
    } // END  for (int i = 0; i < specRate; i++)

    // In the same timestep new species arise.
    // Which individuals become part of the new species is dependent on the incipient age of the
    // individual and the other speciation mechanism (fission or no fission).

    if (m_fission == 0) {
        for (int k = 0; k < m_Xdimensions; k++) {
            for (int j = 0; j < m_Ydimensions; j++) {
                if (m_Individuals[k][j].m_incip_Age == m_protracted) {
                    m_Global_Species_Counter += 1;

                    m_Individuals[k][j].m_Species->m_Children.push_back(m_Global_Species_Counter);
                    m_Individuals[k][j].reportDeath(generation);
                    m_Individuals[k][j].m_Species = new Species(m_Global_Species_Counter,
                                                                m_Individuals[k][j].m_Species->get_species_ID(),
                                                                generation, std::make_pair(k, j), m_SimulationEnd);
                    m_Individuals[k][j].evolveDuringSpeciation();
                    m_Phylogeny.updatePhylogeny(m_Individuals[k][j].m_Species);

                    if (needsDensity()) densityUpdate(k, j);
                }
            }
        }
    }

    if (m_fission != 0) {
        unsigned long long oldspec;
        for (int k = 0; k < m_Xdimensions; k++) {
            for (int j = 0; j < m_Ydimensions; j++) {
                if (m_Individuals[k][j].m_incip_Age == m_protracted) {
                    m_Global_Species_Counter += 1;
                    m_Individuals[k][j].m_Species->m_Children.push_back(m_Global_Species_Counter);
                    m_Individuals[k][j].reportDeath(generation);

                    oldspec = m_Individuals[k][j].m_Species->m_ID;
                    m_Individuals[k][j].m_Species = new Species(m_Global_Species_Counter,
                                                                m_Individuals[k][j].m_Species->get_species_ID(),
                                                                generation, std::make_pair(k, j), m_SimulationEnd);

                    m_Individuals[k][j].evolveDuringSpeciation();
                    m_Phylogeny.updatePhylogeny(m_Individuals[k][j].m_Species);
                    if (needsDensity()) densityUpdate(k, j);
                    m_Individuals[k][j].m_incip_Age = -99999999;

                    for (int z = 0; z < m_Xdimensions; z++) {
                        for (int q = 0; q < m_Ydimensions; q++) {
                            if ((m_Individuals[z][q].m_incip_Age == m_protracted) && (m_Individuals[z][q].m_Species->m_ID == oldspec)) {
                                m_Individuals[z][q].reportDeath(generation);
                                m_Individuals[z][q].m_Species = m_Individuals[k][j].m_Species;
                                m_Individuals[z][q].evolveDuringSpeciation();
                                if (needsDensity()) densityUpdate(z, q);
                                m_Individuals[z][q].m_incip_Age = -99999999;
                            }
                        }
                    }
                }
            }
        }
    }
}
        



