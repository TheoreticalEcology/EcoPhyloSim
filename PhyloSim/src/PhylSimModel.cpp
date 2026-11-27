/*
 * model.c++
 *
 *  Created on: 20.06.2014
 *      Author: Paul
 *              Betim Musa <musab@informatik.uni-freiburg.de>
 */

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include "./RandomGen.h"
#include "Individual.h"
#include "Grid.h"
#include "PhylSimModel.h"
#include "Species.h"


PhylSimModel::PhylSimModel(int X, int Y, int dispersal, int simulationEnd, double specRate, bool dens,
                           bool env, bool neutral, bool mort, int mortalityStrength, bool repro, double dispersalCutoff,
                           int densityCutoff, double envStrength, double compStrength,
                           int fission,
                           double redQueen, double redQueenStrength, int protracted, std::vector<double> airmat,
                           std::vector<double> soilmat, double nicheWidth, DensityOptions ddOptions) {

#ifdef DEBUG
    Rcpp::Rcout << "Running simulation with\n";
    Rcpp::Rcout << "Dispersal " << dispersal << "; cutoff " << dispersalCutoff << "\n";
    Rcpp::Rcout << "Competition " << dens << "; strength " << compStrength << " cutoff " << densityCutoff << "\n";
    Rcpp::Rcout << "Environment " << env << "; strength " << envStrength << "\n";
    Rcpp::Rcout << "\n---- debug message for development purposes, remove debug switch in debug.h for turning this off\n\n";
#endif

    if (dispersal == 1) {
        m_Global = new GlobalEnvironment(X, Y, dispersal, neutral, dens, env, mort, repro, simulationEnd, specRate,
                                         dispersalCutoff, densityCutoff, mortalityStrength, envStrength, compStrength,
                                         fission, redQueen, redQueenStrength, protracted, airmat, soilmat, nicheWidth,
                                         ddOptions);
        m_Local = NULL;
    } else if (dispersal == 2 || dispersal == 3) {
        m_Global = NULL;
        m_Local = new LocalEnvironment(X, Y, dispersal, neutral, dens, env, mort, repro, simulationEnd, specRate,
                                       dispersalCutoff, densityCutoff, mortalityStrength, envStrength, compStrength,
                                       fission, redQueen, redQueenStrength, protracted, airmat, soilmat, nicheWidth,
                                       ddOptions);
    }

    timeStep = 0;
    m_Dispersal = dispersal;
    m_X_coordinate = X;
    m_Y_coordinate = Y;
    m_SimulationEnd = simulationEnd;
    m_CheckEvery = std::max(1, 10000 / (X * Y));
    m_Start = m_NextInterrupt = m_NextPrint = std::chrono::steady_clock::now();
    m_NextPrint += std::chrono::seconds(1);
}

PhylSimModel::~PhylSimModel() {
    delete m_Global;
    delete m_Local;
}

void PhylSimModel::update(unsigned int runs) {
    for (unsigned int generation = 1; generation < runs + 1; generation++) {


#ifdef DEBUG
        if(generation % 1000 == 0){
            Rcpp::Rcout << "generation: " << generation << '/' << runs << "\n";
        }
#endif

        if (m_Dispersal == 1) {
            m_Global->increaseAge(timeStep + generation);
            m_Global->reproduce(timeStep + generation);
            m_Global->speciation(timeStep + generation);
            m_Global->m_Phylogeny.pruneIfLarge();
        } else if (m_Dispersal == 3) {
            m_Local->increaseAge(timeStep + generation);
            m_Local->reproduce(timeStep + generation);
            m_Local->speciation(timeStep + generation);
            m_Local->m_Phylogeny.pruneIfLarge();
        }
        if ((timeStep + generation) % m_CheckEvery == 0) checkProgress(timeStep + generation);
    }
    timeStep += runs;
}

// Checks for user interrupts and time limits every 0.2 s and, with m_Progress, updates the progress bar every
// second. unwindProtect() turns R's jump into a C++ exception, so that destructors free the model and RNGSync
// saves the state of the random number generator; unlike Rcpp::checkUserInterrupt(), it keeps a
// setTimeLimit() error an error instead of turning it into an interrupt.
void PhylSimModel::checkProgress(int generation) {
    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    if (now < m_NextInterrupt) return;
    m_NextInterrupt = now + std::chrono::milliseconds(200);
    try {
        Rcpp::unwindProtect([]() -> SEXP { R_CheckUserInterrupt(); return R_NilValue; });
    } catch (...) {
        if (m_LineLength > 0) Rcpp::Rcout << "\n"; // later output starts below the progress bar
        throw;
    }
    if (m_Progress && now >= m_NextPrint) {
        m_NextPrint = now + std::chrono::seconds(1);
        printProgress(generation, std::chrono::duration<double>(now - m_Start).count());
    }
}

static std::string formatSeconds(double s) {
    long long t = std::llround(s);
    char buf[32];
    if (t >= 3600) std::snprintf(buf, sizeof(buf), "%lld h %02lld min", t / 3600, t % 3600 / 60);
    else if (t >= 60) std::snprintf(buf, sizeof(buf), "%lld min %02lld s", t / 60, t % 60);
    else std::snprintf(buf, sizeof(buf), "%lld s", t);
    return buf;
}

// One console line, overwritten with \r and shorter than 80 characters, so that it does not wrap;
// the remaining time assumes a constant time per generation
void PhylSimModel::printProgress(int generation, double elapsed) {
    double done = m_SimulationEnd > 0 ? (double) generation / m_SimulationEnd : 1;
    int width = 10, bar = (int) (done * width);
    std::string line = "[" + std::string(bar, '=') + std::string(width - bar, ' ') + "] ";
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%3d%%  gen %d/%d, ", (int) (done * 100), generation, m_SimulationEnd);
    line += buf;
    if (generation < m_SimulationEnd)
        line += "elapsed " + formatSeconds(elapsed) + ", ~" +
                formatSeconds(elapsed * (m_SimulationEnd - generation) / generation) + " left";
    else line += "finished after " + formatSeconds(elapsed);
    size_t length = line.size();
    line.resize(std::max(length, m_LineLength), ' '); // overwrites the rest of a longer previous line
    m_LineLength = length;
    Rcpp::Rcout << "\r" << line << std::flush;
}

void PhylSimModel::finishProgress() {
    if (!m_Progress) return;
    printProgress(m_SimulationEnd, std::chrono::duration<double>(std::chrono::steady_clock::now() - m_Start).count());
    Rcpp::Rcout << "\n";
}
