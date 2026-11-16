/*
 * model.c++
 *
 *  Created on: 20.06.2014
 *      Author: Paul
 *              Betim Musa <musab@informatik.uni-freiburg.de>
 */

#include <Rcpp.h>
#include <string>
#include "./RandomGen.h"
//#include "Utils/RandomGen.h"
#include "Individual.h"
#include "Grid.h"
#include "PhylSimModel.h"
#include "Species.h"


PhylSimModel::PhylSimModel(int X, int Y, int dispersal, int simulationEnd, double specRate, bool dens,
                           bool env, bool neutral, bool mort, int mortalityStrength, bool repro, int dispersalCutoff,
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
}

PhylSimModel::~PhylSimModel() {
    delete m_Global;
    delete m_Local;
}

void PhylSimModel::update(unsigned int runs) {
    for (unsigned int generation = 1; generation < runs + 1; generation++) {

        // std::cout << "generation :" <<  generation << '/' << runs << '\n';

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
            //  std::cout << "Disp generation :" <<  generation << '/' << runs << '\n';
        } else if (m_Dispersal == 3) {
            m_Local->increaseAge(timeStep + generation);
            m_Local->reproduce(timeStep + generation);
            m_Local->speciation(timeStep + generation);
            m_Local->m_Phylogeny.pruneIfLarge();
            //   std::cout << "Disp generation :" <<  generation << '/' << runs << '\n';
        }
    }
    timeStep += runs;
}
