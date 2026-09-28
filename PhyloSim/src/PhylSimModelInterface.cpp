/*
 * Author: Betim Musa 
 * This is the main interface between the phylogeny simulation (written in C++) and R.
 */


#include <Rcpp.h>

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include "./RandomGen.h"
#include "Individual.h"
#include "Grid.h"
#include "PhylSimModel.h"
#include "Species.h"
#include "debug.h"

// basic file operations for DEBUG, remove when done with debugging
#include <iostream>
#include <fstream>

using namespace Rcpp;

//' Core phylosim model
//'
//' Low-level engine that runs the actual PhyloSim simulation. This function is called
//' internally by \code{\link{runSimulation}} and is not normally called directly by the
//' user -- most parameters correspond directly to the (user-facing) arguments of
//' \code{\link{createCompletePar}}, but several categorical/continuous parameters there
//' are split here into separate on/off switches and strength values (see Details below).
//'
//' @param x Integer, dimension of the model landscape in x-direction.
//' @param y Integer, dimension of the model landscape in y-direction.
//' @param dispersal Integer, dispersal mode/distance. 0 sets global dispersal (offspring
//'   may colonize any cell in the lattice with equal probability); integers >= 1 set local
//'   dispersal with a characteristic distance (delta) for an exponential dispersal kernel
//'   (see \code{dispersalCutoff}).
//' @param runs Integer vector, the generation(s) at which the simulation state is saved.
//'   The simulation runs until \code{max(runs)}; if \code{runs} has more than one element,
//'   intermediate states are saved at each value as well.
//' @param specRate Numeric, number of new individuals/species introduced to the community
//'   in each generation (speciation rate).
//' @param dens Logical, whether density dependence is active. If \code{FALSE}, density
//'   dependence is disabled (together with \code{env = FALSE} this yields the neutral
//'   model). The strength of density dependence is set by \code{compStrength}.
//' @param env Logical, whether environmental filtering is active. If \code{FALSE},
//'   environmental filtering is disabled (together with \code{dens = FALSE} this yields
//'   the neutral model). The strength of environmental filtering is set by \code{envStrength}.
//' @param neutral Logical, whether the neutral trait (returned as \code{NeutralTrait} /
//'   \code{neutMat}) is tracked and evolved during the simulation. This is independent of
//'   \code{dens}/\code{env} -- it does not switch the model into or out of the neutral
//'   model, it only controls whether this additional, selectively neutral marker trait is
//'   computed.
//' @param mort Logical, whether fitness differences among individuals are expressed via
//'   differential mortality. Corresponds to \code{fitnessActsOn = "mortality"} (or
//'   \code{"both"}, together with \code{repro = TRUE}) in \code{\link{createCompletePar}}.
//' @param mortStrength Integer, the fitness-based mortality ratio (must be >= 1); the
//'   higher the value, the stronger fitness differences translate into mortality
//'   differences. Corresponds to \code{fitnessBaseMortalityRatio} in
//'   \code{\link{createCompletePar}}.
//' @param repro Logical, whether fitness differences among individuals are expressed via
//'   differential reproductive output. Corresponds to \code{fitnessActsOn = "reproduction"}
//'   (or \code{"both"}, together with \code{mort = TRUE}) in \code{\link{createCompletePar}}.
//' @param dispersalCutoff Integer, cutoff distance for the local dispersal kernel (ignored
//'   if \code{dispersal = 0}); individuals only disperse within this distance, typically
//'   set to 2x the dispersal distance, to limit computation.
//' @param densityCutoff Integer, defines the effective spatial range of the competition /
//'   density-dependence calculation (ignored if \code{dens = FALSE}). Corresponds to
//'   \code{densityCut} in \code{\link{createCompletePar}}.
//' @param seed Integer, sets the random seed.
//' @param envStrength Numeric between 0 and 1, strength of the environmental filtering
//'   effect (ignored if \code{env = FALSE}); the higher the value, the stronger the
//'   influence of the environment on individual fitness. Corresponds to \code{environment}
//'   in \code{\link{createCompletePar}}.
//' @param compStrength Numeric between 0 and 1, strength of the density-dependence /
//'   competition effect (ignored if \code{dens = FALSE}); the higher the value, the
//'   stronger the influence of trait-based competition on individual fitness. Corresponds
//'   to \code{density} in \code{\link{createCompletePar}}.
//' @param fission Integer, which fission (speciation) type to use: 0 = none (default),
//'   1 = every second individual of a species becomes part of a new species, 2 = the
//'   population is geographically split into two parts, one of which becomes a new species.
//' @param redQueen Numeric, strength of the fitness decline of an aging species.
//' @param redQueenStrength Numeric, strength of the Red Queen effect; values > 0 give a
//'   newly formed species a fitness boost due to its novelty, which vanishes over time.
//' @param protracted Integer, number of generations a new species stays 'incipient' before
//'   turning into a 'good' species (protracted speciation). Default is 0 (point speciation).
//' @param airmatR Numeric vector, the environmental values of the grid, flattened
//'   column-major (length x*y), scaled between 0 and 1. Corresponds to \code{airmat} in
//'   \code{\link{createCompletePar}}.
//' @param soilmatR Numeric vector, a second, optional environmental grid, flattened
//'   column-major (length x*y), analogous to \code{airmatR}. TODO: confirm exact role --
//'   not currently exposed as a separate argument in \code{createCompletePar}.
//' @param prunePhylogeny Logical, whether the returned phylogeny is pruned to extant
//'   (non-fossil) tips by the internal pruning function, or the full phylogeny (including
//'   extinct lineages) is returned.
//'
//' @return A list with one element per entry in \code{runs}, each containing the raw
//'   simulation output for that generation: \code{Species} (species identity per grid
//'   cell), \code{EnvTrait} (environmental trait per individual), \code{NeutralTrait}
//'   (neutral trait per individual), \code{CompetitionTrait} (competition trait per
//'   individual), \code{Environment} (environmental value per grid cell) and
//'   \code{Phylogeny} (Newick string). \code{\link{runSimulation}} reshapes and renames
//'   these into the final \code{Phylosim} object (\code{specMat}, \code{traitMat},
//'   \code{neutMat}, \code{compMat}, \code{envMat}, \code{phylogeny}, \code{phyloTXT}).
//' @export
// [[Rcpp::export]]
List callModel(int x, int y, int dispersal, IntegerVector runs, double specRate, bool dens, bool env, bool neutral,
               bool mort, int mortStrength, bool repro, int dispersalCutoff, int densityCutoff, int seed,
               double envStrength, double compStrength, int fission, double redQueen, double redQueenStrength,
               int protracted, NumericVector airmatR, NumericVector soilmatR, bool prunePhylogeny) {
#ifdef DEBUG
    std::ofstream debugFile;
    debugFile.open("debug.txt");
    std::cout << "Writing parameters to debug.txt..." << ";" << std::endl;
    debugFile << "int x = " << x << ";" << std::endl;
    debugFile << "int y = " << y << ";" << std::endl;
    debugFile << "int dispersal = " << dispersal << ";" << std::endl;
    debugFile << "IntegerVector runs = " << runs << ";" << std::endl;
    debugFile << "double specRate = " << specRate << ";" << std::endl;
    debugFile << "bool dens = " << dens << ";" << std::endl;
    debugFile << "bool env = " << env << ";" << std::endl;
    debugFile << "bool neutral = " << neutral << ";" << std::endl;
    debugFile << "bool mort = " << mort << ";" << std::endl;
    debugFile << "int mortStrength = " << mortStrength << ";" << std::endl;
    debugFile << "bool repro = " << repro << ";" << std::endl;
    debugFile << "int dispersalCutoff = " << dispersalCutoff << ";" << std::endl;
    debugFile << "int densityCutoff = " << densityCutoff << ";" << std::endl;
    debugFile << "int seed = " << seed << ";" << std::endl;
    debugFile << "double envStrength = " << envStrength << ";" << std::endl;
    debugFile << "double compStrength = " << compStrength << ";" << std::endl;
    debugFile << "int fission = " << fission << ";" << std::endl;
    debugFile << "double redQueen = " << redQueen << ";" << std::endl;
    debugFile << "double redQueenStrength = " << redQueenStrength << ";" << std::endl;
    debugFile << "int protracted = " << protracted << ";" << std::endl;
    debugFile << "NumericVector airmatR = " << airmatR << ";" << std::endl;
    debugFile << "airmatR = " << airmatR.begin() <<", " << airmatR.end() << ";" << std::endl;
    debugFile << "NumericVector airmatR = " << soilmatR << ";" << std::endl;
    debugFile << "soilmatR = " << soilmatR.begin() <<", " << soilmatR.end() << ";" << std::endl;
    debugFile.close();
#endif

    RandomGen ran;
    ran.seedrand(seed); // seed is int while seedrand expects unsigned int

    std::string tempSaveLoc = "./";

    int nRuns = runs.length(); // number of consecutive simulations that are executed sequentially
    int prevNGen = 0;

    std::vector<double> airmat(airmatR.begin(), airmatR.end());
    std::vector<double> soilmat(soilmatR.begin(), soilmatR.end());

    Rcpp::List outList = Rcpp::List::create();

    PhylSimModel phylSimModel(x, y, dispersal, runs[nRuns - 1], specRate, dens, env, neutral, mort, mortStrength, repro,
                              dispersalCutoff, densityCutoff, tempSaveLoc, envStrength, compStrength, fission, redQueen,
                              redQueenStrength, protracted, airmat, soilmat);


    for (int step = 0; step < nRuns; step++) {

        Rcpp::IntegerVector specOut(x * y);
        Rcpp::NumericVector traitOut(x * y);
        Rcpp::NumericVector envOut(x * y);
        Rcpp::NumericVector neutralOut(x * y);
        Rcpp::NumericVector compOut(x * y);
        Rcpp::CharacterVector phyloOut = "";


        int curNGen = runs[step] - prevNGen;

#ifdef DEBUG
        std::cout << "Run model for" <<  curNGen << " generations\n";
#endif
        phylSimModel.update(curNGen);

        prevNGen = runs[step];

        int indCounter = 0;

        if (dispersal == 1) {
            for (int i = 0; i < x; i++) {
                for (int j = 0; j < y; j++) {
                    specOut[indCounter] = phylSimModel.m_Global->m_Individuals[i][j].m_Species->m_ID;
                    traitOut[indCounter] = phylSimModel.m_Global->m_Individuals[i][j].m_Mean;
                    envOut[indCounter] = phylSimModel.m_Global->m_Environment[i * y + j].first;
                    neutralOut[indCounter] = phylSimModel.m_Global->m_Individuals[i][j].m_NeutralMarker;
                    compOut[indCounter] = phylSimModel.m_Global->m_Individuals[i][j].m_CompetitionMarker;
                    indCounter++;
                }
            }
            phylSimModel.m_Global->m_Phylogeny.prunePhylogeny(prevNGen);
            std::string phyloPass("\0");

#ifdef DEBUG
            std::cout << "before writePhylogenyR" << std::endl;
            std::cout << prunePhylogeny << std::endl;
            // consistency checks -> length, etc
#endif
            if (prunePhylogeny) {
                phyloPass = phylSimModel.m_Global->m_Phylogeny.writePhylogenyR(1,
                                                                               phylSimModel.m_Global->m_Phylogeny.m_PrunedPhylo);
            } else {
                phyloPass = phylSimModel.m_Global->m_Phylogeny.writePhylogenyR(1,
                                                                               phylSimModel.m_Global->m_Phylogeny.m_FullPhylogeny);
            }
#ifdef DEBUG
            std::cout << "after writePhylogenyR" << std::endl;
#endif
            char *cstr = new char[phyloPass.length() + 1];
            std::strcpy(cstr, phyloPass.c_str());
            phyloOut[0] = cstr;
            delete[] cstr;

        } else {
            for (int i = 0; i < x; i++) {
                for (int j = 0; j < y; j++) {
                    specOut[indCounter] = phylSimModel.m_Local->m_Individuals[i][j].m_Species->m_ID;
                    traitOut[indCounter] = phylSimModel.m_Local->m_Individuals[i][j].m_Mean;
                    envOut[indCounter] = phylSimModel.m_Local->m_Environment[i * y + j].first;
                    neutralOut[indCounter] = phylSimModel.m_Local->m_Individuals[i][j].m_NeutralMarker;
                    compOut[indCounter] = phylSimModel.m_Local->m_Individuals[i][j].m_CompetitionMarker;
                    indCounter++;
                }
            }
            phylSimModel.m_Local->m_Phylogeny.prunePhylogeny(prevNGen);

            std::string phyloPass("\0");

            if (prunePhylogeny) {
                phyloPass = phylSimModel.m_Local->m_Phylogeny.writePhylogenyR(1,
                                                                              phylSimModel.m_Local->m_Phylogeny.m_PrunedPhylo);
            } else {
                phyloPass = phylSimModel.m_Local->m_Phylogeny.writePhylogenyR(1,
                                                                              phylSimModel.m_Local->m_Phylogeny.m_FullPhylogeny);
            }

            char *cstr = new char[phyloPass.length() + 1];
            std::strcpy(cstr, phyloPass.c_str());
            phyloOut[0] = cstr;
            delete[] cstr;
        }

        Rcpp::List listResults = Rcpp::List::create(Rcpp::Named("Species") = specOut,
                                                    Rcpp::Named("EnvTrait") = traitOut,
                                                    Rcpp::Named("NeutralTrait") = neutralOut,
                                                    Rcpp::Named("CompetitionTrait") = compOut,
                                                    Rcpp::Named("Environment") = envOut,
                                                    Rcpp::Named("Phylogeny") = phyloOut);

        outList.push_back(listResults);
    }


    return outList;
}
