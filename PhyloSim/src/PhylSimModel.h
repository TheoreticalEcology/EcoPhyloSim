/*
 * model.h
 *
 *  Created on: 20.06.2014
 *      Author: Paul
 */

#ifndef PHYLSIMMODEL_H_
#define PHYLSIMMODEL_H_

#include "debug.h"

#include "Grid.h"
#include <chrono>

class PhylSimModel {
public:
    int m_Dispersal;
    int m_X_coordinate;
    int m_Y_coordinate;
    int timeStep;
    // TODO FH I don't see a compelling reason to implement this as pointers, move to normal class fields?
    GlobalEnvironment *m_Global;
    LocalEnvironment *m_Local;

    PhylSimModel(int x, int y, int dispersal, int simulationEnd, double specRate, bool dens,
                 bool env, bool neutral, bool mort, int mortStrength, bool repro, double dispersalCutoff,
                 int densityCutoff, double envStrength, double compStrength,
                 int fission, double redQueen, double redQueenStrength, int protracted, std::vector<double> airmat,
                 std::vector<double> soilmat, double nicheWidth, DensityOptions ddOptions);

    ~PhylSimModel();

    void update(unsigned int runs);

    // progress bar and user interrupts during update() (#221); they draw no random numbers
    bool m_Progress = false;
    void finishProgress();

private:
    int m_SimulationEnd;
    unsigned int m_CheckEvery; // generations between clock reads, about 10^4 cell updates
    std::chrono::steady_clock::time_point m_Start, m_NextInterrupt, m_NextPrint;
    size_t m_LineLength = 0; // length of the last progress line, 0 = none printed
    void checkProgress(int generation);
    void printProgress(int generation, double elapsed);
};


#endif /* PHYLSIMMODEL_H_ */
