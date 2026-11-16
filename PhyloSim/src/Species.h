/*
 * species.h
 *
 *  Created on: 26.06.2014
 *      Author: Paul
 *              Betim Musa <musab@informatik.uni-freiburg.de>
 */

#ifndef SPECIES_H_
#define SPECIES_H_

#define _USE_MATH_DEFINES

#include <utility>
#include  <vector>


class Species {
private:
    double m_MeanSum;
    double m_NeutralSum;
    double m_CompetitionSum;
    double m_Variance;

    void updateMean();

public:
    unsigned long long m_Count;
    unsigned long long m_Ancestor;
    std::pair<int, int> m_Origin;

    unsigned long long m_ID;
    unsigned long long m_Date_of_Emergence;
    unsigned int m_Date_of_Extinction;
    double m_Mean;
    double m_FirstMean, m_FirstNeutral, m_FirstComp;
    double m_NeutralMean;
    double m_CompetitionMean;
    std::vector<unsigned long long> m_Children;


    Species();

    Species(const Species &spec);

    void operator=(const Species &spec);


    Species(unsigned long long id, unsigned long long father,
            unsigned long long date, std::pair<int, int> birthplace, unsigned int runs);

    ~Species();

    void addIndividual(double env, double comp, double neutral);

    void removeIndividual(double env, double comp, double neutral, int generation);

    unsigned long long get_species_ID();
};


#endif /* SPECIES_H_ */
