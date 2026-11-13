/*
 * phylogeny.cpp
 *
 *  Created on: 30.07.2014
 *      Author: Paul
 */

#include "Phylogeny.h"
#include "StringConversion.h"

#include <Rcpp.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include "debug.h"
#include <sstream>
#include <stdexcept>
#include <utility>

Phylogeny::Phylogeny() {
    m_FullPhylogeny = new std::multimap<unsigned long long, Species*>();
    m_PruneDuringRun = false;
    m_PruneThreshold = 1000;
}

Phylogeny::~Phylogeny() {
    for (auto &entry : *m_FullPhylogeny) delete entry.second;
    delete m_FullPhylogeny;
}

void Phylogeny::updatePhylogeny(Species *spec) {
    m_FullPhylogeny->insert(std::make_pair(spec->m_ID, spec));
}


void Phylogeny::writeSpeciesData() {
    std::ofstream specData("..\\species_data.txt");
    specData << "ID, nChildren, environment, competition, neutral, emergence, extinction, ancestor, count" << '\n';
    Species *spec;
    for (unsigned int i = 1; i < m_FullPhylogeny->size(); i++) {
        spec = m_FullPhylogeny->find(i)->second;
        specData << spec->m_ID << ',' << spec->m_Children.size() << ',' << spec->m_Mean << ',';
        specData << spec->m_CompetitionMean << ',' << spec->m_NeutralMean << ',' << spec->m_Date_of_Emergence << ','
                 << spec->m_Date_of_Extinction << ',' << spec->m_Ancestor << ',' << spec->m_Count << '\n';
    }

}


// Removes and deletes all species without living individuals and without remaining
// descendants, i.e. lineages that went extinct completely. Must be called between
// generations, when no individual points to a species with m_Count == 0 (#140).
void Phylogeny::prunePhylogeny() {
    // children have higher IDs than their parents, so one pass in descending ID order suffices
    for (auto it = m_FullPhylogeny->end(); it != m_FullPhylogeny->begin();) {
        --it;
        std::vector<unsigned long long> &children = it->second->m_Children;
        children.erase(std::remove_if(children.begin(), children.end(), [this](unsigned long long id) {
            return m_FullPhylogeny->find(id) == m_FullPhylogeny->end();
        }), children.end());
        if (it->second->m_Count == 0 && children.empty() && it->first != 1) {
            delete it->second;
            it = m_FullPhylogeny->erase(it);
        }
    }
}

// Prunes during the run each time the phylogeny has doubled in size, so that memory
// does not grow with the number of species created between censuses
void Phylogeny::pruneIfLarge() {
    if (!m_PruneDuringRun || m_FullPhylogeny->size() < m_PruneThreshold) return;
    prunePhylogeny();
    m_PruneThreshold = std::max<size_t>(1000, 2 * m_FullPhylogeny->size());
}

void Phylogeny::writePhylogeny(unsigned long long start, std::multimap<unsigned long long, Species *> *phylogenyMap,
                               char suffix) {
//	 std::cout << "write Phylogeny" << std::endl;
    Species *focusSpecies = phylogenyMap->find(start)->second;
    Species *ancestorSpecies = new Species(0, 0, 0, std::make_pair(0, 0), 0);
    phylogenyMap->insert(std::make_pair(ancestorSpecies->m_ID, ancestorSpecies));
    focusSpecies->m_Ancestor = 0;
    //std::map<unsigned long long, int> position; // Initialisieren ??
//	std::vector <unsigned long long> phylchildren;

//	unsigned long long child;
//	unsigned int count = 0;
    int branchLength = 0;
//	int rootLength = 0;
//	double parentBranchlength = 0.0;
//	unsigned long long ancestor = 0;

//	phylchildren = parent->children;
    //position[parent->ID] = 0 ;




    std::stringstream filename;
    filename << "..\\phylo_newick" << suffix << ".txt";
    std::ofstream phylo(filename.str().c_str());

    std::string tree = ";";


    std::vector <size_t> position(phylogenyMap->size(), 0);
    std::vector <size_t> paranthesis(phylogenyMap->size(), 0);


    bool go = true;
    while (go) {

//		if(position[parent->ID]+1 <= parent->children.size() )
//		{

        // todo long term: think about changing children from integer vector + map to pointer


        // If parent is leaf, write down leaf and length and go one up
        if (focusSpecies->m_Children.empty()) {
            if (position[focusSpecies->m_ID] != 0)
                throw std::runtime_error("unexpected value for position value in a leaf");
            branchLength = focusSpecies->m_Date_of_Extinction - focusSpecies->m_Date_of_Emergence;
            tree.insert(0, to_string(branchLength));
            tree.insert(0, ":");
            tree.insert(0, to_string(focusSpecies->m_ID));
            tree.insert(0, "s");
            tree.insert(0, ",");
            focusSpecies = phylogenyMap->find(focusSpecies->m_Ancestor)->second;
            position[focusSpecies->m_ID] += 1;
        }

            // if all children visited write down parent, close and go one up
        else if (position[focusSpecies->m_ID] == focusSpecies->m_Children.size()) {
            branchLength = focusSpecies->m_Date_of_Extinction -
                           phylogenyMap->find(focusSpecies->m_Children.back())->second->m_Date_of_Emergence;

            tree.insert(0, to_string(branchLength));
            tree.insert(0, ":");
            tree.insert(0, to_string(focusSpecies->m_ID));
            tree.insert(0, "s");
            tree.insert(0, paranthesis[focusSpecies->m_ID], '(');
            tree.insert(0, ",");
            focusSpecies = phylogenyMap->find(focusSpecies->m_Ancestor)->second;

            position[focusSpecies->m_ID] += 1;

        }

            // If Parent has children and not yet gone to each child
            // then go to child
        else if (position[focusSpecies->m_ID] < focusSpecies->m_Children.size()) {

            // if two siblings are born in the same generation
            if (position[focusSpecies->m_ID] > 1 &&
                phylogenyMap->find(
                        focusSpecies->m_Children[position[focusSpecies->m_ID]])->second->m_Date_of_Emergence ==
                phylogenyMap->find(
                        focusSpecies->m_Children[position[focusSpecies->m_ID] - 1])->second->m_Date_of_Emergence &&
                phylogenyMap->find(
                        focusSpecies->m_Children[position[focusSpecies->m_ID]])->second->m_Date_of_Emergence ==
                phylogenyMap->find(
                        focusSpecies->m_Children[position[focusSpecies->m_ID] + 1])->second->m_Date_of_Emergence) {
                if (!(phylogenyMap->find(
                        focusSpecies->m_Children[position[focusSpecies->m_ID] - 1])->second->m_Children.empty()) &&
                    (phylogenyMap->find(
                            focusSpecies->m_Children[position[focusSpecies->m_ID]])->second->m_Children.empty())) {

                    ;
//						tree.insert(0, "," );
                }
            } else {
                // if first time going through this node branch lenght is distance to parent node
                if (position[focusSpecies->m_ID] == 0) {
                    // TODO: unsigned long long assigned to int!
                    branchLength = phylogenyMap->find(
                            focusSpecies->m_Children[position[focusSpecies->m_ID]])->second->m_Date_of_Emergence -
                                   focusSpecies->m_Date_of_Emergence;
                }
                    // else branch lenght is distance to previous sibling
                else // ????? or position = children size
                {
                    // TODO: unsigned long long assigned to int!
                    branchLength = phylogenyMap->find(
                            focusSpecies->m_Children[position[focusSpecies->m_ID]])->second->m_Date_of_Emergence -
                                   phylogenyMap->find(
                                           focusSpecies->m_Children[position[focusSpecies->m_ID] -
                                                                    1])->second->m_Date_of_Emergence;
                }

                tree.insert(0, to_string(branchLength));
                tree.insert(0, ":");
                tree.insert(0, to_string(focusSpecies->m_ID));
                tree.insert(0, "s");
                tree.insert(0, ")");
                paranthesis[focusSpecies->m_ID] += 1;
            }

            focusSpecies = phylogenyMap->find(focusSpecies->m_Children[position[focusSpecies->m_ID]])->second;
        } else throw std::runtime_error("Nope");


        if (focusSpecies->m_ID == 0) {
            go = false;
        }
    }
    tree.erase(tree.begin());
 
    phylo << tree;
}


std::string
Phylogeny::writePhylogenyR(unsigned long long start, std::multimap<unsigned long long, Species *> *phylogenyMap,
                           unsigned long long current) {

#ifdef DEBUG
    Rcpp::Rcout << "writePhylogenyR start\n";
#endif

    Species *focusSpecies = phylogenyMap->find(start)->second;
    Species *ancestorSpecies = new Species(0, 0, 0, std::make_pair(0, 0), 0);
    phylogenyMap->insert(std::make_pair(ancestorSpecies->m_ID, ancestorSpecies));
    focusSpecies->m_Ancestor = 0;

    int branchLength = 0;

    std::string tree = ";";
    // indexed by species ID, which has gaps after pruning (#140)
    std::vector <size_t> position(phylogenyMap->rbegin()->first + 1, 0);
    std::vector <size_t> paranthesis(phylogenyMap->rbegin()->first + 1, 0);
    int counter = 1;

    bool go = true;
    while (go) {

        // todo long term: think about changing children from integer vector + map to pointer


        // If parent is leaf, write down leaf and length and go one up
        if (focusSpecies->m_Children.empty()) {

            if (position[focusSpecies->m_ID] != 0)
                throw std::runtime_error("unexpected value for position value in a leaf");
            // TODO: unsigned long long assigned to int!
            branchLength = std::min<unsigned long long>(focusSpecies->m_Date_of_Extinction, current) -
                           focusSpecies->m_Date_of_Emergence;
            tree.insert(0, to_string(branchLength));
            tree.insert(0, ":");
            tree.insert(0, to_string(focusSpecies->m_ID));
            tree.insert(0, "s");
            tree.insert(0, ",");
            focusSpecies = phylogenyMap->find(focusSpecies->m_Ancestor)->second;
            position[focusSpecies->m_ID] += 1;

        }

            // if all children visited write down parent, close and go one up
        else if (position[focusSpecies->m_ID] == focusSpecies->m_Children.size()) {

            // TODO: unsigned long long assigned to int!
            branchLength = std::min<unsigned long long>(focusSpecies->m_Date_of_Extinction, current) -
                           phylogenyMap->find(focusSpecies->m_Children.back())->second->m_Date_of_Emergence;

            tree.insert(0, to_string(branchLength));
            tree.insert(0, ":");
            tree.insert(0, to_string(focusSpecies->m_ID));
            tree.insert(0, "s");
            tree.insert(0, paranthesis[focusSpecies->m_ID], '(');
            tree.insert(0, ",");
            focusSpecies = phylogenyMap->find(focusSpecies->m_Ancestor)->second;
            position[focusSpecies->m_ID] += 1;

        }

            // If Parent has children and not yet gone to each child then go to child
        else if (position[focusSpecies->m_ID] < focusSpecies->m_Children.size()) {

            // if two siblings are born in the same generation
            //if(position[parent->m_ID]==3) throw std::runtime_error("Cannot build phylogeny");
            // TODO the bug arises when evaluating the if statement
            // TODO check
            // Here the program can crash in case of high speciation rates.



            // If already one child visited, but why is this
            if (position[focusSpecies->m_ID] > 1 &&
                phylogenyMap->find(
                        focusSpecies->m_Children[position[focusSpecies->m_ID]])->second->m_Date_of_Emergence ==
                phylogenyMap->find(
                        focusSpecies->m_Children[position[focusSpecies->m_ID] - 1])->second->m_Date_of_Emergence
                // TODO find out why this causes a crash and sort it out
                /*&&
                phylogenyMap->find(focusSpecies->m_Children[position[focusSpecies->m_ID]])->second->m_Date_of_Emergence ==
                phylogenyMap->find(focusSpecies->m_Children[position[focusSpecies->m_ID] + 1])->second->m_Date_of_Emergence*/
                    ) {

                if (!(phylogenyMap->find(
                        focusSpecies->m_Children[position[focusSpecies->m_ID] - 1])->second->m_Children.empty()) &&
                    (phylogenyMap->find(
                            focusSpecies->m_Children[position[focusSpecies->m_ID]])->second->m_Children.empty())) {

                    ;
                    //tree.insert(0, "," );
                }
            } else {
                // if first time going through this node branch lenght is distance to parent node
                if (position[focusSpecies->m_ID] == 0) {
                    // TODO: unsigned long long assigned to int!
                    branchLength = phylogenyMap->find(
                            focusSpecies->m_Children[position[focusSpecies->m_ID]])->second->m_Date_of_Emergence -
                                   focusSpecies->m_Date_of_Emergence;

                }
                    // else branch lenght is distance to previous sibling
                else // ????? or position = children size
                {

                    // TODO: unsigned long long assigned to int!
                    branchLength = phylogenyMap->find(
                            focusSpecies->m_Children[position[focusSpecies->m_ID]])->second->m_Date_of_Emergence -
                                   phylogenyMap->find(
                                           focusSpecies->m_Children[position[focusSpecies->m_ID] -
                                                                    1])->second->m_Date_of_Emergence;

                }
                tree.insert(0, to_string(branchLength));
                tree.insert(0, ":");
                tree.insert(0, to_string(focusSpecies->m_ID));
                tree.insert(0, "s");
                tree.insert(0, ")");
                paranthesis[focusSpecies->m_ID] += 1;
            }

            focusSpecies = phylogenyMap->find(focusSpecies->m_Children[position[focusSpecies->m_ID]])->second;

        } else break;
        //else throw std::runtime_error("Nope");

        counter += 1;
        //std::cout << "Counter: " << counter << " Size: " << phylogenyMap->size() << std::endl;
        //std::cout << "Position max: " << *std::max_element(position.begin(), position.end()) << std::endl;

        // a complete walk takes 2 * N - 1 steps for N species (#105)
        if (counter > 2 * phylogenyMap->size()) Rcpp::stop("writePhylogenyR: tree traversal did not finish");

        if (focusSpecies->m_ID == 0) {
            go = false;
        }
    }
    //std::cout << tree << std::endl;

    phylogenyMap->erase(0); // fixes #105 (remove temporary root again)
    delete ancestorSpecies;
    tree.erase(tree.begin());
#ifdef DEBUG
    Rcpp::Rcout << "writePhylogenyR end\n";
#endif
    return tree;
}
