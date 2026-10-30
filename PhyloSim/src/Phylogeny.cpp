/*
 * phylogeny.cpp
 *
 *  Created on: 30.07.2014
 *      Author: Paul
 */

#include "Phylogeny.h"
#include "StringConversion.h"

#include <fstream>
#include <iostream>
#include <iterator>
#include "debug.h"
#include <sstream>
#include <stdexcept>
#include <utility>

Phylogeny::Phylogeny() {
    m_PrunedPhylo = new std::multimap<unsigned long long, Species*>();
    m_FullPhylogeny = new std::multimap<unsigned long long, Species*>();
}

Phylogeny::~Phylogeny() {

    for (unsigned int i = 1; i < m_FullPhylogeny->size(); i++) {
        delete m_FullPhylogeny->find(i)->second;
    }
    m_FullPhylogeny->clear();
    delete m_FullPhylogeny;

    for (unsigned int i = 1; i < m_PrunedPhylo->size(); i++) {
        delete m_PrunedPhylo->find(i)->second;
    }
    m_PrunedPhylo->clear();
    delete m_PrunedPhylo;
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


void Phylogeny::prunePhylogeny(int current) {
//  std::cout << "Prune Phylogeny" << std::endl;
    m_PrunedPhylo->clear();

    // Create deep copy
    // TODO -> I tried the deep copy of multimap, but it doesn't work ... sort this out!!!
    for (unsigned long long i = 1; i <= m_FullPhylogeny->size(); i++) {
        Species *specp = new Species(*m_FullPhylogeny->find(i)->second);
        if (specp->m_Date_of_Extinction > current) specp->m_Date_of_Extinction = current;
        m_PrunedPhylo->insert(std::pair <unsigned long long, Species*> (specp->m_ID, specp));
    }


    /* TODO - I somehow don't understand 100% why this works .. we go through the phylogeny size, erase things inbetween,
     * why don't we get into trouble for trying to access indices that have already been erase? */

    for (unsigned long long i = 1; i <= m_PrunedPhylo->size(); i++) {
        Species *father = m_PrunedPhylo->find(i)->second;

        for (int i = 0; i < 10; i++) {
            if (!(father->m_Children.empty())) {
                for (unsigned int j = 0; j < father->m_Children.size(); j++) {
                    Species *child = m_PrunedPhylo->find(father->m_Children[j])->second;
                    if (child->m_Children.empty() && child->m_Count == 0) {
                        father->m_Children.erase(father->m_Children.begin() + j);
                    }
                    //				 else if(!(child->children.empty()) && child->count==0)
                    //				 {
                    //					 father->children.erase(father->children.begin()+j);
                    //					 father->children.insert(father->children.end(), child->children.begin(), child->children.end());
                    //				 }
                }
            }
        }
        // std::cout << father->children.size() << '\n';
    }
    // std::cout <<"phylosize : " << prunedPhylo.size() << '\n';
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
Phylogeny::writePhylogenyR(unsigned long long start, std::multimap<unsigned long long, Species *> *phylogenyMap) {

#ifdef DEBUG
   
#endif

    Species *focusSpecies = phylogenyMap->find(start)->second;
    Species *ancestorSpecies = new Species(0, 0, 0, std::make_pair(0, 0), 0);
    phylogenyMap->insert(std::make_pair(ancestorSpecies->m_ID, ancestorSpecies));
    focusSpecies->m_Ancestor = 0;

    int branchLength = 0;

    std::string tree = ";";
    std::vector <size_t> position(phylogenyMap->size(), 0);
    std::vector <size_t> paranthesis(phylogenyMap->size(), 0);
    int counter = 1;

    bool go = true;
    while (go) {

        // todo long term: think about changing children from integer vector + map to pointer


        // If parent is leaf, write down leaf and length and go one up
        if (focusSpecies->m_Children.empty()) {

            if (position[focusSpecies->m_ID] != 0)
                throw std::runtime_error("unexpected value for position value in a leaf");
            // TODO: unsigned long long assigned to int!
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

            // TODO: unsigned long long assigned to int!
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

        if (counter > phylogenyMap->size()) break;

        if (focusSpecies->m_ID == 0) {
            go = false;
        }
    }
    //std::cout << tree << std::endl;

    tree.erase(tree.begin());
#ifdef DEBUG
   
#endif
    return tree;
}
