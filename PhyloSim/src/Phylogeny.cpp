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

// Writes the phylogeny in Newick format. Each species is a lineage from its emergence to
// its extinction (or the present); each speciation event is a node labeled with the parent
// species, and children that emerge in the same generation share one node (polytomy).
// Events of generation g are placed at time g, the census at time current + 0.5, so that tips
// of extant species end 0.5 generations after those of species that went extinct in the
// last generation (#147). With m_PruneDuringRun, prunePhylogeny() must have been called
// before, and only extant species are written as tips: extinct ancestors become internal
// nodes, and nodes with a single descendant are removed. The tree has no root edge.
std::string
Phylogeny::writePhylogenyR(unsigned long long start, std::multimap<unsigned long long, Species *> *phylogenyMap,
                           unsigned long long current) {

    const bool extantOnly = m_PruneDuringRun;
    auto species = [phylogenyMap](unsigned long long id) { return phylogenyMap->find(id)->second; };
    auto hasTip = [extantOnly](const Species *s) { return !extantOnly || s->m_Count > 0; };
    // index after the children of s that emerge in the same generation as child i
    auto groupEnd = [&species](const Species *s, size_t i) {
        unsigned long long date = species(s->m_Children[i])->m_Date_of_Emergence;
        while (i < s->m_Children.size() && species(s->m_Children[i])->m_Date_of_Emergence == date) i++;
        return i;
    };
    auto restNonEmpty = [&hasTip](const Species *s, size_t i) { return i < s->m_Children.size() || hasTip(s); };

    // A segment (s, i) is the lineage of s from the emergence of its child i onwards;
    // (s, 0) is the subtree of s. Frames are kept on a stack, as the tree can be deep.
    struct Frame {
        Species *s;
        size_t begin, end, next;  // children [begin, end) of s form this node, then (s, end)
        unsigned long long time, parentTime;
    };
    std::vector<Frame> stack;
    std::string tree;

    // writes a tip, or opens the node of segment (s, i) after skipping single-child nodes
    auto openSegment = [&](Species *s, size_t i, unsigned long long parentTime, bool isRoot) {
        while (i < s->m_Children.size() && groupEnd(s, i) == i + 1 && !restNonEmpty(s, i + 1)) {
            s = species(s->m_Children[i]);
            i = 0;
        }
        if (i == s->m_Children.size()) {
            if (!hasTip(s)) Rcpp::stop("writePhylogenyR: extinct species without descendants");
            tree += "s" + to_string(s->m_ID);
            if (isRoot) return;
            bool extant = s->m_Count > 0;
            unsigned long long end = extant ? current : s->m_Date_of_Extinction;
            tree += ":" + to_string(end - parentTime) + (extant ? ".5" : "");
        } else {
            tree += "(";
            stack.push_back({s, i, groupEnd(s, i), i,
                             species(s->m_Children[i])->m_Date_of_Emergence, parentTime});
        }
    };

    openSegment(species(start), 0, 0, true);
    while (!stack.empty()) {
        Frame &f = stack.back();
        if (f.next <= f.end && (f.next < f.end || restNonEmpty(f.s, f.end))) {
            if (f.next > f.begin) tree += ",";
            size_t k = f.next++;
            if (k < f.end) openSegment(species(f.s->m_Children[k]), 0, f.time, false);
            else openSegment(f.s, f.end, f.time, false);
        } else {
            tree += ")s" + to_string(f.s->m_ID);
            unsigned long long time = f.time, parentTime = f.parentTime;
            stack.pop_back();
            if (!stack.empty()) tree += ":" + to_string(time - parentTime);
        }
    }
    return tree + ";";
}
