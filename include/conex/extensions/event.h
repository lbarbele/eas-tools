#ifndef _conex_extensions_event_h
#define _conex_extensions_event_h

#include <vector>

#include <conex/extensions/interaction.h>

#include <TTree.h>

namespace conex::extensions {

  class event {
  private:
    std::vector<interaction> m_interactions;

  public:
    // * build an event from the respective trees in the conex::extensions file
    event(TTree* particle, TTree* projectile, TTree* interaction, TTree* seed, const double threshold, const bool check = true);

    // * get the number of interactions in this event
    size_t get_n_interactions() const
    {return m_interactions.size();}

    // * get interaction at index pos
    const interaction& get_interaction(size_t pos) const
    {return m_interactions.at(pos);}

    // - Interaction iterators

    // * standard iterators
    auto begin() {return m_interactions.begin();}
    auto end() {return m_interactions.end();}
    
    // * const iterators
    auto begin() const {return m_interactions.begin();}
    auto cbegin() const {return m_interactions.cbegin();}
    auto end() const {return m_interactions.end();}
    auto cend() const {return m_interactions.cend();}

    // * reverse iterators
    auto rbegin() {return m_interactions.rbegin();}
    auto rend() {return m_interactions.rend();}

    // * const reverse iterators
    auto rbegin() const {return m_interactions.rbegin();}
    auto crbegin() const {return m_interactions.crbegin();}
    auto rend() const {return m_interactions.rend();}
    auto crend() const {return m_interactions.crend();}
  };

}

#endif