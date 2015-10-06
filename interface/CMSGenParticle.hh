#ifndef RIVET_CMSGenParticle_HH
#define RIVET_CMSGenParticle_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Event.hh"
//#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  // @brief Pseudo top finder
  // 
  // Find top quark in the particle level.
  // The definition is based on the agreement at the LHC working group.
  class CMSGenParticle : public FinalState {
  public:
    /// @name Standard constructors and destructors.
    //@{

    /// The default constructor.
    CMSGenParticle()
      : FinalState(-MAXDOUBLE, MAXDOUBLE, 0*GeV),
        _vetoIds({1000022,1000012, 1000014, 1000016,
                  2000012, 2000014, 2000016, 1000039, 5100039,
                  4000012, 4000014, 4000016, 9900012, 9900014, 9900016, 39}),
        _vetoIdsFromResonances({12, 13, 14, 16})
    {
      setName("CMSGenParticle");
    }

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new CMSGenParticle(*this);
    }

    //@}

// partonicFinalState = false
// exclude Resonances  = true
// tauAsjets = false
  protected:
    // Apply the projection to the event
    void project(const Event& e) override;

    bool isParton(int pdgId) const {
      const int iid = abs(pdgId) % 10000;
      return (iid > 0 && iid < 6) || iid == 7 || iid == 9 || iid == 21;
    }
    bool isResonance(const int pdgId) const {
      const int iid = abs(pdgId) % 10000;
      return (iid > 21 && iid <= 42) || iid == 6 || iid == 8;
    }
    bool isIgnored(const unsigned int absId) const {
      auto pos = std::lower_bound(_vetoIds.begin(), _vetoIds.end(), absId);
      return pos != _vetoIds.end() && *pos == absId;
    }
    bool isExcludedFromResonance(const unsigned int absId) const {
      auto pos = std::lower_bound(_vetoIdsFromResonances.begin(), _vetoIdsFromResonances.end(), absId);
      return pos != _vetoIdsFromResonances.end() && *pos == absId;
    }

    int fromResonance(std::set<GenParticle*>& invalid, const std::vector<GenParticle*>& pv, GenParticle* p) const;

  protected:
    std::vector<unsigned int> _vetoIds;
    std::vector<unsigned int> _vetoIdsFromResonances;

  };

}

#endif
