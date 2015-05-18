#ifndef RIVET_PartonTTbarState_HH
#define RIVET_PartonTopFinter_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Event.hh"

namespace Rivet {

  // @brief Parton level top quark finder
  // 
  // Find top quark in the parton level directly tracking particle history.
  // This does not fit with the Rivet philosophy and can be generator dependent,
  // so please use this with your own risks.
  class PartonTTbarState : public FinalState {
  public:
    enum TTbarMode { CH_FULLHADRON = 0, CH_SEMILEPTON, CH_FULLLEPTON };
    enum DecayMode { CH_HADRON = 0, CH_MUON, CH_ELECTRON, CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON };

    /// @name Standard constructors and destructors.
    //@{

    /// The default constructor. May specify the minimum and maximum
    /// pseudorapidity \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
    PartonTTbarState(double mineta = -MAXDOUBLE,
                    double maxeta =  MAXDOUBLE,
                    double minpt = 0.0*GeV)
      : FinalState(mineta, maxeta, minpt)
    {
      setName("PartonTop");
    }

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new PartonTTbarState(*this);
    }

    //@}

  public:
    TTbarMode mode() const { 
      if ( _mode1 == CH_HADRON && _mode2 == CH_HADRON ) return CH_FULLHADRON;
      else if ( _mode1 != CH_HADRON && _mode2 != CH_HADRON ) return CH_FULLLEPTON;
      else return CH_SEMILEPTON;
    }
    DecayMode mode1() const { return _mode1; }
    DecayMode mode2() const { return _mode2; }

    Particle t1() const { return _t1; }
    Particle t2() const { return _t2; }
    Particle b1() const { return _b1; }
    Particle b2() const { return _b2; }
    ParticleVector wDecays1() const { return _wDecays1; }
    ParticleVector wDecays2() const { return _wDecays2; }

  protected:
    // Apply the projection to the event
    void project(const Event& e);

  private:
    DecayMode _mode1, _mode2;
    Particle _t1, _t2;
    Particle _b1, _b2;
    ParticleVector _wDecays1, _wDecays2;
  };

}

#endif
