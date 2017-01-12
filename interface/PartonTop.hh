#ifndef RIVET_PartonTop_HH
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
  class PartonTop : public FinalState {
    public:
      enum TTbarMode { CH_FULLHADRON = 0, CH_SEMILEPTON, CH_FULLLEPTON };
      enum DecayMode { CH_HADRON = 0, CH_MUON, CH_ELECTRON, CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON };

      /// @name Standard constructors and destructors.
      //@{

      /// The default constructor.
      PartonTop() : FinalState(-MAXDOUBLE, MAXDOUBLE, 0.0*GeV)
      {
        setName("PartonTop");
      }
      
      /// Clone on the heap.
      virtual unique_ptr<Projection> clone() const {
        return unique_ptr<Projection>(new PartonTop(*this));
      }

      //@}

    public:
      TTbarMode mode() const {
        const bool isLepton1 = _mode1%3 != 0;
        const bool isLepton2 = _mode2%3 != 0;
        if      (  isLepton1 &&  isLepton2 ) return CH_FULLLEPTON;
        else if ( !isLepton1 && !isLepton2 ) return CH_FULLHADRON;
        return CH_SEMILEPTON;
      }
      DecayMode mode1() const { return _mode1; }
      DecayMode mode2() const { return _mode2; }

      Particle t1() const { return _t1; }
      Particle t2() const { return _t2; }
      Particle b1() const { return _b1; }
      Particle b2() const { return _b2; }
      ParticleVector wDecays1() const { return _wDecays1; }
      ParticleVector wDecays2() const { return _wDecays2; }
      Particle lepton1() const { return findLepton(_wDecays1); };
      Particle lepton2() const { return findLepton(_wDecays2); };

    protected:
      // Apply the projection to the event
      void project(const Event& e);

      Particle findLepton(const ParticleVector& v) const;

    private:
      DecayMode _mode1, _mode2;
      Particle _t1, _t2;
      Particle _b1, _b2;
      ParticleVector _wDecays1, _wDecays2;
  };

}

#endif
