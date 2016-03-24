#ifndef RIVET_PseudoBoostedTop_HH
#define RIVET_PseudoBoostedTop_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/FastJets.hh"

#include <vector>
#include <string>

namespace Rivet {

  // @brief Pseudo boosted top finder
  // 
  // Find boosted top quark in the particle level.
  class PseudoBoostedTop : public FinalState {
  public:

    // Default constructor
    PseudoBoostedTop(double lepMinPt = 45., double lepMaxEta = 2.1, double jetMinPt = 30., double jetMaxEta = 2.4, double bjetR = 0.5, double tjetR = 0.8)
      : FinalState(-MAXDOUBLE, MAXDOUBLE, 0*GeV),
        _lepMinPt(lepMinPt), _lepMaxEta(lepMaxEta), _jetMinPt(jetMinPt), _jetMaxEta(jetMaxEta), _bjetR(bjetR), _tjetR(tjetR)
    {
      setName("PseudoBoostedTop");
    }

    enum DecayMode {CH_HADRON, CH_MUON, CH_ELECTRON, CH_TAU};
    enum TTbarMode {CH_DILEP, CH_ALLHADRON, CH_SEMILEP_MU, CH_SEMILEP_EL, CH_SEMILEP_TAU};

    TTbarMode mode() const {
      if ((_topDecay == CH_MUON && _antitopDecay == CH_HADRON) || (_topDecay == CH_HADRON && _antitopDecay == CH_MUON)) return CH_SEMILEP_MU;
      else if ((_topDecay == CH_ELECTRON && _antitopDecay == CH_HADRON) || (_topDecay == CH_HADRON && _antitopDecay == CH_ELECTRON)) return CH_SEMILEP_EL;
      else if ((_topDecay == CH_TAU && _antitopDecay == CH_HADRON) || (_topDecay == CH_HADRON && _antitopDecay == CH_TAU)) return CH_SEMILEP_TAU;
      else if (_topDecay == CH_HADRON && _antitopDecay == CH_HADRON) return CH_ALLHADRON;
      else return CH_DILEP;
    }

    DecayMode topDecay() const { return _topDecay;}
    DecayMode antitopDecay() const {return _antitopDecay;}
    
    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new PseudoBoostedTop(*this);
    }

    //@}

  public:
    Particle partonTop() const   {
      if (_antitopDecay == CH_HADRON) return _antitop;
      else if (_topDecay == CH_HADRON) return _top;
      else return NULL;
    }
    float mtt() const {return _mtt;}
    Jet particleTop() const      {return _particleTop;}
    
    bool isParticleLep() const {return _isParticleLep;}
    bool isGenBJet() const {return _isGenBJet;}
    bool isGenTopJet() const {return _isGenTopJet;}
    bool passParticle() const    {return _passParticle;}

  protected:
    // Apply the projection to the event
    void project(const Event& e) override;
    bool isFromTop(const GenParticle* p);
    bool fromHardMuon(const GenParticle* p);

  private:
    const double _lepMinPt, _lepMaxEta, _jetMinPt, _jetMaxEta;
    const double _bjetR, _tjetR;

    constexpr static double _halfpi = 3.1415 / 2.;

  private:
    bool _isValid;
    
    DecayMode _topDecay;
    DecayMode _antitopDecay;

    Particle _top;
    Particle _antitop;
    
    Jet _particleTop;

    bool _isParticleLep;
    bool _isGenBJet;
    bool _isGenTopJet;
    bool _passParticle;

    float _mtt;
    
  };

}

#endif
