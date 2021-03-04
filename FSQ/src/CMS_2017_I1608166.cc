// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/ParticleName.hh"

namespace Rivet {


  /// @brief Measurement of charged pion, kaon, and proton production in proton-proton collisions at 13 TeV
  class CMS_2017_I1608166 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2017_I1608166);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs(Cuts::absrap < 1.);
      declare(cfs, "CFS");
      //
      // pt spectra
      book(_h[PID::PIPLUS],  "d01-x01-y01");
      book(_h[PID::KPLUS],   "d01-x01-y02");
      book(_h[PID::PROTON],  "d01-x01-y03");
      book(_h[PID::PIMINUS], "d02-x01-y01");
      book(_h[PID::KMINUS],  "d02-x01-y02");
      book(_h[PID::PBAR],    "d02-x01-y03");
      // negative/positive ratios
      book(_s["pi-/pi+"], "d100-x01-y01");
      book(_s["k-/k+"],   "d100-x01-y02");
      book(_s["p~/p"],    "d100-x01-y03");
      // k/pi and p/pi ratios
      book(_hkpi[PID::PIPLUS], "TMP/hkpi/pi", refData(101, 1, 1));
      book(_hkpi[PID::KPLUS],  "TMP/hkpi/k",  refData(101, 1, 1));
      book(_hppi[PID::PIPLUS], "TMP/hppi/pi", refData(101, 1, 2));
      book(_hppi[PID::PROTON], "TMP/hppi/p",  refData(101, 1, 2));
      book(_s["k/pi"],    "d101-x01-y01");
      book(_s["p/pi"],    "d101-x01-y02");
    }


    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      for (const Particle& p : cfs.particles()) {
        // protections against mc generators decaying long-lived particles
        if (p.hasAncestor(310)  || p.hasAncestor(-310)  ||  // K0s
            p.hasAncestor(130)  || p.hasAncestor(-130)  ||  // K0l
            p.hasAncestor(3322) || p.hasAncestor(-3322) ||  // Xi0
            p.hasAncestor(3122) || p.hasAncestor(-3122) ||  // Lambda
            p.hasAncestor(3222) || p.hasAncestor(-3222) ||  // Sigma+/-
            p.hasAncestor(3312) || p.hasAncestor(-3312) ||  // Xi-/+
            p.hasAncestor(3334) || p.hasAncestor(-3334))    // Omega-/+
          continue;
        
        if (theParticles.find(p.pid()) != theParticles.end()) {
          // fill pt spectra
          _h[p.pid()]->fill(p.pt() / GeV);
          // fill tmp histos for ratios
          if (p.abspid() != PID::PROTON)
            _hkpi[p.abspid()]->fill(p.pt() / GeV);
          if (p.abspid() != PID::KPLUS)
            _hppi[p.abspid()]->fill(p.pt() / GeV);
        }
      }
    }

    void finalize() {

      divide(_h[PID::PIMINUS], _h[PID::PIPLUS], _s["pi-/pi+"]);
      divide(_h[PID::KMINUS],  _h[PID::KPLUS],  _s["k-/k+"]);
      divide(_h[PID::PBAR],    _h[PID::PROTON], _s["p~/p"]);
      
      divide(_hkpi[PID::KPLUS],  _hkpi[PID::PIPLUS], _s["k/pi"]);
      divide(_hppi[PID::PROTON], _hppi[PID::PIPLUS], _s["p/pi"]);
      
      scale({_h[PID::PIPLUS], _h[PID::KPLUS], _h[PID::PROTON], _h[PID::PIMINUS], _h[PID::KMINUS], _h[PID::PBAR], }, 1./2./sumOfWeights());

    }

  private:
  
    set<int> theParticles = {PID::PIPLUS, PID::KPLUS, PID::PROTON, PID::PIMINUS, PID::KMINUS, PID::PBAR};
  
    map<int, Histo1DPtr> _h;
    map<int, Histo1DPtr> _hkpi, _hppi;
    map<string, Scatter2DPtr> _s;


  };


  DECLARE_RIVET_PLUGIN(CMS_2017_I1608166);

}