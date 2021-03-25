// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2020_I1837084 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2020_I1837084);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // Initialise and register projections
      FinalState fs;
      Cut cut = Cuts::pT > 0*GeV;

      ZFinder zmmFind(fs, cut, PID::MUON, 76.1876*GeV, 106.1876*GeV, 0.1, ZFinder::ChargedLeptons::PROMPT, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::YES );
      declare(zmmFind, "ZmmFind");
      
      // Book histograms
      book(_h_Z_pt,      4, 1, 5);
      book(_h_Z_pt_norm, 5, 1, 5);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ZFinder& zmmFS = apply<ZFinder>(event, "ZmmFind");

      const Particles& zmms = zmmFS.bosons();

      if (zmms.size() == 1 && zmms[0].pt() > 200) { 
        _h_Z_pt     ->fill(min(zmms[0].pt(),1499.999));
        _h_Z_pt_norm->fill(min(zmms[0].pt(),1499.999));
      }

    }

    void normalizeToSum(Histo1DPtr hist) {
      double sum = 0.;
      for (size_t i = 0; i < hist->numBins(); ++i) {
        sum += hist->bin(i).height();
      }
      scale(hist, 1./sum);
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      double norm = (sumOfWeights() != 0) ? crossSection()/femtobarn/sumOfWeights() : 1.0;
      
      scale(_h_Z_pt, norm);

      normalizeToSum(_h_Z_pt_norm);

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_Z_pt, _h_Z_pt_norm;
    //@}


  };


  DECLARE_RIVET_PLUGIN(CMS_2020_I1837084);

}
