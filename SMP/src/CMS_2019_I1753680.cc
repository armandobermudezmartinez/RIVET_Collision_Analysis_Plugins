// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2019_I1753680 : public Analysis {
  public:
    
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2019_I1753680);
    

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState fs(Cuts::abseta < 2.4 && Cuts::pT > 25*GeV);

      ZFinder zeeFind(fs, Cuts::abseta < 2.4, PID::ELECTRON, 76.1876*GeV, 106.1876*GeV, 0.1 );
      declare(zeeFind, "ZeeFind");
      ZFinder zmmFind(fs, Cuts::abseta < 2.4, PID::MUON    , 76.1876*GeV, 106.1876*GeV, 0.1 );
      declare(zmmFind, "ZmmFind");
      
      // Book histograms
      // FIXME: use HepData file with new mapping
      book(_h_Zmm_pt            ,  2, 1, 1);
      book(_h_Zee_pt            ,  4, 1, 1);
      book(_h_Zmm_phiStar       ,  6, 1, 1);
      book(_h_Zee_phiStar       ,  8, 1, 1);
      book(_h_Zmm_absY          , 10, 1, 1);
      book(_h_Zee_absY          , 12, 1, 1);
      book(_h_Zmm_pt_Y0         , 14, 1, 1);
      book(_h_Zee_pt_Y0         , 16, 1, 1);
      book(_h_Zmm_pt_Y1         , 18, 1, 1);
      book(_h_Zee_pt_Y1         , 20, 1, 1);
      book(_h_Zmm_pt_Y2         , 22, 1, 1);
      book(_h_Zee_pt_Y2         , 24, 1, 1);
      book(_h_Zmm_pt_Y3         , 26, 1, 1);
      book(_h_Zee_pt_Y3         , 28, 1, 1);
      book(_h_Zmm_pt_Y4         , 30, 1, 1);
      book(_h_Zee_pt_Y4         , 32, 1, 1);
      book(_h_Zmm_pt_norm       , 34, 1, 1);
      book(_h_Zee_pt_norm       , 36, 1, 1);
      book(_h_Zmm_phiStar_norm  , 38, 1, 1);
      book(_h_Zee_phiStar_norm  , 40, 1, 1);
      book(_h_Zmm_absY_norm     , 42, 1, 1);
      book(_h_Zee_absY_norm     , 44, 1, 1);
      book(_h_Zmm_pt_Y0_norm    , 46, 1, 1);
      book(_h_Zee_pt_Y0_norm    , 48, 1, 1);
      book(_h_Zmm_pt_Y1_norm    , 50, 1, 1);
      book(_h_Zee_pt_Y1_norm    , 52, 1, 1);
      book(_h_Zmm_pt_Y2_norm    , 54, 1, 1);
      book(_h_Zee_pt_Y2_norm    , 56, 1, 1);
      book(_h_Zmm_pt_Y3_norm    , 58, 1, 1);
      book(_h_Zee_pt_Y3_norm    , 60, 1, 1);
      book(_h_Zmm_pt_Y4_norm    , 62, 1, 1);
      book(_h_Zee_pt_Y4_norm    , 64, 1, 1);
      
    }
    

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      const ZFinder& zeeFS = apply<ZFinder>(event, "ZeeFind");
      const ZFinder& zmumuFS = apply<ZFinder>(event, "ZmmFind");

      const Particles& zees = zeeFS.bosons();
      const Particles& zmumus = zmumuFS.bosons();

      if (zees.size() + zmumus.size() != 1) {
        MSG_DEBUG("Did not find exactly one good Z candidate");
        vetoEvent;
      }

      //event identification depending on mass window
      bool ee_event=false;
      bool mm_event=false;

      if (zees.size() == 1) { 
        ee_event = true; 
      }
      if (zmumus.size() == 1) { 
        mm_event = true; 
      }

      const Particles& theLeptons = zees.size() ? zeeFS.constituents() : zmumuFS.constituents();
      const Particle& lminus = theLeptons[0].charge() < 0 ? theLeptons[0] : theLeptons[1];
      const Particle& lplus = theLeptons[0].charge() < 0 ? theLeptons[1] : theLeptons[0];

      //calculate phi*
      const double thetaStar = acos(tanh( 0.5 * (lminus.eta() - lplus.eta()) ));
      const double dPhi = M_PI - deltaPhi(lminus, lplus);
      const double phiStar = tan(0.5 * dPhi) * sin(thetaStar);

      if (ee_event) {
        _h_Zee_pt->fill(zees[0].pt());
        _h_Zee_pt_norm->fill(zees[0].pt());
        _h_Zee_phiStar->fill(phiStar);
        _h_Zee_phiStar_norm->fill(phiStar);
        _h_Zee_absY->fill(zees[0].absrap());
        _h_Zee_absY_norm->fill(zees[0].absrap());
        if      (zees[0].absrap()<0.4) {
          _h_Zee_pt_Y0->fill(zees[0].pt());
          _h_Zee_pt_Y0_norm->fill(zees[0].pt());
        }
        else if (zees[0].absrap()<0.8) {
          _h_Zee_pt_Y1->fill(zees[0].pt());
          _h_Zee_pt_Y1_norm->fill(zees[0].pt());
        }
        else if (zees[0].absrap()<1.2) {
          _h_Zee_pt_Y2->fill(zees[0].pt());
          _h_Zee_pt_Y2_norm->fill(zees[0].pt());
        }
        else if (zees[0].absrap()<1.6) {
          _h_Zee_pt_Y3->fill(zees[0].pt());
          _h_Zee_pt_Y3_norm->fill(zees[0].pt());
        }
        else if (zees[0].absrap()<2.4) {
          _h_Zee_pt_Y4->fill(zees[0].pt());
          _h_Zee_pt_Y4_norm->fill(zees[0].pt());
        }

      } 
      else if (mm_event) {
        _h_Zmm_pt->fill(zmumus[0].pt());
        _h_Zmm_pt_norm->fill(zmumus[0].pt());
        _h_Zmm_phiStar->fill(phiStar);
        _h_Zmm_phiStar_norm->fill(phiStar);
        _h_Zmm_absY->fill(zmumus[0].absrap());
        _h_Zmm_absY_norm->fill(zmumus[0].absrap());
        if      (zmumus[0].absrap()<0.4) {
          _h_Zmm_pt_Y0->fill(zmumus[0].pt());
          _h_Zmm_pt_Y0_norm->fill(zmumus[0].pt());
        }
        else if (zmumus[0].absrap()<0.8) {
          _h_Zmm_pt_Y1->fill(zmumus[0].pt());
          _h_Zmm_pt_Y1_norm->fill(zmumus[0].pt());
        }
        else if (zmumus[0].absrap()<1.2) {
          _h_Zmm_pt_Y2->fill(zmumus[0].pt());
          _h_Zmm_pt_Y2_norm->fill(zmumus[0].pt());
        }
        else if (zmumus[0].absrap()<1.6) {
          _h_Zmm_pt_Y3->fill(zmumus[0].pt());
          _h_Zmm_pt_Y3_norm->fill(zmumus[0].pt());
        }
        else if (zmumus[0].absrap()<2.4) {
          _h_Zmm_pt_Y4->fill(zmumus[0].pt());
          _h_Zmm_pt_Y4_norm->fill(zmumus[0].pt());
	}

      }

    }

    /// Normalise histograms etc., after the run
    void finalize() {

      //normalize(_h_Z_pt);
      
      //normalize(_h_YYYY); // normalize to unity
      //scale(_h_ZZZZ, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      //
      scale(_h_Zmm_pt, crossSection()/picobarn/sumOfWeights());
      scale(_h_Zmm_absY, crossSection()/picobarn/sumOfWeights());
      scale(_h_Zmm_phiStar, crossSection()/picobarn/sumOfWeights());
      scale(_h_Zmm_pt_Y0, crossSection()/picobarn/sumOfWeights());
      scale(_h_Zmm_pt_Y1, crossSection()/picobarn/sumOfWeights());
      scale(_h_Zmm_pt_Y2, crossSection()/picobarn/sumOfWeights());
      scale(_h_Zmm_pt_Y3, crossSection()/picobarn/sumOfWeights());
      scale(_h_Zmm_pt_Y4, crossSection()/picobarn/sumOfWeights());
      
      scale(_h_Zee_pt, crossSection()/picobarn/sumOfWeights());
      scale(_h_Zee_absY, crossSection()/picobarn/sumOfWeights());
      scale(_h_Zee_phiStar, crossSection()/picobarn/sumOfWeights());
      scale(_h_Zee_pt_Y0, crossSection()/picobarn/sumOfWeights());
      scale(_h_Zee_pt_Y1, crossSection()/picobarn/sumOfWeights());
      scale(_h_Zee_pt_Y2, crossSection()/picobarn/sumOfWeights());
      scale(_h_Zee_pt_Y3, crossSection()/picobarn/sumOfWeights());
      scale(_h_Zee_pt_Y4, crossSection()/picobarn/sumOfWeights());

      normalize(_h_Zmm_pt_norm);
      normalize(_h_Zmm_absY_norm);
      normalize(_h_Zmm_phiStar_norm);
      normalize(_h_Zmm_pt_Y0_norm);
      normalize(_h_Zmm_pt_Y1_norm);
      normalize(_h_Zmm_pt_Y2_norm);
      normalize(_h_Zmm_pt_Y3_norm);
      normalize(_h_Zmm_pt_Y4_norm);
      
      normalize(_h_Zee_pt_norm);
      normalize(_h_Zee_absY_norm);
      normalize(_h_Zee_phiStar_norm);
      normalize(_h_Zee_pt_Y0_norm);
      normalize(_h_Zee_pt_Y1_norm);
      normalize(_h_Zee_pt_Y2_norm);
      normalize(_h_Zee_pt_Y3_norm);
      normalize(_h_Zee_pt_Y4_norm);

    }

    //@}


    /// @name Histograms

  private:
    
    Histo1DPtr   _h_Zmm_pt, _h_Zmm_phiStar, _h_Zmm_absY;
    Histo1DPtr   _h_Zmm_pt_Y0, _h_Zmm_pt_Y1, _h_Zmm_pt_Y2, _h_Zmm_pt_Y3, _h_Zmm_pt_Y4;
    	
    Histo1DPtr   _h_Zmm_pt_norm, _h_Zmm_phiStar_norm, _h_Zmm_absY_norm;
    Histo1DPtr   _h_Zmm_pt_Y0_norm, _h_Zmm_pt_Y1_norm, _h_Zmm_pt_Y2_norm, _h_Zmm_pt_Y3_norm, _h_Zmm_pt_Y4_norm;
	
    Histo1DPtr   _h_Zee_pt, _h_Zee_phiStar, _h_Zee_absY;
    Histo1DPtr   _h_Zee_pt_Y0, _h_Zee_pt_Y1, _h_Zee_pt_Y2, _h_Zee_pt_Y3, _h_Zee_pt_Y4;
	
    Histo1DPtr   _h_Zee_pt_norm, _h_Zee_phiStar_norm, _h_Zee_absY_norm;
    Histo1DPtr   _h_Zee_pt_Y0_norm, _h_Zee_pt_Y1_norm, _h_Zee_pt_Y2_norm, _h_Zee_pt_Y3_norm, _h_Zee_pt_Y4_norm;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2019_I1753680);


}
