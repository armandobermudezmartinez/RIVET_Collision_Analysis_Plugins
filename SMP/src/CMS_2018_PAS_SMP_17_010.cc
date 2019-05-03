// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2018_PAS_SMP_17_010 : public Analysis {
  public:
    
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2018_PAS_SMP_17_010);
    

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState fs(Cuts::abseta < 2.4 && Cuts::pT > 25*GeV);

      ZFinder zeeFind(fs, Cuts::abseta < 2.4, PID::ELECTRON, 76.1876*GeV, 106.1876*GeV, 0.1 );
      addProjection(zeeFind, "ZeeFind");
      ZFinder zmmFind(fs, Cuts::abseta < 2.4, PID::MUON    , 76.1876*GeV, 106.1876*GeV, 0.1 );
      addProjection(zmmFind, "ZmmFind");
      
      // Book histograms
      //_h_Zmm_pt_cov         = bookHisto2D(1,1,1);
      _h_Zmm_pt             = bookHisto1D(2,1,1);
      //_h_Zee_pt_cov         = bookHisto2D(3,1,1);
      _h_Zee_pt             = bookHisto1D(4,1,1);
      //_h_Zmm_phiStar_cov    = bookHisto2D(5,1,1);
      _h_Zmm_phiStar        = bookHisto1D(6,1,1);
      //_h_Zee_phiStar_cov    = bookHisto2D(7,1,1);
      _h_Zee_phiStar        = bookHisto1D(8,1,1);
      //_h_Zmm_absY_cov       = bookHisto2D(9,1,1);
      _h_Zmm_absY           = bookHisto1D(10,1,1);
      //_h_Zee_absY_cov       = bookHisto2D(11,1,1);
      _h_Zee_absY           = bookHisto1D(12,1,1);
      //_h_Zmm_pt_Y0_cov      = bookHisto2D(13,1,1);
      _h_Zmm_pt_Y0          = bookHisto1D(14,1,1);
      //_h_Zee_pt_Y0_cov      = bookHisto2D(15,1,1);
      _h_Zee_pt_Y0          = bookHisto1D(16,1,1);
      //_h_Zmm_pt_Y1_cov      = bookHisto2D(17,1,1);
      _h_Zmm_pt_Y1          = bookHisto1D(18,1,1);
      //_h_Zee_pt_Y1_cov      = bookHisto2D(19,1,1);
      _h_Zee_pt_Y1          = bookHisto1D(20,1,1);
      //_h_Zmm_pt_Y2_cov      = bookHisto2D(21,1,1);
      _h_Zmm_pt_Y2          = bookHisto1D(22,1,1);
      //_h_Zee_pt_Y2_cov      = bookHisto2D(23,1,1);
      _h_Zee_pt_Y2          = bookHisto1D(24,1,1);
      //_h_Zmm_pt_Y3_cov      = bookHisto2D(25,1,1);
      _h_Zmm_pt_Y3          = bookHisto1D(26,1,1);
      //_h_Zee_pt_Y3_cov      = bookHisto2D(27,1,1);
      _h_Zee_pt_Y3          = bookHisto1D(28,1,1);
      //_h_Zmm_pt_Y4_cov      = bookHisto2D(29,1,1);
      _h_Zmm_pt_Y4          = bookHisto1D(30,1,1);
      //_h_Zee_pt_Y4_cov      = bookHisto2D(31,1,1);
      _h_Zee_pt_Y4          = bookHisto1D(32,1,1);
      //_h_Zmm_pt_norm_cov    = bookHisto2D(33,1,1);
      _h_Zmm_pt_norm        = bookHisto1D(34,1,1);
      //_h_Zee_pt_norm_cov    = bookHisto2D(35,1,1);
      _h_Zee_pt_norm        = bookHisto1D(36,1,1);
      //_h_Zmm_phiStar_norm_cov = bookHisto2D(37,1,1);
      _h_Zmm_phiStar_norm     = bookHisto1D(38,1,1);
      //_h_Zee_phiStar_norm_cov = bookHisto2D(39,1,1);
      _h_Zee_phiStar_norm     = bookHisto1D(40,1,1);
      //_h_Zmm_absY_norm_cov    = bookHisto2D(41,1,1);
      _h_Zmm_absY_norm        = bookHisto1D(42,1,1);
      //_h_Zee_absY_norm_cov    = bookHisto2D(43,1,1);
      _h_Zee_absY_norm        = bookHisto1D(44,1,1);
      //_h_Zmm_pt_Y0_norm_cov   = bookHisto2D(45,1,1);
      _h_Zmm_pt_Y0_norm       = bookHisto1D(46,1,1);
      //_h_Zee_pt_Y0_norm_cov   = bookHisto2D(47,1,1);
      _h_Zee_pt_Y0_norm       = bookHisto1D(48,1,1);
      //_h_Zmm_pt_Y1_norm_cov   = bookHisto2D(49,1,1);
      _h_Zmm_pt_Y1_norm       = bookHisto1D(50,1,1);
      //_h_Zee_pt_Y1_norm_cov   = bookHisto2D(51,1,1);
      _h_Zee_pt_Y1_norm       = bookHisto1D(52,1,1);
      //_h_Zmm_pt_Y2_norm_cov   = bookHisto2D(53,1,1);
      _h_Zmm_pt_Y2_norm       = bookHisto1D(54,1,1);
      //_h_Zee_pt_Y2_norm_cov   = bookHisto2D(55,1,1);
      _h_Zee_pt_Y2_norm       = bookHisto1D(56,1,1);
      //_h_Zmm_pt_Y3_norm_cov   = bookHisto2D(57,1,1);
      _h_Zmm_pt_Y3_norm       = bookHisto1D(58,1,1);
      //_h_Zee_pt_Y3_norm_cov   = bookHisto2D(59,1,1);
      _h_Zee_pt_Y3_norm       = bookHisto1D(60,1,1);
      //_h_Zmm_pt_Y4_norm_cov   = bookHisto2D(61,1,1);
      _h_Zmm_pt_Y4_norm       = bookHisto1D(62,1,1);
      //_h_Zee_pt_Y4_norm_cov   = bookHisto2D(63,1,1);
      _h_Zee_pt_Y4_norm       = bookHisto1D(64,1,1);
      
    }
    

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      /// @todo Do the event by event analysis here
      
      //Event weight
      const double w = event.weight();
      
      const ZFinder& zeeFS = applyProjection<ZFinder>(event, "ZeeFind");
      const ZFinder& zmumuFS = applyProjection<ZFinder>(event, "ZmmFind");

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
      const double thetaStar = tanh( 0.5 * (lminus.eta() - lplus.eta()) );
      const double dPhi = M_PI - deltaPhi(lminus, lplus);
      const double phiStar = tan(0.5 * dPhi) * sin(thetaStar);

      if (ee_event) {
	_h_Zee_pt->fill(zees[0].pt(),w);
	_h_Zee_pt_norm->fill(zees[0].pt(),w);
	_h_Zee_phiStar->fill(phiStar,w);
	_h_Zee_phiStar_norm->fill(phiStar,w);
	_h_Zee_absY->fill(zees[0].absrap(),w);
	_h_Zee_absY_norm->fill(zees[0].absrap(),w);
	if      (zees[0].absrap()<0.4) {
	  _h_Zee_pt_Y0->fill(zees[0].pt(),w);
	  _h_Zee_pt_Y0_norm->fill(zees[0].pt(),w);
	}
	else if (zees[0].absrap()<0.8) {
	  _h_Zee_pt_Y1->fill(zees[0].pt(),w);
	  _h_Zee_pt_Y1_norm->fill(zees[0].pt(),w);
	}
	else if (zees[0].absrap()<1.2) {
	  _h_Zee_pt_Y2->fill(zees[0].pt(),w);
	  _h_Zee_pt_Y2_norm->fill(zees[0].pt(),w);
	}
	else if (zees[0].absrap()<1.6) {
	  _h_Zee_pt_Y3->fill(zees[0].pt(),w);
	  _h_Zee_pt_Y3_norm->fill(zees[0].pt(),w);
	}
	else if (zees[0].absrap()<2.4) {
	  _h_Zee_pt_Y4->fill(zees[0].pt(),w);
	  _h_Zee_pt_Y4_norm->fill(zees[0].pt(),w);
	}

      } 
      else if (mm_event) {
	_h_Zmm_pt->fill(zmumus[0].pt(),w);
	_h_Zmm_pt_norm->fill(zmumus[0].pt(),w);
	_h_Zmm_phiStar->fill(phiStar,w);
	_h_Zmm_phiStar_norm->fill(phiStar,w);
	_h_Zmm_absY->fill(zmumus[0].absrap(),w);
	_h_Zmm_absY_norm->fill(zmumus[0].absrap(),w);
        if      (zmumus[0].absrap()<0.4) {
          _h_Zmm_pt_Y0->fill(zmumus[0].pt(),w);
          _h_Zmm_pt_Y0_norm->fill(zmumus[0].pt(),w);
        }
        else if (zmumus[0].absrap()<0.8) {
          _h_Zmm_pt_Y1->fill(zmumus[0].pt(),w);
          _h_Zmm_pt_Y1_norm->fill(zmumus[0].pt(),w);
        }
        else if (zmumus[0].absrap()<1.2) {
          _h_Zmm_pt_Y2->fill(zmumus[0].pt(),w);
          _h_Zmm_pt_Y2_norm->fill(zmumus[0].pt(),w);
        }
        else if (zmumus[0].absrap()<1.6) {
          _h_Zmm_pt_Y3->fill(zmumus[0].pt(),w);
          _h_Zmm_pt_Y3_norm->fill(zmumus[0].pt(),w);
        }
        else if (zmumus[0].absrap()<2.4) {
          _h_Zmm_pt_Y4->fill(zmumus[0].pt(),w);
          _h_Zmm_pt_Y4_norm->fill(zmumus[0].pt(),w);
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
    Histo1DPtr   _h_Zmm_pt_Y0, _h_Zmm_pt_Y1, _h_Zmm_pt_Y2, _h_Zmm_pt_Y3, _h_Zmm_pt_Y4; //8
    //Histo2DPtr   _h_Zmm_pt_cov, _h_Zmm_phiStar_cov, _h_Zmm_absY_cov;
    //Histo2DPtr   _h_Zmm_pt_Y0_cov, _h_Zmm_pt_Y1_cov, _h_Zmm_pt_Y2_cov, _h_Zmm_pt_Y3_cov, _h_Zmm_pt_Y4_cov; //8
    	
    Histo1DPtr   _h_Zmm_pt_norm, _h_Zmm_phiStar_norm, _h_Zmm_absY_norm;
    Histo1DPtr   _h_Zmm_pt_Y0_norm, _h_Zmm_pt_Y1_norm, _h_Zmm_pt_Y2_norm, _h_Zmm_pt_Y3_norm, _h_Zmm_pt_Y4_norm; //8
    //Histo2DPtr   _h_Zmm_pt_norm_cov, _h_Zmm_phiStar_norm_cov, _h_Zmm_absY_norm_cov;
    //Histo2DPtr   _h_Zmm_pt_Y0_norm_cov, _h_Zmm_pt_Y1_norm_cov, _h_Zmm_pt_Y2_norm_cov, _h_Zmm_pt_Y3_norm_cov, _h_Zmm_pt_Y4_norm_cov; //8
	
    Histo1DPtr   _h_Zee_pt, _h_Zee_phiStar, _h_Zee_absY;
    Histo1DPtr   _h_Zee_pt_Y0, _h_Zee_pt_Y1, _h_Zee_pt_Y2, _h_Zee_pt_Y3, _h_Zee_pt_Y4; //8
    //Histo2DPtr   _h_Zee_pt_cov, _h_Zee_phiStar_cov, _h_Zee_absY_cov;
    //Histo2DPtr   _h_Zee_pt_Y0_cov, _h_Zee_pt_Y1_cov, _h_Zee_pt_Y2_cov, _h_Zee_pt_Y3_cov, _h_Zee_pt_Y4_cov; //8
	
    Histo1DPtr   _h_Zee_pt_norm, _h_Zee_phiStar_norm, _h_Zee_absY_norm;
    Histo1DPtr   _h_Zee_pt_Y0_norm, _h_Zee_pt_Y1_norm, _h_Zee_pt_Y2_norm, _h_Zee_pt_Y3_norm, _h_Zee_pt_Y4_norm; //8
    //Histo2DPtr   _h_Zee_pt_norm_cov, _h_Zee_phiStar_norm_cov, _h_Zee_absY_norm_cov;
    //Histo2DPtr   _h_Zee_pt_Y0_norm_cov, _h_Zee_pt_Y1_norm_cov, _h_Zee_pt_Y2_norm_cov, _h_Zee_pt_Y3_norm_cov, _h_Zee_pt_Y4_norm_cov; //8






  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2018_PAS_SMP_17_010);


}
