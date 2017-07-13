// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// Measurement of the W and Z bosons pT produced in pp collisions at sqrt(s)=8 TeV
  class CMS_2017_I1471281 : public Analysis {
    public:

      /// Constructor
      CMS_2017_I1471281(std::string name="CMS_2017_I1471281")
	: Analysis(name)
	{
	  _mode = 0; // init 
	  //setNeedsCrossSection(true);
	}


      /// @name Analysis methods
      //@{

      /// Book histograms and initialise projections before the run
      void init() {

	// Set up projections
	FinalState fs(-MAXDOUBLE, MAXDOUBLE, 0*GeV);
	addProjection(fs, "FS");

	Cut cut_el = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;
	Cut cut_mu = Cuts::abseta < 2.1 && Cuts::pT > 20*GeV;

	// Dressed Ws ...
	WFinder welnu_Finder(fs, cut_el, PID::ELECTRON, 0*GeV, MAXDOUBLE, 0*GeV, 12, WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);
	addProjection(welnu_Finder, "Welnu_Finder");

	WFinder wmunu_Finder(fs, cut_mu, PID::MUON, 0*GeV, MAXDOUBLE, 0*GeV, 12, WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);
	addProjection(wmunu_Finder, "Wmunu_Finder");

	// Dressed Zs ... 
	ZFinder zmumu_Finder(fs, cut_mu, PID::MUON, 60*GeV, 120*GeV, 12, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
	addProjection(zmumu_Finder, "Zmumu_Finder");

	// Histograms
	if ( _mode == 1) {
	  _hist_WtoElNuPt = bookHisto1D(1, 1, 1);
	  _hist_WtoMuNuPt = bookHisto1D(2, 1, 1);
	  //_hist_tmp_WmToMuNuPt = bookHisto1D("TMP/WmToMuNuPt",refData(4, 1, 1));
	  //_hist_tmp_WpToMuNuPt = bookHisto1D("TMP/WpToMuNuPt",refData(4, 1, 1));
	  //_hist_WmToWpRatio = bookScatter2D(4, 1, 1);
	}
	if( _mode == 2) {
	  _hist_ZtoMuMuPt = bookHisto1D(3, 1, 1);
	}
      }


      /// Perform the per-event analysis
      void analyze(const Event& event) {

	const double weight = event.weight();
	//int    lepton_charge = 0;

	if ( _mode == 1 ) {
	  // Get the W bosons - electron decay channel
	  const WFinder& welnu_Finder = applyProjection<WFinder>(event, "Welnu_Finder");
	  if (!welnu_Finder.bosons().empty()) {
	    const FourMomentum pWelnu = welnu_Finder.bosons()[0].momentum();
	    _hist_WtoElNuPt->fill(pWelnu.pT()/GeV, weight);
	  }

	  // Get the W bosons - muon decay channel
	  const WFinder& wmunu_Finder = applyProjection<WFinder>(event, "Wmunu_Finder");
	  if (!wmunu_Finder.bosons().empty()) {
	    //lepton_charge = wmunu_Finder.constituentLeptons()[0].charge();
	    const FourMomentum pWmunu = wmunu_Finder.bosons()[0].momentum();
	    _hist_WtoMuNuPt->fill(pWmunu.pT()/GeV, weight);
	    /***
	    //Get negative charged the W bosons - muon decay channel
	    if (lepton_charge < 0) {
	      _hist_tmp_WmToMuNuPt->fill(pWmunu.pT()/GeV, weight);
	    }
	    //Get positive charged the W bosons - muon decay channel
	    if (lepton_charge > 0) {
	      _hist_tmp_WpToMuNuPt->fill(pWmunu.pT()/GeV, weight);
	    }
	    ***/
	  }
	}

	else if ( _mode == 2 ) {
	  // Get the Z bosons - muon decay channel
	  const ZFinder& zmumu_Finder = applyProjection<ZFinder>(event, "Zmumu_Finder");
	  if (!zmumu_Finder.bosons().empty()) {
	    const FourMomentum pZmumu = zmumu_Finder.bosons()[0].momentum();
	    _hist_ZtoMuMuPt->fill(pZmumu.pT()/GeV, weight);
	  }
	}

      }


      /// Normalise histograms etc., after the run
      void finalize() {

	MSG_INFO("Cross section = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << crossSection() << " pb");
	MSG_INFO("# Events      = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << numEvents() );
	MSG_INFO("SumW          = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << sumOfWeights());

	if (_mode ==1) {
	  normalize(_hist_WtoElNuPt);
	  normalize(_hist_WtoMuNuPt);
	  //normalize(_hist_tmp_WmToMuNuPt);
	  //normalize(_hist_tmp_WpToMuNuPt);
	
	  // Create ratio histograms
	  //divide(_hist_tmp_WmToMuNuPt,_hist_tmp_WpToMuNuPt, _hist_WmToWpRatio);
	}

	else if (_mode ==2) {
	  normalize(_hist_ZtoMuMuPt);
	}
	      

      }

      //@}

    protected:

      // Data members
      size_t _mode;


    private:


      /// @name Histograms

      Histo1DPtr _hist_WtoElNuPt;
      Histo1DPtr _hist_WtoMuNuPt;
      Histo1DPtr _hist_ZtoMuMuPt;

      Scatter2DPtr _hist_WmToWpRatio;

      // Temporary histograms to be divided for the ratios
      //Histo1DPtr _hist_tmp_WmToMuNuPt, _hist_tmp_WpToMuNuPt;

  };

  /// W part of the CMS_2017_I1471281 analysis
  class CMS_2017_I1471281_W : public CMS_2017_I1471281 {
    public:
      CMS_2017_I1471281_W()
	: CMS_2017_I1471281("CMS_2017_I1471281_W")
      {_mode = 1;}
  };

  /// Z part of the CMS_2017_I1471281 analysis
  class CMS_2017_I1471281_Z : public CMS_2017_I1471281 {
    public:
      CMS_2017_I1471281_Z()
	: CMS_2017_I1471281("CMS_2017_I1471281_Z")
      {_mode = 2;}
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_I1471281);
  DECLARE_RIVET_PLUGIN(CMS_2017_I1471281_W);
  DECLARE_RIVET_PLUGIN(CMS_2017_I1471281_Z);


}
