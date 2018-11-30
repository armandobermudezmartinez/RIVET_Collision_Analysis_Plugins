// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/ProjectionApplier.hh"

namespace Rivet {


  class CMS_2016_I1426517 : public Analysis {
  public:

    /// Constructor
    CMS_2016_I1426517()
      : Analysis("CMS_2016_I1426517")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// @todo Initialise and register projections here
      //FinalState fs(Cuts::abseta < 2.4);
      FinalState fs;
      //declare(fs, "FS");
      addProjection(fs, "FS");

      /// Get W's with muons with |eta| < 2.4, pT > 25 GeV
      WFinder wfinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 25.0*GeV, PID::MUON, 0.0*GeV, MAXDOUBLE, 0.0*GeV, 0.1, 
		      WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);

        addProjection(wfinder, "WFinder");
    //  declare(wfinder, "WFinder");


      /// @todo Book histograms here, e.g.:
      _h_dsigmadeta_Wplus  = bookHisto1D(1,1,1);
      _h_dsigmadeta_Wminus = bookHisto1D(2,1,1);
      _h_asymmetry = bookScatter2D(3,1,1);
      
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      /// @todo Do the event by event analysis here

      // Get W boson 
      const WFinder& wfinder = applyProjection<WFinder>(event, "WFinder");
      if (wfinder.bosons().size() != 1) vetoEvent;
      //const Particle& w = wfinder.bosons()[0];
      
      // Get lepton from W
      const Particle& l = wfinder.constituentLeptons()[0];
      double _lep_abseta = l.abseta();
      int _lep_charge = l.charge();
      if (l.pT() < 25*GeV || _lep_abseta > 2.4) vetoEvent;

      // Fill histograms
      if(_lep_charge == 1){
	_h_dsigmadeta_Wplus->fill(_lep_abseta, weight);
      }

      if(_lep_charge == -1){
	_h_dsigmadeta_Wminus->fill(_lep_abseta, weight);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      /// @todo Normalise, scale and otherwise manipulate histograms here
      double sf = crossSection()/picobarn/sumOfWeights();
      scale(_h_dsigmadeta_Wplus, sf);
      scale(_h_dsigmadeta_Wminus, sf);
      // scale(_h_YYYY, crossSection()/sumOfWeights()); // norm to cross section
      // normalize(_h_YYYY); // normalize to unity
      assert(_h_dsigmadeta_Wplus->numBins() == _h_dsigmadeta_Wminus->numBins());
      for (size_t i = 0; i < _h_dsigmadeta_Wplus->numBins(); ++i) {
	const double num   = _h_dsigmadeta_Wplus->bin(i).sumW() - _h_dsigmadeta_Wminus->bin(i).sumW();
	const double denom = _h_dsigmadeta_Wplus->bin(i).sumW() + _h_dsigmadeta_Wminus->bin(i).sumW();
	const double relerr = _h_dsigmadeta_Wplus->bin(i).relErr()  + _h_dsigmadeta_Wminus->bin(i).relErr();
	const double asym = (num != 0 && denom != 0) ? num / denom : 0;
	const double asym_err = (num != 0 && denom != 0) ? asym*relerr : 0;
	_h_asymmetry->addPoint(_h_dsigmadeta_Wplus->bin(i).xMid(), asym, _h_dsigmadeta_Wplus->bin(i).xWidth()/2.0, asym_err);
      }

    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here

    /// @name Histograms
    //@{
    Histo1DPtr _h_dsigmadeta_Wplus;
    Histo1DPtr _h_dsigmadeta_Wminus;
    Scatter2DPtr _h_asymmetry;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1426517);


}
