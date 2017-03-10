// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include <iostream>
namespace Rivet {

   /// CMS Azimuthal deccorellations at 8TeV
   class CMS_2016_I1421646 : public Analysis {
  public:

    CMS_2016_I1421646() : Analysis("CMS_2016_I1421646") {}


    void init() {
      FinalState fs;
      FastJets akt(fs, FastJets::ANTIKT, 0.7);
      addProjection(akt, "antikT");

      _h_deltaPhi.addHistogram( 200.,  300., bookHisto1D(1, 1, 1));
      _h_deltaPhi.addHistogram( 300.,  400., bookHisto1D(2, 1, 1));
      _h_deltaPhi.addHistogram( 400.,  500., bookHisto1D(3, 1, 1));
      _h_deltaPhi.addHistogram( 500.,  700., bookHisto1D(4, 1, 1));
      _h_deltaPhi.addHistogram( 700.,  900., bookHisto1D(5, 1, 1));
      _h_deltaPhi.addHistogram( 900.,  1100., bookHisto1D(6, 1, 1));
      _h_deltaPhi.addHistogram( 1100., 4000., bookHisto1D(7, 1, 1));

    }


    void analyze(const Event & event) {

      const double weight = event.weight();

      const Jets& jets = applyProjection<JetAlg>(event, "antikT").jetsByPt();
      
      if ( jets.size() < 2 ) vetoEvent;

      if ( fabs(jets[0].rap()) > 2.5 || jets[0].pT() < 200.*GeV ) vetoEvent;
      
      if ( fabs(jets[1].rap()) > 2.5 || jets[1].pT() < 100.*GeV ) vetoEvent;
      
      double dphi = deltaPhi(jets[0].phi(), jets[1].phi());
      
      _h_deltaPhi.fill(jets[0].pT(), dphi, weight);
    }


    void finalize() {
      
      foreach (Histo1DPtr histo, _h_deltaPhi.getHistograms()) {
        normalize(histo); 
      }
    }

  private:

    BinnedHistogram<double> _h_deltaPhi;

  };

  // A hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1421646);

}
