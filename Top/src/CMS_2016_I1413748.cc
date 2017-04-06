// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2016_I1413748 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2016_I1413748);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // Complete final state
      FinalState fs(-MAXDOUBLE, MAXDOUBLE, 0*GeV);

      // Projection for electrons and muons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);
      IdentifiedFinalState l_id(fs);
      l_id.acceptIdPair(PID::ELECTRON);
      l_id.acceptIdPair(PID::MUON);
      PromptFinalState leptons(l_id);
      addProjection(leptons, "Leptons");
      DressedLeptons dressedleptons(photons, leptons, 0.1, Cuts::open(), true, false);
      addProjection(dressedleptons, "DressedLeptons");

      // Weight counter
      _vis_unit_weights = 0.;

      // Booking of histograms
      const std::vector<double> bindphi = {0., 5.*M_PI/60., 10.*M_PI/60., 15.*M_PI/60., 20.*M_PI/60., 25.*M_PI/60., 30.*M_PI/60., 35.*M_PI/60., 40.*M_PI/60., 45.*M_PI/60., 50.*M_PI/60., 55.*M_PI/60., M_PI};
      _h_dphi = bookHisto1D("dphi", bindphi);
      
      const std::vector<double> bindabseta = { -2., -68./60., -48./60., -32./60., -20./60., -8./60., 0., 8./60., 20./60., 32./60., 48./60., 68./60., 2.};
      _h_dabseta = bookHisto1D("dabseta", bindabseta);

    }






    /// Perform the per-event analysis
    void analyze(const Event& event) {

      double EPSILON = 0.00000001;
      
      const double weight = event.weight();
      
      // select ttbar -> lepton+jets without taus
      const DressedLeptons& dressedleptons = applyProjection<DressedLeptons>(event, "DressedLeptons");
      if (dressedleptons.dressedLeptons().size() != 2) vetoEvent;

      // count weights
      if (weight != 0.) _vis_unit_weights += weight/std::abs(weight);

      if (sameSign(dressedleptons.dressedLeptons()[0],dressedleptons.dressedLeptons()[1])) {
        cout<<"error in lepton charge assignment"<<endl;
      }

      FourMomentum lepPlus = dressedleptons.dressedLeptons()[0].charge() > 0 ? dressedleptons.dressedLeptons()[0].momentum() : dressedleptons.dressedLeptons()[1].momentum();
      FourMomentum lepMinus = dressedleptons.dressedLeptons()[0].charge() > 0 ? dressedleptons.dressedLeptons()[1].momentum() : dressedleptons.dressedLeptons()[0].momentum();
      
      double dabseta_unbounded = lepPlus.abseta() - lepMinus.abseta();

      _h_dphi->fill(min(deltaPhi(lepPlus,lepMinus), M_PI-EPSILON), weight);
      _h_dabseta->fill(dabseta_unbounded > 0 ? min(dabseta_unbounded, 2.-EPSILON) : max(dabseta_unbounded, -2.+EPSILON), weight);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_dphi); // normalize to unity
      //scale(_h_dphi, crossSection()/picobarn/sumOfWeights()); // norm to cross section

      normalize(_h_dabseta); // normalize to unity
      //scale(_h_dabseta, crossSection()/picobarn/sumOfWeights()); // norm to cross section

    }

    //@}


  private:


    /// @name Histograms
    //@{
    double _vis_unit_weights;
    Histo1DPtr _h_dphi, _h_dabseta;
    //Histo1DPtr _h_mode0, _h_mode1, _h_mode2;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1413748);


}
