#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/ParticleName.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {

  class CMS_2016_I1473674 : public Analysis {
  public:

    /// Minimal constructor
    CMS_2016_I1473674() : Analysis("CMS_2016_I1473674")
    {
    }


    /// @name Analysis methods
    //@{

    /// Set up projections and book histograms
    void init() {
      // Complete final state
      double MAXRAPIDITY = 1e10;
      FinalState fs(-MAXRAPIDITY, MAXRAPIDITY, 0*GeV);
      
      // Projection for taus
      TauFinder taus(TauFinder::ANY);
      addProjection(taus, "Tau");
      IdentifiedFinalState nu_taus(fs);
      nu_taus.acceptIdPair(PID::NU_TAU);
      addProjection(nu_taus, "nu_tau");

      // Projection for electrons and muons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);
      
      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      addProjection(electrons, "Electrons");
      DressedLeptons dressedelectrons(photons, electrons, 0.1, Cuts::open(), true, false);
      addProjection(dressedelectrons, "DressedElectrons");
      
      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      addProjection(muons, "Muons");
      DressedLeptons dressedmuons(photons, muons, 0.1, Cuts::open(), true, false);
      addProjection(dressedmuons, "DressedMuons");
      
      // Projection for jets
      VetoedFinalState fsForJets(FinalState(-MAXRAPIDITY, MAXRAPIDITY, 0*GeV));
      fsForJets.addVetoOnThisFinalState(dressedmuons);
      addProjection(FastJets(fsForJets, FastJets::ANTIKT, 0.5), "Jets");
      
      // Projections for MET
      addProjection(MissingMomentum(), "MET");
      
      // Weight counter
      _vis_unit_weights = 0.;

      // Booking of histograms      
      _h_met = bookHisto1D(5, 1, 1);
      _h_ht  = bookHisto1D(6, 1, 1);
      _h_st  = bookHisto1D(7, 1, 1);
      _h_wpt = bookHisto1D(8, 1, 1);
    }


    void analyze(const Event& event) {
      const double weight = event.weight();
      
      // select ttbar -> lepton+jets without taus
      const DressedLeptons& dressedelectrons = applyProjection<DressedLeptons>(event, "DressedElectrons");
      const DressedLeptons& dressedmuons = applyProjection<DressedLeptons>(event, "DressedMuons");
      if (dressedelectrons.dressedLeptons().size() + dressedmuons.dressedLeptons().size() != 1) vetoEvent;
      
      FourMomentum lepton;
      if (dressedelectrons.dressedLeptons().size() == 1) lepton = dressedelectrons.dressedLeptons()[0].momentum();
      else lepton = dressedmuons.dressedLeptons()[0].momentum();
      
      const TauFinder& taus = applyProjection<TauFinder>(event, "Tau");
      const IdentifiedFinalState nu_taus = applyProjection<IdentifiedFinalState>(event, "nu_tau");
      foreach(const Particle& tau, taus.taus()) {
        foreach(const Particle& nu, nu_taus.particles()) {
          if (tau.pid() * nu.pid() < 0) continue;
          const FourMomentum wCandidate = tau.momentum() + nu.momentum();
          if (abs(wCandidate.mass() - 80.4) > 5.) vetoEvent;
        }
      }

      // count weights in visible phase space
      if (weight != 0.) _vis_unit_weights += weight/std::abs(weight);
      
      // MET
      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MET");
      _h_met->fill(met.visibleMomentum().pT()/GeV, weight);
      
      // HT and ST
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
      const Jets jets = jetpro.jetsByPt(20*GeV);
      
      double ht = 0.0;
      foreach (const Jet& j, jets) {
        if (deltaR(j.momentum(), lepton) > 0.3) {
          ht += j.pT();
        }
      }
      double st = ht + lepton.pT() + met.visibleMomentum().pT();
      _h_ht->fill(ht/GeV, weight);
      _h_st->fill(st/GeV, weight);
      
      // WPT
      FourMomentum w = lepton - met.visibleMomentum();
      _h_wpt->fill(w.pT()/GeV, weight);
    }


    void finalize() {
      const double s = 1./_vis_unit_weights;
      scale(_h_met, s);
      scale(_h_ht, s);
      scale(_h_st, s);
      scale(_h_wpt, s);
    }

    //@}


  private:

    // @name Histogram data members
    //@{
    
    double _vis_unit_weights;
    
    Histo1DPtr _h_met, _h_ht, _h_st, _h_wpt;

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1473674);

}
