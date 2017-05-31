#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Tools/ParticleName.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {

  class CMS_2017_PAS_TOP_16_014 : public Analysis {
  public:

    // Minimal constructor
    CMS_2017_PAS_TOP_16_014() : Analysis("CMS_2017_PAS_TOP_16_014") {
    }

    // Set up projections and book histograms
    void init() {
      // Complete final state
      FinalState fs( (Cuts::abseta < 5) and (Cuts::pT > 0.0*MeV) );

      // // Projection for electrons and muons
      // IdentifiedFinalState photons(fs);
      // photons.acceptIdPair(PID::PHOTON);

      // IdentifiedFinalState el_id(fs);
      // el_id.acceptIdPair(PID::ELECTRON);
      // PromptFinalState electrons(el_id);
      // addProjection(electrons, "Electrons");
      // Cut electronLooseCuts = Cuts::pt > 15*GeV && Cuts::abseta < 2.4;
      // DressedLeptons dressed_electrons(photons, electrons, 0.1, electronLooseCuts, true, false);
      // addProjection(dressed_electrons, "DressedElectrons");

      // IdentifiedFinalState mu_id(fs);
      // mu_id.acceptIdPair(PID::MUON);
      // PromptFinalState muons(mu_id);
      // addProjection(muons, "Muons");
      // Cut muonLooseCuts = Cuts::pt > 15*GeV && Cuts::abseta < 2.4;
      // DressedLeptons dressed_muons(photons, muons, 0.1, Cuts::open(), true, false);
      // addProjection(dressed_muons, "DressedMuons");

      // Dressed leptons
      ChargedLeptons charged_leptons(fs);
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);
      
      PromptFinalState prompt_leptons(charged_leptons);
      prompt_leptons.acceptMuonDecays(true);
      prompt_leptons.acceptTauDecays(true);
      PromptFinalState prompt_photons(photons);
      prompt_photons.acceptMuonDecays(true);
      prompt_photons.acceptTauDecays(true);
      Cut looseLeptonCuts = Cuts::pt > 15*GeV && Cuts::abseta < 2.4;

      DressedLeptons dressed_leptons(prompt_photons, prompt_leptons, 0.1, 
                     looseLeptonCuts, true, true);
      addProjection(dressed_leptons, "DressedLeptons");

      // Projection for jets
      VetoedFinalState fsForJets(fs);
      fsForJets.addVetoOnThisFinalState(dressed_leptons);
      addProjection(FastJets(fsForJets, FastJets::ANTIKT, 0.4), "Jets");

      // Projections for MET
      addProjection(MissingMomentum(fs), "MET");

      // Booking of histograms
      _hist_norm_met = bookHisto1D(1, 1, 1);
      _hist_norm_ht  = bookHisto1D(2, 1, 1);
      _hist_norm_st  = bookHisto1D(3, 1, 1);
      _hist_norm_wpt = bookHisto1D(4, 1, 1);
      _hist_norm_njets = bookHisto1D(5, 1, 1);
      _hist_norm_lpt = bookHisto1D(6, 1, 1);
      _hist_norm_labseta = bookHisto1D(7, 1, 1);

      _hist_abs_met = bookHisto1D(8, 1, 1);
      _hist_abs_ht  = bookHisto1D(9, 1, 1);
      _hist_abs_st  = bookHisto1D(10, 1, 1);
      _hist_abs_wpt = bookHisto1D(11, 1, 1);
      _hist_abs_njets = bookHisto1D(12, 1, 1);
      _hist_abs_lpt = bookHisto1D(13, 1, 1);
      _hist_abs_labseta = bookHisto1D(14, 1, 1);

    }


    // per event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // select ttbar -> lepton+jets at particle level
      const DressedLeptons& dressed_leptons = applyProjection<DressedLeptons>(event, "DressedLeptons");
      if ( dressed_leptons.dressedLeptons().size() != 1 ) {
        vetoEvent;
      }

      // Lepton selection
      FourMomentum lepton = dressed_leptons.dressedLeptons()[0];

      const double leptonPt = lepton.pT();
      const double leptonAbsEta = std::abs( lepton.eta() );

      if ( leptonPt < 26 or leptonAbsEta > 2.4 ) vetoEvent;


      // Jet selection
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
      const Jets jets = jetpro.jets(Cuts::abseta < 2.4 && Cuts::pT > 30*GeV);
      Jets cleanedJets;
      unsigned int nBJets = 0;
      foreach (const Jet& j, jets) {
        // if (deltaR(j.momentum(), lepton) > 0.4) {
        // std::cout << "Jet pt, eta : " << j.pT() << " " << j.eta() << std::endl;
          cleanedJets.push_back( j );
          if ( j.bTagged() ) ++nBJets;
        // }
      }

      if ( cleanedJets.size() < 4 ) vetoEvent;
      if ( nBJets < 2 ) vetoEvent;

      // MET
      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MET");
      _hist_norm_met->fill(met.visibleMomentum().pT()/GeV, weight);
      _hist_abs_met->fill(met.visibleMomentum().pT()/GeV, weight);

      // HT and ST
      double ht = 0.0;
      foreach (const Jet& j, cleanedJets) {
          ht += j.pT();
      }

      double st = ht + lepton.pT() + met.visibleMomentum().pT();
      _hist_norm_ht->fill(ht/GeV, weight);
      _hist_norm_st->fill(st/GeV, weight);

      _hist_abs_ht->fill(ht/GeV, weight);
      _hist_abs_st->fill(st/GeV, weight);

      // WPT
      FourMomentum w = lepton - met.visibleMomentum();
      _hist_norm_wpt->fill(w.pT()/GeV, weight);
      _hist_abs_wpt->fill(w.pT()/GeV, weight);

      // Lepton pt and eta
      _hist_norm_lpt->fill( leptonPt/GeV, weight);
      _hist_norm_labseta->fill( leptonAbsEta/GeV, weight);

      _hist_abs_lpt->fill( leptonPt/GeV, weight);
      _hist_abs_labseta->fill( leptonAbsEta/GeV, weight);

      // NJets
      _hist_norm_njets->fill( cleanedJets.size(), weight );
      _hist_abs_njets->fill( cleanedJets.size(), weight );

    }

    // scale by 1 over weight
    void finalize() {
      normalize(_hist_norm_met);
      normalize(_hist_norm_ht);
      normalize(_hist_norm_st);
      normalize(_hist_norm_wpt);
      normalize(_hist_norm_njets);
      normalize(_hist_norm_lpt);
      normalize(_hist_norm_labseta);

      scale(_hist_abs_met, crossSection() / sumOfWeights() );
      scale(_hist_abs_ht, crossSection() / sumOfWeights() );
      scale(_hist_abs_st, crossSection() / sumOfWeights() );
      scale(_hist_abs_wpt, crossSection() / sumOfWeights() );
      scale(_hist_abs_njets, crossSection() / sumOfWeights() );
      scale(_hist_abs_lpt, crossSection() / sumOfWeights() );
      scale(_hist_abs_labseta, crossSection() / sumOfWeights() );

    }

  private:
    Histo1DPtr _hist_norm_met, _hist_norm_ht, _hist_norm_st, _hist_norm_wpt, _hist_norm_njets, _hist_norm_lpt, _hist_norm_labseta;
    Histo1DPtr _hist_abs_met, _hist_abs_ht, _hist_abs_st, _hist_abs_wpt, _hist_abs_njets, _hist_abs_lpt, _hist_abs_labseta;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_PAS_TOP_16_014);
}
