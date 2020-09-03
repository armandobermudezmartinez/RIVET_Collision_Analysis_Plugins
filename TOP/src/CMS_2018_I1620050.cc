#include "Rivet/Analysis.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include <vector>

namespace Rivet {
  namespace { //< only visible in this compilation unit

    /// @brief Special dressed lepton finder
    ///
    /// Find dressed leptons by clustering all leptons and photons
    class SpecialDressedLeptons : public FinalState {
      public:
        /// The default constructor. May specify cuts
        SpecialDressedLeptons(const FinalState& fs, const Cut& cut)
          : FinalState(cut)
        {
          setName("SpecialDressedLeptons");
          IdentifiedFinalState ifs(fs);
          ifs.acceptIdPair(PID::PHOTON);
          ifs.acceptIdPair(PID::ELECTRON);
          ifs.acceptIdPair(PID::MUON);
          addProjection(ifs, "IFS");
          addProjection(FastJets(ifs, FastJets::ANTIKT, 0.1), "LeptonJets");
        }

        /// Clone on the heap.
        virtual unique_ptr<Projection> clone() const {
          return unique_ptr<Projection>(new SpecialDressedLeptons(*this));
        }

        /// Retrieve the dressed leptons
        const vector<DressedLepton>& dressedLeptons() const { return _clusteredLeptons; }

      private:
        /// Container which stores the clustered lepton objects
        vector<DressedLepton> _clusteredLeptons;

      public:
        void project(const Event& e) {
          _theParticles.clear();
          _clusteredLeptons.clear();

          vector<DressedLepton> allClusteredLeptons;

          const Jets jets = applyProjection<FastJets>(e, "LeptonJets").jetsByPt(5.*GeV);
          foreach (const Jet& jet, jets) {
            Particle lepCand;
            for (const Particle& cand : jet.particles()) {
              const int absPdgId = abs(cand.pdgId());
              if (absPdgId == PID::ELECTRON || absPdgId == PID::MUON) {
                if (cand.pt() > lepCand.pt()) lepCand = cand;
              }
            }
            //Central lepton must be the major component -> NOTE: not used in this analysis
            if (lepCand.pdgId() == 0) continue;

            DressedLepton lepton = DressedLepton(lepCand);

            for (const Particle& cand : jet.particles()) {
              if (cand == lepCand) continue;
              lepton.addPhoton(cand, true);
            }
            allClusteredLeptons.push_back(lepton);
          }

          for (const DressedLepton& lepton : allClusteredLeptons) {
            if (accept(lepton)) {
              _clusteredLeptons.push_back(lepton);
              _theParticles.push_back(lepton.constituentLepton());
              _theParticles += lepton.constituentPhotons();
            }
          }
        }
    };
  }

  class CMS_2018_I1620050 : public Analysis {
    public:
      CMS_2018_I1620050() : Analysis("CMS_2018_I1620050") {}

      void init() {
        const bool acceptTauDecays = false;

        // Parton level top quark to analyze dilepton channels only
        // Note: 2nd argument of PartonicTops to toggle tau->lepton channel (true to inclusive, false to exclusive)
        declare(PartonicTops(PartonicTops::MUON, acceptTauDecays), "PartonTopsToMuon"); // Partonic top decaying to mu
        declare(PartonicTops(PartonicTops::ELECTRON, acceptTauDecays), "PartonTopsToElectron"); // Partonic top decaying to e

        // Build particle level tops starting from FinalState
        const FinalState fs(Cuts::pT > 0. && Cuts::abseta < 6.);

        // Neutrinos
        IdentifiedFinalState neutrinos(fs);
        neutrinos.acceptNeutrinos();
        PromptFinalState prompt_neutrinos(neutrinos, true, true);
        declare(prompt_neutrinos, "Neutrinos");

        // Projection for electrons and muons
        Cut leptonCuts = Cuts::pt > 20*GeV && Cuts::abseta < 2.4;

        PromptFinalState fsLepton(fs);
        fsLepton.acceptMuonDecays(true);
        fsLepton.acceptTauDecays(true);
        SpecialDressedLeptons dressedLeptons(fsLepton, leptonCuts);
        declare(dressedLeptons, "DressedLeptons");

        // Projection for jets
        VetoedFinalState fs_jets(FinalState(-MAXDOUBLE, MAXDOUBLE, 0*GeV));
        fs_jets.addVetoOnThisFinalState(dressedLeptons);
        fs_jets.vetoNeutrinos();
        declare(FastJets(fs_jets, FastJets::ANTIKT, 0.4), "ak4jets");

        //book hists
        _hist_lep_pt = bookHisto1D("d01-x01-y01");
        _hist_jet_pt = bookHisto1D("d02-x01-y01");
        _hist_top_pt = bookHisto1D("d03-x01-y01");
        _hist_top_y = bookHisto1D("d04-x01-y01");
        _hist_tt_pt = bookHisto1D("d05-x01-y01");
        _hist_tt_y = bookHisto1D("d06-x01-y01");
        _hist_tt_m = bookHisto1D("d07-x01-y01");
        _hist_tt_dphi = bookHisto1D("d08-x01-y01");
      }

      void analyze(const Event& event) {

        const double weight = event.weight();

        // Do the analysis only for the full-dleptonic channel
        const Particles partonTopsToMuon     = apply<ParticleFinder>(event, "PartonTopsToMuon").particles();
        //const Particles partonTopsToElectron = apply<ParticleFinder>(event, "PartonTopsToElectron").particles();
        Particles partonTopsToElectron;
        for ( auto x : apply<ParticleFinder>(event, "PartonTopsToElectron").particles() ) {
          bool isDuplicated = false;
          for ( auto y : partonTopsToMuon ) {
            if ( std::abs(x.pt()-y.pt()) < 0.01 and deltaR(x, y) < 0.01 ) {
              isDuplicated = true;
              break;
            }
          }
          if ( !isDuplicated ) partonTopsToElectron.push_back(x);
        }
        const int nPartonElectrons = partonTopsToElectron.size();
        const int nPartonMuons     = partonTopsToMuon.size();
        if ( nPartonElectrons+nPartonMuons != 2 ) vetoEvent;

        // Select leptons
        const std::vector<DressedLepton>& dressedLeptons = apply<SpecialDressedLeptons>(event, "DressedLeptons").dressedLeptons();
        if ( dressedLeptons.size() < 2 ) vetoEvent;
        sortByPt(dressedLeptons);

        const FourMomentum& lepton1 = dressedLeptons[0].momentum();
        const FourMomentum& lepton2 = dressedLeptons[1].momentum();
        const int channel = std::abs(dressedLeptons[0].pdgId())+std::abs(dressedLeptons[1].pdgId());
        if ( !((channel == 22 and nPartonElectrons == 2) or
               (channel == 24 and nPartonElectrons == 1 and nPartonMuons == 1) or
               (channel == 26 and nPartonMuons == 2)) ) vetoEvent;

        // Select neutrinos
        const Particles neutrinos = apply<PromptFinalState>(event, "Neutrinos").particlesByPt();
        if ( neutrinos.size() < 2 ) vetoEvent;

        // Select bjets
        Jets bJets;
        const FastJets& fjJets = apply<FastJets>(event, "ak4jets");
        const Jets jets = fjJets.jets(Cuts::abseta < 2.4 && Cuts::pT > 30*GeV);
        for ( Jets::const_iterator itjet = jets.begin(); itjet != jets.end() ; ++itjet) {
          if ( itjet->bTagged() ) { // Note: default b tagging algorithm is ghost association (see the Rivet Jet class reference manual)
            bJets.push_back(*itjet);
          }
        }
        // There should at least two b jets.
        if ( bJets.size() < 2 ) vetoEvent;
        sortByPt(bJets);

        // Construct particle level top
        FourMomentum nu1 = neutrinos[0].momentum();
        FourMomentum nu2 = neutrinos[1].momentum();
        if ( std::abs((lepton1+nu1).mass()-80.4) + std::abs((lepton2+nu2).mass()-80.4) >
             std::abs((lepton1+nu2).mass()-80.4) + std::abs((lepton2+nu1).mass()-80.4) ) {
          std::swap(nu1, nu2);
        }
        const FourMomentum w1 = lepton1 + nu1;
        const FourMomentum w2 = lepton2 + nu2;

        FourMomentum bjet1 = bJets[0].momentum();
        FourMomentum bjet2 = bJets[1].momentum();
        if ( std::abs((w1+bjet1).mass()-172.5) + std::abs((w2+bjet2).mass()-172.5) >
             std::abs((w1+bjet2).mass()-172.5) + std::abs((w2+bjet1).mass()-172.5) ) {
          std::swap(bjet1, bjet2);
        }
        const FourMomentum t1 = w1 + bjet1;
        const FourMomentum t2 = w2 + bjet2;
        const FourMomentum tt = t1+t2;

        _hist_lep_pt->fill(lepton1.pt(), weight);
        _hist_lep_pt->fill(lepton2.pt(), weight);
        _hist_jet_pt->fill(bjet1.pt(), weight);
        _hist_jet_pt->fill(bjet2.pt(), weight);

        _hist_top_pt->fill(t1.pt(), weight);
        _hist_top_pt->fill(t2.pt(), weight);
        _hist_top_y->fill(t1.rapidity(), weight);
        _hist_top_y->fill(t2.rapidity(), weight);

        _hist_tt_pt->fill(tt.pt(), weight);
        _hist_tt_y->fill(tt.rapidity(), weight);
        _hist_tt_m->fill(tt.mass(), weight);
        _hist_tt_dphi->fill(deltaPhi(t1.phi(), t2.phi()), weight);
      }

      /// Normalise histograms etc., after the run
      void finalize() {
        normalize(_hist_lep_pt);
        normalize(_hist_jet_pt);
        normalize(_hist_top_pt);
        normalize(_hist_top_y);
        normalize(_hist_tt_pt);
        normalize(_hist_tt_y);
        normalize(_hist_tt_m);
        normalize(_hist_tt_dphi);
      }

    private:

      Histo1DPtr _hist_lep_pt;
      Histo1DPtr _hist_jet_pt;
      Histo1DPtr _hist_top_pt;
      Histo1DPtr _hist_top_y;
      Histo1DPtr _hist_tt_pt;
      Histo1DPtr _hist_tt_y;
      Histo1DPtr _hist_tt_m;
      Histo1DPtr _hist_tt_dphi;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2018_I1620050);
}

