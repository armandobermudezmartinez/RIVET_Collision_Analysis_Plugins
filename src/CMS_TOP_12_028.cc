#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "GeneratorInterface/RivetTop/interface/PartonTTbarState.hh"

namespace Rivet {

class CMS_TOP_12_028 : public Analysis {
public:
  CMS_TOP_12_028() : Analysis("CMS_TOP_12_028") {
  }

  void init() {
    // Parton level top quarks
    PartonTTbarState ttbarState;
    addProjection(ttbarState, "ttbar");

    FinalState fs(-5.0, 5.0, 0*GeV);
    VetoedFinalState fsForJets(fs);
    fsForJets.addDecayProductsVeto(+24);
    fsForJets.addDecayProductsVeto(-24);

    FastJets fj(fsForJets, FastJets::ANTIKT, 0.5);
    fj.useInvisibles();
    addProjection(fj, "Jets");

    _h00_diffXSecTopSemiLepHadronPhaseSpacelepPt    = bookHistogram1D("h00_diffXSecTopSemiLepHadronPhaseSpacelepPt"    );
    _h01_diffXSecTopSemiLepHadronPhaseSpacelepEta   = bookHistogram1D("h01_diffXSecTopSemiLepHadronPhaseSpacelepEta"   );
    _h02_diffXSecTopSemiLepHadronPhaseSpacebqPt     = bookHistogram1D("h02_diffXSecTopSemiLepHadronPhaseSpacebqPt"     );
    _h03_diffXSecTopSemiLepHadronPhaseSpacebqEta    = bookHistogram1D("h03_diffXSecTopSemiLepHadronPhaseSpacebqEta"    );
    _h04_diffXSecTopSemiLepHadronPhaseSpacebbbarPt  = bookHistogram1D("h04_diffXSecTopSemiLepHadronPhaseSpacebbbarPt"  );
    _h05_diffXSecTopSemiLepHadronPhaseSpacebbbarMass= bookHistogram1D("h05_diffXSecTopSemiLepHadronPhaseSpacebbbarMass");
                                                                          
    _h06_diffXSecTopDiLepHadronPhaseSpacelepPt      = bookHistogram1D("h06_diffXSecTopDiLepHadronPhaseSpacelepPt"      );
    _h07_diffXSecTopDiLepHadronPhaseSpacelepEta     = bookHistogram1D("h07_diffXSecTopDiLepHadronPhaseSpacelepEta"     );
    _h08_diffXSecTopDiLepHadronPhaseSpacedilepPt    = bookHistogram1D("h08_diffXSecTopDiLepHadronPhaseSpacedilepPt"    );
    _h09_diffXSecTopDiLepHadronPhaseSpacedilepMass  = bookHistogram1D("h09_diffXSecTopDiLepHadronPhaseSpacedilepMass"  );
    _h10_diffXSecTopDiLepHadronPhaseSpacebqPt       = bookHistogram1D("h10_diffXSecTopDiLepHadronPhaseSpacebqPt"       );
    _h11_diffXSecTopDiLepHadronPhaseSpacebqEta      = bookHistogram1D("h11_diffXSecTopDiLepHadronPhaseSpacebqEta"      );
    _h12_diffXSecTopDiLepHadronPhaseSpacebbbarPt    = bookHistogram1D("h12_diffXSecTopDiLepHadronPhaseSpacebbbarPt"    );
    _h13_diffXSecTopDiLepHadronPhaseSpacebbbarMass  = bookHistogram1D("h13_diffXSecTopDiLepHadronPhaseSpacebbbarMass"  );
  };

  void findBAncestors(const ParticleVector& ps, std::set<int>& bAncestors) {
    bAncestors.clear();

    foreach (const Particle& p, ps) {
      GenVertex* v = p.genParticle().production_vertex();
      if ( !v ) continue;

      foreach (const GenParticle* ancestor, Rivet::particles(v, HepMC::ancestors)) {
        if ( ancestor->status() != 2 ) continue;
        const PdgId pid = ancestor->pdg_id();
        if ( !PID::isHadron(pid) or !PID::hasBottom(pid) ) continue;

        GenVertex* av = ancestor->production_vertex();
        if ( !av ) continue;

        bool isDuplicated = false;
        foreach (const GenParticle* ap, Rivet::particles(av, HepMC::parents)) {
          if ( &(p.genParticle()) != ap && pid == ap->pdg_id() ) {
            isDuplicated = true;
            break;
          }
        }

        if ( !isDuplicated ) {
          bAncestors.insert(ancestor->barcode());
        }
      }
    }
  }
  void analyze(const Event& event) {
    const double weight = event.weight();

    // Get the parton level ttbar candidate
    const PartonTTbarState& ttbarState = applyProjection<PartonTTbarState>(event, "ttbar");

    // Do the analysis only for the ttbar full leptonic or semileptonic channel, without tau decay
    if ( ttbarState.mode() != PartonTTbarState::CH_SEMILEPTON and
         ttbarState.mode() != PartonTTbarState::CH_FULLLEPTON ) vetoEvent;
    if ( ttbarState.mode1() >= PartonTTbarState::CH_TAU_HADRON ||
         ttbarState.mode2() >= PartonTTbarState::CH_TAU_HADRON ) vetoEvent;

    // Build genJets
    const Jets& jets = applyProjection<JetAlg>(event, "Jets").jetsByPt(30*GeV, MAXDOUBLE, -2.4, 2.4);
    if ( jets.size() < 4 ) vetoEvent;

    // Get Leading two jets
    std::set<int> bAncestors;
    const Jet* bjet1 = 0, * bjet2 = 0;
    foreach (const Jet& jet, jets) {
      if ( !bjet1 ) {
        findBAncestors(jet.particles(), bAncestors);
        if ( !bAncestors.empty() ) bjet1 = &jet;
      } else if ( !bjet2 ) {
        std::set<int> bAncestors2;
        findBAncestors(jet.particles(), bAncestors2);
        bool hasUniqueB = false;
        foreach (int barCode, bAncestors2) {
          if ( bAncestors.find(barCode) == bAncestors.end() ) {
            hasUniqueB = true;
            break;
          }
        }
        if ( hasUniqueB ) {
          bjet2 = &jet;
          break;
        }
      } else break;
    }
    if ( !bjet1 or !bjet2 ) vetoEvent;

    const FourMomentum& b1P4 = bjet1->momentum();
    const FourMomentum& b2P4 = bjet2->momentum();
    const FourMomentum bbP4 = b1P4+b2P4;
    const double b1Pt = b1P4.pT(), b1Eta = b1P4.eta();
    const double b2Pt = b2P4.pT(), b2Eta = b2P4.eta();
    const double bbPt = bbP4.pT(), bbMass = bbP4.mass();

    // Find leptons, apply channel dependent phase space cuts, fill histograms
    if ( ttbarState.mode() == PartonTTbarState::CH_SEMILEPTON ) {
      // Do the semileptonic channel
      Particle lCand;
      foreach (const Particle& p, ttbarState.wDecays1()) {
        if ( PID::isLepton(p.pdgId()) ) { lCand = p; break; }
      }
      foreach (const Particle& p, ttbarState.wDecays2()) {
        if ( PID::isLepton(p.pdgId()) ) { lCand = p; break; }
      }

      // Build lepton
      const FourMomentum& lP4 = lCand.momentum();
      const double lPt = lP4.pT(), lEta = lP4.eta();

      // Leptons and b jet observables are considered only within the particle level phase space
      if ( lPt <= 33 or std::abs(lEta) >= 2.1 ) vetoEvent;
      if ( b1Pt <= 30 or std::abs(b1Eta) >= 2.4 or b2Pt <= 30 or std::abs(b2Eta) >= 2.4 ) vetoEvent;

      _h00_diffXSecTopSemiLepHadronPhaseSpacelepPt->fill(lPt, weight);
      _h01_diffXSecTopSemiLepHadronPhaseSpacelepEta->fill(lEta, weight);
      _h02_diffXSecTopSemiLepHadronPhaseSpacebqPt->fill(b1Pt, weight);
      _h02_diffXSecTopSemiLepHadronPhaseSpacebqPt->fill(b2Pt, weight);
      _h03_diffXSecTopSemiLepHadronPhaseSpacebqEta->fill(b1Eta, weight);
      _h03_diffXSecTopSemiLepHadronPhaseSpacebqEta->fill(b2Eta, weight);
      _h04_diffXSecTopSemiLepHadronPhaseSpacebbbarPt->fill(bbPt, weight);
      _h05_diffXSecTopSemiLepHadronPhaseSpacebbbarMass->fill(bbMass, weight);
    }
    else if ( ttbarState.mode() == PartonTTbarState::CH_FULLLEPTON ) {
      // Do the dileptonic channel
      Particle l1Cand, l2Cand;
      foreach (const Particle& p, ttbarState.wDecays1()) {
        if ( PID::isLepton(p.pdgId()) ) { l1Cand = p; break; }
      }
      foreach (const Particle& p, ttbarState.wDecays2()) {
        if ( PID::isLepton(p.pdgId()) ) { l2Cand = p; break; }
      }

      // Build lepton
      const FourMomentum& l1P4 = l1Cand.momentum();
      const FourMomentum& l2P4 = l2Cand.momentum();
      const FourMomentum dilP4 = l1P4+l2P4;
      const double l1Pt = l1P4.pT(), l1Eta = l1P4.eta();
      const double l2Pt = l2P4.pT(), l2Eta = l2P4.eta();

      if ( l1Pt <= 20 or std::abs(l1Eta) >= 2.4 or l2Pt <= 20 or std::abs(l2Eta) >= 2.4 ) vetoEvent;
      if ( b1Pt <= 30 or std::abs(b1Eta) >= 2.4 or b2Pt <= 30 or std::abs(b2Eta) >= 2.4 ) vetoEvent;

      _h06_diffXSecTopDiLepHadronPhaseSpacelepPt->fill(l1Pt, weight);
      _h06_diffXSecTopDiLepHadronPhaseSpacelepPt->fill(l2Pt, weight);
      _h07_diffXSecTopDiLepHadronPhaseSpacelepEta->fill(l1Eta, weight);
      _h07_diffXSecTopDiLepHadronPhaseSpacelepEta->fill(l2Eta, weight);
      _h08_diffXSecTopDiLepHadronPhaseSpacedilepPt->fill(dilP4.pT(), weight);
      _h09_diffXSecTopDiLepHadronPhaseSpacedilepMass->fill(dilP4.mass(), weight);

      _h10_diffXSecTopDiLepHadronPhaseSpacebqPt->fill(b1Pt, weight);
      _h10_diffXSecTopDiLepHadronPhaseSpacebqPt->fill(b2Pt, weight);
      _h11_diffXSecTopDiLepHadronPhaseSpacebqEta->fill(b1Eta, weight);
      _h11_diffXSecTopDiLepHadronPhaseSpacebqEta->fill(b2Eta, weight);
      _h12_diffXSecTopDiLepHadronPhaseSpacebbbarPt->fill(bbPt, weight);
      _h13_diffXSecTopDiLepHadronPhaseSpacebbbarMass->fill(bbMass, weight);
    }
  };

  void finalize() {
    normalize(_h00_diffXSecTopSemiLepHadronPhaseSpacelepPt    );
    normalize(_h01_diffXSecTopSemiLepHadronPhaseSpacelepEta   );
    normalize(_h02_diffXSecTopSemiLepHadronPhaseSpacebqPt     );
    normalize(_h03_diffXSecTopSemiLepHadronPhaseSpacebqEta    );
    normalize(_h04_diffXSecTopSemiLepHadronPhaseSpacebbbarPt  );
    normalize(_h05_diffXSecTopSemiLepHadronPhaseSpacebbbarMass);
                  
    normalize(_h06_diffXSecTopDiLepHadronPhaseSpacelepPt      );
    normalize(_h07_diffXSecTopDiLepHadronPhaseSpacelepEta     );
    normalize(_h08_diffXSecTopDiLepHadronPhaseSpacedilepPt    );
    normalize(_h09_diffXSecTopDiLepHadronPhaseSpacedilepMass  );
    normalize(_h10_diffXSecTopDiLepHadronPhaseSpacebqPt       );
    normalize(_h11_diffXSecTopDiLepHadronPhaseSpacebqEta      );
    normalize(_h12_diffXSecTopDiLepHadronPhaseSpacebbbarPt    );
    normalize(_h13_diffXSecTopDiLepHadronPhaseSpacebbbarMass  );
  };

private:
  AIDA::IHistogram1D* _h00_diffXSecTopSemiLepHadronPhaseSpacelepPt;
  AIDA::IHistogram1D* _h01_diffXSecTopSemiLepHadronPhaseSpacelepEta;
  AIDA::IHistogram1D* _h02_diffXSecTopSemiLepHadronPhaseSpacebqPt;
  AIDA::IHistogram1D* _h03_diffXSecTopSemiLepHadronPhaseSpacebqEta;
  AIDA::IHistogram1D* _h04_diffXSecTopSemiLepHadronPhaseSpacebbbarPt;
  AIDA::IHistogram1D* _h05_diffXSecTopSemiLepHadronPhaseSpacebbbarMass;
                          
  AIDA::IHistogram1D* _h06_diffXSecTopDiLepHadronPhaseSpacelepPt;
  AIDA::IHistogram1D* _h07_diffXSecTopDiLepHadronPhaseSpacelepEta;
  AIDA::IHistogram1D* _h08_diffXSecTopDiLepHadronPhaseSpacedilepPt;
  AIDA::IHistogram1D* _h09_diffXSecTopDiLepHadronPhaseSpacedilepMass;
  AIDA::IHistogram1D* _h10_diffXSecTopDiLepHadronPhaseSpacebqPt;
  AIDA::IHistogram1D* _h11_diffXSecTopDiLepHadronPhaseSpacebqEta;
  AIDA::IHistogram1D* _h12_diffXSecTopDiLepHadronPhaseSpacebbbarPt;
  AIDA::IHistogram1D* _h13_diffXSecTopDiLepHadronPhaseSpacebbbarMass;
};

// This global object acts as a hook for the plugin system
AnalysisBuilder<CMS_TOP_12_028> plugin_CMS_TOP_12_028;

}
