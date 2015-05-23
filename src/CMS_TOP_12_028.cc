#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
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

    _h00_diffXSecTopSemiLepHadronPhaseSpacelepPt    = bookHisto1D("h00_diffXSecTopSemiLepHadronPhaseSpacelepPt"    );
    _h01_diffXSecTopSemiLepHadronPhaseSpacelepEta   = bookHisto1D("h01_diffXSecTopSemiLepHadronPhaseSpacelepEta"   );
    _h02_diffXSecTopSemiLepHadronPhaseSpacebqPt     = bookHisto1D("h02_diffXSecTopSemiLepHadronPhaseSpacebqPt"     );
    _h03_diffXSecTopSemiLepHadronPhaseSpacebqEta    = bookHisto1D("h03_diffXSecTopSemiLepHadronPhaseSpacebqEta"    );
    _h04_diffXSecTopSemiLepHadronPhaseSpacebbbarPt  = bookHisto1D("h04_diffXSecTopSemiLepHadronPhaseSpacebbbarPt"  );
    _h05_diffXSecTopSemiLepHadronPhaseSpacebbbarMass= bookHisto1D("h05_diffXSecTopSemiLepHadronPhaseSpacebbbarMass");
                                                                          
    _h06_diffXSecTopDiLepHadronPhaseSpacelepPt      = bookHisto1D("h06_diffXSecTopDiLepHadronPhaseSpacelepPt"      );
    _h07_diffXSecTopDiLepHadronPhaseSpacelepEta     = bookHisto1D("h07_diffXSecTopDiLepHadronPhaseSpacelepEta"     );
    _h08_diffXSecTopDiLepHadronPhaseSpacedilepPt    = bookHisto1D("h08_diffXSecTopDiLepHadronPhaseSpacedilepPt"    );
    _h09_diffXSecTopDiLepHadronPhaseSpacedilepMass  = bookHisto1D("h09_diffXSecTopDiLepHadronPhaseSpacedilepMass"  );
    _h10_diffXSecTopDiLepHadronPhaseSpacebqPt       = bookHisto1D("h10_diffXSecTopDiLepHadronPhaseSpacebqPt"       );
    _h11_diffXSecTopDiLepHadronPhaseSpacebqEta      = bookHisto1D("h11_diffXSecTopDiLepHadronPhaseSpacebqEta"      );
    _h12_diffXSecTopDiLepHadronPhaseSpacebbbarPt    = bookHisto1D("h12_diffXSecTopDiLepHadronPhaseSpacebbbarPt"    );
    _h13_diffXSecTopDiLepHadronPhaseSpacebbbarMass  = bookHisto1D("h13_diffXSecTopDiLepHadronPhaseSpacebbbarMass"  );
  };

  void findBAncestors(const ParticleVector& ps, std::set<int>& bAncestors) {
    bAncestors.clear();

    foreach (const Particle& p, ps) {
      GenVertex* v = p.genParticle()->production_vertex();
      if ( !v ) continue;

      foreach (const GenParticle* ancestor, Rivet::particles(v, HepMC::ancestors)) {
        if ( ancestor->status() != 2 ) continue;
        const PdgId pid = ancestor->pdg_id();
        if ( !PID::isHadron(pid) or !PID::hasBottom(pid) ) continue;

        GenVertex* av = ancestor->production_vertex();
        if ( !av ) continue;

        bool isDuplicated = false;
        foreach (const GenParticle* ap, Rivet::particles(av, HepMC::parents)) {
          if ( p.genParticle() != ap && pid == ap->pdg_id() ) {
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

    // Find leptons
    Particles lCands;
    if ( ttbarState.mode() == PartonTTbarState::CH_SEMILEPTON ) {
      lCands.push_back(Particle());
      foreach (const Particle& p, ttbarState.wDecays1()) {
        const int absId = std::abs(p.pdgId());
        if ( absId == 11 or absId == 13 ) { lCands[0] = p; break; }
      }
      foreach (const Particle& p, ttbarState.wDecays2()) {
        const int absId = std::abs(p.pdgId());
        if ( absId == 11 or absId == 13 ) { lCands[0] = p; break; }
      }
      // Apply the particle level phase space cut
      if ( lCands[0].pT() <= 33 or std::abs(lCands[0].eta()) >= 2.1 ) vetoEvent;
    }
    else if ( ttbarState.mode() == PartonTTbarState::CH_FULLLEPTON ) {
      lCands.push_back(Particle());
      foreach (const Particle& p, ttbarState.wDecays1()) {
        const int absId = std::abs(p.pdgId());
        if ( absId == 11 or absId == 13 ) { lCands[0] = p; break; }
      }
      lCands.push_back(Particle());
      foreach (const Particle& p, ttbarState.wDecays2()) {
        const int absId = std::abs(p.pdgId());
        if ( absId == 11 or absId == 13 ) { lCands[1] = p; break; }
      }
      if ( lCands[0].pT() < lCands[1].pT() ) std::swap(lCands[0], lCands[1]);
      const double l1Pt = lCands[0].pT(), l1Abseta = std::abs(lCands[0].eta());
      const double l2Pt = lCands[1].pT(), l2Abseta = std::abs(lCands[1].eta());

      // Apply the particle level phase space cut
      if ( l1Pt <= 20 or l1Abseta >= 2.4 or l2Pt <= 20 or l2Abseta >= 2.4 ) vetoEvent;
    }

    // Build genJets
    const Jets& jetsIn = applyProjection<JetAlg>(event, "Jets").jetsByPt(30*GeV, MAXDOUBLE, -2.4, 2.4);
    Jets jets;
    foreach ( const Jet& jet, jetsIn ) {
      foreach ( const Particle& lCand, lCands ) {
        if ( deltaR(lCand.momentum(), jet.momentum()) < 0.3 ) continue;
        jets.push_back(jet);
      }
    }
    if ( ttbarState.mode() != PartonTTbarState::CH_SEMILEPTON and jets.size() < 4 ) vetoEvent;
    else if ( ttbarState.mode() != PartonTTbarState::CH_FULLLEPTON and jets.size() < 2 ) vetoEvent;

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
    // Require at least 2 b jets.
    // It is safe with checking bjet2 only. But let's check bjet1 for clearity.
    if ( !bjet2 or !bjet1 ) vetoEvent;

    const FourMomentum& b1P4 = bjet1->momentum();
    const FourMomentum& b2P4 = bjet2->momentum();
    const FourMomentum bbP4 = b1P4+b2P4;
    const double b1Pt = b1P4.pT(), b1Eta = b1P4.eta();
    const double b2Pt = b2P4.pT(), b2Eta = b2P4.eta();
    const double bbPt = bbP4.pT(), bbMass = bbP4.mass();

    // Find leptons, apply channel dependent phase space cuts, fill histograms
    if ( ttbarState.mode() == PartonTTbarState::CH_SEMILEPTON ) {
      // Do the semileptonic channel
      const FourMomentum& lP4 = lCands[0].momentum();
      const double lPt = lP4.pT(), lEta = lP4.eta();

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
      const FourMomentum& l1P4 = lCands[0].momentum();
      const FourMomentum& l2P4 = lCands[1].momentum();
      const FourMomentum dilP4 = l1P4+l2P4;
      const double dilMass = dilP4.mass();
      if ( dilMass < 20 ) vetoEvent;
      const double l1Pt = l1P4.pT(), l1Eta = l1P4.eta();
      const double l2Pt = l2P4.pT(), l2Eta = l2P4.eta();

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
  Histo1DPtr _h00_diffXSecTopSemiLepHadronPhaseSpacelepPt;
  Histo1DPtr _h01_diffXSecTopSemiLepHadronPhaseSpacelepEta;
  Histo1DPtr _h02_diffXSecTopSemiLepHadronPhaseSpacebqPt;
  Histo1DPtr _h03_diffXSecTopSemiLepHadronPhaseSpacebqEta;
  Histo1DPtr _h04_diffXSecTopSemiLepHadronPhaseSpacebbbarPt;
  Histo1DPtr _h05_diffXSecTopSemiLepHadronPhaseSpacebbbarMass;
                          
  Histo1DPtr _h06_diffXSecTopDiLepHadronPhaseSpacelepPt;
  Histo1DPtr _h07_diffXSecTopDiLepHadronPhaseSpacelepEta;
  Histo1DPtr _h08_diffXSecTopDiLepHadronPhaseSpacedilepPt;
  Histo1DPtr _h09_diffXSecTopDiLepHadronPhaseSpacedilepMass;
  Histo1DPtr _h10_diffXSecTopDiLepHadronPhaseSpacebqPt;
  Histo1DPtr _h11_diffXSecTopDiLepHadronPhaseSpacebqEta;
  Histo1DPtr _h12_diffXSecTopDiLepHadronPhaseSpacebbbarPt;
  Histo1DPtr _h13_diffXSecTopDiLepHadronPhaseSpacebbbarMass;
};

// This global object acts as a hook for the plugin system
AnalysisBuilder<CMS_TOP_12_028> plugin_CMS_TOP_12_028;

}
