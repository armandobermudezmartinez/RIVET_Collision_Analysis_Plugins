#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "GeneratorInterface/RivetTop/interface/CMSGenParticle.hh"
#include "GeneratorInterface/RivetTop/interface/PartonTop.hh"

namespace Rivet {

class CMS_TOP_12_028 : public Analysis {
public:
  CMS_TOP_12_028() : Analysis("CMS_TOP_12_028") {
  }

  void init() {
    // Parton level top quarks
    PartonTop ttbarState;
    addProjection(ttbarState, "ttbar");

    //FinalState fs(-5.0, 5.0, 0*GeV);
    //VetoedFinalState fsForJets(fs);
    //fsForJets.addDecayProductsVeto(+24);
    //fsForJets.addDecayProductsVeto(-24);

    CMSGenParticle fsForJets;
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

    _h14_diffXSecTopSemiLepPartontopPt         = bookHisto1D("h14_diffXSecTopSemiLepPartontopPt"        );
    _h15_diffXSecTopSemiLepPartontopPtTtbarSys = bookHisto1D("h15_diffXSecTopSemiLepPartontopPtTtbarSys");
    _h16_diffXSecTopSemiLepPartontopY          = bookHisto1D("h16_diffXSecTopSemiLepPartontopY"         );
    _h17_diffXSecTopSemiLepPartonttbarDelPhi   = bookHisto1D("h17_diffXSecTopSemiLepPartonttbarDelPhi"  );
    _h18_diffXSecTopSemiLepPartontopPtLead     = bookHisto1D("h18_diffXSecTopSemiLepPartontopPtLead"    );
    _h19_diffXSecTopSemiLepPartontopPtSubLead  = bookHisto1D("h19_diffXSecTopSemiLepPartontopPtSubLead" );
    _h20_diffXSecTopSemiLepPartonttbarPt       = bookHisto1D("h20_diffXSecTopSemiLepPartonttbarPt"      );
    _h21_diffXSecTopSemiLepPartonttbarY        = bookHisto1D("h21_diffXSecTopSemiLepPartonttbarY"       );
    _h22_diffXSecTopSemiLepPartonttbarMass     = bookHisto1D("h22_diffXSecTopSemiLepPartonttbarMass"    );

    _h23_diffXSecTopDiLepPartontopPt           = bookHisto1D("h23_diffXSecTopDiLepPartontopPt"          );
    _h24_diffXSecTopDiLepPartontopPtTtbarSys   = bookHisto1D("h24_diffXSecTopDiLepPartontopPtTtbarSys"  );
    _h25_diffXSecTopDiLepPartontopY            = bookHisto1D("h25_diffXSecTopDiLepPartontopY"           );
    _h26_diffXSecTopDiLepPartonttbarDelPhi     = bookHisto1D("h26_diffXSecTopDiLepPartonttbarDelPhi"    );
    _h27_diffXSecTopDiLepPartontopPtLead       = bookHisto1D("h27_diffXSecTopDiLepPartontopPtLead"      );
    _h28_diffXSecTopDiLepPartontopPtSubLead    = bookHisto1D("h28_diffXSecTopDiLepPartontopPtSubLead"   );
    _h29_diffXSecTopDiLepPartonttbarPt         = bookHisto1D("h29_diffXSecTopDiLepPartonttbarPt"        );
    _h30_diffXSecTopDiLepPartonttbarY          = bookHisto1D("h30_diffXSecTopDiLepPartonttbarY"         );
    _h31_diffXSecTopDiLepPartonttbarMass       = bookHisto1D("h31_diffXSecTopDiLepPartonttbarMass"      );

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
    const PartonTop& ttbarState = applyProjection<PartonTop>(event, "ttbar");
    // Do the analysis only for the ttbar full leptonic or semileptonic channel, without tau decay
    if ( ttbarState.mode() != PartonTop::CH_SEMILEPTON and
         ttbarState.mode() != PartonTop::CH_FULLLEPTON ) vetoEvent;
    if ( ttbarState.mode1() >= PartonTop::CH_TAU_HADRON ||
         ttbarState.mode2() >= PartonTop::CH_TAU_HADRON ) vetoEvent;

    // Parton level at full phase space
    // Fill top quarks they are defined in the parton level, full phase space
    const FourMomentum& t1P4 = ttbarState.t1().momentum();
    const FourMomentum& t2P4 = ttbarState.t2().momentum();
    const double t1Pt = t1P4.pT(), t2Pt = t2P4.pT();
    const FourMomentum ttbarP4 = t1P4+t2P4;
    const FourMomentum t1P4AtCM = LorentzTransform(-ttbarP4.boostVector()).transform(t1P4);
    const double dPhi = deltaPhi(t1P4.phi(), t2P4.phi());

    if ( ttbarState.mode() == PartonTop::CH_SEMILEPTON ) {
      _h14_diffXSecTopSemiLepPartontopPt->fill(t1Pt, weight);
      _h14_diffXSecTopSemiLepPartontopPt->fill(t2Pt, weight);
      _h15_diffXSecTopSemiLepPartontopPtTtbarSys->fill(t1P4AtCM.pT(), weight);
      _h16_diffXSecTopSemiLepPartontopY->fill(t1P4.rapidity(), weight);
      _h16_diffXSecTopSemiLepPartontopY->fill(t2P4.rapidity(), weight);
      _h17_diffXSecTopSemiLepPartonttbarDelPhi->fill(dPhi, weight);
      _h18_diffXSecTopSemiLepPartontopPtLead->fill(std::max(t1Pt, t2Pt), weight);
      _h19_diffXSecTopSemiLepPartontopPtSubLead->fill(std::min(t1Pt, t2Pt), weight);
      _h20_diffXSecTopSemiLepPartonttbarPt->fill(ttbarP4.pT(), weight);
      _h21_diffXSecTopSemiLepPartonttbarY->fill(ttbarP4.rapidity(), weight);
      _h22_diffXSecTopSemiLepPartonttbarMass->fill(ttbarP4.mass(), weight);
    }
    else if ( ttbarState.mode() == PartonTop::CH_FULLLEPTON ) {
      _h23_diffXSecTopDiLepPartontopPt->fill(t1Pt, weight);
      _h23_diffXSecTopDiLepPartontopPt->fill(t2Pt, weight);
      _h24_diffXSecTopDiLepPartontopPtTtbarSys->fill(t1P4AtCM.pT(), weight);
      _h25_diffXSecTopDiLepPartontopY->fill(t1P4.rapidity(), weight);
      _h25_diffXSecTopDiLepPartontopY->fill(t2P4.rapidity(), weight);
      _h26_diffXSecTopDiLepPartonttbarDelPhi->fill(dPhi, weight);
      _h27_diffXSecTopDiLepPartontopPtLead->fill(std::max(t1Pt, t2Pt), weight);
      _h28_diffXSecTopDiLepPartontopPtSubLead->fill(std::min(t1Pt, t2Pt), weight);
      _h29_diffXSecTopDiLepPartonttbarPt->fill(ttbarP4.pT(), weight);
      _h30_diffXSecTopDiLepPartonttbarY->fill(ttbarP4.rapidity(), weight);
      _h31_diffXSecTopDiLepPartonttbarMass->fill(ttbarP4.mass(), weight);
    }

    // Find leptons
    Particles lCands;
    if ( ttbarState.mode() == PartonTop::CH_SEMILEPTON ) {
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
    else if ( ttbarState.mode() == PartonTop::CH_FULLLEPTON ) {
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
    const Jets& jetsIn = applyProjection<JetAlg>(event, "Jets").jetsByPt(30*GeV);
    Jets jets;
    foreach ( const Jet& jet, jetsIn ) {
      if ( std::abs(jet.eta()) > 2.4 ) continue;
      foreach ( const Particle& lCand, lCands ) {
        if ( deltaR(lCand.momentum(), jet.momentum()) < 0.3 ) continue;
        jets.push_back(jet);
      }
    }
    if ( ttbarState.mode() != PartonTop::CH_SEMILEPTON and jets.size() < 4 ) vetoEvent;
    else if ( ttbarState.mode() != PartonTop::CH_FULLLEPTON and jets.size() < 2 ) vetoEvent;

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
    if ( ttbarState.mode() == PartonTop::CH_SEMILEPTON ) {
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
    else if ( ttbarState.mode() == PartonTop::CH_FULLLEPTON ) {
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

    normalize(_h14_diffXSecTopSemiLepPartontopPt        );
    normalize(_h15_diffXSecTopSemiLepPartontopPtTtbarSys);
    normalize(_h16_diffXSecTopSemiLepPartontopY         );
    normalize(_h17_diffXSecTopSemiLepPartonttbarDelPhi  );
    normalize(_h18_diffXSecTopSemiLepPartontopPtLead    );
    normalize(_h19_diffXSecTopSemiLepPartontopPtSubLead );
    normalize(_h20_diffXSecTopSemiLepPartonttbarPt      );
    normalize(_h21_diffXSecTopSemiLepPartonttbarY       );
    normalize(_h22_diffXSecTopSemiLepPartonttbarMass    );

    normalize(_h23_diffXSecTopDiLepPartontopPt        );
    normalize(_h24_diffXSecTopDiLepPartontopPtTtbarSys);
    normalize(_h25_diffXSecTopDiLepPartontopY         );
    normalize(_h26_diffXSecTopDiLepPartonttbarDelPhi  );
    normalize(_h27_diffXSecTopDiLepPartontopPtLead    );
    normalize(_h28_diffXSecTopDiLepPartontopPtSubLead );
    normalize(_h29_diffXSecTopDiLepPartonttbarPt      );
    normalize(_h30_diffXSecTopDiLepPartonttbarY       );
    normalize(_h31_diffXSecTopDiLepPartonttbarMass    );

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

  Histo1DPtr _h14_diffXSecTopSemiLepPartontopPt        ;
  Histo1DPtr _h15_diffXSecTopSemiLepPartontopPtTtbarSys;
  Histo1DPtr _h16_diffXSecTopSemiLepPartontopY         ;
  Histo1DPtr _h17_diffXSecTopSemiLepPartonttbarDelPhi  ;
  Histo1DPtr _h18_diffXSecTopSemiLepPartontopPtLead    ;
  Histo1DPtr _h19_diffXSecTopSemiLepPartontopPtSubLead ;
  Histo1DPtr _h20_diffXSecTopSemiLepPartonttbarPt      ;
  Histo1DPtr _h21_diffXSecTopSemiLepPartonttbarY       ;
  Histo1DPtr _h22_diffXSecTopSemiLepPartonttbarMass    ;

  Histo1DPtr _h23_diffXSecTopDiLepPartontopPt        ;
  Histo1DPtr _h24_diffXSecTopDiLepPartontopPtTtbarSys;
  Histo1DPtr _h25_diffXSecTopDiLepPartontopY         ;
  Histo1DPtr _h26_diffXSecTopDiLepPartonttbarDelPhi  ;
  Histo1DPtr _h27_diffXSecTopDiLepPartontopPtLead    ;
  Histo1DPtr _h28_diffXSecTopDiLepPartontopPtSubLead ;
  Histo1DPtr _h29_diffXSecTopDiLepPartonttbarPt      ;
  Histo1DPtr _h30_diffXSecTopDiLepPartonttbarY       ;
  Histo1DPtr _h31_diffXSecTopDiLepPartonttbarMass    ;

};

// This global object acts as a hook for the plugin system
AnalysisBuilder<CMS_TOP_12_028> plugin_CMS_TOP_12_028;

}
