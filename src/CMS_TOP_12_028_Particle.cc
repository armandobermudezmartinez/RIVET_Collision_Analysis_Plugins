#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Cuts.hh"

// This is an approximate version of TOP-12-028

namespace Rivet {

class CMS_TOP_12_028_Particle : public Analysis {
public:
  CMS_TOP_12_028_Particle() : Analysis("CMS_TOP_12_028_Particle") {
  }

  void init() {
    FinalState fs(-5.0, 5.0, 0*GeV);

    WFinder wefs(fs, Cuts::open(), PID::ELECTRON, 0*GeV, 160*GeV, 0);
    WFinder wmfs(fs, Cuts::open(), PID::MUON, 0*GeV, 160*GeV, 0);
    addProjection(wefs, "wefs");
    addProjection(wmfs, "wmfs");

    VetoedFinalState fsForJets(fs);
    fsForJets.addDecayProductsVeto(+24);
    fsForJets.addDecayProductsVeto(-24);

    FastJets fj(fsForJets, FastJets::ANTIKT, 0.5);
    fj.useInvisibles();
    addProjection(fj, "JetsParticle");

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

  void analyze(const Event& event) {
    const double weight = event.weight();

    const Particles elFromW = applyProjection<WFinder>(event, "wefs").constituentLeptons();
    const Particles muFromW = applyProjection<WFinder>(event, "wmfs").constituentLeptons();

    Particles leptons;
    if ( true ) {
      Particles leptonsForSemiLep, leptonsForDiLep;
      foreach ( const Particle& p, elFromW ) {
        const double abseta = std::abs(p.eta());
        const double pt = p.pT();
        if ( abseta < 2.1 && pt > 33*GeV ) leptonsForSemiLep.push_back(p);
        if ( abseta < 2.4 && pt > 20*GeV ) leptonsForDiLep.push_back(p);
      }
      foreach ( const Particle& p, muFromW ) {
        const double abseta = std::abs(p.eta());
        const double pt = p.pT();
        if ( abseta < 2.1 && pt > 33*GeV ) leptonsForSemiLep.push_back(p);
        if ( abseta < 2.4 && pt > 20*GeV ) leptonsForDiLep.push_back(p);
      }

      if ( leptonsForDiLep.size() >= 2 ) {
        leptons.push_back(Particle());
        leptons.push_back(Particle());
        foreach ( Particle& lepton, leptonsForDiLep ) {
          if ( leptons[0].pT() < lepton.pT() ) {
            leptons[1] = leptons[0];
            leptons[0] = lepton;
          }
          else if ( leptons[1].pT() < lepton.pT() ) {
            leptons[1] = lepton;
          }
        }
      }
      else if ( leptonsForSemiLep.size() == 1 ) {
        leptons.push_back(leptonsForSemiLep[0]);
      }
    }
    if ( leptons.empty() ) vetoEvent; // dilepton or semilepton channel
    const int nLepton = leptons.size();

    // Build genJets
    const Jets& jetsIn = applyProjection<JetAlg>(event, "JetsParticle").jetsByPt(30*GeV, MAXDOUBLE, -2.4, 2.4);
    Jets jets;
    // Check overlap with leptons
    foreach ( const Jet& jet, jetsIn ) {
      bool isOverlap = false;
      foreach ( const Particle& lep, leptons ) {
        if ( deltaR(lep.momentum(), jet.momentum()) < 0.3 ) {
          isOverlap = true;
          break;
        }
      }
      if ( !isOverlap ) jets.push_back(jet);
    }
    const int nJet = jets.size();
    if ( nLepton == 1 && nJet < 4 ) vetoEvent;
    if ( nLepton == 2 && nJet < 2 ) vetoEvent;

    // Get leading two b-jets. Note that jets are already sorted by pT.
    const Jet* bjet1 = 0, * bjet2 = 0;
    foreach (const Jet& jet, jets) {
      if ( !jet.containsBottom() ) continue;
      if ( !bjet1 ) {
        bjet1 = &jet;
      } else if ( !bjet2 ) {
        bjet2 = &jet;
        break;
      }
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
    if ( leptons.size() == 1 ) {
      // Do the semileptonic channel
      const FourMomentum& lP4 = leptons[0].momentum();
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
    else {
      // Do the dileptonic channel
      const FourMomentum l1P4 = leptons[0].momentum();
      const FourMomentum l2P4 = leptons[1].momentum();

      const FourMomentum dilP4 = l1P4+l2P4;
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
AnalysisBuilder<CMS_TOP_12_028_Particle> plugin_CMS_TOP_12_028_Particle;

}
