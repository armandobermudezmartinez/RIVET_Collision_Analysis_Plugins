#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"


namespace Rivet {

  class CMS_2015_I1370682_parton : public Analysis {
    public:
      CMS_2015_I1370682_parton() : Analysis("CMS_2015_I1370682_parton") {}

      void init() {
        // Parton level top quarks
        declare(PartonicTops(PartonicTops::E_MU, false), "LeptonicPartonTops");
        declare(PartonicTops(PartonicTops::HADRONIC),    "HadronicPartonTops");

        // Find jets not related to the top/W decays
        VetoedFinalState vfs;
        vfs.addDecayProductsVeto(PID::WPLUSBOSON);
        vfs.addDecayProductsVeto(PID::WMINUSBOSON);
        FastJets fj(vfs, FastJets::ANTIKT, 0.5, JetAlg::ALL_MUONS, JetAlg::ALL_INVISIBLES);
        declare(fj, "Jets");

        _hSL_lepPt    = bookHisto1D("d01-x01-y01");
        _hSL_lepEta   = bookHisto1D("d02-x01-y01");
        _hSL_bqPt     = bookHisto1D("d03-x01-y01");
        _hSL_bqEta    = bookHisto1D("d04-x01-y01");
        _hSL_bbbarPt  = bookHisto1D("d05-x01-y01");
        _hSL_bbbarMass= bookHisto1D("d06-x01-y01");

        _hDL_lepPt      = bookHisto1D("d07-x01-y01");
        _hDL_lepEta     = bookHisto1D("d08-x01-y01");
        _hDL_dilepPt    = bookHisto1D("d09-x01-y01");
        _hDL_dilepMass  = bookHisto1D("d10-x01-y01");
        _hDL_bqPt       = bookHisto1D("d11-x01-y01");
        _hDL_bqEta      = bookHisto1D("d12-x01-y01");
        _hDL_bbbarPt    = bookHisto1D("d13-x01-y01");
        _hDL_bbbarMass  = bookHisto1D("d14-x01-y01");

        _hSL_topPt         = bookHisto1D("d15-x01-y01");
        _hSL_topPtTtbarSys = bookHisto1D("d16-x01-y01");
        _hSL_topY          = bookHisto1D("d17-x01-y01");
        _hSL_ttbarDelPhi   = bookHisto1D("d18-x01-y01");
        _hSL_topPtLead     = bookHisto1D("d19-x01-y01");
        _hSL_topPtSubLead  = bookHisto1D("d20-x01-y01");
        _hSL_ttbarPt       = bookHisto1D("d21-x01-y01");
        _hSL_ttbarY        = bookHisto1D("d22-x01-y01");
        _hSL_ttbarMass     = bookHisto1D("d23-x01-y01");

        _hDL_topPt           = bookHisto1D("d24-x01-y01");
        _hDL_topPtTtbarSys   = bookHisto1D("d25-x01-y01");
        _hDL_topY            = bookHisto1D("d26-x01-y01");
        _hDL_ttbarDelPhi     = bookHisto1D("d27-x01-y01");
        _hDL_topPtLead       = bookHisto1D("d28-x01-y01");
        _hDL_topPtSubLead    = bookHisto1D("d29-x01-y01");
        _hDL_ttbarPt         = bookHisto1D("d30-x01-y01");
        _hDL_ttbarY          = bookHisto1D("d31-x01-y01");
        _hDL_ttbarMass       = bookHisto1D("d32-x01-y01");

      };

      void analyze(const Event& event) {
        const double weight = event.weight();
        
        bool isDilepton   = false;
        bool isSemilepton = false;
        
        // Do the analysis only for the ttbar full leptonic or semileptonic channel, without tau decay
        const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();
        const Particles hadronicpartontops = apply<ParticleFinder>(event, "HadronicPartonTops").particlesByPt();
        if (leptonicpartontops.size() == 1 and hadronicpartontops.size() == 1)
          isSemilepton = true;
        else if (leptonicpartontops.size() == 2)
          isDilepton = true;
        else
          vetoEvent;

        // Parton level at full phase space
        // Fill top quarks defined in the parton level, full phase space
        FourMomentum t1P4;
        FourMomentum t2P4;
        if (isSemilepton) {
          t1P4 = leptonicpartontops[0].momentum();
          t2P4 = hadronicpartontops[0].momentum();
        }
        else if (isDilepton) {
          t1P4 = leptonicpartontops[0].momentum();
          t2P4 = leptonicpartontops[1].momentum();
        }
        double t1Pt = t1P4.pT(), t2Pt = t2P4.pT();
        FourMomentum ttbarP4 = t1P4+t2P4;
        FourMomentum t1P4AtCM = LorentzTransform::mkFrameTransformFromBeta(ttbarP4.betaVec()).transform(t1P4);
        double dPhi = deltaPhi(t1P4.phi(), t2P4.phi());

        if (isSemilepton) {
          _hSL_topPt->fill(t1Pt, weight);
          _hSL_topPt->fill(t2Pt, weight);
          _hSL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
          _hSL_topY->fill(t1P4.rapidity(), weight);
          _hSL_topY->fill(t2P4.rapidity(), weight);
          _hSL_ttbarDelPhi->fill(dPhi, weight);
          _hSL_topPtLead->fill(std::max(t1Pt, t2Pt), weight);
          _hSL_topPtSubLead->fill(std::min(t1Pt, t2Pt), weight);
          _hSL_ttbarPt->fill(ttbarP4.pT(), weight);
          _hSL_ttbarY->fill(ttbarP4.rapidity(), weight);
          _hSL_ttbarMass->fill(ttbarP4.mass(), weight);
        }
        else if (isDilepton) {
          _hDL_topPt->fill(t1Pt, weight);
          _hDL_topPt->fill(t2Pt, weight);
          _hDL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
          _hDL_topY->fill(t1P4.rapidity(), weight);
          _hDL_topY->fill(t2P4.rapidity(), weight);
          _hDL_ttbarDelPhi->fill(dPhi, weight);
          _hDL_topPtLead->fill(std::max(t1Pt, t2Pt), weight);
          _hDL_topPtSubLead->fill(std::min(t1Pt, t2Pt), weight);
          _hDL_ttbarPt->fill(ttbarP4.pT(), weight);
          _hDL_ttbarY->fill(ttbarP4.rapidity(), weight);
          _hDL_ttbarMass->fill(ttbarP4.mass(), weight);
        }

        // Find leptons
        Particles lCands;
        const auto isPromptChLepton = [](const Particle& p){return isChargedLepton(p) && !fromDecay(p);};
        for (const Particle& top : leptonicpartontops) {
          lCands.push_back(top.allDescendants(lastParticleWith(isPromptChLepton)).front());
        }
        if (isSemilepton) {
          // Apply the particle level phase space cut
          if ( lCands[0].pT() <= 33 or std::abs(lCands[0].eta()) >= 2.1 ) vetoEvent;
        }
        else if (isDilepton) {
          if ( lCands[0].pT() < lCands[1].pT() ) std::swap(lCands[0], lCands[1]);
          const double l1Pt = lCands[0].pT(), l1Abseta = std::abs(lCands[0].eta());
          const double l2Pt = lCands[1].pT(), l2Abseta = std::abs(lCands[1].eta());

          // Apply the particle level phase space cut
          if ( l1Pt <= 20 or l1Abseta >= 2.4 or l2Pt <= 20 or l2Abseta >= 2.4 ) vetoEvent;
          if ( (lCands[0].momentum()+lCands[1].momentum()).mass() < 20 ) vetoEvent;
        }

        // Build genJets
        const Jets& jetsIn = apply<JetAlg>(event, "Jets").jetsByPt(30*GeV);
        Jets jets;
        for (const Jet& jet : jetsIn) {
          if ( std::abs(jet.eta()) > 2.4 ) continue;
          bool isOverlapped = false;
          foreach ( const Particle& lCand, lCands ) {
            double lepJetMinDR = isSemilepton ? 0.4 : 0.3;
            if ( deltaR(lCand.momentum(), jet.momentum()) < lepJetMinDR ) {
              isOverlapped = true;
              break;
            }
          }
          if ( isOverlapped ) continue;
          jets.push_back(jet);
        }
        if (isSemilepton and jets.size() < 4) vetoEvent;
        else if (isDilepton and jets.size() < 2) vetoEvent;

        // Require at least 2 b jets
        std::vector<FourMomentum> bjets;
        for (const Jet& jet : jets) {
          if (jet.bTagged()) bjets.push_back(jet.momentum());
        }
        if (bjets.size() < 2) vetoEvent;

        const FourMomentum& b1P4 = bjets[0];
        const FourMomentum& b2P4 = bjets[1];
        const FourMomentum bbP4 = b1P4+b2P4;
        const double b1Pt = b1P4.pT(), b1Eta = b1P4.eta();
        const double b2Pt = b2P4.pT(), b2Eta = b2P4.eta();
        const double bbPt = bbP4.pT(), bbMass = bbP4.mass();

        // Find leptons, apply channel dependent phase space cuts, fill histograms
        if (isSemilepton) {
          // Do the semileptonic channel
          const FourMomentum& lP4 = lCands[0].momentum();
          const double lPt = lP4.pT(), lEta = lP4.eta();

          _hSL_lepPt->fill(lPt, weight);
          _hSL_lepEta->fill(lEta, weight);
          _hSL_bqPt->fill(b1Pt, weight);
          _hSL_bqPt->fill(b2Pt, weight);
          _hSL_bqEta->fill(b1Eta, weight);
          _hSL_bqEta->fill(b2Eta, weight);
          _hSL_bbbarPt->fill(bbPt, weight);
          _hSL_bbbarMass->fill(bbMass, weight);
        }
        else if (isDilepton) {
          // Do the dileptonic channel
          const FourMomentum& l1P4 = lCands[0].momentum();
          const FourMomentum& l2P4 = lCands[1].momentum();
          const FourMomentum dilP4 = l1P4+l2P4;
          const double l1Pt = l1P4.pT(), l1Eta = l1P4.eta();
          const double l2Pt = l2P4.pT(), l2Eta = l2P4.eta();

          _hDL_lepPt->fill(l1Pt, weight);
          _hDL_lepPt->fill(l2Pt, weight);
          _hDL_lepEta->fill(l1Eta, weight);
          _hDL_lepEta->fill(l2Eta, weight);
          _hDL_dilepPt->fill(dilP4.pT(), weight);
          _hDL_dilepMass->fill(dilP4.mass(), weight);

          _hDL_bqPt->fill(b1Pt, weight);
          _hDL_bqPt->fill(b2Pt, weight);
          _hDL_bqEta->fill(b1Eta, weight);
          _hDL_bqEta->fill(b2Eta, weight);
          _hDL_bbbarPt->fill(bbPt, weight);
          _hDL_bbbarMass->fill(bbMass, weight);
        }
      };

      void finalize() {
        normalize(_hSL_lepPt    );
        normalize(_hSL_lepEta   );
        normalize(_hSL_bqPt     );
        normalize(_hSL_bqEta    );
        normalize(_hSL_bbbarPt  );
        normalize(_hSL_bbbarMass);

        normalize(_hDL_lepPt      );
        normalize(_hDL_lepEta     );
        normalize(_hDL_dilepPt    );
        normalize(_hDL_dilepMass  );
        normalize(_hDL_bqPt       );
        normalize(_hDL_bqEta      );
        normalize(_hDL_bbbarPt    );
        normalize(_hDL_bbbarMass  );

        normalize(_hSL_topPt        );
        normalize(_hSL_topPtTtbarSys);
        normalize(_hSL_topY         );
        normalize(_hSL_ttbarDelPhi  );
        normalize(_hSL_topPtLead    );
        normalize(_hSL_topPtSubLead );
        normalize(_hSL_ttbarPt      );
        normalize(_hSL_ttbarY       );
        normalize(_hSL_ttbarMass    );

        normalize(_hDL_topPt        );
        normalize(_hDL_topPtTtbarSys);
        normalize(_hDL_topY         );
        normalize(_hDL_ttbarDelPhi  );
        normalize(_hDL_topPtLead    );
        normalize(_hDL_topPtSubLead );
        normalize(_hDL_ttbarPt      );
        normalize(_hDL_ttbarY       );
        normalize(_hDL_ttbarMass    );

      };

    private:
      Histo1DPtr _hSL_lepPt;
      Histo1DPtr _hSL_lepEta;
      Histo1DPtr _hSL_bqPt;
      Histo1DPtr _hSL_bqEta;
      Histo1DPtr _hSL_bbbarPt;
      Histo1DPtr _hSL_bbbarMass;

      Histo1DPtr _hDL_lepPt;
      Histo1DPtr _hDL_lepEta;
      Histo1DPtr _hDL_dilepPt;
      Histo1DPtr _hDL_dilepMass;
      Histo1DPtr _hDL_bqPt;
      Histo1DPtr _hDL_bqEta;
      Histo1DPtr _hDL_bbbarPt;
      Histo1DPtr _hDL_bbbarMass;

      Histo1DPtr _hSL_topPt        ;
      Histo1DPtr _hSL_topPtTtbarSys;
      Histo1DPtr _hSL_topY         ;
      Histo1DPtr _hSL_ttbarDelPhi  ;
      Histo1DPtr _hSL_topPtLead    ;
      Histo1DPtr _hSL_topPtSubLead ;
      Histo1DPtr _hSL_ttbarPt      ;
      Histo1DPtr _hSL_ttbarY       ;
      Histo1DPtr _hSL_ttbarMass    ;

      Histo1DPtr _hDL_topPt        ;
      Histo1DPtr _hDL_topPtTtbarSys;
      Histo1DPtr _hDL_topY         ;
      Histo1DPtr _hDL_ttbarDelPhi  ;
      Histo1DPtr _hDL_topPtLead    ;
      Histo1DPtr _hDL_topPtSubLead ;
      Histo1DPtr _hDL_ttbarPt      ;
      Histo1DPtr _hDL_ttbarY       ;
      Histo1DPtr _hDL_ttbarMass    ;

  };

  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CMS_2015_I1370682_parton> plugin_CMS_2015_I1370682_parton;

}
