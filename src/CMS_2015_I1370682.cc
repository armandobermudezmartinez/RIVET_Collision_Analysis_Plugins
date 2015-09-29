#include "Rivet/Analysis.hh"
//#include "Rivet/AnalysisLoader.hh"
//#include "Rivet/Particle.fhh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "TopMonteCarlo/RivetTop/interface/PseudoTop.hh"

namespace Rivet {

class CMS_2015_I1370682 : public Analysis {
public:
  CMS_2015_I1370682() : Analysis("CMS_2015_I1370682") {
  }

  void init() {
    addProjection(PseudoTop(0.1, 20, 2.4, 0.5, 30, 2.4), "ttbar");

    // Lepton + Jet channel
    _hSL_topPt         = bookHisto1D("d15-x01-y01"); // 1/sigma dsigma/dpt(top)
    _hSL_topPtTtbarSys = bookHisto1D("d16-x01-y01"); // 1/sigma dsigma/dpt*(top)
    _hSL_topY          = bookHisto1D("d17-x01-y01"); // 1/sigma dsigma/dy(top)
    _hSL_ttbarDelPhi   = bookHisto1D("d18-x01-y01"); // 1/sigma dsigma/ddeltaphi(t,tbar)
    _hSL_topPtLead     = bookHisto1D("d19-x01-y01"); // 1/sigma dsigma/dpt(t1)
    _hSL_topPtSubLead  = bookHisto1D("d20-x01-y01"); // 1/sigma dsigma/dpt(t2)
    _hSL_ttbarPt       = bookHisto1D("d21-x01-y01"); // 1/sigma dsigma/dpt(ttbar)
    _hSL_ttbarY        = bookHisto1D("d22-x01-y01"); // 1/sigma dsigma/dy(ttbar)
    _hSL_ttbarMass     = bookHisto1D("d23-x01-y01"); // 1/sigma dsigma/dm(ttbar)

    // Dilepton channel
    _hDL_topPt         = bookHisto1D("d24-x01-y01"); // 1/sigma dsigma/dpt(top)
    _hDL_topPtTtbarSys = bookHisto1D("d25-x01-y01"); // 1/sigma dsigma/dpt*(top)
    _hDL_topY          = bookHisto1D("d26-x01-y01"); // 1/sigma dsigma/dy(top)
    _hDL_ttbarDelPhi   = bookHisto1D("d27-x01-y01"); // 1/sigma dsigma/ddeltaphi(t,tbar)
    _hDL_topPtLead     = bookHisto1D("d28-x01-y01"); // 1/sigma dsigma/dpt(t1)
    _hDL_topPtSubLead  = bookHisto1D("d29-x01-y01"); // 1/sigma dsigma/dpt(t2)
    _hDL_ttbarPt       = bookHisto1D("d30-x01-y01"); // 1/sigma dsigma/dpt(ttbar)
    _hDL_ttbarY        = bookHisto1D("d31-x01-y01"); // 1/sigma dsigma/dy(ttbar)
    _hDL_ttbarMass     = bookHisto1D("d32-x01-y01"); // 1/sigma dsigma/dm(ttbar)

  };

  void analyze(const Event& event) {
    const double weight = event.weight();

    // Get the parton level ttbar candidate
    const PseudoTop& ttbar = applyProjection<PseudoTop>(event, "ttbar");
    if ( ttbar.mode() == PseudoTop::CH_NONE ) vetoEvent;

    const FourMomentum& t1P4 = ttbar.t1().momentum();
    const FourMomentum& t2P4 = ttbar.t2().momentum();
    const double pt1 = std::max(t1P4.pT(), t2P4.pT());
    const double pt2 = std::min(t1P4.pT(), t2P4.pT());
    const double dPhi = deltaPhi(t1P4, t2P4);
    const FourMomentum ttP4 = t1P4+t2P4;
    const FourMomentum t1P4AtCM = LorentzTransform(-ttP4.boostVector()).transform(t1P4);

    if ( ttbar.mode() == PseudoTop::CH_SEMILEPTON ) {
      const Particle lCand1 = ttbar.wDecays1()[0]; // w1 dau0 is the lepton in the PseudoTop
      if ( lCand1.pt() < 30 or std::abs(lCand1.eta()) > 2.4 ) vetoEvent;

      _hSL_topPt->fill(t1P4.pT(), weight);
      _hSL_topPt->fill(t2P4.pT(), weight);
      _hSL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
      _hSL_topY->fill(t1P4.rapidity(), weight);
      _hSL_topY->fill(t2P4.rapidity(), weight);
      _hSL_ttbarDelPhi->fill(dPhi, weight);
      _hSL_topPtLead->fill(pt1, weight);
      _hSL_topPtSubLead->fill(pt2, weight);
      _hSL_ttbarPt->fill(ttP4.pT(), weight);
      _hSL_ttbarY->fill(ttP4.rapidity(), weight);
      _hSL_ttbarMass->fill(ttP4.mass(), weight);
    }
    else if ( ttbar.mode() == PseudoTop::CH_FULLLEPTON ) {
      const Particle lCand1 = ttbar.wDecays1()[0]; // dau0 are the lepton in the PseudoTop
      const Particle lCand2 = ttbar.wDecays2()[0]; // dau0 are the lepton in the PseudoTop
      if ( lCand1.pt() < 20 or std::abs(lCand1.eta()) > 2.4 or
           lCand2.pt() < 20 or std::abs(lCand2.eta()) > 2.4 ) vetoEvent;

      _hDL_topPt->fill(t1P4.pT(), weight);
      _hDL_topPt->fill(t2P4.pT(), weight);
      _hDL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
      _hDL_topY->fill(t1P4.rapidity(), weight);
      _hDL_topY->fill(t2P4.rapidity(), weight);
      _hDL_ttbarDelPhi->fill(dPhi, weight);
      _hDL_topPtLead->fill(pt1, weight);
      _hDL_topPtSubLead->fill(pt2, weight);
      _hDL_ttbarPt->fill(ttP4.pT(), weight);
      _hDL_ttbarY->fill(ttP4.rapidity(), weight);
      _hDL_ttbarMass->fill(ttP4.mass(), weight);
    }
/*
    const FourMomentum& l1P4 = lCand1.momentum();
    const FourMomentum& l2P4 = lCand2.momentum();
    const FourMomentum& b1P4 = ttbar.b1().momentum();
    const FourMomentum& b2P4 = ttbar.b2().momentum();

    const double l1Pt = l1P4.pT(), l1Eta = l1P4.eta();
    const double l2Pt = l2P4.pT(), l2Eta = l2P4.eta();
    const double b1Pt = b1P4.pT(), b1Eta = b1P4.eta();
    const double b2Pt = b2P4.pT(), b2Eta = b2P4.eta();
    // Leptons and b jet observables are considered only within the particle level phase space
    if ( l1Pt > 20 and l2Pt > 20 and std::abs(l1Eta) < 2.4 and std::abs(l2Eta) < 2.4 and
         b1Pt > 30 and b2Pt > 30 and std::abs(b1Eta) < 2.4 and std::abs(b2Eta) < 2.4 ) {
      _h_lepton_pt->fill(l1Pt, weight);
      _h_lepton_pt->fill(l2Pt, weight);
      _h_lepton_eta->fill(l1Eta, weight);
      _h_lepton_eta->fill(l2Eta, weight);

      _h_bjet_pt->fill(b1Pt, weight);
      _h_bjet_pt->fill(b2Pt, weight);
      _h_bjet_eta->fill(b1Eta, weight);
      _h_bjet_eta->fill(b2Eta, weight);

      const FourMomentum dileptonP4 = l1P4+l2P4;
      _h_dilepton_mass->fill(dileptonP4.mass(), weight);
      _h_dilepton_pt->fill(dileptonP4.pT(), weight);

      const FourMomentum lb11P4 = l1P4 + b1P4;
      const FourMomentum lb22P4 = l2P4 + b2P4;
      const FourMomentum lb12P4 = l1P4 + b2P4;
      const FourMomentum lb21P4 = l2P4 + b1P4;
      _h_lb_mass->fill(lb11P4.mass(), weight);
      _h_lb_mass->fill(lb22P4.mass(), weight);
      _h_lb_mass->fill(lb12P4.mass(), weight);
      _h_lb_mass->fill(lb21P4.mass(), weight);

      const FourMomentum dijetP4 = b1P4 + b2P4;
      _h_dijet_mass->fill(dijetP4.mass(), weight);
      _h_dijet_pt->fill(dijetP4.pT(), weight);
    }
*/

  };

  void finalize() {
    // Correction functions for TOP-12-028 paper, (parton bin height)/(pseudotop bin height)
    const double ch15[] = {1.161660, 1.136480, 1.020996, 0.895649, 0.772136, 0.685911, 0.559711, 0.566430};
    const double ch17[] = {2.101211, 1.099831, 0.937698, 0.883005, 0.868135, 0.882153, 0.878180, 0.941096, 1.095958, 2.056497};
    const double ch21[] = {1.602612, 0.913407, 0.816876, 0.849766, 0.889415, 0.857082};
    const double ch22[] = {2.461665, 1.147150, 0.908031, 0.848166, 0.814687, 0.803214, 0.824948, 0.947269, 1.122359, 2.428979};
    const double ch23[] = {1.498358, 1.362128, 1.024490, 0.819021, 0.646227, 0.475925, 0.372441};

    const double ch24[] = {0.933825, 1.069645, 1.051336, 0.919932, 0.774565};
    const double ch26[] = {1.682022, 1.002849, 0.925246, 0.924734, 0.880097, 0.901330, 1.042041, 1.733911};
    const double ch30[] = {1.129278, 0.908123, 0.933110, 0.963850};
    const double ch31[] = {2.401265, 1.140515, 0.937143, 0.889803, 0.833903, 0.946386, 1.179555, 2.445021};
    const double ch32[] = {0.803342, 1.136017, 1.206834, 1.037619, 1.081579, 0.741247};

    applyCorrection(_hSL_topPt, ch15);
    applyCorrection(_hSL_topY, ch17);
    applyCorrection(_hSL_ttbarPt, ch21);
    applyCorrection(_hSL_ttbarY, ch22);
    applyCorrection(_hSL_ttbarMass, ch23);

    applyCorrection(_hDL_topPt, ch24);
    applyCorrection(_hDL_topY, ch26);
    applyCorrection(_hDL_ttbarPt, ch30);
    applyCorrection(_hDL_ttbarY, ch31);
    applyCorrection(_hDL_ttbarMass, ch32);

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

  void applyCorrection(Histo1DPtr h, const double* cf) {
    std::vector<YODA::HistoBin1D>& bins = h->bins();
    for ( int i=0, n=bins.size(); i<n; ++i ) {
      const double s = cf[i];
      YODA::HistoBin1D& bin = bins[i];
      bin.scaleW(s);
    }
  };

private:
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
AnalysisBuilder<CMS_2015_I1370682> plugin_CMS_2015_I1370682;

}
