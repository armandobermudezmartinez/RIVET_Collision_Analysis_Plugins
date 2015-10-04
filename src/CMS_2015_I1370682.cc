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
  CMS_2015_I1370682() : Analysis("CMS_2015_I1370682"),
    _applyCorrection(true) {
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
      if ( lCand1.pt() < 33 or std::abs(lCand1.eta()) > 2.1 ) vetoEvent;

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
    if ( _applyCorrection ) {
      // Correction functions for TOP-12-028 paper, (parton bin height)/(pseudotop bin height)
      const double ch15[] = { 1.332751, 1.202828, 1.016043, 0.825949, 0.678531, 0.578008, 0.534732, 0.535978, };
      const double ch16[] = { 1.331849, 1.204749, 0.993864, 0.785609, 0.637532, 0.546045, 0.518751, 0.533375, };
      const double ch17[] = { 2.436009, 1.107317, 0.932163, 0.876795, 0.857605, 0.858273, 0.876685, 0.927404, 1.103424, 2.434089, };
      const double ch18[] = { 1.072753, 0.987128, 0.946092, 1.026064, };
      const double ch19[] = { 1.505408, 1.279948, 1.076546, 0.869024, 0.703783, 0.588815, 0.533738, 0.530076, };
      const double ch20[] = { 1.266092, 1.142511, 0.950296, 0.765820, 0.635268, 0.555882, 0.537098, 0.551669, };
      const double ch21[] = { 1.474220, 0.919568, 0.867454, 0.877010, 0.869179, 0.830976, };
      const double ch22[] = { 2.907781, 1.169675, 0.920911, 0.825742, 0.785750, 0.783991, 0.823403, 0.918792, 1.165880, 2.898277, };
      const double ch23[] = { 1.739471, 1.372242, 0.986298, 0.736899, 0.566952, 0.432300, 0.301161, };
      const double ch24[] = { 0.997789, 1.044846, 1.024955, 0.899288, 0.769006, };
      const double ch25[] = { 0.981824, 1.049661, 1.032651, 0.897714, 0.766686, 0.659779, };
      const double ch26[] = { 1.758026, 1.017009, 0.914441, 0.898547, 0.897409, 0.916156, 1.012132, 1.767733, };
      const double ch27[] = { 0.971406, 0.958558, 0.976420, 1.090009, };
      const double ch28[] = { 1.020914, 1.044463, 1.031316, 0.909401, 0.776720, };
      const double ch29[] = { 0.978018, 1.045177, 1.015406, 0.880289, 0.751794, };
      const double ch30[] = { 1.144800, 0.893456, 0.943555, 0.982827, };
      const double ch31[] = { 2.551711, 1.160696, 0.926956, 0.867168, 0.865276, 0.924450, 1.162228, 2.546409, };
      const double ch32[] = { 0.882339, 1.120747, 1.092742, 1.032253, 0.971692, 0.918227, };

      applyCorrection(_hSL_topPt, ch15);
      applyCorrection(_hSL_topPtTtbarSys, ch16);
      applyCorrection(_hSL_topY, ch17);
      applyCorrection(_hSL_ttbarDelPhi, ch18);
      applyCorrection(_hSL_topPtLead, ch19);
      applyCorrection(_hSL_topPtSubLead, ch20);
      applyCorrection(_hSL_ttbarPt, ch21);
      applyCorrection(_hSL_ttbarY, ch22);
      applyCorrection(_hSL_ttbarMass, ch23);

      applyCorrection(_hDL_topPt, ch24);
      applyCorrection(_hDL_topPtTtbarSys, ch25);
      applyCorrection(_hDL_topY, ch26);
      applyCorrection(_hDL_ttbarDelPhi, ch27);
      applyCorrection(_hDL_topPtLead, ch28);
      applyCorrection(_hDL_topPtSubLead, ch29);
      applyCorrection(_hDL_ttbarPt, ch30);
      applyCorrection(_hDL_ttbarY, ch31);
      applyCorrection(_hDL_ttbarMass, ch32);
    }

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
  const bool _applyCorrection;

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
