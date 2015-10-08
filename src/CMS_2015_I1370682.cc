#include "Rivet/Analysis.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "TopMonteCarlo/RivetTop/interface/PseudoTop.hh"

namespace Rivet {

class CMS_2015_I1370682 : public Analysis {
public:
  CMS_2015_I1370682() : Analysis("CMS_2015_I1370682"),
    _applyCorrection(false), _doShapeOnly(false) {
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

    // Get the ttbar candidate
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
      const double ch15[] = { 1.332861, 1.202524, 1.015941, 0.825966, 0.679034, 0.578274, 0.534191, 0.535590, };
      const double ch16[] = { 1.330970, 1.204910, 0.994211, 0.785767, 0.637293, 0.546148, 0.518764, 0.531851, };
      const double ch17[] = { 2.439632, 1.107419, 0.932594, 0.877266, 0.857611, 0.858057, 0.876822, 0.927184, 1.101292, 2.431395, };
      const double ch18[] = { 1.071972, 0.987473, 0.946021, 1.026136, };
      const double ch19[] = { 1.505596, 1.279310, 1.076561, 0.868829, 0.704581, 0.589095, 0.533348, 0.529244, };
      const double ch20[] = { 1.266209, 1.142433, 0.950088, 0.766105, 0.635270, 0.556113, 0.536189, 0.552520, };
      const double ch21[] = { 1.473991, 0.919882, 0.867136, 0.877298, 0.868842, 0.830410, };
      const double ch22[] = { 2.910018, 1.170504, 0.920888, 0.826197, 0.786376, 0.784372, 0.822743, 0.917418, 1.164525, 2.893186, };
      const double ch23[] = { 1.739299, 1.371855, 0.985983, 0.737965, 0.567306, 0.431933, 0.300742, };

      const double ch24[] = { 0.999101, 1.045283, 1.024105, 0.897331, 0.768020, };
      const double ch25[] = { 0.983285, 1.050363, 1.031595, 0.894776, 0.766134, 0.661792, };
      const double ch26[] = { 1.757030, 1.017485, 0.913990, 0.898459, 0.897156, 0.916592, 1.012075, 1.771837, };
      const double ch27[] = { 0.971477, 0.958043, 0.976221, 1.090951, };
      const double ch28[] = { 1.022223, 1.045561, 1.030319, 0.907276, 0.776254, };
      const double ch29[] = { 0.979248, 1.045432, 1.014236, 0.878637, 0.749696, };
      const double ch30[] = { 1.146521, 0.892631, 0.942240, 0.983815, };
      const double ch31[] = { 2.547908, 1.161046, 0.927098, 0.866168, 0.865782, 0.925177, 1.162697, 2.549990, };
      const double ch32[] = { 0.884560, 1.121415, 1.091602, 1.029943, 0.970465, 0.913087, };

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

    if ( _doShapeOnly ) {
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
    }
    else {
      const double s = 1./sumOfWeights();
      scale(_hSL_topPt        , s);
      scale(_hSL_topPtTtbarSys, s);
      scale(_hSL_topY         , s);
      scale(_hSL_ttbarDelPhi  , s);
      scale(_hSL_topPtLead    , s);
      scale(_hSL_topPtSubLead , s);
      scale(_hSL_ttbarPt      , s);
      scale(_hSL_ttbarY       , s);
      scale(_hSL_ttbarMass    , s);

      scale(_hDL_topPt        , s);
      scale(_hDL_topPtTtbarSys, s);
      scale(_hDL_topY         , s);
      scale(_hDL_ttbarDelPhi  , s);
      scale(_hDL_topPtLead    , s);
      scale(_hDL_topPtSubLead , s);
      scale(_hDL_ttbarPt      , s);
      scale(_hDL_ttbarY       , s);
      scale(_hDL_ttbarMass    , s);
    }

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
  const bool _applyCorrection, _doShapeOnly;

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
