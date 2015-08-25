#include "Rivet/Analysis.hh"
//#include "Rivet/AnalysisLoader.hh"
//#include "Rivet/Particle.fhh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "TopMonteCarlo/RivetTop/interface/PseudoTopLesHouches.hh"

namespace Rivet {

class CMS_LesHouches2015 : public Analysis {
public:
  CMS_LesHouches2015() : Analysis("CMS_LesHouches2015") {
  }

  void init() {
    addProjection(PseudoTopLesHouches(), "ttbar");

    // Booking of histograms
    _h_njets  = bookHisto1D("jet_mult", 11, -0.5, 10.5);
    _h_nbjets = bookHisto1D("bjet_mult", 6, -0.5, 5.5);
    _h_nljets = bookHisto1D("ljet_mult", 6, -0.5, 5.5);
 
    _h_jet_1_pT = bookHisto1D("jet_1_pT", logspace(50, 20., 1000.));
    _h_jet_2_pT = bookHisto1D("jet_2_pT", logspace(50, 20., 700.));
    _h_jet_3_pT = bookHisto1D("jet_3_pT", logspace(50, 20., 500.));
    _h_jet_4_pT = bookHisto1D("jet_4_pT", logspace(50, 20., 400.));
    _h_bjet_1_pT = bookHisto1D("jetb_1_pT", logspace(50, 20., 1000.));
    _h_bjet_2_pT = bookHisto1D("jetb_2_pT", logspace(50, 20., 700.));
    _h_bjet_3_pT = bookHisto1D("jetb_3_pT", logspace(50, 20., 500.));
    _h_bjet_4_pT = bookHisto1D("jetb_4_pT", logspace(50, 20., 400.));
    _h_ljet_1_pT = bookHisto1D("jetl_1_pT", logspace(50, 20., 1000.));
    _h_ljet_2_pT = bookHisto1D("jetl_2_pT", logspace(50, 20., 700.));
    _h_ljet_3_pT = bookHisto1D("jetl_3_pT", logspace(50, 20., 500.));
    _h_ljet_4_pT = bookHisto1D("jetl_4_pT", logspace(50, 20., 400.));

    _h_top_pT = bookHisto1D("top_pT", logspace(50, 10., 500));
    _h_top_mass = bookHisto1D("top_mass", 100, 0., 500.);
    _h_ttbar_pT = bookHisto1D("ttbar_pT", logspace(50, 10., 1000.));
    _h_ttbar_mass = bookHisto1D("ttbar_mass", 200, 0., 2000.);
    _h_ttjet_pT = bookHisto1D("ttjet_pT", logspace(50, 10., 1000.));
    _h_b1b2_mass = bookHisto1D("b1b2_mass", 100, 0, 1000);

    _h_lep_1_lep_2_dphi = bookHisto1D("lep_1_lep_2_dphi", 32, 0.0, 3.2);
    _h_jet_1_jet_2_dphi = bookHisto1D("jet_1_jet_2_dphi", 32, 0.0, 3.2);
  
    _h_nLepton_deta = bookHisto1D("n_lepton_deta", 2, -1., 1.);
    _h_nTtbar_dy = bookHisto1D("n_ttbar_dy", 2, -1., 1.);
    _h_nNoAddJet = bookHisto1D("n_noaddjet", 2, -0.5, 1.5);

    _h_eta_lplus = bookHisto1D("eta_lplus",50,-2.5,2.5);
    _h_eta_lminus = bookHisto1D("eta_lminus",50,-2.5,2.5);
    _h_y_tplus = bookHisto1D("eta_tplus",50,-5.,5.);
    _h_y_tminus = bookHisto1D("eta_tminus",50,-5.,5.);

    _h_ttbar_mass_deta_pos = bookHisto1D("ttbar_mass_deta_pos", 20, 0., 2000.);
    _h_ttbar_mass_deta_neg = bookHisto1D("ttbar_mass_deta_neg", 20, 0., 2000.);
    _h_ttbar_mass_dy_pos = bookHisto1D("ttbar_mass_dy_pos", 20, 0., 2000.);
    _h_ttbar_mass_dy_neg = bookHisto1D("ttbar_mass_dy_neg", 20, 0., 2000.);
    _h_ttbar_mass_noaddjet = bookHisto1D("ttbar_mass_noaddjet", 20, 0., 2000.);
    _h_ttbar_mass_addjet = bookHisto1D("ttbar_mass_addjet", 20, 0., 2000.);
    
    _s_lepton_asym_ttbar_mass = bookScatter2D("lepton_asym_ttbar_mass",20,0.,2000.);
    _s_ttbar_asym_ttbar_mass = bookScatter2D("ttbar_asym_ttbar_mass",20,0.,2000.);
    _s_gap_fraction_ttbar_mass = bookScatter2D("gap_fraction_ttbar_mass",20,0.,2000.);
  };
  
  void analyze(const Event& event) {
    const double weight = event.weight();

    // Get the parton level ttbar candidate
    const PseudoTopLesHouches& ttbar = applyProjection<PseudoTopLesHouches>(event, "ttbar");
    if (ttbar.mode() != PseudoTopLesHouches::CH_FULLLEPTON) {
      MSG_DEBUG("Event fail channel topology cuts");
      vetoEvent;
    }

    _h_nNoAddJet->fill(1., weight);
    _h_ttbar_mass_addjet->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).mass(), weight);
    if (ttbar.ljets().size() < 1) {
      _h_nNoAddJet->fill(0., weight);
      _h_ttbar_mass_noaddjet->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).mass(), weight);
    }
      
    _h_njets->fill(ttbar.jets().size(), weight);
    _h_nbjets->fill(ttbar.bjets().size(), weight);
    _h_nljets->fill(ttbar.ljets().size(), weight);
    
    if (ttbar.jets().size() > 0) {
      _h_jet_1_pT->fill(ttbar.jets()[0].momentum().pT(), weight);
      if (ttbar.jets().size() > 1) {
        _h_jet_2_pT->fill(ttbar.jets()[1].momentum().pT(), weight);
        if (ttbar.jets().size() > 2) {
          _h_jet_3_pT->fill(ttbar.jets()[2].momentum().pT(), weight);
          if (ttbar.jets().size() > 3) {
            _h_jet_4_pT->fill(ttbar.jets()[3].momentum().pT(), weight);
          }
        }
      }
    }
    if (ttbar.bjets().size() > 0) {
      _h_bjet_1_pT->fill(ttbar.bjets()[0].momentum().pT(), weight);
      if (ttbar.bjets().size() > 1) {
        _h_bjet_2_pT->fill(ttbar.bjets()[1].momentum().pT(), weight);
        if (ttbar.bjets().size() > 2) {
          _h_bjet_3_pT->fill(ttbar.bjets()[2].momentum().pT(), weight);
          if (ttbar.bjets().size() > 3) {
            _h_bjet_4_pT->fill(ttbar.bjets()[3].momentum().pT(), weight);
          }
        }
      }
    }
    if (ttbar.ljets().size() > 0) {
      _h_ljet_1_pT->fill(ttbar.ljets()[0].momentum().pT(), weight);
      if (ttbar.ljets().size() > 1) {
        _h_ljet_2_pT->fill(ttbar.ljets()[1].momentum().pT(), weight);
        if (ttbar.ljets().size() > 2) {
          _h_ljet_3_pT->fill(ttbar.ljets()[2].momentum().pT(), weight);
          if (ttbar.ljets().size() > 3) {
            _h_ljet_4_pT->fill(ttbar.ljets()[3].momentum().pT(), weight);
          }
        }
      }
    }
    
    _h_top_pT->fill(ttbar.t1().momentum().pT(), weight);
    _h_top_pT->fill(ttbar.t2().momentum().pT(), weight);
    _h_top_mass->fill(ttbar.t1().momentum().mass(), weight);
    _h_top_mass->fill(ttbar.t2().momentum().mass(), weight);
    _h_ttbar_pT->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).pT(), weight);
    _h_ttbar_mass->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).mass(), weight);
    if (ttbar.ljets().size() > 0) 
      _h_ttjet_pT->fill((ttbar.t1().momentum()+ttbar.t2().momentum()+ttbar.ljets()[0].momentum()).pT(), weight);
    _h_b1b2_mass->fill((ttbar.b1().momentum()+ttbar.b2().momentum()).mass(), weight);

    _h_lep_1_lep_2_dphi->fill(deltaPhi(ttbar.wDecays1()[0].momentum(),ttbar.wDecays2()[0].momentum()), weight);
    if (ttbar.jets().size() > 1)
      _h_jet_1_jet_2_dphi->fill(deltaPhi(ttbar.jets()[0].momentum(),ttbar.jets()[1].momentum()), weight);

    double lepton_deta = 0.;
    if (ttbar.wDecays1()[0].charge() > 0) { 
      lepton_deta = ttbar.wDecays1()[0].momentum().abseta() - ttbar.wDecays2()[0].momentum().abseta();
      _h_eta_lplus->fill(ttbar.wDecays1()[0].momentum().eta());
      _h_eta_lminus->fill(ttbar.wDecays2()[0].momentum().eta());
    }
    else if (ttbar.wDecays2()[0].charge() > 0) { 
      lepton_deta = ttbar.wDecays2()[0].momentum().abseta() - ttbar.wDecays1()[0].momentum().abseta();
      _h_eta_lplus->fill(ttbar.wDecays2()[0].momentum().eta());
      _h_eta_lminus->fill(ttbar.wDecays1()[0].momentum().eta());
    }
    if (lepton_deta > 1e-6) {
      _h_nLepton_deta->fill(0.5, weight);
      _h_ttbar_mass_deta_pos->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).mass(), weight);
    }
    if (lepton_deta < -1e-6) {
      _h_nLepton_deta->fill(-0.5, weight);
      _h_ttbar_mass_deta_neg->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).mass(), weight);
    }

    double ttbar_dy = 0.;
    if (ttbar.t1().charge() > 0) {
      ttbar_dy = ttbar.t1().momentum().absrapidity() - ttbar.t2().momentum().absrapidity();
      _h_y_tplus->fill(ttbar.t1().momentum().rapidity());
      _h_y_tminus->fill(ttbar.t2().momentum().rapidity());
      }
    else if (ttbar.t2().charge() > 0) {
      ttbar_dy = ttbar.t2().momentum().absrapidity() - ttbar.t1().momentum().absrapidity();
      _h_y_tplus->fill(ttbar.t2().momentum().rapidity());
      _h_y_tminus->fill(ttbar.t1().momentum().rapidity());
      }
    if (ttbar_dy > 1e-6) {
      _h_nTtbar_dy->fill(0.5, weight);
      _h_ttbar_mass_dy_pos->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).mass(), weight);
    }
    if (ttbar_dy < -1e-6) {
      _h_nTtbar_dy->fill(-0.5, weight);
      _h_ttbar_mass_dy_neg->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).mass(), weight);
    }
  };
  
  void finalize() {
    asymm(_h_ttbar_mass_deta_pos, _h_ttbar_mass_deta_neg, _s_lepton_asym_ttbar_mass);
    asymm(_h_ttbar_mass_dy_pos, _h_ttbar_mass_dy_neg, _s_ttbar_asym_ttbar_mass);
    divide(_h_ttbar_mass_noaddjet, _h_ttbar_mass_addjet, _s_gap_fraction_ttbar_mass);
   
    normalize(_h_njets);
    normalize(_h_nbjets);
    normalize(_h_nljets);
    normalize(_h_jet_1_pT);
    normalize(_h_jet_2_pT);
    normalize(_h_jet_3_pT);
    normalize(_h_jet_4_pT);
    normalize(_h_bjet_1_pT);
    normalize(_h_bjet_2_pT);
    normalize(_h_bjet_3_pT);
    normalize(_h_bjet_4_pT);
    normalize(_h_ljet_1_pT);
    normalize(_h_ljet_2_pT);
    normalize(_h_ljet_3_pT);
    normalize(_h_ljet_4_pT);
    normalize(_h_eta_lplus);
    normalize(_h_eta_lminus);
    normalize(_h_y_tplus);
    normalize(_h_y_tminus);
    normalize(_h_top_pT);
    normalize(_h_top_mass);
    normalize(_h_ttbar_pT);
    normalize(_h_ttbar_mass);
    normalize(_h_ttjet_pT);
    normalize(_h_b1b2_mass);
    normalize(_h_lep_1_lep_2_dphi);
    normalize(_h_jet_1_jet_2_dphi);
  };

private:
  
  Histo1DPtr _h_njets, _h_nbjets, _h_nljets;

  Histo1DPtr _h_jet_1_pT, _h_jet_2_pT, _h_jet_3_pT, _h_jet_4_pT;
  Histo1DPtr _h_ljet_1_pT, _h_ljet_2_pT, _h_ljet_3_pT, _h_ljet_4_pT;
  Histo1DPtr _h_bjet_1_pT, _h_bjet_2_pT, _h_bjet_3_pT, _h_bjet_4_pT;
  
  Histo1DPtr _h_top_pT, _h_top_mass, _h_ttbar_pT, _h_ttbar_mass, _h_ttjet_pT, _h_b1b2_mass;
 
  Histo1DPtr _h_lep_1_lep_2_dphi, _h_jet_1_jet_2_dphi;
 
  Histo1DPtr _h_nLepton_deta, _h_eta_lplus, _h_eta_lminus, _h_ttbar_mass_deta_pos, _h_ttbar_mass_deta_neg;
  Histo1DPtr _h_nTtbar_dy, _h_y_tplus, _h_y_tminus, _h_ttbar_mass_dy_pos, _h_ttbar_mass_dy_neg;
  Histo1DPtr _h_nNoAddJet, _h_ttbar_mass_noaddjet, _h_ttbar_mass_addjet;

  Scatter2DPtr _s_lepton_asym_ttbar_mass, _s_ttbar_asym_ttbar_mass, _s_gap_fraction_ttbar_mass; 

};

// This global object acts as a hook for the plugin system
AnalysisBuilder<CMS_LesHouches2015> plugin_CMS_LesHouches2015;

}
