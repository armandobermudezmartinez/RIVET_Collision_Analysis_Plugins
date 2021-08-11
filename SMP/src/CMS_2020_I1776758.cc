// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2020_I1776758 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2020_I1776758);


    // based on SMP_19_004
    
    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState fs; ///< @todo No cuts?
      VisibleFinalState visfs(fs);

      ZFinder zeeFinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 25*GeV, PID::ELECTRON, 71.0*GeV, 111.0*GeV, 0.1 );
      declare(zeeFinder, "ZeeFinder");

      ZFinder zmumuFinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 25*GeV, PID::MUON, 71.0*GeV, 111.0*GeV, 0.1 );
      declare(zmumuFinder, "ZmumuFinder");

      VetoedFinalState jetConstits(visfs);
      jetConstits.addVetoOnThisFinalState(zeeFinder);
      jetConstits.addVetoOnThisFinalState(zmumuFinder);

      FastJets akt04Jets(jetConstits, FastJets::ANTIKT, 0.4);
      declare(akt04Jets, "AntiKt04Jets");
      
      book(_h_jet_pt_combined,"TMP/jet_pt_combined", refData(1,1,1));
      book(_h_jet_pt_cjet_combined,"TMP/jet_pt_cjet_combined", refData(1,1,1));
      book(_h_jet_pt_bjet_combined,"TMP/jet_pt_bjet_combined", refData(1,1,1));
      
      book(_h_Z_pt_combined,"TMP/Z_pt_combined", refData(2, 1, 1));
      book(_h_Z_pt_cjet_combined,"TMP/Z_pt_cjet_combined", refData(2, 1, 1));
      book(_h_Z_pt_bjet_combined,"TMP/Z_pt_bjet_combined", refData(2, 1, 1));
      
      // book ratio histos      

      book(_h_R_jet_pt_cjet_combined, 1, 1, 1);
      book(_h_R_jet_pt_bjet_combined, 3, 1, 1);
      book(_h_R_jet_pt_cb_combined, 5, 1, 1);

      book(_h_R_Z_pt_cjet_combined, 2, 1, 1);
      book(_h_R_Z_pt_bjet_combined, 4, 1, 1);
      book(_h_R_Z_pt_cb_combined, 6, 1, 1);


    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ZFinder& zeeFS = applyProjection<ZFinder>(event, "ZeeFinder");
      const ZFinder& zmumuFS = applyProjection<ZFinder>(event, "ZmumuFinder");

      const Particles& zees = zeeFS.bosons();
      const Particles& zmumus = zmumuFS.bosons();

      // We did not find exactly one Z. No good.
      if (zees.size() + zmumus.size() != 1) {
        MSG_DEBUG("Did not find exactly one good Z candidate");
        vetoEvent;
      }

      //event identification depending on mass window
      bool ee_event=false;
      bool mm_event=false;
            
      if (zees.size() == 1) { ee_event = true; }
      if (zmumus.size() == 1) { mm_event = true; }
      // const Particles& theLeptons = zees.size() ? zeeFS.constituents() : zmumuFS.constituents();

      // Cluster jets
      // NB. Veto has already been applied on leptons and photons used for dressing
      const FastJets& fj = applyProjection<FastJets>(event, "AntiKt04Jets");
      const Jets& jets = fj.jetsByPt(Cuts::abseta < 2.4 && Cuts::pT > 30*GeV);
      
      // Perform lepton-jet overlap 
      Jets goodjets;
      for (const Jet& j : jets) {
        if ( ee_event ) { if(deltaR(j,zees[0]) > 0.4 )  goodjets.push_back(j);}
        if ( mm_event ) { if(deltaR(j,zmumus[0]) > 0.4 )  goodjets.push_back(j);}
      }

      // We don't care about events with no isolated jets
      if (goodjets.empty()) {
        MSG_DEBUG("No jets in event");
        vetoEvent;
      }

      Jets jc_final;
      Jets jb_final;
            
      //identification of bjets
      int n_ctag = 0; 
      int n_btag = 0; 
       for (const Jet& j : goodjets) {
        if ( j.cTagged() ) { n_ctag = n_ctag + 1; }
        if ( j.bTagged() ) { n_btag = n_btag + 1; }
      }
            
      for (const Jet& j : goodjets) {
        if ( j.cTagged() && n_btag == 0 ) { jc_final.push_back(j); }
        if ( j.bTagged() ) { jb_final.push_back(j); }
      }
      //histogram filling

      if ((ee_event || mm_event) && goodjets.size() > 0) {
        
        FourMomentum j1(goodjets[0].momentum());

        if ( ee_event ) {
           _h_Z_pt_combined->fill(zees[0].pt()); 
           _h_jet_pt_combined -> fill(j1.pt()) ;
           }
        if ( mm_event ) {
           _h_Z_pt_combined->fill(zmumus[0].pt()); 
           _h_jet_pt_combined -> fill(j1.pt()) ;
           }
           
     
        if ( jc_final.size() > 0 ) { 
          FourMomentum c1(jc_final[0].momentum());
          _h_jet_pt_cjet_combined -> fill(c1.pt());
          if ( ee_event ) { 
             _h_Z_pt_cjet_combined->fill(zees[0].pt());
             }
          if ( mm_event ) { 
             _h_Z_pt_cjet_combined->fill(zmumus[0].pt());
             }
          
        }
        if ( jb_final.size() > 0 ) { 
          FourMomentum b1(jb_final[0].momentum());
          _h_jet_pt_bjet_combined -> fill(b1.pt());
          if ( ee_event ) { 
             _h_Z_pt_bjet_combined->fill(zees[0].pt());
             }
          if ( mm_event ) { 
             _h_Z_pt_bjet_combined->fill(zmumus[0].pt());
             }
          
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      MSG_INFO("Cross section = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << crossSection() << " pb");
 
      divide( _h_jet_pt_cjet_combined , _h_jet_pt_combined , _h_R_jet_pt_cjet_combined );
      divide( _h_jet_pt_bjet_combined , _h_jet_pt_combined , _h_R_jet_pt_bjet_combined );
      divide( _h_jet_pt_cjet_combined , _h_jet_pt_bjet_combined , _h_R_jet_pt_cb_combined );
      
      divide ( _h_Z_pt_cjet_combined , _h_Z_pt_combined, _h_R_Z_pt_cjet_combined);
      divide ( _h_Z_pt_bjet_combined , _h_Z_pt_combined, _h_R_Z_pt_bjet_combined);
      divide ( _h_Z_pt_cjet_combined , _h_Z_pt_bjet_combined, _h_R_Z_pt_cb_combined);

    }

    ///@}


    /// @name Histograms
     
     Histo1DPtr _h_jet_pt_combined;
     Histo1DPtr _h_jet_pt_cjet_combined;
     Histo1DPtr _h_jet_pt_bjet_combined;
     Histo1DPtr _h_Z_pt_combined, _h_Z_pt_cjet_combined, _h_Z_pt_bjet_combined;


     Scatter2DPtr _h_R_jet_pt_cjet_combined, _h_R_jet_pt_bjet_combined, _h_R_jet_pt_cb_combined;

 
     Scatter2DPtr _h_R_Z_pt_cjet_combined, _h_R_Z_pt_bjet_combined, _h_R_Z_pt_cb_combined ;

  };


  DECLARE_RIVET_PLUGIN(CMS_2020_I1776758);

}
