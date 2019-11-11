// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include <iostream>

namespace Rivet {

  class CMS_2019_I1705068 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2019_I1705068);


    void init() {

      FinalState fs;
      WFinder wfinder_mu(fs, Cuts::abseta < 2.4 && Cuts::pT > 0*GeV, PID::MUON, 0*GeV, 1000000*GeV, 0*GeV, 0.1, WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);
      addProjection(wfinder_mu, "WFinder_mu");
      
      UnstableParticles dst(Cuts::pT > 5*GeV && Cuts::abseta < 2.4); 
      addProjection(dst, "Dstar");
      
      // Particle-Level Histograms form the paper: 
      _hist_WplusMinus_MuAbseta = bookHisto1D("d04-x01-y01");
      _hist_Wplus_MuAbseta = bookHisto1D("d05-x01-y01");
      _hist_Wminus_MuAbseta = bookHisto1D("d06-x01-y01");   
           
    }

    void analyze(const Event& event) {
        
        const double weight = event.weight();        
        const WFinder& wfinder_mu = applyProjection<WFinder>(event, "WFinder_mu"); 
        
        if (wfinder_mu.bosons().size() != 1) vetoEvent;  
        // No Missing Energy or MT Cut at generator level:
        const FourMomentum& lepton0 = wfinder_mu.constituentLeptons()[0].momentum();
         
        double pt0 = lepton0.pT();
        double eta0 = fabs( lepton0.eta() );
        if ( (eta0 > 2.4) || (pt0 < 26.0*GeV) ) vetoEvent;
        
        int muID = wfinder_mu.constituentLeptons()[0].pid(); 

        
        // D* selection: 
        // OS = W boson and D* Meson have Opposite (charge) Signs
        // SS = W Boson and D* Meson have Same (charge) Signs
        // Associated W+c only has OS contributions, 
        // W+ccbar (ccbar from gluon splitting) has equal probability to be OS or SS
        // OS-SS to remove the gluon splitting background
        
        const UnstableParticles& dst = applyProjection<UnstableFinalState>(event, "Dstar");
        for(auto p: dst.particles()) {
          if(muID == -13 && p.pid() == -413){ // OS
            _hist_Wplus_MuAbseta->fill(eta0, weight); 
            _hist_WplusMinus_MuAbseta->fill(eta0, weight); 
          }
          else if(muID == 13 && p.pid() == 413){ // OS
            _hist_Wminus_MuAbseta->fill(eta0, weight);
            _hist_WplusMinus_MuAbseta->fill(eta0, weight);
          }
          else if (muID == -13 && p.pid() == 413) { // SS
            _hist_Wplus_MuAbseta->fill(eta0, weight*-1); 
            _hist_WplusMinus_MuAbseta->fill(eta0, weight*-1); 
          }
          else if (muID == 13 && p.pid() == -413) { // SS
            _hist_Wminus_MuAbseta->fill(eta0, weight*-1); 
            _hist_WplusMinus_MuAbseta->fill(eta0, weight*-1); 
          }
        }
        
        
    }


    void finalize() {

      scale(_hist_Wplus_MuAbseta, crossSection()/picobarn/sumOfWeights());
      scale(_hist_Wminus_MuAbseta, crossSection()/picobarn/sumOfWeights());
      scale(_hist_WplusMinus_MuAbseta, crossSection()/picobarn/sumOfWeights());
    }


    Histo1DPtr _hist_Wplus_MuAbseta;
    Histo1DPtr _hist_Wminus_MuAbseta;
    Histo1DPtr _hist_WplusMinus_MuAbseta;
    
    
  };

  DECLARE_RIVET_PLUGIN(CMS_2019_I1705068);

}
