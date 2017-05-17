// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/InitialQuarks.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include <cmath>

/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {
  
  
  class CMS_2016_I1486238 : public Analysis {
  public:
    
    /// @name Constructors etc.
    //@{
    
    /// Constructor
    CMS_2016_I1486238()
      : Analysis("CMS_2016_I1486238")
    {
      /// @todo Set whether your finalize method needs the generator cross section
    }
    
    //@}
    
    
  public:
    
    /// @name Analysis methods
    //@{
    
    /// Book histograms and initialise projections before the run
    void init() {
      
      /// @todo Initialise and register projections here
      
      FinalState fs;
      FastJets akt(fs, FastJets::ANTIKT, 0.5);
      addProjection(akt, "antikT");
      
      /// @todo Book histograms here, e.g.:
      
      _h_Deltaphi_newway = bookHisto1D(1,1,1);
      _h_deltaphiafterlight = bookHisto1D(9,1,1);
      _h_SumPLight = bookHisto1D(5,1,1);
      
      _h_LeadingBJetpt = bookHisto1D(11,1,1);
      _h_SubleadingBJetpt = bookHisto1D(15,1,1);
      _h_LeadingLightJetpt = bookHisto1D(13,1,1);
      _h_SubleadingLightJetpt = bookHisto1D(17,1,1);

      _h_LeadingBJeteta = bookHisto1D(10,1,1);
      _h_SubleadingBJeteta = bookHisto1D(14,1,1);
      _h_LeadingLightJeteta = bookHisto1D(12,1,1);
      _h_SubleadingLightJeteta = bookHisto1D(16,1,1);
      
    }
    
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      
      const Jets& jets = applyProjection<JetAlg>(event, "antikT").jetsByPt( Cuts::absrap < 4.7 && Cuts::pT > 20.*GeV );
      
      //Initial quarks
      
      if(jets.size()<4) vetoEvent;

      Jets bjets;
      Jets ljets;
      
      foreach (const Jet& j, jets) {
        
        if(j.momentum().pT()>20) {
          
          bool btag=false;
          
          foreach (const GenParticle* p, particles(event.genEvent())) {
            
            const PdgId pid = p->pdg_id();
            
            if (abs(pid) == 5) { 
              
              if (deltaR(j.momentum(), (FourMomentum) p->momentum()) < 0.3) {
                btag=true;
              }
            }
          }
          
          if(btag && fabs(j.momentum().eta()) < 2.4) bjets.push_back(j);
          else ljets.push_back(j);
        }
      }
      
      if (bjets.size() >= 2 && ljets.size() >= 2){
        
        _h_LeadingBJetpt->fill(bjets[0].pT(),weight);
        _h_SubleadingBJetpt->fill(bjets[1].pT(),weight);
        _h_LeadingLightJetpt->fill(ljets[0].pT(),weight);
        _h_SubleadingLightJetpt->fill(ljets[1].pT(),weight);
        
        _h_LeadingBJeteta->fill(bjets[0].eta(),weight);
        _h_SubleadingBJeteta->fill(bjets[1].eta(),weight);
        _h_LeadingLightJeteta->fill(ljets[0].eta(),weight);
        _h_SubleadingLightJeteta->fill(ljets[1].eta(),weight);
        
        double lightdphi=deltaPhi(ljets[0].phi(),ljets[1].phi());
        _h_deltaphiafterlight->fill(lightdphi,weight);
        
        const double vecsumlightjets=sqrt(pow(ljets[0].px()+ljets[1].px(),2)+pow(ljets[0].py()+ljets[1].py(),2));
        
        const double term2=vecsumlightjets/(sqrt(ljets[0].px()*ljets[0].px()+ljets[0].py()*ljets[0].py())+sqrt(ljets[1].px()*ljets[1].px()+ljets[1].py()*ljets[1].py()));
        
        _h_SumPLight->fill(term2,weight);
        
        const double pxBsyst2=bjets[0].px()+bjets[1].px();
        const double pyBsyst2=bjets[0].py()+bjets[1].py();
        const double pxJetssyst2=ljets[0].px()+ljets[1].px();
        const double pyJetssyst2=ljets[0].py()+ljets[1].py();
        const double modulusB2=sqrt(pow(pxBsyst2,2)+pow(pyBsyst2,2));
        const double modulusJets2=sqrt(pow(pxJetssyst2,2)+pow(pyJetssyst2,2));
        const double cosphiBsyst2=pxBsyst2/modulusB2;
        const double cosphiJetssyst2=pxJetssyst2/modulusJets2;
        double phiBsyst2=0;
        double phiJetssyst2=0;
        if(pyBsyst2>0) {phiBsyst2=acos(cosphiBsyst2);}
        if(pyBsyst2<0) {phiBsyst2=-acos(cosphiBsyst2);}
        if(pyJetssyst2>0) {phiJetssyst2=acos(cosphiJetssyst2);}
        if(pyJetssyst2<0) {phiJetssyst2=-acos(cosphiJetssyst2);}
        
        const double Dphi2=deltaPhi(phiBsyst2,phiJetssyst2);

        _h_Deltaphi_newway->fill(Dphi2,weight);
        
      }   
      
      
    }
    
    /// Normalise histograms etc., after the run
    void finalize() {
      
      /// @todo Normalise, scale and otherwise manipulate histograms here
      
      // scale(_h_YYYY, crossSection()/sumOfWeights()); # norm to cross section
      // normalize(_h_YYYY); # normalize to unity
      //sumOfWeights() returns the weights of the observations, that is the number of generated events
      
      double invlumi = crossSection()/picobarn/sumOfWeights(); //norm to cross section
      
      normalize(_h_SumPLight);
      normalize(_h_deltaphiafterlight);
      normalize(_h_Deltaphi_newway);
      
      scale(_h_LeadingLightJetpt, invlumi);
      scale(_h_SubleadingLightJetpt, invlumi);
      scale(_h_LeadingBJetpt, invlumi);
      scale(_h_SubleadingBJetpt, invlumi);      
      
      scale(_h_LeadingLightJeteta, invlumi);
      scale(_h_SubleadingLightJeteta, invlumi);
      scale(_h_LeadingBJeteta, invlumi);
      scale(_h_SubleadingBJeteta, invlumi);  
      
    }
    
    //@}
    

  private:
    
    // Data members like post-cuts event weight counters go here
    
    
  private:
    
    /// @name Histograms
    //@{
    
    Histo1DPtr _h_deltaphiafterlight;
    Histo1DPtr _h_Deltaphi_newway;
    Histo1DPtr _h_SumPLight;
    
    Histo1DPtr _h_LeadingBJetpt;
    Histo1DPtr _h_SubleadingBJetpt;
    Histo1DPtr _h_LeadingLightJetpt;
    Histo1DPtr _h_SubleadingLightJetpt;
    
    Histo1DPtr _h_LeadingBJeteta;
    Histo1DPtr _h_SubleadingBJeteta;
    Histo1DPtr _h_LeadingLightJeteta;
    Histo1DPtr _h_SubleadingLightJeteta;
    
  };

  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CMS_2016_I1486238> plugin_CMS_2016_I1486238;
  

}
