#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Tools/ParticleName.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {
  
  class CMS_2016_I1454211 : public Analysis {
  public:
    
    // Minimal constructor
    CMS_2016_I1454211() : Analysis("CMS_2016_I1454211") {
    }
    
    // Set up projections and book histograms
    void init() {
      
      // Complete final state
      FinalState fs(-MAXDOUBLE, MAXDOUBLE, 0*GeV);
      
      // Partonic tops
      declare(PartonicTops(PartonicTops::ELECTRON, false), "ElectronPartonTops");
      declare(PartonicTops(PartonicTops::MUON, false),     "MuonPartonTops");
      declare(PartonicTops(PartonicTops::HADRONIC),        "HadronicPartonTops");
      
      // Projection for electrons and muons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);
      
      Cut leptonCuts = Cuts::pt > 45*GeV && Cuts::abseta < 2.1;
      
      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      DressedLeptons dressed_electrons(photons, electrons, 0.1, leptonCuts, true, false);
      addProjection(dressed_electrons, "DressedElectrons");
      
      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      DressedLeptons dressed_muons(photons, muons, 0.1, leptonCuts, true, false);
      addProjection(dressed_muons, "DressedMuons");
      
      // Projection for jets
      VetoedFinalState fs_jets(FinalState(-MAXDOUBLE, MAXDOUBLE, 0*GeV));
      fs_jets.addVetoOnThisFinalState(dressed_muons);
      fs_jets.addVetoOnThisFinalState(dressed_electrons);
      fs_jets.vetoNeutrinos();
      addProjection(FastJets(fs_jets, FastJets::ANTIKT, 0.5), "ak5jets");
      addProjection(FastJets(fs_jets, FastJets::CAM, 0.8), "ca8jets");
      
      _hEl_topPt_parton          = bookHisto1D("d01-x01-y01"); // dsigma/dpt(top quark), el ch
      _hEl_topPt_particle        = bookHisto1D("d02-x01-y01"); // dsigma/dpt(top jet), el ch
      _hEl_topY_parton           = bookHisto1D("d03-x01-y01"); // dsigma/dy(top quark), el ch 
      _hEl_topY_particle         = bookHisto1D("d04-x01-y01"); // dsigma/dy(top jet), el ch
      _hMu_topPt_parton          = bookHisto1D("d05-x01-y01"); // dsigma/dpt(top quark), mu ch
      _hMu_topPt_particle        = bookHisto1D("d06-x01-y01"); // dsigma/dpt(top jet), mu ch
      _hMu_topY_parton           = bookHisto1D("d07-x01-y01"); // dsigma/dy(top quark), mu ch
      _hMu_topY_particle         = bookHisto1D("d08-x01-y01"); // dsigma/dy(top jet), mu ch
      _hComb_topPt_parton        = bookHisto1D("d09-x01-y01"); // dsigma/dpt(top quark), comb ch
      _hComb_topPt_particle      = bookHisto1D("d10-x01-y01"); // dsigma/dpt(top jet), comb ch
      _hComb_topY_parton         = bookHisto1D("d11-x01-y01"); // dsigma/dy(top quark), comb ch
      _hComb_topY_particle       = bookHisto1D("d12-x01-y01"); // dsigma/dy(top jet), comb ch
      
      _hEl_topPt_parton_norm     = bookHisto1D("d13-x01-y01"); // 1/sigma dsigma/dpt(top quark), el ch
      _hEl_topPt_particle_norm   = bookHisto1D("d14-x01-y01"); // 1/sigma dsigma/dpt(top jet), el ch
      _hEl_topY_parton_norm      = bookHisto1D("d15-x01-y01"); // 1/sigma dsigma/dy(top quark), el ch
      _hEl_topY_particle_norm    = bookHisto1D("d16-x01-y01"); // 1/sigma dsigma/dy(top jet), el ch
      _hMu_topPt_parton_norm     = bookHisto1D("d17-x01-y01"); // 1/sigma dsigma/dpt(top quark), mu ch
      _hMu_topPt_particle_norm   = bookHisto1D("d18-x01-y01"); // 1/sigma dsigma/dpt(top jet), mu ch
      _hMu_topY_parton_norm      = bookHisto1D("d19-x01-y01"); // 1/sigma dsigma/dy(top quark), mu ch
      _hMu_topY_particle_norm    = bookHisto1D("d20-x01-y01"); // 1/sigma dsigma/dy(top jet), mu ch
      _hComb_topPt_parton_norm   = bookHisto1D("d21-x01-y01"); // 1/sigma dsigma/dpt(top quark), comb ch
      _hComb_topPt_particle_norm = bookHisto1D("d22-x01-y01"); // 1/sigma dsigma/dpt(top jet), comb ch
      _hComb_topY_parton_norm    = bookHisto1D("d23-x01-y01"); // 1/sigma dsigma/dy(top quark), comb ch
      _hComb_topY_particle_norm  = bookHisto1D("d24-x01-y01"); // 1/sigma dsigma/dy(top jet), comb ch
      
      _hMu_cutflow = bookHisto1D("mu_cutflow",7,-0.5,6.5);
      _hEl_cutflow = bookHisto1D("el_cutflow",7,-0.5,6.5);
      
      nMu = 0.;
      nEl = 0.;
      nPassParton_mu = 0.;
      nPassParton_el = 0.;
      nPassParticle_mu = 0.;
      nPassParticle_el = 0.;
    }
    
    
    // per event analysis
    void analyze(const Event& event) {
      
      _hMu_cutflow->fill(0.); // total events
      _hEl_cutflow->fill(0.);
      
      const double weight = event.weight();
      
      // Do parton-level selection and channel determination
      int partonCh = 0; //0 non-semi-lep, 1 muon, 2 electron
      const Particles muonpartontops = apply<ParticleFinder>(event, "MuonPartonTops").particlesByPt();
      const Particles electronpartontops = apply<ParticleFinder>(event, "ElectronPartonTops").particlesByPt();
      if (electronpartontops.size() == 0 && muonpartontops.size() == 1) partonCh = 1;
      else if (electronpartontops.size() == 1 && muonpartontops.size() == 0) partonCh = 2;
      else vetoEvent;
      const Particles hadronicpartontops = apply<ParticleFinder>(event, "HadronicPartonTops").particlesByPt();
      if (hadronicpartontops.size() != 1) vetoEvent;
      
      if (partonCh == 1) _hMu_cutflow->fill(1.); // muon at parton level
      if (partonCh == 2) _hEl_cutflow->fill(1.); // electron at parton level
      
      // Get hadronic parton-level top
      const FourMomentum& partonTopP4 = hadronicpartontops.at(0).momentum();
      
      // Do particle-level selection and channel determination
      const DressedLeptons& dressed_electrons = applyProjection<DressedLeptons>(event, "DressedElectrons");
      const DressedLeptons& dressed_muons = applyProjection<DressedLeptons>(event, "DressedMuons");
      
      bool passParticleLep = false;
      bool passParticleTop = false;
      FourMomentum lepton;
      FourMomentum particleTopP4;
      
      if (partonCh == 1 && dressed_muons.dressedLeptons().size() == 1 && dressed_electrons.dressedLeptons().size() == 0){
        passParticleLep = true;
        _hMu_cutflow->fill(3.); //muon at particle level
        lepton = dressed_muons.dressedLeptons()[0].momentum();
      }
      if (partonCh == 2 && dressed_muons.dressedLeptons().size() == 0 && dressed_electrons.dressedLeptons().size() == 1){
        passParticleLep = true;
        _hEl_cutflow->fill(3.); //electron at particle level
        lepton = dressed_electrons.dressedLeptons()[0].momentum();
      }
      
      if (passParticleLep){
        
        // Jet cuts
        Cut jetCuts = Cuts::pt > 30*GeV && Cuts::abseta < 2.4;
        Jets genBjets;
        Jets genTjets;
        int nGenBjets = 0;
        int nGenTjets = 0;
        
        const FastJets& AK5jets = applyProjection<FastJets>(event, "ak5jets");
        
        foreach (const Jet& jet, AK5jets.jetsByPt(jetCuts)) {
          if (deltaR(jet.momentum(),lepton) > 3.1415 / 2.0) continue;
          if (deltaR(jet.momentum(),lepton) < 0.1) continue;
          genBjets.push_back(jet);
          nGenBjets += 1;
        }
        
        const FastJets& CA8jets = applyProjection<FastJets>(event, "ca8jets");
        
        foreach (const Jet& jet, CA8jets.jetsByPt(jetCuts)) {
          if (deltaR(jet.momentum(), lepton) < 3.1415 / 2.0) continue;
          if (jet.momentum().mass() < 140.) continue;
          if (jet.momentum().mass() > 250.) continue;
          genTjets.push_back(jet);
          nGenTjets += 1;
        }
        
        if (nGenBjets >=1){
          if (partonCh == 1) _hMu_cutflow->fill(4.); // muon at parton level
          if (partonCh == 2) _hEl_cutflow->fill(4.); // electron at parton level
          
          if (nGenTjets >= 1){
            passParticleTop = true;
            if (partonCh == 1) _hMu_cutflow->fill(5.); // muon at parton level
            if (partonCh == 2) _hEl_cutflow->fill(5.); // electron at parton level 
            
            particleTopP4 = genTjets[0].momentum();
          }
        }
      }
      
      if (partonCh == 1){
        nMu += 1;
        _hMu_topPt_parton->fill(partonTopP4.pT(), weight);
        _hMu_topPt_parton_norm->fill(partonTopP4.pT(), weight);
        _hComb_topPt_parton->fill(partonTopP4.pT(), weight);
        _hComb_topPt_parton_norm->fill(partonTopP4.pT(), weight);
        
        if (partonTopP4.pT() >= 400.){
          nPassParton_mu += 1;
          _hMu_cutflow->fill(2.);
          _hMu_topY_parton->fill(partonTopP4.rapidity(), weight);
          _hMu_topY_parton_norm->fill(partonTopP4.rapidity(), weight);
          _hComb_topY_parton->fill(partonTopP4.rapidity(), weight);
          _hComb_topY_parton_norm->fill(partonTopP4.rapidity(), weight);
        }
        
        if (passParticleTop){
          _hMu_topPt_particle->fill(particleTopP4.pT(), weight);
          _hMu_topPt_particle_norm->fill(particleTopP4.pT(), weight);
          _hComb_topPt_particle->fill(particleTopP4.pT(), weight);
          _hComb_topPt_particle_norm->fill(particleTopP4.pT(), weight);
          
          if (particleTopP4.pT() >= 400.){
            nPassParticle_mu += 1;
            _hMu_cutflow->fill(6.);
            _hMu_topY_particle->fill(particleTopP4.rapidity(), weight);
            _hMu_topY_particle_norm->fill(particleTopP4.rapidity(), weight);
            _hComb_topY_particle->fill(particleTopP4.rapidity(), weight);
            _hComb_topY_particle_norm->fill(particleTopP4.rapidity(), weight);
          }
        }
      }
      
      if (partonCh == 2){
        nEl += 1;
        _hEl_topPt_parton->fill(partonTopP4.pT(), weight);
        _hEl_topPt_parton_norm->fill(partonTopP4.pT(), weight);
        _hComb_topPt_parton->fill(partonTopP4.pT(), weight);
        _hComb_topPt_parton_norm->fill(partonTopP4.pT(), weight);
        
        if (partonTopP4.pT() >= 400.){
          nPassParton_el += 1;
          _hEl_cutflow->fill(2.);
          _hEl_topY_parton->fill(partonTopP4.rapidity(), weight);
          _hEl_topY_parton_norm->fill(partonTopP4.rapidity(), weight);
          _hComb_topY_parton->fill(partonTopP4.rapidity(), weight);
          _hComb_topY_parton_norm->fill(partonTopP4.rapidity(), weight);
        }
        
        if (passParticleTop){
          _hEl_topPt_particle->fill(particleTopP4.pT(), weight);
          _hEl_topPt_particle_norm->fill(particleTopP4.pT(), weight);
          _hComb_topPt_particle->fill(particleTopP4.pT(), weight);
          _hComb_topPt_particle_norm->fill(particleTopP4.pT(), weight);
          
          if (particleTopP4.pT() >= 400.){
            nPassParticle_el += 1;
            _hEl_cutflow->fill(6.);
            _hEl_topY_particle->fill(particleTopP4.rapidity(), weight);
            _hEl_topY_particle_norm->fill(particleTopP4.rapidity(), weight);
            _hComb_topY_particle->fill(particleTopP4.rapidity(), weight);
            _hComb_topY_particle_norm->fill(particleTopP4.rapidity(), weight);
          }
        }
      }
    }
    
    void finalize() {
      
      double xs_mu_parton = 252.89 * 1000. * nPassParton_mu / nMu;
      double xs_mu_particle = 252.89 * 1000. * nPassParticle_mu / nMu;
      double xs_el_parton = 252.89 * 1000. * nPassParton_el / nEl;
      double xs_el_particle = 252.89 * 1000. * nPassParticle_el / nEl;
      double xs_comb_parton = 252.89 * 1000. * (nPassParton_el + nPassParton_mu) / (nEl + nMu);
      double xs_comb_particle = 252.89 * 1000. * (nPassParticle_el + nPassParticle_mu) / (nEl + nMu);
      cout << "Inclusive xsec (pt>400), mu channel, parton-level: " << xs_mu_parton << endl;
      cout << "Inclusive xsec (pt>400), mu channel, particle-level: " << xs_mu_particle << endl;
      cout << "Inclusive xsec (pt>400), el channel, parton-level: " << xs_el_parton << endl;
      cout << "Inclusive xsec (pt>400), el channel, particle-level: " << xs_el_particle << endl;
      cout << "Inclusive xsec (pt>400), comb channel, parton-level: " << xs_comb_parton << endl;
      cout << "Inclusive xsec (pt>400), comb channel, particle-level: " << xs_comb_particle << endl;

      normalize(_hMu_topPt_parton,xs_mu_parton,false);
      normalize(_hMu_topPt_particle,xs_mu_particle,false);
      normalize(_hMu_topY_parton,xs_mu_parton,false);
      normalize(_hMu_topY_particle,xs_mu_particle,false);
      normalize(_hEl_topPt_parton,xs_el_parton,false);
      normalize(_hEl_topPt_particle,xs_el_particle,false);
      normalize(_hEl_topY_parton,xs_el_parton,false);
      normalize(_hEl_topY_particle,xs_el_particle,false);
      normalize(_hComb_topPt_parton,xs_comb_parton,false);
      normalize(_hComb_topPt_particle,xs_comb_particle,false);
      normalize(_hComb_topY_parton,xs_comb_parton,false);
      normalize(_hComb_topY_particle,xs_comb_particle,false);
      
      normalize(_hMu_topPt_parton_norm,1.0,false);
      normalize(_hMu_topPt_particle_norm,1.0,false);
      normalize(_hMu_topY_parton_norm,1.0,false);
      normalize(_hMu_topY_particle_norm,1.0,false);
      normalize(_hEl_topPt_parton_norm,1.0,false);
      normalize(_hEl_topPt_particle_norm,1.0,false);
      normalize(_hEl_topY_parton_norm,1.0,false);
      normalize(_hEl_topY_particle_norm,1.0,false);
      normalize(_hComb_topPt_parton_norm,1.0,false);
      normalize(_hComb_topPt_particle_norm,1.0,false);
      normalize(_hComb_topY_parton_norm,1.0,false);
      normalize(_hComb_topY_particle_norm,1.0,false);
    }
    
  private:
    
    Histo1DPtr _hMu_topPt_parton;
    Histo1DPtr _hMu_topPt_particle;
    Histo1DPtr _hMu_topY_parton;
    Histo1DPtr _hMu_topY_particle;
    Histo1DPtr _hEl_topPt_parton;
    Histo1DPtr _hEl_topPt_particle;
    Histo1DPtr _hEl_topY_parton;
    Histo1DPtr _hEl_topY_particle;
    Histo1DPtr _hComb_topPt_parton;
    Histo1DPtr _hComb_topPt_particle;
    Histo1DPtr _hComb_topY_parton;
    Histo1DPtr _hComb_topY_particle;
    Histo1DPtr _hMu_topPt_parton_norm;
    Histo1DPtr _hMu_topPt_particle_norm;
    Histo1DPtr _hMu_topY_parton_norm;
    Histo1DPtr _hMu_topY_particle_norm;
    Histo1DPtr _hEl_topPt_parton_norm;
    Histo1DPtr _hEl_topPt_particle_norm;
    Histo1DPtr _hEl_topY_parton_norm;
    Histo1DPtr _hEl_topY_particle_norm;
    Histo1DPtr _hComb_topPt_parton_norm;
    Histo1DPtr _hComb_topPt_particle_norm;
    Histo1DPtr _hComb_topY_parton_norm;
    Histo1DPtr _hComb_topY_particle_norm;
    Histo1DPtr _hMu_cutflow;
    Histo1DPtr _hEl_cutflow;
    
    int nEvents;
    int nMu;
    int nEl;
    int nPassParton_mu;
    int nPassParton_el;
    int nPassParticle_mu;
    int nPassParticle_el;
  };
  
  // The hook for the plugin system                                                                                                                                     
  DECLARE_RIVET_PLUGIN(CMS_2016_I1454211);

}
