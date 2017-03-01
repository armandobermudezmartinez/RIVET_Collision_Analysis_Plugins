#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/AnalysisLoader.hh"

namespace Rivet {

  class MC_TOPMASS_LJETS : public Analysis {
    public:

      /// Minimal constructor
      MC_TOPMASS_LJETS() : Analysis("MC_TOPMASS_LJETS")
    {
    }


      /// @name Analysis methods
      //@{

      /// Set up projections and book histograms
      void init() {

        // Cuts
        Cut particle_cut = (Cuts::abseta < 5.0) and (Cuts::pT > 0.0*MeV);
        Cut lepton_cut   = (Cuts::abseta < 2.1) and (Cuts::pT > 30.*GeV);
        // Jet cut needs to be applied in analyze() function
        // Cut jet_cut = (Cuts::abseta < 2.4) and (Cuts::pT > 30.*GeV);
        
        // Generic final state
        FinalState fs(particle_cut);
        
        // Dressed leptons
        ChargedLeptons charged_leptons(fs);
        PromptFinalState prompt_leptons(charged_leptons);
        IdentifiedFinalState photons(fs);
        photons.acceptIdPair(PID::PHOTON);
        PromptFinalState prompt_photons(photons); //photons not from hadrons
        DressedLeptons dressed_leptons(prompt_photons, prompt_leptons, 0.1, lepton_cut, true, true);
        declare(dressed_leptons, "DressedLeptons");
        
        // Jets
        VetoedFinalState jet_fs(fs);
        jet_fs.addVetoOnThisFinalState(dressed_leptons);
        declare(FastJets(jet_fs, FastJets::ANTIKT, 0.4), "Jets");
        
        // Neutrinos
        IdentifiedFinalState neutrinos(fs);
        neutrinos.acceptNeutrinos();
        PromptFinalState prompt_neutrinos(neutrinos);
        declare(neutrinos, "Neutrinos");
        
        // MET
        declare(MissingMomentum(fs), "MET");

        // Booking of histograms
        _h["weight"] = bookHisto1D("weight", 200, -2., 2.);
        _h["scale"]  = bookHisto1D("scale", logspace(50, 100.0, 1000.0));
        _h["stage"]  = bookHisto1D("stage", 10, 0, 10);
        //
        _h["nbjet"]    = bookHisto1D("nbjet", 10, 0, 10);
        _h["nljet"]    = bookHisto1D("nljet", 10, 0, 10);
        _h["delta_mt"] = bookHisto1D("delta_mt", 100, -100., 100.);
        //
        _h["wlep_mass"]       = bookHisto1D("wlep_mass", 40, 60.4, 100.4);
        _h["wlep_pt"]         = bookHisto1D("wlep_pt", 50, 0., 400.);
        _h["wlep_pt_scaled"]  = bookHisto1D("wlep_pt_scaled", 50, 0., 400.);
        _h["whad_mass"]       = bookHisto1D("whad_mass", 75, 30., 180.);
        _h["whad_pt"]         = bookHisto1D("whad_pt", 50, 0., 400.);
        _h["whad_pt_scaled"]  = bookHisto1D("whad_pt_scaled", 50, 0., 400.);
        _h["tlep_mass"]       = bookHisto1D("tlep_mass", 100, 100., 300.);
        _h["tlep_pt"]         = bookHisto1D("tlep_pt", 50, 0., 400.);
        _h["thad_mass"]       = bookHisto1D("thad_mass", 100, 100., 300.);
        _h["thad_pt"]         = bookHisto1D("thad_pt", 50, 0., 400.);
        _h["ttbar_mass"]      = bookHisto1D("ttbar_mass", 50, 0., 1000.);
        _h["ttbar_pt"]        = bookHisto1D("ttbar_pt", 50, 0., 400.);
      }


      void analyze(const Event& event) {
        const double wmass = 80.4;
        
        const double weight = event.weight();
        _h["weight"]->fill(weight);
        _h["scale"]->fill(event.genEvent()->event_scale(), weight);
        _h["stage"]->fill(0, weight);

        // Get analysis objects from projections
        
        // Exactly 1 lepton for lepton+jets channel
        const std::vector<DressedLepton>& leptons = applyProjection<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
        if (leptons.size() != 1)
          vetoEvent;
        _h["stage"]->fill(1, weight);
        const DressedLepton lepton = leptons[0];

        const Particle neutrino = apply<IdentifiedFinalState>(event, "Neutrinos").particlesByPt()[0];
        if (not(neutrino.abspid() == lepton.abspid()+1 and neutrino.pid()*lepton.pid()<0))
          vetoEvent;
        _h["stage"]->fill(2, weight);
        
        FourMomentum wlep(lepton.momentum() + neutrino.momentum());
        _h["wlep_mass"]->fill(wlep.mass(), weight);
        if (not(inRange(wlep.mass(), wmass-5., wmass+5.)))
          vetoEvent;
        _h["stage"]->fill(3, weight);
        
        // Select events with at least 4 jets (2 b-tagged, 2 untagged)
        Cut jet_cut = (Cuts::abseta < 2.4) and (Cuts::pT > 30.*GeV);
        const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(jet_cut);
        Jets bjets, ljets;
        for (const Jet& jet : jets) {
          if (jet.bTagged()) bjets.push_back(jet);
          else               ljets.push_back(jet);
        }
        _h["nbjet"]->fill(bjets.size(), weight);
        _h["nljet"]->fill(ljets.size(), weight);
        if (bjets.size() < 2 or ljets.size() < 2)
          vetoEvent;
        _h["stage"]->fill(4, weight);
        
        // Reconstruct top quarks
        // Constrain W candidate masses to 80.4 GeV
        // Top candidates should have similar masses
        FourMomentum whad, thad, tlep;
        double kmin = numeric_limits<double>::max();
        for (const Jet& ljet1 : ljets) {
          for (const Jet& ljet2 : ljets) {
            if (&ljet1 == &ljet2 or ljet1.pt() < ljet2.pt())
              continue;
            FourMomentum whad_cand(ljet1.momentum() + ljet2.momentum());
            for (const Jet& bjethad : bjets) {
              FourMomentum thad_cand(whad_cand * wmass/whad_cand.mass() + bjethad.momentum());
              for (const Jet& bjetlep : bjets) {
                if (&bjethad == &bjetlep)
                  continue;
                FourMomentum tlep_cand(wlep * wmass/wlep.mass() + bjetlep.momentum());

                double K = pow(whad_cand.mass() - wmass, 2) + pow(thad_cand.mass() - tlep_cand.mass(), 2);
                if(K < kmin) {
                  kmin = K;
                  whad = whad_cand;
                  thad = thad_cand;
                  tlep = tlep_cand;
                }
              }
            }
          }
        }
        FourMomentum ttbar = thad + tlep;
        
        _h["whad_mass"]->fill(whad.mass(), weight);
        double delta_mt = thad.mass() - tlep.mass();
        _h["delta_mt"] ->fill(delta_mt, weight);
        if (not(inRange(whad.mass(), wmass-10., wmass+10.) and abs(delta_mt)<20.))
          vetoEvent;
        _h["stage"]->fill(5, weight);
        
        // Fill plots
        _h["wlep_pt"]->fill(wlep.pt(), weight);
        _h["wlep_pt_scaled"]->fill(wlep.pt() * wmass/wlep.mass(), weight);
        _h["whad_pt"]->fill(whad.pt(), weight);
        _h["whad_pt_scaled"]->fill(whad.pt() * wmass/whad.mass(), weight);
        _h["tlep_mass"]->fill(tlep.mass(), weight);
        _h["tlep_pt"]->fill(tlep.pt(), weight);
        _h["thad_mass"]->fill(thad.mass(), weight);
        _h["thad_pt"]->fill(thad.pt(), weight);
        _h["ttbar_mass"]->fill(ttbar.mass(), weight);
        _h["ttbar_pt"]->fill(ttbar.pt(), weight);

      }


      void finalize() {
      }

      //@}


    private:

      // @name Histogram data members
      //@{

      map<string, Histo1DPtr> _h;

      //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TOPMASS_LJETS);

}
