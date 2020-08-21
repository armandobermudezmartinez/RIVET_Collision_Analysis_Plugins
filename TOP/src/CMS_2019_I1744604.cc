// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PartonicTops.hh"

namespace Rivet {
  class CMS_2019_I1744604 : public Analysis {
  private:

    static std::vector<double> EquationSolve(double a, double b, double c, double d)
    {
      std::vector<double> result;
      
      std::complex<double> x1;
      std::complex<double> x2;
      std::complex<double> x3;

      if (a != 0) {
        
        double q = (3*a*c-b*b)/(9*a*a);
        double r = (9*a*b*c - 27*a*a*d - 2*b*b*b)/(54*a*a*a);
        double Delta = q*q*q + r*r;

        std::complex<double> s;
        std::complex<double> t;

        double rho=0;
        double theta=0;
        
        if( Delta<=0) {
          rho = sqrt(-(q*q*q));

          theta = acos(r/rho);

          s = std::polar<double>(sqrt(-q),theta/3.0);
          t = std::polar<double>(sqrt(-q),-theta/3.0);
        }
        
        if(Delta>0) {
          s = std::complex<double>(cbrt(r+sqrt(Delta)),0);
          t = std::complex<double>(cbrt(r-sqrt(Delta)),0);
        }
      
        std::complex<double> i(0,1.0);
        
        
        x1 = s+t+std::complex<double>(-b/(3.0*a),0);
        x2 = (s+t)*std::complex<double>(-0.5,0)-std::complex<double>(b/(3.0*a),0)+(s-t)*i*std::complex<double>(sqrt(3)/2.0,0);
        x3 = (s+t)*std::complex<double>(-0.5,0)-std::complex<double>(b/(3.0*a),0)-(s-t)*i*std::complex<double>(sqrt(3)/2.0,0);

        if (fabs(x1.imag())<0.0001) result.push_back(x1.real());
        if (fabs(x2.imag())<0.0001) result.push_back(x2.real());
        if (fabs(x3.imag())<0.0001) result.push_back(x3.real());

        return result;
        
      } else {
      
        return result;
      }
      
      return result;
    }

    static std::pair<FourMomentum,FourMomentum> NuMomentum(
      double leptonPx, double leptonPy, double leptonPz, double leptonPt, 
      double leptonE, double metPx, double metPy
    )
    {
      double  mW = 80.399;

      FourMomentum result(0,0,0,0);
      FourMomentum result2(0,0,0,0);


      double MisET2 = (metPx * metPx + metPy * metPy);
      double mu = (mW * mW) / 2 + metPx * leptonPx + metPy * leptonPy;
      double a  = (mu * leptonPz) / (leptonE * leptonE - leptonPz * leptonPz);
      double a2 = std::pow(a, 2);
      double b  = (std::pow(leptonE, 2.) * (MisET2) - std::pow(mu, 2.)) / (std::pow(leptonE, 2) - std::pow(leptonPz, 2));
      double pz1(0), pz2(0), pznu(0), pznu2(0);

      FourMomentum p4W_rec;
      FourMomentum p4b_rec;
      FourMomentum p4Top_rec;
      FourMomentum p4lep_rec;

      p4lep_rec.setXYZE(leptonPx, leptonPy, leptonPz, leptonE);

      FourMomentum p40_rec(0, 0, 0, 0);

      if (a2 - b > 0 )
      {

        double root = sqrt(a2 - b);
        pz1 = a + root;
        pz2 = a - root;

        pznu = pz1;
        pznu2 = pz2;
        
        if (fabs(pz1) > fabs(pz2)) {
          pznu = pz2;
          pznu2 = pz1;
        }


        double Enu = sqrt(MisET2 + pznu * pznu);
        double Enu2 = sqrt(MisET2 + pznu2 * pznu2);

        result.setXYZE(metPx, metPy, pznu, Enu);
        result2.setXYZE(metPx, metPy, pznu2, Enu2);
   

      } else {


        double ptlep = leptonPt, pxlep = leptonPx, pylep = leptonPy, metpx = metPx, metpy = metPy;

        double EquationA = 1;
        double EquationB = -3 * pylep * mW / (ptlep);
        double EquationC = mW * mW * (2 * pylep * pylep) / (ptlep * ptlep) + mW * mW - 4 * pxlep * pxlep * pxlep * metpx / (ptlep * ptlep) - 4 * pxlep * pxlep * pylep * metpy / (ptlep * ptlep);
        double EquationD = 4 * pxlep * pxlep * mW * metpy / (ptlep) - pylep * mW * mW * mW / ptlep;

        std::vector<double> solutions = EquationSolve(EquationA, EquationB, EquationC, EquationD);

        std::vector<double> solutions2 = EquationSolve(EquationA, -EquationB, EquationC, -EquationD);


        double deltaMin = 14000 * 14000;
        double zeroValue = -mW * mW / (4 * pxlep);
        double minPx = 0;
        double minPy = 0;

        for ( size_t i = 0; i < solutions.size(); ++i) {
          if (solutions[i] < 0 ) continue;
          double p_x = (solutions[i] * solutions[i] - mW * mW) / (4 * pxlep);
          double p_y = ( mW * mW * pylep + 2 * pxlep * pylep * p_x - mW * ptlep * solutions[i]) / (2 * pxlep * pxlep);
          double Delta2 = (p_x - metpx) * (p_x - metpx) + (p_y - metpy) * (p_y - metpy);


          if (Delta2 < deltaMin && Delta2 > 0) {
            deltaMin = Delta2;
            minPx = p_x;
            minPy = p_y;
          }
            
        }

        for ( size_t i = 0; i < solutions2.size(); ++i) {
          if (solutions2[i] < 0 ) continue;
          double p_x = (solutions2[i] * solutions2[i] - mW * mW) / (4 * pxlep);
          double p_y = ( mW * mW * pylep + 2 * pxlep * pylep * p_x + mW * ptlep * solutions2[i]) / (2 * pxlep * pxlep);
          double Delta2 = (p_x - metpx) * (p_x - metpx) + (p_y - metpy) * (p_y - metpy);
          if (Delta2 < deltaMin && Delta2 > 0) {
            deltaMin = Delta2;
            minPx = p_x;
            minPy = p_y;
          }
        }

        double pyZeroValue = ( mW * mW * pxlep + 2 * pxlep * pylep * zeroValue);
        double delta2ZeroValue = (zeroValue - metpx) * (zeroValue - metpx) + (pyZeroValue - metpy) * (pyZeroValue - metpy);

        if (deltaMin == 14000 * 14000)return std::make_pair(result,result2);

        if (delta2ZeroValue < deltaMin) {
          deltaMin = delta2ZeroValue;
          minPx = zeroValue;
          minPy = pyZeroValue;
        }


        double mu_Minimum = (mW * mW) / 2 + minPx * pxlep + minPy * pylep;
        double a_Minimum  = (mu_Minimum * leptonPz) / (leptonE * leptonE - leptonPz * leptonPz);
        pznu = a_Minimum;

        double Enu = sqrt(minPx * minPx + minPy * minPy + pznu * pznu);
        result.setXYZE(minPx, minPy, pznu , Enu);


      }
      return std::make_pair(result,result2);
    }
  
  
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2019_I1744604);


    void init() {

      Cut particle_cut = (Cuts::abseta < 5.0) and (Cuts::pT > 0.*MeV);
      Cut lepton_cut   = (Cuts::abseta < 2.4) and (Cuts::pT > 26.*GeV);
      
      // Generic final state
      FinalState fs(particle_cut);

      // Dressed leptons
      ChargedLeptons charged_leptons(fs);
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      PromptFinalState prompt_leptons(charged_leptons);
      prompt_leptons.acceptMuonDecays(true);
      prompt_leptons.acceptTauDecays(true);
      PromptFinalState prompt_photons(photons);
      prompt_photons.acceptMuonDecays(true);
      prompt_photons.acceptTauDecays(true);

      DressedLeptons dressed_leptons(
        prompt_photons, prompt_leptons, 0.1,
        lepton_cut, true, true
      );
      declare(dressed_leptons, "DressedLeptons");


      // Jets
      VetoedFinalState fsForJets(fs);
      fsForJets.addVetoOnThisFinalState(dressed_leptons);
      declare(
        //excludes all neutrinos by default
        FastJets(fsForJets, FastJets::ANTIKT, 0.4), 
        "Jets"
      );

      // Neutrinos
      IdentifiedFinalState neutrinos(fs);
      neutrinos.acceptNeutrinos();
      PromptFinalState prompt_neutrinos(neutrinos);
      prompt_neutrinos.acceptMuonDecays(true);
      prompt_neutrinos.acceptTauDecays(true);
      declare(prompt_neutrinos, "Neutrinos");
      
      //Partonic top (for differentiating between t and tbar events only)
      declare(PartonicTops(),"TopQuarks");
      

      book(_hist_abs_top_pt,"d13-x01-y01",_binning_top_pt);
      book(_hist_norm_top_pt,"d37-x01-y01",_binning_top_pt);
      book(_hist_ratio_top_pt,"d59-x01-y01",_binning_top_pt);


      book(_hist_t_top_pt,"t_top_pt",_binning_top_pt);
      book(_hist_tbar_top_pt,"tbar_top_pt",_binning_top_pt);

      

      //book(_hist_abs_top_pt,"d13-x01-y01");
      //top y 0.,0.2,0.5,0.8,1.3,2.6
      //l pt 26.,35.,45.,60.,85.,200.
      //l y 0.0,0.4,0.8,1.5,1.9,2.4
      //w pt 0.,35.,55.,80.,140.,250.
      //cos -1.0,-0.6,-0.3,0.0,0.3,0.6,1.0

      /*
      // Book histograms
      // specify custom binning
      book(_h["XXXX"], "myh1", 20, 0.0, 100.0);
      book(_h["YYYY"], "myh2", logspace(20, 1e-2, 1e3));
      book(_h["ZZZZ"], "myh3", {0.0, 1.0, 2.0, 4.0, 8.0, 16.0});
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_h["AAAA"], 1, 1, 1);
      book(_p["BBBB"], 2, 1, 1);
      book(_c["CCCC"], 3, 1, 1);
      */

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      vector<DressedLepton> leptons = applyProjection<DressedLeptons>(
        event,
        "DressedLeptons"
      ).dressedLeptons();
      
      if (leptons.size() != 1) {
        return;
      }
      
      Cut jet_cut((Cuts::abseta < 4.7) and (Cuts::pT > 40.*GeV));
      vector<Jet> jets = applyProjection<FastJets>(
        event,
        "Jets"
      ).jets(jet_cut);
      
      if (jets.size() != 2) {
        return;
      }
      
      vector<Particle> neutrinos = applyProjection<PromptFinalState>(
        event,
        "Neutrinos"
      ).particles();
      
      if (neutrinos.size() == 0) {
        return;
      }

      vector<Particle> topQuarks = applyProjection<PartonicTops>(
        event,
        "TopQuarks"
      ).tops();
      
      if (topQuarks.size() != 1) {
        return;
      }

      std::cout<<"Top quarks: "<<topQuarks.size()<<std::endl;
      for (size_t i = 0; i < topQuarks.size(); ++i)
      {
        std::cout<<" - "<<topQuarks[i].pid()<<", "<<topQuarks[i].pt()<<std::endl;
      }
      

      _hist_t_top_pt->fill(leptons[0].pt()/GeV);
      _hist_tbar_top_pt->fill(leptons[0].pt()/GeV);

      /*
      Cut jet_cut = (Cuts::abseta < _jetMaxEta) and (Cuts::pT > _jetMinPt*GeV);
      _jets = applyProjection<FastJets>(event, "Jets").jetsByPt(jet_cut);
      for (const Jet& jet : _jets) {
        if (jet.bTagged()) _bjets.push_back(jet);
        else               _ljets.push_back(jet);
      }

      _nujets = applyProjection<FastJets>(event, "NuJets").jetsByPt(jet_cut);

      Cut fatjet_cut = (Cuts::abseta < _fatJetMaxEta) and (Cuts::pT > _fatJetMinPt*GeV);
      _fatjets = applyProjection<FastJets>(event, "FatJets").jetsByPt(fatjet_cut);

      _photons = applyProjection<FinalState>(event, "Photons").particlesByPt();

      _neutrinos = applyProjection<FinalState>(event, "Neutrinos").particlesByPt();

      _met = -applyProjection<MissingMomentum>(event, "MET").vectorEt();




      // the final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      // FinalState of prompt photons and bare muons and electrons in the event
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);

      // dress the prompt bare leptons with prompt photons within dR < 0.1
      // apply some fiducial cuts on the dressed leptons
      Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);

      /// @todo Do the event by event analysis here

      // retrieve dressed leptons, sorted by pT
      vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons();

      // retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);

      // remove all jets within dR < 0.2 of a dressed lepton
      idiscardIfAnyDeltaRLess(jets, leptons, 0.2);

      // select jets ghost-associated to B-hadrons with a certain fiducial selection
      Jets bjets = filter_select(jets, [](const Jet& jet) {
        return  jet.bTagged(Cuts::pT > 5*GeV && Cuts::abseta < 2.5);
      });

      // veto event if there are no b-jets
      if (bjets.empty())  vetoEvent;

      // apply a missing-momentum cut
      if (apply<MissingMomentum>(event, "MET").missingPt() < 30*GeV)  vetoEvent;
      */
      // fill histogram with leading b-jet pT
      //_h["XXXX"]->fill(bjets[0].pT()/GeV);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //normalize(_h["YYYY"]); // normalize to unity
      scale(_hist_t_top_pt, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_tbar_top_pt, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      
      for (size_t ipt = 0; ipt < _hist_t_top_pt->numBins(); ++ ipt) {
        const auto& t_bin = _hist_t_top_pt->bin(ipt);
        const auto& tbar_bin = _hist_tbar_top_pt->bin(ipt);
        auto& abs_bin = _hist_abs_top_pt->bin(ipt);
        auto& norm_bin = _hist_norm_top_pt->bin(ipt);
        auto& ratio_bin = _hist_ratio_top_pt->bin(ipt);
        
        abs_bin.fillBin(t_bin.height()+tbar_bin.height());
        norm_bin.fillBin(t_bin.height()+tbar_bin.height());
        ratio_bin.fillBin(t_bin.height()/(t_bin.height()+tbar_bin.height()));
      }
    }

    static const std::vector<double> _binning_top_pt; 
     
    Histo1DPtr _hist_t_top_pt;
    Histo1DPtr _hist_tbar_top_pt;
    Histo1DPtr _hist_abs_top_pt;
    Histo1DPtr _hist_norm_top_pt;
    Histo1DPtr _hist_ratio_top_pt;
  };

  const std::vector<double> CMS_2019_I1744604::_binning_top_pt{{
    0.,50.,80.,120.,180.,300.
  }};

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2019_I1744604);


}
