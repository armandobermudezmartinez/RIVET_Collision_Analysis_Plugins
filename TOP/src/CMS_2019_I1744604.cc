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

    static std::vector<double> EquationSolve(double a, double b, double c, double d) {
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
    ) {
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
    
    static void fillAbsHist(Histo1DPtr& hist_abs, const Histo1DPtr& hist_t, const Histo1DPtr& hist_tbar) {
      (*hist_abs)+=(*hist_t);
      (*hist_abs)+=(*hist_tbar);
      /*
      for (size_t i = 0; i < hist_abs->numBins(); ++i) {
        const auto& t_bin = hist_t->bin(i);
        const auto& tbar_bin = hist_tbar->bin(i);
        hist_abs->bin(i).fillBin(t_bin.height()+tbar_bin.height());
      }
      */
    }
    
    static void fillNormHist(Histo1DPtr& hist_norm, const Histo1DPtr& hist_t, const Histo1DPtr& hist_tbar) {
      for (size_t i = 0; i < hist_norm->numBins(); ++i) {
        const auto& t_bin = hist_t->bin(i);
        const auto& tbar_bin = hist_tbar->bin(i);
        hist_norm->bin(i).fillBin(t_bin.height()+tbar_bin.height());
      }
    }
    
    static void fillRatioHist(Histo1DPtr& hist_ratio, const Histo1DPtr& hist_t, const Histo1DPtr& hist_tbar) {
      for (size_t i = 0; i < hist_ratio->numBins(); ++i) {
        const auto& t_bin = hist_t->bin(i);
        const auto& tbar_bin = hist_tbar->bin(i);
        hist_ratio->bin(i).fillBin(t_bin.height()/(t_bin.height()+tbar_bin.height()));
      }
    }
    
    
    static double calcXsec(const Histo1DPtr& hist) {
      double xsec = 0.;
      for (size_t i = 0; i < hist->numBins(); ++i) {
        auto bin = hist->bin(i);
        xsec+=bin.height()*bin.width();
      }
      return xsec;
    }
  
  
  public:
    CMS_2019_I1744604(): 
      Analysis("CMS_2019_I1744604"),
      _sum_t_weights(0.),
      _sum_tbar_weights(0.) {
    }

    


    void init() {

      Cut particle_cut = (Cuts::abseta < 5.0) and (Cuts::pT > 0.*MeV);
      Cut lepton_cut   = (Cuts::abseta < 2.4) and (Cuts::pT > 26.*GeV);
      

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
        // excludes all neutrinos by default
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
      
      // Partonic top (for differentiating between t and tbar events only)
      declare(PartonicTops(),"TopQuarks");

      book(_hist_abs_top_pt,"d13-x01-y01",_binning_top_pt);
      book(_hist_norm_top_pt,"d37-x01-y01",_binning_top_pt);
      book(_hist_ratio_top_pt,"d59-x01-y01",_binning_top_pt);
      book(_hist_t_top_pt,"t_top_pt",_binning_top_pt);
      book(_hist_tbar_top_pt,"tbar_top_pt",_binning_top_pt);
      
      book(_hist_abs_top_y,"d15-x01-y01",_binning_top_y);
      book(_hist_norm_top_y,"d39-x01-y01",_binning_top_y);
      book(_hist_ratio_top_y,"d61-x01-y01",_binning_top_y);
      book(_hist_t_top_y,"t_top_y",_binning_top_y);
      book(_hist_tbar_top_y,"tbar_top_y",_binning_top_y);
            
      book(_hist_abs_lepton_pt,"d17-x01-y01",_binning_lepton_pt);
      book(_hist_norm_lepton_pt,"d41-x01-y01",_binning_lepton_pt);
      book(_hist_ratio_lepton_pt,"d63-x01-y01",_binning_lepton_pt);
      book(_hist_t_lepton_pt,"t_lepton_pt",_binning_lepton_pt);
      book(_hist_tbar_lepton_pt,"tbar_lepton_pt",_binning_lepton_pt);
      
      book(_hist_abs_lepton_y,"d19-x01-y01",_binning_lepton_y);
      book(_hist_norm_lepton_y,"d43-x01-y01",_binning_lepton_y);
      book(_hist_ratio_lepton_y,"d65-x01-y01",_binning_lepton_y);
      book(_hist_t_lepton_y,"t_lepton_y",_binning_lepton_y);
      book(_hist_tbar_lepton_y,"tbar_lepton_y",_binning_lepton_y);
      
      book(_hist_abs_w_pt,"d21-x01-y01",_binning_w_pt);
      book(_hist_norm_w_pt,"d45-x01-y01",_binning_w_pt);
      book(_hist_ratio_w_pt,"d67-x01-y01",_binning_w_pt);
      book(_hist_t_w_pt,"t_w_pt",_binning_w_pt);
      book(_hist_tbar_w_pt,"tbar_w_pt",_binning_w_pt);
      
      book(_hist_abs_top_cos,"d23-x01-y01",_binning_top_cos);
      book(_hist_norm_top_cos,"d47-x01-y01",_binning_top_cos);
      book(_hist_t_top_cos,"t_top_cos",_binning_top_cos);
      book(_hist_tbar_top_cos,"tbar_top_cos",_binning_top_cos);
      
    }


    void analyze(const Event& event) {
      vector<Particle> topQuarks = applyProjection<PartonicTops>(
        event,
        "TopQuarks"
      ).tops();
      
      if (topQuarks.size() != 1) {
        return;
      }
      
      if (topQuarks[0].charge() > 0) {
        _sum_t_weights += 1.;
      } else {
        _sum_tbar_weights += 1.;
      }
      
      
      vector<DressedLepton> leptons = applyProjection<DressedLeptons>(
        event,
        "DressedLeptons"
      ).dressedLeptons();
      
      
      if (leptons.size() != 1) {
        return;
      }
      std::cout<<leptons[0].pid()<<std::endl;
      
      Cut jet_cut((Cuts::abseta < 4.7) and (Cuts::pT > 40.*GeV));
      vector<Jet> jets = applyProjection<FastJets>(
        event,
        "Jets"
      ).jets(jet_cut);
      
      std::vector<Jet> cleanedJets;
      DeltaRLess dRFct(leptons[0],0.4);
      for (const auto& jet: jets) {
        if (not dRFct(jet)) {
          cleanedJets.push_back(jet);
        }
      }
      
      if (cleanedJets.size() != 2) {
        return;
      }
      
      vector<Particle> neutrinos = applyProjection<PromptFinalState>(
        event,
        "Neutrinos"
      ).particles();
      
      if (neutrinos.size() == 0) {
        return;
      }

      
      if (topQuarks[0].charge() > 0) {
        _hist_t_lepton_pt->fill(leptons[0].pt()/GeV);
        _hist_t_lepton_y->fill(leptons[0].absrapidity());
        
      } else {
        _hist_tbar_lepton_pt->fill(leptons[0].pt()/GeV);
        _hist_tbar_lepton_y->fill(leptons[0].absrapidity());
      }
    }

    void finalize() {
      std::cout<<"weights="<<sumOfWeights()<<", t="<<_sum_t_weights<<", tbar="<<_sum_tbar_weights<<", xsec="<<(crossSection()/picobarn)<<std::endl;
      std::cout<<"int (w/o overflow)="<<_hist_t_lepton_pt->integral(false)<<", int="<<_hist_t_lepton_pt->integral(true)<<", xsec="<<(0.5*_hist_t_lepton_pt->integral(false)*_t_xsec_fraction*crossSection()/picobarn/_sum_t_weights)<<std::endl;
      std::cout<<"int (w/o overflow)="<<_hist_tbar_lepton_pt->integral(false)<<", int="<<_hist_tbar_lepton_pt->integral(true)<<", xsec="<<(0.5*_hist_tbar_lepton_pt->integral(false)*_tbar_xsec_fraction*crossSection()/picobarn/_sum_tbar_weights)<<std::endl;
      

      scale(_hist_t_top_pt, 0.5*_t_xsec_fraction*crossSection()/picobarn/_sum_t_weights);
      scale(_hist_tbar_top_pt, 0.5*_tbar_xsec_fraction*crossSection()/picobarn/sumOfWeights());
      
      scale(_hist_t_top_y, 0.5*_t_xsec_fraction*crossSection()/picobarn/_sum_t_weights);
      scale(_hist_tbar_top_y, 0.5*_tbar_xsec_fraction*crossSection()/picobarn/_sum_tbar_weights);
      
      scale(_hist_t_lepton_pt, 0.5*_t_xsec_fraction*crossSection()/picobarn/_sum_t_weights);
      scale(_hist_tbar_lepton_pt, 0.5*_tbar_xsec_fraction*crossSection()/picobarn/_sum_tbar_weights);
      
      scale(_hist_t_lepton_y, 0.5*_t_xsec_fraction*crossSection()/picobarn/_sum_t_weights);
      scale(_hist_tbar_lepton_y, 0.5*_tbar_xsec_fraction*crossSection()/picobarn/_sum_tbar_weights);

      scale(_hist_t_w_pt,0.5*_t_xsec_fraction*crossSection()/picobarn/_sum_t_weights);
      scale(_hist_tbar_w_pt, 0.5*_tbar_xsec_fraction*crossSection()/picobarn/_sum_tbar_weights);
      
      scale(_hist_t_top_cos, 0.5*_t_xsec_fraction*crossSection()/picobarn/_sum_t_weights);
      scale(_hist_tbar_top_cos, 0.5*_tbar_xsec_fraction*crossSection()/picobarn/_sum_tbar_weights);
      
      
      fillAbsHist(_hist_abs_top_pt,_hist_t_top_pt,_hist_tbar_top_pt);
      fillNormHist(_hist_norm_top_pt,_hist_t_top_pt,_hist_tbar_top_pt);
      fillRatioHist(_hist_ratio_top_pt,_hist_t_top_pt,_hist_tbar_top_pt);
      
      fillAbsHist(_hist_abs_top_y,_hist_t_top_y,_hist_tbar_top_y);
      fillNormHist(_hist_norm_top_y,_hist_t_top_y,_hist_tbar_top_y);
      fillRatioHist(_hist_ratio_top_y,_hist_t_top_y,_hist_tbar_top_y);
      
      fillAbsHist(_hist_abs_lepton_pt,_hist_t_lepton_pt,_hist_tbar_lepton_pt);
      fillNormHist(_hist_norm_lepton_pt,_hist_t_lepton_pt,_hist_tbar_lepton_pt);
      fillRatioHist(_hist_ratio_lepton_pt,_hist_t_lepton_pt,_hist_tbar_lepton_pt);
      
      fillAbsHist(_hist_abs_lepton_y,_hist_t_lepton_y,_hist_tbar_lepton_y);
      fillNormHist(_hist_norm_lepton_y,_hist_t_lepton_y,_hist_tbar_lepton_y);
      fillRatioHist(_hist_ratio_lepton_y,_hist_t_lepton_y,_hist_tbar_lepton_y);
      
      fillAbsHist(_hist_abs_w_pt,_hist_t_w_pt,_hist_tbar_w_pt);
      fillNormHist(_hist_norm_w_pt,_hist_t_w_pt,_hist_tbar_w_pt);
      fillRatioHist(_hist_ratio_w_pt,_hist_t_w_pt,_hist_tbar_w_pt);
      
      fillAbsHist(_hist_abs_top_cos,_hist_t_top_cos,_hist_tbar_top_cos);
      fillNormHist(_hist_norm_top_cos,_hist_t_top_cos,_hist_tbar_top_cos);
      
      std::cout<<"int xsec: lpt="<<calcXsec(_hist_abs_lepton_pt)<<", ly="<<calcXsec(_hist_abs_lepton_y)<<std::endl;
    }

    //calculated with Hathor v2.1 for a top quark mass of 172.5 GeV at NLO in QCD 
    static constexpr double _t_xsec_fraction = 136.02/216.99;
    static constexpr double _tbar_xsec_fraction = 80.95/216.99;
    
    static const std::vector<double> _binning_top_pt;
    static const std::vector<double> _binning_top_y;
    static const std::vector<double> _binning_lepton_pt;
    static const std::vector<double> _binning_lepton_y;
    static const std::vector<double> _binning_w_pt;
    static const std::vector<double> _binning_top_cos;
    
    double _sum_t_weights;
    double _sum_tbar_weights;
    
    Histo1DPtr _hist_abs_top_pt;
    Histo1DPtr _hist_norm_top_pt;
    Histo1DPtr _hist_ratio_top_pt;
    Histo1DPtr _hist_t_top_pt;
    Histo1DPtr _hist_tbar_top_pt;
    
    Histo1DPtr _hist_abs_top_y;
    Histo1DPtr _hist_norm_top_y;
    Histo1DPtr _hist_ratio_top_y;
    Histo1DPtr _hist_t_top_y;
    Histo1DPtr _hist_tbar_top_y;
    
    Histo1DPtr _hist_abs_lepton_pt;
    Histo1DPtr _hist_norm_lepton_pt;
    Histo1DPtr _hist_ratio_lepton_pt;
    Histo1DPtr _hist_t_lepton_pt;
    Histo1DPtr _hist_tbar_lepton_pt;
    
    Histo1DPtr _hist_abs_lepton_y;
    Histo1DPtr _hist_norm_lepton_y;
    Histo1DPtr _hist_ratio_lepton_y;
    Histo1DPtr _hist_t_lepton_y;
    Histo1DPtr _hist_tbar_lepton_y;
    
    Histo1DPtr _hist_abs_w_pt;
    Histo1DPtr _hist_norm_w_pt;
    Histo1DPtr _hist_ratio_w_pt;
    Histo1DPtr _hist_t_w_pt;
    Histo1DPtr _hist_tbar_w_pt;
    
    Histo1DPtr _hist_abs_top_cos;
    Histo1DPtr _hist_norm_top_cos;
    Histo1DPtr _hist_t_top_cos;
    Histo1DPtr _hist_tbar_top_cos;
  };


  const std::vector<double> CMS_2019_I1744604::_binning_top_pt{{
    0.,50.,80.,120.,180.,300.
  }};
  const std::vector<double> CMS_2019_I1744604::_binning_top_y{{
    0.,0.2,0.5,0.8,1.3,2.6
  }};
  const std::vector<double> CMS_2019_I1744604::_binning_lepton_pt{{
    26.,35.,45.,60.,85.,200.
  }};
  const std::vector<double> CMS_2019_I1744604::_binning_lepton_y{{
    0.0,0.4,0.8,1.5,1.9,2.4
  }};
  const std::vector<double> CMS_2019_I1744604::_binning_w_pt{{
    0.,35.,55.,80.,140.,250.
  }};
  const std::vector<double> CMS_2019_I1744604::_binning_top_cos{{
    -1.0,-0.6,-0.3,0.0,0.3,0.6,1.0
  }};

  DECLARE_RIVET_PLUGIN(CMS_2019_I1744604);
}

