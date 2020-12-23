// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2020_PAS_SMP_20_007 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2020_PAS_SMP_20_007);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // the basic final-state projection: 
      // all final-state particles within 
      // the given eta acceptance: 4.7 (jet range) + 0.4 (cone size)
      const FinalState fs(Cuts::abseta < 5.1); 

      // the final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      // Book histograms
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_h["JetPt1"], 1, 1, 1);
      book(_h["JetPt2"], 2, 1, 1);
      book(_h["JetPt3"], 3, 1, 1);
      book(_h["JetPt4"], 4, 1, 1);
      book(_h["JetEta1"], 5, 1, 1);
      book(_h["JetEta2"], 6, 1, 1);
      book(_h["JetEta3"], 7, 1, 1);
      book(_h["JetEta4"], 8, 1, 1);
      
      book(_h["DeltaPhiSoft"], 9, 1, 1);
      book(_h["DeltaPhi3"], 10, 1, 1);
      book(_h["DeltaY"], 11, 1, 1);
      book(_h["DeltaPhiY"], 12, 1, 1);
      book(_h["DeltaPtSoft"], 13, 1, 1);
      book(_h["DeltaS"], 14, 1, 1);
      
      book(_h["DeltaPhiSoft_binNorm"], 15, 1, 1);
      book(_h["DeltaPhi3_binNorm"], 16, 1, 1);
      book(_h["DeltaY_binNorm"], 17, 1, 1);
      book(_h["DeltaPhiY_binNorm"], 18, 1, 1);
      book(_h["DeltaPtSoft_binNorm"], 19, 1, 1);
      book(_h["DeltaS_binNorm"], 20, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// @todo Do the event by event analysis here

      // retrieve clustered jets, sorted by pT, with a minimum pT cut 10 GeV and eta range 4.7 (similar to PFJet collection)
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::abseta < 4.7 && Cuts::pT > 10*GeV);
      
      // fill only if there are at least 4 jets
      if (jets.size() >= 4) {
        if (jets[0].momentum().pt() > 35.0 && jets[1].momentum().pt() > 30.0 && jets[2].momentum().pt() > 25.0 && jets[3].momentum().pt() > 20.0) {
		
		double pt0 = jets[0].momentum().pt();
		double pt1 = jets[1].momentum().pt();
		double pt2 = jets[2].momentum().pt();
		double pt3 = jets[3].momentum().pt();
		double phi0 = jets[0].momentum().phi();
		double phi1 = jets[1].momentum().phi();
		double phi2 = jets[2].momentum().phi();
		double phi3 = jets[3].momentum().phi();
				
		_h["JetPt1"]->fill(pt0);
		_h["JetPt2"]->fill(pt1);
		_h["JetPt3"]->fill(pt2);
		_h["JetPt4"]->fill(pt3);
	
		_h["JetEta1"]->fill(jets[0].momentum().eta());
		_h["JetEta2"]->fill(jets[1].momentum().eta());
		_h["JetEta3"]->fill(jets[2].momentum().eta());
		_h["JetEta4"]->fill(jets[3].momentum().eta());
	
		// delta phi and eta of the 2 soft jets
		_h["DeltaPhiSoft"]->fill(abs(deltaPhi(phi2,phi3)));
		_h["DeltaPhiSoft_binNorm"]->fill(abs(deltaPhi(phi2,phi3)));
		
		// delta pt between soft jets
        	double DptSoft = sqrt(pow(pt2*cos(phi2) + pt3*cos(phi3), 2) + pow(pt2*sin(phi2) + pt3*sin(phi3), 2))/(pt2 + pt3);
		_h["DeltaPtSoft"]->fill(DptSoft);
		_h["DeltaPtSoft_binNorm"]->fill(DptSoft);
		
    		// delta S
		if (pt0 > 50.0 && pt1 > 30.0 && pt2 > 30.0 && pt3 > 30.0) {
                	double phiH = atan2(pt0*sin(phi0) + pt1*sin(phi1) , pt0*cos(phi0) + pt1*cos(phi1));
                	double phiS = atan2(pt2*sin(phi2) + pt3*sin(phi3) , pt2*cos(phi2) + pt3*cos(phi3));
                	double DS = abs(deltaPhi(phiH,phiS));
			_h["DeltaS"]->fill(DS);
			_h["DeltaS_binNorm"]->fill(DS);
		}
		
		// delta Y: most remote jets in rapidity, find min & max eta
		double mineta = 99999;
		double maxeta = -99999;
		int minetapos = -1;
		int maxetapos = -1;
		for (int i=0;i<4;i++) {
			if (jets[i].momentum().eta() < mineta) {
				mineta = jets[i].momentum().eta();
				minetapos = i;
			}
			if (jets[i].momentum().eta() > maxeta) {
				maxeta = jets[i].momentum().eta();
				maxetapos = i;
			}
		}
		_h["DeltaY"]->fill(abs(jets[minetapos].momentum().eta()-jets[maxetapos].momentum().eta()));
		_h["DeltaY_binNorm"]->fill(abs(jets[minetapos].momentum().eta()-jets[maxetapos].momentum().eta()));
		
		// Delta phi Y: azimuthal angle between most remote jets in eta
		_h["DeltaPhiY"]->fill(abs(deltaPhi(jets[minetapos].momentum().phi(),jets[maxetapos].momentum().phi())));
		_h["DeltaPhiY_binNorm"]->fill(abs(deltaPhi(jets[minetapos].momentum().phi(),jets[maxetapos].momentum().phi())));
		
		// delta phi3
        	double minphi3 = 999;
        	for (int iphi1=0;iphi1<4;iphi1++) {
		    for (int iphi2=0;iphi2<4;iphi2++) {
		        for (int iphi3=0;iphi3<4;iphi3++) {
            		    if (!(iphi1 == iphi2 || iphi2 == iphi3 || iphi1 == iphi3)) {
            		        double tempphi3 = abs(deltaPhi(jets[iphi1].momentum().phi(),jets[iphi2].momentum().phi())) + abs(deltaPhi(jets[iphi2].momentum().phi(),jets[iphi3].momentum().phi()));
				if (tempphi3 < minphi3) minphi3 = tempphi3;
			    }
			}
		    }
		}
        	_h["DeltaPhi3"]->fill(minphi3);
		_h["DeltaPhi3_binNorm"]->fill(minphi3);
		
		
	}
      }
      

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h["JetPt1"], crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_h["JetPt2"], crossSection()/picobarn/sumOfWeights());
      scale(_h["JetPt3"], crossSection()/picobarn/sumOfWeights());
      scale(_h["JetPt4"], crossSection()/picobarn/sumOfWeights());
      scale(_h["JetEta1"], crossSection()/picobarn/sumOfWeights());
      scale(_h["JetEta2"], crossSection()/picobarn/sumOfWeights());
      scale(_h["JetEta3"], crossSection()/picobarn/sumOfWeights());
      scale(_h["JetEta4"], crossSection()/picobarn/sumOfWeights());
      scale(_h["DeltaPhiSoft"], crossSection()/picobarn/sumOfWeights());
      scale(_h["DeltaPhi3"], crossSection()/picobarn/sumOfWeights());
      scale(_h["DeltaY"], crossSection()/picobarn/sumOfWeights());
      scale(_h["DeltaPhiY"], crossSection()/picobarn/sumOfWeights());
      scale(_h["DeltaPtSoft"], crossSection()/picobarn/sumOfWeights());
      scale(_h["DeltaS"], crossSection()/picobarn/sumOfWeights());
      
      // create bin normalised histograms
      
      // Correct for binwidths: Rivet automatically normalises histograms to binwidth when plotting AFTER the normalisation executed here. 
      // So we must calculate an extra correction here so that finally our bin-normalised histograms end up around 1 as in the paper.
      // in YODA bin index starts from 0
      
      // For DeltaY, which has variable binwidths we need to do following steps
      // divide histograms by binwidth
      for (unsigned int i = 0; i < _h["DeltaY_binNorm"]->numBins(); i++) {
        _h["DeltaY_binNorm"]->bin(i).scaleW(1.0/_h["DeltaY_binNorm"]->bin(i).xWidth());
      }
      
      // normalise to average of first 4 bins
      scale(_h["DeltaY_binNorm"], 1.0/(_h["DeltaY_binNorm"]->integralRange(0,3)/4.0));
      
      // multiply again with binwidth
      for (unsigned int i = 0; i < _h["DeltaY_binNorm"]->numBins(); i++) {
        _h["DeltaY_binNorm"]->bin(i).scaleW(_h["DeltaY_binNorm"]->bin(i).xWidth());
      }
      
      // DeltaPhiSoft and DeltaPhi3 histograms have uniform binwidths, so multiply with first binwidth is sufficient
      scale(_h["DeltaPhiSoft_binNorm"], _h["DeltaPhiSoft_binNorm"]->bin(0).xWidth()/(_h["DeltaPhiSoft_binNorm"]->integralRange(0,4)/5.0));
      scale(_h["DeltaPhi3_binNorm"], _h["DeltaPhi3_binNorm"]->bin(0).xWidth()/(_h["DeltaPhi3_binNorm"]->integralRange(0,3)/4.0));
      
      // DeltaPhiY, DeltaPtSoft and DeltaS are normalised to last bin
      scale(_h["DeltaPhiY_binNorm"], _h["DeltaPhiY_binNorm"]->bin(11).xWidth()/_h["DeltaPhiY_binNorm"]->bin(11).sumW() );
      scale(_h["DeltaPtSoft_binNorm"], _h["DeltaPtSoft_binNorm"]->bin(7).xWidth()/_h["DeltaPtSoft_binNorm"]->bin(7).sumW() );
      scale(_h["DeltaS_binNorm"], _h["DeltaS_binNorm"]->bin(6).xWidth()/_h["DeltaS_binNorm"]->bin(6).sumW() );

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2020_PAS_SMP_20_007);


}
