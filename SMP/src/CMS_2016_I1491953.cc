// -*- C++ -*-

//********** Bhawandeep & Apichart **********
//******** WJets@8TeV ************

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/WFinder.hh"

#include "Rivet/AnalysisLoader.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Tools/RivetYODA.hh"

#include <iostream>

namespace Rivet {

  //--- define bool used in sorting
  bool orderByIncRap(const FourMomentum& a, const FourMomentum& b) {
		return (a.rapidity() < b.rapidity());
  }

  /// @brief Add a short analysis description here
  class CMS_2016_I1491953 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2016_I1491953);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

        // Initialise and register projections
		FinalState fs;
		WFinder wfinder_mu(fs, Cuts::abseta < 2.4 && Cuts::pT > 0*GeV, PID::MUON, 0*GeV, 1000000*GeV, 0*GeV, 0.1, WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);
		addProjection(wfinder_mu, "WFinder_mu");
		
		// Define veto FS
		VetoedFinalState vfs;
		vfs.addVetoOnThisFinalState(wfinder_mu);
		vfs.addVetoPairId(PID::MUON);
		vfs.vetoNeutrinos();
		
		FastJets fastjets(vfs, FastJets::ANTIKT, 0.5);
		addProjection(fastjets, "Jets");

		//----- define variable binning
		double bins_addjetPt_Winc1jet[] = {30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 500, 590, 700, 1000};
		double bins_addjetPt_Winc2jet[] = {30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 500, 590, 800};
		double bins_addjetPt_Winc3jet[] = {30, 39, 49, 62, 78, 105, 142, 185, 235, 300};
		double bins_addjetPt_Winc4jet[] = {30, 39, 49, 62, 78, 96, 150};
		
		double bins_addjetHT_Winc1jet[] = {30, 39, 49, 62, 78, 96, 118, 150, 190, 240, 300, 370, 450, 540, 650, 800, 1000, 1500};
		double bins_addjetHT_Winc2jet[] = {60, 78, 96, 118, 150, 190, 240, 300, 370, 450, 540, 650, 800, 1200};
		double bins_addjetHT_Winc3jet[] = {90, 105, 125, 151, 185, 230, 290, 366, 466, 586, 767, 990};
		double bins_addjetHT_Winc4jet[] = {120, 140, 167, 203, 253, 320, 410, 530, 690, 910};
		
		double bins_dijetPt_Winc2jet[] = {20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 500, 590, 800};
		double bins_dijetPt_Winc3jet[] = {20, 24, 30, 39, 49, 62, 78, 105, 142, 185, 235, 300};
		double bins_dijetPt_Winc4jet[] = {20, 24, 30, 39, 49, 62, 78, 96, 150};

		double bins_dijetM_Winc2jet[] = {0, 25, 52, 81, 112, 145, 180, 217, 256, 297, 340, 385, 432, 481, 532, 585, 640, 700};
		double bins_dijetM_Winc3jet[] = {0, 30, 62, 96, 132, 170, 210, 252, 296, 342, 390, 440, 492, 546, 602, 660};
		double bins_dijetM_Winc4jet[] = {0, 34, 70, 108, 148, 190, 234, 280, 328, 378, 430, 484, 540, 598, 660};
		
		double bins_njets_meanNj[] = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5};

		std::vector<double> addjetPt_Winc1jet(bins_addjetPt_Winc1jet,bins_addjetPt_Winc1jet+sizeof(bins_addjetPt_Winc1jet)/sizeof(double));
		std::vector<double> addjetPt_Winc2jet(bins_addjetPt_Winc2jet,bins_addjetPt_Winc2jet+sizeof(bins_addjetPt_Winc2jet)/sizeof(double));
		std::vector<double> addjetPt_Winc3jet(bins_addjetPt_Winc3jet,bins_addjetPt_Winc3jet+sizeof(bins_addjetPt_Winc3jet)/sizeof(double));
		std::vector<double> addjetPt_Winc4jet(bins_addjetPt_Winc4jet,bins_addjetPt_Winc4jet+sizeof(bins_addjetPt_Winc4jet)/sizeof(double));

		std::vector<double> addjetHT_Winc1jet(bins_addjetHT_Winc1jet,bins_addjetHT_Winc1jet+sizeof(bins_addjetHT_Winc1jet)/sizeof(double));
		std::vector<double> addjetHT_Winc2jet(bins_addjetHT_Winc2jet,bins_addjetHT_Winc2jet+sizeof(bins_addjetHT_Winc2jet)/sizeof(double));
		std::vector<double> addjetHT_Winc3jet(bins_addjetHT_Winc3jet,bins_addjetHT_Winc3jet+sizeof(bins_addjetHT_Winc3jet)/sizeof(double));
		std::vector<double> addjetHT_Winc4jet(bins_addjetHT_Winc4jet,bins_addjetHT_Winc4jet+sizeof(bins_addjetHT_Winc4jet)/sizeof(double));
		
		std::vector<double> dijetPt_Winc2jet(bins_dijetPt_Winc2jet,bins_dijetPt_Winc2jet+sizeof(bins_dijetPt_Winc2jet)/sizeof(double));
		std::vector<double> dijetPt_Winc3jet(bins_dijetPt_Winc3jet,bins_dijetPt_Winc3jet+sizeof(bins_dijetPt_Winc3jet)/sizeof(double));
		std::vector<double> dijetPt_Winc4jet(bins_dijetPt_Winc4jet,bins_dijetPt_Winc4jet+sizeof(bins_dijetPt_Winc4jet)/sizeof(double));
		
		std::vector<double> dijetM_Winc2jet(bins_dijetM_Winc2jet,bins_dijetM_Winc2jet+sizeof(bins_dijetM_Winc2jet)/sizeof(double));
		std::vector<double> dijetM_Winc3jet(bins_dijetM_Winc3jet,bins_dijetM_Winc3jet+sizeof(bins_dijetM_Winc3jet)/sizeof(double));
		std::vector<double> dijetM_Winc4jet(bins_dijetM_Winc4jet,bins_dijetM_Winc4jet+sizeof(bins_dijetM_Winc4jet)/sizeof(double));
		
		std::vector<double> njets_meanNj(bins_njets_meanNj,bins_njets_meanNj+sizeof(bins_njets_meanNj)/sizeof(double));
		//-------------
		
      // Book histograms
      //_h_XXXX = bookProfile1D(1, 1, 1);
      //_h_YYYY = bookHisto1D(2, 1, 1);
      //_h_ZZZZ = bookCounter(3, 1, 1);
		
		//-------------
		_hist_Mult_exc      = bookHisto1D("d01-x01-y01", 8, -0.5, 7.5);
		_hist_inc_WJetMult  = bookHisto1D("d02-x01-y01", 8, -0.5, 7.5);
		
		//-------------
		_hist_addJetPt1j = bookHisto1D("d03-x01-y01", addjetPt_Winc1jet);
		_hist_addJetPt2j = bookHisto1D("d04-x01-y01", addjetPt_Winc2jet);
		_hist_addJetPt3j = bookHisto1D("d05-x01-y01", addjetPt_Winc3jet);
		_hist_addJetPt4j = bookHisto1D("d06-x01-y01", addjetPt_Winc4jet);
		
		//-------------
		_hist_addHt_1j = bookHisto1D("d07-x01-y01", addjetHT_Winc1jet);
		_hist_addHt_2j = bookHisto1D("d08-x01-y01", addjetHT_Winc2jet);
		_hist_addHt_3j = bookHisto1D("d09-x01-y01", addjetHT_Winc3jet);
		_hist_addHt_4j = bookHisto1D("d10-x01-y01", addjetHT_Winc4jet);
		
		//-------------
		_hist_diJetPt_2j = bookHisto1D("d11-x01-y01", dijetPt_Winc2jet);
		_hist_diJetPt_3j = bookHisto1D("d12-x01-y01", dijetPt_Winc3jet);
		_hist_diJetPt_4j = bookHisto1D("d13-x01-y01", dijetPt_Winc4jet);

		//-------------
		_hist_dijetM_2j = bookHisto1D("d14-x01-y01", dijetM_Winc2jet);
		_hist_dijetM_3j = bookHisto1D("d15-x01-y01", dijetM_Winc3jet);
		_hist_dijetM_4j = bookHisto1D("d16-x01-y01", dijetM_Winc4jet);
		
		//-------------
		_hist_Jeteta1j = bookHisto1D("d17-x01-y01", 32, 0, 2.4);
		_hist_Jeteta2j = bookHisto1D("d18-x01-y01", 32, 0, 2.4);
		_hist_Jeteta3j = bookHisto1D("d19-x01-y01", 12, 0, 2.4);
		_hist_Jeteta4j = bookHisto1D("d20-x01-y01", 12, 0, 2.4);
		
		//-------------
		_hist_dyj1j2_2j = bookHisto1D("d21-x01-y01", 20, 0, 4.8);
		_hist_dyj1j2_3j = bookHisto1D("d22-x01-y01", 20, 0, 4.8);
		_hist_dyj1j2_4j = bookHisto1D("d23-x01-y01", 16, 0, 4.8);
		
		_hist_dyj1j3_3j = bookHisto1D("d24-x01-y01", 20, 0, 4.8);
		_hist_dyj2j3_3j = bookHisto1D("d25-x01-y01", 20, 0, 4.8);

		_hist_dyjFjB_2j = bookHisto1D("d26-x01-y01", 20, 0, 4.8);
		_hist_dyjFjB_3j = bookHisto1D("d27-x01-y01", 20, 0, 4.8);
		_hist_dyjFjB_4j = bookHisto1D("d28-x01-y01", 16, 0, 4.8);
		
		_hist_dphij1j2_2j = bookHisto1D("d29-x01-y01", 20, 0, 3.14159265359);
		_hist_dphijFjB_2j = bookHisto1D("d30-x01-y01", 20, 0, 3.14159265359);
		_hist_dRj1j2_2j = bookHisto1D("d31-x01-y01", 30, 0, 6.);
		
		//-------------
		_hist_dphij1mu_1j = bookHisto1D("d32-x01-y01", 20, 0, 3.14159265359);
		_hist_dphij2mu_2j = bookHisto1D("d33-x01-y01", 20, 0, 3.14159265359);
		_hist_dphij3mu_3j = bookHisto1D("d34-x01-y01", 16, 0, 3.14159265359);
		_hist_dphij4mu_4j = bookHisto1D("d35-x01-y01", 16, 0, 3.14159265359);

		//------- MeanNJ ------
		_hist_MeanNJht_1j     = bookProfile1D("d36-x01-y01", addjetHT_Winc1jet);
		_hist_MeanNJht_2j     = bookProfile1D("d37-x01-y01", addjetHT_Winc2jet);
		_hist_MeanNJdyj1j2_2j = bookProfile1D("d38-x01-y01", 20, 0, 4.8);
		_hist_MeanNJdyjFjB_2j = bookProfile1D("d39-x01-y01", 20, 0, 4.8);

    }

	  
	// define function used for filiing inc Njets histo
	void Fill(Histo1DPtr& _histJetMult, const double& weight, std::vector<FourMomentum>& finaljet_list){

		  _histJetMult->fill(0, weight);
		  for (size_t i=0 ; i<finaljet_list.size() ; ++i) {
			  if (i==7) break;
			  _histJetMult->fill(i+1, weight);  // inclusive multiplicity
		  }
	}

	  
    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// @todo Do the event by event analysis here
		
		const double weight = event.weight();
		const WFinder& wfinder_mu = applyProjection<WFinder>(event, "WFinder_mu");
		
		if (wfinder_mu.bosons().size() != 1) {
			vetoEvent;
		}
		
		if (wfinder_mu.bosons().size() == 1) {
			
			const FourMomentum& lepton0 = wfinder_mu.constituentLeptons()[0].momentum();
			const FourMomentum& neutrino = wfinder_mu.constituentNeutrinos()[0].momentum();
			double WmT = sqrt( 2 * lepton0.pT() * neutrino.pT() * (1 - cos(deltaPhi(lepton0, neutrino))) );
			
			if (WmT < 50.0*GeV) vetoEvent;
			
			double pt0 = lepton0.pT();
			double eta0 = lepton0.eta();
			
			if ( (fabs(eta0) > 2.1) || (pt0 < 25.0*GeV) ) vetoEvent;
			
			//--- Obtain the jets::::::::::::::
			vector<FourMomentum> finaljet_list;
			double HT = 0.0;
			
			//--- loop over jets in a event, pushback in finaljet_list collection
			foreach (const Jet& j, applyProjection<FastJets>(event, "Jets").jetsByPt(30.0*GeV)) {
				
				const double jeta = j.momentum().eta();
				const double jpt = j.momentum().pT();
				
				if ( (fabs(jeta) < 2.4) && (deltaR(lepton0, j.momentum()) > 0.5) ) {
					
					if(jpt > 30.0*GeV) {
						finaljet_list.push_back(j.momentum());
						HT += j.momentum().pT();
					}
				}
			} // end looping over jets
			
			//--- new jet_list sorted by increasing rapidity
			vector<FourMomentum> jListRap = finaljet_list;
			std::sort(jListRap.begin(), jListRap.end(), orderByIncRap);
			
			//--- Multiplicity exc plot.
			if(finaljet_list.size()<=7) {
				_hist_Mult_exc->fill(finaljet_list.size(), weight);
			}
			else if (finaljet_list.size()>7){
				_hist_Mult_exc->fill(7., weight);
			}
			
			//--- Multiplicity inc plot.
			Fill(_hist_inc_WJetMult, weight, finaljet_list);
			
			if(finaljet_list.size()>=1) {
				_hist_addJetPt1j->fill(finaljet_list[0].pT(), weight);
				_hist_Jeteta1j->fill(fabs(finaljet_list[0].eta()), weight);
				_hist_addHt_1j->fill(HT, weight);
				_hist_dphij1mu_1j->fill( deltaPhi(finaljet_list[0].phi(), lepton0.phi()), weight );
				_hist_MeanNJht_1j->fill( HT, finaljet_list.size(), weight);
			}
			
			if(finaljet_list.size()>=2) {
				_hist_addJetPt2j->fill(finaljet_list[1].pT(), weight);
				_hist_Jeteta2j->fill(fabs(finaljet_list[1].eta()), weight);
				_hist_addHt_2j->fill(HT, weight);
				
				_hist_dyj1j2_2j   ->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
				_hist_dyjFjB_2j   ->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
				_hist_dphij1j2_2j ->fill( deltaPhi(finaljet_list[0].phi(), finaljet_list[1].phi()), weight);
				_hist_dphijFjB_2j ->fill( deltaPhi(jListRap[0].phi(), jListRap[jListRap.size()-1].phi()) , weight);
				
				_hist_dijetM_2j   ->fill( (add(finaljet_list[0], finaljet_list[1])).mass(), weight);
				_hist_diJetPt_2j  ->fill( (add(finaljet_list[0], finaljet_list[1])).pT(), weight);
				_hist_dRj1j2_2j   ->fill( deltaR(finaljet_list[0].rapidity(), finaljet_list[0].phi(), finaljet_list[1].rapidity(), finaljet_list[1].phi()), weight);
				
				_hist_dphij2mu_2j ->fill( deltaPhi(finaljet_list[1].phi(), lepton0.phi()), weight );
				
				_hist_MeanNJht_2j->fill( HT, finaljet_list.size(), weight);
				_hist_MeanNJdyj1j2_2j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), finaljet_list.size(), weight);
				_hist_MeanNJdyjFjB_2j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), finaljet_list.size(), weight);
			}
			
			if(finaljet_list.size()>=3) {
				_hist_addJetPt3j->fill(finaljet_list[2].pT(), weight);
				_hist_Jeteta3j->fill(fabs(finaljet_list[2].eta()), weight);
				_hist_addHt_3j->fill(HT, weight);
				
				_hist_dyj1j2_3j     ->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
				_hist_dyj1j3_3j     ->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[2].rapidity()), weight);
				_hist_dyj2j3_3j     ->fill( fabs(finaljet_list[1].rapidity() - finaljet_list[2].rapidity()), weight);
				_hist_dyjFjB_3j     ->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
				
				_hist_dijetM_3j  ->fill( (add(finaljet_list[0], finaljet_list[1])).mass(), weight);
				_hist_diJetPt_3j ->fill( (add(finaljet_list[0], finaljet_list[1])).pT(), weight);
				
				_hist_dphij3mu_3j->fill( deltaPhi(finaljet_list[2].phi(), lepton0.phi()), weight );
			}
			
			if(finaljet_list.size()>=4) {
				_hist_addJetPt4j->fill(finaljet_list[3].pT(), weight);
				_hist_Jeteta4j->fill(fabs(finaljet_list[3].eta()), weight);
				_hist_addHt_4j->fill(HT, weight);
				
				_hist_dyj1j2_4j     ->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
				_hist_dyjFjB_4j     ->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
				
				_hist_dijetM_4j  ->fill( (add(finaljet_list[0], finaljet_list[1])).mass(), weight);
				_hist_diJetPt_4j ->fill( (add(finaljet_list[0], finaljet_list[1])).pT(), weight);
				_hist_dphij4mu_4j->fill( deltaPhi(finaljet_list[3].phi(), lepton0.phi()), weight );
			}
		} //////

    } //void loop


    /// Normalise histograms etc., after the run
    void finalize() {

      //normalize(_h_YYYY); // normalize to unity
      //scale(_h_ZZZZ, crossSection()/picobarn/sumOfWeights()); // norm to cross section
		
		double crossec = 36703;    //inclusive nnlo xsec calculated by FEWZ
		//double crossec = 6440.58;  // W1jets exclusive
		//double crossec = 2087.225; // W2jets exclusive
		//double crossec = 619.0113; // W3jets exclusive
		//double crossec = 255.2378; // W4jets exclusive
		
		scale(_hist_Mult_exc, crossec/sumOfWeights());
		scale(_hist_inc_WJetMult, crossec/sumOfWeights());
		
		scale(_hist_addJetPt1j, crossec/sumOfWeights());
		scale(_hist_addJetPt2j, crossec/sumOfWeights());
		scale(_hist_addJetPt3j, crossec/sumOfWeights());
		scale(_hist_addJetPt4j, crossec/sumOfWeights());
		
		scale(_hist_Jeteta1j, crossec/sumOfWeights());
		scale(_hist_Jeteta2j, crossec/sumOfWeights());
		scale(_hist_Jeteta3j, crossec/sumOfWeights());
		scale(_hist_Jeteta4j, crossec/sumOfWeights());
		
		scale(_hist_addHt_1j, crossec/sumOfWeights());
		scale(_hist_addHt_2j, crossec/sumOfWeights());
		scale(_hist_addHt_3j, crossec/sumOfWeights());
		scale(_hist_addHt_4j, crossec/sumOfWeights());
		
		//-------------------------------------
		scale(_hist_dyj1j2_2j, crossec/sumOfWeights());
		scale(_hist_dyj1j2_3j, crossec/sumOfWeights());
		scale(_hist_dyj1j2_4j, crossec/sumOfWeights());
		
		scale(_hist_dyjFjB_2j, crossec/sumOfWeights());
		scale(_hist_dyjFjB_3j, crossec/sumOfWeights());
		scale(_hist_dyjFjB_4j, crossec/sumOfWeights());
		
		scale(_hist_dyj1j3_3j, crossec/sumOfWeights());
		scale(_hist_dyj2j3_3j, crossec/sumOfWeights());
		
		scale(_hist_dphij1j2_2j, crossec/sumOfWeights());
		scale(_hist_dphijFjB_2j, crossec/sumOfWeights());
		
		scale(_hist_dRj1j2_2j, crossec/sumOfWeights());
		
		scale(_hist_dijetM_2j, crossec/sumOfWeights());
		scale(_hist_dijetM_3j, crossec/sumOfWeights());
		scale(_hist_dijetM_4j, crossec/sumOfWeights());
		
		scale(_hist_diJetPt_2j, crossec/sumOfWeights());
		scale(_hist_diJetPt_3j, crossec/sumOfWeights());
		scale(_hist_diJetPt_4j, crossec/sumOfWeights());
		
		scale(_hist_dphij1mu_1j, crossec/sumOfWeights());
		scale(_hist_dphij2mu_2j, crossec/sumOfWeights());
		scale(_hist_dphij3mu_3j, crossec/sumOfWeights());
		scale(_hist_dphij4mu_4j, crossec/sumOfWeights());
		
		//-------------------------------------
		//scale(_hist_MeanNJht_1j, crossec/sumOfWeights());
		//scale(_hist_MeanNJht_2j, crossec/sumOfWeights());
		//scale(_hist_MeanNJdyj1j2_2j, crossec/sumOfWeights());
		//scale(_hist_MeanNJdyjFjB_2j, crossec/sumOfWeights());

    }

    //@}


  private:


    /// @name Histograms
    //@{
    //Profile1DPtr _h_XXXX;
    //Histo1DPtr _h_YYYY;
    //CounterPtr _h_ZZZZ;
	  
	  Histo1DPtr _hist_inc_WJetMult;
	  Histo1DPtr _hist_Mult_exc;
	  
	  Histo1DPtr _hist_addJetPt1j;
	  Histo1DPtr _hist_addJetPt2j;
	  Histo1DPtr _hist_addJetPt3j;
	  Histo1DPtr _hist_addJetPt4j;
	  
	  Histo1DPtr _hist_Jeteta1j;
	  Histo1DPtr _hist_Jeteta2j;
	  Histo1DPtr _hist_Jeteta3j;
	  Histo1DPtr _hist_Jeteta4j;
	  
	  Histo1DPtr _hist_addHt_1j;
	  Histo1DPtr _hist_addHt_2j;
	  Histo1DPtr _hist_addHt_3j;
	  Histo1DPtr _hist_addHt_4j;
	  
	  //-------------------------------------
	  Histo1DPtr _hist_dyj1j2_2j;
	  Histo1DPtr _hist_dyj1j2_3j;
	  Histo1DPtr _hist_dyj1j2_4j;
	  
	  Histo1DPtr _hist_dyjFjB_2j;
	  Histo1DPtr _hist_dyjFjB_3j;
	  Histo1DPtr _hist_dyjFjB_4j;
	  
	  Histo1DPtr _hist_dyj1j3_3j;
	  Histo1DPtr _hist_dyj2j3_3j;
	  
	  Histo1DPtr _hist_dphij1j2_2j;
	  Histo1DPtr _hist_dphijFjB_2j;
	  
	  Histo1DPtr _hist_dRj1j2_2j;
	  
	  Histo1DPtr _hist_dijetM_2j;
	  Histo1DPtr _hist_dijetM_3j;
	  Histo1DPtr _hist_dijetM_4j;
	  
	  Histo1DPtr _hist_diJetPt_2j;
	  Histo1DPtr _hist_diJetPt_3j;
	  Histo1DPtr _hist_diJetPt_4j;
	  
	  Histo1DPtr _hist_dphij1mu_1j;
	  Histo1DPtr _hist_dphij2mu_2j;
	  Histo1DPtr _hist_dphij3mu_3j;
	  Histo1DPtr _hist_dphij4mu_4j;
	  
	  //-------------------------------------	  
	  Profile1DPtr _hist_MeanNJht_1j;
	  Profile1DPtr _hist_MeanNJht_2j;
	  Profile1DPtr _hist_MeanNJdyj1j2_2j;
	  Profile1DPtr _hist_MeanNJdyjFjB_2j;

	  
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1491953);


}
