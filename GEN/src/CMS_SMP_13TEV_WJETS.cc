#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...
#include "Rivet/Particle.fhh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/ParticleName.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/AnalysisLoader.hh"
//#include "Rivet/RivetAIDA.hh"

#include <iostream>

namespace Rivet {
    
    
    class CMS_SMP_13TEV_WJETS: public Analysis {
        public:
        
            // Constructors
            CMS_SMP_13TEV_WJETS()
                : Analysis("CMS_SMP_13TEV_WJETS")
            {
                //setBeams(PID::PROTON, PID::PROTON);
                setNeedsCrossSection(true);
            }
        
        public:
        
            // Book histograms and initialise projections before the run
            void init() {
                
                FinalState fs;
                WFinder wfinder_mu(fs, Cuts::abseta < 2.4 && Cuts::pT > 0*GeV, PID::MUON, 0*GeV, 1000000*GeV, 0*GeV, 0.1, WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);
               
                addProjection(wfinder_mu, "WFinder_mu");
                
                
                // Define veto FS
                VetoedFinalState vfs;
                vfs.addVetoOnThisFinalState(wfinder_mu);
                vfs.addVetoPairId(PID::MUON);
                vfs.vetoNeutrinos();
                
                FastJets fastjets(vfs, FastJets::ANTIKT, 0.4);
                addProjection(fastjets, "Jets");
                
                // ---- book all the required histograms with correct binning ----
                // define binning and book Histos
                vector<double> jetPt_Winc1jet ;
                vector<double> jetPt_Winc2jet ;
                vector<double> jetPt_Winc3jet ;
                
                vector<double> jetHT_Winc1jet ;
                vector<double> jetHT_Winc2jet ;
                vector<double> jetHT_Winc3jet ;

                jetPt_Winc1jet += 20, 24, 30, 39, 49, 62, 79, 105, 138, 181, 231, 294, 375, 494, 800;
                jetPt_Winc2jet += 20, 24, 30, 39, 49, 62, 78, 105, 142, 185, 235, 300, 380, 500;
                jetPt_Winc3jet += 20, 24, 30, 41, 59, 81, 110, 152, 200, 300;

                jetHT_Winc1jet += 30, 39, 49, 62, 78, 96, 118, 150, 190, 240, 300, 370, 450, 540, 650, 800, 1000, 1500;
                jetHT_Winc2jet += 60, 78, 96, 118, 150, 190, 240, 300, 370, 450, 540, 650, 800, 1200;
                jetHT_Winc3jet += 90, 118, 168, 220, 300, 400, 550, 780, 1100;

                _hist_JetPt1j =bookHisto1D("jet_Pt1jetcase", jetPt_Winc1jet);
                _hist_JetPt2j =bookHisto1D("jet_Pt2jetcase", jetPt_Winc2jet);
                _hist_JetPt3j =bookHisto1D("jet_Pt3jetcase", jetPt_Winc3jet);
               
                _hist_Ht_1j =bookHisto1D("JetsHT_inc1jet", jetHT_Winc1jet);
                _hist_Ht_2j =bookHisto1D("JetsHT_inc2jet", jetHT_Winc2jet);
                _hist_Ht_3j =bookHisto1D("JetsHT_inc3jet", jetHT_Winc3jet);
                
                //-------------
                _hist_inc_WJetMult = bookHisto1D("njetWJet_incl", 7, -0.5, 6.5);
                _hist_excl_WJetMult= bookHisto1D("njetWJet_excl", 7, -0.5, 6.5);
                _hist_Mult_exc = bookHisto1D("njet_exc_fbin", 7, -0.5, 6.5);
                
                //-------------
                _hist_JetRap1j =bookHisto1D("jet_Rap1jetcase", 12, 0, 2.4);
                _hist_JetRap2j =bookHisto1D("jet_Rap2jetcase", 12, 0, 2.4);
                _hist_JetRap3j =bookHisto1D("jet_Rap3jetcase", 8, 0., 2.4);
               
            }
        
            // define function used for filiing inc Njets histo
            void Fill(Histo1DPtr& _histJetMult, const double& weight, std::vector<FourMomentum>& finaljet_list){
                _histJetMult->fill(0, weight);
                for (size_t i=0 ; i<finaljet_list.size() ; ++i) {
                    if (i==6) break;
                    _histJetMult->fill(i+1, weight);  // inclusive multiplicity
                }
            }
            
            
            /// Perform the per-event analysis
            void analyze(const Event& event) {
                
                //cout << " beamIdpairs " << beamIds().first << " " << beamIds().second << " energy " << sqrtS() << endl;
                
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
                    
                    if ( (fabs(eta0) > 2.4) || (pt0 < 25.0*GeV) ) vetoEvent;
                    
                    // Obtain the jets.
                    vector<FourMomentum> finaljet_list;
                    double HT = 0.0;
                    
                    // loop over jets in an event, pushback in finaljet_list collection
                    foreach (const Jet& j, applyProjection<FastJets>(event, "Jets").jetsByPt(30.0*GeV)) {
                        const double jrap = j.momentum().rap();
                        const double jpt = j.momentum().pT();
                        
                        if ( (fabs(jrap) < 2.4) && (deltaR(lepton0, j.momentum()) > 0.4) ) {
                            
                            if(jpt > 30.0*GeV) {
                                finaljet_list.push_back(j.momentum());
                                HT += j.momentum().pT();
                            }
                        }
                    } // end looping over jets
                    
                    // Multiplicity exc plot.
                    _hist_excl_WJetMult->fill(finaljet_list.size(), weight);
                    
                    if(finaljet_list.size()<=7) {
                        _hist_Mult_exc->fill(finaljet_list.size(), weight);
                    }
                    else if (finaljet_list.size()>7){
                        _hist_Mult_exc->fill(7., weight);
                    }
                    
                    // Multiplicity inc plot.
                    Fill(_hist_inc_WJetMult, weight, finaljet_list);
                    
                    if(finaljet_list.size()>=1) {
                        _hist_JetPt1j->fill(finaljet_list[0].pT(), weight);
                        _hist_JetRap1j->fill(fabs(finaljet_list[0].rap()), weight);
                        _hist_Ht_1j->fill(HT, weight);
                    }
                    
                    if(finaljet_list.size()>=2) {
                        _hist_JetPt2j->fill(finaljet_list[1].pT(), weight);
                        _hist_JetRap2j->fill(fabs(finaljet_list[1].rap()), weight);
                        _hist_Ht_2j->fill(HT, weight);
                    }
                    
                    if(finaljet_list.size()>=3) {
                        _hist_JetPt3j->fill(finaljet_list[2].pT(), weight);
                        _hist_JetRap3j->fill(fabs(finaljet_list[2].rap()), weight);
                        _hist_Ht_3j->fill(HT, weight);
                    }
                } // W loop
                
            } // void loop
        
        
            /// Normalise histograms etc., after the run
            void finalize() {
                
                double crossSec=60290.0;
                
                scale(_hist_inc_WJetMult, crossSec/sumOfWeights());
                scale(_hist_excl_WJetMult, crossSec/sumOfWeights());
                scale(_hist_Mult_exc, crossSec/sumOfWeights());
                
                scale(_hist_JetPt1j, crossSec/sumOfWeights());
                scale(_hist_JetPt2j, crossSec/sumOfWeights());
                scale(_hist_JetPt3j, crossSec/sumOfWeights());
                
                scale(_hist_JetRap1j, crossSec/sumOfWeights());
                scale(_hist_JetRap2j, crossSec/sumOfWeights());
                scale(_hist_JetRap3j, crossSec/sumOfWeights());
                
                scale(_hist_Ht_1j, crossSec/sumOfWeights());
                scale(_hist_Ht_2j, crossSec/sumOfWeights());
                scale(_hist_Ht_3j, crossSec/sumOfWeights());
                
            }
        
        private:
            
            // Data members like post-cuts event weight counters go here
        
        private:
            
            Histo1DPtr _hist_inc_WJetMult;
            Histo1DPtr _hist_excl_WJetMult;
            Histo1DPtr _hist_Mult_exc;
            
            Histo1DPtr _hist_JetPt1j;
            Histo1DPtr _hist_JetPt2j;
            Histo1DPtr _hist_JetPt3j;
            
            Histo1DPtr _hist_JetRap1j;
            Histo1DPtr _hist_JetRap2j;
            Histo1DPtr _hist_JetRap3j;
        
            Histo1DPtr _hist_Ht_1j;
            Histo1DPtr _hist_Ht_2j;
            Histo1DPtr _hist_Ht_3j;
            //-------------------------------------

    };
    
    // This global object acts as a hook for the plugin system
    AnalysisBuilder<CMS_SMP_13TEV_WJETS> plugin_CMS_SMP_13TEV_WJETS;
    
}
