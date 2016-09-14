// -*- C++ -*-
// Author: F. Schaaf
// Created: 02Sep2014
// Modified by A. Descroix: 20Nov2014
// Modified by M.A. Harrendorf: January2015
// Modified by A. Descroix: Feb2015
// Implemented ghost b-tagging and dressed leptons, by A. Descroix: March2015

#include "Rivet/Analysis.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "Rivet/Tools/Logging.hh"


#include "HepMC/GenEvent.h"
#include <fstream>
#include <iostream>
#include <sstream>

namespace Rivet {


////////////
// Config //
////////////

// Lepton Cuts
const double MIN_PT_GOOD_ELECTRON = 30*GeV;
const double MIN_PT_GOOD_MUON = 30*GeV;
const double MIN_PT_LOOSE_ELECTRON = 15*GeV;
const double MIN_PT_LOOSE_MUON = 15*GeV;

const double MAX_ETA_GOOD_ELECTRON = 2.4;
const double MAX_ETA_GOOD_MUON = 2.4;
const double MAX_ETA_LOOSE_ELECTRON = 2.5;
const double MAX_ETA_LOOSE_MUON = 2.5;

// Jet Cuts
const double JET_MIN_PT = 30 * GeV;
const double JET_MAX_ETA = 2.5;
const double JET_MIN_DELTA_R = 0.5;

//Values up-to-date on the 25th of February 2015. Please check for updates!
const double TTbarXS = 252.89;						//From Top++ at 8 TeV with m_t=172.5 GeV/c^2. It's the official value from TOPLHCWG-> https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
const double MG5_BanchingRationCorrection = (0.108*9)*(0.676*1.5);	//Correction of the branching ratio (pdg: 0.108) which is wrongly modeled in MadGraph (taking 1/9). See page 147-148 of http://escholarship.org/uc/item/5f78t6s6

class CMS_2014_ttbarPlusXJets : public Analysis
{
public:

  CMS_2014_ttbarPlusXJets() : Analysis( "CMS_2014_ttbarPlusXJets" ), _nameLength(18), _nJets(20), _nJets_err(20), _nJetsE(20), _nJetsE_err(20), _nJetsMu(20), _nJetsMu_err(20)
  {
  }

private:

  enum Cuts
  {
    TotalEvents=0,
    GoodLepton,
    GoodLepton_Electron,
    GoodLepton_Muon,
    LooseLeptonVeto,
    LooseLeptonVeto_Electron,
    LooseLeptonVeto_Muon,
    NJetsGE4,
    NJetsGE4_Electron,
    NJetsGE4_Muon,
    NBJetsGE2,
    NBJetsGE2_Electron,
    NBJetsGE2_Muon,
    Size
  };

  char** _cutNames;

  size_t* _cuts;

#if DEBUG
  Log& _log = getLog();
#endif
  AIDA::IHistogram1D* _eventHisto;

  AIDA::IHistogram1D* _electronHisto;
  AIDA::IHistogram1D* _muonHisto;
  AIDA::IHistogram1D* _electronMuonHisto;
  
  AIDA::IHistogram1D* _normedElectronHisto;
  AIDA::IHistogram1D* _normedMuonHisto;
  AIDA::IHistogram1D* _normedElectronMuonHisto;
  
  AIDA::IHistogram1D* _absXSElectronHisto;
  AIDA::IHistogram1D* _absXSMuonHisto;
  AIDA::IHistogram1D* _absXSElectronMuonHisto;

private:

  const size_t _nameLength;

  void init_Names()
  {
	fjLepDef_ = std::shared_ptr<JetDef>(new JetDef(fastjet::antikt_algorithm, 0.1));
	fjJetDef_ = std::shared_ptr<JetDef>(new JetDef(fastjet::antikt_algorithm, 0.5));
	
    _nonOneWeights=false;
    _cuts = new size_t[Size];

        _cutNames = new char* [Size];

    for (size_t i = 0; i<Size; i++)
    {
            _cutNames[i] = new char[_nameLength];
    }
    sprintf(_cutNames[0], "TotalEvents       ");
    sprintf(_cutNames[1], "GoodLepton        ");
    sprintf(_cutNames[2], "GoodLepton|e      ");
    sprintf(_cutNames[3], "GoodLepton|mu     ");
    sprintf(_cutNames[4], "LooseLeptonVeto   ");
    sprintf(_cutNames[5], "LooseLeptonVeto|e ");
    sprintf(_cutNames[6], "LooseLeptonVeto|mu");
    sprintf(_cutNames[7], "nJets>=4          ");
    sprintf(_cutNames[8], "nJets>=4|e        ");
    sprintf(_cutNames[9], "nJets>=4|mu       ");
    sprintf(_cutNames[10],"nBJets>=2         ");
    sprintf(_cutNames[11],"nBJets>=2|e       ");
    sprintf(_cutNames[12],"nBJets>=2|mu      ");
  }

  void init_Projections()
  {
	const FinalState fState;
	const UnstableFinalState unstableFState;
	addProjection(FastJets(fState, FastJets::ANTIKT, 0.5), "Jets");
	addProjection( fState, "stableParticles" );
	addProjection( unstableFState, "unstableParticles" );
  }

  void init_Output()
  {
    for (size_t i = 0; i<Size; i++)
    {
            _cuts[i] = 0;
    }
  }
  

public:

  void init()
  {
    init_Names();
    init_Projections();
    init_Output();

    // Initialize histograms
    
    _eventHisto = bookHistogram1D("eventHisto", 1, 0, 1, "Inclusive Event Counter", "Event Count", "Number of events");
    
    _electronHisto = bookHistogram1D("electronHisto", 7, 4, 11, "Electron Jet Multiplicity", "Jet Multiplicity", "Number of events");
    _muonHisto = bookHistogram1D("muonHisto", 7, 4, 11, "Muon Jet Multiplicity", "Jet Multiplicity", "Number of events");
    _electronMuonHisto = bookHistogram1D("electronMuonHisto", 7, 4, 11, "Electron+Muon Jet Multiplicity", "Jet Multiplicity", "Number of events");
    
    _normedElectronHisto = bookHistogram1D("normedElectronHisto", 7, 4, 11, "Normalized Differential Cross Section in Electron+Jets Channel", "Jet Multiplicity", "Normed units");
    _normedMuonHisto = bookHistogram1D("normedMuonHisto", 7, 4, 11, "Normalized Differential Cross Section in Muon+Jets Channel", "Jet Multiplicity", "Normed units");
    _normedElectronMuonHisto = bookHistogram1D("normedElectronMuonHisto", 7, 4, 11, "Normalized Differential Cross Section in Lepton+Jets Channel", "Jet Multiplicity", "Normed units");
    
    _absXSElectronHisto = bookHistogram1D("absXSElectronHisto", 7, 4, 11, "Differential Cross Section in Electron+Jets Channel", "Jet Multiplicity", "pb");
    _absXSMuonHisto = bookHistogram1D("absXSMuonHisto", 7, 4, 11, "Differential Cross Section in Muon+Jets Channel", "Jet Multiplicity", "pb");
    _absXSElectronMuonHisto = bookHistogram1D("absXSElectronMuonHisto", 7, 4, 11, "Differential Cross Section in Lepton+Jets Channel", "Jet Multiplicity", "pb");

  }

private:
	
  std::vector<double> _nJets;
  std::vector<double> _nJets_err;
  std::vector<double> _nJetsE;
  std::vector<double> _nJetsE_err;
  std::vector<double> _nJetsMu;
  std::vector<double> _nJetsMu_err;
  
  typedef fastjet::JetDefinition JetDef;
  std::shared_ptr<JetDef> fjLepDef_;
  std::shared_ptr<JetDef> fjJetDef_;
  bool _nonOneWeights;

  bool isLepton( Particle const& p )
  {
    int const& absid = abs( p.pdgId() );
    return absid == ELECTRON || absid == MUON || absid == NU_E || absid == NU_MU;
  }

  enum ParticleTypeMatch
  {
    EXACT,
    ALLOW_SYMMETRIC_PARTNER
  };

  bool isParticleType( Particle const& p, int type, const ParticleTypeMatch particleTypeMatch = ALLOW_SYMMETRIC_PARTNER )
  {
    int pid = p.pdgId();
    if( particleTypeMatch == ALLOW_SYMMETRIC_PARTNER )
    {
      pid = abs( pid );
      type = abs( type );
    }
    return pid == type;
  }

  bool isE;
  bool isMu;

  bool analyze_Lepton( const Event& event, std::vector<Particle> leptons, Particle& chosenLepton, const double wgt)
  {
    vector<Particle> goodElectrons;
    vector<Particle> goodMuons;
    vector<Particle> looseElectrons;
    vector<Particle> looseMuons;
    foreach( const Particle p, leptons )
    {
      double pT = p.momentum().pT();
      double absEta = abs( p.momentum().eta() );
      if( isParticleType( p, ELECTRON ) && pT > MIN_PT_GOOD_ELECTRON && absEta < MAX_ETA_GOOD_ELECTRON )
      {
        goodElectrons.push_back( p );
      }
      else if( isParticleType( p, MUON ) && pT > MIN_PT_GOOD_MUON && absEta < MAX_ETA_GOOD_MUON )
      {
        goodMuons.push_back( p );
      }
      if( isParticleType( p, ELECTRON ) && pT > MIN_PT_LOOSE_ELECTRON && absEta < MAX_ETA_LOOSE_ELECTRON )
      {
        looseElectrons.push_back( p );
      }
      else if( isParticleType( p, MUON ) && pT > MIN_PT_LOOSE_MUON && absEta < MAX_ETA_LOOSE_MUON )
      {
        looseMuons.push_back( p );
      }
    }

    // Good Lepton Cut
    bool goodElectron = goodElectrons.size()==1 && goodMuons.size() == 0;
    bool goodMuon = goodMuons.size()==1 && goodElectrons.size() == 0;
    bool goodLepton = goodElectron | goodMuon;
    if (goodLepton)
    {
      _cuts[GoodLepton] += wgt;
    }

    if (goodElectron)
    {
      _cuts[GoodLepton_Electron] += wgt;
    }

    if (goodMuon)
    {
      _cuts[GoodLepton_Muon] += wgt;
    }

    // loose Lepton Veto
    isE = goodElectrons.size() == 1 && looseElectrons.size() == 1 && looseMuons.size() == 0;
    isMu = goodMuons.size() == 1 && looseElectrons.size() == 0 && looseMuons.size() == 1;
    bool isSelected = isE ^ isMu;

    if (!isSelected)
    {
      return false;
    }

    if( isE )
    {
      _cuts[LooseLeptonVeto_Electron] += wgt;
      chosenLepton = goodElectrons[0];
    }
    if( isMu )
    {
      _cuts[LooseLeptonVeto_Muon] += wgt;
      chosenLepton = goodMuons[0];
    }
    return true;
  }

  bool analyze_CountJets( const Event& event, std::vector<Particle> inputJets, std::vector<bool> bTags, const Particle chosenLepton, const double wgt,  const double sqWgt)
  {
    vector<FourMomentum> jets;
    int nBJets( 0 );
    int nJetsWrongEta( 0 );
    int nJetsTooCloseToGoodLepton(0);
    for( unsigned int i=0 ; i<inputJets.size() ; i++ )
    {
      Particle j=inputJets.at(i);
      //The pT cut is missing
      const double pT = j.momentum().pT();
      if( pT < JET_MIN_PT ) continue;
      
      const double absEta = abs( j.momentum().eta() );
      if( absEta > JET_MAX_ETA )
      {
    nJetsWrongEta++;
    continue;
      }
      if ( deltaR( chosenLepton.momentum(), j.momentum() ) < JET_MIN_DELTA_R )
      {
    nJetsTooCloseToGoodLepton++;
    continue;
      }

      jets.push_back( j.momentum() );
      if( bTags.at(i) )
      {
    nBJets++;
      }
    }

    int nJets = jets.size();

    if (nJets < 4)
    {
      return false;
    }

    _cuts[NJetsGE4] += wgt;
    if (isE) _cuts[NJetsGE4_Electron]  += wgt;
    if (isMu) _cuts[NJetsGE4_Muon]  += wgt;

    if (nBJets < 2)
    {
      return false;
    }
    
    _nJets[min(10,nJets)] += wgt;
    _electronMuonHisto->fill(min(10,nJets), wgt);
    _normedElectronMuonHisto->fill(min(10,nJets), wgt);
    _absXSElectronMuonHisto->fill(min(10,nJets), wgt);
    _nJets_err[min(10,nJets)] += sqWgt;
    
    if (isE) {
        _nJetsE[min(10,nJets)] += wgt;
        _electronHisto->fill(min(10,nJets), wgt);
        _normedElectronHisto->fill(min(10,nJets), wgt);
	_absXSElectronHisto->fill(min(10,nJets), wgt);
        _nJetsE_err[min(10,nJets)] += sqWgt;
    }
    
    if (isMu) {
        _nJetsMu[min(10,nJets)] += wgt;
        _muonHisto->fill(min(10,nJets), wgt);
        _normedMuonHisto->fill(min(10,nJets), wgt);
        _absXSMuonHisto->fill(min(10,nJets), wgt);
        _nJetsMu_err[min(10,nJets)] += sqWgt;
    }

    return true;
  }

public:

  void analyze( const Event& event )
  {
    const double wgt = event.weight();
    const double sqWgt = wgt*wgt;
    
    //Event counter, so it is accessible in the aida output
    _eventHisto->fill(0, wgt);
    
    if((!_nonOneWeights) && (wgt!=1)) _nonOneWeights=true;
    
    const FinalState fs = applyProjection<FinalState>( event, "stableParticles" );
    const UnstableFinalState ufs = applyProjection<UnstableFinalState>( event, "unstableParticles" );

    ParticleVector pVec = fs.particlesByPt();//stable particles
    ParticleVector uVec = ufs.particlesByPt();//unstable particles

    std::vector<Particle> leptons;
    std::vector<Particle> jets;
    
    int pID;
    std::vector<fastjet::PseudoJet> fjLepInputs;
    std::vector<size_t> neutrinoIdxs; // keep lepton constituents to remove from GenJet construction

    //Leptons:
    //Collect the lepton-clustering input objects and store the neutrino indexes to remove them from the jet clustering
    for( unsigned int i=0 ; i<pVec.size() ; i++ ) {
	    //Cleaning of the inputs as done by John seems to be not possible with RIVET_1.8.2 ...

	    Particle p = pVec.at(i);
	    if ( std::isnan(p.momentum().pT()) or p.momentum().pT() <= 0 ) continue;

	    pID=fabs(p.pdgId());

	    if((pID==11) || (pID==13) || (pID==22)){
		    fjLepInputs.push_back(fastjet::PseudoJet(p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().E())); //photons are also included for the dressing
		    fjLepInputs.back().set_user_index(i);
	    }					

	    //Storing neutrino indexes
	    if((pID==12) || (pID==14) || (pID==16)){
	    	//cout<<"Neutrino pT: "<<p.momentum().pT()<<" eta: "<<p.momentum().eta()<<endl;
	    	neutrinoIdxs.push_back(i);
	    }
    }

    // Run the jet algorithm for the leptons
    fastjet::ClusterSequence fjLepClusterSeq(fjLepInputs, *fjLepDef_);
    std::vector<fastjet::PseudoJet> fjLepJets = fastjet::sorted_by_pt(fjLepClusterSeq.inclusive_jets( 5.*GeV ));

    //// Build dressed lepton objects from the FJ output
    std::vector<size_t> lepDauIdxs;
    for ( auto& fjJet : fjLepJets )
    {
	    //Get jet constituents from fastJet
	    Particle lepCand;
	    if(lepCand.pdgId()!=0) cout<<"Warning, the pdgId of a new particle should be default 0. If not the program is faulty."<<endl;

	    foreach( const fastjet::PseudoJet& pJet, fastjet::sorted_by_pt(fjJet.constituents()) )
	    {
		    const size_t index = pJet.user_index();
		    Particle cand = pVec.at(index);

		    const int absPdgId = abs(cand.pdgId());
		    if ( absPdgId == 11 || absPdgId == 13 )
		    {
			    if ( lepCand.momentum().pT() > cand.momentum().pT() ) continue; // Choose one with highest pt
			    lepCand = cand;
		    }
	    }

	    //Central lepton must be the major component
	    if (( (lepCand.momentum().pT()>0.0) && (lepCand.momentum().pT() < fjJet.pt()/2) ) || (lepCand.pdgId()==0)) continue;

	    //Keep constituent indexes in order to cancel them from the RIVET-jet clustering input
	    foreach( const fastjet::PseudoJet& pJet, fastjet::sorted_by_pt(fjJet.constituents()) ) lepDauIdxs.push_back(pJet.user_index());

	    //Storing the lepton candidates found
	    FourMomentum mom(fjJet.E(),fjJet.px(),fjJet.py(),fjJet.pz());
	    Particle part(lepCand.pdgId(),mom);
	    leptons.push_back(part);
    }
	
    //for(unsigned int i=0 ; i<leptons.size() ; i++) cout<<"Lepton pT: "<<leptons.at(i).momentum().pT()<<" eta: "<<leptons.at(i).momentum().eta()<<endl;
    
    //Jets:
    ////Prepare input particle list.
    std::vector<fastjet::PseudoJet> fjJetInputs;
    std::vector<size_t> bHadronIdxs;
    for ( unsigned int i=0 ; i<pVec.size() ; i++ )
    {
    	const Particle p = pVec.at(i);
	if ( (std::isnan(p.momentum().pT())) || (p.momentum().pT() <= 0) ) continue;
	
	//Some quality cuts from John cannot be included
			
	//Remove neutrinos
	bool toSkip=false;
	for( unsigned int j=0 ; j<neutrinoIdxs.size() ; j++ ){
		if ( neutrinoIdxs.at(j) == i ){
			toSkip=true;
			continue;
		}
	}
	if(toSkip) continue;

	//Remove particles used in lepton clusters, 
	for( unsigned int j=0 ; j<lepDauIdxs.size() ; j++ ){
		if ( lepDauIdxs.at(j) == i ){
			toSkip=true;
			continue;
		}
	}
	if(toSkip) continue;
	
	fjJetInputs.push_back(fastjet::PseudoJet(p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().E()));
	fjJetInputs.back().set_user_index(-1);
    }
		
    //// Also don't forget to put B hadrons
    for ( unsigned int i=0 ; i<uVec.size() ; i++ )
    {
	const Particle p = uVec.at(i);
	int pid = p.pdgId();
	const double scale = 1.0E-20/p.momentum().E();

	//This is a reproduction of RIVET_2_1_2::HeavyHadrons.cc
	if( (!PID::isHadron(pid)) || (!PID::hasBottom(pid)) || (p.momentum().pT()<5.0*GeV) ) continue;

	// "An unbound, or undecayed status 2 hadron: this is weird, but I guess is allowed..."
	if (!p.hasGenParticle() || !p.genParticle().end_vertex()) {
		MSG_DEBUG("Heavy hadron " << pid << " with no GenParticle or decay found");
		fjJetInputs.push_back(fastjet::PseudoJet(p.momentum().px()*scale, p.momentum().py()*scale, p.momentum().pz()*scale, p.momentum().E()*scale));
		fjJetInputs.back().set_user_index(i);
		bHadronIdxs.push_back(i);
		continue;
	}

	//Alexis: The test whether the particle also decayed into a bottom hadron isn't available in 1.8.2 ... I guess it is not essential
	/*const vector<GenParticle*> children = particles_out((&p.genParticle()), HepMC::children);
	bool has_b_child = false;
	foreach (const GenParticle* p2, children) {
		if (p2.isHadron() && p2.hasBottom()) {
			has_b_child = true;
			break;
		}
	}
	if (!has_b_child) {*/
		fjJetInputs.push_back(fastjet::PseudoJet(p.momentum().px()*scale, p.momentum().py()*scale, p.momentum().pz()*scale, p.momentum().E()*scale));
		fjJetInputs.back().set_user_index(i);
		bHadronIdxs.push_back(i);
	/*}*/
    }
		
    //// Run the jet algorithm
    fastjet::ClusterSequence fjJetClusterSeq(fjJetInputs, *fjJetDef_);
    std::vector<fastjet::PseudoJet> fjJets = fastjet::sorted_by_pt(fjJetClusterSeq.inclusive_jets( 5.*GeV ));
    
    //// Build jets
    std::vector<bool> bTags;
    for ( auto& fjJet : fjJets )
    {
    	FourMomentum mom(fjJet.E(),fjJet.px(),fjJet.py(),fjJet.pz());
	Particle j(0,mom);
	jets.push_back(j);
	
	//cout<<"Jet pT: "<<fjJet.pt()<<" eta: "<<fjJet.eta()<<endl;
	
	//Check jet constituents
	bool hasBHadron = false;
	foreach( const fastjet::PseudoJet& pJet, fastjet::sorted_by_pt(fjJet.constituents()) )
	{
		if(pJet.user_index()!=-1){
			hasBHadron=true;
			continue;
		}
	}

	if( hasBHadron ){
		bTags.push_back(true);
		//cout<<"-> Tagged!"<<endl;
	}
	
	else bTags.push_back(false);
    }
    
    //cout<<endl;
    
    //Start the analysis itself
    
    _cuts[TotalEvents] += wgt;
    
    //leptonic cuts
    Particle chosenLepton;
    if( !analyze_Lepton( event, leptons, chosenLepton, wgt ) )
    {
      vetoEvent;
    }
    _cuts[LooseLeptonVeto] += wgt;

    //jet cuts
    if (!analyze_CountJets( event, jets, bTags, chosenLepton, wgt, sqWgt ))
    {
      vetoEvent;
    }
    _cuts[NBJetsGE2] += wgt;
    if (isE) _cuts[NBJetsGE2_Electron] += wgt;
    if (isMu) _cuts[NBJetsGE2_Muon] += wgt;
  }

  void finalize_Names()
  {

    for (size_t i = 0; i<Size; i++)
    {
            delete [] _cutNames[i];
    }
    delete [] _cutNames;
    delete [] _cuts;

}

  void finalize_DisplayCuts()
  {

    int count (-1);
    int countIndex( 0 );
    int childCount (0);

#if DEBUG
    _log << Log::ALWAYS << "[Cut flow with event weights]" << endl;
    _log << Log::ALWAYS << "--------------------------------------------------------" << endl;
    _log << Log::ALWAYS << "[#.#] (Cut Descr.|Branch):nEvents ( Parent% /   Total% )" << endl;
    _log << Log::ALWAYS << "--------------------------------------------------------" << endl;
#endif

    for( size_t i=0; i<Size; i++ )
    {
      double percentage = 100.*_cuts[i]/_cuts[TotalEvents];
      double ofParentPercentage = 100. * _cuts[i]/_cuts[countIndex];
      //_log.precision( 3 );

      bool isChild = false;
      if (std::string(_cutNames[i]).find('|') != std::string::npos)
      {
    isChild = true;
      }

      if (!isChild)
      {
    count++;
    countIndex = i;
      }

      stringstream ss;
      ss << "[" << count;

      if (isChild)
      {
        ss << "." << childCount;
      }
      else
      {
        ss << "  ";
      }
        ss << "] (" << _cutNames[i] << "): " << setw(6) << _cuts[i] << " (" << setw(7) << ofParentPercentage << "% / " << setw(7) << percentage << "% )" << endl;
      if (!isChild)
      {
    childCount=0;
      }
      else
      {
    childCount++;
      }
      #if DEBUG
      _log << Log::ALWAYS << ss.str();
      #endif
    }
    #if DEBUG
    _log << Log::ALWAYS << "--------------------------------------------------------" << endl;
    #endif
  }

  void finalize_DisplayJetMultiplicy()
  {
//    #if DEBUG
    cout <<  "[Jet Multiplicy in Merged Channel]" << endl;
    for (int i=4;i<=10;i++)
    {
      cout << "[" << i << "]: " << (long)_nJets[i] << " +- " << (long)sqrt(_nJets_err[i]) << endl;
    }
    
    cout << "[Jet Multiplicy in E+Jets Channel]" << endl;
    for (int i=4;i<=10;i++)
    {
      cout << "[" << i << "]: " << (long)_nJetsE[i] << " +- " << (long)sqrt(_nJetsE_err[i]) << endl;
    }
    
    cout << "[Jet Multiplicy in Mu+Jets Channel]" << endl;
    for (int i=4;i<=10;i++)
    {
      cout << "[" << i << "]: " <<(long)_nJetsMu[i] << " +- " << (long)sqrt(_nJetsMu_err[i]) << endl;
    }
//    #endif
  }

  void finalize_DisplayResults()
  {
    #if DEBUG
    _log << Log::ALWAYS << endl;
    _log << Log::ALWAYS << "[Displaying Results]" << endl;
    _log << Log::ALWAYS << endl;
    
    #endif
    finalize_DisplayCuts();
    finalize_DisplayJetMultiplicy();
  }

  void finalize()
  {
    #if DEBUG
    if(_nonOneWeights) _log<< Log::ALWAYS <<  endl<<"Warning! This sample has event weights different from one! Are you taking it into account?"<<endl<<endl;
    #endif
    
    finalize_DisplayResults();

    // Normalize histograms
    normalize(_normedElectronHisto, 1.);
    normalize(_normedMuonHisto, 1.);
    normalize(_normedElectronMuonHisto, 1.);

    const double conversionFactor = TTbarXS*MG5_BanchingRationCorrection/_cuts[TotalEvents];

    scale(_absXSElectronHisto, conversionFactor);
    scale(_absXSMuonHisto, conversionFactor);
    scale(_absXSElectronMuonHisto, conversionFactor);

  }



};

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2014_ttbarPlusXJets);

}
