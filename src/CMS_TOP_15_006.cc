// -*- C++ -*-
// Author: F. Schaaf
// Created: 02Sep2014
// Modified by A. Descroix: 20Nov2014
// Modified by M.A. Harrendorf: January2015
// Modified by A. Descroix: Feb2015
// Implemented ghost b-tagging and dressed leptons, by A. Descroix: March2015

#include "Rivet/Analysis.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/ParticleName.hh"
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

class CMS_TOP_15_006 : public Analysis
{
public:

  CMS_TOP_15_006() : Analysis( "CMS_TOP_15_006" )
  {
  }

private:

  Histo1DPtr _normedElectronMuonHisto;
  Histo1DPtr _absXSElectronMuonHisto;

public:

  void init()
  {
    fjLepDef_ = std::shared_ptr<JetDef>(new JetDef(fastjet::antikt_algorithm, 0.1));
    fjJetDef_ = std::shared_ptr<JetDef>(new JetDef(fastjet::antikt_algorithm, 0.5));

    const FinalState fState;
    const UnstableFinalState unstableFState;
    addProjection( fState, "stableParticles" );
    addProjection( unstableFState, "unstableParticles" );

    // Initialize histograms
    _normedElectronMuonHisto = bookHisto1D("normedElectronMuonHisto", 7, 3.5, 10.5, "Normalized Differential Cross Section in Lepton+Jets Channel", "Jet Multiplicity", "Normed units");
    _absXSElectronMuonHisto = bookHisto1D("absXSElectronMuonHisto", 7, 3.5, 10.5, "Differential Cross Section in Lepton+Jets Channel", "Jet Multiplicity", "pb");

  }

private:
  
  typedef fastjet::JetDefinition JetDef;
  std::shared_ptr<JetDef> fjLepDef_;
  std::shared_ptr<JetDef> fjJetDef_;

  bool isLepton( Particle const& p )
  {
    int const& absid = abs( p.pdgId() );
    return absid == PID::ELECTRON || absid == PID::MUON || absid == PID::NU_E || absid == PID::NU_MU;
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
      if( isParticleType( p, PID::ELECTRON ) && pT > MIN_PT_GOOD_ELECTRON && absEta < MAX_ETA_GOOD_ELECTRON )
      {
        goodElectrons.push_back( p );
      }
      else if( isParticleType( p, PID::MUON ) && pT > MIN_PT_GOOD_MUON && absEta < MAX_ETA_GOOD_MUON )
      {
        goodMuons.push_back( p );
      }
      if( isParticleType( p, PID::ELECTRON ) && pT > MIN_PT_LOOSE_ELECTRON && absEta < MAX_ETA_LOOSE_ELECTRON )
      {
        looseElectrons.push_back( p );
      }
      else if( isParticleType( p, PID::MUON ) && pT > MIN_PT_LOOSE_MUON && absEta < MAX_ETA_LOOSE_MUON )
      {
        looseMuons.push_back( p );
      }
    }

    // Good Lepton Cut and loose Lepton Veto
    isE = goodElectrons.size() == 1 && looseElectrons.size() == 1 && looseMuons.size() == 0;
    isMu = goodMuons.size() == 1 && looseElectrons.size() == 0 && looseMuons.size() == 1;
    bool isSelected = isE ^ isMu;

    if (!isSelected)
    {
      return false;
    }

    if( isE )
    {
      chosenLepton = goodElectrons[0];
    }
    if( isMu )
    {
      chosenLepton = goodMuons[0];
    }
    return true;
  }

  bool analyze_CountJets( const Event& event, std::vector<Particle> inputJets, std::vector<bool> bTags, const Particle chosenLepton, const double wgt)
  {
    vector<FourMomentum> jets;
    int nBJets( 0 );
    for( unsigned int i=0 ; i<inputJets.size() ; i++ )
    {
      Particle j=inputJets.at(i);
      //The pT cut is missing
      const double pT = j.momentum().pT();
      if( pT < JET_MIN_PT ) continue;
      
      const double absEta = abs( j.momentum().eta() );
      if( absEta > JET_MAX_ETA )
      {
        continue;
      }
      if ( deltaR( chosenLepton.momentum(), j.momentum() ) < JET_MIN_DELTA_R )
      {
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

    if (nBJets < 2)
    {
      return false;
    }
    
    _normedElectronMuonHisto->fill(min(10,nJets), wgt);
    _absXSElectronMuonHisto->fill(min(10,nJets), wgt);

    return true;
  }

public:

  void analyze( const Event& event )
  {
    const double wgt = event.weight();
    
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
	    //if (fjJet.pt() < MIN_PT_LOOSE_ELECTRON || abs(fjJet.eta()) > MAX_ETA_LOOSE_ELECTRON) continue; //TODO: Need feedback from Alexis
	    
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
    
    //std::cout << "-- Alexis --" << std::endl;
    //for(unsigned int i=0 ; i<leptons.size() ; i++) cout<<"lepton pT: "<<leptons.at(i).momentum().pT()<<" eta: "<<leptons.at(i).momentum().eta()<<" id: "<<leptons.at(i).pdgId()<<endl;
    
    //Jets:
    ////Prepare input particle list.
    std::vector<fastjet::PseudoJet> fjJetInputs;
    std::vector<size_t> bHadronIdxs;
    for ( unsigned int i=0 ; i<pVec.size() ; i++ )
    {
    	const Particle p = pVec.at(i);
	if ( (std::isnan(p.momentum().pT())) || (p.momentum().pT() <= 0) ) continue;
			
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

	fjJetInputs.push_back(fastjet::PseudoJet(p.momentum().px()*scale, p.momentum().py()*scale, p.momentum().pz()*scale, p.momentum().E()*scale));
	fjJetInputs.back().set_user_index(i);
	bHadronIdxs.push_back(i);
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
	}
	
	else bTags.push_back(false);
    }
    
    //Start the analysis itself
    
    //leptonic cuts
    Particle chosenLepton;
    if( !analyze_Lepton( event, leptons, chosenLepton, wgt ) )
    {
      vetoEvent;
    }

    //jet cuts
    if (!analyze_CountJets( event, jets, bTags, chosenLepton, wgt ))
    {
      vetoEvent;
    }
  }


  void finalize()
  {
    const double ttbarXS = !isnan(crossSectionPerEvent()) ? crossSection() : TTbarXS*picobarn;
    if (isnan(crossSectionPerEvent()))
      MSG_INFO("No valid cross-section given, using NNLO (arXiv:1303.6254; sqrt(s)=8 TeV, m_t=172.5 GeV): " << ttbarXS/picobarn << " pb");

    // Normalize histograms
    normalize(_normedElectronMuonHisto, 1.);

    const double xsPerWeight = ttbarXS/picobarn / sumOfWeights();
    scale(_absXSElectronMuonHisto, xsPerWeight);

  }



};

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_TOP_15_006);

}
