// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PartonicTops.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2016_I1413748 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2016_I1413748);


    /// Book histograms and initialise projections
    void init() {

      // Complete final state
      FinalState fs(-MAXDOUBLE, MAXDOUBLE, 0*GeV);

      // Projection for electrons and muons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      addProjection(electrons, "Electrons");
      DressedLeptons dressed_electrons(photons, electrons, 0.1, Cuts::open(), true, false);
      addProjection(dressed_electrons, "DressedElectrons");
  
      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      addProjection(muons, "Muons");
      DressedLeptons dressed_muons(photons, muons, 0.1, Cuts::open(), true, false);
      addProjection(dressed_muons, "DressedMuons");

      // Parton level top quarks
      declare(PartonicTops(PartonicTops::E_MU, false), "LeptonicPartonTops");

      // Booking of histograms



      const std::vector<double> bindphi = {0., 5.*M_PI/60., 10.*M_PI/60., 15.*M_PI/60., 20.*M_PI/60., 25.*M_PI/60., 30.*M_PI/60., 35.*M_PI/60., 40.*M_PI/60., 45.*M_PI/60., 50.*M_PI/60., 55.*M_PI/60., M_PI};
      _h_dphi = bookHisto1D("dphi", bindphi);
      dphi_max = bindphi[12];
      
      const std::vector<double> bindabseta = { -2., -68./60., -48./60., -32./60., -20./60., -8./60., 0., 8./60., 20./60., 32./60., 48./60., 68./60., 2.};
      _h_dabseta = bookHisto1D("dabseta", bindabseta);
      dabseta_max = bindabseta[12];


      //these two histograms are independent of parton-level information
      _h_dphidressedleptons = bookHisto1D("dphidressedleptons", bindphi);
      _h_dabsetadressedleptons = bookHisto1D("dabsetadressedleptons", bindabseta);


      const std::vector<double> bintt_mass = {300., 430., 530., 1200.};
      _h_tt_mass = bookHisto1D("tt_mass", bintt_mass);
      tt_mass_min = bintt_mass[0];
      tt_mass_max = bintt_mass[3];

      const std::vector<double> bintt_rapidity = {0., 0.34, 0.75, 1.5};
      _h_tt_rapidity = bookHisto1D("tt_rapidity", bintt_rapidity);
      tt_rapidity_max = bintt_rapidity[3];

      const std::vector<double> bintt_pT = {0., 41., 92., 300.};
      _h_tt_pT = bookHisto1D("tt_pT", bintt_pT);
      tt_pT_max = bintt_pT[3];

      //for each inclusive histogram, there are 9 further histograms in 3 bins each of ttbar invariant mass, pT, and absolute rapidity
      _h_dphi_bin[0][0] = bookHisto1D("dphi_mttbin1", bindphi);
      _h_dphi_bin[1][0] = bookHisto1D("dphi_ttptbin1", bindphi);
      _h_dphi_bin[2][0] = bookHisto1D("dphi_ttrapbin1", bindphi);
      _h_dphi_bin[0][1] = bookHisto1D("dphi_mttbin2", bindphi);
      _h_dphi_bin[1][1] = bookHisto1D("dphi_ttptbin2", bindphi);
      _h_dphi_bin[2][1] = bookHisto1D("dphi_ttrapbin2", bindphi);
      _h_dphi_bin[0][2] = bookHisto1D("dphi_mttbin3", bindphi);
      _h_dphi_bin[1][2] = bookHisto1D("dphi_ttptbin3", bindphi);
      _h_dphi_bin[2][2] = bookHisto1D("dphi_ttrapbin3", bindphi);

      _h_dabseta_bin[0][0] = bookHisto1D("dabseta_mttbin1", bindabseta);
      _h_dabseta_bin[1][0] = bookHisto1D("dabseta_ttptbin1", bindabseta);
      _h_dabseta_bin[2][0] = bookHisto1D("dabseta_ttrapbin1", bindabseta);
      _h_dabseta_bin[0][1] = bookHisto1D("dabseta_mttbin2", bindabseta);
      _h_dabseta_bin[1][1] = bookHisto1D("dabseta_ttptbin2", bindabseta);
      _h_dabseta_bin[2][1] = bookHisto1D("dabseta_ttrapbin2", bindabseta);
      _h_dabseta_bin[0][2] = bookHisto1D("dabseta_mttbin3", bindabseta);
      _h_dabseta_bin[1][2] = bookHisto1D("dabseta_ttptbin3", bindabseta);
      _h_dabseta_bin[2][2] = bookHisto1D("dabseta_ttrapbin3", bindabseta);

      const std::vector<double> bindabsrapidity = {-2., -44./60., -20./60., 0., 20./60., 44./60., 2.};
      _h_dabsrapidity = bookHisto1D("dabsrapidity", bindabsrapidity);
      _h_dabsrapidity_bin[0][0] = bookHisto1D("dabsrapidity_mttbin1", bindabsrapidity);
      _h_dabsrapidity_bin[1][0] = bookHisto1D("dabsrapidity_ttptbin1", bindabsrapidity);
      _h_dabsrapidity_bin[2][0] = bookHisto1D("dabsrapidity_ttrapbin1", bindabsrapidity);
      _h_dabsrapidity_bin[0][1] = bookHisto1D("dabsrapidity_mttbin2", bindabsrapidity);
      _h_dabsrapidity_bin[1][1] = bookHisto1D("dabsrapidity_ttptbin2", bindabsrapidity);
      _h_dabsrapidity_bin[2][1] = bookHisto1D("dabsrapidity_ttrapbin2", bindabsrapidity);
      _h_dabsrapidity_bin[0][2] = bookHisto1D("dabsrapidity_mttbin3", bindabsrapidity);
      _h_dabsrapidity_bin[1][2] = bookHisto1D("dabsrapidity_ttptbin3", bindabsrapidity);
      _h_dabsrapidity_bin[2][2] = bookHisto1D("dabsrapidity_ttrapbin3", bindabsrapidity);
      dabsrapidity_max = bindabsrapidity[6];

      const std::vector<double> binlep_costheta = {-1., -2./3., -1./3., 0., 1./3., 2./3., 1.};
      _h_lep_costheta = bookHisto1D("lep_costheta", binlep_costheta);
      _h_lep_costheta_bin[0][0] = bookHisto1D("lep_costheta_mttbin1", binlep_costheta);
      _h_lep_costheta_bin[1][0] = bookHisto1D("lep_costheta_ttptbin1", binlep_costheta);
      _h_lep_costheta_bin[2][0] = bookHisto1D("lep_costheta_ttrapbin1", binlep_costheta);
      _h_lep_costheta_bin[0][1] = bookHisto1D("lep_costheta_mttbin2", binlep_costheta);
      _h_lep_costheta_bin[1][1] = bookHisto1D("lep_costheta_ttptbin2", binlep_costheta);
      _h_lep_costheta_bin[2][1] = bookHisto1D("lep_costheta_ttrapbin2", binlep_costheta);
      _h_lep_costheta_bin[0][2] = bookHisto1D("lep_costheta_mttbin3", binlep_costheta);
      _h_lep_costheta_bin[1][2] = bookHisto1D("lep_costheta_ttptbin3", binlep_costheta);
      _h_lep_costheta_bin[2][2] = bookHisto1D("lep_costheta_ttrapbin3", binlep_costheta);
      lep_costheta_max = binlep_costheta[6];

      const std::vector<double> binlep_costheta_CPV = {-1., -2./3., -1./3., 0., 1./3., 2./3., 1.};
      _h_lep_costheta_CPV = bookHisto1D("lep_costheta_CPV", binlep_costheta_CPV);
      _h_lep_costheta_CPV_bin[0][0] = bookHisto1D("lep_costheta_CPV_mttbin1", binlep_costheta_CPV);
      _h_lep_costheta_CPV_bin[1][0] = bookHisto1D("lep_costheta_CPV_ttptbin1", binlep_costheta_CPV);
      _h_lep_costheta_CPV_bin[2][0] = bookHisto1D("lep_costheta_CPV_ttrapbin1", binlep_costheta_CPV);
      _h_lep_costheta_CPV_bin[0][1] = bookHisto1D("lep_costheta_CPV_mttbin2", binlep_costheta_CPV);
      _h_lep_costheta_CPV_bin[1][1] = bookHisto1D("lep_costheta_CPV_ttptbin2", binlep_costheta_CPV);
      _h_lep_costheta_CPV_bin[2][1] = bookHisto1D("lep_costheta_CPV_ttrapbin2", binlep_costheta_CPV);
      _h_lep_costheta_CPV_bin[0][2] = bookHisto1D("lep_costheta_CPV_mttbin3", binlep_costheta_CPV);
      _h_lep_costheta_CPV_bin[1][2] = bookHisto1D("lep_costheta_CPV_ttptbin3", binlep_costheta_CPV);
      _h_lep_costheta_CPV_bin[2][2] = bookHisto1D("lep_costheta_CPV_ttrapbin3", binlep_costheta_CPV);

      const std::vector<double> binc1c2 = {-1., -0.4, -10./60., 0., 10./60., 0.4, 1.};
      _h_c1c2 = bookHisto1D("c1c2", binc1c2);
      _h_c1c2_bin[0][0] = bookHisto1D("c1c2_mttbin1", binc1c2);
      _h_c1c2_bin[1][0] = bookHisto1D("c1c2_ttptbin1", binc1c2);
      _h_c1c2_bin[2][0] = bookHisto1D("c1c2_ttrapbin1", binc1c2);
      _h_c1c2_bin[0][1] = bookHisto1D("c1c2_mttbin2", binc1c2);
      _h_c1c2_bin[1][1] = bookHisto1D("c1c2_ttptbin2", binc1c2);
      _h_c1c2_bin[2][1] = bookHisto1D("c1c2_ttrapbin2", binc1c2);
      _h_c1c2_bin[0][2] = bookHisto1D("c1c2_mttbin3", binc1c2);
      _h_c1c2_bin[1][2] = bookHisto1D("c1c2_ttptbin3", binc1c2);
      _h_c1c2_bin[2][2] = bookHisto1D("c1c2_ttrapbin3", binc1c2);
      c1c2_max = binc1c2[6];

      const std::vector<double> bincos_opening_angle = {-1., -2./3., -1./3., 0., 1./3., 2./3., 1.};
      _h_cos_opening_angle = bookHisto1D("cos_opening_angle", bincos_opening_angle);
      _h_cos_opening_angle_bin[0][0] = bookHisto1D("cos_opening_angle_mttbin1", bincos_opening_angle);
      _h_cos_opening_angle_bin[1][0] = bookHisto1D("cos_opening_angle_ttptbin1", bincos_opening_angle);
      _h_cos_opening_angle_bin[2][0] = bookHisto1D("cos_opening_angle_ttrapbin1", bincos_opening_angle);
      _h_cos_opening_angle_bin[0][1] = bookHisto1D("cos_opening_angle_mttbin2", bincos_opening_angle);
      _h_cos_opening_angle_bin[1][1] = bookHisto1D("cos_opening_angle_ttptbin2", bincos_opening_angle);
      _h_cos_opening_angle_bin[2][1] = bookHisto1D("cos_opening_angle_ttrapbin2", bincos_opening_angle);
      _h_cos_opening_angle_bin[0][2] = bookHisto1D("cos_opening_angle_mttbin3", bincos_opening_angle);
      _h_cos_opening_angle_bin[1][2] = bookHisto1D("cos_opening_angle_ttptbin3", bincos_opening_angle);
      _h_cos_opening_angle_bin[2][2] = bookHisto1D("cos_opening_angle_ttrapbin3", bincos_opening_angle);
      cos_opening_angle_max = bincos_opening_angle[6];

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      double EPSILON = 0.00000001;
      
      const double weight = event.weight();

      // use particle-level leptons for the first 2 histos
      const DressedLeptons& dressed_electrons = applyProjection<DressedLeptons>(event, "DressedElectrons");
      const DressedLeptons& dressed_muons = applyProjection<DressedLeptons>(event, "DressedMuons");

      const std::vector<DressedLepton> delectrons = dressed_electrons.dressedLeptons();
      const std::vector<DressedLepton> dmuons = dressed_muons.dressedLeptons();

      int ndressede = delectrons.size();
      int ndressedm = dmuons.size();

      //For the particle-level histos, require an odd number of electrons and an odd number of muons, to select ttbar->emu channel. This means we can easily identify additional dilepton pairs from the shower (which will be same-flavour with invariant mass M_ll~0 GeV ), which distort the distributions relative to those used in the analysis where the leptons are prompt leptons from top decay only.
      if ( ndressede % 2 == 1 && ndressedm % 2 == 1  ) {

        int electrontouse = 0;
        int muontouse = 0;

        if( ndressede > 1 ) electrontouse = returnPrimaryLepton(delectrons);
        if( ndressedm > 1 ) muontouse = returnPrimaryLepton(dmuons);
  
        //Get the four-momenta of the positively- and negatively-charged leptons
        FourMomentum lepPlus = delectrons[electrontouse].charge() > 0 ? delectrons[electrontouse].momentum() : dmuons[muontouse].momentum();
        FourMomentum lepMinus = delectrons[electrontouse].charge() > 0 ? dmuons[muontouse].momentum() : delectrons[electrontouse].momentum();

        //now calculate the variables

        double dabseta_unbounded = lepPlus.abseta() - lepMinus.abseta();

        double _dphi = min(deltaPhi(lepPlus,lepMinus), dphi_max-EPSILON);
        double _dabseta = dabseta_unbounded > 0 ? min(dabseta_unbounded, dabseta_max-EPSILON) : max(dabseta_unbounded, -dabseta_max+EPSILON);

        //fill particle-level histos for opposite-charge leptons only
        if (sameSign(delectrons[electrontouse],dmuons[muontouse])) cout<<"error, emu same charge"<<endl;
        else{
          _h_dphidressedleptons->fill(_dphi, weight);
          _h_dabsetadressedleptons->fill(_dabseta, weight);
        }
      }

      //the remaining variables use parton-level information

      // Get the leptonically decaying tops
      const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();
      Particles chargedleptons;

      unsigned int ntrueleptonictops = 0;
      bool oppositesign = false;
  
      if (leptonicpartontops.size() == 2) {
        for (unsigned int k = 0; k < leptonicpartontops.size(); ++k) {
    
          //get the lepton
          const Particle lepTop = leptonicpartontops[k];
          //const auto isPromptChargedLepton = [](const Particle& p){return (isChargedLepton(p) && !fromDecay(p));}; //this works for PYTHIA but not MC@NLO+herwig
          //const auto isPromptChargedLepton = [](const Particle& p){return (isChargedLepton(p) && !fromHadron(p));}; //this works for PYTHIA but not MC@NLO+herwig
          //const auto isPromptChargedLepton = [](const Particle& p){return (isChargedLepton(p) && p.genParticle()->status() == 3 );}; //this works for MC@NLO+herwig in combination with allDescendants,false
          const auto isPromptChargedLepton = [](const Particle& p){return (isChargedLepton(p) && isPrompt(p, false, false));}; //this works for MC@NLO+herwig in combination with allDescendants,false
          Particles lepton_candidates = lepTop.allDescendants(firstParticleWith(isPromptChargedLepton), false); 
          //Particles lepton_candidates = lepTop.allDescendants(lastParticleWith(isPromptChargedLepton), false); //not all leptons with gamma parents are identified in this case
          if( lepton_candidates.size() < 1 ) cout<<"error, leptonicpartontop gives no lepton_candidates: "<< lepton_candidates.size() <<endl;
          bool istrueleptonictop = false;
          
          // In some cases there is no lepton from the W decay but only leptons from a decay of a radiated gamma. 
          // These are mistakenly identified by LeptonicPartonTops (as of April 2017), and need to be rejected.
          // LeptonicPartonTops is being fixed in Rivet, and when it is the veto below should do nothing.
          for (unsigned int i = 0; i < lepton_candidates.size(); ++i) {
            Particle lepton_candidate = lepton_candidates.at(i);
            if ( hasParentWith(lepton_candidate, Cuts::abspid == 22)) {cout<<"found gamma parent, top: "<<k+1<<" of "<<leptonicpartontops.size()<<" , lepton: "<<i+1<<" of "<<lepton_candidates.size()<<endl; continue;}
            if(!istrueleptonictop && sameSign(lepTop,lepton_candidate) ) {
              chargedleptons.push_back(lepton_candidate);
              istrueleptonictop = true;
            }
            else cout<<"error, found extra lepton"<<endl;
          }
          if(istrueleptonictop) ++ntrueleptonictops;
        }
      }

      if(ntrueleptonictops != chargedleptons.size()) cout<<"error in lepton count"<<endl;

      if( ntrueleptonictops == 2 ) {
        oppositesign = !( sameSign(chargedleptons[0],chargedleptons[1]) );
        if(!oppositesign) cout<<"error, same charge tops"<<endl;
      }

      if ( ntrueleptonictops == 2 && oppositesign ) {
 
        //Get the four-momenta of the positively- and negatively-charged leptons
        FourMomentum lepPlus = chargedleptons[0].charge() > 0 ? chargedleptons[0].momentum() : chargedleptons[1].momentum();
        FourMomentum lepMinus = chargedleptons[0].charge() > 0 ? chargedleptons[1].momentum() : chargedleptons[0].momentum();
 
        double dabseta_unbounded = lepPlus.abseta() - lepMinus.abseta();

        double _dphi = min(deltaPhi(lepPlus,lepMinus), dphi_max-EPSILON);
        double _dabseta = dabseta_unbounded > 0 ? min(dabseta_unbounded, dabseta_max-EPSILON) : max(dabseta_unbounded, -dabseta_max+EPSILON);

        // Get the top quark 4-momenta
        FourMomentum topPlus_p4 = leptonicpartontops[0].pdgId() > 0 ? leptonicpartontops[0].momentum() : leptonicpartontops[1].momentum();
        FourMomentum topMinus_p4 = leptonicpartontops[0].pdgId() > 0 ? leptonicpartontops[1].momentum() : leptonicpartontops[0].momentum();

        FourMomentum ttbar_p4 = topPlus_p4 + topMinus_p4;

        double tt_mass = (topPlus_p4 + topMinus_p4).mass();
        double tt_rapidity = (topPlus_p4 + topMinus_p4).rapidity();
        double tt_pT = ttbar_p4.pT();

        double dabsrapidity_unbounded = abs(topPlus_p4.rapidity()) - abs(topMinus_p4.rapidity());

        LorentzTransform ttCM;
        ttCM.setBetaVec(-ttbar_p4.boostVector());

        topPlus_p4 = ttCM.transform(topPlus_p4);
        topMinus_p4 = ttCM.transform(topMinus_p4);

        lepPlus = ttCM.transform(lepPlus);
        lepMinus = ttCM.transform(lepMinus);

        LorentzTransform topPlus, topMinus;
        topPlus.setBetaVec(-topPlus_p4.boostVector());
        topMinus.setBetaVec(-topMinus_p4.boostVector());

        lepPlus = topPlus.transform(lepPlus);
        lepMinus = topMinus.transform(lepMinus);

        double lepPlus_costheta = lepPlus.vector3().dot(topPlus_p4.vector3()) / (lepPlus.vector3().mod() * topPlus_p4.vector3().mod());
        double lepMinus_costheta = lepMinus.vector3().dot(topMinus_p4.vector3()) / (lepMinus.vector3().mod() * topMinus_p4.vector3().mod());
        double c1c2 = lepPlus_costheta * lepMinus_costheta;
        double cos_opening_angle = lepPlus.vector3().dot(lepMinus.vector3()) / (lepPlus.vector3().mod() * lepMinus.vector3().mod());

        //fill parton-level histos
        double _tt_mass = max(tt_mass_min,min(tt_mass,tt_mass_max-EPSILON));
        double _tt_rapidity = min(abs(tt_rapidity),tt_rapidity_max-EPSILON);
        double _tt_pT = min(tt_pT,tt_pT_max-EPSILON);
        double _dabsrapidity = dabsrapidity_unbounded > 0 ? min(dabsrapidity_unbounded, dabsrapidity_max-EPSILON) : max(dabsrapidity_unbounded, -dabsrapidity_max+EPSILON);
        double _lepPlus_costheta = min(lepPlus_costheta,lep_costheta_max-EPSILON);
        double _lepMinus_costheta = min(lepMinus_costheta,lep_costheta_max-EPSILON);
        double _c1c2 = min(c1c2,c1c2_max-EPSILON);
        double _cos_opening_angle = min(cos_opening_angle,cos_opening_angle_max-EPSILON);

        _h_tt_mass->fill(_tt_mass, weight);
        _h_tt_rapidity->fill(_tt_rapidity, weight);
        _h_tt_pT->fill(_tt_pT, weight);
        _h_dphi->fill(_dphi, weight);
        _h_dabseta->fill(_dabseta, weight);
        _h_dabsrapidity->fill(_dabsrapidity, weight);
        _h_lep_costheta->fill(_lepPlus_costheta, weight);
        _h_lep_costheta->fill(_lepMinus_costheta, weight);
        _h_lep_costheta_CPV->fill(_lepPlus_costheta, weight);
        _h_lep_costheta_CPV->fill(-_lepMinus_costheta, weight);
        _h_c1c2->fill(_c1c2, weight);
        _h_cos_opening_angle->fill(_cos_opening_angle, weight);

        //now fill the same variables in each of their 3 bins of ttbar invariant mass, pT, and absolute rapidity
        for (int i_var = 0; i_var < 3; ++i_var)
        {
          double var;
          std::vector<double> bins_var;

          if(i_var == 0) {
            var = _tt_mass;
            bins_var = {300., 430., 530., 1200.};
          }
          else if(i_var == 1){
            var = _tt_pT;
            bins_var = {0., 41., 92., 300.};
          }
          else {
            var = _tt_rapidity;
            bins_var = {0., 0.34, 0.75, 1.5};
          }

          int j_bin = -1;

          if(var < bins_var[1]) j_bin = 0;
          else if(var < bins_var[2]) j_bin = 1;
          else j_bin = 2;

          _h_dphi_bin[i_var][j_bin]->fill(_dphi, weight);
          _h_dabseta_bin[i_var][j_bin]->fill(_dabseta, weight);
          _h_dabsrapidity_bin[i_var][j_bin]->fill(_dabsrapidity, weight);
          _h_lep_costheta_bin[i_var][j_bin]->fill(_lepPlus_costheta, weight);
          _h_lep_costheta_bin[i_var][j_bin]->fill(_lepMinus_costheta, weight);
          _h_lep_costheta_CPV_bin[i_var][j_bin]->fill(_lepPlus_costheta, weight);
          _h_lep_costheta_CPV_bin[i_var][j_bin]->fill(-_lepMinus_costheta, weight);
          _h_c1c2_bin[i_var][j_bin]->fill(_c1c2, weight);
          _h_cos_opening_angle_bin[i_var][j_bin]->fill(_cos_opening_angle, weight);

        }

      }

    }


    /// Normalise histograms to unit area
    void finalize() {

      normalize(_h_dphidressedleptons);
      normalize(_h_dabsetadressedleptons);
      normalize(_h_dphi);
      normalize(_h_dabseta);
      normalize(_h_tt_mass);
      normalize(_h_tt_rapidity);
      normalize(_h_tt_pT);
      normalize(_h_dabsrapidity);
      normalize(_h_lep_costheta);
      normalize(_h_lep_costheta_CPV);
      normalize(_h_c1c2);
      normalize(_h_cos_opening_angle);

      for (int i_var = 0; i_var < 3; ++i_var)
      {
        for (int j_bin = 0; j_bin < 3; ++j_bin)
        {
          normalize(_h_dphi_bin[i_var][j_bin]);
          normalize(_h_dabseta_bin[i_var][j_bin]);
          normalize(_h_dabsrapidity_bin[i_var][j_bin]);
          normalize(_h_lep_costheta_bin[i_var][j_bin]);
          normalize(_h_lep_costheta_CPV_bin[i_var][j_bin]);
          normalize(_h_c1c2_bin[i_var][j_bin]);
          normalize(_h_cos_opening_angle_bin[i_var][j_bin]);
        }
      }

    }


  private:
    double dphi_max, dabseta_max, tt_mass_min, tt_mass_max, tt_rapidity_max, tt_pT_max, dabsrapidity_max, lep_costheta_max, c1c2_max, cos_opening_angle_max;
    Histo1DPtr _h_dphidressedleptons, _h_dabsetadressedleptons, _h_dphi, _h_dabseta, _h_tt_mass, _h_tt_rapidity, _h_tt_pT, _h_dabsrapidity, _h_lep_costheta, _h_lep_costheta_CPV, _h_c1c2, _h_cos_opening_angle;
    Histo1DPtr _h_dphi_bin[3][3], _h_dabseta_bin[3][3], _h_dabsrapidity_bin[3][3], _h_lep_costheta_bin[3][3], _h_lep_costheta_CPV_bin[3][3], _h_c1c2_bin[3][3], _h_cos_opening_angle_bin[3][3];

    struct ilepsmll{
        double mll;
        int ilep1;
        int ilep2;
    };

    typedef vector< ilepsmll > Vofilepsmll;

    //inline bool sortbymll(ilepsmll ilepsmll1, ilepsmll ilepsmll2) {
    //    return ilepsmll1.mll > ilepsmll2.mll;
    //}

    int returnPrimaryLepton(const std::vector<DressedLepton> dleptons) {
      //when there are additional lepton pair(s) from the shower, remove the combination of same-flavour opposite-charge lepton pair(s) that has the smallest maximum invariant mass m_ll

      int ndressedlep = dleptons.size();

      if(ndressedlep < 3 || ndressedlep % 2 == 0) {
        cout<<"error, returnPrimaryLepton only works with an odd number of 3 or more leptons, returning 0."<<endl;
        return 0;
      }

      int leptontouse = 0; //variable that will be updated with the ordinal number of the primary lepton

      int nuniquesets = 1; //counter of unique sets of lepton pairs. For 2M+1 dressed leptons, there will be nuniquesets = M!(M+1)! unique sets of pairs. Probability of M>3 negligible (i.e. >3 extra lepton pairs from radiation), so create arrays of size 3!*4! = 144.
      vector<int> leptonsused[144]; //keep track of lepton numbers used forming each unique set of pairs
      Vofilepsmll vimlls[144]; //vector of information about each pair for each unique set of lepton pairs
      //clear the vector for the first unique set of lepton pairs. The others will be cleared later if they are used.
      vimlls[0].clear();
      leptonsused[0].clear();

      if(ndressedlep > 7) {
        cout<<"error, found "<<ndressedlep<<" leptons. returnPrimaryLepton is only configured to use up to 7 leptons. Returning 0."<<endl;
        return 0;
      }

      for (int i = 0; i < ndressedlep; ++i)
      {
        for (int j = i+1; j < ndressedlep; ++j)
        {
          if (!sameSign(dleptons[i],dleptons[j]))
          {
            double mll = (dleptons[i].momentum() + dleptons[j].momentum()).mass();
            ilepsmll imll = { mll, dleptons[i].charge() > 0 ? i : j, dleptons[i].charge() > 0 ? j : i };

            //fill unique sets of lepton pairs
            bool leptonsfilled = false;
            for (int uniqueset = 0; uniqueset < nuniquesets; ++uniqueset)
            {
              //fill leptons into this set if neither have already been filled
              if ( !(std::find(leptonsused[uniqueset].begin(), leptonsused[uniqueset].end(), i) != leptonsused[uniqueset].end()) && !(std::find(leptonsused[uniqueset].begin(), leptonsused[uniqueset].end(), j) != leptonsused[uniqueset].end()) )
              {
                leptonsused[uniqueset].push_back(i);
                leptonsused[uniqueset].push_back(j);
                vimlls[uniqueset].push_back(imll);
                leptonsfilled = true;
              }

            }
            //if the leptons can't be filled in an existing set create a new one
            if (!leptonsfilled){
              leptonsused[nuniquesets].clear();
              leptonsused[nuniquesets].push_back(i);
              leptonsused[nuniquesets].push_back(j);
              vimlls[nuniquesets].clear();
              vimlls[nuniquesets].push_back(imll);
              leptonsfilled = true;  
              ++nuniquesets;
            } 

          }
        }
      }
      if ( vimlls[0].size() == 0 ) cout<<"error, all same charge"<<endl;

      cout<<"nuniquesets count check: "<<nuniquesets<<" ndressedlep: "<<ndressedlep<<endl;


      for (int uniqueset = 0; uniqueset < nuniquesets; ++uniqueset)
      {
        sort(vimlls[uniqueset].begin(), vimlls[uniqueset].end(), [](ilepsmll ilepsmll1, ilepsmll ilepsmll2) { return ilepsmll1.mll > ilepsmll2.mll; } );
      }

      //which set has the lowest mass for the last pair to be removed?
      int lowestmassuniqueset = -1;
      double lowestmass = 99999;
      for (int uniqueset = 0; uniqueset < nuniquesets; ++uniqueset)
      {
        if (vimlls[uniqueset][0].mll < lowestmass) {
          lowestmass = vimlls[uniqueset][0].mll;
          lowestmassuniqueset = uniqueset;
        }
      }



      cout<<"lowest mass unique set is "<<lowestmassuniqueset<<" :"<<endl;
      for (unsigned int i = 0; i < vimlls[lowestmassuniqueset].size(); ++i)
      {
        cout<<vimlls[lowestmassuniqueset][i].mll<<" "<<vimlls[lowestmassuniqueset][i].ilep1<<" "<<vimlls[lowestmassuniqueset][i].ilep2<<endl;
      }
      cout<<leptonsused[lowestmassuniqueset]<<endl;

      while ( (std::find(leptonsused[lowestmassuniqueset].begin(), leptonsused[lowestmassuniqueset].end(), leptontouse) != leptonsused[lowestmassuniqueset].end())  ) ++leptontouse;
      cout<<"Using lepton "<<leptontouse<<endl;
      cout<<endl;

      if(leptontouse >= ndressedlep) {
        cout<<"error, returnPrimaryLepton couldn't find an unpaired lepton, returning 0."<<endl;
        return 0;
      }

      return leptontouse;

    }






  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1413748);

}
