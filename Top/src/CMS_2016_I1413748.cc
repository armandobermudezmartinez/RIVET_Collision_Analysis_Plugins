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
      //IdentifiedFinalState l_id(fs);
      //l_id.acceptIdPair(PID::ELECTRON);
      //l_id.acceptIdPair(PID::MUON);
      //PromptFinalState leptons(l_id);
      //addProjection(leptons, "Leptons");
      //DressedLeptons dressedleptons(photons, leptons, 0.1, Cuts::open(), true, false);
      //addProjection(dressedleptons, "DressedLeptons");


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

      const std::vector<double> binntops = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5};
      _h_ntops = bookHisto1D("ntops", binntops);
      ntops_max = binntops[5];

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

      const std::vector<double> bindabsrapidity = {-2., -44./60., -20./60., 0., 20./60., 44./60., 2.};
      _h_dabsrapidity = bookHisto1D("dabsrapidity", bindabsrapidity);
      dabsrapidity_max = bindabsrapidity[6];

      const std::vector<double> bintop_costheta = {-1., -2./3., -1./3., 0., 1./3., 2./3., 1.};
      _h_top_costheta = bookHisto1D("top_costheta", bintop_costheta);
      top_costheta_max = bintop_costheta[6];

      const std::vector<double> binlepPlus_costheta = {-1., -2./3., -1./3., 0., 1./3., 2./3., 1.};
      _h_lepPlus_costheta = bookHisto1D("lepPlus_costheta", binlepPlus_costheta);
      lepPlus_costheta_max = binlepPlus_costheta[6];

      const std::vector<double> binlepMinus_costheta = {-1., -2./3., -1./3., 0., 1./3., 2./3., 1.};
      _h_lepMinus_costheta = bookHisto1D("lepMinus_costheta", binlepMinus_costheta);
      lepMinus_costheta_max = binlepMinus_costheta[6];

      const std::vector<double> binlep_costheta = {-1., -2./3., -1./3., 0., 1./3., 2./3., 1.};
      _h_lep_costheta = bookHisto1D("lep_costheta", binlep_costheta);
      lep_costheta_max = binlep_costheta[6];

      const std::vector<double> binc1c2 = {-1., -0.4, -10./60., 0., 10./60., 0.4, 1.};
      _h_c1c2 = bookHisto1D("c1c2", binc1c2);
      c1c2_max = binc1c2[6];

      const std::vector<double> bincos_opening_angle = {-1., -2./3., -1./3., 0., 1./3., 2./3., 1.};
      _h_cos_opening_angle = bookHisto1D("cos_opening_angle", bincos_opening_angle);
      cos_opening_angle_max = bincos_opening_angle[6];

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      double EPSILON = 0.00000001;
      
      const double weight = event.weight();
      
      // select ttbar -> lepton+jets without taus
      //const DressedLeptons& dressedleptons = applyProjection<DressedLeptons>(event, "DressedLeptons");
      //if (dressedleptons.dressedLeptons().size() != 2) vetoEvent;

      // select ttbar -> lepton+jets at particle level
      const DressedLeptons& dressed_electrons = applyProjection<DressedLeptons>(event, "DressedElectrons");
      const DressedLeptons& dressed_muons = applyProjection<DressedLeptons>(event, "DressedMuons");
      if (dressed_electrons.dressedLeptons().size() +
	  dressed_muons.dressedLeptons().size() != 2) {
	vetoEvent;
      }
  
      //emu only
      if (dressed_electrons.dressedLeptons().size() != 1 || dressed_muons.dressedLeptons().size() != 1) {
	vetoEvent;
      }


      if (sameSign(dressed_electrons.dressedLeptons()[0],dressed_muons.dressedLeptons()[0])) {                                                                                                            
        cout<<"error in lepton charge assignment"<<endl;    
      }    

      //FourMomentum lepPlus = dressedleptons.dressedLeptons()[0].charge() > 0 ? dressedleptons.dressedLeptons()[0].momentum() : dressedleptons.dressedLeptons()[1].momentum();
      //FourMomentum lepMinus = dressedleptons.dressedLeptons()[0].charge() > 0 ? dressedleptons.dressedLeptons()[1].momentum() : dressedleptons.dressedLeptons()[0].momentum();
     

      FourMomentum lepPlus = dressed_electrons.dressedLeptons()[0].charge() > 0 ? dressed_electrons.dressedLeptons()[0].momentum() : dressed_muons.dressedLeptons()[0].momentum();
      FourMomentum lepMinus = dressed_electrons.dressedLeptons()[0].charge() > 0 ? dressed_muons.dressedLeptons()[0].momentum() : dressed_electrons.dressedLeptons()[0].momentum();


 
      double dabseta_unbounded = lepPlus.abseta() - lepMinus.abseta();

      _h_dphi->fill(min(deltaPhi(lepPlus,lepMinus), dphi_max-EPSILON), weight);
      _h_dabseta->fill(dabseta_unbounded > 0 ? min(dabseta_unbounded, dabseta_max-EPSILON) : max(dabseta_unbounded, -dabseta_max+EPSILON), weight);


      // Get the lepton+jets ttbar candidate
      const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();
      if (leptonicpartontops.size() != 2) cout<<"error, leptonicpartontops.size() = "<<leptonicpartontops.size()<<endl;

      int ntops = leptonicpartontops.size();
      _h_ntops->fill(min(double(ntops),ntops_max-EPSILON), weight);


      // for the emu final state leptonicpartontops.size() always = 2, but explicitly requiring it here protects against crashes (just in case)
      if (leptonicpartontops.size() == 2) {
        // Get top quark 4-momenta
        FourMomentum topPlus_p4 = leptonicpartontops[0].pdgId() > 0 ? leptonicpartontops[0].momentum() : leptonicpartontops[1].momentum();
        FourMomentum topMinus_p4 = leptonicpartontops[0].pdgId() > 0 ? leptonicpartontops[1].momentum() : leptonicpartontops[0].momentum();

        if ( leptonicpartontops[0].pdgId()*leptonicpartontops[1].pdgId() > 0 ) {
          cout<<"error in top charge assignment"<<endl;
          _h_ntops->fill(-2., weight);
        }
        else {

          FourMomentum ttbar_p4 = topPlus_p4 + topMinus_p4;

          double tt_mass = (topPlus_p4 + topMinus_p4).mass();
          double tt_rapidity = (topPlus_p4 + topMinus_p4).rapidity();
          double tt_pT = ttbar_p4.pT();

          double dabsrapidity_unbounded = abs(topPlus_p4.rapidity()) - abs(topMinus_p4.rapidity());

          LorentzTransform ttCM;
          ttCM.setBetaVec(-ttbar_p4.boostVector());

          topPlus_p4 = ttCM.transform(topPlus_p4);
          topMinus_p4 = ttCM.transform(topMinus_p4);

          double top_costheta = topPlus_p4.vector3().dot(ttbar_p4.vector3()) / (topPlus_p4.vector3().mod() * ttbar_p4.vector3().mod());

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

          _h_tt_mass->fill(max(tt_mass_min,min(tt_mass,tt_mass_max-EPSILON)), weight);
          _h_tt_rapidity->fill(min(abs(tt_rapidity),tt_rapidity_max-EPSILON), weight);
          _h_tt_pT->fill(min(tt_pT,tt_pT_max-EPSILON), weight);
          _h_dabsrapidity->fill(dabsrapidity_unbounded > 0 ? min(dabsrapidity_unbounded, dabsrapidity_max-EPSILON) : max(dabsrapidity_unbounded, -dabsrapidity_max+EPSILON), weight);
          _h_top_costheta->fill(min(top_costheta,top_costheta_max-EPSILON), weight);
          _h_lepPlus_costheta->fill(min(lepPlus_costheta,lepPlus_costheta_max-EPSILON), weight);
          _h_lepMinus_costheta->fill(min(lepMinus_costheta,lepMinus_costheta_max-EPSILON), weight);
          _h_lep_costheta->fill(min(lepPlus_costheta,lepPlus_costheta_max-EPSILON), weight);
          _h_lep_costheta->fill(min(lepMinus_costheta,lepMinus_costheta_max-EPSILON), weight);
          _h_c1c2->fill(min(c1c2,c1c2_max-EPSILON), weight);
          _h_cos_opening_angle->fill(min(cos_opening_angle,cos_opening_angle_max-EPSILON), weight);

        }
      }

    }


    /// Normalise histograms to unit area
    void finalize() {

      normalize(_h_dphi);
      normalize(_h_dabseta);
      normalize(_h_ntops);
      normalize(_h_tt_mass);
      normalize(_h_tt_rapidity);
      normalize(_h_tt_pT);
      normalize(_h_dabsrapidity);
      normalize(_h_top_costheta);
      normalize(_h_lepPlus_costheta);
      normalize(_h_lepMinus_costheta);
      normalize(_h_lep_costheta);
      normalize(_h_c1c2);
      normalize(_h_cos_opening_angle);

    }


  private:
    double dphi_max, dabseta_max, ntops_max, tt_mass_min, tt_mass_max, tt_rapidity_max, tt_pT_max, dabsrapidity_max, top_costheta_max, lepPlus_costheta_max, lepMinus_costheta_max, lep_costheta_max, c1c2_max, cos_opening_angle_max;
    Histo1DPtr _h_dphi, _h_dabseta, _h_ntops, _h_tt_mass, _h_tt_rapidity, _h_tt_pT, _h_dabsrapidity, _h_top_costheta, _h_lepPlus_costheta, _h_lepMinus_costheta, _h_lep_costheta, _h_c1c2, _h_cos_opening_angle;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1413748);

}
