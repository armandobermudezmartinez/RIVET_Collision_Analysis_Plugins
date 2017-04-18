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

      //the first two histograms are independent of parton-level information

      const std::vector<double> bindphi = {0., 5.*M_PI/60., 10.*M_PI/60., 15.*M_PI/60., 20.*M_PI/60., 25.*M_PI/60., 30.*M_PI/60., 35.*M_PI/60., 40.*M_PI/60., 45.*M_PI/60., 50.*M_PI/60., 55.*M_PI/60., M_PI};
      _h_dphi = bookHisto1D("dphi", bindphi);
      dphi_max = bindphi[12];
      
      const std::vector<double> bindabseta = { -2., -68./60., -48./60., -32./60., -20./60., -8./60., 0., 8./60., 20./60., 32./60., 48./60., 68./60., 2.};
      _h_dabseta = bookHisto1D("dabseta", bindabseta);
      dabseta_max = bindabseta[12];

      //the remaining histograms use parton-level information

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

      // select ttbar -> dileptons at particle level
      const DressedLeptons& dressed_electrons = applyProjection<DressedLeptons>(event, "DressedElectrons");
      const DressedLeptons& dressed_muons = applyProjection<DressedLeptons>(event, "DressedMuons");

      if (dressed_electrons.dressedLeptons().size() + dressed_muons.dressedLeptons().size() != 2) vetoEvent;
  
      //emu channel only. This additional selection is made because in same-flavour ttbar events (particularly ee) generated by PYTHIA the leptons can come from the decay of low-mass resonances rather than from ttbar, so the leptons have invariant mass M_ll~0 which distorts the distributions relative to those used in the analysis where the leptons are prompt leptons from top decay only.
      if (dressed_electrons.dressedLeptons().size() != 1 || dressed_muons.dressedLeptons().size() != 1) vetoEvent;

      //opposite-charge leptons only
      if (sameSign(dressed_electrons.dressedLeptons()[0],dressed_muons.dressedLeptons()[0])) vetoEvent;

      //Get the four-momenta of the positively- and negatively-charged leptons
      FourMomentum lepPlus = dressed_electrons.dressedLeptons()[0].charge() > 0 ? dressed_electrons.dressedLeptons()[0].momentum() : dressed_muons.dressedLeptons()[0].momentum();
      FourMomentum lepMinus = dressed_electrons.dressedLeptons()[0].charge() > 0 ? dressed_muons.dressedLeptons()[0].momentum() : dressed_electrons.dressedLeptons()[0].momentum();

      //now calculate the variables

      //the first two variables are independent of parton-level information

      double dabseta_unbounded = lepPlus.abseta() - lepMinus.abseta();

      //fill particle-level histos
      double _dphi = min(deltaPhi(lepPlus,lepMinus), dphi_max-EPSILON);
      double _dabseta = dabseta_unbounded > 0 ? min(dabseta_unbounded, dabseta_max-EPSILON) : max(dabseta_unbounded, -dabseta_max+EPSILON);

      _h_dphi->fill(_dphi, weight);
      _h_dabseta->fill(_dabseta, weight);

      //the remaining variables use parton-level information

      // Get the leptonically decaying tops
      const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();

      // for the emu final state leptonicpartontops.size() always = 2, but explicitly requiring it here protects against crashes (just in case)
      if (leptonicpartontops.size() == 2 && leptonicpartontops[0].pdgId()*leptonicpartontops[1].pdgId() < 0) {
        
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
    Histo1DPtr _h_dphi, _h_dabseta, _h_tt_mass, _h_tt_rapidity, _h_tt_pT, _h_dabsrapidity, _h_lep_costheta, _h_lep_costheta_CPV, _h_c1c2, _h_cos_opening_angle;
    Histo1DPtr _h_dphi_bin[3][3], _h_dabseta_bin[3][3], _h_dabsrapidity_bin[3][3], _h_lep_costheta_bin[3][3], _h_lep_costheta_CPV_bin[3][3], _h_c1c2_bin[3][3], _h_cos_opening_angle_bin[3][3];

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1413748);

}
