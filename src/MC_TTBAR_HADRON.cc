#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/AnalysisLoader.hh"

namespace Rivet {

  class MC_TTBAR_HADRON : public Analysis {
  public:

    /// Minimal constructor
    MC_TTBAR_HADRON() : Analysis("MC_TTBAR_HADRON")
    {
    }


    /// @name Analysis methods
    //@{

    /// Set up projections and book histograms
    void init() {

      // A FinalState is used to select particles within |eta| < 4.2 and with pT
      // > 30 GeV, out of which the ChargedLeptons projection picks only the
      // electrons and muons, to be accessed later as "LFS".
      ChargedLeptons lfs(FinalState(-2.5, 2.5, 30*GeV));
      addProjection(lfs, "LFS");
      // A second FinalState is used to select all particles in |eta| < 4.2,
      // with no pT cut. This is used to construct jets and measure missing
      // transverse energy.
      VetoedFinalState fs(FinalState(-2.5, 2.5, 0*GeV));
      fs.addVetoOnThisFinalState(lfs);
      addProjection(FastJets(fs, FastJets::ANTIKT, 0.5), "Jets");
      addProjection(MissingMomentum(fs), "MissingET");

      // Booking of histograms
      _h_weight = bookHisto1D("weight", 200, -2, 2);
      _h_scale  = bookHisto1D("scale", logspace(50, 100.0, 2000.0));
      //
      _h_njets  = bookHisto1D("jet_mult", 11, -0.5, 10.5);
      _h_njets60  = bookHisto1D("jet_mult60", 11, -0.5, 10.5);
      _h_njets90  = bookHisto1D("jet_mult90", 11, -0.5, 10.5);
      _h_nbjets = bookHisto1D("bjet_mult", 6, -0.5, 5.5);
      _h_nleps  = bookHisto1D("lep_mult", 6, -0.5, 5.5);
      //
      _h_jet_1_pT = bookHisto1D("jet_1_pT", logspace(50, 20.0, 500.0));
      _h_jet_2_pT = bookHisto1D("jet_2_pT", logspace(50, 20.0, 400.0));
      _h_jet_3_pT = bookHisto1D("jet_3_pT", logspace(50, 20.0, 300.0));
      _h_jet_4_pT = bookHisto1D("jet_4_pT", logspace(50, 20.0, 200.0));
      _h_jet_5_pT = bookHisto1D("jet_5_pT", logspace(50, 20.0, 200.0));
      _h_jet_HT   = bookHisto1D("jet_HT", logspace(50, 100.0, 2000.0));
      //
      _h_bjet_1_pT = bookHisto1D("jetb_1_pT", logspace(50, 20.0, 400.0));
      _h_bjet_2_pT = bookHisto1D("jetb_2_pT", logspace(50, 20.0, 300.0));
      //
      _h_ljet_1_pT = bookHisto1D("jetl_1_pT", logspace(50, 20.0, 400.0));
      _h_ljet_2_pT = bookHisto1D("jetl_2_pT", logspace(50, 20.0, 300.0));
      //
      _h_Rbq       = bookHisto1D("Rbq", 50, 0.0, 5.0);
      _h_Rbq_W_cut = bookHisto1D("Rbq_W_cut", 50, 0.0, 5.0);
      //
      _h_jet_1_eta  = bookHisto1D( "jet_1_eta", 30, -3., 3.);
      _h_jet_2_eta  = bookHisto1D( "jet_2_eta", 30, -3., 3.);
      _h_jet_3_eta  = bookHisto1D( "jet_3_eta", 30, -3., 3.);
      _h_jet_4_eta  = bookHisto1D( "jet_4_eta", 30, -3., 3.);
      _h_jet_5_eta  = bookHisto1D( "jet_5_eta", 30, -3., 3.);
      _h_bjet_1_eta = bookHisto1D("jetb_1_eta", 30, -3., 3.);
      _h_bjet_2_eta = bookHisto1D("jetb_2_eta", 30, -3., 3.);
      _h_ljet_1_eta = bookHisto1D("jetl_1_eta", 30, -3., 3.);
      _h_ljet_2_eta = bookHisto1D("jetl_2_eta", 30, -3., 3.);
      //
      _h_bjet_1_mass = bookHisto1D("jetb_1_mass", 50, 0., 50.);
      _h_bjet_2_mass = bookHisto1D("jetb_2_mass", 50, 0., 50.);
      _h_ljet_1_mass = bookHisto1D("jetl_1_mass", 50, 0., 50.);
      _h_ljet_2_mass = bookHisto1D("jetl_2_mass", 50, 0., 50.);
      //
      _h_W_mass       = bookHisto1D("W_mass",        75,  30, 180);
      _h_t_mass       = bookHisto1D("t_mass",       100, 100, 300);
      _h_t_mass_W_cut = bookHisto1D("t_mass_W_cut", 100, 100, 300);
      _h_t_mass_W_scaled = bookHisto1D("t_mass_W_scaled", 100, 100, 300);
      //
      _h_W_pT         = bookHisto1D("W_pT",       50, 0, 400);
      _h_W_pT_W_cut   = bookHisto1D("W_pT_W_cut", 50, 0, 400);
      _h_t_pT_W_cut   = bookHisto1D("t_pT_W_cut", 50, 0, 400);
      //
      _h_jetb_1_jetb_2_dR   = bookHisto1D("jetb_1_jetb_2_dR", 20, 0.0, 7.0);
      _h_jetb_1_jetb_2_deta = bookHisto1D("jetb_1_jetb_2_deta", 20, 0.0, 7.0);
      _h_jetb_1_jetb_2_dphi = bookHisto1D("jetb_1_jetb_2_dphi", 20, 0.0, M_PI);
      _h_jetb_1_jetl_1_dR   = bookHisto1D("jetb_1_jetl_1_dR", 20, 0.0, 7.0);
      _h_jetb_1_jetl_1_deta = bookHisto1D("jetb_1_jetl_1_deta", 20, 0.0, 7.0);
      _h_jetb_1_jetl_1_dphi = bookHisto1D("jetb_1_jetl_1_dphi", 20, 0.0, M_PI);
      _h_jetl_1_jetl_2_dR   = bookHisto1D("jetl_1_jetl_2_dR", 20, 0.0, 7.0);
      _h_jetl_1_jetl_2_deta = bookHisto1D("jetl_1_jetl_2_deta", 20, 0.0, 7.0);
      _h_jetl_1_jetl_2_dphi = bookHisto1D("jetl_1_jetl_2_dphi", 20, 0.0, M_PI);
      _h_jetb_1_W_dR        = bookHisto1D("jetb_1_W_dR", 20, 0.0, 7.0);
      _h_jetb_1_W_deta      = bookHisto1D("jetb_1_W_deta", 20, 0.0, 7.0);
      _h_jetb_1_W_dphi      = bookHisto1D("jetb_1_W_dphi", 20, 0.0, M_PI);
      _h_jetb_1_l_dR        = bookHisto1D("jetb_1_l_dR", 20, 0.0, 7.0);
      _h_jetb_1_l_deta      = bookHisto1D("jetb_1_l_deta", 20, 0.0, 7.0);
      _h_jetb_1_l_dphi      = bookHisto1D("jetb_1_l_dphi", 20, 0.0, M_PI);
      _h_jetb_1_l_mass      = bookHisto1D("jetb_1_l_mass", 40, 0.0, 500.0);
    }


    void analyze(const Event& event) {
      const double weight = event.weight();
      _h_weight->fill(weight);
      _h_scale->fill(event.genEvent()->event_scale(), weight);

      // Use the "LFS" projection to require at least one hard charged
      // lepton. This is an experimental signature for the leptonically decaying
      // W. This helps to reduce pure QCD backgrounds.
      const ChargedLeptons& lfs = applyProjection<ChargedLeptons>(event, "LFS");
      MSG_DEBUG("Charged lepton multiplicity = " << lfs.chargedLeptons().size());
      _h_nleps->fill(lfs.chargedLeptons().size(), weight);
      foreach (const Particle& lepton, lfs.chargedLeptons()) {
        MSG_DEBUG("Lepton pT = " << lepton.pT());
      }
      if (lfs.chargedLeptons().size() != 1) {
        MSG_DEBUG("Event failed lepton multiplicity cut");
        vetoEvent;
      }

      // Use a missing ET cut to bias toward events with a hard neutrino from
      // the leptonically decaying W. This helps to reduce pure QCD backgrounds.
      //const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MissingET");
      //MSG_DEBUG("Vector ET = " << met.vectorEt().mod() << " GeV");
      //if (met.vectorEt().mod() < 30*GeV) {
      //  MSG_DEBUG("Event failed missing ET cut");
      //  vetoEvent;
      //}

      // Use the "Jets" projection to check that there are at least 4 jets of
      // any pT. Getting the jets sorted by pT ensures that the first jet is the
      // hardest, and so on. We apply no pT cut here only because we want to
      // plot all jet pTs to help optimise our jet pT cut.
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
      const Jets alljets = jetpro.jetsByPt();
      if (alljets.size() < 4) {
        MSG_DEBUG("Event failed jet multiplicity cut");
        vetoEvent;
      }

      // Update passed-cuts counter and fill all-jets histograms
      _h_jet_1_pT->fill(alljets[0].pT()/GeV, weight);
      _h_jet_2_pT->fill(alljets[1].pT()/GeV, weight);
      _h_jet_3_pT->fill(alljets[2].pT()/GeV, weight);
      _h_jet_4_pT->fill(alljets[3].pT()/GeV, weight);
      _h_jet_5_pT->fill(alljets[4].pT()/GeV, weight);
      
      _h_jet_1_eta->fill(alljets[0].eta(), weight);
      _h_jet_2_eta->fill(alljets[1].eta(), weight);
      _h_jet_3_eta->fill(alljets[2].eta(), weight);
      _h_jet_4_eta->fill(alljets[3].eta(), weight);
      _h_jet_5_eta->fill(alljets[4].eta(), weight);

      // Insist that the hardest 4 jets pass pT hardness cuts. If we don't find
      // at least 4 such jets, we abandon this event.
      //Jets jets;
      //foreach (const Jet& j, jetpro) {
      //  /// @todo Use direct kinematics access
      //  const double pt = j->momentum().pT();
      //  const double eta = j->momentum().eta();
      //  if (pt > 30*GeV && fabs(eta) < 2.4) jets.push_back(j);
      //}
      
      const Jets jets = jetpro.jetsByPt(30*GeV);
      _h_njets->fill(jets.size(), weight);
      const Jets jets60 = jetpro.jetsByPt(60*GeV);
      _h_njets60->fill(jets60.size(), weight);
      const Jets jets90 = jetpro.jetsByPt(90*GeV);
      _h_njets90->fill(jets90.size(), weight);
      
      double ht = 0.0;
      foreach (const Jet& j, jets) { ht += j.pT(); }
      _h_jet_HT->fill(ht/GeV, weight);
      if (jets.size() < 4) {
        MSG_DEBUG("Event failed jet cuts");
        vetoEvent;
      }

      // Sort the jets into b-jets and light jets. We expect one hard b-jet from
      // each top decay, so our 4 hardest jets should include two b-jets. The
      // Jet::containsBottom() method is equivalent to perfect experimental
      // b-tagging, in a generator-independent way.
      
      // Get b hadrons with pT > 5 GeV
      /// @todo This is a hack -- replace with UnstableFinalState
      vector<HepMC::GenParticle*> B_hadrons;
      vector<HepMC::GenParticle*> allParticles = particles(event.genEvent());
      for (size_t i = 0; i < allParticles.size(); i++) {
        GenParticle* p = allParticles[i];
        if (!PID::isHadron(p->pdg_id()) || !PID::hasBottom(p->pdg_id())) continue;
        if (p->momentum().perp() < 5*GeV) continue;
        B_hadrons.push_back(p);
      }
      
      // Check energy/momentum conservation
      double energy = 0.;
      double px = 0;
      double py = 0;
      double pz = 0;
      for (size_t i = 0; i < allParticles.size(); i++) {
        GenParticle* p = allParticles[i];
        if (p->status() == 1) {
          energy += p->momentum().e();
          px += p->momentum().px();
          py += p->momentum().py();
          pz += p->momentum().pz();
        }
      }
      std::cout << "Total energy/momentum: " << energy << "/" << px << "/" << py << "/" << pz << std::endl;
      
      Jets bjets, bbarjets, bmixjets, ljets;
      foreach (const Jet& jet, jets) {
        // // Don't count jets that overlap with the hard leptons
        bool isolated = true;
        foreach (const Particle& lepton, lfs.chargedLeptons()) {
          if (deltaR(jet.momentum(), lepton.momentum()) < 0.4) {
            isolated = false;
            break;
          }
        }
        if (!isolated) {
          MSG_DEBUG("Jet failed lepton isolation cut");
          break;
        }
        bool isbJet    = false;
        bool isbbarJet = false;
        //std::cout << std::endl;
        foreach(HepMC::GenParticle* b, B_hadrons) {
          if (deltaR(jet.momentum(), FourMomentum(b->momentum())) < 0.4) {
            if (b->pdg_id() < 0) isbJet    = true;
            if (b->pdg_id() > 0) isbbarJet = true;
            //std::cout << "id " << b->pdg_id() << " e " << b->momentum().e() << " status " << b->status() << std::endl;
          }
        }
        if (isbJet && !isbbarJet) {
          bjets.push_back(jet);
        } else if (isbbarJet && !isbJet) {
          bbarjets.push_back(jet);
        } else if (!isbJet && !isbbarJet) {
          ljets.push_back(jet);
        }
        else bmixjets.push_back(jet);
      }
      MSG_DEBUG("Number of b-jets = " << bjets.size()+bbarjets.size()+bmixjets.size());
      //std::cout << "b  " <<    bjets.size() << std::endl;
      //std::cout << "b~ " << bbarjets.size() << std::endl;
      _h_nbjets->fill(bjets.size()+bbarjets.size()+bmixjets.size(), weight);
      MSG_DEBUG("Number of l-jets = " << ljets.size());
      if (bjets.size() != 1 || bbarjets.size() != 1) {
        MSG_DEBUG("Event failed post-lepton-isolation b-tagging cut");
        vetoEvent;
      }
      if (ljets.size() < 2) {
        MSG_DEBUG("Event failed since not enough light jets remaining after lepton-isolation");
        vetoEvent;
      }

      // Plot the pTs of the identified jets.
      _h_bjet_1_pT->fill(bjets[0].pT(), weight);
      _h_bjet_2_pT->fill(bbarjets[0].pT(), weight);
      _h_ljet_1_pT->fill(ljets[0].pT(), weight);
      _h_ljet_2_pT->fill(ljets[1].pT(), weight);
      
      _h_Rbq->fill((bjets[0].pT()+bbarjets[0].pT())/(ljets[0].pT()+ljets[1].pT()), weight);
      
      _h_bjet_1_eta->fill(bjets[0].eta(), weight);
      _h_bjet_2_eta->fill(bbarjets[0].eta(), weight);
      _h_ljet_1_eta->fill(ljets[0].eta(), weight);
      _h_ljet_2_eta->fill(ljets[1].eta(), weight);
      
      _h_bjet_1_mass->fill(bjets[0].mass(), weight);
      _h_bjet_2_mass->fill(bbarjets[0].mass(), weight);
      _h_ljet_1_mass->fill(ljets[0].mass(), weight);
      _h_ljet_2_mass->fill(ljets[1].mass(), weight);

      // Construct the hadronically decaying W momentum 4-vector from pairs of
      // non-b-tagged jets. The pair which best matches the W mass is used. We start
      // with an always terrible 4-vector estimate which should always be "beaten" by
      // a real jet pair.
      FourMomentum W(10*sqrtS(), 0, 0, 0);
      for (size_t i = 0; i < ljets.size()-1; ++i) {
        for (size_t j = i + 1; j < ljets.size(); ++j) {
          const FourMomentum Wcand = ljets[i].momentum() + ljets[j].momentum();
          MSG_TRACE(i << "," << j << ": candidate W mass = " << Wcand.mass()/GeV
                    << " GeV, vs. incumbent candidate with " << W.mass()/GeV << " GeV");
          if (fabs(Wcand.mass() - 80.4*GeV) < fabs(W.mass() - 80.4*GeV)) {
            W = Wcand;
          }
        }
      }
      MSG_DEBUG("Candidate W mass = " << W.mass() << " GeV");

      // There are two b-jets with which this can be combined to make the
      // hadronically decaying top, one of which is correct and the other is
      // not... but we have no way to identify which is which, so we construct
      // both possible top momenta and fill the histograms with both.
      FourMomentum t1;
      if (lfs.chargedLeptons()[0].charge() < 0) t1 = W +    bjets[0].momentum();
      else                                      t1 = W + bbarjets[0].momentum();
      //const FourMomentum t2 = W + bbarjets[0].momentum();
      _h_W_mass->fill(W.mass(),  weight);
      _h_W_pT  ->fill(W.pT(),    weight);
      _h_t_mass->fill(t1.mass(), weight);
      //_h_t_mass->fill(t2.mass(), weight);
      
      FourMomentum W_scaled = W * 80.4 / W.mass();      
      FourMomentum t1_W_scaled;
      if (lfs.chargedLeptons()[0].charge() < 0) t1_W_scaled = W_scaled +    bjets[0].momentum();
      else                                      t1_W_scaled = W_scaled + bbarjets[0].momentum();

      // Placing a cut on the well-known W mass helps to reduce backgrounds
      if (inRange(W.mass()/GeV, 75.0, 85.0)) {
        MSG_DEBUG("W found with mass " << W.mass()/GeV << " GeV");
        _h_W_pT_W_cut->fill(W.pT(), weight);
        
        _h_t_mass_W_cut->fill(t1.mass(), weight);
        _h_t_mass_W_scaled->fill(t1_W_scaled.mass(), weight);
        //_h_t_mass_W_cut->fill(t2.mass(), weight);
        
        _h_t_pT_W_cut->fill(t1.pT(), weight);
        //_h_t_pT_W_cut->fill(t2.pT(), weight);
        
        _h_Rbq->fill((bjets[0].pT()+bbarjets[0].pT())/(ljets[0].pT()+ljets[1].pT()), weight);

        _h_jetb_1_jetb_2_dR->fill(deltaR(bjets[0].momentum(), bbarjets[0].momentum()),weight);
        _h_jetb_1_jetb_2_deta->fill(fabs(bjets[0].eta()-bbarjets[0].eta()),weight);
        _h_jetb_1_jetb_2_dphi->fill(deltaPhi(bjets[0].momentum(),bbarjets[0].momentum()),weight);

        _h_jetb_1_jetl_1_dR->fill(deltaR(bjets[0].momentum(), ljets[0].momentum()),weight);
        _h_jetb_1_jetl_1_deta->fill(fabs(bjets[0].eta()-ljets[0].eta()),weight);
        _h_jetb_1_jetl_1_dphi->fill(deltaPhi(bjets[0].momentum(),ljets[0].momentum()),weight);

        _h_jetl_1_jetl_2_dR->fill(deltaR(ljets[0].momentum(), ljets[1].momentum()),weight);
        _h_jetl_1_jetl_2_deta->fill(fabs(ljets[0].eta()-ljets[1].eta()),weight);
        _h_jetl_1_jetl_2_dphi->fill(deltaPhi(ljets[0].momentum(),ljets[1].momentum()),weight);

        _h_jetb_1_W_dR->fill(deltaR(bjets[0].momentum(), W),weight);
        _h_jetb_1_W_deta->fill(fabs(bjets[0].eta()-W.eta()),weight);
        _h_jetb_1_W_dphi->fill(deltaPhi(bjets[0].momentum(),W),weight);

        FourMomentum l=lfs.chargedLeptons()[0].momentum();
        _h_jetb_1_l_dR->fill(deltaR(bjets[0].momentum(), l),weight);
        _h_jetb_1_l_deta->fill(fabs(bjets[0].eta()-l.eta()),weight);
        _h_jetb_1_l_dphi->fill(deltaPhi(bjets[0].momentum(),l),weight);
        _h_jetb_1_l_mass->fill(FourMomentum(bjets[0].momentum()+l).mass(), weight);
      }

    }


    void finalize() {
      // Normalization here tends to break job merging
      /*
      normalize(_h_weight);
      normalize(_h_scale);
      normalize(_h_njets);
      normalize(_h_nbjets);
      normalize(_h_nleps);
      normalize(_h_jet_1_pT);
      normalize(_h_jet_2_pT);
      normalize(_h_jet_3_pT);
      normalize(_h_jet_4_pT);
      normalize(_h_jet_5_pT);
      normalize(_h_jet_1_eta);
      normalize(_h_jet_2_eta);
      normalize(_h_jet_3_eta);
      normalize(_h_jet_4_eta);
      normalize(_h_jet_5_eta);
      normalize(_h_jet_HT);
      normalize(_h_bjet_1_pT);
      normalize(_h_bjet_2_pT);
      normalize(_h_ljet_1_pT);
      normalize(_h_ljet_2_pT);
      normalize(_h_bjet_1_eta);
      normalize(_h_bjet_2_eta);
      normalize(_h_ljet_1_eta);
      normalize(_h_ljet_2_eta);
      normalize(_h_bjet_1_mass);
      normalize(_h_bjet_2_mass);
      normalize(_h_ljet_1_mass);
      normalize(_h_ljet_2_mass);
      normalize(_h_W_mass);
      normalize(_h_t_mass);
      normalize(_h_t_mass_W_cut);
      normalize(_h_t_mass_W_scaled);
      normalize(_h_t_pT_W_cut);
      normalize(_h_jetb_1_jetb_2_dR);
      normalize(_h_jetb_1_jetb_2_deta);
      normalize(_h_jetb_1_jetb_2_dphi);
      normalize(_h_jetb_1_jetl_1_dR);
      normalize(_h_jetb_1_jetl_1_deta);
      normalize(_h_jetb_1_jetl_1_dphi);
      normalize(_h_jetl_1_jetl_2_dR);
      normalize(_h_jetl_1_jetl_2_deta);
      normalize(_h_jetl_1_jetl_2_dphi);
      normalize(_h_jetb_1_W_dR);
      normalize(_h_jetb_1_W_deta);
      normalize(_h_jetb_1_W_dphi);
      normalize(_h_jetb_1_l_dR);
      normalize(_h_jetb_1_l_deta);
      normalize(_h_jetb_1_l_dphi);
      normalize(_h_jetb_1_l_mass);
      */
    }

    //@}


  private:

    // @name Histogram data members
    //@{

    Histo1DPtr _h_weight, _h_scale;
    Histo1DPtr _h_njets, _h_njets60, _h_njets90, _h_nbjets, _h_nleps;
    Histo1DPtr _h_jet_1_pT, _h_jet_2_pT, _h_jet_3_pT, _h_jet_4_pT, _h_jet_5_pT;
    Histo1DPtr _h_jet_1_eta, _h_jet_2_eta, _h_jet_3_eta, _h_jet_4_eta, _h_jet_5_eta;
    Histo1DPtr _h_jet_HT;
    Histo1DPtr _h_bjet_1_pT, _h_bjet_2_pT;
    Histo1DPtr _h_ljet_1_pT, _h_ljet_2_pT;
    Histo1DPtr _h_Rbq, _h_Rbq_W_cut;
    Histo1DPtr _h_bjet_1_eta, _h_bjet_2_eta;
    Histo1DPtr _h_ljet_1_eta, _h_ljet_2_eta;
    Histo1DPtr _h_bjet_1_mass, _h_bjet_2_mass;
    Histo1DPtr _h_ljet_1_mass, _h_ljet_2_mass;
    Histo1DPtr _h_W_mass, _h_W_pT, _h_W_pT_W_cut;
    Histo1DPtr _h_t_mass, _h_t_mass_W_cut, _h_t_mass_W_scaled, _h_t_pT_W_cut;
    Histo1DPtr _h_jetb_1_jetb_2_dR, _h_jetb_1_jetb_2_deta, _h_jetb_1_jetb_2_dphi;
    Histo1DPtr _h_jetb_1_jetl_1_dR, _h_jetb_1_jetl_1_deta, _h_jetb_1_jetl_1_dphi;
    Histo1DPtr _h_jetl_1_jetl_2_dR, _h_jetl_1_jetl_2_deta, _h_jetl_1_jetl_2_dphi;
    Histo1DPtr _h_jetb_1_W_dR, _h_jetb_1_W_deta, _h_jetb_1_W_dphi;
    Histo1DPtr _h_jetb_1_l_dR, _h_jetb_1_l_deta, _h_jetb_1_l_dphi,_h_jetb_1_l_mass;


    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TTBAR_HADRON);

}
