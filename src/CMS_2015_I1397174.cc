#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"


namespace Rivet {

  namespace { //< only visible in this compilation unit

    // @brief Parton level top quark finder
    // 
    // Find top quark in the parton level directly tracking particle history.
    // This does not fit with the Rivet philosophy and can be generator dependent,
    // so please use this with your own risks.
    class PartonTop : public FinalState {
      public:
        enum TTbarMode { CH_FULLHADRON = 0, CH_SEMILEPTON, CH_FULLLEPTON };
        enum DecayMode { CH_HADRON = 0, CH_MUON, CH_ELECTRON, CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON };

        /// @name Standard constructors and destructors.
        //@{

        /// The default constructor.
        PartonTop() : FinalState(-MAXDOUBLE, MAXDOUBLE, 0.0*GeV)
      {
        setName("PartonTop");
      }

        /// Clone on the heap.
        virtual const Projection* clone() const {
          return new PartonTop(*this);
        }

        //@}

        TTbarMode mode() const {
          const bool isLepton1 = _mode1%3 != 0;
          const bool isLepton2 = _mode2%3 != 0;
          if      (  isLepton1 &&  isLepton2 ) return CH_FULLLEPTON;
          else if ( !isLepton1 && !isLepton2 ) return CH_FULLHADRON;
          return CH_SEMILEPTON;
        }
        DecayMode mode1() const { return _mode1; }
        DecayMode mode2() const { return _mode2; }

        Particle t1() const { return _t1; }
        Particle t2() const { return _t2; }
        Particle b1() const { return _b1; }
        Particle b2() const { return _b2; }
        ParticleVector wDecays1() const { return _wDecays1; }
        ParticleVector wDecays2() const { return _wDecays2; }
        Particle lepton1() const { return findLepton(_wDecays1); };
        Particle lepton2() const { return findLepton(_wDecays2); };

      protected:
        // Apply the projection to the event
        void project(const Event& e) {
          _theParticles.clear();
          _wDecays1.clear();
          _wDecays2.clear();
          _mode1 = _mode2 = CH_HADRON; // Set default decay mode to full-hadronic
          _t1 = _t2 = _b1 = _b2 = Particle();

          const double ptmin = 0;
          const double etamin = -MAXDOUBLE, etamax = MAXDOUBLE;

          int nTop = 0;
          bool isTau1 = false, isTau2 = false;
          foreach (const GenParticle* p, Rivet::particles(e.genEvent())) {
            const int pdgId = p->pdg_id();
            const int absId = abs(pdgId);
            if ( absId > 20 ) continue; // We are only interested in quarks and leptons
            //if ( PID::isHadron(pdgId) ) continue; // skip hadrons
            //if ( pdgId == 22 ) continue; // skip photons
            //if ( pdgId == 91 or pdgId == 92 ) continue; // Skip cluster, strings

            if ( isZero(p->momentum().perp()) || p->momentum().perp() < ptmin ) continue;
            if ( !inRange(p->momentum().eta(), etamin, etamax) ) continue;

            // Avoid double counting by skipping if particle ID == parent ID
            std::vector<GenParticle*> pps;
            if ( absId == 6 and p->end_vertex() != 0 ) {
              pps = Rivet::particles(p->end_vertex(), HepMC::children);
            }
            else if ( absId != 6 and p->production_vertex() != 0 )
            {
              pps = Rivet::particles(p->production_vertex(), HepMC::parents);
            }
            else continue;

            bool isDuplicated = false;
            foreach (GenParticle* pp, pps) {
              if ( p != pp && p->pdg_id() == pp->pdg_id() ) {
                isDuplicated = true;
                break;
              }
            }
            if ( isDuplicated ) continue;

            // Build Rivet::Particle
            Particle rp(*p);
            // Skip particles from hadronization (and keep tau decay)
            if ( rp.fromDecay() and !rp.hasAncestor(15) and !rp.hasAncestor(-15) ) continue;

            if      ( pdgId ==  6 ) { nTop++; _t1 = rp; }
            else if ( pdgId == -6 ) { nTop++; _t2 = rp; }
            else if ( pdgId ==  5 and rp.pT() > _b1.pT() ) _b1 = rp;
            else if ( pdgId == -5 and rp.pT() > _b2.pT() ) _b2 = rp;
            else if ( absId <= 16 && rp.hasAncestor( 24) ) {
              if ( pdgId == -15 ) isTau1 = true;
              else if ( pdgId == -11 ) _mode1 = CH_ELECTRON;
              else if ( pdgId == -13 ) _mode1 = CH_MUON;
              _wDecays1.push_back(rp);
            }
            else if ( absId <= 16 && rp.hasAncestor(-24) ) {
              if ( pdgId == 15 ) isTau2 = true;
              else if ( pdgId == 11 ) _mode2 = CH_ELECTRON;
              else if ( pdgId == 13 ) _mode2 = CH_MUON;
              _wDecays2.push_back(rp);
            }
          }
          if ( isTau1 ) _mode1 = static_cast<DecayMode>(_mode1+3);
          if ( isTau2 ) _mode2 = static_cast<DecayMode>(_mode2+3);

        }

        Particle findLepton(const ParticleVector& v) const {
          Particle pp = Particle();
          foreach (const Particle& p, v) {
            const int aid = std::abs(p.pdgId());
            if ( (aid == 11 or aid == 13) and  p.pT() > pp.pT() ) pp = p;
          }
          return pp;
        }

      private:
        DecayMode _mode1, _mode2;
        Particle _t1, _t2;
        Particle _b1, _b2;
        ParticleVector _wDecays1, _wDecays2;
    };

  }

  class CMS_2015_I1397174 : public Analysis {
    public:
      /// Minimal constructor
      CMS_2015_I1397174() : Analysis("CMS_2015_I1397174") {}


      /// @name Analysis methods
      //@{

      /// Set up projections and book histograms
      void init() {
        // Parton level top quarks
        addProjection(PartonTop(), "partonTop");

        //FastJets fj(CMSGenParticle(), FastJets::ANTIKT, 0.5);
        VetoedFinalState vfs;
        vfs.addDecayProductsVeto( 24);
        vfs.addDecayProductsVeto(-24);
        FastJets fj(vfs, FastJets::ANTIKT, 0.5);
        fj.useInvisibles();
        addProjection(fj, "Jets");

        _vis_unit_weights = 0.;
        _full_unit_weights = 0.;
        // Book histograms
        _hVis_nJet30_abs       = bookHisto1D( 1, 1, 1);
        _hVis_nJet30           = bookHisto1D( 2, 1, 1);
        _hVis_nJet60_abs       = bookHisto1D( 3, 1, 1);
        _hVis_nJet60           = bookHisto1D( 4, 1, 1);
        _hVis_nJet100_abs      = bookHisto1D( 5, 1, 1);
        _hVis_nJet100          = bookHisto1D( 6, 1, 1);

        _hVis_addJet1Pt_abs    = bookHisto1D( 7, 1, 1);
        _hVis_addJet1Pt        = bookHisto1D( 8, 1, 1);
        _hVis_addJet1Eta_abs   = bookHisto1D( 9, 1, 1);
        _hVis_addJet1Eta       = bookHisto1D(10, 1, 1);
        _hVis_addJet2Pt_abs    = bookHisto1D(11, 1, 1);
        _hVis_addJet2Pt        = bookHisto1D(12, 1, 1);
        _hVis_addJet2Eta_abs   = bookHisto1D(13, 1, 1);
        _hVis_addJet2Eta       = bookHisto1D(14, 1, 1);
        _hVis_addJJMass_abs    = bookHisto1D(15, 1, 1);
        _hVis_addJJMass        = bookHisto1D(16, 1, 1);
        _hVis_addJJDR_abs      = bookHisto1D(17, 1, 1);
        _hVis_addJJDR          = bookHisto1D(18, 1, 1);
        _hVis_addJJHT_abs      = bookHisto1D(19, 1, 1);
        _hVis_addJJHT          = bookHisto1D(20, 1, 1);

        _hFull_addJet1Pt_abs   = bookHisto1D(21, 1, 1);
        _hFull_addJet1Pt       = bookHisto1D(22, 1, 1);
        _hFull_addJet1Eta_abs  = bookHisto1D(23, 1, 1);
        _hFull_addJet1Eta      = bookHisto1D(24, 1, 1);
        _hFull_addJet2Pt_abs   = bookHisto1D(25, 1, 1);
        _hFull_addJet2Pt       = bookHisto1D(26, 1, 1);
        _hFull_addJet2Eta_abs  = bookHisto1D(27, 1, 1);
        _hFull_addJet2Eta      = bookHisto1D(28, 1, 1);
        _hFull_addJJMass_abs   = bookHisto1D(29, 1, 1);
        _hFull_addJJMass       = bookHisto1D(30, 1, 1);
        _hFull_addJJDR_abs     = bookHisto1D(31, 1, 1);
        _hFull_addJJDR         = bookHisto1D(32, 1, 1);
        _hFull_addJJHT_abs     = bookHisto1D(33, 1, 1);
        _hFull_addJJHT         = bookHisto1D(34, 1, 1);

        _hVis_addBJet1Pt_abs   = bookHisto1D(35, 1, 1);
        _hVis_addBJet1Pt       = bookHisto1D(36, 1, 1);
        _hVis_addBJet1Eta_abs  = bookHisto1D(37, 1, 1);
        _hVis_addBJet1Eta      = bookHisto1D(38, 1, 1);
        _hVis_addBJet2Pt_abs   = bookHisto1D(39, 1, 1);
        _hVis_addBJet2Pt       = bookHisto1D(40, 1, 1);
        _hVis_addBJet2Eta_abs  = bookHisto1D(41, 1, 1);
        _hVis_addBJet2Eta      = bookHisto1D(42, 1, 1);
        _hVis_addBBMass_abs    = bookHisto1D(43, 1, 1);
        _hVis_addBBMass        = bookHisto1D(44, 1, 1);
        _hVis_addBBDR_abs      = bookHisto1D(45, 1, 1);
        _hVis_addBBDR          = bookHisto1D(46, 1, 1);

        _hFull_addBJet1Pt_abs  = bookHisto1D(47, 1, 1);
        _hFull_addBJet1Pt      = bookHisto1D(48, 1, 1);
        _hFull_addBJet1Eta_abs = bookHisto1D(49, 1, 1);
        _hFull_addBJet1Eta     = bookHisto1D(50, 1, 1);
        _hFull_addBJet2Pt_abs  = bookHisto1D(51, 1, 1);
        _hFull_addBJet2Pt      = bookHisto1D(52, 1, 1);
        _hFull_addBJet2Eta_abs = bookHisto1D(53, 1, 1);
        _hFull_addBJet2Eta     = bookHisto1D(54, 1, 1);
        _hFull_addBBMass_abs   = bookHisto1D(55, 1, 1);
        _hFull_addBBMass       = bookHisto1D(56, 1, 1);
        _hFull_addBBDR_abs     = bookHisto1D(57, 1, 1);
        _hFull_addBBDR         = bookHisto1D(58, 1, 1);

        _gap_weights = 0.;
        _h_gap_addJet1Pt       = bookHisto1D(59, 1, 1);
        _h_gap_addJet1Pt_eta0  = bookHisto1D(60, 1, 1);
        _h_gap_addJet1Pt_eta1  = bookHisto1D(61, 1, 1);
        _h_gap_addJet1Pt_eta2  = bookHisto1D(62, 1, 1);
        _h_gap_addJet2Pt       = bookHisto1D(63, 1, 1);
        _h_gap_addJet2Pt_eta0  = bookHisto1D(64, 1, 1);
        _h_gap_addJet2Pt_eta1  = bookHisto1D(65, 1, 1);
        _h_gap_addJet2Pt_eta2  = bookHisto1D(66, 1, 1);
        _h_gap_addJetHT        = bookHisto1D(67, 1, 1);
        _h_gap_addJetHT_eta0   = bookHisto1D(68, 1, 1);
        _h_gap_addJetHT_eta1   = bookHisto1D(69, 1, 1);
        _h_gap_addJetHT_eta2   = bookHisto1D(70, 1, 1);
      }


      void analyze(const Event& event) {
        const double weight = event.weight();
        if (weight != 0.) _vis_unit_weights += weight/std::abs(weight);
        // The objects used in the PAPER 12-041 is defined as follows (see p.16 for details):
        // 
        //   * Leptons    : from the W boson decays after FSR
        //   * Jets       : anti-kT R=0.5 to all stable particles 
        //                               exclude W->enu, munu, taunu
        //   * B jet      : B-Ghost matched
        //   * B from top : B hadron from top->b decay
        //
        // Visible phase space definition:
        //
        //   * Leptons         : pT > 20, |eta| < 2.4
        //   * B jets from top : pT > 30, |eta| < 2.4
        //     Additional jets : pT > 20, |eta| < 2.4
        //   * 
        // Full phase space definition:
        //
        //   * Correction to dilepton BR from W boson BR
        //   * No cut on top decay products
        //   * Additional jets : pT > 20, |eta| < 2.4
        //

        const PartonTop& partonTop = applyProjection<PartonTop>(event, "partonTop");
        // Do the analysis only for the ttbar full leptonic channel, removing tau decays
        if ( partonTop.mode() != PartonTop::CH_FULLLEPTON ) vetoEvent;
        if ( partonTop.mode1() >= PartonTop::CH_TAU_HADRON or
            partonTop.mode2() >= PartonTop::CH_TAU_HADRON ) vetoEvent;

        // Apply acceptance cut on muon or electrons
        // Get the lepton of the ttbar decay
        const Particle lep1 = getLast(partonTop.lepton1());
        const Particle lep2 = getLast(partonTop.lepton2());
        if ( lep1.pT() <= 1e-9 or lep2.pT() <= 1e-9 ) vetoEvent; // Just for sanity check

        const Jets& jets = applyProjection<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.4);

        if (weight != 0.) _full_unit_weights += weight/std::abs(weight);

        int nJet30 = 0, nJet60 = 0, nJet100 = 0;
        Jets topBJets, addJets, addBJets, addJets_eta0, addJets_eta1, addJets_eta2;
        foreach ( const Jet& jet, jets ) {
          //if ( deltaR(lep1.momentum(), jet.momentum()) < 0.4 ) continue;
          //if ( deltaR(lep2.momentum(), jet.momentum()) < 0.4 ) continue;

          const double pt = jet.pT();
          if ( pt >  30*GeV ) ++nJet30;
          if ( pt >  60*GeV ) ++nJet60;
          if ( pt > 100*GeV ) ++nJet100;

          bool isBtagged = false, isBFromTop = false;
          foreach(const Particle& p, jet.bTags()) {
            isBtagged = true;
            if ( isFromTop(p) ) {
              isBFromTop = true;
              break;
            }
          }
          if ( isBFromTop ) {
            if ( jet.pT() > 30*GeV ) topBJets.push_back(jet);
          } else {
            addJets.push_back(jet);
            if      ( std::abs(jet.eta()) < 0.8 ) addJets_eta0.push_back(jet);
            else if ( std::abs(jet.eta()) < 1.5 ) addJets_eta1.push_back(jet);
            else if ( std::abs(jet.eta()) < 2.4 ) addJets_eta2.push_back(jet);
            if ( isBtagged ) addBJets.push_back(jet);
          }
        }
        const bool isVisiblePS = lep1.pT() > 20*GeV and std::abs(lep1.eta()) < 2.4 and
          lep2.pT() > 20*GeV and std::abs(lep2.eta()) < 2.4 and
          topBJets.size() >= 2;
        if ( isVisiblePS ) {
          fillWithOF(_hVis_nJet30_abs,  nJet30, weight);
          fillWithOF(_hVis_nJet30,      nJet30, weight);
          fillWithOF(_hVis_nJet60_abs,  nJet60, weight);
          fillWithOF(_hVis_nJet60,      nJet60, weight);
          fillWithOF(_hVis_nJet100_abs, nJet100, weight);
          fillWithOF(_hVis_nJet100,     nJet100, weight);

          _gap_weights += weight;
          fillGapFractions(addJets, _h_gap_addJet1Pt, _h_gap_addJet2Pt, _h_gap_addJetHT, weight);
          fillGapFractions(addJets_eta0, _h_gap_addJet1Pt_eta0, _h_gap_addJet2Pt_eta0, _h_gap_addJetHT_eta0, weight);
          fillGapFractions(addJets_eta1, _h_gap_addJet1Pt_eta1, _h_gap_addJet2Pt_eta1, _h_gap_addJetHT_eta1, weight);
          fillGapFractions(addJets_eta2, _h_gap_addJet1Pt_eta2, _h_gap_addJet2Pt_eta2, _h_gap_addJetHT_eta2, weight);
        }

        // Plots with at least two additional jets
        do {
          if ( addJets.size() < 1 ) break;
          const double ht = std::accumulate(addJets.begin(), addJets.end(),
              0., [](double x, const Jet& jj){return x+jj.pT();});
          _hFull_addJJHT_abs->fill(ht, weight);
          _hFull_addJJHT    ->fill(ht, weight);
          if ( isVisiblePS ) {
            _hVis_addJJHT_abs->fill(ht, weight);
            _hVis_addJJHT    ->fill(ht, weight);
          }

          const double j1pt = addJets[0].pT(), j1aeta = std::abs(addJets[0].eta());
          _hFull_addJet1Pt_abs ->fill(j1pt  , weight);
          _hFull_addJet1Pt     ->fill(j1pt  , weight);
          _hFull_addJet1Eta_abs->fill(j1aeta, weight);
          _hFull_addJet1Eta    ->fill(j1aeta, weight);
          if ( isVisiblePS ) {
            _hVis_addJet1Pt_abs ->fill(j1pt  , weight);
            _hVis_addJet1Pt     ->fill(j1pt  , weight);
            _hVis_addJet1Eta_abs->fill(j1aeta, weight);
            _hVis_addJet1Eta    ->fill(j1aeta, weight);
          }

          if ( addJets.size() < 2 ) break;

          const double j2pt = addJets[1].pT(), j2aeta = std::abs(addJets[1].eta());
          const double jjmass = (addJets[0].momentum()+addJets[1].momentum()).mass();
          const double jjdR = deltaR(addJets[0], addJets[1]);

          _hFull_addJet2Pt_abs ->fill(j2pt  , weight);
          _hFull_addJet2Pt     ->fill(j2pt  , weight);
          _hFull_addJet2Eta_abs->fill(j2aeta, weight);
          _hFull_addJet2Eta    ->fill(j2aeta, weight);
          if ( isVisiblePS ) {
            _hVis_addJet2Pt_abs ->fill(j2pt  , weight);
            _hVis_addJet2Pt     ->fill(j2pt  , weight);
            _hVis_addJet2Eta_abs->fill(j2aeta, weight);
            _hVis_addJet2Eta    ->fill(j2aeta, weight);
          }

          _hFull_addJJMass_abs->fill(jjmass, weight);
          _hFull_addJJMass    ->fill(jjmass, weight);
          _hFull_addJJDR_abs  ->fill(jjdR, weight);
          _hFull_addJJDR      ->fill(jjdR, weight);
          if ( isVisiblePS ) {
            _hVis_addJJMass_abs->fill(jjmass, weight);
            _hVis_addJJMass    ->fill(jjmass, weight);
            _hVis_addJJDR_abs  ->fill(jjdR, weight);
            _hVis_addJJDR      ->fill(jjdR, weight);
          }
        } while ( false );

        // Same set of plots if there are additional b-jets
        do {
          if ( addBJets.size() < 1 ) break;
          const double b1pt = addBJets[0].pT(), b1aeta = std::abs(addBJets[0].eta());
          _hFull_addBJet1Pt_abs ->fill(b1pt  , weight);
          _hFull_addBJet1Pt     ->fill(b1pt  , weight);
          _hFull_addBJet1Eta_abs->fill(b1aeta, weight);
          _hFull_addBJet1Eta    ->fill(b1aeta, weight);
          if ( isVisiblePS ) {
            _hVis_addBJet1Pt_abs ->fill(b1pt  , weight);
            _hVis_addBJet1Pt     ->fill(b1pt  , weight);
            _hVis_addBJet1Eta_abs->fill(b1aeta, weight);
            _hVis_addBJet1Eta    ->fill(b1aeta, weight);
          }

          if ( addBJets.size() < 2 ) break;
          const double b2pt = addBJets[1].pT(), b2aeta = std::abs(addBJets[1].eta());
          const double bbmass = (addBJets[0].momentum()+addBJets[1].momentum()).mass();
          const double bbdR = deltaR(addBJets[0], addBJets[1]);

          _hFull_addBJet2Pt_abs ->fill(b2pt  , weight);
          _hFull_addBJet2Pt     ->fill(b2pt  , weight);
          _hFull_addBJet2Eta_abs->fill(b2aeta, weight);
          _hFull_addBJet2Eta    ->fill(b2aeta, weight);
          if ( isVisiblePS ) {
            _hVis_addBJet2Pt_abs ->fill(b2pt  , weight);
            _hVis_addBJet2Pt     ->fill(b2pt  , weight);
            _hVis_addBJet2Eta_abs->fill(b2aeta, weight);
            _hVis_addBJet2Eta    ->fill(b2aeta, weight);
          }

          _hFull_addBBMass_abs->fill(bbmass, weight);
          _hFull_addBBMass    ->fill(bbmass, weight);
          _hFull_addBBDR_abs  ->fill(bbdR, weight);
          _hFull_addBBDR      ->fill(bbdR, weight);

          // Fill plots in visible phase space
          if ( isVisiblePS ) {
            _hVis_addBBMass_abs->fill(bbmass, weight);
            _hVis_addBBMass    ->fill(bbmass, weight);
            _hVis_addBBDR_abs  ->fill(bbdR, weight);
            _hVis_addBBDR      ->fill(bbdR, weight);
          }
        } while ( false );

      }

      void finalize() {
        double ttbarXS = 0.;
        if (!isnan(crossSectionPerEvent())) {
          std::cout << "Using generator cross section: " << crossSection() << " pb" << std::endl;
          ttbarXS = crossSection();
        }
        else {
          std::cout << "No valid cross section given, using NNLO value: 252.89 pb" << std::endl;
          ttbarXS = 252.89; // NNLO (arXiv:1303.6254; sqrt(s)=8 TeV, m_t=172.5 GeV)
          // see also https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
        }

        std::cout << "Sum vis unit weights: " << _vis_unit_weights
          << ", sum full unit weights: " << _full_unit_weights << std::endl;

        scale(_hVis_nJet30_abs,       ttbarXS / _vis_unit_weights);
        scale(_hVis_nJet60_abs,       ttbarXS / _vis_unit_weights);
        scale(_hVis_nJet100_abs,      ttbarXS / _vis_unit_weights);
        scale(_hVis_addJet1Pt_abs,    ttbarXS / _vis_unit_weights);
        scale(_hVis_addJet1Eta_abs,   ttbarXS / _vis_unit_weights);
        scale(_hVis_addJet2Pt_abs,    ttbarXS / _vis_unit_weights);
        scale(_hVis_addJet2Eta_abs,   ttbarXS / _vis_unit_weights);
        scale(_hVis_addJJMass_abs,    ttbarXS / _vis_unit_weights);
        scale(_hVis_addJJDR_abs,      ttbarXS / _vis_unit_weights);
        scale(_hVis_addJJHT_abs,      ttbarXS / _vis_unit_weights);
        scale(_hFull_addJet1Pt_abs,   ttbarXS / _full_unit_weights);
        scale(_hFull_addJet1Eta_abs,  ttbarXS / _full_unit_weights);
        scale(_hFull_addJet2Pt_abs,   ttbarXS / _full_unit_weights);
        scale(_hFull_addJet2Eta_abs,  ttbarXS / _full_unit_weights);
        scale(_hFull_addJJMass_abs,   ttbarXS / _full_unit_weights);
        scale(_hFull_addJJDR_abs,     ttbarXS / _full_unit_weights);
        scale(_hFull_addJJHT_abs,     ttbarXS / _full_unit_weights);
        scale(_hVis_addBJet1Pt_abs,   ttbarXS / _vis_unit_weights);
        scale(_hVis_addBJet1Eta_abs,  ttbarXS / _vis_unit_weights);
        scale(_hVis_addBJet2Pt_abs,   ttbarXS / _vis_unit_weights);
        scale(_hVis_addBJet2Eta_abs,  ttbarXS / _vis_unit_weights);
        scale(_hVis_addBBMass_abs,    ttbarXS / _vis_unit_weights);
        scale(_hVis_addBBDR_abs,      ttbarXS / _vis_unit_weights);
        scale(_hFull_addBJet1Pt_abs,  ttbarXS / _full_unit_weights);
        scale(_hFull_addBJet1Eta_abs, ttbarXS / _full_unit_weights);
        scale(_hFull_addBJet2Pt_abs,  ttbarXS / _full_unit_weights);
        scale(_hFull_addBJet2Eta_abs, ttbarXS / _full_unit_weights);
        scale(_hFull_addBBMass_abs,   ttbarXS / _full_unit_weights);
        scale(_hFull_addBBDR_abs,     ttbarXS / _full_unit_weights);

        normalize(_hVis_nJet30);
        normalize(_hVis_nJet60);
        normalize(_hVis_nJet100);
        normalize(_hVis_addJet1Pt);
        normalize(_hVis_addJet1Eta );
        normalize(_hVis_addJet2Pt);
        normalize(_hVis_addJet2Eta);
        normalize(_hVis_addJJMass);
        normalize(_hVis_addJJDR);
        normalize(_hVis_addJJHT);
        normalize(_hFull_addJet1Pt);
        normalize(_hFull_addJet1Eta);
        normalize(_hFull_addJet2Pt);
        normalize(_hFull_addJet2Eta);
        normalize(_hFull_addJJMass);
        normalize(_hFull_addJJDR);
        normalize(_hFull_addJJHT);
        normalize(_hVis_addBJet1Pt);
        normalize(_hVis_addBJet1Eta);
        normalize(_hVis_addBJet2Pt);
        normalize(_hVis_addBJet2Eta);
        normalize(_hVis_addBBMass);
        normalize(_hVis_addBBDR);
        normalize(_hFull_addBJet1Pt);
        normalize(_hFull_addBJet1Eta);
        normalize(_hFull_addBJet2Pt);
        normalize(_hFull_addBJet2Eta);
        normalize(_hFull_addBBMass);
        normalize(_hFull_addBBDR);

        const double s = 1./_gap_weights;
        scale(_h_gap_addJet1Pt     , s);
        scale(_h_gap_addJet1Pt_eta0, s);
        scale(_h_gap_addJet1Pt_eta1, s);
        scale(_h_gap_addJet1Pt_eta2, s);
        scale(_h_gap_addJet2Pt     , s);
        scale(_h_gap_addJet2Pt_eta0, s);
        scale(_h_gap_addJet2Pt_eta1, s);
        scale(_h_gap_addJet2Pt_eta2, s);
        scale(_h_gap_addJetHT      , s);
        scale(_h_gap_addJetHT_eta0 , s);
        scale(_h_gap_addJetHT_eta1 , s);
        scale(_h_gap_addJetHT_eta2 , s);
      }

      //@}


    private:
      inline void fillWithOF(Histo1DPtr h, const double x, const double w) const {
        h->fill(std::min(x, h->xMax()-1e-9), w);
      }

      const Particle getLast(const Particle& p) const {
        if ( !p.genParticle() or !p.genParticle()->end_vertex() ) return p;

        foreach (const Particle& dp, Rivet::particles(p.genParticle()->end_vertex(), HepMC::children)) {
          if ( dp.pdgId() == p.pdgId() ) return getLast(dp);;
        }

        return p;
      }

      bool isFromTop(const Particle& p) const {
        if ( !p.genParticle() or !p.genParticle()->production_vertex() ) return false;

        foreach ( const Particle& ap, Rivet::particles(p.genParticle()->production_vertex(), HepMC::ancestors)) {
          if ( abs(ap.pdgId()) == 6 ) return true;
        }

        return false;
      }

      void fillGapFractions(Jets addJets, Histo1DPtr h_gap_addJet1Pt, Histo1DPtr h_gap_addJet2Pt, Histo1DPtr h_gap_addJetHT, double weight) const {
        double j1pt   = 0.;      
        if  ( addJets.size() > 0 ) {
          j1pt   = addJets[0].pT();
        }
        for (unsigned int i = 0; i < h_gap_addJet1Pt->numBins(); ++i) {
          double binCenter = h_gap_addJet1Pt->bin(i).xMid();
          double binWidth  = h_gap_addJet1Pt->bin(i).xWidth();
          if (j1pt < binCenter) {
            h_gap_addJet1Pt->fillBin(i, binWidth*weight);
          }
        }

        double j2pt   = 0.;
        if  ( addJets.size() > 1 ) {
          j2pt   = addJets[1].pT();
        }
        for (unsigned int i = 0; i < h_gap_addJet2Pt->numBins(); ++i) {
          double binCenter = h_gap_addJet2Pt->bin(i).xMid();
          double binWidth  = h_gap_addJet2Pt->bin(i).xWidth();
          if (j2pt < binCenter) {
            h_gap_addJet2Pt->fillBin(i, binWidth*weight);
          }
        }

        double ht = std::accumulate(addJets.begin(), addJets.end(),
            0., [](double x, const Jet& jj){return x+jj.pT();});
        for (unsigned int i = 0; i < h_gap_addJetHT->numBins(); ++i) {
          double binCenter = h_gap_addJetHT->bin(i).xMid();
          double binWidth  = h_gap_addJetHT->bin(i).xWidth();
          if (ht < binCenter) {
            h_gap_addJetHT->fillBin(i, binWidth*weight);
          }
        }
      }

      // @name Histogram data members
      //@{

      double _vis_unit_weights;
      double _full_unit_weights;

      Histo1DPtr _hVis_nJet30_abs, _hVis_nJet60_abs, _hVis_nJet100_abs;
      Histo1DPtr _hVis_addJet1Pt_abs, _hVis_addJet1Eta_abs, _hVis_addJet2Pt_abs, _hVis_addJet2Eta_abs;
      Histo1DPtr _hVis_addJJMass_abs, _hVis_addJJDR_abs, _hVis_addJJHT_abs;
      Histo1DPtr _hFull_addJet1Pt_abs, _hFull_addJet1Eta_abs, _hFull_addJet2Pt_abs, _hFull_addJet2Eta_abs;
      Histo1DPtr _hFull_addJJMass_abs, _hFull_addJJDR_abs, _hFull_addJJHT_abs;
      Histo1DPtr _hVis_addBJet1Pt_abs, _hVis_addBJet1Eta_abs, _hVis_addBJet2Pt_abs, _hVis_addBJet2Eta_abs;
      Histo1DPtr _hVis_addBBMass_abs, _hVis_addBBDR_abs;
      Histo1DPtr _hFull_addBJet1Pt_abs, _hFull_addBJet1Eta_abs, _hFull_addBJet2Pt_abs, _hFull_addBJet2Eta_abs;
      Histo1DPtr _hFull_addBBMass_abs, _hFull_addBBDR_abs;

      Histo1DPtr _hVis_nJet30, _hVis_nJet60, _hVis_nJet100;
      Histo1DPtr _hVis_addJet1Pt, _hVis_addJet1Eta, _hVis_addJet2Pt, _hVis_addJet2Eta;
      Histo1DPtr _hVis_addJJMass, _hVis_addJJDR, _hVis_addJJHT;
      Histo1DPtr _hFull_addJet1Pt, _hFull_addJet1Eta, _hFull_addJet2Pt, _hFull_addJet2Eta;
      Histo1DPtr _hFull_addJJMass, _hFull_addJJDR, _hFull_addJJHT;
      Histo1DPtr _hVis_addBJet1Pt, _hVis_addBJet1Eta, _hVis_addBJet2Pt, _hVis_addBJet2Eta;
      Histo1DPtr _hVis_addBBMass, _hVis_addBBDR;
      Histo1DPtr _hFull_addBJet1Pt, _hFull_addBJet1Eta, _hFull_addBJet2Pt, _hFull_addBJet2Eta;
      Histo1DPtr _hFull_addBBMass, _hFull_addBBDR;

      double _gap_weights;
      Histo1DPtr _h_gap_addJet1Pt, _h_gap_addJet1Pt_eta0, _h_gap_addJet1Pt_eta1, _h_gap_addJet1Pt_eta2;
      Histo1DPtr _h_gap_addJet2Pt, _h_gap_addJet2Pt_eta0, _h_gap_addJet2Pt_eta1, _h_gap_addJet2Pt_eta2;
      Histo1DPtr _h_gap_addJetHT, _h_gap_addJetHT_eta0, _h_gap_addJetHT_eta1, _h_gap_addJetHT_eta2;

      //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2015_I1397174);

}

