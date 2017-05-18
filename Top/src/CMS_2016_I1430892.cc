// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PartonicTops.hh"

namespace Rivet {


  /// Dilepton channel ttbar charge asymmetry analysis
  class CMS_2016_I1430892 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2016_I1430892);


    /// Book histograms and initialise projections
    void init() {

      // Complete final state
      FinalState fs(-MAXDOUBLE, MAXDOUBLE, 0*GeV);

      // Projection for dressed electrons and muons
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

      //this histogram is independent of the parton-level information
      _h_dabsetadressedleptons = bookHisto1D("d01-x01-y01", _bins_dabseta);

      //the remaining histos use parton-level information

      _h_dabseta = bookHisto1D("d05-x01-y01", _bins_dabseta);
      _h_dabsrapidity = bookHisto1D("d02-x01-y01", _bins_dabsrapidity);

      //2D histos
      _h_dabsrapidity_var[0] = bookHisto2D("d11-x01-y01", _bins_dabsrapidity, _bins_tt_mass);
      _h_dabseta_var[0] = bookHisto2D("d17-x01-y01", _bins_dabseta, _bins_tt_mass);

      _h_dabsrapidity_var[1] = bookHisto2D("d23-x01-y01", _bins_dabsrapidity, _bins_tt_pT);
      _h_dabseta_var[1] = bookHisto2D("d29-x01-y01", _bins_dabseta, _bins_tt_pT);

      _h_dabsrapidity_var[2] = bookHisto2D("d35-x01-y01", _bins_dabsrapidity, _bins_tt_absrapidity);
      _h_dabseta_var[2] = bookHisto2D("d41-x01-y01", _bins_dabseta, _bins_tt_absrapidity);

      //profile histos for asymmetries
      _h_dabsrapidity_profile[0] = bookProfile1D("d08-x01-y01", _bins_tt_mass);
      _h_dabseta_profile[0] = bookProfile1D("d14-x01-y01", _bins_tt_mass);

      _h_dabsrapidity_profile[1] = bookProfile1D("d20-x01-y01", _bins_tt_pT);
      _h_dabseta_profile[1] = bookProfile1D("d26-x01-y01", _bins_tt_pT);

      _h_dabsrapidity_profile[2] = bookProfile1D("d32-x01-y01", _bins_tt_absrapidity);
      _h_dabseta_profile[2] = bookProfile1D("d38-x01-y01", _bins_tt_absrapidity);
      
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      // use particle-level leptons for the first 2 histos
      const DressedLeptons& dressed_electrons = applyProjection<DressedLeptons>(event, "DressedElectrons");
      const DressedLeptons& dressed_muons = applyProjection<DressedLeptons>(event, "DressedMuons");

      const std::vector<DressedLepton> dressedels = dressed_electrons.dressedLeptons();
      const std::vector<DressedLepton> dressedmus = dressed_muons.dressedLeptons();

      int ndressedel = dressedels.size();
      int ndressedmu = dressedmus.size();

      // For the particle-level histos, require an odd number of electrons and an odd number of muons, to select ttbar->emu channel. This means we can easily identify additional dilepton pairs from the shower (which will be same-flavour with invariant mass m_ll~0 GeV ), which distort the distributions relative to those used in the analysis where the leptons are prompt leptons from top decay only.
      if ( ndressedel % 2 == 1 && ndressedmu % 2 == 1  ) {

        int electrontouse = 0;
        int muontouse = 0;

        if ( ndressedel > 1 ) electrontouse = returnPrimaryLepton(dressedels);
        if ( ndressedmu > 1 ) muontouse = returnPrimaryLepton(dressedmus);
  
        if ( electrontouse != -1 && muontouse != -1 ) {
          //fill particle-level histos for opposite-charge leptons only
          if ( sameSign(dressedels[electrontouse],dressedmus[muontouse]) ) {
            MSG_WARNING("error, e and mu have same charge, skipping event");
          }
          else {
            //Get the four-momenta of the positively- and negatively-charged leptons
            FourMomentum lepPlus = dressedels[electrontouse].charge() > 0 ? dressedels[electrontouse].momentum() : dressedmus[muontouse].momentum();
            FourMomentum lepMinus = dressedels[electrontouse].charge() > 0 ? dressedmus[muontouse].momentum() : dressedels[electrontouse].momentum();

            //now calculate the variable
            double dabseta_temp = lepPlus.abseta() - lepMinus.abseta();

            fillWithUFOF( _h_dabsetadressedleptons, dabseta_temp, weight );
          }
        }

      }

      // The remaining variables use parton-level information.

      // Get the leptonically decaying tops
      const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();
      Particles chargedleptons;

      unsigned int ntrueleptonictops = 0;
      bool oppositesign = false;
  
      if ( leptonicpartontops.size() == 2 ) {
        for (unsigned int k = 0; k < leptonicpartontops.size(); ++k) {
    
          //get the lepton
          const Particle lepTop = leptonicpartontops[k];
          const auto isPromptChargedLepton = [](const Particle& p){return (isChargedLepton(p) && isPrompt(p, false, false));};
          Particles lepton_candidates = lepTop.allDescendants(firstParticleWith(isPromptChargedLepton), false); 
          if ( lepton_candidates.size() < 1 ) MSG_WARNING("error, PartonicTops::E_MU top quark had no daughter lepton candidate, skipping event.");
          bool istrueleptonictop = false;
          
          // In some cases there is no lepton from the W decay but only leptons from the decay of a radiated gamma. 
          // These hadronic PartonicTops are currently being mistakenly selected by PartonicTops::E_MU (as of April 2017), and need to be rejected.
          // PartonicTops::E_MU is being fixed in Rivet, and when it is the veto below should do nothing.
          for (unsigned int i = 0; i < lepton_candidates.size(); ++i) {
            Particle lepton_candidate = lepton_candidates[i];
            if ( hasParentWith(lepton_candidate, Cuts::abspid == 22) ) {
              MSG_DEBUG("found gamma parent, top: "<<k+1<<" of "<<leptonicpartontops.size()<<" , lepton: "<<i+1<<" of "<<lepton_candidates.size());
              continue;
            }
            if ( !istrueleptonictop && sameSign(lepTop,lepton_candidate) ) {
              chargedleptons.push_back(lepton_candidate);
              istrueleptonictop = true;
            }
            else MSG_WARNING("error, found extra prompt charged lepton from top decay (and without gamma parent), ignoring it.");
          }
          if ( istrueleptonictop ) ++ntrueleptonictops;
        }
      }

      if ( ntrueleptonictops == 2 ) {
        oppositesign = !( sameSign(chargedleptons[0],chargedleptons[1]) );
        if ( !oppositesign ) MSG_WARNING("error, same charge tops, skipping event.");
      }

      if ( ntrueleptonictops == 2 && oppositesign ) {
 
        //Get the four-momenta of the positively- and negatively-charged leptons
        FourMomentum lepPlus = chargedleptons[0].charge() > 0 ? chargedleptons[0].momentum() : chargedleptons[1].momentum();
        FourMomentum lepMinus = chargedleptons[0].charge() > 0 ? chargedleptons[1].momentum() : chargedleptons[0].momentum();
 
        double dabseta_temp = lepPlus.abseta() - lepMinus.abseta();

        //Get the four-momenta of the positively- and negatively-charged tops
        FourMomentum topPlus_p4 = leptonicpartontops[0].pdgId() > 0 ? leptonicpartontops[0].momentum() : leptonicpartontops[1].momentum();
        FourMomentum topMinus_p4 = leptonicpartontops[0].pdgId() > 0 ? leptonicpartontops[1].momentum() : leptonicpartontops[0].momentum();

        FourMomentum ttbar_p4 = topPlus_p4 + topMinus_p4;

        double tt_mass_temp = ttbar_p4.mass();
        double tt_absrapidity_temp = ttbar_p4.absrapidity();
        double tt_pT_temp = ttbar_p4.pT();

        double dabsrapidity_temp = topPlus_p4.absrapidity() - topMinus_p4.absrapidity();

        //fill parton-level histos
        fillWithUFOF( _h_dabseta, dabseta_temp, weight );
        fillWithUFOF( _h_dabsrapidity, dabsrapidity_temp, weight );

        //now fill the same variables in each of their 3 bins of ttbar invariant mass, pT, and absolute rapidity
        for (int i_var = 0; i_var < 3; ++i_var) {
          double var;
          std::vector<double> bins_var;

          if ( i_var == 0 ) {
            var = tt_mass_temp;
            bins_var = _bins_tt_mass;
          }
          else if ( i_var == 1 ) {
            var = tt_pT_temp;
            bins_var = _bins_tt_pT;
          }
          else {
            var = tt_absrapidity_temp;
            bins_var = _bins_tt_absrapidity;
          }

          fillWithUFOF( _h_dabsrapidity_var[i_var], dabsrapidity_temp, var, weight );
          fillWithUFOF( _h_dabseta_var[i_var], dabseta_temp, var, weight );

          fillWithUFOF( _h_dabsrapidity_profile[i_var], dabsrapidity_temp, var, weight, (_h_dabsrapidity->xMax() + _h_dabsrapidity->xMin())/2. );
          fillWithUFOF( _h_dabseta_profile[i_var], dabseta_temp, var, weight, (_h_dabseta->xMax() + _h_dabseta->xMin())/2. );
        }

      }

    }


    /// Normalise histograms to unit area
    void finalize() {

      normalize(_h_dabsetadressedleptons);

      normalize(_h_dabseta);
      normalize(_h_dabsrapidity);

      for (int i_var = 0; i_var < 3; ++i_var) {
        normalize(_h_dabsrapidity_var[i_var]);
        normalize(_h_dabseta_var[i_var]);
      }

    }


  private:
    Histo1DPtr _h_dabsetadressedleptons, _h_dabseta, _h_dabsrapidity;
    Histo2DPtr _h_dabseta_var[3], _h_dabsrapidity_var[3];
    Profile1DPtr _h_dabseta_profile[3], _h_dabsrapidity_profile[3];

    const std::vector<double> _bins_tt_mass = {300., 430., 530., 1200.};
    const std::vector<double> _bins_tt_pT = {0., 41., 92., 300.};
    const std::vector<double> _bins_tt_absrapidity = {0., 0.34, 0.75, 1.5};
    const std::vector<double> _bins_dabseta = { -2., -68./60., -48./60., -32./60., -20./60., -8./60., 0., 8./60., 20./60., 32./60., 48./60., 68./60., 2.};
    const std::vector<double> _bins_dabsrapidity = {-2., -44./60., -20./60., 0., 20./60., 44./60., 2.};

    struct ilepsmll {
      double mll;
      int ilep1;
      int ilep2;
    };

    typedef vector< ilepsmll > Vofilepsmll;

    void fillWithUFOF(Histo1DPtr h, double x, double w) {
      h->fill(std::max(std::min(x, h->xMax()-1e-9),h->xMin()+1e-9), w);
    }

    void fillWithUFOF(Histo2DPtr h, double x, double y, double w) {
      h->fill(std::max(std::min(x, h->xMax()-1e-9),h->xMin()+1e-9), std::max(std::min(y, h->yMax()-1e-9),h->yMin()+1e-9), w);
    }

    void fillWithUFOF(Profile1DPtr h, double x, double y, double w, double c) {
      h->fill(std::max(std::min(y, h->xMax()-1e-9),h->xMin()+1e-9), float(x > c) - float(x < c), w);
    }

    int returnPrimaryLepton(const std::vector<DressedLepton> dressedleptons) {
      //when there are additional lepton pair(s) from the shower, remove the combination of same-flavour opposite-charge lepton pair(s) that has the smallest maximum invariant mass m_ll

      int ndressedlep = dressedleptons.size();

      if ( ndressedlep < 3 || ndressedlep % 2 == 0 ) {
        MSG_WARNING("error, returnPrimaryLepton only works with an odd number of 3 or more leptons, returning -1.");
        return -1;
      }

      int leptontouse = 0; //variable that will be updated with the ordinal number of the primary lepton

      int nuniquesets = 1; //counter of unique sets of lepton pairs. For 2M+1 dressed leptons, there will be nuniquesets = M!(M+1)! unique sets of pairs. Probability of M>3 negligible (i.e. >3 extra lepton pairs from radiation), so create arrays of size 3!*4! = 144.
      vector<int> leptonsused[144]; //keep track of lepton numbers used forming each unique set of pairs
      Vofilepsmll vimlls[144]; //vector of information about each pair for each unique set of lepton pairs
      //clear the vector for the first unique set of lepton pairs. The others will be cleared later if they are used.
      vimlls[0].clear();
      leptonsused[0].clear();

      if ( ndressedlep > 7 ) {
        MSG_WARNING("error, found "<<ndressedlep<<" leptons. returnPrimaryLepton is only configured to use up to 7 leptons. Returning -1.");
        return -1;
      }

      for (int i = 0; i < ndressedlep; ++i) {
        for (int j = i+1; j < ndressedlep; ++j) {

          if ( !sameSign(dressedleptons[i],dressedleptons[j]) ) {

            double mll = (dressedleptons[i].momentum() + dressedleptons[j].momentum()).mass();
            ilepsmll imll = { mll, dressedleptons[i].charge() > 0 ? i : j, dressedleptons[i].charge() > 0 ? j : i };

            //fill unique sets of lepton pairs
            bool leptonsfilled = false;
            for (int uniqueset = 0; uniqueset < nuniquesets; ++uniqueset) {
              //fill leptons into this set if neither have already been filled
              if ( !(std::find(leptonsused[uniqueset].begin(), leptonsused[uniqueset].end(), i) != leptonsused[uniqueset].end()) && !(std::find(leptonsused[uniqueset].begin(), leptonsused[uniqueset].end(), j) != leptonsused[uniqueset].end()) ) {
                leptonsused[uniqueset].push_back(i);
                leptonsused[uniqueset].push_back(j);
                vimlls[uniqueset].push_back(imll);
                leptonsfilled = true;
              }
            }
            //if the leptons can't be filled in an existing set create a new one
            if ( !leptonsfilled ) {
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
      if ( vimlls[0].size() == 0 ) {
        MSG_WARNING("error, returnPrimaryLepton found all the same-flavour leptons have the same charge. Returning -1.");
        return -1;
      }

      MSG_DEBUG("nuniquesets count check: "<<nuniquesets<<" ndressedlep: "<<ndressedlep);


      for (int uniqueset = 0; uniqueset < nuniquesets; ++uniqueset) {
        sort(vimlls[uniqueset].begin(), vimlls[uniqueset].end(), [](ilepsmll ilepsmll1, ilepsmll ilepsmll2) { return ilepsmll1.mll > ilepsmll2.mll; } );
      }

      //Find the set that has the lowest mass for the last pair to be removed.
      int lowestmassuniqueset = -1;
      double lowestmass = 99999;
      for (int uniqueset = 0; uniqueset < nuniquesets; ++uniqueset) {
        if ( vimlls[uniqueset][0].mll < lowestmass ) {
          lowestmass = vimlls[uniqueset][0].mll;
          lowestmassuniqueset = uniqueset;
        }
      }


      MSG_DEBUG("lowest mass unique set is "<<lowestmassuniqueset<<" :");
      for (unsigned int i = 0; i < vimlls[lowestmassuniqueset].size(); ++i) {
        MSG_DEBUG(vimlls[lowestmassuniqueset][i].mll<<" "<<vimlls[lowestmassuniqueset][i].ilep1<<" "<<vimlls[lowestmassuniqueset][i].ilep2);
      }

      while ( (std::find(leptonsused[lowestmassuniqueset].begin(), leptonsused[lowestmassuniqueset].end(), leptontouse) != leptonsused[lowestmassuniqueset].end()) ) ++leptontouse;
      MSG_DEBUG("Using lepton "<<leptontouse);

      if ( leptontouse >= ndressedlep ) {
        MSG_WARNING("error, returnPrimaryLepton couldn't find an unpaired lepton, returning -1.");
        return -1;
      }

      return leptontouse;

    }






  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1430892);


}
