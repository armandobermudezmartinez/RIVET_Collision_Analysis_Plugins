// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

#include "HepMC/GenParticle.h"

#include <vector>

namespace Rivet
{


 class CMS_2016_I1491950 : public Analysis
 {
  private:
   Particles _muelphs;
   Particles _neutrinos;
   Particles _jetparticles;
   Jets _dressLeps;
   vector<const Jet*> _bJets;
   const Jet*_bl;
   const Jet* _bh;
   const Jet* _wja;
   const Jet* _wjb;
   FourMomentum _tl;
   FourMomentum _th;
   FourMomentum _wl;
   FourMomentum _wh;
   FourMomentum _nusum;

   Histo1DPtr _hist_thadpt;
   Histo1DPtr _hist_thady;
   Histo1DPtr _hist_tleppt;
   Histo1DPtr _hist_tlepy;
   Histo1DPtr _hist_ttpt;
   Histo1DPtr _hist_tty;
   Histo1DPtr _hist_ttm;
   Histo1DPtr _hist_njet;
   Histo1DPtr _hist_njets_thadpt_1;
   Histo1DPtr _hist_njets_thadpt_2;
   Histo1DPtr _hist_njets_thadpt_3;
   Histo1DPtr _hist_njets_thadpt_4;
   Histo1DPtr _hist_njets_ttpt_1;
   Histo1DPtr _hist_njets_ttpt_2;
   Histo1DPtr _hist_njets_ttpt_3;
   Histo1DPtr _hist_njets_ttpt_4;
   Histo1DPtr _hist_thady_thadpt_1;
   Histo1DPtr _hist_thady_thadpt_2;
   Histo1DPtr _hist_thady_thadpt_3;
   Histo1DPtr _hist_thady_thadpt_4;
   Histo1DPtr _hist_ttm_tty_1;
   Histo1DPtr _hist_ttm_tty_2;
   Histo1DPtr _hist_ttm_tty_3;
   Histo1DPtr _hist_ttm_tty_4;
   Histo1DPtr _hist_ttpt_ttm_1;
   Histo1DPtr _hist_ttpt_ttm_2;
   Histo1DPtr _hist_ttpt_ttm_3;
   Histo1DPtr _hist_ttpt_ttm_4;

   Histo1DPtr _histnorm_thadpt;
   Histo1DPtr _histnorm_thady;
   Histo1DPtr _histnorm_tleppt;
   Histo1DPtr _histnorm_tlepy;
   Histo1DPtr _histnorm_ttpt;
   Histo1DPtr _histnorm_tty;
   Histo1DPtr _histnorm_ttm;
   Histo1DPtr _histnorm_njet;
   Histo1DPtr _histnorm_njets_thadpt_1;
   Histo1DPtr _histnorm_njets_thadpt_2;
   Histo1DPtr _histnorm_njets_thadpt_3;
   Histo1DPtr _histnorm_njets_thadpt_4;
   Histo1DPtr _histnorm_njets_ttpt_1;
   Histo1DPtr _histnorm_njets_ttpt_2;
   Histo1DPtr _histnorm_njets_ttpt_3;
   Histo1DPtr _histnorm_njets_ttpt_4;
   Histo1DPtr _histnorm_thady_thadpt_1;
   Histo1DPtr _histnorm_thady_thadpt_2;
   Histo1DPtr _histnorm_thady_thadpt_3;
   Histo1DPtr _histnorm_thady_thadpt_4;
   Histo1DPtr _histnorm_ttm_tty_1;
   Histo1DPtr _histnorm_ttm_tty_2;
   Histo1DPtr _histnorm_ttm_tty_3;
   Histo1DPtr _histnorm_ttm_tty_4;
   Histo1DPtr _histnorm_ttpt_ttm_1;
   Histo1DPtr _histnorm_ttpt_ttm_2;
   Histo1DPtr _histnorm_ttpt_ttm_3;
   Histo1DPtr _histnorm_ttpt_ttm_4;

  public:

   double _jetetamax;
   double _jetptmin;
   double _lepetamax;
   double _lepptmin;


   /// Constructor
   CMS_2016_I1491950()
    : Analysis("CMS_2016_I1491950"), _jetetamax(2.5), _jetptmin(25.), _lepetamax(2.5), _lepptmin(30.)
   {    }

   /// Book histograms and initialise projections before the run
   void init()
   {
    const FinalState fs(Cuts::pT > 0. && Cuts::abseta < 6.);
    addProjection(fs, "FS");

    //book hists
    _hist_thadpt = bookHisto1D("d01-x02-y01");
    _hist_thady = bookHisto1D("d03-x02-y01");
    _hist_tleppt = bookHisto1D("d05-x02-y01");
    _hist_tlepy = bookHisto1D("d07-x02-y01");
    _hist_ttpt = bookHisto1D("d09-x02-y01");
    _hist_tty = bookHisto1D("d13-x02-y01");
    _hist_ttm = bookHisto1D("d11-x02-y01");
    _hist_njet = bookHisto1D("d15-x02-y01");
    _hist_njets_thadpt_1 = bookHisto1D("d17-x02-y01");
    _hist_njets_thadpt_2 = bookHisto1D("d18-x02-y01");
    _hist_njets_thadpt_3 = bookHisto1D("d19-x02-y01");
    _hist_njets_thadpt_4 = bookHisto1D("d20-x02-y01");
    _hist_njets_ttpt_1 = bookHisto1D("d22-x02-y01");
    _hist_njets_ttpt_2 = bookHisto1D("d23-x02-y01");
    _hist_njets_ttpt_3 = bookHisto1D("d24-x02-y01");
    _hist_njets_ttpt_4 = bookHisto1D("d25-x02-y01");
    _hist_thady_thadpt_1 = bookHisto1D("d27-x02-y01");
    _hist_thady_thadpt_2 = bookHisto1D("d28-x02-y01");
    _hist_thady_thadpt_3 = bookHisto1D("d29-x02-y01");
    _hist_thady_thadpt_4 = bookHisto1D("d30-x02-y01");
    _hist_ttm_tty_1 = bookHisto1D("d32-x02-y01");
    _hist_ttm_tty_2 = bookHisto1D("d33-x02-y01");
    _hist_ttm_tty_3 = bookHisto1D("d34-x02-y01");
    _hist_ttm_tty_4 = bookHisto1D("d35-x02-y01");
    _hist_ttpt_ttm_1 = bookHisto1D("d37-x02-y01");
    _hist_ttpt_ttm_2 = bookHisto1D("d38-x02-y01");
    _hist_ttpt_ttm_3 = bookHisto1D("d39-x02-y01");
    _hist_ttpt_ttm_4 = bookHisto1D("d40-x02-y01");

    _histnorm_thadpt = bookHisto1D("d42-x02-y01");
    _histnorm_thady = bookHisto1D("d44-x02-y01");
    _histnorm_tleppt = bookHisto1D("d46-x02-y01");
    _histnorm_tlepy = bookHisto1D("d48-x02-y01");
    _histnorm_ttpt = bookHisto1D("d50-x02-y01");
    _histnorm_tty = bookHisto1D("d54-x02-y01");
    _histnorm_ttm = bookHisto1D("d52-x02-y01");
    _histnorm_njet = bookHisto1D("d56-x02-y01");
    _histnorm_njets_thadpt_1 = bookHisto1D("d58-x02-y01");
    _histnorm_njets_thadpt_2 = bookHisto1D("d59-x02-y01");
    _histnorm_njets_thadpt_3 = bookHisto1D("d60-x02-y01");
    _histnorm_njets_thadpt_4 = bookHisto1D("d61-x02-y01");
    _histnorm_njets_ttpt_1 = bookHisto1D("d63-x02-y01");
    _histnorm_njets_ttpt_2 = bookHisto1D("d64-x02-y01");
    _histnorm_njets_ttpt_3 = bookHisto1D("d65-x02-y01");
    _histnorm_njets_ttpt_4 = bookHisto1D("d66-x02-y01");
    _histnorm_thady_thadpt_1 = bookHisto1D("d68-x02-y01");
    _histnorm_thady_thadpt_2 = bookHisto1D("d69-x02-y01");
    _histnorm_thady_thadpt_3 = bookHisto1D("d70-x02-y01");
    _histnorm_thady_thadpt_4 = bookHisto1D("d71-x02-y01");
    _histnorm_ttm_tty_1 = bookHisto1D("d73-x02-y01");
    _histnorm_ttm_tty_2 = bookHisto1D("d74-x02-y01");
    _histnorm_ttm_tty_3 = bookHisto1D("d75-x02-y01");
    _histnorm_ttm_tty_4 = bookHisto1D("d76-x02-y01");
    _histnorm_ttpt_ttm_1 = bookHisto1D("d78-x02-y01");
    _histnorm_ttpt_ttm_2 = bookHisto1D("d79-x02-y01");
    _histnorm_ttpt_ttm_3 = bookHisto1D("d80-x02-y01");
    _histnorm_ttpt_ttm_4 = bookHisto1D("d81-x02-y01");

   }


   /// Perform the per-event analysis
   void analyze(const Event& event)
   {
    const double weight = event.weight();
    _muelphs.clear();
    _neutrinos.clear();
    _dressLeps.clear();
    _bJets.clear();
    _jetparticles.clear();

    const ParticleVector& FSpars = applyProjection<FinalState>(event, "FS").particles();
    for(ParticleVector::const_iterator itpar = FSpars.begin() ; itpar != FSpars.end() ; ++itpar)
    {
     const Particle& par = *itpar;
     if(isFromHadron(par.genParticle())) continue;
     int pdgid = abs(par.pdgId());
     if(pdgid == 12 || pdgid == 14 || pdgid == 16)
     {
      _neutrinos.push_back(par);
     }
     else if(pdgid == 11 || pdgid == 13 || pdgid == 22)
     {
      _muelphs.push_back(par);
     }
    }

    FastJets _fjDressLep(FinalState(), FastJets::ANTIKT, 0.1);
    _fjDressLep.calc(_muelphs);
    const Jets dressLeps = _fjDressLep.jets(Cuts::abseta < _lepetamax && Cuts::pT > _lepptmin*GeV);
    for(Jets::const_iterator itdl = dressLeps.begin() ; itdl != dressLeps.end() ; ++itdl)
    {
     Particles jetconst = itdl->particles();
     for(Particles::const_iterator itpar = jetconst.begin() ; itpar != jetconst.end(); ++itpar)
     {
      if((abs(itpar->pdgId()) == 11 || abs(itpar->pdgId()) == 13) && itpar->pt()/itdl->pt() > 0.5)
      {
       _dressLeps.push_back(*itdl);
      } 
     }
    }

    if(_dressLeps.size() != 1) return;

    const vector<const HepMC::GenParticle*>& allgenpars = Rivet::particles(event.genEvent());
    for(vector<const HepMC::GenParticle*>::const_iterator itgenpar = allgenpars.begin() ; itgenpar != allgenpars.end() ; ++itgenpar)
    {
     const HepMC::GenParticle* genpar = *itgenpar;
     const int pdgid = genpar->pdg_id();
     if(isFinalBHadron(genpar))
     {
      _jetparticles.push_back(Particle(555555555, 1E-20*FourMomentum(genpar->momentum())));
      continue;
     }
     if(genpar->status() != 1) continue;

     //remove selected neutrinos
     bool isneutrino = false; 
     for(Particles::const_iterator itnu = _neutrinos.begin() ; itnu != _neutrinos.end(); ++itnu)
     {
      if(itnu->genParticle() == genpar)
      {
       isneutrino = true;
       break;  
      } 
     }
     if(isneutrino) continue;

     //remove dressed lepton constituents
     bool islepton = false; 
     for(Jets::const_iterator itdl = _dressLeps.begin() ; itdl != _dressLeps.end() ; ++itdl)
     {
      Particles jetconst = itdl->particles();
      for(Particles::const_iterator itpar = jetconst.begin() ; itpar != jetconst.end(); ++itpar)
      {
       if(itpar->genParticle() == genpar)
       {
        islepton = true;   
        break;
       } 
      }
     }
     if(islepton) continue;

     _jetparticles.push_back(Particle(pdgid, FourMomentum(genpar->momentum())));
    }

    FastJets _fjJets(FinalState(), FastJets::ANTIKT, 0.4);
    _fjJets.calc(_jetparticles);
    const Jets allJets = _fjJets.jets(Cuts::abseta < _jetetamax && Cuts::pT > _jetptmin*GeV);
    for(Jets::const_iterator itjet = allJets.begin() ; itjet != allJets.end() ; ++itjet)
    {
     const Particles jetconst = itjet->particles();
     for(Particles::const_iterator itpar = jetconst.begin() ; itpar != jetconst.end(); ++itpar)
     {
      if(itpar->pdgId() == 555555555)
      {
       _bJets.push_back(&(*itjet));
       break;
      } 
     }
    }

    if(_bJets.size() < 2 || allJets.size() < 4) return;


    _nusum = FourMomentum(0., 0., 0., 0.);
    for(Particles::const_iterator itnu = _neutrinos.begin() ; itnu != _neutrinos.end(); ++itnu)
    {
     _nusum += itnu->momentum();
    }
    _wl = _nusum + _dressLeps[0].momentum();

    //construct top quark proxies
    double Kmin = numeric_limits<double>::max();
    for(Jets::const_iterator itaj = allJets.begin() ; itaj != allJets.end() ; ++itaj)
    {
     for(Jets::const_iterator itbj = allJets.begin() ; itbj != itaj ; ++itbj)
     {
      FourMomentum wh(itaj->momentum() + itbj->momentum());
      for(vector<const Jet*>::const_iterator ithbj = _bJets.begin() ; ithbj != _bJets.end() ; ++ithbj)
      {
       FourMomentum th(wh + (*ithbj)->momentum());
       if(&(*itaj) == *ithbj || &(*itbj) == *ithbj) continue;
       for(vector<const Jet*>::const_iterator itlbj = _bJets.begin() ; itlbj != _bJets.end() ; ++itlbj)
       {
        if(&(*itaj) == *itlbj || &(*itbj) == *itlbj || ithbj == itlbj) continue;
        FourMomentum tl(_wl + (*itlbj)->momentum());

        double K = pow(wh.mass() - 80.4, 2) + pow(th.mass() - 172.5, 2) + pow(tl.mass() - 172.5, 2);
        if(K < Kmin)
        {
         Kmin = K;
         _bl = *itlbj;
         _bh = *ithbj;
         _wja = &(*itaj); 
         _wjb = &(*itbj);
         _tl = tl;
         _th = th;
         _wh = wh;
        }
       }
      }
     }
    }

    _hist_thadpt->fill(_th.pt(), weight); 
    _hist_thady->fill(abs(_th.rapidity()) , weight);
    _hist_tleppt->fill(_tl.pt() , weight);
    _hist_tlepy->fill(abs(_tl.rapidity()) , weight);
    _histnorm_thadpt->fill(_th.pt(), weight); 
    _histnorm_thady->fill(abs(_th.rapidity()) , weight);
    _histnorm_tleppt->fill(_tl.pt() , weight);
    _histnorm_tlepy->fill(abs(_tl.rapidity()) , weight);
    FourMomentum tt(_tl+_th);
    _hist_ttpt->fill(tt.pt() , weight);
    _hist_tty->fill(abs(tt.rapidity()) , weight);
    _hist_ttm->fill(tt.mass() , weight);
    _hist_njet->fill(min(allJets.size()-4., 4.), weight);
    _histnorm_ttpt->fill(tt.pt() , weight);
    _histnorm_tty->fill(abs(tt.rapidity()) , weight);
    _histnorm_ttm->fill(tt.mass() , weight);
    _histnorm_njet->fill(min(allJets.size()-4., 4.), weight);
    if(allJets.size() == 4)
    {
     _hist_njets_thadpt_1->fill(_th.pt(), weight);
     _hist_njets_ttpt_1->fill(tt.pt(), weight);
     _histnorm_njets_thadpt_1->fill(_th.pt(), weight);
     _histnorm_njets_ttpt_1->fill(tt.pt(), weight);
    }
    else if(allJets.size() == 5)
    {
     _hist_njets_thadpt_2->fill(_th.pt(), weight);
     _hist_njets_ttpt_2->fill(tt.pt(), weight);
     _histnorm_njets_thadpt_2->fill(_th.pt(), weight);
     _histnorm_njets_ttpt_2->fill(tt.pt(), weight);
    }
    else if(allJets.size() == 6)
    {
     _hist_njets_thadpt_3->fill(_th.pt(), weight);
     _hist_njets_ttpt_3->fill(tt.pt(), weight);
     _histnorm_njets_thadpt_3->fill(_th.pt(), weight);
     _histnorm_njets_ttpt_3->fill(tt.pt(), weight);
    }
    else //>= 4 jets
    {
     _hist_njets_thadpt_4->fill(_th.pt(), weight);
     _hist_njets_ttpt_4->fill(tt.pt(), weight);
     _histnorm_njets_thadpt_4->fill(_th.pt(), weight);
     _histnorm_njets_ttpt_4->fill(tt.pt(), weight);
    }

    if(abs(_th.rapidity()) < 0.5)
    {
     _hist_thady_thadpt_1->fill(_th.pt(), weight);
     _histnorm_thady_thadpt_1->fill(_th.pt(), weight);
    }
    else if(abs(_th.rapidity()) < 1.0)
    {
     _hist_thady_thadpt_2->fill(_th.pt(), weight);
     _histnorm_thady_thadpt_2->fill(_th.pt(), weight);
    }
    else if(abs(_th.rapidity()) < 1.5)
    {
     _hist_thady_thadpt_3->fill(_th.pt(), weight);
     _histnorm_thady_thadpt_3->fill(_th.pt(), weight);
    }
    else if(abs(_th.rapidity()) < 2.5)
    {
     _hist_thady_thadpt_4->fill(_th.pt(), weight);
     _histnorm_thady_thadpt_4->fill(_th.pt(), weight);
    }

    if(tt.mass() >= 300. && tt.mass() < 450.)
    {
     _hist_ttm_tty_1->fill(abs(tt.rapidity()), weight);
     _histnorm_ttm_tty_1->fill(abs(tt.rapidity()), weight);
    }
    else if(tt.mass() >= 450. && tt.mass() < 625.)
    {
     _hist_ttm_tty_2->fill(abs(tt.rapidity()), weight);
     _histnorm_ttm_tty_2->fill(abs(tt.rapidity()), weight);
    }
    else if(tt.mass() >= 625. && tt.mass() < 850.)
    {
     _hist_ttm_tty_3->fill(abs(tt.rapidity()), weight);
     _histnorm_ttm_tty_3->fill(abs(tt.rapidity()), weight);
    }
    else if(tt.mass() >= 850. && tt.mass() < 2000.)
    {
     _hist_ttm_tty_4->fill(abs(tt.rapidity()), weight);
     _histnorm_ttm_tty_4->fill(abs(tt.rapidity()), weight);
    }

    if(tt.pt() < 35.)
    {
     _hist_ttpt_ttm_1->fill(tt.mass(), weight);
     _histnorm_ttpt_ttm_1->fill(tt.mass(), weight);
    }
    else if(tt.pt() < 80.)
    {
     _hist_ttpt_ttm_2->fill(tt.mass(), weight);
     _histnorm_ttpt_ttm_2->fill(tt.mass(), weight);
    }
    else if(tt.pt() < 140.)
    {
     _hist_ttpt_ttm_3->fill(tt.mass(), weight);
     _histnorm_ttpt_ttm_3->fill(tt.mass(), weight);
    }
    else if(tt.pt() < 500.)
    {
     _hist_ttpt_ttm_4->fill(tt.mass(), weight);
     _histnorm_ttpt_ttm_4->fill(tt.mass(), weight);
    }

   }


   /// Normalise histograms etc., after the run
   void finalize()
   {
    scale(_hist_thadpt, crossSection()/sumOfWeights());
    scale(_hist_thady, crossSection()/sumOfWeights());
    scale(_hist_tleppt, crossSection()/sumOfWeights());
    scale(_hist_tlepy, crossSection()/sumOfWeights());
    scale(_hist_ttpt, crossSection()/sumOfWeights());
    scale(_hist_tty, crossSection()/sumOfWeights());
    scale(_hist_ttm, crossSection()/sumOfWeights());
    scale(_hist_njet, crossSection()/sumOfWeights());
    scale(_hist_njets_thadpt_1, crossSection()/sumOfWeights());
    scale(_hist_njets_thadpt_2, crossSection()/sumOfWeights());
    scale(_hist_njets_thadpt_3, crossSection()/sumOfWeights());
    scale(_hist_njets_thadpt_4, crossSection()/sumOfWeights());
    scale(_hist_njets_ttpt_1, crossSection()/sumOfWeights());
    scale(_hist_njets_ttpt_2, crossSection()/sumOfWeights());
    scale(_hist_njets_ttpt_3, crossSection()/sumOfWeights());
    scale(_hist_njets_ttpt_4, crossSection()/sumOfWeights());
    scale(_hist_thady_thadpt_1, crossSection()/sumOfWeights()/0.5);
    scale(_hist_thady_thadpt_2, crossSection()/sumOfWeights()/0.5);
    scale(_hist_thady_thadpt_3, crossSection()/sumOfWeights()/0.5);
    scale(_hist_thady_thadpt_4, crossSection()/sumOfWeights()/1.0);
    scale(_hist_ttm_tty_1, crossSection()/sumOfWeights()/150.);
    scale(_hist_ttm_tty_2, crossSection()/sumOfWeights()/175.);
    scale(_hist_ttm_tty_3, crossSection()/sumOfWeights()/225.);
    scale(_hist_ttm_tty_4, crossSection()/sumOfWeights()/1150.);
    scale(_hist_ttpt_ttm_1, crossSection()/sumOfWeights()/35.);
    scale(_hist_ttpt_ttm_2, crossSection()/sumOfWeights()/45.);
    scale(_hist_ttpt_ttm_3, crossSection()/sumOfWeights()/60.);
    scale(_hist_ttpt_ttm_4, crossSection()/sumOfWeights()/360.);

    scale(_histnorm_thadpt, 1./_histnorm_thadpt->sumW(false));
    scale(_histnorm_thady, 1./_histnorm_thady->sumW(false));
    scale(_histnorm_tleppt, 1./_histnorm_tleppt->sumW(false));
    scale(_histnorm_tlepy, 1./_histnorm_tlepy->sumW(false));
    scale(_histnorm_ttpt, 1./_histnorm_ttpt->sumW(false));
    scale(_histnorm_tty, 1./_histnorm_tty->sumW(false));
    scale(_histnorm_ttm, 1./_histnorm_ttm->sumW(false));
    scale(_histnorm_njet, 1./_histnorm_njet->sumW(false));
    double sum_njets_thadpt = _histnorm_njets_thadpt_1->sumW(false) + _histnorm_njets_thadpt_2->sumW(false) + _histnorm_njets_thadpt_3->sumW(false) + _histnorm_njets_thadpt_4->sumW(false);
    scale(_histnorm_njets_thadpt_1, 1./sum_njets_thadpt);
    scale(_histnorm_njets_thadpt_2, 1./sum_njets_thadpt);
    scale(_histnorm_njets_thadpt_3, 1./sum_njets_thadpt);
    scale(_histnorm_njets_thadpt_4, 1./sum_njets_thadpt);
    double sum_njets_ttpt = _histnorm_njets_ttpt_1->sumW(false) + _histnorm_njets_ttpt_2->sumW(false) + _histnorm_njets_ttpt_3->sumW(false) + _histnorm_njets_ttpt_4->sumW(false);
    scale(_histnorm_njets_ttpt_1, 1./sum_njets_ttpt);
    scale(_histnorm_njets_ttpt_2, 1./sum_njets_ttpt);
    scale(_histnorm_njets_ttpt_3, 1./sum_njets_ttpt);
    scale(_histnorm_njets_ttpt_4, 1./sum_njets_ttpt);
    double sum_thady_thadpt = _histnorm_thady_thadpt_1->sumW(false) + _histnorm_thady_thadpt_2->sumW(false) + _histnorm_thady_thadpt_3->sumW(false) + _histnorm_thady_thadpt_4->sumW(false);
    scale(_histnorm_thady_thadpt_1, 1./sum_thady_thadpt/0.5);
    scale(_histnorm_thady_thadpt_2, 1./sum_thady_thadpt/0.5);
    scale(_histnorm_thady_thadpt_3, 1./sum_thady_thadpt/0.5);
    scale(_histnorm_thady_thadpt_4, 1./sum_thady_thadpt/1.0);
    double sum_ttm_tty = _histnorm_ttm_tty_1->sumW(false) + _histnorm_ttm_tty_2->sumW(false) + _histnorm_ttm_tty_3->sumW(false) + _histnorm_ttm_tty_4->sumW(false);
    scale(_histnorm_ttm_tty_1, 1./sum_ttm_tty/150.);
    scale(_histnorm_ttm_tty_2, 1./sum_ttm_tty/175.);
    scale(_histnorm_ttm_tty_3, 1./sum_ttm_tty/225.);
    scale(_histnorm_ttm_tty_4, 1./sum_ttm_tty/1150.);
    double sum_ttpt_ttm = _histnorm_ttpt_ttm_1->sumW(false) + _histnorm_ttpt_ttm_2->sumW(false) + _histnorm_ttpt_ttm_3->sumW(false) + _histnorm_ttpt_ttm_4->sumW(false);
    scale(_histnorm_ttpt_ttm_1, 1./sum_ttpt_ttm/35.);
    scale(_histnorm_ttpt_ttm_2, 1./sum_ttpt_ttm/45.);
    scale(_histnorm_ttpt_ttm_3, 1./sum_ttpt_ttm/60.);
    scale(_histnorm_ttpt_ttm_4, 1./sum_ttpt_ttm/360.);

   }

   bool isFinalBHadron(const HepMC::GenParticle* genp)
   {
    if(isBHadron(genp->pdg_id()) == false) return false;

    const GenVertex* prodVtx = genp->production_vertex();
    if (prodVtx == nullptr) return false;
    const vector<const GenParticle*> ancestors(particles(prodVtx, HepMC::ancestors));
    for(vector<const GenParticle*>::const_iterator itan = ancestors.begin() ; itan != ancestors.end() ; ++itan)
    {
     if(isBHadron((*itan)->pdg_id()) == true) return false;
    }
    return true;
   }

   bool isBHadron(int pdgid)
   {
    int abspdgid = abs(pdgid);
    if(abspdgid <= 100) return false;

    int nq3 = (abspdgid / 10) % 10;
    int nq2 = (abspdgid / 100) % 10;
    int nq1 = (abspdgid / 1000) % 10;

    if ( nq3 == 0 ) return false;
    if ( nq1 == 5 or nq2 == 5 ) return true;

    return false;
   }

   bool isFromHadron(const HepMC::GenParticle* genp)
   {
    const GenVertex* prodVtx = genp->production_vertex();
    const vector<const GenParticle*> ancestors(particles(prodVtx, HepMC::ancestors));
    for(vector<const GenParticle*>::const_iterator itan = ancestors.begin() ; itan != ancestors.end() ; ++itan)
    {
     if((*itan)->status() == 2 && abs((*itan)->pdg_id()) > 100) return true;
    }
    return false;
   }


 };



 // The hook for the plugin system
 DECLARE_RIVET_PLUGIN(CMS_2016_I1491950);


}

