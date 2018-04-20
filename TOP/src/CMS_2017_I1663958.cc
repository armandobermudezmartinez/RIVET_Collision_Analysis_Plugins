#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Rivet/Analysis.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"

#include "Rivet/Event.hh"
#include "Rivet/Tools/RivetHepMC.hh"
#include "HepMC/GenParticle.h"

#include <vector>
#include <algorithm>

using namespace std;
using namespace Rivet;

namespace Rivet
{
  class CMS_2017_I1663958;
}

namespace Rivet_CMS_2017_I1663958
{
  class Histo1DGroup
  {
    private:
      CMS_2017_I1663958* m_an;
      vector<double> m_xbins;
      vector<Histo1DPtr> m_histos;
    public:
      Histo1DGroup(CMS_2017_I1663958* an, const vector<string>& hnames, const vector<double>& xbinranges);

      void fill(double x, double y, double w = 1.);
      void scale(double s, bool diff = true);
      void norm(bool diff = true);

      void gapfractionfromjetpt(Scatter2DPtr hgap, int njet);
  };
}

namespace Rivet
{
  class CMS_2017_I1663958 : public Analysis
  {
    friend class Rivet_CMS_2017_I1663958::Histo1DGroup;
    private:
    Histo1DPtr m_hist_thadpt;
    Histo1DPtr m_hist_thady;
    Histo1DPtr m_hist_tleppt;
    Histo1DPtr m_hist_tlepy;
    Histo1DPtr m_hist_ttpt;
    Histo1DPtr m_hist_tty;
    Histo1DPtr m_hist_ttm;
    Histo1DPtr m_hist_njet;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_hist_njet_ttm;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_hist_njet_thadpt;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_hist_njet_ttpt;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_hist_thady_thadpt;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_hist_ttm_tty;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_hist_thadpt_ttm;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_hist_jetspt;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_hist_jetseta;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_hist_jetsdrjets;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_hist_jetsdrtops;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_hist_njetspt;

    Histo1DPtr m_nhist_thadpt;
    Histo1DPtr m_nhist_thady;
    Histo1DPtr m_nhist_tleppt;
    Histo1DPtr m_nhist_tlepy;
    Histo1DPtr m_nhist_ttm;
    Histo1DPtr m_nhist_ttpt;
    Histo1DPtr m_nhist_tty;
    Histo1DPtr m_nhist_njet;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_nhist_njet_ttm;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_nhist_njet_thadpt;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_nhist_njet_ttpt;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_nhist_thady_thadpt;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_nhist_ttm_tty;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_nhist_thadpt_ttm;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_nhist_jetspt;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_nhist_jetseta;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_nhist_jetsdrjets;
    Rivet_CMS_2017_I1663958::Histo1DGroup* m_nhist_jetsdrtops;
    Scatter2DPtr m_hist_gap1;
    Scatter2DPtr m_hist_gap2;

    Particles m_neutrinos;
    Particles m_leptons;
    Particles m_vetoleptons;
    Jets m_bjets;
    Jets m_ljets;

    Particle m_thad;
    Particles m_thad_decay;
    Particle m_tlep;
    Particles m_tlep_decay;
    Particles m_tt_jets;
    Particles m_additionalobjects;
    Particles m_additionaljets;

    double m_jetptmin = 25.;
    double m_jetetamax = 2.4;
    double m_addjetptmin = 30.;
    double m_addjetetamax = 2.4;
    double m_jetdr = 0.4;
    double m_lepptmin = 30.;
    double m_lepetamax = 2.4;
    double m_vetolepptmin = 15.;
    double m_vetolepetamax = 2.4;
    double m_lepisomax = 0.35;
    double m_lepisodr = 0.4;
    double m_lepdressdr = 0.1;
    double m_phptmin = 15.;
    double m_phetamax = 2.4;
    double m_phisomax = 0.25;
    double m_phisodr = 0.4;

    public:

    CMS_2017_I1663958() : Analysis("CMS_2017_I1663958"), m_thad_decay(3), m_tlep_decay(3), m_tt_jets(4)
    {}
    virtual ~CMS_2017_I1663958()
    {}

    virtual void init()
    {
      const FinalState fs(Cuts::abseta < 6.);
      addProjection(fs, "FS");
      const VisibleFinalState vfs(Cuts::abseta < 6.);
      addProjection(vfs, "vFS");

      VetoedFinalState invisibles(fs);
      invisibles.addVetoOnThisFinalState(vfs);
      addProjection(invisibles, "Invisibles");

      IdentifiedFinalState all_photons(vfs);
      all_photons.acceptId(22);
      IdentifiedFinalState leptons(vfs);
      leptons.acceptIds({11,-11,13,-13});

      DressedLeptons dressed_leptons(all_photons, leptons, m_lepdressdr, Cuts::abseta < m_lepetamax && Cuts::pT > m_vetolepptmin*GeV, true, true);
      addProjection(dressed_leptons, "MyLeptons");

      VetoedFinalState photons(all_photons);
      photons.addVetoOnThisFinalState(dressed_leptons);
      addProjection(photons, "MyPhotons");

      VetoedFinalState isolationparticles(vfs);
      isolationparticles.addVetoOnThisFinalState(dressed_leptons);
      addProjection(isolationparticles, "IsoParticles");

      addProjection(FastJets(vfs, FastJets::ANTIKT, m_jetdr), "Jets");


      m_hist_thadpt = bookHisto1D("d01-x01-y01");
      m_hist_thady = bookHisto1D("d03-x01-y01");
      m_hist_tleppt = bookHisto1D("d05-x01-y01");
      m_hist_tlepy = bookHisto1D("d07-x01-y01");
      m_hist_ttm = bookHisto1D("d09-x01-y01");
      m_hist_ttpt = bookHisto1D("d11-x01-y01");
      m_hist_tty = bookHisto1D("d13-x01-y01");
      m_hist_njet = bookHisto1D("d15-x01-y01");
      m_hist_njet_ttm = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d17-x01-y01", "d18-x01-y01", "d19-x01-y01", "d20-x01-y01"}, {-0.5, 0.5, 1.5, 2.5, 3.5});
      m_hist_njet_thadpt = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d22-x01-y01", "d23-x01-y01", "d24-x01-y01", "d25-x01-y01"}, {-0.5, 0.5, 1.5, 2.5, 3.5});
      m_hist_njet_ttpt = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d27-x01-y01", "d28-x01-y01", "d29-x01-y01", "d30-x01-y01"}, {-0.5, 0.5, 1.5, 2.5, 3.5});
      m_hist_thady_thadpt = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d32-x01-y01", "d33-x01-y01", "d34-x01-y01", "d35-x01-y01"}, {0.0,0.5, 1.0, 1.5, 2.5});
      m_hist_ttm_tty = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d37-x01-y01", "d38-x01-y01", "d39-x01-y01", "d40-x01-y01"}, {300., 450., 625., 850., 2000.});
      m_hist_thadpt_ttm = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d42-x01-y01", "d43-x01-y01", "d44-x01-y01", "d45-x01-y01"}, {0., 90., 180., 270., 800.});
      m_hist_jetspt = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d47-x01-y01", "d48-x01-y01", "d49-x01-y01", "d50-x01-y01", "d51-x01-y01", "d52-x01-y01", "d53-x01-y01", "d54-x01-y01"}, {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5});
      m_hist_jetseta = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d56-x01-y01", "d57-x01-y01", "d58-x01-y01", "d59-x01-y01", "d60-x01-y01", "d61-x01-y01", "d62-x01-y01", "d63-x01-y01"}, {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5});
      m_hist_jetsdrjets = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d65-x01-y01", "d66-x01-y01", "d67-x01-y01", "d68-x01-y01", "d69-x01-y01", "d70-x01-y01", "d71-x01-y01", "d72-x01-y01"}, {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5});
      m_hist_jetsdrtops = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d74-x01-y01", "d75-x01-y01", "d76-x01-y01", "d77-x01-y01", "d78-x01-y01", "d79-x01-y01", "d80-x01-y01", "d81-x01-y01"}, {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5});
      m_hist_njetspt = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d169-x01-y01", "d170-x01-y01", "d171-x01-y01", "d172-x01-y01"}, {0., 40., 60., 80., 120.});

      m_nhist_thadpt = bookHisto1D("d83-x01-y01");
      m_nhist_thady = bookHisto1D("d85-x01-y01");
      m_nhist_tleppt = bookHisto1D("d87-x01-y01");
      m_nhist_tlepy = bookHisto1D("d89-x01-y01");
      m_nhist_ttm = bookHisto1D("d91-x01-y01");
      m_nhist_ttpt = bookHisto1D("d93-x01-y01");
      m_nhist_tty = bookHisto1D("d95-x01-y01");
      m_nhist_njet = bookHisto1D("d97-x01-y01");
      m_nhist_njet_ttm = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d99-x01-y01", "d100-x01-y01", "d101-x01-y01", "d102-x01-y01"}, {-0.5, 0.5, 1.5, 2.5, 3.5});
      m_nhist_njet_thadpt = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d104-x01-y01", "d105-x01-y01", "d106-x01-y01", "d107-x01-y01"}, {-0.5, 0.5, 1.5, 2.5, 3.5});
      m_nhist_njet_ttpt = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d109-x01-y01", "d110-x01-y01", "d111-x01-y01", "d112-x01-y01"}, {-0.5, 0.5, 1.5, 2.5, 3.5});
      m_nhist_thady_thadpt = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d114-x01-y01", "d115-x01-y01", "d116-x01-y01", "d117-x01-y01"}, {0.0,0.5, 1.0, 1.5, 2.5});
      m_nhist_ttm_tty = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d119-x01-y01", "d120-x01-y01", "d121-x01-y01", "d122-x01-y01"}, {300., 450., 625., 850., 2000.});
      m_nhist_thadpt_ttm = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d124-x01-y01", "d125-x01-y01", "d126-x01-y01", "d127-x01-y01"}, {0., 90., 180., 270., 800.});
      m_nhist_jetspt = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d129-x01-y01", "d130-x01-y01", "d131-x01-y01", "d132-x01-y01", "d133-x01-y01", "d134-x01-y01", "d135-x01-y01", "d136-x01-y01"}, {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5});
      m_nhist_jetseta = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d138-x01-y01", "d139-x01-y01", "d140-x01-y01", "d141-x01-y01", "d142-x01-y01", "d143-x01-y01", "d144-x01-y01", "d145-x01-y01"}, {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5});
      m_nhist_jetsdrjets = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d147-x01-y01", "d148-x01-y01", "d149-x01-y01", "d150-x01-y01", "d151-x01-y01", "d152-x01-y01", "d153-x01-y01", "d154-x01-y01"}, {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5});
      m_nhist_jetsdrtops = new Rivet_CMS_2017_I1663958::Histo1DGroup(this, {"d156-x01-y01", "d157-x01-y01", "d158-x01-y01", "d159-x01-y01", "d160-x01-y01", "d161-x01-y01", "d162-x01-y01", "d163-x01-y01"}, {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5});
      m_hist_gap1 = bookScatter2D("d165-x01-y01");
      m_hist_gap2 = bookScatter2D("d167-x01-y01");

    }

    void analyze(const Event& event)
    {
      const double weight = event.weight();
      m_leptons.clear();
      m_vetoleptons.clear();
      m_neutrinos.clear();
      m_bjets.clear();
      m_ljets.clear();
      m_additionalobjects.clear();
      m_additionaljets.clear();

      const ParticleVector& isopars = applyProjection<VetoedFinalState>(event, "IsoParticles").particles();

      const vector<DressedLepton>& dressedleptons = applyProjection<DressedLeptons>(event, "MyLeptons").dressedLeptons();
      for(const DressedLepton& lep : dressedleptons)
      {
        double isolation = accumulate(isopars.begin(), isopars.end(), 0., [&](double iso, const Particle& par) {if(deltaR(lep, par) < m_lepisodr){return iso + par.pt();} else {return iso;}});

        isolation = isolation/lep.pt();

        if(isolation > m_lepisomax) continue;
        if(lep.pt() > m_lepptmin && lep.abseta() < m_lepetamax)
        {
          m_leptons.push_back(lep);
        }
        else if(lep.pt() > m_vetolepptmin && lep.abseta() < m_vetolepetamax)
        {
          m_vetoleptons.push_back(lep);
        }
      }

      const ParticleVector& photons = applyProjection<VetoedFinalState>(event, "MyPhotons").particles(Cuts::abseta < m_phetamax && Cuts::pT > m_phptmin*GeV);
      for(const Particle& ph : photons)
      {
        double isolation = accumulate(isopars.begin(), isopars.end(), 0., [&](double iso, const Particle& par) {if(deltaR(ph, par) < m_phisodr){return iso + par.pt();} else {return iso;}});

        isolation = isolation/ph.pt() - 1.;

        if(isolation > m_phisomax) continue;
        m_additionalobjects.push_back(ph);
      }

      const ParticleVector& invfspars = applyProjection<FinalState>(event, "Invisibles").particles();
      for(const Particle& par : invfspars)
      {
        m_neutrinos.push_back(par);
      }

      const Jets& allJets = applyProjection<FastJets>(event, "Jets").jetsByPt(Cuts::abseta < m_jetetamax && Cuts::pT > m_jetptmin*GeV);
      for(const Jet& jet : allJets)
      {
        //clean jets from leptons
        if(find_if(m_leptons.begin(), m_leptons.end(), [&](const Particle& par){return deltaR(jet, par) < m_jetdr;}) != m_leptons.end()) continue;
        //clean jets from veto leptons
        if(find_if(m_vetoleptons.begin(), m_vetoleptons.end(), [&](const Particle& par){return deltaR(jet, par) < m_jetdr;}) != m_vetoleptons.end()) continue;
        //clean other objects (photons)
        if(find_if(m_additionalobjects.begin(), m_additionalobjects.end(), [&](const Particle& par){return deltaR(jet, par) < m_jetdr;}) != m_additionalobjects.end()) continue;

        if(jet.bTagged())
        {
          m_bjets.push_back(jet);
        }
        else
        {
          m_ljets.push_back(jet);
        }
      }

      //Semi-leptonic reconstruction
      if(m_leptons.size() != 1 || m_vetoleptons.size() != 0 || m_bjets.size() < 2 || m_ljets.size() < 2) {return;}

      FourMomentum nusum = accumulate(m_neutrinos.begin(), m_neutrinos.end(), FourMomentum(0.,0.,0.,0.), [&](FourMomentum& invmom, const Particle& par) {return invmom += par.momentum();});

      FourMomentum wl = nusum + m_leptons[0].momentum();

      double Kmin = numeric_limits<double>::max();
      for(size_t a = 0 ; a <  m_ljets.size() ; ++a)
      {
        const Jet& lja = m_ljets[a];
        for(size_t b = 0 ; b < a ; ++b)
        {
          const Jet& ljb = m_ljets[b];
          FourMomentum wh(lja.momentum() + ljb.momentum());
          for(const Jet& bjh : m_bjets)
          {
            FourMomentum th(wh + bjh.momentum());
            for(const Jet& bjl : m_bjets)
            {
              if(&bjh == &bjl) continue;
              FourMomentum tl(wl + bjl.momentum());

              double K = pow(wh.mass() - 80.4, 2) + pow(th.mass() - 172.5, 2) + pow(tl.mass() - 172.5, 2);
              if(K < Kmin)
              {
                Kmin = K;
                m_thad = Particle(6, th);
                m_thad_decay[0] = Particle(5, bjh);
                m_thad_decay[1] = lja.pt() > ljb.pt() ? Particle(1, lja) : Particle(1, ljb);
                m_thad_decay[2] = lja.pt() <= ljb.pt() ? Particle(1, lja) : Particle(1, ljb);
                m_tlep = Particle(-6, tl);
                m_tlep_decay[0] = Particle(5, bjl);
                m_tlep_decay[1] = m_leptons[0];
                m_tlep_decay[2] = Particle(-1*(m_leptons[0].pdgId()+1), nusum);
              }
            }
          }
        }
      }

      m_tt_jets[0] = m_tlep_decay[0];
      m_tt_jets[1] = m_thad_decay[0];
      m_tt_jets[2] = m_thad_decay[1];
      m_tt_jets[3] = m_thad_decay[2];

      const double eps = 1E-5;
      for(const Jet& jet : m_bjets)
      {
        if(jet.pt() < m_addjetptmin || jet.abseta() > m_addjetetamax) continue;
        if(find_if(m_tt_jets.begin(), m_tt_jets.end(), [&](const Particle& par){return deltaR(jet, par) < eps;}) != m_tt_jets.end()) continue;
        m_additionaljets.push_back(Particle(5, jet.momentum()));
      }
      for(const Jet& jet : m_ljets)
      {
        if(jet.pt() < m_addjetptmin || jet.abseta() > m_addjetetamax) continue;
        if(find_if(m_tt_jets.begin(), m_tt_jets.end(), [&](const Particle& par){return deltaR(jet, par) < eps;}) != m_tt_jets.end()) continue;
        if(jet.cTagged())
        {
          m_additionaljets.push_back(Particle(4, jet.momentum()));
        }
        else
        {
          m_additionaljets.push_back(Particle(1, jet.momentum()));
        }
      }

      sort(m_additionaljets.begin(), m_additionaljets.end(), [](const Particle& ja, const Particle& jb) {return ja.pt() > jb.pt();});

      FourMomentum tt(m_thad.momentum() + m_tlep.momentum());

      m_hist_thadpt->fill(m_thad.pt(), weight);
      m_nhist_thadpt->fill(m_thad.pt(), weight);
      m_hist_thady->fill(abs(m_thad.rapidity()), weight);
      m_nhist_thady->fill(abs(m_thad.rapidity()), weight);
      m_hist_tleppt->fill(m_tlep.pt(), weight);
      m_nhist_tleppt->fill(m_tlep.pt(), weight);
      m_hist_tlepy->fill(abs(m_tlep.rapidity()), weight);
      m_nhist_tlepy->fill(abs(m_tlep.rapidity()), weight);
      m_hist_ttm->fill(tt.mass(), weight);
      m_nhist_ttm->fill(tt.mass(), weight);
      m_hist_ttpt->fill(tt.pt(), weight);
      m_nhist_ttpt->fill(tt.pt(), weight);
      m_hist_tty->fill(abs(tt.rapidity()), weight);
      m_nhist_tty->fill(abs(tt.rapidity()), weight);
      m_hist_njet->fill(min(m_additionaljets.size(), (size_t)5), weight);
      m_nhist_njet->fill(min(m_additionaljets.size(), (size_t)5), weight);
      int njet = min((size_t)3, m_additionaljets.size());
      m_hist_njet_ttm->fill(njet, tt.mass(), weight);
      m_nhist_njet_ttm->fill(njet, tt.mass(), weight);
      m_hist_njet_thadpt->fill(njet, m_thad.pt(), weight);
      m_nhist_njet_thadpt->fill(njet, m_thad.pt(), weight);
      m_hist_njet_ttpt->fill(njet, tt.pt(), weight);
      m_nhist_njet_ttpt->fill(njet, tt.pt(), weight);
      m_hist_thady_thadpt->fill(abs(m_thad.rapidity()), m_thad.pt(), weight);
      m_nhist_thady_thadpt->fill(abs(m_thad.rapidity()), m_thad.pt(), weight);
      m_hist_ttm_tty->fill(tt.mass(), abs(tt.rapidity()), weight);
      m_nhist_ttm_tty->fill(tt.mass(), abs(tt.rapidity()), weight);
      m_hist_thadpt_ttm->fill(m_thad.pt(), tt.mass(), weight);
      m_nhist_thadpt_ttm->fill(m_thad.pt(), tt.mass(), weight);
      int jpos = -1;
      for(const Particles& jets : {m_tt_jets, m_additionaljets})
      {
        for(const Particle& jet : jets)
        {
          jpos++;
          m_hist_jetspt->fill(jpos, jet.pt(), weight);
          m_nhist_jetspt->fill(jpos, jet.pt(), weight);
          m_hist_jetseta->fill(jpos, abs(jet.eta()), weight);
          m_nhist_jetseta->fill(jpos, abs(jet.eta()), weight);
          double drmin = 1E10;
          for(const Particle& tjet : m_tt_jets)
          {
            double dr = deltaR(jet, tjet);
            if(dr > eps && dr < drmin)
            {
              drmin = dr;
            }
          }
          m_hist_jetsdrjets->fill(jpos, drmin, weight);
          m_nhist_jetsdrjets->fill(jpos, drmin, weight);
          m_hist_jetsdrtops->fill(jpos, min(deltaR(jet, m_thad), deltaR(jet, m_tlep)), weight);
          m_nhist_jetsdrtops->fill(jpos, min(deltaR(jet, m_thad), deltaR(jet, m_tlep)), weight);
        }
      }
      for(double ptcut : {30, 50, 75, 100})
      {
        m_hist_njetspt->fill(ptcut , count_if(m_additionaljets.begin(), m_additionaljets.end(), [&ptcut](const Particle& j) {return j.pt() > ptcut;}) , weight);
      }
    }


    virtual void finalize()
    {
      m_hist_jetspt->gapfractionfromjetpt(m_hist_gap1, 1);
      m_hist_jetspt->gapfractionfromjetpt(m_hist_gap2, 2);
      scale(m_hist_thadpt, crossSection()/sumOfWeights());
      scale(m_hist_thady, crossSection()/sumOfWeights());
      scale(m_hist_tleppt, crossSection()/sumOfWeights());
      scale(m_hist_tlepy, crossSection()/sumOfWeights());
      scale(m_hist_ttpt, crossSection()/sumOfWeights());
      scale(m_hist_tty, crossSection()/sumOfWeights());
      scale(m_hist_ttm, crossSection()/sumOfWeights());
      scale(m_hist_njet, crossSection()/sumOfWeights());
      m_hist_njet_ttm->scale(crossSection()/sumOfWeights(), false);
      m_hist_njet_thadpt->scale(crossSection()/sumOfWeights(), false);
      m_hist_njet_ttpt->scale(crossSection()/sumOfWeights(), false);
      m_hist_thady_thadpt->scale(crossSection()/sumOfWeights());
      m_hist_ttm_tty->scale(crossSection()/sumOfWeights());
      m_hist_thadpt_ttm->scale(crossSection()/sumOfWeights());
      m_hist_jetspt->scale(crossSection()/sumOfWeights(), false);
      m_hist_jetseta->scale(crossSection()/sumOfWeights(), false);
      m_hist_jetsdrjets->scale(crossSection()/sumOfWeights(), false);
      m_hist_jetsdrtops->scale(crossSection()/sumOfWeights(), false);
      m_hist_njetspt->scale(crossSection()/sumOfWeights(), false);


      scale(m_nhist_thadpt, 1./m_nhist_thadpt->sumW(false));
      scale(m_nhist_thady, 1./m_nhist_thady->sumW(false));
      scale(m_nhist_tleppt, 1./m_nhist_tleppt->sumW(false));
      scale(m_nhist_tlepy, 1./m_nhist_tlepy->sumW(false));
      scale(m_nhist_ttpt, 1./m_nhist_ttpt->sumW(false));
      scale(m_nhist_tty, 1./m_nhist_tty->sumW(false));
      scale(m_nhist_ttm, 1./m_nhist_ttm->sumW(false));
      scale(m_nhist_njet, 1./m_nhist_njet->sumW(false));
      m_nhist_njet_ttm->norm(false);
      m_nhist_njet_thadpt->norm(false);
      m_nhist_njet_ttpt->norm(false);
      m_nhist_thady_thadpt->norm(true);
      m_nhist_ttm_tty->norm(true);
      m_nhist_thadpt_ttm->norm(true);
      m_nhist_jetspt->norm(false);
      m_nhist_jetseta->norm(false);
      m_nhist_jetsdrjets->norm(false);
      m_nhist_jetsdrtops->norm(false);

    }
    //  Histo2DPtr bookHisto2D(const string& hname)
    //  {
    //    const Scatter3D& refscatter = refData<Scatter3D>(hname);
    //    Histo2DPtr hist( new Histo2D(refscatter, hname));
    //    addAnalysisObject(hist);
    //    return hist;
    //  }



  };
  DECLARE_RIVET_PLUGIN(CMS_2017_I1663958);
}


namespace Rivet_CMS_2017_I1663958
{
  Histo1DGroup::Histo1DGroup(CMS_2017_I1663958* an, const vector<string>& hnames, const vector<double>& xbinranges) : m_an(an), m_xbins(xbinranges)
  {
    for(const string& hname : hnames)
    {
      m_histos.push_back(an->bookHisto1D(hname));
    }
  }

  void Histo1DGroup::fill(double x, double y, double w)
  {
    if(x < m_xbins[0]) {return;}    
    if(x >= m_xbins[m_xbins.size()-1]) {return;}
    int xbin = upper_bound(m_xbins.begin(), m_xbins.end(), x) - m_xbins.begin()-1;
    m_histos[xbin]->fill(y, w);
  }

  void Histo1DGroup::scale(double s, bool diff)
  {
    for(size_t h = 0 ; h < m_histos.size() ; ++h)
    {
      double sc = s;
      if(diff) {sc /= m_xbins[h+1]-m_xbins[h];}
      m_an->scale(m_histos[h], sc);
    }

  }

  void Histo1DGroup::norm(bool diff)
  {
    double sum = 0.;
    for(Histo1DPtr hist : m_histos) {sum+=hist->sumW(false);}
    scale(1./sum, diff);
  }

  void Histo1DGroup::gapfractionfromjetpt(Scatter2DPtr hgap, int njet)
  {
    double total = m_histos[0]->integral();
    int hn = njet+3;
    double totalj = m_histos[hn]->integral();
    double acc = total-totalj;
    for(size_t nb = 0 ; nb < m_histos[hn]->numBins() ; ++nb)
    {
      double gf = acc/total;
      double bl = m_histos[njet+3]->bin(nb).xMin();
      double bh = m_histos[njet+3]->bin(nb).xMax();
      double bc = 0.5*(bh+bl);
      hgap->addPoint(bc, gf, bc-bl, bh-bc, 0., 0.);
      acc += m_histos[njet+3]->bin(nb).area();
    }

  }

}
