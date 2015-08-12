#include "TopMonteCarlo/RivetTop/interface/PartonTop.hh"
//#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

using namespace Rivet;

void PartonTop::project(const Event& e) {
  _theParticles.clear();
  _wDecays1.clear();
  _wDecays2.clear();
  _mode1 = _mode2 = CH_HADRON; // Set default decay mode to full-hadronic

  const double ptmin = 0;
  const double etamin = -MAXDOUBLE, etamax = MAXDOUBLE;

  int nTop = 0;
  bool isTau1 = false, isTau2 = false;
  foreach (GenParticle* p, Rivet::particles(e.genEvent())) {
    const int status = p->status();
    if ( status == 1 ) continue;

    const int pdgId = p->pdg_id();
    if ( PID::isHadron(pdgId) ) continue; // skip hadrons
    if ( pdgId == 22 ) continue; // skip photons

    if ( isZero(p->momentum().perp()) || p->momentum().perp() < ptmin ) continue;
    if ( !inRange(p->momentum().eta(), etamin, etamax) ) continue;

    // Avoid double counting by skipping if particle ID == parent ID
    std::vector<GenParticle*> pps;
    if ( abs(pdgId) == 6 and p->end_vertex() != 0 ) {
      pps = Rivet::particles(p->end_vertex(), HepMC::children);
    }
    else if ( abs(pdgId) != 6 and p->production_vertex() != 0 )
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

    if      ( pdgId ==  6 ) { nTop++; _t1 = Particle(rp); }
    else if ( pdgId == -6 ) { nTop++; _t2 = Particle(rp); }
    else if ( pdgId ==  5 ) _b1 = Particle(rp);
    else if ( pdgId == -5 ) _b2 = Particle(rp);
    else if ( pdgId != 24 && rp.hasAncestor( 24) ) {
      if ( pdgId == -15 ) isTau1 = true;
      else if ( pdgId == -11 ) _mode1 = CH_ELECTRON;
      else if ( pdgId == -13 ) _mode1 = CH_MUON;
      _wDecays1.push_back(rp);
    }
    else if ( pdgId != -24 && rp.hasAncestor(-24) ) {
      if ( pdgId == 15 ) isTau2 = true;
      else if ( pdgId == 11 ) _mode2 = CH_ELECTRON;
      else if ( pdgId == 13 ) _mode2 = CH_MUON;
      _wDecays2.push_back(rp);
    }
  }
  if ( isTau1 ) _mode1 = static_cast<DecayMode>(_mode1+3);
  if ( isTau2 ) _mode2 = static_cast<DecayMode>(_mode2+3);
}
