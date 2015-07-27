#include "GeneratorInterface/RivetTop/interface/CMSGenParticle.hh"
//#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

using namespace Rivet;

bool CMSGenParticle::isFromResonance(GenParticle* p) const {
  const int aid = abs(p->pdg_id());
  if ( _vetoIds.find(aid) == _vetoIds.end() and isParton(aid) ) return false;

  if ( p->production_vertex() == 0 ) return false;
  foreach ( GenParticle* pp,  Rivet::particles(p->production_vertex(), HepMC::parents) ) {
    if ( isResonance(pp->pdg_id()) ) return true;
    if ( isFromResonance(pp) ) return true;
  }

  return false;
}

void CMSGenParticle::project(const Event& e) {
  _theParticles.clear();

  foreach (GenParticle* p, Rivet::particles(e.genEvent())) {
    const int status = p->status();
    if ( status == 1 ) continue;

    const int pdgId = p->pdg_id();
    if ( _vetoIds.find(pdgId) != _vetoIds.end() ) continue;
    if ( isResonance(pdgId) ) continue;
    // if (isResonance(id) && (particle->status() == 3 || particle->status() == 22) ) <- status code in CMSSW, but this is not allowed in RIVET
    if ( isFromResonance(p) ) continue;

    //if ( isZero(p->momentum().perp()) || p->momentum().perp() < ptmin ) continue;
    //if ( !inRange(p->momentum().eta(), etamin, etamax) ) continue;
    
    // Build Rivet::Particle
    _theParticles.push_back(Particle(*p));
  }

}
