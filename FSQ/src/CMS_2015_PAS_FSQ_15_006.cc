/// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {
    
    ///forward energy flow at 13 TeV with CMS
    
    class CMS_FSQ_15_006 : public Analysis {
    public:
        
    CMS_FSQ_15_006()
        : Analysis("CMS_FSQ_15_006")
       {    }
        
//-----------------------------------------------------------------
        
    void init() {

        _noe_inel = 0.;
        _noe_nsd = 0.;
        _noe_bsc = 0.;
        _noe_sd = 0.;
        _noe_nsd_sd = 0.;
        
        addProjection(FinalState(),"FS");
        
        const ChargedFinalState cfsBSCplus(3.9, 4.4, 0*MeV);
        addProjection(cfsBSCplus, "cfsBSCplus");
        
        const ChargedFinalState cfsBSCminus(-4.4, -3.9, 0*MeV);
        addProjection(cfsBSCminus, "cfsBSCminus");
        
        _h_inel = bookHisto1D(1, 1, 1);
        _h_nsd  = bookHisto1D(2, 1, 1);
        _h_et   = bookHisto1D(3, 1, 1);
        _h_sd   = bookHisto1D(4, 1, 1);
        
     }
        
//-----------------------------------------------------------

    void analyze(const Event& event) {
        
    const double weight = event.weight();

    double ybeam = 9.54;
    bool bscplus  = true;
    bool bscminus = true;

     const ChargedFinalState& cfsBSCplus = applyProjection<ChargedFinalState>(event, "cfsBSCplus");
    if (cfsBSCplus.empty()) bscplus = false;
        const ChargedFinalState& cfsBSCminus = applyProjection<ChargedFinalState>(event, "cfsBSCminus");
        if (cfsBSCminus.empty()) bscminus = false;

        const FinalState& fs = applyProjection<FinalState>(event, "FS");

        const ParticleVector particlesByRapidity = fs.particlesByRapidity();
        const size_t num_particles = particlesByRapidity.size();

         vector<double> gaps;
         vector<double> midpoints;
        for (size_t ip = 1; ip < num_particles; ++ip) {
            const Particle& p1 = particlesByRapidity[ip-1];
            const Particle& p2 = particlesByRapidity[ip];
            const double gap = p2.momentum().rapidity() - p1.momentum().rapidity();
            const double mid = (p2.momentum().rapidity() + p1.momentum().rapidity()) / 2.;
            gaps.push_back(gap);
            midpoints.push_back(mid);
            
        }

        int imid = std::distance(gaps.begin(), max_element(gaps.begin(), gaps.end()));
        double gapcenter = midpoints[imid];

        FourMomentum MxFourVector(0.,0.,0.,0.);
        FourMomentum MyFourVector(0.,0.,0.,0.);
        
        foreach(const Particle& p, fs.particlesByRapidity()) {
            if (p.momentum().rapidity() > gapcenter) {
                MxFourVector += p.momentum();
            } else {
                MyFourVector += p.momentum();
            }
        }

        double Mx = MxFourVector.mass();
        double My = MyFourVector.mass();

        double xix = (Mx * Mx) / (sqrtS()/GeV * sqrtS()/GeV);
        double xiy = (My * My) / (sqrtS()/GeV * sqrtS()/GeV);
        double xi_sd = std::max(xix, xiy);

        bool inel = false;

        if (xi_sd > 1E-6) {
            inel = true;
            ++_noe_inel;
        }

        bool nsd = false;
        bool bsc = false;
        bool sd = false;
        int nplus  = 0;
        int nminus = 0;
        
        if (bscplus && bscminus) {
            ++_noe_bsc;
            bsc = true;
        }
        
        foreach(const Particle& p, fs.particlesByRapidity()) {
            double eta = p.momentum().eta();
            //double energy = p.momentum().E();
            double tenergy = p.momentum().Et();
            if (abs(p.pid()) >= 12 && abs(p.pid()) <= 16) continue;
            if (eta > 2.866 && eta < 5.205) ++nplus;
            if (eta < -2.866 && eta > -5.205) ++nminus;
            
            if (bsc) {
                _h_et->fill(eta-ybeam, tenergy*weight);
            }
        }

        if (nminus > 0 && nplus > 0) {
            nsd = true;
            ++_noe_nsd;
        }
        
        foreach(const Particle& p, fs.particlesByRapidity()) {
            double eta = p.momentum().eta();
            double energy = p.momentum().E();
            if (inel) {
                _h_inel->fill(abs(eta), energy*weight);
            }
            if (nsd) {
                _h_nsd->fill(abs(eta), energy*weight);
            }
        }

        //////////////
        //sd selection
        //////////////
        float etamaxminus = -3.152;
        float etaminminus = -5.205;
        float etamin   = 3.152;
        float etamax   = 5.205;
        float emin     = 5.;
        bool StableParticleEnergyCutMinus = false;
        bool StableParticleEnergyCutPlus  = false;
        
        foreach(const Particle& p, fs.particlesByRapidity()) {
            double eta    = p.momentum().eta();
            double energy = p.momentum().E();
            //check whether this particle is stable according to the generator
            if (p.genParticle()->status() == 1) {
                if (eta >= etaminminus && eta <= etamaxminus && energy > emin)  StableParticleEnergyCutMinus=true;
                if (eta >= etamin && eta <= etamax && energy > emin)            StableParticleEnergyCutPlus=true;
            }
        }//decision on SD-event selection
        
        //in order to select SD-enhanced events use the following condition
        if((StableParticleEnergyCutPlus && !StableParticleEnergyCutMinus) ||
           (!StableParticleEnergyCutPlus && StableParticleEnergyCutMinus)) {
            foreach(const Particle& p, fs.particlesByRapidity()) {
                 double eta    = p.momentum().eta();
                 double energy = p.momentum().E();
                 if(StableParticleEnergyCutPlus && !StableParticleEnergyCutMinus){
                     if (abs(eta) >= etamin && abs(eta) <= etamax) {
                         if (eta > 0) _h_sd->fill(abs(eta), energy*weight);
                     }
                     //for CASTOR
                     else  _h_sd->fill(abs(eta), energy*weight*0.5);
                 }
                 if(!StableParticleEnergyCutPlus && StableParticleEnergyCutMinus ){
                     if (abs(eta) >= etamin && abs(eta) <= etamax) {
                         if (eta < 0) _h_sd->fill(abs(eta), energy*weight);
                     }
                     //for CASTOR
                     else  _h_sd->fill(abs(eta), energy*weight*0.5);
                 }
             }
            sd = true;
            ++_noe_sd;
        }//SD-selection
        if (nsd && sd ) ++_noe_nsd_sd;
    }//event

//---------------------------------------------------------------

        void finalize() {

            scale(_h_inel, (1./(2.*_noe_inel)));
            scale(_h_nsd, (1./(2.*_noe_nsd)));
            scale(_h_et, 1./_noe_bsc);
            scale(_h_sd, 1./_noe_sd);
            MSG_INFO( "Number of events of INEL : " << _noe_inel );
            MSG_INFO( "Number of events of NSD : " << _noe_nsd );
            MSG_INFO( "Number of events of SD : " << _noe_sd );
            MSG_INFO( "Number of events of NSD and SD contribution :" << _noe_nsd_sd );

        }

    private:
        
        Histo1DPtr _h_inel;
        Histo1DPtr _h_nsd;
        Histo1DPtr _h_et;
        Histo1DPtr _h_sd;
        double _noe_inel;
        double _noe_nsd;
        double _noe_bsc;
        double _noe_sd;
        double _noe_nsd_sd;
    };
    
    DECLARE_RIVET_PLUGIN(CMS_FSQ_15_006);

 }
