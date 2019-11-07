#include "Rivet/Analysis.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include <fstream>
#include <iostream>


namespace Rivet {

  class CMS_2019_I1747892 : public Analysis {
  public:

    CMS_2019_I1747892()
      : Analysis("CMS_2019_I1747892")
    { }

    void init() {

        // all the final states we need:
        declare(FinalState(), "FS"); // for event selection
        declare(VisibleFinalState(Cuts::charge != 0 && Cuts::pT > 200*MeV && Cuts::abseta < 2.), "CentralChargedParticles"); // central charged particles
        declare(VisibleFinalState(Cuts::eta > -6.6 && Cuts::eta < -5.2 && Cuts::abspid != PID::MUON), "ForwardParticles"); // forward stable particles

        // three 2D histograms will be filled on generator level, the rest is done in the finalize method
        _hist_nCh_E_tot = new Histo2D(_multiplicityBinEdges, _energyBinEdges,"totalEnergy","totalEnergy");
        _hist_nCh_E_em = new Histo2D(_multiplicityBinEdges, _energyBinEdges,"electromagneticEnergy","electromagneticEnergy");
        _hist_nCh_E_had = new Histo2D(_multiplicityBinEdges, _energyBinEdges,"hadronicEnergy","hadronicEnergy");

        // find the dataFile with the smearing matrices
        std::string foldingFileName = Rivet::findAnalysisRefFile("CMS_2019_I1747892_ForwardFolding.dat");
        if (foldingFileName.empty()) {
            throw std::runtime_error( "reference file for forward folding not found!" );
        }

        // load the folding matrices and bring them to a format to work with
        std::ifstream foldingFile(foldingFileName);
        std::string line;
        while(std::getline(foldingFile, line)) {
            std::istringstream iss(line);
            std::string composition_in, variation_in;
            iss >> composition_in >> variation_in;
            for(int i = 0; i<_nMultiplicityBins; ++i) { 
                for(int j = 0; j<_nEnergyBins; ++j) {
                    for(int k = 0; k<_nMultiplicityBins; ++k) {
                        for(int l = 0; l<_nEnergyBins; ++l) {
                            double nextValue;
                            iss >> nextValue;
                            _forwardFolding[composition_in][variation_in][i][j][k][l] = nextValue;
                        }
                    }
                }
            }
        }
        foldingFile.close();

    }

    void analyze(const Event& event) {

        const double weight = event.weight();

        // ----------------------------
        // Step 1: event selection: reject low mass diffractive events based on largest rapidity gap
        // ----------------------------
        const Particles& AllParticles = apply<FinalState>(event, "FS").particles(cmpMomByRap); // sorted by rapidity
        if (AllParticles.size() < 2) vetoEvent; // need at least two particles to calculate gaps

        double gapCenter = 0.;
        double largestGap = 0.;
        double previousRapidity = 0.;
        bool first = true;

        foreach(const Particle& p, AllParticles) {
            if (first) { // First particle
                first = false;
                previousRapidity = p.rapidity();
            } else {
                double gap = fabs(p.rapidity()-previousRapidity);
                if (gap > largestGap) {
                   largestGap = gap; // largest gap
                   gapCenter = (p.rapidity()+previousRapidity)/2.; // find the center of the gap to separate the X and Y systems.
                }
                previousRapidity = p.rapidity();
            }
        }

        // masses of the two subsystems
        FourMomentum mxFourVector, myFourVector;
        foreach(const Particle& p, AllParticles) {
            ((p.rapidity() > gapCenter) ? mxFourVector : myFourVector) += p.momentum();
        }
        const double xiX = mxFourVector.mass2()/sqr(sqrtS());
        const double xiY = myFourVector.mass2()/sqr(sqrtS());
        const double xi = max(xiX,xiY);

        //need minimal mass to be visible to the detector
        if (xi < 1e-6) vetoEvent;

        // ----------------------------
        // Step 2: central multiplicity
        // ----------------------------

        const Particles& CentralChargedParticles = apply<FinalState>(event, "CentralChargedParticles").particles();
        int nCharged_abseta2 = CentralChargedParticles.size();
            if (nCharged_abseta2==0) vetoEvent; // at least one charged particle in central acceptance

        // ----------------------------
        // Step 3: calculate energies in forward acceptance
        // ----------------------------
        const Particles& ForwardParticles = apply<FinalState>(event, "ForwardParticles").particles();

        // we look at three energy compositions
        double totEnergy = 0.; // all visible particles except muons
        double emEnergy = 0.; // only e and gamma
        double hadEnergy = 0.; // only hadrons

        foreach(const Particle& p, ForwardParticles) {
            totEnergy += p.energy();
            if ( p.abspid() == 11 || p.abspid() == 22 || p.abspid() == 111){
                emEnergy += p.energy();
            }
            if ( p.abspid() != 11 && p.abspid() != 22 && p.abspid() != 111){
                hadEnergy += p.energy();
            }
        }
	

        // ----------------------------
        // Fill histograms
        // ----------------------------
        _hist_nCh_E_tot->fill(nCharged_abseta2, totEnergy, weight);
        _hist_nCh_E_em->fill(nCharged_abseta2, emEnergy, weight);
        _hist_nCh_E_had->fill(nCharged_abseta2, hadEnergy, weight);

        // ----------------------------
        // Done with generator level
        // ----------------------------

    }

    void finalize() {


        // ----------------------------
        // post-processing of the generator-level histograms
        // ----------------------------

        // ----------------------------
        // Step 1: Apply 4-dimensional smearing matrix to the generator-level distribution
        //     for the central value as well as four systematic variations
        // ----------------------------

        std::map <std::string , std::map < std::string , TwoDimMatrix > > foldedHistogram;
        for(const auto& variation : _variations) { // we do this for all systematic variations
            for (unsigned int iNreco=0; iNreco<_nMultiplicityBins; ++iNreco) {
                for (unsigned int iEreco=0; iEreco<_nEnergyBins; ++iEreco) {
                    // we do this for all three energy compositions
                    double smearedValue_tot = 0.;
                    double smearedValue_em = 0.;
                    double smearedValue_had = 0.;
                    for (unsigned int iNtrue=0; iNtrue<_nMultiplicityBins; ++iNtrue) {
                        for (unsigned int iEtrue=0; iEtrue<_nEnergyBins; ++iEtrue) {
                            int binIndex = _hist_nCh_E_tot->binIndexAt((_multiplicityBinEdges[iNtrue]+_multiplicityBinEdges[iNtrue+1])/2., (_energyBinEdges[iEtrue]+_energyBinEdges[iEtrue+1])/2.);
                            double generatorData_tot = _hist_nCh_E_tot->bin(binIndex).volume();
                            double generatorData_em = _hist_nCh_E_em->bin(binIndex).volume();
                            double generatorData_had = _hist_nCh_E_had->bin(binIndex).volume();
                            smearedValue_tot += generatorData_tot*_forwardFolding["total"][variation][iNtrue][iEtrue][iNreco][iEreco];
                            smearedValue_em += generatorData_em*_forwardFolding["electromagnetic"][variation][iNtrue][iEtrue][iNreco][iEreco];
                            smearedValue_had += generatorData_had*_forwardFolding["hadronic"][variation][iNtrue][iEtrue][iNreco][iEreco];
                        }
                    }
                    // these are the folded distributions on detector level
                    foldedHistogram["total"][variation][iNreco][iEreco] = smearedValue_tot;
                    foldedHistogram["electromagnetic"][variation][iNreco][iEreco] = smearedValue_em;
                    foldedHistogram["hadronic"][variation][iNreco][iEreco] = smearedValue_had;
                }
            }
        }


        // ----------------------------
        // Step 2: Calculate average energies for each distribution 
        // ----------------------------

        std::map < std::string, std::map < std::string , double[_nReducedMultiplicityBins]  > > AverageEnergy, ErrorOnTheMean;

        for(const auto& composition : _compositions) { // we do this for all energy compositions
            for(const auto& variation : _variations) { // we do this for all systematic variations

                // only run over limited range to compare to data
                for (unsigned int iBin=0; iBin<_nReducedMultiplicityBins; ++iBin) {

                    double nEvents=0.;
                    double averageEnergy=0.;
                    double sigma2=0.;
                    double sigma=0.;

                    for (unsigned int j=1; j<_nEnergyBins; ++j){
                        double binCenter = (_energyBinEdges[j]+_energyBinEdges[j+1])/2.;
                        double binContent = foldedHistogram[composition][variation][iBin+1][j];
                        averageEnergy += binContent*binCenter;
                        nEvents += binContent;
                    }
                    if (nEvents != 0) averageEnergy /= nEvents;

                    for (unsigned int j=1; j<_nEnergyBins; ++j){
                        double binCenter = (_energyBinEdges[j]+_energyBinEdges[j+1])/2.;
                        double binContent = foldedHistogram[composition][variation][iBin+1][j];
                        sigma2 += binContent * (binCenter-averageEnergy)*(binCenter-averageEnergy);
                    }
                    if (nEvents >1) sigma = std::sqrt(sigma2/nEvents/(nEvents-1));

                    AverageEnergy[composition][variation][iBin] = averageEnergy;
                    ErrorOnTheMean[composition][variation][iBin] = sigma;
                }
            }
        }


        // ----------------------------
        // Step 3: Calculate shape normalized distribution of total energy and the ratio of em and had energies
        // ----------------------------

        std::map < std::string , double[_nReducedMultiplicityBins] > shapeAnalysis, shapeAnalysis_uncertainty, emHadRatio, emHadRatio_uncertainty;

        for(const auto& variation : _variations) { // we do this for all systematic variations
            for (unsigned int iBin=0; iBin<_nReducedMultiplicityBins; ++iBin) {

                shapeAnalysis[variation][iBin] = AverageEnergy["total"][variation][iBin] / AverageEnergy["total"][variation][0];
                shapeAnalysis_uncertainty[variation][iBin] = std::sqrt( std::pow(ErrorOnTheMean["total"][variation][iBin]/AverageEnergy["total"][variation][0],2)
                        +std::pow(ErrorOnTheMean["total"][variation][0]*AverageEnergy["total"][variation][iBin]/AverageEnergy["total"][variation][0]/AverageEnergy["total"][variation][0],2) ) ;

                emHadRatio[variation][iBin] = AverageEnergy["electromagnetic"][variation][iBin] / AverageEnergy["hadronic"][variation][iBin];
                emHadRatio_uncertainty[variation][iBin] = std::sqrt( std::pow(ErrorOnTheMean["electromagnetic"][variation][iBin]/AverageEnergy["hadronic"][variation][iBin],2)
                        +std::pow(ErrorOnTheMean["hadronic"][variation][iBin]*AverageEnergy["electromagnetic"][variation][iBin]/AverageEnergy["hadronic"][variation][iBin]/AverageEnergy["hadronic"][variation][0],2) ) ;
            }
        }


        // ----------------------------
        // Step 4: Calculate range of systematic variations and add this as uncertainty to the central distribution
        // ----------------------------

        std::map < std::string, double[_nReducedMultiplicityBins] > SystErrorOnTheMeanUp, SystErrorOnTheMeanDown;

        for(const auto& composition : _compositions) { // we do this for all energy compositions
            for (unsigned int iBin=0; iBin<_nReducedMultiplicityBins; ++iBin) {
                std::vector <double> tmp;
                for(const auto& variation : _variations) {
                    tmp.push_back(AverageEnergy[composition][variation][iBin]);
                }
                double systValueMax = *max_element(tmp.begin(), tmp.end()); // max value of systematic variations
                double systValueMin = *min_element(tmp.begin(), tmp.end()); // min value of systematic variations

                double centralValue = AverageEnergy[composition]["central"][iBin];
                double statError = ErrorOnTheMean[composition]["central"][iBin];
            
                // final uncertainty on the central value: square sum of statistical uncertainty, systematic variation, and additional 1% due to the binning
                if (statError>0. || (systValueMax != centralValue && systValueMin != centralValue)){
                    SystErrorOnTheMeanUp[composition][iBin] = std::sqrt( std::pow(statError,2) + std::pow((systValueMax-centralValue),2)+std::pow(0.01*centralValue,2));
                    SystErrorOnTheMeanDown[composition][iBin] = std::sqrt( std::pow(statError,2) + std::pow((centralValue-systValueMin),2) + std::pow(0.01*centralValue,2) );
                } else {
                    SystErrorOnTheMeanUp[composition][iBin] = 0.;
                    SystErrorOnTheMeanDown[composition][iBin] = 0.;
                }
            }
        }


        // do the same for the shape analysis and em/had ratio
        for (unsigned int iBin=0; iBin<_nReducedMultiplicityBins; iBin++) {
            std::vector <double> tmp;
            for(const auto& variation : _variations) {
                tmp.push_back(shapeAnalysis[variation][iBin]);
            }
            double systValueMax = *max_element(tmp.begin(), tmp.end()); // max value of systematic variations
            double systValueMin = *min_element(tmp.begin(), tmp.end()); // min value of systematic variations

            double centralValue = shapeAnalysis["central"][iBin];
            double statError = shapeAnalysis_uncertainty["central"][iBin];
            // final uncertainty on the central value: square sum of statistical uncertainty, systematic variation, and additional 1% due to the binning
            if (statError>0. || (systValueMax != centralValue && systValueMin != centralValue)){
                SystErrorOnTheMeanUp["shape"][iBin] = std::sqrt(std::pow(statError,2) + std::pow((systValueMax-centralValue),2) + std::pow(0.01*centralValue,2) );
                SystErrorOnTheMeanDown["shape"][iBin] = std::sqrt(std::pow(statError,2) + std::pow((centralValue-systValueMin),2) + std::pow(0.01*centralValue,2) );
            } else {
                SystErrorOnTheMeanUp["shape"][iBin] = 0.;
                SystErrorOnTheMeanDown["shape"][iBin] = 0.;
            }
        }

        for (unsigned int iBin=0; iBin<_nReducedMultiplicityBins; iBin++) {
            std::vector <double> tmp;
            for(const auto& variation : _variations) {
                tmp.push_back(emHadRatio[variation][iBin]);
            }
            double systValueMax = *max_element(tmp.begin(), tmp.end()); // max value of systematic variations
            double systValueMin = *min_element(tmp.begin(), tmp.end()); // min value of systematic variations

            double centralValue = emHadRatio["central"][iBin];
            double statError = emHadRatio_uncertainty["central"][iBin];
            // final uncertainty on the central value: square sum of statistical uncertainty, systematic variation, and additional 1% due to the binning
            if (statError>0. || (systValueMax != centralValue && systValueMin != centralValue)){
                SystErrorOnTheMeanUp["emhadratio"][iBin] = std::sqrt(std::pow(statError,2) + std::pow((systValueMax-centralValue),2)+std::pow(0.01*centralValue,2) );
                SystErrorOnTheMeanDown["emhadratio"][iBin] = std::sqrt(std::pow(statError,2) + std::pow((centralValue-systValueMin),2)+std::pow(0.01*centralValue,2) );
            } else {
                SystErrorOnTheMeanUp["emhadratio"][iBin] = 0.;
                SystErrorOnTheMeanDown["emhadratio"][iBin] = 0.;
            }
        }

        // ----------------------------
        // Step 5: Create and save distributions as Scatter2Ds
        // ----------------------------

	// d01-x01-y01: average total energy, central value and variations
        Scatter2DPtr _h_nCh_E_tot_smeared_profile = bookScatter2D(1,1,1);
        Scatter2DPtr _h_nCh_E_tot_smeared_profile_var1 = bookScatter2D("d01-x01-y01_var1");
        Scatter2DPtr _h_nCh_E_tot_smeared_profile_var2 = bookScatter2D("d01-x01-y01_var2");
        Scatter2DPtr _h_nCh_E_tot_smeared_profile_var3 = bookScatter2D("d01-x01-y01_var3");
        Scatter2DPtr _h_nCh_E_tot_smeared_profile_var4 = bookScatter2D("d01-x01-y01_var4");

	// d01-x01-y02: shape normalized average total energy, central value and variations
        Scatter2DPtr _h_nCh_E_tot_normalized_smeared_profile = bookScatter2D(1,1,2);
        Scatter2DPtr _h_nCh_E_tot_normalized_smeared_profile_var1 = bookScatter2D("d01-x01-y02_var1");
        Scatter2DPtr _h_nCh_E_tot_normalized_smeared_profile_var2 = bookScatter2D("d01-x01-y02_var2");
        Scatter2DPtr _h_nCh_E_tot_normalized_smeared_profile_var3 = bookScatter2D("d01-x01-y02_var3");
        Scatter2DPtr _h_nCh_E_tot_normalized_smeared_profile_var4 = bookScatter2D("d01-x01-y02_var4");

	// d02-x01-y01: average electromagnetic energy, central value and variations
        Scatter2DPtr _h_nCh_E_em_smeared_profile = bookScatter2D(2,1,1);
        Scatter2DPtr _h_nCh_E_em_smeared_profile_var1 = bookScatter2D("d02-x01-y01_var1");
        Scatter2DPtr _h_nCh_E_em_smeared_profile_var2 = bookScatter2D("d02-x01-y01_var2");
        Scatter2DPtr _h_nCh_E_em_smeared_profile_var3 = bookScatter2D("d02-x01-y01_var3");
        Scatter2DPtr _h_nCh_E_em_smeared_profile_var4 = bookScatter2D("d02-x01-y01_var4");

	// d03-x01-y01: average hadronic energy, central value and variations
        Scatter2DPtr _h_nCh_E_had_smeared_profile = bookScatter2D(3,1,1);
        Scatter2DPtr _h_nCh_E_had_smeared_profile_var1 = bookScatter2D("d03-x01-y01_var1");
        Scatter2DPtr _h_nCh_E_had_smeared_profile_var2 = bookScatter2D("d03-x01-y01_var2");
        Scatter2DPtr _h_nCh_E_had_smeared_profile_var3 = bookScatter2D("d03-x01-y01_var3");
        Scatter2DPtr _h_nCh_E_had_smeared_profile_var4 = bookScatter2D("d03-x01-y01_var4");

	// d04-x01-y01: ratio of average electromagnetic to hadronic energy, central value and variations
        Scatter2DPtr _h_nCh_E_emhadRatio_smeared_profile = bookScatter2D(4,1,1);
        Scatter2DPtr _h_nCh_E_emhadRatio_smeared_profile_var1 = bookScatter2D("d04-x01-y01_var1");
        Scatter2DPtr _h_nCh_E_emhadRatio_smeared_profile_var2 = bookScatter2D("d04-x01-y01_var2");
        Scatter2DPtr _h_nCh_E_emhadRatio_smeared_profile_var3 = bookScatter2D("d04-x01-y01_var3");
        Scatter2DPtr _h_nCh_E_emhadRatio_smeared_profile_var4 = bookScatter2D("d04-x01-y01_var4");


        for (unsigned int iBin=0; iBin<_nReducedMultiplicityBins; iBin++) {

            double multiplicityCenter = (_multiplicityBinEdges[iBin+1]+_multiplicityBinEdges[iBin+2])/2.;
            double multiplicityHalfWidth = multiplicityCenter-_multiplicityBinEdges[iBin+1];


	    // d01-x01-y01: average total energy, central value and variations
            _h_nCh_E_tot_smeared_profile->addPoint(multiplicityCenter, AverageEnergy["total"]["central"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                SystErrorOnTheMeanDown["total"][iBin], SystErrorOnTheMeanUp["total"][iBin]);
	    _h_nCh_E_tot_smeared_profile_var1->addPoint(multiplicityCenter, AverageEnergy["total"]["variation1"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                ErrorOnTheMean["total"]["variation1"][iBin], ErrorOnTheMean["total"]["variation1"][iBin]);
	    _h_nCh_E_tot_smeared_profile_var2->addPoint(multiplicityCenter, AverageEnergy["total"]["variation2"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                ErrorOnTheMean["total"]["variation2"][iBin], ErrorOnTheMean["total"]["variation2"][iBin]);
	    _h_nCh_E_tot_smeared_profile_var3->addPoint(multiplicityCenter, AverageEnergy["total"]["variation3"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                ErrorOnTheMean["total"]["variation3"][iBin], ErrorOnTheMean["total"]["variation3"][iBin]);
	    _h_nCh_E_tot_smeared_profile_var4->addPoint(multiplicityCenter, AverageEnergy["total"]["variation4"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                ErrorOnTheMean["total"]["variation4"][iBin], ErrorOnTheMean["total"]["variation4"][iBin]);


            // d01-x01-y02: shape normalized average total energy, central value and variations
            _h_nCh_E_tot_normalized_smeared_profile->addPoint(multiplicityCenter, shapeAnalysis["central"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                SystErrorOnTheMeanDown["shape"][iBin], SystErrorOnTheMeanUp["shape"][iBin]);
           _h_nCh_E_tot_normalized_smeared_profile_var1->addPoint(multiplicityCenter, shapeAnalysis["variation1"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                shapeAnalysis_uncertainty["variation1"][iBin], shapeAnalysis_uncertainty["variation1"][iBin]);
           _h_nCh_E_tot_normalized_smeared_profile_var2->addPoint(multiplicityCenter, shapeAnalysis["variation3"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                shapeAnalysis_uncertainty["variation2"][iBin], shapeAnalysis_uncertainty["variation2"][iBin]);
           _h_nCh_E_tot_normalized_smeared_profile_var3->addPoint(multiplicityCenter, shapeAnalysis["variation3"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                shapeAnalysis_uncertainty["variation3"][iBin], shapeAnalysis_uncertainty["variation3"][iBin]);
           _h_nCh_E_tot_normalized_smeared_profile_var4->addPoint(multiplicityCenter, shapeAnalysis["variation4"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                shapeAnalysis_uncertainty["variation4"][iBin], shapeAnalysis_uncertainty["variation4"][iBin]);


            // d02-x01-y01: average electromagnetic energy, central value and variations
            _h_nCh_E_em_smeared_profile->addPoint(multiplicityCenter, AverageEnergy["electromagnetic"]["central"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                SystErrorOnTheMeanDown["electromagnetic"][iBin], SystErrorOnTheMeanUp["electromagnetic"][iBin]);
	    _h_nCh_E_em_smeared_profile_var1->addPoint(multiplicityCenter, AverageEnergy["electromagnetic"]["variation1"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                ErrorOnTheMean["electromagnetic"]["variation1"][iBin], ErrorOnTheMean["electromagnetic"]["variation1"][iBin]);
	    _h_nCh_E_em_smeared_profile_var2->addPoint(multiplicityCenter, AverageEnergy["electromagnetic"]["variation2"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                ErrorOnTheMean["electromagnetic"]["variation2"][iBin], ErrorOnTheMean["electromagnetic"]["variation2"][iBin]);
	    _h_nCh_E_em_smeared_profile_var3->addPoint(multiplicityCenter, AverageEnergy["electromagnetic"]["variation3"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                ErrorOnTheMean["electromagnetic"]["variation3"][iBin], ErrorOnTheMean["electromagnetic"]["variation3"][iBin]);
	    _h_nCh_E_em_smeared_profile_var4->addPoint(multiplicityCenter, AverageEnergy["electromagnetic"]["variation4"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                ErrorOnTheMean["electromagnetic"]["variation4"][iBin], ErrorOnTheMean["electromagnetic"]["variation4"][iBin]);


            // d03-x01-y01: average hadronic energy, central value and variations
            _h_nCh_E_had_smeared_profile->addPoint(multiplicityCenter, AverageEnergy["hadronic"]["central"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                SystErrorOnTheMeanDown["hadronic"][iBin], SystErrorOnTheMeanUp["hadronic"][iBin]);
	    _h_nCh_E_had_smeared_profile_var1->addPoint(multiplicityCenter, AverageEnergy["hadronic"]["variation1"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                ErrorOnTheMean["hadronic"]["variation1"][iBin], ErrorOnTheMean["hadronic"]["variation1"][iBin]);
	    _h_nCh_E_had_smeared_profile_var2->addPoint(multiplicityCenter, AverageEnergy["hadronic"]["variation2"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                ErrorOnTheMean["hadronic"]["variation2"][iBin], ErrorOnTheMean["hadronic"]["variation2"][iBin]);
	    _h_nCh_E_had_smeared_profile_var3->addPoint(multiplicityCenter, AverageEnergy["hadronic"]["variation3"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                ErrorOnTheMean["hadronic"]["variation3"][iBin], ErrorOnTheMean["hadronic"]["variation3"][iBin]);
	    _h_nCh_E_had_smeared_profile_var4->addPoint(multiplicityCenter, AverageEnergy["hadronic"]["variation4"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                ErrorOnTheMean["hadronic"]["variation4"][iBin], ErrorOnTheMean["hadronic"]["variation4"][iBin]);


	    // d04-x01-y01: ratio of average electromagnetic to hadronic energy, central value and variations
            _h_nCh_E_emhadRatio_smeared_profile->addPoint(multiplicityCenter, emHadRatio["central"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                SystErrorOnTheMeanDown["emhadratio"][iBin], SystErrorOnTheMeanUp["emhadratio"][iBin]);
	    _h_nCh_E_emhadRatio_smeared_profile_var1->addPoint(multiplicityCenter, emHadRatio["variation1"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                emHadRatio_uncertainty["variation1"][iBin], emHadRatio_uncertainty["variation1"][iBin]);
	    _h_nCh_E_emhadRatio_smeared_profile_var2->addPoint(multiplicityCenter, emHadRatio["variation2"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                emHadRatio_uncertainty["variation2"][iBin], emHadRatio_uncertainty["variation2"][iBin]);
	    _h_nCh_E_emhadRatio_smeared_profile_var3->addPoint(multiplicityCenter, emHadRatio["variation3"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                emHadRatio_uncertainty["variation3"][iBin], emHadRatio_uncertainty["variation3"][iBin]);
	    _h_nCh_E_emhadRatio_smeared_profile_var4->addPoint(multiplicityCenter, emHadRatio["variation4"][iBin],
                                                multiplicityHalfWidth, multiplicityHalfWidth,
                                                emHadRatio_uncertainty["variation4"][iBin], emHadRatio_uncertainty["variation4"][iBin]);
        }

    }

  private:

    // histograms needed globally
    Histo2D* _hist_nCh_E_tot;
    Histo2D* _hist_nCh_E_em;
    Histo2D* _hist_nCh_E_had;

    // bins for the analysis, including over- and underflow bins
    const static int _nMultiplicityBins = 21;
    const static int _nReducedMultiplicityBins = _nMultiplicityBins-6;
    const static int _nEnergyBins = 47;
    const std::vector<double> _multiplicityBinEdges = {-0.5,0.5,9.5,19.5,29.5,39.5,49.5,59.5,69.5,79.5,89.5,99.5,
                                                            109.5,119.5,129.5,139.5,149.5,159.5,169.5,179.5,189.5,199.5};
    const std::vector<double> _energyBinEdges = {-25,0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,
                                                    375,400,425,450,475,500,525,550,575,600,650,700,750,800,850,900,950,
                                                    1000,1100,1200,1300,1400,1500,1750,2000,2250,2500,3000,3500,5000,7000,10000};

    // a couple of definitions to make life easier
    typedef double TwoDimMatrix[_nMultiplicityBins][_nEnergyBins];
    typedef double FourDimMatrix[_nMultiplicityBins][_nEnergyBins][_nMultiplicityBins][_nEnergyBins];

    const std::array<std::string, 5> _variations = {{"central","variation1","variation2","variation3","variation4"}};
    const std::array<std::string, 3> _compositions = {{"total","electromagnetic","hadronic"}};

    // the structure to hold the forward folding information
    std::map < std::string, std::map < std::string , FourDimMatrix > > _forwardFolding;

  };

  DECLARE_RIVET_PLUGIN(CMS_2019_I1747892);
}
