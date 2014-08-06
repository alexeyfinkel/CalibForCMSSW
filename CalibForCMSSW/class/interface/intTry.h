//this is going to be a first draft/brainstorming area for a class that interfaces MN calibration with CMSSW.



// Standard Library
#include <map>  // std::map
#include <string>  // std::string
#include <utility>  // std::pair
#include <vector>  // std::vector
#include <stdio>	//FILE

//Root
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>




typedef std::map< std::pair<int,int>, TH1D* > xtlHistMap;

class EECalibration 
{
	public:
		//default constructor
		EECalibration(){};
		
		//Normal constructor
		EECalibration(TChain* inputChain, std::string outputDir);
			/*	this needs to take in the following:
			path to output directory, output file name
			?possibly names of the variables in the provided ntuples?
			?possibly optional numbner of iterations?
			Apparently, I need to use a TChain for this.	*/
			
		//destructor
		virtual ~EECalibration();
		
	private:
		//list of variables
		std::vector<std::string> branchNames = {"PtEle", "etaSCEle", "phiSCEle", "PtEle", "etaSCEle", "phiSCEle"}; 
		//may also need for NT: EIso, HIso, H/EM, sigmaIeIe--depdnds if selection is already made
		struct zVars
		{
			//to eta-sort electrons:
			int tag;
			int proble;
			//arrays for getting values out of NTuple
			float pts[2];//float because that's what Shevin's NTuples use
			float etas[2];
			float phis[2];
			std::vector<int> hix[2]; 
			std::vector<int> hiy[2];
			std::vector<float> hitEnergies[2];
			std::vector<float> rawEnergies[2];
			//calculated values:
			double expectedPt;
			double observedPt;
			double energyRatio;
			double correctedE;
			double correction;
			std::vector<float> hitFractions;
			
			//may also need shower-shape vars
		} Z;
		//safety checks:
		bool histsReady = false;
		bool mapsReady = false;
		bool treeReady = false;
		
		//number of iterations
		const int nIterations; //default should be ~15
		//nominal Z mass
		const double Mz =  91.1876;
		//maps associating a ratio histogram to each Xtal
		xtlHistMap ratioMapP;
		xtlHistMap ratioMapN;
		//hists of final calibration constants
		TH2D *calConstsN;
		TH2D *calConstsP;
		//a Gaussian to fit ratio plots
		TF1 *gfit = new TF1("gfit","gaus",0.2,1.8);
		
		//methods:
		//once per run
		void initHists();//sets all cal consts to 1 with uncerts of 999
		void initMaps();//sets up the maps of ratio hists
		void updateConsts();//updates consts by multiplying with residual mean of E_exp/E_obs--it should call ztlFit() to get those
		void calibrate();//holds actual calibration steps
		void pruneHists();//go through and set any problematic xtls to "-1+/-999"
		void printConsts();//printsconsts to a file in the specified directory
		
		//once per iteration:
		void initTree();//initialize a TTree, needs to take ROOT file name as argument
						//an possibly a vector of needed variable names
		void fillRatHists();//fill the ratio hists
		void resetRatHists(); //reset ratio hists for next iteration
		
		//once per event, or many times per iteration
		void initVars(); //initializes variables to safe values at each new event... probably redundant!
		void get electrons(int evt);//retrieve electrons, eta-ordered, passing selection
		void xtlFit(TH1D* xtlHist);//fits agiven histogram with a gaussian
		
}
		
		
		
		
		
		
		
		
		
		
		
		
