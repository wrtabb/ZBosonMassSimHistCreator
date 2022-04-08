#include <TMath.h>
#include <TBranch.h>
#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TH2F.h>
#include <TH1F.h>
#include <iostream>
#include "file_list.h"

// Functions
bool PassDileptonSelection(double eta1,double eta2,double pt1,double pt2);
bool PassMuonAngle(double pt1,double eta1,double phi1,double mass1,
                   double pt2,double eta2,double phi2,double mass2);
vector<double> GetVariables(double eta1,double eta2,double pt1,double pt2,double phi1,
                            double phi2);
bool GetRecoLeptons(int &idxRecoLead, int &idxRecoSub);
bool GetHardLeptons(int &idxHard1,int &idxHard2);
std::vector<TLorentzVector> GetDressedLeptons(int &idxHardLead,int &idxHardSub);
void Counter(Long64_t event,Long64_t total);
double GetVertexChi2(double pt1,double pt2);

TString base_directory = "/mnt/t2ceph/cms/store/user/wtabb/DrellYan_13TeV_2016/v2p6/DYLL_M50toInf/base/";

TString treeName = "recoTree/DYTree";
const TString muonTrigger1 = "HLT_IsoMu24_v*";
const TString muonTrigger2 = "HLT_IsoTkMu24_v*";
const double etaGapLow = 1.4442;
const double etaGapHigh = 1.566;
const double etaHigh = 2.4;
const double ptLow = 17;
const double ptHigh = 28;
const float dRMinCut = 0.3;
const double ptBinHigh = 499;
const double ptBinLow = 26;
const double etaBinLow = -2.5;
const double etaBinHigh = 2.5;
const double pi = TMath::Pi();
const double xSec = 2009.41; // https://twiki.cern.ch/twiki/bin/view/CMS/SNUCMSYooDYntuple
const double dataLumi = 35867.0; // https://twiki.cern.ch/twiki/bin/view/CMS/SNUCMSYooDYntuple
const int MPSIZE = 5000;
int GENnPair;//number of gen leptons per event
double GENEvt_weight;
double GENLepton_phi[MPSIZE];
double GENLepton_eta[MPSIZE];
double GENLepton_pT[MPSIZE];
double GENLepton_Px[MPSIZE];
double GENLepton_Py[MPSIZE];
double GENLepton_Pz[MPSIZE];
double GENLepton_E[MPSIZE];
int GENLepton_ID[MPSIZE];
int GENLepton_isHardProcess[MPSIZE];
int GENLepton_fromHardProcessFinalState[MPSIZE];
int nGenOthers;
double GenOthers_phi[MPSIZE];
double GenOthers_eta[MPSIZE];
double GenOthers_pT[MPSIZE];
double GenOthers_Px[MPSIZE];
double GenOthers_Py[MPSIZE];
double GenOthers_Pz[MPSIZE];
double GenOthers_E[MPSIZE];
int GenOthers_ID[MPSIZE];
int GenOthers_isHardProcess[MPSIZE];
int GenOthers_isPromptFinalState[MPSIZE];
int nMuon;
int Nmuons;
double Muon_pT[MPSIZE];
double Muon_Px[MPSIZE];
double Muon_Py[MPSIZE];
double Muon_Pz[MPSIZE];
double Muon_eta[MPSIZE];
double Muon_phi[MPSIZE];
double Muon_Inner_pT[MPSIZE];
int Muon_charge[MPSIZE];
double Muon_dxy[MPSIZE];
double Muon_dz[MPSIZE];
bool Muon_passTightID[MPSIZE];
double Muon_PfChargedHadronIsoR04[MPSIZE];
double Muon_PfNeutralHadronIsoR04[MPSIZE];
double Muon_PfGammaIsoR04[MPSIZE];
double Muon_PFSumPUIsoR04[MPSIZE];
double Muon_trkiso[MPSIZE];
double _prefiringweight;
double muMass = 0.1056583715;
int HLT_ntrig;
int HLT_trigType[MPSIZE];
int HLT_trigFired[MPSIZE];
std::vector<std::string> HLT_trigName;
std::vector<std::string> *pHLT_trigName = &HLT_trigName;
int nVertices;
int nPileUp;
std::vector<double> vtxTrkCkt1Pt;
std::vector<double>*pvtxTrkCkt1Pt = &vtxTrkCkt1Pt;

std::vector<double> vtxTrkCkt2Pt;
std::vector<double>*pvtxTrkCkt2Pt = &vtxTrkCkt2Pt;

std::vector<double> vtxTrkChi2;
std::vector<double>*pvtxTrkChi2 = &vtxTrkChi2;

std::vector<double> vtxTrkNdof;
std::vector<double>*pvtxTrkNdof = &vtxTrkNdof;
TBranch*b_runNum;
TBranch*b_evtNum;
TBranch*b_lumiBlock;
TBranch*b_PUweight;
TBranch*b_Nelectrons;
TBranch*b_nVertices;
TBranch*b_nPileUp;

TBranch*b__prefiringweight;
TBranch*b__prefiringweightup;
TBranch*b__prefiringweightdown;
TBranch*b_HLT_ntrig;
TBranch*b_HLT_trigType;
TBranch*b_HLT_trigFired;
TBranch*b_nMuon;
TBranch*b_Nmuons;
TBranch*b_Muon_pT;
TBranch*b_Muon_Px;
TBranch*b_Muon_Py;
TBranch*b_Muon_Pz;
TBranch*b_Muon_eta;
TBranch*b_Muon_phi;
TBranch*b_Muon_Inner_pT;
TBranch*b_Muon_charge;
TBranch*b_Muon_dxy;
TBranch*b_Muon_dz;
TBranch*b_Muon_passTightID;
TBranch*b_Muon_PfChargedHadronIsoR04;
TBranch*b_Muon_PfNeutralHadronIsoR04;
TBranch*b_Muon_PfGammaIsoR04;
TBranch*b_Muon_PFSumPUIsoR04;
TBranch*b_Muon_trkiso;
TBranch*b_GENnPair;
TBranch*b_GENLepton_phi;
TBranch*b_GENLepton_eta;
TBranch*b_GENLepton_pT;
TBranch*b_GENLepton_Px;
TBranch*b_GENLepton_Py;
TBranch*b_GENLepton_Pz;
TBranch*b_GENLepton_E;
TBranch*b_GENLepton_mother;
TBranch*b_GENLepton_mother_pT;
TBranch*b_GENLepton_charge;
TBranch*b_GENLepton_status;
TBranch*b_GENLepton_ID;
TBranch*b_GENLepton_isPrompt;
TBranch*b_GENLepton_isPromptFinalState;
TBranch*b_GENLepton_isTauDecayProduct;
TBranch*b_GENLepton_isPromptTauDecayProduct;
TBranch*b_GENLepton_isDirectPromptTauDecayProductFinalState;
TBranch*b_GENLepton_isHardProcess;
TBranch*b_GENLepton_isLastCopy;
TBranch*b_GENLepton_isLastCopyBeforeFSR;
TBranch*b_GENLepton_isPromptDecayed;
TBranch*b_GENLepton_isDecayedLeptonHadron;
TBranch*b_GENLepton_fromHardProcessBeforeFSR;
TBranch*b_GENLepton_fromHardProcessDecayed;
TBranch*b_GENLepton_fromHardProcessFinalState;
TBranch*b_GENLepton_isMostlyLikePythia6Status3;
TBranch*b_GENEvt_weight;
TBranch*b_GENEvt_QScale;
TBranch*b_GENEvt_x1;
TBranch*b_GENEvt_x2;
TBranch*b_GENEvt_alphaQCD;
TBranch*b_GENEvt_alphaQED;
TBranch*b_nGenOthers;
TBranch*b_GenOthers_phi;
TBranch*b_GenOthers_eta;
TBranch*b_GenOthers_pT;
TBranch*b_GenOthers_Px;
TBranch*b_GenOthers_Py;
TBranch*b_GenOthers_Pz;
TBranch*b_GenOthers_E ;
TBranch*b_GenOthers_ID;
TBranch*b_GenOthers_isHardProcess;
TBranch*b_GenOthers_isPromptFinalState;

void produceHistograms()
{
	TH1::SetDefaultSumw2();

	// Load trees
    int nFiles = file_list.size();
	TChain*chain = new TChain(treeName);
    cout << "Loading trees" << endl;
    cout << "This may take a few minutes" << endl;
    cout << "..." << endl;
    for(int i=0;i<nFiles;i++){
        TString loadFile = base_directory+file_list.at(i);
	    chain->Add(loadFile);
    }

	Long64_t nEntries = chain->GetEntries();
    cout << nEntries << " events loaded" << endl;	
    cout << endl;

	// Define histograms
	TH1D*hInvMassReco    = new TH1D("hInvMassReco","",150,50,200);
	TH1D*hInvMassHard    = new TH1D("hInvMassHard","",150,50,200);
	TH1D*hInvMassDressed = new TH1D("hInvMassDressed","",150,50,200);
    TH2D*hMatrixHard     = new TH2D("hMatrixHard","",150,50,200,150,50,200);
    TH2D*hMatrixDressed  = new TH2D("hMatrixDressed","",150,50,200,150,50,200);

	// Define branches
	chain->SetBranchAddress("nMuon",&nMuon,&b_nMuon);
	chain->SetBranchAddress("Muon_pT",&Muon_pT,&b_Muon_pT);
	chain->SetBranchAddress("Muon_Px",&Muon_Px,&b_Muon_Px);
	chain->SetBranchAddress("Muon_Py",&Muon_Py,&b_Muon_Py);
	chain->SetBranchAddress("Muon_Pz",&Muon_Pz,&b_Muon_Pz);
	chain->SetBranchAddress("Muon_eta",&Muon_eta,&b_Muon_eta);
	chain->SetBranchAddress("Muon_phi",&Muon_phi,&b_Muon_phi);
	chain->SetBranchAddress("Muon_Inner_pT",&Muon_Inner_pT,&b_Muon_Inner_pT);
	chain->SetBranchAddress("Muon_passTightID",&Muon_passTightID,
				&b_Muon_passTightID);
	chain->SetBranchAddress("Muon_charge",&Muon_charge,&b_Muon_charge);
	chain->SetBranchAddress("Muon_PfChargedHadronIsoR04",
				&Muon_PfChargedHadronIsoR04,
				&b_Muon_PfChargedHadronIsoR04);
	chain->SetBranchAddress("Muon_PfNeutralHadronIsoR04",
				&Muon_PfNeutralHadronIsoR04,
				&b_Muon_PfNeutralHadronIsoR04);
	chain->SetBranchAddress("Muon_PfGammaIsoR04",
				&Muon_PfGammaIsoR04,
				&b_Muon_PfGammaIsoR04);
	chain->SetBranchAddress("Muon_PFSumPUIsoR04",
				&Muon_PFSumPUIsoR04,
				&b_Muon_PFSumPUIsoR04);
	chain->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
	chain->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
	chain->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
	chain->SetBranchAddress("HLT_trigName",&pHLT_trigName);
	chain->SetBranchAddress("nVertices",&nVertices,&b_nVertices);
	chain->SetBranchAddress("nPileUp",&nPileUp,&b_nPileUp);
	chain->SetBranchAddress("vtxTrkCkt1Pt",&pvtxTrkCkt1Pt);
	chain->SetBranchAddress("vtxTrkCkt2Pt",&pvtxTrkCkt2Pt);
	chain->SetBranchAddress("vtxTrkChi2",&pvtxTrkChi2);
	chain->SetBranchAddress("vtxTrkNdof",&pvtxTrkNdof);

    chain->SetBranchAddress("_prefiringweight", &_prefiringweight,
                &b__prefiringweight);
    chain->SetBranchAddress("GENnPair", &GENnPair, &b_GENnPair);
    chain->SetBranchAddress("GENLepton_eta", &GENLepton_eta, 
                &b_GENLepton_eta);
    chain->SetBranchAddress("GENLepton_phi",&GENLepton_phi, 
                &b_GENLepton_phi);
    chain->SetBranchAddress("GENLepton_pT",&GENLepton_pT, 
                &b_GENLepton_pT);
    chain->SetBranchAddress("GENLepton_Px",&GENLepton_Px,
                &b_GENLepton_Px);
    chain->SetBranchAddress("GENLepton_Py",&GENLepton_Py,
                &b_GENLepton_Py);
    chain->SetBranchAddress("GENLepton_Pz",&GENLepton_Pz,
                &b_GENLepton_Pz);
    chain->SetBranchAddress("GENLepton_E",&GENLepton_E,
                &b_GENLepton_E);
    chain->SetBranchAddress("GENLepton_ID",&GENLepton_ID, 
                &b_GENLepton_ID);
    chain->SetBranchAddress("GENLepton_isHardProcess",
                &GENLepton_isHardProcess,
                &b_GENLepton_isHardProcess);
    chain->SetBranchAddress("GENLepton_fromHardProcessFinalState",
                &GENLepton_fromHardProcessFinalState,
                &b_GENLepton_fromHardProcessFinalState);
    chain->SetBranchAddress("nGenOthers",&nGenOthers,&b_nGenOthers);
    chain->SetBranchAddress("GenOthers_eta",&GenOthers_eta,
                &b_GenOthers_eta);
    chain->SetBranchAddress("GenOthers_phi",&GenOthers_phi,
                &b_GenOthers_phi);
    chain->SetBranchAddress("GenOthers_pT",&GenOthers_pT,
                &b_GenOthers_pT);
    chain->SetBranchAddress("GenOthers_Px",&GenOthers_Px,
                &b_GenOthers_Px);
    chain->SetBranchAddress("GenOthers_Py",&GenOthers_Py,
                &b_GenOthers_Py);
    chain->SetBranchAddress("GenOthers_Pz",&GenOthers_Pz,
                &b_GenOthers_Pz);
    chain->SetBranchAddress("GenOthers_E",&GenOthers_E,
                &b_GenOthers_E);
    chain->SetBranchAddress("GenOthers_ID",&GenOthers_ID,
                &b_GenOthers_ID);
    chain->SetBranchAddress("GenOthers_isHardProcess",
                &GenOthers_isHardProcess,
                &b_GenOthers_isHardProcess);
    chain->SetBranchAddress("GenOthers_isPromptFinalState",
                &GenOthers_isPromptFinalState,
                &b_GenOthers_isPromptFinalState);
    chain->SetBranchAddress("GENEvt_weight",&GENEvt_weight,
                &b_GENEvt_weight);

	// Find the gen weight sum
//	double sumGenWeight = 5.11429e+07;
	double sumGenWeight = 5.11429e+07;
	double genWeight;
/*
    // Calculate generator weights
    cout << endl;
    cout << "Starting gen weight calculation" << endl;
    Long64_t localEntry;
    for(Long64_t iGen=0;iGen<nEntries;iGen++){
		Counter(iGen,nEntries);
        localEntry = chain->LoadTree(iGen);
        b_GENEvt_weight->GetEntry(localEntry);
        genWeight = GENEvt_weight/fabs(GENEvt_weight);
        sumGenWeight += genWeight;
    }
*/
    cout << "gen weight sum is " << sumGenWeight << endl;
	// Loop over events
	cout << "*********************************" << endl;
    cout << "Starting event loop" << endl;
	for(Long64_t iEntry=0;iEntry<nEntries;iEntry++){
		chain->GetEntry(iEntry);
		Counter(iEntry,nEntries);

		// Check if event passes HLT cut
		TString trigName;
		int trigNameSize = pHLT_trigName->size();
		bool passHLT = false;
		for(int iHLT=0;iHLT<trigNameSize;iHLT++){
			trigName = pHLT_trigName->at(iHLT);
			if(((trigName.CompareTo(muonTrigger1)==0)||
			    (trigName.CompareTo(muonTrigger2)==0)) && 
			     HLT_trigFired[iHLT]==1){
				passHLT = true;
				break;
			} // end if trigName
		}// end loop over triggers

		//-----Get Reconstructed Quantities-----//
		double invMassReco      = -1000;
        double rapidityReco     = -1000;
        double leadPtReco       = -1000;
        double subPtReco        = -1000;

		double ptRecoLead  = -1000;
		double ptRecoSub   = -1000;
		double etaRecoLead = -1000;
		double etaRecoSub  = -1000;
		double phiRecoLead = -1000;
		double phiRecoSub  = -1000;

		int idxRecoLead = -1;
		int idxRecoSub = -1;

		if(passHLT){
			bool recoLep = GetRecoLeptons(idxRecoLead,idxRecoSub);
			if(recoLep){
				ptRecoLead  = Muon_pT[idxRecoLead]; 
				ptRecoSub   = Muon_pT[idxRecoSub];
				etaRecoLead = Muon_eta[idxRecoLead];
				etaRecoSub  = Muon_eta[idxRecoSub];
				phiRecoLead = Muon_phi[idxRecoLead];
				phiRecoSub  = Muon_phi[idxRecoSub];
			}// end if recoLep
		}// end if passHLT

		bool passRecoSelection = PassDileptonSelection(etaRecoLead,etaRecoSub,
							       ptRecoLead,ptRecoSub);
		vector<double> recoVariables;
        recoVariables = GetVariables(etaRecoLead,etaRecoSub,ptRecoLead,ptRecoSub,
                                     phiRecoLead,phiRecoSub);
		if(passRecoSelection){
            invMassReco     = recoVariables.at(0);
            rapidityReco    = recoVariables.at(1);
            leadPtReco      = recoVariables.at(2);
            subPtReco       = recoVariables.at(3);
        }

		//-----Get Hard Process Quantities-----//
		double invMassHard      = -1000;
        double rapidityHard     = -1000;
        double leadPtHard       = -1000;
        double subPtHard        = -1000;

        double ptHardLead  = -1000;
        double ptHardSub   = -1000;
        double etaHardLead = -1000;
        double etaHardSub  = -1000;
        double phiHardLead = -1000;
        double phiHardSub  = -1000;

        int idxHardLead = -1;
        int idxHardSub  = -1;

        bool hardLep = GetHardLeptons(idxHardLead,idxHardSub);
        if(hardLep){
            ptHardLead  = GENLepton_pT[idxHardLead];
            ptHardSub   = GENLepton_pT[idxHardSub];
            etaHardLead = GENLepton_eta[idxHardLead];
            etaHardSub  = GENLepton_eta[idxHardSub];
            phiHardLead = GENLepton_phi[idxHardLead];
            phiHardSub  = GENLepton_phi[idxHardSub];
        }// end if hardLep

        bool passHardSelection = false;
        if(ptHardLead >=0 && ptHardSub >= 0 && etaHardLead > -3 &&
           etaHardSub > -3 && phiHardLead > -100 && phiHardSub > -100 &&
           idxHardLead > -1 && idxHardSub > -1){
                passHardSelection = PassDileptonSelection(etaHardLead,etaHardSub,
                                                          ptHardLead,ptHardSub);
        }

		// Get Hard Variables
		vector<double> hardVariables;
        hardVariables = GetVariables(etaHardLead,etaHardSub,ptHardLead,ptHardSub,
                                     phiHardLead,phiHardSub);

        if(passHardSelection){
            invMassHard     = hardVariables.at(0);
            rapidityHard    = hardVariables.at(1);
            leadPtHard      = hardVariables.at(2);
            subPtHard       = hardVariables.at(3);
        }

		//-----Get Dressed Quantities-----//
		double invMassDressed   = -1000;
        double rapidityDressed  = -1000;
        double leadPtDressed    = -1000;
        double subPtDressed     = -1000;

        int idxDressedLead = -1;
        int idxDressedSub  = -1;

        vector<TLorentzVector> dressedLeptons =
            GetDressedLeptons(idxDressedLead,idxDressedSub);

        double ptDressedLead  = -1000;
        double ptDressedSub   = -1000;
        double etaDressedLead = -1000;
        double etaDressedSub  = -1000;
        double phiDressedLead = -1000;
        double phiDressedSub  = -1000;

        ptDressedLead  = dressedLeptons.at(0).Pt();
        ptDressedSub   = dressedLeptons.at(1).Pt();
        etaDressedLead = dressedLeptons.at(0).Eta();
        etaDressedSub  = dressedLeptons.at(1).Eta();
        phiDressedLead = dressedLeptons.at(0).Phi();
        phiDressedSub  = dressedLeptons.at(1).Phi();

        bool passDressedSelection = PassDileptonSelection(etaDressedLead,
                                                          etaDressedSub,
                                                          ptDressedLead,
                                                          ptDressedSub);

		vector<double> dressedVariables;
                dressedVariables = GetVariables(etaDressedLead,etaDressedSub,ptDressedLead,
                                                ptDressedSub,phiDressedLead,phiDressedSub);


        if(passDressedSelection){
                invMassDressed  = dressedVariables.at(0);
                rapidityDressed = dressedVariables.at(1);
                leadPtDressed   = dressedVariables.at(2);
                subPtDressed    = dressedVariables.at(3);
        }

        // Get gen weight
        double genWeight = (GENEvt_weight/fabs(GENEvt_weight))/sumGenWeight;

		double weight = xSec*dataLumi*genWeight;

		// Fill histograms
		hInvMassReco->Fill(invMassReco,weight);
		hInvMassHard->Fill(invMassHard,weight);
		hInvMassDressed->Fill(invMassDressed,weight);
		hMatrixHard->Fill(invMassHard,invMassReco,weight);
		hMatrixDressed->Fill(invMassDressed,invMassReco,weight);
	}// end loop over entries

	// Save results to output file
	TString saveName = "output_data/unfolding_histograms.root";
	TFile*file = new TFile(saveName,"recreate");
    hInvMassReco->Write();
    hInvMassHard->Write();
    hInvMassDressed->Write();
    hMatrixHard->Write();
    hMatrixDressed->Write();
	file->Close();
}// end analyze

bool PassDileptonSelection(double eta1,double eta2,double pt1,double pt2)
{
        if(abs(eta1)>etaGapLow && abs(eta1)<etaGapHigh) return false;
        if(abs(eta2)>etaGapLow && abs(eta2)<etaGapHigh) return false;
        if(abs(eta1)>etaHigh||abs(eta2)>etaHigh) return false;
        if(pt1>pt2 && (pt1<ptHigh || pt2<ptLow)) return false;
        if(pt2>pt1 && (pt2<ptHigh || pt1<ptLow)) return false;

        return true;
}// end PassDileptonSelection()

vector<double> GetVariables(double eta1,double eta2,double pt1,double pt2,double phi1,
                            double phi2)
{
        TLorentzVector v1;
        TLorentzVector v2;

        v1.SetPtEtaPhiM(pt1,eta1,phi1,muMass);
        v2.SetPtEtaPhiM(pt2,eta2,phi2,muMass);

	// Dimuon invariant mass
        double invMass = (v1+v2).M();

	// Dimuon rapidity
        double rapidity = (v1+v2).Rapidity();

        double ptLead;
        double ptSub;

        if(pt1>pt2){
                ptLead = pt1;
                ptSub  = pt2;
        }
        else{
                ptLead = pt2;
                ptSub  = pt1;
        }

        vector<double> variableReturn = {
                invMass,        // 0
                rapidity,       // 1
                ptLead,         // 2
                ptSub           // 3
        };

        return variableReturn;
}// end GetVariables()

bool GetRecoLeptons(int &idxRecoLead, int &idxRecoSub)
{
	double chargedIso;
	double neutralIso;
	double gammaIso;
	double sumPUPt;
	double pt1,pt2;  
	double eta1,eta2;
	double phi1,phi2;
	double iso_dBeta;
	double charge1,charge2;
	double chi2Old = 1000000000;

	int nDileptons = 0;

	// NOTE: Need to choose two muons by smallest vertex chi2
	// Choosing highest two pT is a temporary placeholder 
	// Add angular cut for muons
	for(int iMu=0;iMu<nMuon;iMu++){
		if(!Muon_passTightID[iMu]) continue;
		chargedIso = Muon_PfChargedHadronIsoR04[iMu];
		neutralIso = Muon_PfNeutralHadronIsoR04[iMu];
		gammaIso = Muon_PfGammaIsoR04[iMu];
		sumPUPt = Muon_PFSumPUIsoR04[iMu];
		pt1 = Muon_pT[iMu];
		eta1 = Muon_eta[iMu];
		phi1 = Muon_phi[iMu];
		iso_dBeta = 
			(chargedIso+max(0.0,neutralIso+gammaIso-0.5*sumPUPt))/pt1;
		if(iso_dBeta > 0.15) continue;
		charge1 = Muon_charge[iMu]; 

		for(int jMu=iMu+1;jMu<nMuon;jMu++){
			if(!Muon_passTightID[jMu]) continue;
			chargedIso = Muon_PfChargedHadronIsoR04[jMu];
			neutralIso = Muon_PfNeutralHadronIsoR04[jMu];
			gammaIso = Muon_PfGammaIsoR04[jMu];
			sumPUPt = Muon_PFSumPUIsoR04[jMu];
			pt2 = Muon_pT[jMu];
			eta2 = Muon_eta[jMu];
			phi2 = Muon_phi[jMu];
			iso_dBeta = 
				(chargedIso+max(0.0,neutralIso+gammaIso-0.5*sumPUPt))/pt2;
			if(iso_dBeta > 0.15) continue;
			charge2 = Muon_charge[jMu]; 
			
			// ensure muons with opposite charge
			if(charge1*charge2 > 0) continue;

			if(!PassMuonAngle(pt1,eta1,phi1,muMass,pt2,eta2,phi2,muMass)) 
				continue;

			double ptIn1 = Muon_Inner_pT[iMu];
			double ptIn2 = Muon_Inner_pT[jMu];
			double chi2 = GetVertexChi2(ptIn1,ptIn2);

			if(chi2>0 && chi2<chi2Old){
				chi2Old = chi2;
				
				if(pt1 > pt2){
					idxRecoLead = iMu;
					idxRecoSub  = jMu;
				}
				else{
					idxRecoLead = jMu;
					idxRecoSub  = iMu;
				} 
				nDileptons++;
			}// end if chi2>0 and chi2 < chi2old
		}//end inner muon loop
	}// end outer muon loop 

	if(nDileptons==1) return true; 
	else return false;
}// end GetRecoLeptons()

bool GetHardLeptons(int &idxHardLead,int &idxHardSub)
{
	int nDileptons = 0;
        for(int iLep=0;iLep<GENnPair;iLep++){
                for(int jLep=iLep+1;jLep<GENnPair;jLep++){
                        if(!(abs(GENLepton_ID[iLep])==13 && abs(GENLepton_ID[jLep])==13))
                                continue;
                        if(GENLepton_ID[iLep]*GENLepton_ID[jLep]>0) continue;
                        if(GENLepton_isHardProcess[iLep]==1 &&
                           GENLepton_isHardProcess[jLep]==1){
                                if(GENLepton_pT[iLep] > GENLepton_pT[jLep]){
                                        idxHardLead = iLep;
                                        idxHardSub = jLep;
                                }// end if iLep is leading electron
                                else{
                                        idxHardLead = jLep;
                                        idxHardSub = iLep;
                                }// end if jLep is leading electron
				nDileptons++;
                        }// end if hard process
                }//end inner loop over gen leptons
        }//end outer loop over gen leptons

        if(nDileptons==1) return true;
        else return false;
}// end GetHardLeptons()

void Counter(Long64_t event,Long64_t total)
{
        int P = 100*(event)/(total);
        if(event%(total/100)==0) 
                cout << P << "%" << endl;
         return;
}

std::vector<TLorentzVector> GetDressedLeptons(int &idxDressedLead,int &idxDressedSub)
{
        TLorentzVector dressed1;
        TLorentzVector dressed2;
	int nDileptons = 0;

	// Loop over muons and select two post-fsr gen-level muons 
        for(int iLep=0;iLep<GENnPair;iLep++){
                for(int jLep=iLep+1;jLep<GENnPair;jLep++){

			// require both leptons be muons
                        if(!(abs(GENLepton_ID[iLep])==13 && abs(GENLepton_ID[jLep])==13))
                                continue;

			// require that they be opposite sign
                        if(GENLepton_ID[iLep]*GENLepton_ID[jLep]>0) continue;
			
			// require that they be post-fsr
                        if(!(GENLepton_fromHardProcessFinalState[iLep]==1 &&
                           GENLepton_fromHardProcessFinalState[jLep]==1)) continue;

			// determine which is lead and which is sub-lead
			if(GENLepton_pT[iLep] > GENLepton_pT[jLep]){
				idxDressedLead = iLep;
				idxDressedSub = jLep;
			}// end if iLep is leading electron
			else{
				idxDressedLead = jLep;
				idxDressedSub = iLep;
			}// end if jLep is leading electron
			nDileptons++;
                }//end inner loop over gen leptons
        }//end outer loop over gen leptons

        double px1 = GENLepton_Px[idxDressedLead];
        double px2 = GENLepton_Px[idxDressedSub];
        double py1 = GENLepton_Py[idxDressedLead];
        double py2 = GENLepton_Py[idxDressedSub];
        double pz1 = GENLepton_Pz[idxDressedLead];
        double pz2 = GENLepton_Pz[idxDressedSub];
        double E1 =  GENLepton_E[idxDressedLead];
        double E2 =  GENLepton_E[idxDressedSub];

        dressed1.SetPxPyPzE(px1,py1,pz1,E1);
        dressed2.SetPxPyPzE(px2,py2,pz2,E2);

        double eta1 = GENLepton_eta[idxDressedLead];
        double eta2 = GENLepton_eta[idxDressedSub];
        double phi1 = GENLepton_phi[idxDressedLead];
        double phi2 = GENLepton_phi[idxDressedSub];

        double dRMin = 0.1;
        double etaPho,phiPho;
        double etaDiff1,phiDiff1;

        double etaDiff2,phiDiff2;
        double dR1Squared,dR1;
        double dR2Squared,dR2;
        double pxPho,pyPho,pzPho,EPho;

	// Loop over photons
	if(nDileptons==1){
		for(int iPho=0;iPho<nGenOthers;iPho++){
			// require that they be photons and are prompt final state
			if(abs(GenOthers_ID[iPho])!=22 ||
			   GenOthers_isPromptFinalState[iPho]!=1) continue;

			// define location of photon in the eta-phi plane
			etaPho = GenOthers_eta[iPho];
			phiPho = GenOthers_phi[iPho];

			// find distance, dR1, of photon from muon1
			etaDiff1 = eta1-etaPho;
			phiDiff1 = phi1-phiPho;
			dR1Squared = etaDiff1*etaDiff1+phiDiff1*phiDiff1;
			dR1 = sqrt(dR1Squared);

			// find distance, dR2, of photon from muon2
			etaDiff2 = eta2-etaPho;
			phiDiff2 = phi2-phiPho;
			dR2Squared = etaDiff2*etaDiff2+phiDiff2*phiDiff2;
			dR2 = sqrt(dR2Squared);

			// create lorentz vector for photon
			pxPho = GenOthers_Px[iPho];
			pyPho = GenOthers_Py[iPho];
			pzPho = GenOthers_Pz[iPho];
			EPho  = GenOthers_E[iPho];
			TLorentzVector phoVec;
			phoVec.SetPxPyPzE(pxPho,pyPho,pzPho,EPho);

			// only keep photons which are closer than 0.1 from one of the muons
			if(dR1>dRMin && dR2>dRMin) continue;

			// add its four momentum to the four momentum of the muon
			// it is closest to
			if(dR1<dR2){
				dressed1 += phoVec;
			}
			else{
				dressed2 += phoVec;
			}
		}// end loop over photons       
	}// end if nDileptons==1

	// after looping over each photon, the dressed vectors should contain the 
	// original post-fsr muons with the associated photons added to them
	// place these two lorentz vectors into a vector and return
        vector<TLorentzVector> returnVector;
        returnVector.push_back(dressed1);
        returnVector.push_back(dressed2);

        return returnVector;
}// end GetDressedLeptons()

bool PassMuonAngle(double pt1,double eta1,double phi1,double mass1,
                   double pt2,double eta2,double phi2,double mass2)
{
	double angle;
	double limit = pi - 0.005;
	TLorentzVector v1;
	TLorentzVector v2;
	
	v1.SetPtEtaPhiM(pt1,eta1,phi1,mass1);
	v2.SetPtEtaPhiM(pt2,eta2,phi2,mass2);
	
	angle = v1.Angle(v2.Vect());

	if(angle>limit) return false;
	else return true;
}

double GetVertexChi2(double pt1,double pt2)
{
	int nPt1 = vtxTrkCkt1Pt.size();
	int nPt2 = vtxTrkCkt2Pt.size();
	double chi2_dof = -1000;
	if(nPt1 != nPt2) cout << "nPt1 = " << nPt1 << ", nPt2 = " << nPt2 << endl;

	for(int i=0;i<nPt1;i++){
		if( (vtxTrkCkt1Pt.at(i) == pt1 && vtxTrkCkt2Pt.at(i) == pt2) ||
		    (vtxTrkCkt1Pt.at(i) == pt2 && vtxTrkCkt2Pt.at(i) == pt1))
			chi2_dof = vtxTrkChi2[i]/vtxTrkNdof[i];
	}
	return chi2_dof;	
}