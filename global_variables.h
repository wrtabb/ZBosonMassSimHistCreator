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
void SetBranches(TChain*chain);
double GetGenWeights(Long64_t nEntries,TChain*chain);

// Where data is saves
TString base_directory = "/mnt/t2ceph/cms/store/user/wtabb/DrellYan_13TeV_2016/v2p6/DYLL_M50toInf/base/";

// Global parameters
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

// Branch variables
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

// Branches
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

