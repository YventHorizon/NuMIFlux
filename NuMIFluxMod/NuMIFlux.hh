#ifndef NuMIFlux_h
#define NuMIFlux_h

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TSystem.h"
#include "TVector3.h"


#include "FluggNtuple/FluxNtuple.h"
//#include "calcLocationWeights.h"

class NuMIFlux {
public :

  int Nfiles = 0;

  static const int numu  = 56;
  static const int anumu = 55;
  static const int nue   = 53;
  static const int anue  = 52;
    
    static const int decaymode_one = 1;
    static const int decaymode_two = 2;
    static const int decaymode_three = 3;
    static const int decaymode_four = 4;
    static const int decaymode_five = 5;
    static const int decaymode_six = 6;
    static const int decaymode_seven = 7;
    static const int decaymode_eight = 8;
    static const int decaymode_nine = 9;
    static const int decaymode_ten = 10;
    static const int decaymode_eleven = 11;
    static const int decaymode_twelve = 12;
    static const int decaymode_thirteen = 13;
    static const int decaymode_fourteen = 14;
    
//   static const int kpdg_pionplus   =   8;  // Geant  8

  int highest_evtno = 0;
  double NominalPOT = 6e20;
  bool debug = false;
  double fDefaultWeightCorrection = 1./(10000. * TMath::Pi());
  double Ntarget = 4.76e31/56.41e6*256.35*233*1036.8; //TPC active!!!
  double AccumulatedPOT=0.;
  int treeNumber = -1;

  double histMin = 0;
  double histMax = 20;
  int histNbins = 400;

  TChain *cflux;

  FluxNtuple *fluxNtuple;

  TH1D* numuFluxHisto;
  TH1D* anumuFluxHisto;
  TH1D* nueFluxHisto;
  TH1D* anueFluxHisto;
    
    TH1D* nueFluxHisto_orig_MIPP_con;
    TH1D* nueFluxHisto_orig_NA49_con;
    TH1D* nueFluxHisto_orig_NA49_proton_con;
    TH1D* nueFluxHisto_orig_NA49_pion_con;
    TH1D* nueFluxHisto_orig_NA49_kaon_con;
    TH1D* nueFluxHisto_orig_NA49_neutron_con;
    TH1D* nueFluxHisto_orig;
    
    TH1D* numuFluxHisto_up;
    TH1D* anumuFluxHisto_up;
    TH1D* nueFluxHisto_up;
    TH1D* nueFluxHisto_tgp;
    TH1D* nueFluxHisto_up_tgp;
    TH1D* nueFluxHisto_down_tgp;
    TH1D* anueFluxHisto_up;
    
    TH1D* numuFluxHisto_down;
    TH1D* anumuFluxHisto_down;
    TH1D* nueFluxHisto_down;
    TH1D* anueFluxHisto_down;

  TH1D* nueFluxHisto_noDAR; //added new hists
  TH1D* anueFluxHisto_noDAR; //added new hists
  TH1D* nueCCHisto;
    TH1D* nueCCHisto_tgp;
    TH1D* nueCCHisto_down;
    TH1D* nueCCHisto_up;
    TH1D* nueCCHisto_down_tgp;
    TH1D* nueCCHisto_up_tgp;
  TGraph *genieXsecNueCC;
    
  TH1D* pionplusAngleHisto;//added new hists
  TH1D* pionplusMomentumHisto;//added new hists
  TH2D* pionplus2DHisto;//added new hists
    TH2D* pionplus_nue;
    TH2D* pionplus_numu;
    TH2D* pionplus_anue;
    TH2D* pionplus_anumu;
    
    
TH1D* pionminusAngleHisto;//added new hists
TH1D* pionminusMomentumHisto;//added new hists
TH2D* pionminus2DHisto;//added new hists
    TH2D* pionminus_nue;
    TH2D* pionminus_numu;
    TH2D* pionminus_anue;
    TH2D* pionminus_anumu;

TH1D* KlongAngleHisto;//added new hists
TH1D* KlongMomentumHisto;//added new hists
TH2D* Klong2DHisto;//added new hists
    TH2D* Klong_nue;
    TH2D* Klong_numu;
    TH2D* Klong_anue;
    TH2D* Klong_anumu;

TH1D* KminusAngleHisto;//added new hists
TH1D* KminusMomentumHisto;//added new hists
TH2D* Kminus2DHisto;//added new hists
    TH2D* Kminus_nue;
    TH2D* Kminus_numu;
    TH2D* Kminus_anue;
    TH2D* Kminus_anumu;

TH1D* KplusAngleHisto;//added new hists
TH1D* KplusMomentumHisto;//added new hists
TH2D* Kplus2DHisto;//added new hists
    TH2D* Kplus_nue;
    TH2D* Kplus_numu;
    TH2D* Kplus_anue;
    TH2D* Kplus_anumu;
    
    TH1D* muonplusAngleHisto;//added new hists
    TH1D* muonplusMomentumHisto;//added new hists
    TH2D* muonplus2DHisto;//added new hists
    
    TH1D* muonminusAngleHisto;//added new hists
    TH1D* muonminusMomentumHisto;//added new hists
    TH2D* muonminus2DHisto;//added new hists
    
    TH1D* protonAngleHisto;//added new hists
    TH1D* protonMomentumHisto;//added new hists
    TH2D* proton2DHisto;//added new hists
    TH2D* proton_nue;
    TH2D* proton_numu;
    TH2D* proton_anue;
    TH2D* proton_anumu;
    
    TH1D* neutronAngleHisto;//added new hists
    TH1D* neutronMomentumHisto;//added new hists
    TH2D* neutron2DHisto;//added new hists
    TH2D* neutron_nue;
    TH2D* neutron_numu;
    TH2D* neutron_anue;
    TH2D* neutron_anumu;
    
    TH1D* lambdaAngleHisto;//added new hists
    TH1D* lambdaMomentumHisto;//added new hists
    TH2D* lambda2DHisto;//added new hists
    
    TH1D* sigmaplusAngleHisto;//added new hists
    TH1D* sigmaplusMomentumHisto;//added new hists
    TH2D* sigmaplus2DHisto;//added new hists
    
    
    ////////////// nue_parent //////////////
    TH1D* nue_parent_muonplus;
    TH1D* nue_parent_kaonplus;
    TH1D* nue_parent_kaonlong;
    
    TH1D* nue_tptype_muonplus;
    TH1D* nue_tptype_kaonplus;
    TH1D* nue_tptype_pionplus;
    TH1D* nue_tptype_kaonlong;
    
    TH1D* nue_tgptype_muonplus_kaonplus;
    TH1D* nue_tgptype_muonplus_kaonlong;
    
    ////////////// numu_parent //////////////
    TH1D* numu_parent_pionplus;
    TH1D* numu_parent_kaonplus;
    TH1D* numu_parent_kaonlong;
    TH1D* numu_parent_muonplus;
    TH1D* numu_parent_muonminus;
    
    ////////////// anue_parent //////////////
    TH1D* anue_parent_kaonplus;
    
    ////////////// anumu_parent //////////////
    TH1D* anumu_parent_kaonplus;
    
    /////////// nu Ndecay ////////////
    TH1D* nue_Ndecay_one;
    TH1D* nue_Ndecay_six;
    TH1D* nue_Ndecay_eleven;
    TH1D* nue_Ndecay_five;
    
    TH1D* numu_Ndecay_thirteen;
    TH1D* numu_Ndecay_seven;
    TH1D* numu_Ndecay_three;
    TH1D* numu_Ndecay_five;
    TH1D* numu_Ndecay_twelve;
    
    ////// nu tgptype //////////
    TH1D* nue_tgptype_pionplus;
    TH1D* nue_tgptype_pionminus;
    TH1D* nue_tgptype_kaonlong;
    TH1D* nue_tgptype_kaonplus;
    TH1D* nue_tgptype_kaonminus;
    TH1D* nue_tgptype_neutron;
    TH1D* nue_tgptype_proton;
    TH1D* nue_tgptype_lambda;
    TH1D* nue_tgptype_sigmaplus;
    
    TH1D* numu_tgptype_pionplus;
    TH1D* numu_tgptype_pionminus;
    TH1D* numu_tgptype_kaonlong;
    TH1D* numu_tgptype_kaonplus;
    TH1D* numu_tgptype_kaonminus;
    TH1D* numu_tgptype_neutron;
    TH1D* numu_tgptype_proton;
    TH1D* numu_tgptype_lambda;
    TH1D* numu_tgptype_sigmaplus;
    
    //////////// nu tgen ////////////
    TH1D* nue_tgen_two;
    TH1D* nue_tgen_three;
    TH1D* nue_tgen_four;
    TH1D* nue_tgen_five;
    TH1D* nue_tgen_six;
    TH1D* nue_tgen_seven;
    TH1D* nue_tgen_eight;
    TH1D* nue_tgen_nine;
    TH1D* nue_tgen_ten;
    TH1D* nue_tgen_eleven;
    TH1D* nue_tgen_twelve;
    
    TH1D* numu_tgen_two;
    TH1D* numu_tgen_three;
    TH1D* numu_tgen_four;
    TH1D* numu_tgen_five;
    TH1D* numu_tgen_six;
    TH1D* numu_tgen_seven;
    TH1D* numu_tgen_eight;
    TH1D* numu_tgen_nine;
    TH1D* numu_tgen_ten;
    TH1D* numu_tgen_eleven;
    TH1D* numu_tgen_twelve;
    
    
    ////////////// decay mode hists /////////////////////////////////////
    // nue
    TH2D* ppnue_dm_one;
    TH2D* ppnue_dm_two;
    TH2D* ppnue_dm_three;
    TH2D* ppnue_dm_four;
    TH2D* ppnue_dm_five;
    TH2D* ppnue_dm_six;
    TH2D* ppnue_dm_seven;
    TH2D* ppnue_dm_eight;
    TH2D* ppnue_dm_nine;
    TH2D* ppnue_dm_ten;
    TH2D* ppnue_dm_eleven;
    TH2D* ppnue_dm_twelve;
    TH2D* ppnue_dm_thirteen;
    TH2D* ppnue_dm_fourteen;
    
    TH2D* pmnue_dm_one;
    TH2D* pmnue_dm_two;
    TH2D* pmnue_dm_three;
    TH2D* pmnue_dm_four;
    TH2D* pmnue_dm_five;
    TH2D* pmnue_dm_six;
    TH2D* pmnue_dm_seven;
    TH2D* pmnue_dm_eight;
    TH2D* pmnue_dm_nine;
    TH2D* pmnue_dm_ten;
    TH2D* pmnue_dm_eleven;
    TH2D* pmnue_dm_twelve;
    TH2D* pmnue_dm_thirteen;
    TH2D* pmnue_dm_fourteen;
    
    TH2D* klnue_dm_one;
    TH2D* klnue_dm_two;
    TH2D* klnue_dm_three;
    TH2D* klnue_dm_four;
    TH2D* klnue_dm_five;
    TH2D* klnue_dm_six;
    TH2D* klnue_dm_seven;
    TH2D* klnue_dm_eight;
    TH2D* klnue_dm_nine;
    TH2D* klnue_dm_ten;
    TH2D* klnue_dm_eleven;
    TH2D* klnue_dm_twelve;
    TH2D* klnue_dm_thirteen;
    TH2D* klnue_dm_fourteen;
    
    TH2D* kpnue_dm_one;
    TH2D* kpnue_dm_two;
    TH2D* kpnue_dm_three;
    TH2D* kpnue_dm_four;
    TH2D* kpnue_dm_five;
    TH2D* kpnue_dm_six;
    TH2D* kpnue_dm_seven;
    TH2D* kpnue_dm_eight;
    TH2D* kpnue_dm_nine;
    TH2D* kpnue_dm_ten;
    TH2D* kpnue_dm_eleven;
    TH2D* kpnue_dm_twelve;
    TH2D* kpnue_dm_thirteen;
    TH2D* kpnue_dm_fourteen;
    
    TH2D* kmnue_dm_one;
    TH2D* kmnue_dm_two;
    TH2D* kmnue_dm_three;
    TH2D* kmnue_dm_four;
    TH2D* kmnue_dm_five;
    TH2D* kmnue_dm_six;
    TH2D* kmnue_dm_seven;
    TH2D* kmnue_dm_eight;
    TH2D* kmnue_dm_nine;
    TH2D* kmnue_dm_ten;
    TH2D* kmnue_dm_eleven;
    TH2D* kmnue_dm_twelve;
    TH2D* kmnue_dm_thirteen;
    TH2D* kmnue_dm_fourteen;
    
    TH2D* pnue_dm_one;
    TH2D* pnue_dm_two;
    TH2D* pnue_dm_three;
    TH2D* pnue_dm_four;
    TH2D* pnue_dm_five;
    TH2D* pnue_dm_six;
    TH2D* pnue_dm_seven;
    TH2D* pnue_dm_eight;
    TH2D* pnue_dm_nine;
    TH2D* pnue_dm_ten;
    TH2D* pnue_dm_eleven;
    TH2D* pnue_dm_twelve;
    TH2D* pnue_dm_thirteen;
    TH2D* pnue_dm_fourteen;
    
    TH2D* nnue_dm_one;
    TH2D* nnue_dm_two;
    TH2D* nnue_dm_three;
    TH2D* nnue_dm_four;
    TH2D* nnue_dm_five;
    TH2D* nnue_dm_six;
    TH2D* nnue_dm_seven;
    TH2D* nnue_dm_eight;
    TH2D* nnue_dm_nine;
    TH2D* nnue_dm_ten;
    TH2D* nnue_dm_eleven;
    TH2D* nnue_dm_twelve;
    TH2D* nnue_dm_thirteen;
    TH2D* nnue_dm_fourteen;
    
    // numu
    TH2D* ppnumu_dm_one;
    TH2D* ppnumu_dm_two;
    TH2D* ppnumu_dm_three;
    TH2D* ppnumu_dm_four;
    TH2D* ppnumu_dm_five;
    TH2D* ppnumu_dm_six;
    TH2D* ppnumu_dm_seven;
    TH2D* ppnumu_dm_eight;
    TH2D* ppnumu_dm_nine;
    TH2D* ppnumu_dm_ten;
    TH2D* ppnumu_dm_eleven;
    TH2D* ppnumu_dm_twelve;
    TH2D* ppnumu_dm_thirteen;
    TH2D* ppnumu_dm_fourteen;
    
    TH2D* pmnumu_dm_one;
    TH2D* pmnumu_dm_two;
    TH2D* pmnumu_dm_three;
    TH2D* pmnumu_dm_four;
    TH2D* pmnumu_dm_five;
    TH2D* pmnumu_dm_six;
    TH2D* pmnumu_dm_seven;
    TH2D* pmnumu_dm_eight;
    TH2D* pmnumu_dm_nine;
    TH2D* pmnumu_dm_ten;
    TH2D* pmnumu_dm_eleven;
    TH2D* pmnumu_dm_twelve;
    TH2D* pmnumu_dm_thirteen;
    TH2D* pmnumu_dm_fourteen;
    
    TH2D* klnumu_dm_one;
    TH2D* klnumu_dm_two;
    TH2D* klnumu_dm_three;
    TH2D* klnumu_dm_four;
    TH2D* klnumu_dm_five;
    TH2D* klnumu_dm_six;
    TH2D* klnumu_dm_seven;
    TH2D* klnumu_dm_eight;
    TH2D* klnumu_dm_nine;
    TH2D* klnumu_dm_ten;
    TH2D* klnumu_dm_eleven;
    TH2D* klnumu_dm_twelve;
    TH2D* klnumu_dm_thirteen;
    TH2D* klnumu_dm_fourteen;
    
    TH2D* kpnumu_dm_one;
    TH2D* kpnumu_dm_two;
    TH2D* kpnumu_dm_three;
    TH2D* kpnumu_dm_four;
    TH2D* kpnumu_dm_five;
    TH2D* kpnumu_dm_six;
    TH2D* kpnumu_dm_seven;
    TH2D* kpnumu_dm_eight;
    TH2D* kpnumu_dm_nine;
    TH2D* kpnumu_dm_ten;
    TH2D* kpnumu_dm_eleven;
    TH2D* kpnumu_dm_twelve;
    TH2D* kpnumu_dm_thirteen;
    TH2D* kpnumu_dm_fourteen;
    
    TH2D* kmnumu_dm_one;
    TH2D* kmnumu_dm_two;
    TH2D* kmnumu_dm_three;
    TH2D* kmnumu_dm_four;
    TH2D* kmnumu_dm_five;
    TH2D* kmnumu_dm_six;
    TH2D* kmnumu_dm_seven;
    TH2D* kmnumu_dm_eight;
    TH2D* kmnumu_dm_nine;
    TH2D* kmnumu_dm_ten;
    TH2D* kmnumu_dm_eleven;
    TH2D* kmnumu_dm_twelve;
    TH2D* kmnumu_dm_thirteen;
    TH2D* kmnumu_dm_fourteen;
    
    TH2D* pnumu_dm_one;
    TH2D* pnumu_dm_two;
    TH2D* pnumu_dm_three;
    TH2D* pnumu_dm_four;
    TH2D* pnumu_dm_five;
    TH2D* pnumu_dm_six;
    TH2D* pnumu_dm_seven;
    TH2D* pnumu_dm_eight;
    TH2D* pnumu_dm_nine;
    TH2D* pnumu_dm_ten;
    TH2D* pnumu_dm_eleven;
    TH2D* pnumu_dm_twelve;
    TH2D* pnumu_dm_thirteen;
    TH2D* pnumu_dm_fourteen;
    
    TH2D* nnumu_dm_one;
    TH2D* nnumu_dm_two;
    TH2D* nnumu_dm_three;
    TH2D* nnumu_dm_four;
    TH2D* nnumu_dm_five;
    TH2D* nnumu_dm_six;
    TH2D* nnumu_dm_seven;
    TH2D* nnumu_dm_eight;
    TH2D* nnumu_dm_nine;
    TH2D* nnumu_dm_ten;
    TH2D* nnumu_dm_eleven;
    TH2D* nnumu_dm_twelve;
    TH2D* nnumu_dm_thirteen;
    TH2D* nnumu_dm_fourteen;
    
    
    ////////////// MIPP hists /////////////////////////////////////
    
    TH2D* pionplus_MIPP;
    TH2D* pionplus_MIPP_nue;
    TH2D* pionplus_MIPP_anue;
    TH2D* pionplus_MIPP_numu;
    TH2D* pionplus_MIPP_anumu;
    
    TH2D* pionminus_MIPP;
    TH2D* pionminus_MIPP_nue;
    TH2D* pionminus_MIPP_anue;
    TH2D* pionminus_MIPP_numu;
    TH2D* pionminus_MIPP_anumu;
    
    TH2D* Kplus_MIPP;
    TH2D* Kplus_MIPP_nue;
    TH2D* Kplus_MIPP_anue;
    TH2D* Kplus_MIPP_numu;
    TH2D* Kplus_MIPP_anumu;
    
    TH2D* Kminus_MIPP;
    TH2D* Kminus_MIPP_nue;
    TH2D* Kminus_MIPP_nue_dm_one;
    TH2D* Kminus_MIPP_nue_dm_eleven;
    TH2D* Kminus_MIPP_nue_dm_six;
    
    TH2D* Kminus_MIPP_anue;
    TH2D* Kminus_MIPP_numu;
    TH2D* Kminus_MIPP_anumu;
    
    TH2D* pionplus_NA49;
    TH2D* pionplus_NA49_nue;
    TH2D* pionplus_MIPP_nue_dm1;
    TH2D* pionplus_MIPP_nue_dm6;
    TH2D* pionplus_MIPP_nue_dm5;
    TH2D* pionplus_MIPP_nue_dm11;
    TH2D* pionminus_NA49;
    TH2D* pionminus_NA49_nue;
    TH2D* Kplus_NA49;
    TH2D* Kplus_NA49_second;
    TH2D* Kplus_NA49_nue;
    TH2D* Kminus_NA49;
    TH2D* Kminus_NA49_nue;
    
    TH2D* proton_NA49;
    TH2D* proton_NA49_nue;
    
    TH2D* neutron_NA49;
    TH2D* neutron_NA49_nue;
    
    TH2D* antiproton_NA49;
    TH2D* antiproton_NA49_nue;
    
    
  TFile* f = new TFile("NuMIFlux.root", "RECREATE");


  NuMIFlux(string pattern="/uboone/data/flux/numi/v2/flugg_mn000z200i_rp11_lowth_pnut_f112c0f093bbird/flugg_mn000z200i_rp11_bs1.1_pnut_lowth_f112c0f093bbird_0*.root");
//"/uboone/data/flux/numi/current/flugg_mn000z200i_20101117.gpcfgrid_lowth/flugg_mn000z200i_20101117.gpcfgrid_lowth_00*.root"
//
  virtual ~NuMIFlux();

  void CalculateFlux();
  TVector3 RandomInTPC();
  TVector3 FromDetToBeam(const TVector3& det);
  double estimate_pots(int highest_potnum);
  int calcEnuWgt( FluxNtuple* decay, const TVector3& xyz, double& enu, double& wgt_xy);
    
};

#endif

#ifdef NuMIFlux_cxx

NuMIFlux::NuMIFlux(string pattern) {

  const char* path = "/uboone/app/users/mdeltutt/NuMIFlux";
  if ( path ) {
    TString libs = gSystem->GetDynamicPath();
    libs += ":";
    libs += path;
    //libs += "/lib";
    gSystem->SetDynamicPath(libs.Data());       
    gSystem->Load("FluxNtuple_C.so");
  }

  cflux = new TChain("h10");
  cflux->Add(pattern.c_str());

  Nfiles = cflux->GetNtrees();
  cout << "Number of files: " << Nfiles << endl;

  //Inizialise histos
  TString titleBase1 = "Neutrino Flux;";
  TString titleBase2 = " Energy [GeV];";
  TString titleBase3 = " / cm^{2} / 6e20 POT";
  // numu
  numuFluxHisto = new TH1D("numuFluxHisto", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
  // anumu
  anumuFluxHisto = new TH1D("anumuFluxHisto", (titleBase1 + "#bar{#nu}_{#mu}" + titleBase2 +"#bar{#nu}_{#mu}" + titleBase3),histNbins,histMin,histMax);
  // nue
  nueFluxHisto = new TH1D("nueFluxHisto", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    
     nueFluxHisto_tgp = new TH1D("nueFluxHisto_tgp", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    
    nueFluxHisto_orig = new TH1D("nueFluxHisto_orig", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    
    
    nueFluxHisto_orig_MIPP_con = new TH1D("nueFluxHisto_orig_MIPP_con", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    
    nueFluxHisto_orig_NA49_con = new TH1D("nueFluxHisto_orig_NA49_con", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    
       nueFluxHisto_orig_NA49_proton_con = new TH1D("nueFluxHisto_orig_NA49_proton_con", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    
       nueFluxHisto_orig_NA49_kaon_con = new TH1D("nueFluxHisto_orig_NA49_kaon_con", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    
       nueFluxHisto_orig_NA49_pion_con = new TH1D("nueFluxHisto_orig_NA49_pion_con", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    
    nueFluxHisto_orig_NA49_neutron_con = new TH1D("nueFluxHisto_orig_NA49_neutron_con", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);

    anueFluxHisto = new TH1D("anueFluxHisto", (titleBase1 + "#bar{#nu}_{e}" + titleBase2 + "#bar{#nu}_{e}" + titleBase3),histNbins,histMin,histMax);

    // weight reweighted up histos
  anueFluxHisto_up = new TH1D("anueFluxHisto_up", (titleBase1 + "#bar{#nu}_{e}" + titleBase2 + "#bar{#nu}_{e}" + titleBase3),histNbins,histMin,histMax);
    
    numuFluxHisto_up = new TH1D("numuFluxHisto_up", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    
    anumuFluxHisto_up = new TH1D("anumuFluxHisto_up", (titleBase1 + "#bar{#nu}_{#mu}" + titleBase2 +"#bar{#nu}_{#mu}" + titleBase3),histNbins,histMin,histMax);

    nueFluxHisto_up = new TH1D("nueFluxHisto_up", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    
    nueFluxHisto_up_tgp = new TH1D("nueFluxHisto_up_tgp", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3 + "up"),histNbins,histMin,histMax);
    
    
    // weight reweighted down histos
    anueFluxHisto_down = new TH1D("anueFluxHisto_down", (titleBase1 + "#bar{#nu}_{e}" + titleBase2 + "#bar{#nu}_{e}" + titleBase3),histNbins,histMin,histMax);
    
    numuFluxHisto_down = new TH1D("numuFluxHisto_down", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    
    anumuFluxHisto_down = new TH1D("anumuFluxHisto_down", (titleBase1 + "#bar{#nu}_{#mu}" + titleBase2 +"#bar{#nu}_{#mu}" + titleBase3),histNbins,histMin,histMax);
    
    nueFluxHisto_down = new TH1D("nueFluxHisto_down", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nueFluxHisto_down_tgp = new TH1D("nueFluxHisto_down_tgp", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3 + "down"),histNbins,histMin,histMax);
 
    // no DAR histos
    nueFluxHisto_noDAR = new TH1D("nueFluxHisto_noDAR", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);

  anueFluxHisto_noDAR = new TH1D("anueFluxHisto_noDAR", (titleBase1 + "#bar{#nu}_{e}" + titleBase2 + "#bar{#nu}_{e}" + titleBase3),histNbins,histMin,histMax);
    
  nueCCHisto = new TH1D("nueCCHisto", "nue CC; #nu_{e} Energy [GeV]; #nu_{e} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
    
    nueCCHisto_tgp = new TH1D("nueCCHisto_tgp", "nue CC; #nu_{e} Energy [GeV]; #nu_{e} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
    
    nueCCHisto_up = new TH1D("nueCCHisto_up", "nue CC; #nu_{e} Energy [GeV]; #nu_{e} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
    
    nueCCHisto_down = new TH1D("nueCCHisto_down", "nue CC; #nu_{e} Energy [GeV]; #nu_{e} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
    
    
    nueCCHisto_up_tgp = new TH1D("nueCCHisto_up_tgp", "nue CC; #nu_{e} Energy [GeV]; #nu_{e} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
    
    nueCCHisto_down_tgp = new TH1D("nueCCHisto_down_tgp", "nue CC; #nu_{e} Energy [GeV]; #nu_{e} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
    
    //////////////////// nue_parent /////////////////////
    nue_parent_muonplus = new TH1D("nue_parent_muonplus", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_parent_kaonlong = new TH1D("nue_parent_kaonlong", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_parent_kaonplus = new TH1D("nue_parent_kaonplus", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    
    nue_tptype_muonplus = new TH1D("nue_tptype_muonplus", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tptype_kaonlong = new TH1D("nue_tptype_kaonlong", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tptype_kaonplus = new TH1D("nue_tptype_kaonplus", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tptype_pionplus = new TH1D("nue_tptype_pionplus", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    
    nue_tgptype_muonplus_kaonplus = new TH1D("nue_tgptype_muonplus_kaonplus", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgptype_muonplus_kaonlong = new TH1D("nue_tgptype_muonplus_kaonlong", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    
    //////////////////// numu_parent /////////////////////
    numu_parent_pionplus = new TH1D("numu_parent_pionplus", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_parent_kaonlong = new TH1D("numu_parent_kaonlong", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_parent_kaonplus = new TH1D("numu_parent_kaonplus", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_parent_muonplus = new TH1D("numu_parent_muonplus", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_parent_muonminus = new TH1D("numu_parent_muonminus", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    
    
    anue_parent_kaonplus = new TH1D("anue_parent_kaonplus", (titleBase1 + "#bar{#nu}_{e}" + titleBase2 +"#bar{#nu}_{e}" + titleBase3),histNbins,histMin,histMax);
    
    anumu_parent_kaonplus = new TH1D("anumu_parent_kaonplus", (titleBase1 + "#bar{#nu}_{#mu}" + titleBase2 +"#bar{#nu}_{#mu}" + titleBase3),histNbins,histMin,histMax);
   
    /////////////// nu Ndecay ////////////////////
    nue_Ndecay_one = new TH1D("nue_Ndecay_one", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_Ndecay_six = new TH1D("nue_Ndecay_six", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_Ndecay_eleven = new TH1D("nue_Ndecay_eleven", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_Ndecay_five = new TH1D("nue_Ndecay_five", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    
    numu_Ndecay_thirteen = new TH1D("numu_Ndecay_thirteen", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_Ndecay_seven = new TH1D("numu_Ndecay_seven", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_Ndecay_three = new TH1D("numu_Ndecay_three", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_Ndecay_five = new TH1D("numu_Ndecay_five", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_Ndecay_twelve = new TH1D("numu_Ndecay_twelve", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    
    ///////////// nu tgptype /////////////////
    nue_tgptype_pionplus = new TH1D("nue_tgptype_pionplus", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgptype_pionminus = new TH1D("nue_tgptype_pionminus", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgptype_kaonlong = new TH1D("nue_tgptype_kaonlong", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgptype_kaonplus = new TH1D("nue_tgptype_kaonplus", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgptype_kaonminus = new TH1D("nue_tgptype_kaonminus", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgptype_neutron = new TH1D("nue_tgptype_neutron", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgptype_proton = new TH1D("nue_tgptype_proton", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgptype_lambda = new TH1D("nue_tgptype_lambda", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgptype_sigmaplus = new TH1D("nue_tgptype_sigmaplus", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    
    numu_tgptype_pionplus = new TH1D("numu_tgptype_pionplus", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgptype_pionminus = new TH1D("numu_tgptype_pionminus", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgptype_kaonlong = new TH1D("numu_tgptype_kaonlong", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgptype_kaonplus = new TH1D("numu_tgptype_kaonplus", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgptype_kaonminus = new TH1D("numu_tgptype_kaonminus", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgptype_neutron = new TH1D("numu_tgptype_neutron", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgptype_proton = new TH1D("numu_tgptype_proton", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgptype_lambda = new TH1D("numu_tgptype_lambda", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgptype_sigmaplus = new TH1D("numu_tgptype_sigmaplus", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    
    ///////////// nu tgen /////////////////
    
    nue_tgen_two = new TH1D("nue_tgen_two", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgen_three = new TH1D("nue_tgen_three", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgen_four = new TH1D("nue_tgen_four", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgen_five = new TH1D("nue_tgen_five", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgen_six = new TH1D("nue_tgen_six", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgen_seven = new TH1D("nue_tgen_seven", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgen_eight = new TH1D("nue_tgen_eight", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgen_nine = new TH1D("nue_tgen_nine", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgen_ten = new TH1D("nue_tgen_ten", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgen_eleven = new TH1D("nue_tgen_eleven", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    nue_tgen_twelve = new TH1D("nue_tgen_twelve", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
    
    numu_tgen_two = new TH1D("numu_tgen_two", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgen_three = new TH1D("numu_tgen_three", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgen_four = new TH1D("numu_tgen_four", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgen_five = new TH1D("numu_tgen_five", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgen_six = new TH1D("numu_tgen_six", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgen_seven = new TH1D("numu_tgen_seven", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgen_eight = new TH1D("numu_tgen_eight", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgen_nine = new TH1D("numu_tgen_nine", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgen_ten = new TH1D("numu_tgen_ten", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgen_eleven = new TH1D("numu_tgen_eleven", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    numu_tgen_twelve = new TH1D("numu_tgen_twelve", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
    
    
    
    // new (particles)
    pionplus2DHisto = new TH2D("pionplus2DHisto",  "pionplus2DHisto ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    pionplusAngleHisto = new TH1D("pionplusAngleHisto",  "pionplusAngleHisto ; theta [mrad] ;  no. of entries ", 600, 0, 3000 );
    pionplusMomentumHisto = new TH1D("pionplusMomentumHisto",  "pionplusMomentumHisto ; p [GeV/c]; no. of entries ", 600, 0, 400 );
    pionplus_nue= new TH2D("pionplus_nue",  "pionplus_nue ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    pionplus_numu = new TH2D("pionplus_numu",  "pionplus_numu ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    pionplus_anue= new TH2D("pionplus_anue",  "pionplus_anue ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    pionplus_anumu = new TH2D("pionplus_anumu",  "pionplus_anumu ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    
    pionminus2DHisto = new TH2D("pionminus2DHisto",  "pionminus2DHisto ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    pionminusAngleHisto = new TH1D("pionminusAngleHisto",  "pionminusAngleHisto ; theta [mrad] ;  no. of entries ", 600, 0, 3000 );
    pionminusMomentumHisto = new TH1D("pionminusMomentumHisto",  "pionminusMomentumHisto ; p [GeV/c]; no. of entries ", 600, 0, 400 );
    pionminus_nue= new TH2D("pionminus_nue",  "pionminus_nue ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    pionminus_numu = new TH2D("pionminus_numu",  "pionminus_numu ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    pionminus_anue= new TH2D("pionminus_anue",  "pionminus_anue ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    pionminus_anumu = new TH2D("pionminus_anumu",  "pionminus_anumu ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    
    Klong2DHisto = new TH2D("Klong2DHisto",  "Klong2DHisto ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    KlongAngleHisto = new TH1D("KlongAngleHisto",  "KlongAngleHisto ; theta [mrad] ;  no. of entries ", 600, 0, 3000 );
    KlongMomentumHisto = new TH1D("KlongMomentumHisto",  "KlongMomentumHisto ; p [GeV/c]; no. of entries ", 600, 0, 400 );
    Klong_nue= new TH2D("Klong_nue",  "Klong_nue ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    Klong_numu = new TH2D("Klong_numu",  "Klong_numu ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    Klong_anue= new TH2D("Klong_anue",  "Klong_anue ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    Klong_anumu = new TH2D("Klong_anumu",  "Klong_anumu ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    
    Kplus2DHisto = new TH2D("Kplus2DHisto",  "KplusHisto ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    KplusAngleHisto = new TH1D("KplusAngleHisto",  "KplusAngleHisto ; theta [mrad] ;  no. of entries ", 600, 0, 3000 );
    KplusMomentumHisto = new TH1D("KplusMomentumHisto",  "KplusMomentumHisto ; p [GeV/c]; no. of entries ", 600, 0, 400 );
    Kplus_nue= new TH2D("Kplus_nue",  "Kplus_nue ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    Kplus_numu = new TH2D("Kplus_numu",  "Kplus_numu ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    Kplus_anue= new TH2D("Kplus_anue",  "Kplus_anue ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    Kplus_anumu = new TH2D("Kplus_anumu",  "Kplus_anumu ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    
    Kminus2DHisto = new TH2D("Kminus2DHisto",  "Kminus2DHisto ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    KminusAngleHisto = new TH1D("KminusAngleHisto",  "KminusAngleHisto ; theta [mrad] ;  no. of entries ", 600, 0, 3000 );
    KminusMomentumHisto = new TH1D("KminusMomentumHisto",  "KminusMomentumHisto ; p [GeV/c]; no. of entries ", 600, 0, 400 );
    Kminus_nue= new TH2D("Kminus_nue",  "Kminus_nue ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    Kminus_numu = new TH2D("Kminus_numu",  "Kminus_numu ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    Kminus_anue= new TH2D("Kminus_anue",  "Kminus_anue ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    Kminus_anumu = new TH2D("Kminus_anumu",  "Kminus_anumu ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    
    muonminus2DHisto = new TH2D("muonminus2DHisto",  "muonminus2DHisto ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    muonminusAngleHisto = new TH1D("muonminusAngleHisto",  "muonminusAngleHisto ; theta [mrad] ;  no. of entries ", 600, 0, 3000 );
    muonminusMomentumHisto = new TH1D("muonminusMomentumHisto",  "muonminusMomentumHisto ; p [GeV/c]; no. of entries ", 600, 0, 400 );
    
    muonplus2DHisto = new TH2D("muonplus2DHisto",  "muonplus2DHisto ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    muonplusAngleHisto = new TH1D("muonplusAngleHisto",  "muonplusAngleHisto ; theta [mrad] ;  no. of entries ", 600, 0, 3000 );
    muonplusMomentumHisto = new TH1D("muonplusMomentumHisto",  "muonplusMomentumHisto ; p [GeV/c]; no. of entries ", 600, 0, 400 );
    
    proton2DHisto = new TH2D("proton2DHisto",  "proton2DHisto ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    protonAngleHisto = new TH1D("protonAngleHisto",  "protonAngleHisto ; theta [mrad] ;  no. of entries ", 600, 0, 3000 );
    protonMomentumHisto = new TH1D("protonMomentumHisto",  "protonMomentumHisto ; p [GeV/c]; no. of entries ", 600, 0, 400 );
    proton_nue= new TH2D("proton_nue",  "proton_nue ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    proton_numu = new TH2D("proton_numu",  "proton_numu ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    proton_anue= new TH2D("proton_anue",  "proton_anue ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    proton_anumu = new TH2D("proton_anumu",  "proton_anumu ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    
    neutron2DHisto = new TH2D("neutron2DHisto",  "neutron2DHisto ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    neutronAngleHisto = new TH1D("neutronAngleHisto",  "neutronAngleHisto ; theta [mrad] ;  no. of entries ", 600, 0, 3000 );
    neutronMomentumHisto = new TH1D("neutronMomentumHisto",  "neutronMomentumHisto ; p [GeV/c]; no. of entries ", 600, 0, 400 );
    neutron_nue= new TH2D("neutron_nue",  "neutron_nue ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    neutron_numu = new TH2D("neutron_numu",  "neutron_numu ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    neutron_anue= new TH2D("neutron_anue",  "neutron_anue ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    neutron_anumu = new TH2D("neutron_anumu",  "neutron_anumu ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    
    lambda2DHisto = new TH2D("lambda2DHisto",  "lambda2DHisto ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    lambdaAngleHisto = new TH1D("lambdaAngleHisto",  "lambdaAngleHisto ; theta [mrad] ;  no. of entries ", 600, 0, 3000 );
    lambdaMomentumHisto = new TH1D("lambdaMomentumHisto",  "lambdaMomentumHisto ; p [GeV/c]; no. of entries ", 600, 0, 400 );
    
    sigmaplus2DHisto = new TH2D("sigmaplus2DHisto",  "sigmaplus2DHisto ;p [GeV/c] ;  theta [mrad]", 600, 0, 400, 600, 0, 3000 );
    sigmaplusAngleHisto = new TH1D("sigmaplusAngleHisto",  "sigmaplusAngleHisto ; theta [mrad] ;  no. of entries ", 600, 0, 3000 );
    sigmaplusMomentumHisto = new TH1D("sigmaplusMomentumHisto",  "sigmaplusMomentumHisto ; p [GeV/c]; no. of entries ", 600, 0, 400 );
    
    ///////////////////////////////// decays /////////////////////////////
    ppnue_dm_one = new TH2D("ppnue_dm_one",  "ppnue_dm_one ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnue_dm_two = new TH2D("ppnue_dm_two",  "ppnue_dm_two ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnue_dm_three = new TH2D("ppnue_dm_three",  "ppnue_dm_three ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnue_dm_four = new TH2D("ppnue_dm_four",  "ppnue_dm_four ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnue_dm_five = new TH2D("ppnue_dm_five",  "ppnue_dm_five ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnue_dm_six = new TH2D("ppnue_dm_six",  "ppnue_dm_six ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnue_dm_seven = new TH2D("ppnue_dm_seven",  "ppnue_dm_seven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnue_dm_eight = new TH2D("ppnue_dm_eight",  "ppnue_dm_eight ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnue_dm_nine = new TH2D("ppnue_dm_nine",  "ppnue_dm_nine ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnue_dm_ten = new TH2D("ppnue_dm_ten",  "ppnue_dm_ten ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnue_dm_eleven = new TH2D("ppnue_dm_eleven",  "ppnue_dm_eleven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnue_dm_twelve = new TH2D("ppnue_dm_twelve",  "ppnue_dm_twelve ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnue_dm_thirteen = new TH2D("ppnue_dm_thirteen",  "ppnue_dm_thirteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnue_dm_fourteen = new TH2D("ppnue_dm_fourteen",  "ppnue_dm_fourteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    
    pmnue_dm_one = new TH2D("pmnue_dm_one",  "pmnue_dm_one ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnue_dm_two = new TH2D("pmnue_dm_two",  "pmnue_dm_two ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnue_dm_three = new TH2D("pmnue_dm_three",  "pmnue_dm_three ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnue_dm_four = new TH2D("pmnue_dm_four",  "pmnue_dm_four ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnue_dm_five = new TH2D("pmnue_dm_five",  "pmnue_dm_five ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnue_dm_six = new TH2D("pmnue_dm_six",  "pmnue_dm_six ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnue_dm_seven = new TH2D("pmnue_dm_seven",  "pmnue_dm_seven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnue_dm_eight = new TH2D("pmnue_dm_eight",  "pmnue_dm_eight ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnue_dm_nine = new TH2D("pmnue_dm_nine",  "pmnue_dm_nine ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnue_dm_ten = new TH2D("pmnue_dm_ten",  "pmnue_dm_ten ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnue_dm_eleven = new TH2D("pmnue_dm_eleven",  "pmnue_dm_eleven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnue_dm_twelve = new TH2D("pmnue_dm_twelve",  "pmnue_dm_twelve ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnue_dm_thirteen = new TH2D("pmnue_dm_thirteen",  "pmnue_dm_thirteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnue_dm_fourteen = new TH2D("pmnue_dm_fourteen",  "pmnue_dm_fourteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    
    klnue_dm_one = new TH2D("klnue_dm_one",  "klnue_dm_one ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnue_dm_two = new TH2D("klnue_dm_two",  "klnue_dm_two ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnue_dm_three = new TH2D("klnue_dm_three",  "klnue_dm_three ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnue_dm_four = new TH2D("klnue_dm_four",  "klnue_dm_four ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnue_dm_five = new TH2D("klnue_dm_five",  "klnue_dm_five ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnue_dm_six = new TH2D("klnue_dm_six",  "klnue_dm_six ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnue_dm_seven = new TH2D("klnue_dm_seven",  "klnue_dm_seven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnue_dm_eight = new TH2D("klnue_dm_eight",  "klnue_dm_eight ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnue_dm_nine = new TH2D("klnue_dm_nine",  "klnue_dm_nine ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnue_dm_ten = new TH2D("klnue_dm_ten",  "klnue_dm_ten ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnue_dm_eleven = new TH2D("klnue_dm_eleven",  "klnue_dm_eleven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnue_dm_twelve = new TH2D("klnue_dm_twelve",  "klnue_dm_twelve ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnue_dm_thirteen = new TH2D("klnue_dm_thirteen",  "klnue_dm_thirteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnue_dm_fourteen = new TH2D("klnue_dm_fourteen",  "klnue_dm_fourteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    
    kpnue_dm_one = new TH2D("kpnue_dm_one",  "kpnue_dm_one ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnue_dm_two = new TH2D("kpnue_dm_two",  "kpnue_dm_two ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnue_dm_three = new TH2D("kpnue_dm_three",  "kpnue_dm_three ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnue_dm_four = new TH2D("kpnue_dm_four",  "kpnue_dm_four ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnue_dm_five = new TH2D("kpnue_dm_five",  "kpnue_dm_five ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnue_dm_six = new TH2D("kpnue_dm_six",  "kpnue_dm_six ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnue_dm_seven = new TH2D("kpnue_dm_seven",  "kpnue_dm_seven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnue_dm_eight = new TH2D("kpnue_dm_eight",  "kpnue_dm_eight ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnue_dm_nine = new TH2D("kpnue_dm_nine",  "kpnue_dm_nine ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnue_dm_ten = new TH2D("kpnue_dm_ten",  "kpnue_dm_ten ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnue_dm_eleven = new TH2D("kpnue_dm_eleven",  "kpnue_dm_eleven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnue_dm_twelve = new TH2D("kpnue_dm_twelve",  "kpnue_dm_twelve ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnue_dm_thirteen = new TH2D("kpnue_dm_thirteen",  "kpnue_dm_thirteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnue_dm_fourteen = new TH2D("kpnue_dm_fourteen",  "kpnue_dm_fourteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );

    kmnue_dm_one = new TH2D("kmnue_dm_one",  "kmnue_dm_one ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnue_dm_two = new TH2D("kmnue_dm_two",  "kmnue_dm_two ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnue_dm_three = new TH2D("kmnue_dm_three",  "kmnue_dm_three ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnue_dm_four = new TH2D("kmnue_dm_four",  "kmnue_dm_four ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnue_dm_five = new TH2D("kmnue_dm_five",  "kmnue_dm_five ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnue_dm_six = new TH2D("kmnue_dm_six",  "kmnue_dm_six ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnue_dm_seven = new TH2D("kmnue_dm_seven",  "kmnue_dm_seven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnue_dm_eight = new TH2D("kmnue_dm_eight",  "kmnue_dm_eight ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnue_dm_nine = new TH2D("kmnue_dm_nine",  "kmnue_dm_nine ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnue_dm_ten = new TH2D("kmnue_dm_ten",  "kmnue_dm_ten ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnue_dm_eleven = new TH2D("kmnue_dm_eleven",  "kmnue_dm_eleven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnue_dm_twelve = new TH2D("kmnue_dm_twelve",  "kmnue_dm_twelve ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnue_dm_thirteen = new TH2D("kmnue_dm_thirteen",  "kmnue_dm_thirteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnue_dm_fourteen = new TH2D("kmnue_dm_fourteen",  "kmnue_dm_fourteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    
    pnue_dm_one = new TH2D("pnue_dm_one",  "pnue_dm_one ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnue_dm_two = new TH2D("pnue_dm_two",  "pnue_dm_two ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnue_dm_three = new TH2D("pnue_dm_three",  "pnue_dm_three ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnue_dm_four = new TH2D("pnue_dm_four",  "pnue_dm_four ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnue_dm_five = new TH2D("pnue_dm_five",  "pnue_dm_five ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnue_dm_six = new TH2D("pnue_dm_six",  "pnue_dm_six ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnue_dm_seven = new TH2D("pnue_dm_seven",  "pnue_dm_seven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnue_dm_eight = new TH2D("pnue_dm_eight",  "pnue_dm_eight ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnue_dm_nine = new TH2D("pnue_dm_nine",  "pnue_dm_nine ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnue_dm_ten = new TH2D("pnue_dm_ten",  "pnue_dm_ten ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnue_dm_eleven = new TH2D("pnue_dm_eleven",  "pnue_dm_eleven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnue_dm_twelve = new TH2D("pnue_dm_twelve",  "pnue_dm_twelve ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnue_dm_thirteen = new TH2D("pnue_dm_thirteen",  "pnue_dm_thirteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnue_dm_fourteen = new TH2D("pnue_dm_fourteen",  "pnue_dm_fourteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    
    nnue_dm_one = new TH2D("nnue_dm_one",  "nnue_dm_one ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnue_dm_two = new TH2D("nnue_dm_two",  "nnue_dm_two ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnue_dm_three = new TH2D("nnue_dm_three",  "nnue_dm_three ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnue_dm_four = new TH2D("nnue_dm_four",  "nnue_dm_four ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnue_dm_five = new TH2D("nnue_dm_five",  "nnue_dm_five ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnue_dm_six = new TH2D("nnue_dm_six",  "nnue_dm_six ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnue_dm_seven = new TH2D("nnue_dm_seven",  "nnue_dm_seven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnue_dm_eight = new TH2D("nnue_dm_eight",  "nnue_dm_eight ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnue_dm_nine = new TH2D("nnue_dm_nine",  "nnue_dm_nine ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnue_dm_ten = new TH2D("nnue_dm_ten",  "nnue_dm_ten ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnue_dm_eleven = new TH2D("nnue_dm_eleven",  "nnue_dm_eleven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnue_dm_twelve = new TH2D("nnue_dm_twelve",  "nnue_dm_twelve ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnue_dm_thirteen = new TH2D("nnue_dm_thirteen",  "nnue_dm_thirteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnue_dm_fourteen = new TH2D("nnue_dm_fourteen",  "nnue_dm_fourteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    
    ppnumu_dm_one = new TH2D("ppnumu_dm_one",  "ppnumu_dm_one ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnumu_dm_two = new TH2D("ppnumu_dm_two",  "ppnumu_dm_two ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnumu_dm_three = new TH2D("ppnumu_dm_three",  "ppnumu_dm_three ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnumu_dm_four = new TH2D("ppnumu_dm_four",  "ppnumu_dm_four ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnumu_dm_five = new TH2D("ppnumu_dm_five",  "ppnumu_dm_five ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnumu_dm_six = new TH2D("ppnumu_dm_six",  "ppnumu_dm_six ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnumu_dm_seven = new TH2D("ppnumu_dm_seven",  "ppnumu_dm_seven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnumu_dm_eight = new TH2D("ppnumu_dm_eight",  "ppnumu_dm_eight ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnumu_dm_nine = new TH2D("ppnumu_dm_nine",  "ppnumu_dm_nine ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnumu_dm_ten = new TH2D("ppnumu_dm_ten",  "ppnumu_dm_ten ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnumu_dm_eleven = new TH2D("ppnumu_dm_eleven",  "ppnumu_dm_eleven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnumu_dm_twelve = new TH2D("ppnumu_dm_twelve",  "ppnumu_dm_twelve ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnumu_dm_thirteen = new TH2D("ppnumu_dm_thirteen",  "ppnumu_dm_thirteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    ppnumu_dm_fourteen = new TH2D("ppnumu_dm_fourteen",  "ppnumu_dm_fourteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    
    pmnumu_dm_one = new TH2D("pmnumu_dm_one",  "pmnumu_dm_one ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnumu_dm_two = new TH2D("pmnumu_dm_two",  "pmnumu_dm_two ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnumu_dm_three = new TH2D("pmnumu_dm_three",  "pmnumu_dm_three ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnumu_dm_four = new TH2D("pmnumu_dm_four",  "pmnumu_dm_four ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnumu_dm_five = new TH2D("pmnumu_dm_five",  "pmnumu_dm_five ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnumu_dm_six = new TH2D("pmnumu_dm_six",  "pmnumu_dm_six ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnumu_dm_seven = new TH2D("pmnumu_dm_seven",  "pmnumu_dm_seven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnumu_dm_eight = new TH2D("pmnumu_dm_eight",  "pmnumu_dm_eight ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnumu_dm_nine = new TH2D("pmnumu_dm_nine",  "pmnumu_dm_nine ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnumu_dm_ten = new TH2D("pmnumu_dm_ten",  "pmnumu_dm_ten ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnumu_dm_eleven = new TH2D("pmnumu_dm_eleven",  "pmnumu_dm_eleven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnumu_dm_twelve = new TH2D("pmnumu_dm_twelve",  "pmnumu_dm_twelve ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnumu_dm_thirteen = new TH2D("pmnumu_dm_thirteen",  "pmnumu_dm_thirteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pmnumu_dm_fourteen = new TH2D("pmnumu_dm_fourteen",  "pmnumu_dm_fourteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    
    klnumu_dm_one = new TH2D("klnumu_dm_one",  "klnumu_dm_one ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnumu_dm_two = new TH2D("klnumu_dm_two",  "klnumu_dm_two ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnumu_dm_three = new TH2D("klnumu_dm_three",  "klnumu_dm_three ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnumu_dm_four = new TH2D("klnumu_dm_four",  "klnumu_dm_four ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnumu_dm_five = new TH2D("klnumu_dm_five",  "klnumu_dm_five ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnumu_dm_six = new TH2D("klnumu_dm_six",  "klnumu_dm_six ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnumu_dm_seven = new TH2D("klnumu_dm_seven",  "klnumu_dm_seven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnumu_dm_eight = new TH2D("klnumu_dm_eight",  "klnumu_dm_eight ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnumu_dm_nine = new TH2D("klnumu_dm_nine",  "klnumu_dm_nine ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnumu_dm_ten = new TH2D("klnumu_dm_ten",  "klnumu_dm_ten ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnumu_dm_eleven = new TH2D("klnumu_dm_eleven",  "klnumu_dm_eleven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnumu_dm_twelve = new TH2D("klnumu_dm_twelve",  "klnumu_dm_twelve ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnumu_dm_thirteen = new TH2D("klnumu_dm_thirteen",  "klnumu_dm_thirteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    klnumu_dm_fourteen = new TH2D("klnumu_dm_fourteen",  "klnumu_dm_fourteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    
    kpnumu_dm_one = new TH2D("kpnumu_dm_one",  "kpnumu_dm_one ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnumu_dm_two = new TH2D("kpnumu_dm_two",  "kpnumu_dm_two ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnumu_dm_three = new TH2D("kpnumu_dm_three",  "kpnumu_dm_three ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnumu_dm_four = new TH2D("kpnumu_dm_four",  "kpnumu_dm_four ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnumu_dm_five = new TH2D("kpnumu_dm_five",  "kpnumu_dm_five ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnumu_dm_six = new TH2D("kpnumu_dm_six",  "kpnumu_dm_six ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnumu_dm_seven = new TH2D("kpnumu_dm_seven",  "kpnumu_dm_seven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnumu_dm_eight = new TH2D("kpnumu_dm_eight",  "kpnumu_dm_eight ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnumu_dm_nine = new TH2D("kpnumu_dm_nine",  "kpnumu_dm_nine ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnumu_dm_ten = new TH2D("kpnumu_dm_ten",  "kpnumu_dm_ten ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnumu_dm_eleven = new TH2D("kpnumu_dm_eleven",  "kpnumu_dm_eleven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnumu_dm_twelve = new TH2D("kpnumu_dm_twelve",  "kpnumu_dm_twelve ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnumu_dm_thirteen = new TH2D("kpnumu_dm_thirteen",  "kpnumu_dm_thirteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kpnumu_dm_fourteen = new TH2D("kpnumu_dm_fourteen",  "kpnumu_dm_fourteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    
    kmnumu_dm_one = new TH2D("kmnumu_dm_one",  "kmnumu_dm_one ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnumu_dm_two = new TH2D("kmnumu_dm_two",  "kmnumu_dm_two ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnumu_dm_three = new TH2D("kmnumu_dm_three",  "kmnumu_dm_three ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnumu_dm_four = new TH2D("kmnumu_dm_four",  "kmnumu_dm_four ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnumu_dm_five = new TH2D("kmnumu_dm_five",  "kmnumu_dm_five ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnumu_dm_six = new TH2D("kmnumu_dm_six",  "kmnumu_dm_six ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnumu_dm_seven = new TH2D("kmnumu_dm_seven",  "kmnumu_dm_seven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnumu_dm_eight = new TH2D("kmnumu_dm_eight",  "kmnumu_dm_eight ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnumu_dm_nine = new TH2D("kmnumu_dm_nine",  "kmnumu_dm_nine ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnumu_dm_ten = new TH2D("kmnumu_dm_ten",  "kmnumu_dm_ten ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnumu_dm_eleven = new TH2D("kmnumu_dm_eleven",  "kmnumu_dm_eleven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnumu_dm_twelve = new TH2D("kmnumu_dm_twelve",  "kmnumu_dm_twelve ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnumu_dm_thirteen = new TH2D("kmnumu_dm_thirteen",  "kmnumu_dm_thirteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    kmnumu_dm_fourteen = new TH2D("kmnumu_dm_fourteen",  "kmnumu_dm_fourteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    
    pnumu_dm_one = new TH2D("pnumu_dm_one",  "pnumu_dm_one ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnumu_dm_two = new TH2D("pnumu_dm_two",  "pnumu_dm_two ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnumu_dm_three = new TH2D("pnumu_dm_three",  "pnumu_dm_three ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnumu_dm_four = new TH2D("pnumu_dm_four",  "pnumu_dm_four ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnumu_dm_five = new TH2D("pnumu_dm_five",  "pnumu_dm_five ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnumu_dm_six = new TH2D("pnumu_dm_six",  "pnumu_dm_six ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnumu_dm_seven = new TH2D("pnumu_dm_seven",  "pnumu_dm_seven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnumu_dm_eight = new TH2D("pnumu_dm_eight",  "pnumu_dm_eight ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnumu_dm_nine = new TH2D("pnumu_dm_nine",  "pnumu_dm_nine ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnumu_dm_ten = new TH2D("pnumu_dm_ten",  "pnumu_dm_ten ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnumu_dm_eleven = new TH2D("pnumu_dm_eleven",  "pnumu_dm_eleven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnumu_dm_twelve = new TH2D("pnumu_dm_twelve",  "pnumu_dm_twelve ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnumu_dm_thirteen = new TH2D("pnumu_dm_thirteen",  "pnumu_dm_thirteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    pnumu_dm_fourteen = new TH2D("pnumu_dm_fourteen",  "pnumu_dm_fourteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    
    nnumu_dm_one = new TH2D("nnumu_dm_one",  "nnumu_dm_one ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnumu_dm_two = new TH2D("nnumu_dm_two",  "nnumu_dm_two ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnumu_dm_three = new TH2D("nnumu_dm_three",  "nnumu_dm_three ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnumu_dm_four = new TH2D("nnumu_dm_four",  "nnumu_dm_four ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnumu_dm_five = new TH2D("nnumu_dm_five",  "nnumu_dm_five ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnumu_dm_six = new TH2D("nnumu_dm_six",  "nnumu_dm_six ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnumu_dm_seven = new TH2D("nnumu_dm_seven",  "nnumu_dm_seven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnumu_dm_eight = new TH2D("nnumu_dm_eight",  "nnumu_dm_eight ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnumu_dm_nine = new TH2D("nnumu_dm_nine",  "nnumu_dm_nine ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnumu_dm_ten = new TH2D("nnumu_dm_ten",  "nnumu_dm_ten ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnumu_dm_eleven = new TH2D("nnumu_dm_eleven",  "nnumu_dm_eleven ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnumu_dm_twelve = new TH2D("nnumu_dm_twelve",  "nnumu_dm_twelve ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnumu_dm_thirteen = new TH2D("nnumu_dm_thirteen",  "nnumu_dm_thirteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    nnumu_dm_fourteen = new TH2D("nnumu_dm_fourteen",  "nnumu_dm_fourteen ;p [GeV/c] ;  theta [mrad]", 600, 0, 200, 600, 0, 3000 );
    
    
    
    ////////////////////////////// MIPP ///////////////////////////////////
    pionplus_MIPP =  new TH2D("pionplus_MIPP",  "pionplus_MIPP ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
     pionplus_MIPP_nue =  new TH2D("pionplus_MIPP_nue",  "pionplus_MIPP_nue ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
     pionplus_MIPP_numu =  new TH2D("pionplus_MIPP_numu",  "pionplus_MIPP_numu ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
     pionplus_MIPP_anue =  new TH2D("pionplus_MIPP_anue",  "pionplus_MIPP_anue ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
     pionplus_MIPP_anumu =  new TH2D("pionplus_MIPP_anumu",  "pionplus_MIPP_anumu ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    
    pionminus_MIPP =  new TH2D("pionminus_MIPP",  "pionminus_MIPP ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    pionminus_MIPP_nue =  new TH2D("pionminus_MIPP_nue",  "pionminus_MIPP_nue ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    pionminus_MIPP_numu =  new TH2D("pionminus_MIPP_numu",  "pionminus_MIPP_numu ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    pionminus_MIPP_anue =  new TH2D("pionminus_MIPP_anue",  "pionminus_MIPP_anue ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    pionminus_MIPP_anumu =  new TH2D("pionminus_MIPP_anumu",  "pionminus_MIPP_anumu ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    
    Kplus_MIPP =  new TH2D("Kplus_MIPP",  "Kplus_MIPP ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    Kplus_MIPP_nue =  new TH2D("Kplus_MIPP_nue",  "Kplus_MIPP_nue ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    Kplus_MIPP_numu =  new TH2D("Kplus_MIPP_numu",  "Kplus_MIPP_numu ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    Kplus_MIPP_anue =  new TH2D("Kplus_MIPP_anue",  "Kplus_MIPP_anue ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    Kplus_MIPP_anumu =  new TH2D("Kplus_MIPP_anumu",  "Kplus_MIPP_anumu ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    
    Kminus_MIPP =  new TH2D("Kminus_MIPP",  "Kminus_MIPP ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    Kminus_MIPP_nue =  new TH2D("Kminus_MIPP_nue",  "Kminus_MIPP_nue ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    Kminus_MIPP_numu =  new TH2D("Kminus_MIPP_numu",  "Kminus_MIPP_numu ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    Kminus_MIPP_anue =  new TH2D("Kminus_MIPP_anue",  "Kminus_MIPP_anue ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    Kminus_MIPP_anumu =  new TH2D("Kminus_MIPP_anumu",  "Kminus_MIPP_anumu ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
     Kminus_MIPP_nue_dm_one =  new TH2D("Kminus_MIPP_nue_dm_one",  "Kminus_MIPP_nue ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    Kminus_MIPP_nue_dm_six =  new TH2D("Kminus_MIPP_nue_dm_six",  "Kminus_MIPP_nue ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    
    Kminus_MIPP_nue_dm_eleven =  new TH2D("Kminus_MIPP_nue_dm_eleven",  "Kminus_MIPP_nue ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    
    
    // feynman x histos
    
    pionplus_NA49 =  new TH2D("pionplus_NA49",  "pionplus_NA49 ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 200, -1, 1, 100, 0, 10 );
    pionplus_NA49_nue =  new TH2D("pionplus_NA49_nue",  "pionplus_NA49_nue ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 200, -1, 1, 100, 0, 10 );
    pionplus_MIPP_nue_dm1 =  new TH2D("pionplus_MIPP_nue_dm1",  "pionplus_MIPP_nue_dm1 ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    pionplus_MIPP_nue_dm5 =  new TH2D("pionplus_MIPP_nue_dm5",  "pionplus_MIPP_nue_dm5 ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    pionplus_MIPP_nue_dm6 =  new TH2D("pionplus_MIPP_nue_dm6",  "pionplus_MIPP_nue_dm6 ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
     pionplus_MIPP_nue_dm11 =  new TH2D("pionplus_MIPP_nue_dm11",  "pionplus_MIPP_nue_dm11 ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 600, 0, 200, 200, 0, 200 );
    
    pionminus_NA49 =  new TH2D("pionminus_NA49",  "pionminus_NA49 ; xF ;  p_{T} [GeV/c]", 200, -1, 1, 100, 0, 10 );
    pionminus_NA49_nue =  new TH2D("pionminus_NA49_nue",  "pionminus_NA49_nue ; xF ;  p_{T} [GeV/c]", 200, -1, 1, 100, 0, 10 );
    
    Kplus_NA49 =  new TH2D("Kplus_NA49",  "Kplus_NA49 ;xF ;  p_{T} [GeV/c]", 200, -1, 1, 100, 0, 10 );
    Kplus_NA49_second =  new TH2D("Kplus_NA49_second ",  "Kplus_NA49_second ;xF ;  p_{T} [GeV/c]", 200, -1, 1, 100, 0, 10);
    Kplus_NA49_nue =  new TH2D("Kplus_NA49_nue",  "Kplus_NA49_nue ;p_{z} [GeV/c] ;  p_{T} [GeV/c]", 200, -1, 1, 100, 0, 10 );
    
    Kminus_NA49 =  new TH2D("Kminus_NA49",  "Kminus_NA49 ; xF;  p_{T} [GeV/c]", 200, -1, 1, 100, 0, 10 );
    Kminus_NA49_nue =  new TH2D("Kminus_NA49_nue",  "Kminus_NA49_nue ; xF ;  p_{T} [GeV/c]", 200, -1, 1, 100, 0, 10 );
    
    
    proton_NA49 =  new TH2D("proton_NA49",  "proton_NA49 ; xF;  p_{T} [GeV/c]", 200, -1, 1, 100, 0, 10 );
    proton_NA49_nue =  new TH2D("proton_NA49_nue",  "proton_NA49_nue ; xF ;  p_{T} [GeV/c]", 200, -1, 1, 100, 0, 10 );
    
    neutron_NA49 =  new TH2D("neutron_NA49",  "neutron_NA49 ; xF;  p_{T} [GeV/c]", 200, -1, 1, 100, 0, 10 );
    neutron_NA49_nue =  new TH2D("neutron_NA49_nue",  "neutron_NA49_nue ; xF ;  p_{T} [GeV/c]", 200, -1, 1, 100, 0, 10 );
    
    antiproton_NA49 =  new TH2D("antiproton_NA49",  "antiproton_NA49 ; xF;  p_{T} [GeV/c]", 200, -1, 1, 100, 0, 10 );
    antiproton_NA49_nue =  new TH2D("antiproton_NA49_nue",  "antiproton_NA49_nue ; xF ;  p_{T} [GeV/c]", 200, -1, 1, 100, 0, 10 );
    
    
    
    
    
}

NuMIFlux::~NuMIFlux() {

}

#endif // #ifdef NuMIFlux_cxx
