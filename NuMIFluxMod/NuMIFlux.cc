#define NuMIFlux_cxx
#include <iostream>	
#include <iomanip>
#include <string>
	
using namespace std;
	
#include "TChain.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TRotation.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLorentzVector.h"

#include "NuMIFlux.hh"
#include "FluggNtuple/FluxNtuple.h"

//#include "dk2nu.h"
//#include "dkmeta.h"
//#include "calcLocationWeights.cxx"
//#include "FluxNtuple_C.so"

double calculateXf(double sqrts,double mass,double pz,double pt, double ZBoost)
{
    
    if(!sqrts)
        return 0;
    
    double ptot=TMath::Sqrt(pz*pz+pt*pt);
    
    double Energy=mass*sqrt(1+(ptot/mass)*(ptot/mass));
    
    TLorentzVector Particle;
    
    Particle.SetPz(pz);
    Particle.SetPy(pt);
    Particle.SetE(Energy);
    
    Particle.Boost(0,0,-ZBoost);
    
    
    return 2*Particle.Pz()/sqrts; ///////////////
    
}

double fixedTargetBoost(double incidentPz, double incidentmass, double targetmass)
{
    
    double incidentPt=0;  //in case we need it in the future.
    
    double ptot=TMath::Sqrt(incidentPz*incidentPz+incidentPt*incidentPt);
    
    double Energy=incidentmass*sqrt(1+(ptot/incidentmass)*(ptot/incidentmass));
    
    return  incidentPz/(Energy+targetmass);
}

void NuMIFlux::CalculateFlux() {

  fluxNtuple = new FluxNtuple(cflux);  // 

  //***************************************
  //
  //  Loop over the entries.
  //
  //***************************************

  Long64_t nflux = cflux->GetEntries();
  std::cout << "Total number of entries: " << nflux << std::endl;
  for (Long64_t i=0; i < 10000000; ++i ) {

    // Get entry i. fluxNtuple is now filled with entry i info.
    cflux->GetEntry(i);

    // Alert the user
    if (i % 10000000 == 0) cout << "On entry " << i << endl;
    if(treeNumber != cflux->GetTreeNumber()) {
      treeNumber = cflux->GetTreeNumber();
      std::cout << "Moving to tree number " << treeNumber << "." << std::endl;	
      AccumulatedPOT += estimate_pots(highest_evtno);
      cout << "AccumulatedPOT: " << AccumulatedPOT << endl;
    }

    double wgt_xy = 0.;  // neutrino weight
    double enu    = 0.;  // neutrino energy in lab frame

    // Pick a random point in the TPC (in detector coordinates)
    TVector3 xyz_det = RandomInTPC();
    if (debug) cout << "xyz_det = [" << xyz_det.X() << ", " << xyz_det.Y() << ", " << xyz_det.Z() << "]" << endl;

    // From detector to beam coordinates
    TVector3 xyz_beam = FromDetToBeam(xyz_det);
    if (debug) cout << "xyz_beam = [" << xyz_beam.X() << ", " << xyz_beam.Y() << ", " << xyz_beam.Z() << "]" << endl;   
 
    // Calculate the weight
    int ret = calcEnuWgt(fluxNtuple, xyz_beam, enu, wgt_xy);     
    if (ret != 0) cout << "Error with calcEnuWgt. Return " << ret << endl;
    if (debug) cout << "wgt_xy " << wgt_xy << endl;

    // Calculate the total weight
    double weight = wgt_xy * fluxNtuple->Nimpwt * fDefaultWeightCorrection;
      
      
      // rescaling the weight based on MIPP constrains
      double pZ = fluxNtuple->tpz;
      double pT = sqrt( pow(fluxNtuple->tpx,2) + pow(fluxNtuple->tpy,2));
   
          
     // double weight_up = 0
    double parent_mom2 = ( fluxNtuple->pdPx*fluxNtuple->pdPx +
                        fluxNtuple->pdPy*fluxNtuple->pdPy +
                       fluxNtuple->pdPz*fluxNtuple->pdPz );
      
      
      
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     // std::cout<<"particle type: "<<fluxNtuple->tgptype<<endl;
      // momentum of particles following the initial proton interaction; the ones that produce the neutrino parent
      double mc_flux_px = fluxNtuple->tgppx;
      double mc_flux_py = fluxNtuple->tgppy;
      double mc_flux_pz = fluxNtuple-> tgppz;
      double momentum = sqrt((mc_flux_px * mc_flux_px) + (mc_flux_py * mc_flux_py) + (mc_flux_pz * mc_flux_pz));
     // std::cout<<"momentum: "<<momentum<<endl;
      
      //create a unit vector in the Z direction - assuming z is the beam
      std::vector<double> mc_flux_unit;
      mc_flux_unit.push_back(0);
      mc_flux_unit.push_back(0);
      mc_flux_unit.push_back(1);
      
      //now calculate the angle between the unit vector and the momentum vector
      double dot_prod = (mc_flux_unit.at(0) * mc_flux_px) + (mc_flux_unit.at(1) * mc_flux_py) + (mc_flux_unit.at(2) * mc_flux_pz);
      double angle = acos(dot_prod / (1 * std::abs(momentum) ) ) *1000;
     // std::cout<<"angle: "<<angle <<"\n"<<std::endl;
      //double momentum_transverse =sqrt( pow(mc_flux_px,2) + pow(mc_flux_py , 2) );
      
      
      

      
      //// FEYNMAN SCALING VARIABLE (x_F) calcs /////
      
      /////// Andrzej's xF Calcs ///////
     
      ///////////////// needed setting part
      double Mass_tgp =0.0 ;
      switch (fluxNtuple -> tgptype){
          case 8: //pionplus
          case 9: //pion minus
              Mass_tgp = 0.13957;
              break;
              
          case 11: //kaonplus
          case 12: //kaonminus
              Mass_tgp = 0.49368;
              break;
      }
  
      const double pmass=0.93827;
      const double cmass=11.118;
      
      const double BeamEnergy=120;
      
      
      TLorentzVector C1;   //original C1 vector
      TLorentzVector p1;   //original proton vector
      
      //max energy is 149.1
      
      // proton momentum is 120 GeV/c, all in pZ
      
      double protonPz=BeamEnergy;
      double protonE=pmass*sqrt(1+(protonPz/pmass)*(protonPz/pmass));
      p1.SetPz(protonPz);
      p1.SetE(protonE);
      
      C1.SetE(cmass);
      
      TLorentzVector sumbefore=p1+C1;
      
      const double InvMass=sumbefore.M();           //////////// about 60.3 - sqrt( 2 * mass_carbon^2 + 2 * Energy_beam * mass_carbon)
      
      /////////////////// end of setting part.

      
      double Boost=protonPz/(protonE+cmass);
      
     // std::cout << "Boost " << Boost << " protonE "<< protonE << std::endl;
      
      p1.Print();
      C1.Print();
    
      p1.Boost(0,0,-Boost);
      C1.Boost(0,0,-Boost);
      
      p1.Print();
      C1.Print();
      
      TLorentzVector sumafter=p1+C1;
      
      //std::cout << " xf?? " << 2*p1.Pz()/17.3 << " " << p1.Dot(C1) << " " << sumafter.Mag() << " " << sumbefore.Mag() << " .m "<< sumafter.M() << " " << sumbefore.M() <<  std::endl;
  
    
      double Momentum_tgp = sqrt( pow(fluxNtuple -> tgppx ,2) + pow(fluxNtuple -> tgppy ,2) + pow(fluxNtuple -> tgppz ,2) );
      
      double pL_tgp = fluxNtuple -> tgppz;
      double pT_tgp = sqrt( pow(fluxNtuple -> tgppx ,2) + pow(fluxNtuple -> tgppy ,2));
      
      double Energy_tgp = sqrt( pow(Mass_tgp,2) + pow(Momentum_tgp,2));
      
      
              
      TLorentzVector FinalVector;
              
      FinalVector.SetPz(pL_tgp);
      FinalVector.SetE(Energy_tgp);
      
      FinalVector.Boost(0,0,-Boost);
      double Xf= 2*(FinalVector.Pz())/InvMass;    //////////////////////
      
    //  if(fluxNtuple->tgptype == 11 ||fluxNtuple->tgptype == 12 ||fluxNtuple->tgptype == 8 ||fluxNtuple->tgptype == 9 ){
      //    std::cout << "tgptype: "<<fluxNtuple->tgptype<<" E: " << Energy_tgp << " xF: " << Xf << " alt method " << calculateXf(InvMass,Mass_tgp,pL_tgp,pT_tgp,fixedTargetBoost(BeamEnergy,pmass,cmass)) <<  std::endl;
    
              
          //}
      
      
      ////// MIPP (Seun's) /////////////////////////////////////////
      
      double momentum_L = fluxNtuple -> tgppz;
      
      double momentum_T = sqrt( pow(fluxNtuple->tgppx,2) + pow(fluxNtuple->tgppy , 2) );
    
      

      switch( fluxNtuple -> tgptype)
      {
          case 8: //pion plus
              pionplus_MIPP->Fill(momentum_L , momentum_T);
              
              if (fluxNtuple -> Ntype == nue){
                  pionplus_MIPP_nue -> Fill(momentum_L , momentum_T);
                  nue_tptype_pionplus->Fill(enu,weight);
                  
                  if (fluxNtuple -> Ndecay == 6) //K^{+} -> #nu_{e} #pi^{0} e^{+}
                  {pionplus_MIPP_nue_dm6->Fill(momentum_L, momentum_T);}
                  if (fluxNtuple -> Ndecay == 1)// K^{0}_{L} -> #nu_{e} #pi^{-} e^{+}
                  {pionplus_MIPP_nue_dm1->Fill(momentum_L, momentum_T);}
                  if (fluxNtuple -> Ndecay == 11) //#mu^{+} -> #bar{#nu}_{#mu} #nu_{e} e^{+}
                  {pionplus_MIPP_nue_dm11->Fill(momentum_L, momentum_T);}
                  if (fluxNtuple -> Ndecay == 5)// K^{+} -> #nu_{#mu} #mu^{+}
                  {pionplus_MIPP_nue_dm5->Fill(momentum_L, momentum_T);}
              
              
              }
                          
              if (fluxNtuple -> Ntype == numu)
              { pionplus_MIPP_numu -> Fill(momentum_L , momentum_T);}
              if (fluxNtuple -> Ntype == anue)
              { pionplus_MIPP_anue -> Fill(momentum_L , momentum_T);}
              if (fluxNtuple -> Ntype == anumu)
              { pionplus_MIPP_anumu -> Fill(momentum_L , momentum_T);}
              break;
              
          case 9: //pion minus
              pionminus_MIPP->Fill(momentum_L , momentum_T);
              if (fluxNtuple -> Ntype == nue)
              { pionminus_MIPP_nue -> Fill(momentum_L , momentum_T);}
              if (fluxNtuple -> Ntype == numu)
              { pionminus_MIPP_numu -> Fill(momentum_L , momentum_T);}
              if (fluxNtuple -> Ntype == anue)
              { pionminus_MIPP_anue -> Fill(momentum_L , momentum_T);}
              if (fluxNtuple -> Ntype == anumu)
              { pionminus_MIPP_anumu -> Fill(momentum_L , momentum_T);}
              break;
              
          case 11: // k plus
              Kplus_MIPP->Fill(momentum_L , momentum_T);
              if (fluxNtuple -> Ntype == nue)
              { Kplus_MIPP_nue -> Fill(momentum_L , momentum_T);}
              if (fluxNtuple -> Ntype == numu)
              { Kplus_MIPP_numu -> Fill(momentum_L , momentum_T);}
              if (fluxNtuple -> Ntype == anue)
              { Kplus_MIPP_anue -> Fill(momentum_L , momentum_T);}
              if (fluxNtuple -> Ntype == anumu)
              { Kplus_MIPP_anumu -> Fill(momentum_L , momentum_T);}
              break;
              
          case 12: // k minus
              Kminus_MIPP->Fill(momentum_L , momentum_T);
              if (fluxNtuple -> Ntype == nue)
              { Kminus_MIPP_nue -> Fill(momentum_L , momentum_T);
                  if(fluxNtuple->Ndecay == 1){
                      Kminus_MIPP_nue_dm_one -> Fill(momentum_L , momentum_T);
                  }
                  if(fluxNtuple->Ndecay == 6){
                      Kminus_MIPP_nue_dm_six -> Fill(momentum_L , momentum_T);
                  }
                  if(fluxNtuple->Ndecay == 11){
                      Kminus_MIPP_nue_dm_eleven -> Fill(momentum_L , momentum_T);
                  }
              }
              if (fluxNtuple -> Ntype == numu)
              { Kminus_MIPP_numu -> Fill(momentum_L , momentum_T);}
              if (fluxNtuple -> Ntype == anue)
              { Kminus_MIPP_anue -> Fill(momentum_L , momentum_T);}
              if (fluxNtuple -> Ntype == anumu)
              { Kminus_MIPP_anumu -> Fill(momentum_L , momentum_T);}
              break;
      }
      
      //////////////////////////////////////////////////////////////
      
      
      
      switch ( fluxNtuple-> tgptype ){
          case 8: // for pion plus
              
              //std::cout<<"case 8 (pi plus): "<<" x_F "<<Xf<<" pT_tgp: "<<pT_tgp<<"\n"<<endl;
              pionplusAngleHisto->Fill(angle);
              pionplusMomentumHisto->Fill(momentum);
              pionplus2DHisto->Fill(momentum,angle);
              pionplus_NA49->Fill(Xf , pT_tgp);
              
              // if(fluxNtuple->ptype == 53)
             // {
              //    pionplus_nue->Fill(momentum, angle);
             // }
              
              if(fluxNtuple->Ntype == 53)
            {
                pionplus_nue->Fill(momentum, angle);
                pionplus_NA49_nue -> Fill(Xf ,pT_tgp);
            }
              if(fluxNtuple->Ntype == 56)
          {
              pionplus_numu->Fill(momentum, angle);
          }
              if(fluxNtuple->Ntype == 55)
              {
                  pionplus_anumu->Fill(momentum, angle);
              }
              if(fluxNtuple->Ntype == 52)
              {
                  pionplus_anue->Fill(momentum, angle);
              }
              
              if(fluxNtuple->Ntype == 56)
              {
                  switch (fluxNtuple-> Ndecay)
                  {
                      case decaymode_one:
                          ppnumu_dm_one->Fill(momentum,angle);
                          break;
                          
                      case decaymode_two:
                          ppnumu_dm_two->Fill(momentum,angle);
                          break;
                          
                      case decaymode_three:
                          ppnumu_dm_three->Fill(momentum,angle);
                          break;
                          
                      case decaymode_four:
                          ppnumu_dm_four->Fill(momentum,angle);
                          break;
                          
                      case decaymode_five:
                          ppnumu_dm_five->Fill(momentum,angle);
                          break;
                          
                      case decaymode_six:
                          ppnumu_dm_six->Fill(momentum,angle);
                          break;
                          
                      case decaymode_seven:
                          ppnumu_dm_seven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eight:
                          ppnumu_dm_eight->Fill(momentum,angle);
                          break;
                          
                      case decaymode_nine:
                          ppnumu_dm_nine->Fill(momentum,angle);
                          break;
                          
                      case decaymode_ten:
                          ppnumu_dm_ten->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eleven:
                          ppnumu_dm_eleven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_twelve:
                          ppnumu_dm_twelve->Fill(momentum,angle);
                          break;
                          
                      case decaymode_thirteen:
                          ppnumu_dm_thirteen->Fill(momentum,angle);
                          break;
                          
                      case decaymode_fourteen:
                          ppnumu_dm_fourteen->Fill(momentum,angle);
                          break;
                  }
              }
              
              if(fluxNtuple->Ntype == 53)
              {
                  switch (fluxNtuple-> Ndecay)
                  {
                      case decaymode_one:
                          ppnue_dm_one->Fill(momentum,angle);
                          break;
                          
                      case decaymode_two:
                          ppnue_dm_two->Fill(momentum,angle);
                          break;
                      
                      case decaymode_three:
                          ppnue_dm_three->Fill(momentum,angle);
                          break;
                          
                      case decaymode_four:
                          ppnue_dm_four->Fill(momentum,angle);
                          break;
                          
                      case decaymode_five:
                          ppnue_dm_five->Fill(momentum,angle);
                          break;
                          
                      case decaymode_six:
                          ppnue_dm_six->Fill(momentum,angle);
                          break;
                          
                      case decaymode_seven:
                          ppnue_dm_seven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eight:
                          ppnue_dm_eight->Fill(momentum,angle);
                          break;
                          
                      case decaymode_nine:
                          ppnue_dm_nine->Fill(momentum,angle);
                          break;
                          
                      case decaymode_ten:
                          ppnue_dm_ten->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eleven:
                          ppnue_dm_eleven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_twelve:
                          ppnue_dm_twelve->Fill(momentum,angle);
                          break;
                          
                      case decaymode_thirteen:
                          ppnue_dm_thirteen->Fill(momentum,angle);
                          break;
                          
                      case decaymode_fourteen:
                          ppnue_dm_fourteen->Fill(momentum,angle);
                          break;
                  }
                  
              }
              break;
              
          case 9: // for pion minus
              //std::cout<<"case 9 (pi minus): "<<" x_F "<<Xf<<" pT_tgp: "<<pT_tgp<<"\n"<<endl;
              
              pionminusAngleHisto->Fill(angle);
              pionminusMomentumHisto->Fill(momentum);
              pionminus2DHisto->Fill(momentum,angle);
              pionminus_NA49->Fill(Xf , pT_tgp);
             // pionminus_NA49->Fill(mc_flux_pz , momentum_transverse);
              if(fluxNtuple->Ntype == 53)
          {
              pionminus_nue->Fill(momentum, angle);
              pionminus_NA49_nue -> Fill(Xf ,pT_tgp);
          }
              if(fluxNtuple->Ntype == 56)
          {
              pionminus_numu->Fill(momentum, angle);
          }
              if(fluxNtuple->Ntype == 55)
              {
                  pionminus_anumu->Fill(momentum, angle);
              }
              if(fluxNtuple->Ntype == 52)
              {
                  pionminus_anue->Fill(momentum, angle);
              }
              
              if(fluxNtuple->Ntype == 53)
              {
                  switch (fluxNtuple-> Ndecay)
                  {
                      case decaymode_one:
                          pmnue_dm_one->Fill(momentum,angle);
                          break;
                          
                      case decaymode_two:
                          pmnue_dm_two->Fill(momentum,angle);
                          break;
                          
                      case decaymode_three:
                          pmnue_dm_three->Fill(momentum,angle);
                          break;
                          
                      case decaymode_four:
                          pmnue_dm_four->Fill(momentum,angle);
                          break;
                          
                      case decaymode_five:
                          pmnue_dm_five->Fill(momentum,angle);
                          break;
                          
                      case decaymode_six:
                          pmnue_dm_six->Fill(momentum,angle);
                          break;
                          
                      case decaymode_seven:
                          pmnue_dm_seven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eight:
                          pmnue_dm_eight->Fill(momentum,angle);
                          break;
                          
                      case decaymode_nine:
                          pmnue_dm_nine->Fill(momentum,angle);
                          break;
                          
                      case decaymode_ten:
                          pmnue_dm_ten->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eleven:
                          pmnue_dm_eleven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_twelve:
                          pmnue_dm_twelve->Fill(momentum,angle);
                          break;
                          
                      case decaymode_thirteen:
                          pmnue_dm_thirteen->Fill(momentum,angle);
                          break;
                          
                      case decaymode_fourteen:
                          pmnue_dm_fourteen->Fill(momentum,angle);
                          break;
                  }
                  
              }
              
              if(fluxNtuple->Ntype == 56)
              {
                  switch (fluxNtuple-> Ndecay)
                  {
                      case decaymode_one:
                          pmnumu_dm_one->Fill(momentum,angle);
                          break;
                          
                      case decaymode_two:
                          pmnumu_dm_two->Fill(momentum,angle);
                          break;
                          
                      case decaymode_three:
                          pmnumu_dm_three->Fill(momentum,angle);
                          break;
                          
                      case decaymode_four:
                          pmnumu_dm_four->Fill(momentum,angle);
                          break;
                          
                      case decaymode_five:
                          pmnumu_dm_five->Fill(momentum,angle);
                          break;
                          
                      case decaymode_six:
                          pmnumu_dm_six->Fill(momentum,angle);
                          break;
                          
                      case decaymode_seven:
                          pmnumu_dm_seven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eight:
                          pmnumu_dm_eight->Fill(momentum,angle);
                          break;
                          
                      case decaymode_nine:
                          pmnumu_dm_nine->Fill(momentum,angle);
                          break;
                          
                      case decaymode_ten:
                          pmnumu_dm_ten->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eleven:
                          pmnumu_dm_eleven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_twelve:
                          pmnumu_dm_twelve->Fill(momentum,angle);
                          break;
                          
                      case decaymode_thirteen:
                          pmnumu_dm_thirteen->Fill(momentum,angle);
                          break;
                          
                      case decaymode_fourteen:
                          pmnumu_dm_fourteen->Fill(momentum,angle);
                          break;
                  }
              }
              
              break;
              
          case 10: // for K long
             // std::cout<<"case 10 (K long): "<<" angle: "<<angle<<" momentum: "<<momentum<<"\n"<<endl;
              KlongAngleHisto->Fill(angle);
              KlongMomentumHisto->Fill(momentum);
              Klong2DHisto->Fill(momentum,angle);
             // Klong_NA49->Fill(mc_flux_pz , momentum_transverse);
              
              if(fluxNtuple->Ntype == 53)
          {
              Klong_nue->Fill(momentum, angle);
          }
              if(fluxNtuple->Ntype == 56)
          {
              Klong_numu->Fill(momentum, angle);
          }
              if(fluxNtuple->Ntype == 55)
              {
                  Klong_anumu->Fill(momentum, angle);
              }
              if(fluxNtuple->Ntype == 52)
              {
                  Klong_anue->Fill(momentum, angle);
              }
              
              if(fluxNtuple->Ntype == 53)
              {
                  switch (fluxNtuple-> Ndecay)
                  {
                      case decaymode_one:
                          klnue_dm_one->Fill(momentum,angle);
                          break;
                          
                      case decaymode_two:
                          klnue_dm_two->Fill(momentum,angle);
                          break;
                          
                      case decaymode_three:
                          klnue_dm_three->Fill(momentum,angle);
                          break;
                          
                      case decaymode_four:
                          klnue_dm_four->Fill(momentum,angle);
                          break;
                          
                      case decaymode_five:
                          klnue_dm_five->Fill(momentum,angle);
                          break;
                          
                      case decaymode_six:
                          klnue_dm_six->Fill(momentum,angle);
                          break;
                          
                      case decaymode_seven:
                          klnue_dm_seven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eight:
                          klnue_dm_eight->Fill(momentum,angle);
                          break;
                          
                      case decaymode_nine:
                          klnue_dm_nine->Fill(momentum,angle);
                          break;
                          
                      case decaymode_ten:
                          klnue_dm_ten->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eleven:
                          klnue_dm_eleven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_twelve:
                          klnue_dm_twelve->Fill(momentum,angle);
                          break;
                          
                      case decaymode_thirteen:
                          klnue_dm_thirteen->Fill(momentum,angle);
                          break;
                          
                      case decaymode_fourteen:
                          klnue_dm_fourteen->Fill(momentum,angle);
                          break;
                  }
                  
              }
              
              if(fluxNtuple->Ntype == 56)
              {
                  switch (fluxNtuple-> Ndecay)
                  {
                      case decaymode_one:
                          klnumu_dm_one->Fill(momentum,angle);
                          break;
                          
                      case decaymode_two:
                          klnumu_dm_two->Fill(momentum,angle);
                          break;
                          
                      case decaymode_three:
                          klnumu_dm_three->Fill(momentum,angle);
                          break;
                          
                      case decaymode_four:
                          klnumu_dm_four->Fill(momentum,angle);
                          break;
                          
                      case decaymode_five:
                          klnumu_dm_five->Fill(momentum,angle);
                          break;
                          
                      case decaymode_six:
                          klnumu_dm_six->Fill(momentum,angle);
                          break;
                          
                      case decaymode_seven:
                          klnumu_dm_seven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eight:
                          klnumu_dm_eight->Fill(momentum,angle);
                          break;
                          
                      case decaymode_nine:
                          klnumu_dm_nine->Fill(momentum,angle);
                          break;
                          
                      case decaymode_ten:
                          klnumu_dm_ten->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eleven:
                          klnumu_dm_eleven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_twelve:
                          klnumu_dm_twelve->Fill(momentum,angle);
                          break;
                          
                      case decaymode_thirteen:
                          klnumu_dm_thirteen->Fill(momentum,angle);
                          break;
                          
                      case decaymode_fourteen:
                          klnumu_dm_fourteen->Fill(momentum,angle);
                          break;
                  }
              }
              break;
              
          case 11: // for K plus
              //std::cout<<"case 11 (K plus): "<<" x_F "<<Xf<<" pT_tgp: "<<pT_tgp<<"\n"<<endl;
              KplusAngleHisto->Fill(angle);
              KplusMomentumHisto->Fill(momentum);
              Kplus2DHisto->Fill(momentum,angle);
              Kplus_NA49->Fill(Xf , pT_tgp);
              Kplus_NA49_second->Fill(calculateXf(15,0.49368 ,pL_tgp, pT_tgp, fixedTargetBoost(120, pmass,cmass)) , pT_tgp);
              
            //  Kplus_NA49->Fill(mc_flux_pz , momentum_transverse);
              
              if(fluxNtuple->Ntype == 53)
          {
              Kplus_nue->Fill(momentum, angle);
              Kplus_NA49_nue -> Fill(Xf ,pT_tgp);
          }
              if(fluxNtuple->Ntype == 56)
          {
              Kplus_numu->Fill(momentum, angle);
          }
              if(fluxNtuple->Ntype == 55)
              {
                  Kplus_anumu->Fill(momentum, angle);
              }
              if(fluxNtuple->Ntype == 52)
              {
                  Kplus_anue->Fill(momentum, angle);
              }
              
              if(fluxNtuple->Ntype == 53)
              {
                  switch (fluxNtuple-> Ndecay)
                  {
                      case decaymode_one:
                          kpnue_dm_one->Fill(momentum,angle);
                          break;
                          
                      case decaymode_two:
                          kpnue_dm_two->Fill(momentum,angle);
                          break;
                          
                      case decaymode_three:
                          kpnue_dm_three->Fill(momentum,angle);
                          break;
                          
                      case decaymode_four:
                          kpnue_dm_four->Fill(momentum,angle);
                          break;
                          
                      case decaymode_five:
                          kpnue_dm_five->Fill(momentum,angle);
                          break;
                          
                      case decaymode_six:
                          kpnue_dm_six->Fill(momentum,angle);
                          break;
                          
                      case decaymode_seven:
                          kpnue_dm_seven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eight:
                          kpnue_dm_eight->Fill(momentum,angle);
                          break;
                          
                      case decaymode_nine:
                          kpnue_dm_nine->Fill(momentum,angle);

                          break;
                          
                      case decaymode_ten:
                          kpnue_dm_ten->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eleven:
                          kpnue_dm_eleven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_twelve:
                          kpnue_dm_twelve->Fill(momentum,angle);
                          break;
                          
                      case decaymode_thirteen:
                          kpnue_dm_thirteen->Fill(momentum,angle);
                          break;
                          
                      case decaymode_fourteen:
                          kpnue_dm_fourteen->Fill(momentum,angle);
                          break;
                  }
                  
              }
              
              if(fluxNtuple->Ntype == 56)
              {
                  switch (fluxNtuple-> Ndecay)
                  {
                      case decaymode_one:
                          kpnumu_dm_one->Fill(momentum,angle);
                          break;
                          
                      case decaymode_two:
                          kpnumu_dm_two->Fill(momentum,angle);
                          break;
                          
                      case decaymode_three:
                          kpnumu_dm_three->Fill(momentum,angle);
                          break;
                          
                      case decaymode_four:
                          kpnumu_dm_four->Fill(momentum,angle);
                          break;
                          
                      case decaymode_five:
                          kpnumu_dm_five->Fill(momentum,angle);
                          break;
                          
                      case decaymode_six:
                          kpnumu_dm_six->Fill(momentum,angle);
                          break;
                          
                      case decaymode_seven:
                          kpnumu_dm_seven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eight:
                          kpnumu_dm_eight->Fill(momentum,angle);
                          break;
                          
                      case decaymode_nine:
                          kpnumu_dm_nine->Fill(momentum,angle);
                          break;
                          
                      case decaymode_ten:
                          kpnumu_dm_ten->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eleven:
                          kpnumu_dm_eleven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_twelve:
                          kpnumu_dm_twelve->Fill(momentum,angle);
                          break;
                          
                      case decaymode_thirteen:
                          kpnumu_dm_thirteen->Fill(momentum,angle);
                          break;
                          
                      case decaymode_fourteen:
                          kpnumu_dm_fourteen->Fill(momentum,angle);
                          break;
                  }
              }
              break;
              
          case 12: // for K minus
              //std::cout<<"case 12 (K minus): "<<" x_F "<<Xf<<" pT_tgp: "<<pT_tgp<<"\n"<<endl;
              KminusAngleHisto->Fill(angle);
              KminusMomentumHisto->Fill(momentum);
              Kminus2DHisto->Fill(momentum,angle);
              Kminus_NA49->Fill(Xf , pT_tgp);
            //  Kminus_NA49->Fill(mc_flux_pz , momentum_L);
              
              
              if(fluxNtuple->Ntype == 53)
          {
              Kminus_nue->Fill(momentum, angle);
              Kminus_NA49_nue -> Fill(Xf ,pT_tgp);
          }
              if(fluxNtuple->Ntype == 56)
          {
              Kminus_numu->Fill(momentum, angle);
          }
              if(fluxNtuple->Ntype == 55)
              {
                  Kminus_anumu->Fill(momentum, angle);
              }
              if(fluxNtuple->Ntype == 52)
              {
                  Kminus_anue->Fill(momentum, angle);
              }
              
              if(fluxNtuple->Ntype == 53)
              {
                  switch (fluxNtuple-> Ndecay)
                  {
                      case decaymode_one:
                          kmnue_dm_one->Fill(momentum,angle);
                          break;
                          
                      case decaymode_two:
                          kmnue_dm_two->Fill(momentum,angle);
                          break;
                          
                      case decaymode_three:
                          kmnue_dm_three->Fill(momentum,angle);
                          break;
                          
                      case decaymode_four:
                          kmnue_dm_four->Fill(momentum,angle);
                          break;
                          
                      case decaymode_five:
                          kmnue_dm_five->Fill(momentum,angle);
                          break;
                          
                      case decaymode_six:
                          kmnue_dm_six->Fill(momentum,angle);
                          break;
                          
                      case decaymode_seven:
                          kmnue_dm_seven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eight:
                          kmnue_dm_eight->Fill(momentum,angle);
                          break;
                          
                      case decaymode_nine:
                          kmnue_dm_nine->Fill(momentum,angle);
                          break;
                          
                      case decaymode_ten:
                          kmnue_dm_ten->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eleven:
                          kmnue_dm_eleven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_twelve:
                          kmnue_dm_twelve->Fill(momentum,angle);
                          break;
                          
                      case decaymode_thirteen:
                          kmnue_dm_thirteen->Fill(momentum,angle);
                          break;
                          
                      case decaymode_fourteen:
                          kmnue_dm_fourteen->Fill(momentum,angle);
                          break;
                  }
                  
              }
              
              if(fluxNtuple->Ntype == 56)
              {
                  switch (fluxNtuple-> Ndecay)
                  {
                      case decaymode_one:
                          kmnumu_dm_one->Fill(momentum,angle);
                          break;
                          
                      case decaymode_two:
                          kmnumu_dm_two->Fill(momentum,angle);
                          break;
                          
                      case decaymode_three:
                          kmnumu_dm_three->Fill(momentum,angle);
                          break;
                          
                      case decaymode_four:
                          kmnumu_dm_four->Fill(momentum,angle);
                          break;
                          
                      case decaymode_five:
                          kmnumu_dm_five->Fill(momentum,angle);
                          break;
                          
                      case decaymode_six:
                          kmnumu_dm_six->Fill(momentum,angle);
                          break;
                          
                      case decaymode_seven:
                          kmnumu_dm_seven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eight:
                          kmnumu_dm_eight->Fill(momentum,angle);
                          break;
                          
                      case decaymode_nine:
                          kmnumu_dm_nine->Fill(momentum,angle);
                          break;
                          
                      case decaymode_ten:
                          kmnumu_dm_ten->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eleven:
                          kmnumu_dm_eleven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_twelve:
                          kmnumu_dm_twelve->Fill(momentum,angle);
                          break;
                          
                      case decaymode_thirteen:
                          kmnumu_dm_thirteen->Fill(momentum,angle);
                          break;
                          
                      case decaymode_fourteen:
                          kmnumu_dm_fourteen->Fill(momentum,angle);
                          break;
                  }
              }
              break;
              
          case 5: // for mu plus
              //std::cout<<"case 5 (mu plus): "<<" angle: "<<angle<<" momentum: "<<momentum<<"\n"<<endl;
              muonplusAngleHisto->Fill(angle);
              muonplusMomentumHisto->Fill(momentum);
              muonplus2DHisto->Fill(momentum,angle);
              //muonplusNA49->Fill(momentum , momentum_T);
              break;
              
          case 6: // for mu minus
             // std::cout<<"case 6 (mu minus): "<<" angle: "<<angle<<" momentum: "<<momentum<<"\n"<<endl;
              muonminusAngleHisto->Fill(angle);
              muonminusMomentumHisto->Fill(momentum);
              muonminus2DHisto->Fill(momentum,angle);
              break;
              
          case 13: // for neutron
             // std::cout<<"case 13 (neutron): "<<" angle: "<<angle<<" momentum: "<<momentum<<"\n"<<endl;
              neutronAngleHisto->Fill(angle);
              neutronMomentumHisto->Fill(momentum);
              neutron2DHisto->Fill(momentum,angle);
              neutron_NA49->Fill(Xf , pT_tgp);
              
              if(fluxNtuple->Ntype == 53)
          {
              neutron_nue->Fill(momentum, angle);
              neutron_NA49_nue->Fill(Xf , pT_tgp);
          }
              if(fluxNtuple->Ntype == 56)
          {
              neutron_numu->Fill(momentum, angle);
          }
              if(fluxNtuple->Ntype == 55)
              {
                  neutron_anumu->Fill(momentum, angle);
              }
              if(fluxNtuple->Ntype == 52)
              {
                  neutron_anue->Fill(momentum, angle);
              }
              if(fluxNtuple->Ntype == 53)
              {
                  switch (fluxNtuple-> Ndecay)
                  {
                      case decaymode_one:
                          nnue_dm_one->Fill(momentum,angle);
                          break;
                          
                      case decaymode_two:
                          nnue_dm_two->Fill(momentum,angle);
                          break;
                          
                      case decaymode_three:
                          nnue_dm_three->Fill(momentum,angle);
                          break;
                          
                      case decaymode_four:
                          nnue_dm_four->Fill(momentum,angle);
                          break;
                          
                      case decaymode_five:
                          nnue_dm_five->Fill(momentum,angle);
                          break;
                          
                      case decaymode_six:
                          nnue_dm_six->Fill(momentum,angle);
                          break;
                          
                      case decaymode_seven:
                          nnue_dm_seven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eight:
                          nnue_dm_eight->Fill(momentum,angle);
                          break;
                          
                      case decaymode_nine:
                          nnue_dm_nine->Fill(momentum,angle);
                          break;
                          
                      case decaymode_ten:
                          nnue_dm_ten->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eleven:
                          nnue_dm_eleven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_twelve:
                          nnue_dm_twelve->Fill(momentum,angle);
                          break;
                          
                      case decaymode_thirteen:
                          nnue_dm_thirteen->Fill(momentum,angle);
                          break;
                          
                      case decaymode_fourteen:
                          nnue_dm_fourteen->Fill(momentum,angle);
                          break;
                  }
                  
              }
              
              if(fluxNtuple->Ntype == 56)
              {
                  switch (fluxNtuple-> Ndecay)
                  {
                      case decaymode_one:
                          nnumu_dm_one->Fill(momentum,angle);
                          break;
                          
                      case decaymode_two:
                          nnumu_dm_two->Fill(momentum,angle);
                          break;
                          
                      case decaymode_three:
                          nnumu_dm_three->Fill(momentum,angle);
                          break;
                          
                      case decaymode_four:
                          nnumu_dm_four->Fill(momentum,angle);
                          break;
                          
                      case decaymode_five:
                          nnumu_dm_five->Fill(momentum,angle);
                          break;
                          
                      case decaymode_six:
                          nnumu_dm_six->Fill(momentum,angle);
                          break;
                          
                      case decaymode_seven:
                          nnumu_dm_seven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eight:
                          nnumu_dm_eight->Fill(momentum,angle);
                          break;
                          
                      case decaymode_nine:
                          nnumu_dm_nine->Fill(momentum,angle);
                          break;
                          
                      case decaymode_ten:
                          nnumu_dm_ten->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eleven:
                          nnumu_dm_eleven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_twelve:
                          nnumu_dm_twelve->Fill(momentum,angle);
                          break;
                          
                      case decaymode_thirteen:
                          nnumu_dm_thirteen->Fill(momentum,angle);
                          break;
                          
                      case decaymode_fourteen:
                          nnumu_dm_fourteen->Fill(momentum,angle);
                          break;
                  }
              }
              break;
              
            
              
          case 14: // for proton
           //   std::cout<<"case 14 (proton): "<<" angle: "<<angle<<" momentum: "<<momentum<<"\n"<<endl;
              protonAngleHisto->Fill(angle);
              protonMomentumHisto->Fill(momentum);
              proton2DHisto->Fill(momentum,angle);
          
              proton_NA49->Fill(Xf , pT_tgp);
              
              if(fluxNtuple->Ntype == 53)
          {
              proton_nue->Fill(momentum, angle);
              proton_NA49_nue -> Fill(Xf ,pT_tgp);
          }
              if(fluxNtuple->Ntype == 56)
          {
              proton_numu->Fill(momentum, angle);
          }
              
              if(fluxNtuple->Ntype == 55)
              {
                  proton_anumu->Fill(momentum, angle);
              }
              if(fluxNtuple->Ntype == 52)
              {
                  proton_anue->Fill(momentum, angle);
              }
              if(fluxNtuple->Ntype == 53)
              {
                  switch (fluxNtuple-> Ndecay)
                  {
                      case decaymode_one:
                          pnue_dm_one->Fill(momentum,angle);
                          break;
                          
                      case decaymode_two:
                          pnue_dm_two->Fill(momentum,angle);
                          break;
                          
                      case decaymode_three:
                          pnue_dm_three->Fill(momentum,angle);
                          break;
                          
                      case decaymode_four:
                          pnue_dm_four->Fill(momentum,angle);
                          break;
                          
                      case decaymode_five:
                          pnue_dm_five->Fill(momentum,angle);
                          break;
                          
                      case decaymode_six:
                          pnue_dm_six->Fill(momentum,angle);
                          break;
                          
                      case decaymode_seven:
                          pnue_dm_seven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eight:
                          pnue_dm_eight->Fill(momentum,angle);
                          break;
                          
                      case decaymode_nine:
                          pnue_dm_nine->Fill(momentum,angle);
                          break;
                          
                      case decaymode_ten:
                          pnue_dm_ten->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eleven:
                          pnue_dm_eleven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_twelve:
                          pnue_dm_twelve->Fill(momentum,angle);
                          break;
                          
                      case decaymode_thirteen:
                          pnue_dm_thirteen->Fill(momentum,angle);
                          break;
                          
                      case decaymode_fourteen:
                          pnue_dm_fourteen->Fill(momentum,angle);
                          break;
                  }
                  
              }
              if(fluxNtuple->Ntype == 56)
              {
                  switch (fluxNtuple-> Ndecay)
                  {
                      case decaymode_one:
                          pnumu_dm_one->Fill(momentum,angle);
                          break;
                          
                      case decaymode_two:
                          pnumu_dm_two->Fill(momentum,angle);
                          break;
                          
                      case decaymode_three:
                          pnumu_dm_three->Fill(momentum,angle);
                          break;
                          
                      case decaymode_four:
                          pnumu_dm_four->Fill(momentum,angle);
                          break;
                          
                      case decaymode_five:
                          pnumu_dm_five->Fill(momentum,angle);
                          break;
                          
                      case decaymode_six:
                          pnumu_dm_six->Fill(momentum,angle);
                          break;
                          
                      case decaymode_seven:
                          pnumu_dm_seven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eight:
                          pnumu_dm_eight->Fill(momentum,angle);
                          break;
                          
                      case decaymode_nine:
                          pnumu_dm_nine->Fill(momentum,angle);
                          break;
                          
                      case decaymode_ten:
                          pnumu_dm_ten->Fill(momentum,angle);
                          break;
                          
                      case decaymode_eleven:
                          pnumu_dm_eleven->Fill(momentum,angle);
                          break;
                          
                      case decaymode_twelve:
                          pnumu_dm_twelve->Fill(momentum,angle);
                          break;
                          
                      case decaymode_thirteen:
                          pnumu_dm_thirteen->Fill(momentum,angle);
                          break;
                          
                      case decaymode_fourteen:
                          pnumu_dm_fourteen->Fill(momentum,angle);
                          break;
                  }
              }
              
              break;
              
          case 18: // for lambda
           //   std::cout<<"case 18 (lambda): "<<" angle: "<<angle<<" momentum: "<<momentum<<"\n"<<endl;
              lambdaAngleHisto->Fill(angle);
              lambdaMomentumHisto->Fill(momentum);
              lambda2DHisto->Fill(momentum,angle);
             // lambda_NA49->Fill(mc_flux_pz , momentum_transverse);
              
              
          case 19: // for sigma plus
            //  std::cout<<"case 19 (sigmaplus): "<<" angle: "<<angle<<" momentum: "<<momentum<<"\n"<<endl;
              sigmaplusAngleHisto->Fill(angle);
              sigmaplusMomentumHisto->Fill(momentum);
              sigmaplus2DHisto->Fill(momentum,angle);
           //   sigmaplus_NA49->Fill(mc_flux_pz , momentum_transverse);
      }
      
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Fill the histograms
    switch (fluxNtuple->Ntype) {
      case numu:
        numuFluxHisto->Fill(enu, weight);
            //std::cout<< "Normal: "<<"enu: "<<enu<<" weight: "<<weight<<std::endl;
            if (fluxNtuple -> tptype==11 || fluxNtuple -> tptype==12 ){ //Kplus or Kminus
                if (pZ < 20 || (pZ < 24 && pT > 1) || (pZ<31 && pT>1.2) || (pZ<42 && pT>1.55) || pZ>90 || pT>2 ){
               //     std::cout<<"tptype: "<<fluxNtuple->tptype << " pZ: " <<pZ<<" pT: "<<pT<<std::endl;
                //    std::cout<<"enu: "<<enu<<" weight: "<<weight<<std::endl;
                    numuFluxHisto_up -> Fill(enu, weight*1.75);
                    numuFluxHisto_down -> Fill(enu, weight*0.75);
                }
            }
            else {
                numuFluxHisto_up -> Fill(enu, weight);
                numuFluxHisto_down -> Fill(enu, weight);
            }
            
          //  numuFluxHisto_up->Fill(enu,weight);
            
        if (fluxNtuple->ptype == 8)                 // numu_parent
        {numu_parent_pionplus->Fill(enu,weight);
        }
        if (fluxNtuple->ptype == 11)
        {numu_parent_kaonplus->Fill(enu,weight);
        }
        if (fluxNtuple->ptype == 10)
        {numu_parent_kaonlong->Fill(enu,weight);
        }
        if (fluxNtuple->ptype == 6)  //muon minus
        {numu_parent_muonminus->Fill(enu,weight);
        }
        if (fluxNtuple->ptype == 5)  //muon plus
        {numu_parent_muonplus->Fill(enu,weight);
        }
            
        // plot Ndecays 13, 5, 7, 3, 12
        
        if (fluxNtuple->Ndecay==13)
        {
            numu_Ndecay_thirteen->Fill(enu,weight);
        }
        
        if (fluxNtuple->Ndecay==5)
        {
            numu_Ndecay_five->Fill(enu,weight);
        }
        
        if (fluxNtuple->Ndecay==7)
        {
            numu_Ndecay_seven->Fill(enu,weight);
        }
        
        if (fluxNtuple->Ndecay==12)
        {
            numu_Ndecay_twelve->Fill(enu,weight);
        }
        
                // plot tgptype s
    
        if (fluxNtuple -> tgptype == 8)
        {numu_tgptype_pionplus -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 9)
        {numu_tgptype_pionminus -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 10)
        {numu_tgptype_kaonlong -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 11)
        {numu_tgptype_kaonplus -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 12)
        {numu_tgptype_kaonminus -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 13)
        {numu_tgptype_neutron -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 14)
        {numu_tgptype_proton -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 18)
        {numu_tgptype_lambda -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 19)
        {numu_tgptype_sigmaplus -> Fill(enu,weight);
        }
            
            // plot tgen s
        
        if (fluxNtuple -> tgen ==2)
        {numu_tgen_two -> Fill(enu,weight);
        }
        if (fluxNtuple -> tgen ==3)
        {numu_tgen_three -> Fill(enu,weight);
        }
        if (fluxNtuple -> tgen ==4)
        {numu_tgen_four -> Fill(enu,weight);
        }
        if (fluxNtuple -> tgen ==5)
        {numu_tgen_five -> Fill(enu,weight);
        }
        if (fluxNtuple -> tgen ==6)
        {numu_tgen_six -> Fill(enu,weight);
        }
        if (fluxNtuple -> tgen ==7)
        {numu_tgen_seven -> Fill(enu,weight);
        }
        if (fluxNtuple -> tgen ==8)
        {numu_tgen_eight -> Fill(enu,weight);
        }
        if (fluxNtuple -> tgen ==9)
        {numu_tgen_nine -> Fill(enu,weight);
        }
        if (fluxNtuple -> tgen ==10)
        {numu_tgen_ten -> Fill(enu,weight);
        }
        if (fluxNtuple -> tgen ==11)
        {numu_tgen_eleven -> Fill(enu,weight);
        }
        if (fluxNtuple -> tgen ==12)
        {numu_tgen_twelve -> Fill(enu,weight);
        }
            
        
        break;
            
      case anumu:
        anumuFluxHisto->Fill(enu, weight);
            
            if (fluxNtuple -> tptype==11 || fluxNtuple -> tptype==12 ){ //Kplus or Kminus
                if (pZ < 20 || (pZ < 24 && pT > 1) || (pZ<31 && pT>1.2) || (pZ<42 && pT>1.55) || pZ>90 || pT>2 ){
                    
                    anumuFluxHisto_up -> Fill(enu, weight*1.75);
                    anumuFluxHisto_down -> Fill(enu, weight*0.75);
                }
            }
            else {
                anumuFluxHisto_up -> Fill(enu, weight);
                anumuFluxHisto_down -> Fill(enu, weight);
            }
            
            
            if (fluxNtuple->ptype == 11)
            {anumu_parent_kaonplus->Fill(enu,weight);
            }
        break;
            
      case nue:
        //nueFluxHisto->Fill(enu, weight);
            
        //    Colton's Efficiency
        //    0    0.25    0.00380228
        //    0.25    0.5    0.0666667
        //    0.5    0.75    0.104741
        //    0.75    1    0.121728
        //    1    1.5    0.114269
        //    1.5    4    0.0788979
            
        //    now just multiplying the flux histograms by the efficiencies in the ranges given
            nueFluxHisto_orig->Fill(enu,weight);
            // with tgptype, pL_tgp and pT_tgp, since thin target. /////////////////////////////////////////////
            
            
            if ( 0<enu && enu<0.25){
                if (fluxNtuple -> tgptype==11 || fluxNtuple -> tgptype==12 || fluxNtuple ->tgptype==8 || fluxNtuple ->tgptype==9){ //Kplus or Kminus or pionplus or pionminus
                    if (pL_tgp < 20 || (pL_tgp < 24 && pT_tgp > 1) || (pL_tgp<31 && pT_tgp>1.2) || (pL_tgp<42 && pT_tgp>1.55) || pL_tgp>90 || pT_tgp>2 ){
                        nueFluxHisto_up_tgp -> Fill(enu, weight*(1.25)*0.00380228);
                        nueFluxHisto_down_tgp -> Fill(enu, weight*(0.75)*0.00380228);
                       // std::cout<<"up: "<<weight*(1.25)*0.00380228<<" down: "<<weight*(0.75)*0.00380228<<" normal: "<<weight*0.00380228<<std::endl;
                    }
                
                else{
                    nueFluxHisto_up_tgp -> Fill(enu, weight*0.00380228);
                    nueFluxHisto_down_tgp -> Fill(enu, weight*0.00380228);
                }
                }
                else{
                    nueFluxHisto_up_tgp -> Fill(enu, weight*(1.25)*0.00380228);
                    nueFluxHisto_down_tgp -> Fill(enu, weight*(0.75)*0.00380228);
                }
                
            
            nueFluxHisto_tgp->Fill(enu, weight*0.00380228);
            }
            
            //////
            
            if ( 0.25<enu && enu<0.5){
                if (fluxNtuple -> tgptype==11 || fluxNtuple -> tgptype==12 || fluxNtuple ->tgptype==8 || fluxNtuple ->tgptype==9){ //Kplus or Kminus or pionplus or pionminus
                    if (pL_tgp < 20 || (pL_tgp < 24 && pT_tgp > 1) || (pL_tgp<31 && pT_tgp>1.2) || (pL_tgp<42 && pT_tgp>1.55) || pL_tgp>90 || pT_tgp>2 ){
                        nueFluxHisto_up_tgp -> Fill(enu, weight*(1.25)*0.0666667);
                        nueFluxHisto_down_tgp -> Fill(enu, weight*0.75*0.0666667);
                    }
                
                else{
                    nueFluxHisto_up_tgp -> Fill(enu, weight*0.0666667);
                    nueFluxHisto_down_tgp -> Fill(enu, weight*0.0666667);
                }
                }
                else{
                    nueFluxHisto_up_tgp -> Fill(enu, weight*(1.25)*0.0666667);
                    nueFluxHisto_down_tgp -> Fill(enu, weight*(0.75)*0.0666667);
                }
            
            nueFluxHisto_tgp->Fill(enu, weight*0.0666667);

            }
            
            //////
            
            if ( 0.5<enu && enu<0.75){
                if (fluxNtuple -> tgptype==11 || fluxNtuple -> tgptype==12 || fluxNtuple ->tgptype==8 || fluxNtuple ->tgptype==9){ //Kplus or Kminus or pionplus or pionminus
                    if (pL_tgp < 20 || (pL_tgp < 24 && pT_tgp > 1) || (pL_tgp<31 && pT_tgp>1.2) || (pL_tgp<42 && pT_tgp>1.55) || pL_tgp>90 || pT_tgp>2 ){
                        nueFluxHisto_up_tgp -> Fill(enu, weight*(1.25)*0.104741);
                        nueFluxHisto_down_tgp -> Fill(enu, weight*(0.75)*0.104741);
                    }
                
                else{
                    nueFluxHisto_up_tgp -> Fill(enu, weight*0.104741);
                    nueFluxHisto_down_tgp -> Fill(enu, weight*0.104741);
                }
                }
                else{
                    nueFluxHisto_up_tgp -> Fill(enu, weight*(1.25)*0.104741);
                    nueFluxHisto_down_tgp -> Fill(enu, weight*(0.75)*0.104741);
                }
            
            nueFluxHisto_tgp->Fill(enu, weight*0.104741);
            }
            
            //////
            
            if ( 0.75<enu && enu<1){
                if (fluxNtuple -> tgptype==11 || fluxNtuple -> tgptype==12 || fluxNtuple ->tgptype==8 || fluxNtuple ->tgptype==9){ //Kplus or Kminus or pionplus or pionminus
                    if (pL_tgp < 20 || (pL_tgp < 24 && pT_tgp > 1) || (pL_tgp<31 && pT_tgp>1.2) || (pL_tgp<42 && pT_tgp>1.55) || pL_tgp>90 || pT_tgp>2 ){
                        nueFluxHisto_up_tgp -> Fill(enu, weight*(1.25)*0.121728);
                        nueFluxHisto_down_tgp -> Fill(enu, weight*(0.75)*0.121728);
                    }
                
                else{
                    nueFluxHisto_up_tgp -> Fill(enu, weight*0.121728);
                    nueFluxHisto_down_tgp -> Fill(enu, weight*0.121728);
                }
                }
                else{
                    nueFluxHisto_up_tgp -> Fill(enu, weight*(1.25)*0.121728);
                    nueFluxHisto_down_tgp -> Fill(enu, weight*(0.75)*0.121728);
                }
                
            
            nueFluxHisto_tgp->Fill(enu, weight*0.121728);
            }
            
            ////
            
            if ( 1<enu && enu<1.5){
                
              if (fluxNtuple -> tgptype==11 || fluxNtuple -> tgptype==12 || fluxNtuple ->tgptype==8 || fluxNtuple ->tgptype==9){ //Kplus or Kminus or pionplus or pionminus
                    if (pL_tgp < 20 || (pL_tgp < 24 && pT_tgp > 1) || (pL_tgp<31 && pT_tgp>1.2) || (pL_tgp<42 && pT_tgp>1.55) || pL_tgp>90 || pT_tgp>2 ){
                        nueFluxHisto_up_tgp -> Fill(enu, weight*(1.25)*0.114269);
                        nueFluxHisto_down_tgp -> Fill(enu, weight*(0.75)*0.114269);
                    }
                
                else{
                    nueFluxHisto_up_tgp -> Fill(enu, weight*0.114269);
                    nueFluxHisto_down_tgp -> Fill(enu, weight*0.114269);
                }
              }
              else{
                  nueFluxHisto_up_tgp -> Fill(enu, weight*(1.25)*0.114269);
                  nueFluxHisto_down_tgp -> Fill(enu, weight*(0.75)*0.114269);
              }
                
            
            nueFluxHisto_tgp->Fill(enu, weight*0.114269);
            }
            
            /////
            
            if ( 1.5<enu && enu<4){
                if (fluxNtuple -> tgptype==11 || fluxNtuple -> tgptype==12 || fluxNtuple ->tgptype==8 || fluxNtuple ->tgptype==9){ //Kplus or Kminus or pionplus or pionminus
                    if (pL_tgp < 20 || (pL_tgp < 24 && pT_tgp > 1) || (pL_tgp<31 && pT_tgp>1.2) || (pL_tgp<42 && pT_tgp>1.55) || pL_tgp>90 || pT_tgp>2 ){
                        nueFluxHisto_up_tgp -> Fill(enu, weight*(1.25)*0.0788979);
                        nueFluxHisto_down_tgp -> Fill(enu, weight*(0.75)*0.0788979);
                    }
                else{
                    nueFluxHisto_up_tgp -> Fill(enu, weight*0.0788979);
                    nueFluxHisto_down_tgp -> Fill(enu, weight*0.0788979);
                }
                }
                else{
                    nueFluxHisto_up_tgp -> Fill(enu, weight*(1.25)*0.0788979);
                    nueFluxHisto_down_tgp -> Fill(enu, weight*(0.75)*0.0788979);
                }
                nueFluxHisto_tgp->Fill(enu, weight*0.0788979);
            
            }
            
            
            
            if ( enu>4){   //no efficiency
                if (fluxNtuple -> tgptype==11 || fluxNtuple -> tgptype==12 || fluxNtuple ->tgptype==8 || fluxNtuple ->tgptype==9){ //Kplus or Kminus or pionplus or pionminus
                    if (pL_tgp < 20 || (pL_tgp < 24 && pT_tgp > 1) || (pL_tgp<31 && pT_tgp>1.2) || (pL_tgp<42 && pT_tgp>1.55) || pL_tgp>90 || pT_tgp>2 ){
                        nueFluxHisto_up_tgp -> Fill(enu, weight*(1.25));
                        nueFluxHisto_down_tgp -> Fill(enu, weight*(0.75));
                    }
                    else{
                        nueFluxHisto_up_tgp -> Fill(enu, weight);
                        nueFluxHisto_down_tgp -> Fill(enu, weight);
                    }
                }
                    else{
                        nueFluxHisto_up_tgp -> Fill(enu, (1.25)*weight);
                        nueFluxHisto_down_tgp -> Fill(enu, (0.75)*weight);
                        
                    }
                nueFluxHisto_tgp->Fill(enu, weight);
            }
            
            
            ////// filling in MIPP constrained region ////////
                if (fluxNtuple -> tgptype==11 || fluxNtuple -> tgptype==12 || fluxNtuple ->tgptype==8 || fluxNtuple ->tgptype==9){ //Kplus or Kminus or pionplus or pionminus
                    if ( (pL_tgp > 20 && pL_tgp<24 && pT_tgp<1)||
                        (pL_tgp > 24 && pL_tgp<31 && pT_tgp<1.2) ||
                        (pL_tgp > 31 && pL_tgp<42 && pT_tgp<1.55) ||
                        (pL_tgp > 42 && pL_tgp<90 && pT_tgp<2)
                    ){
                        nueFluxHisto_orig_MIPP_con->Fill(enu,weight);
                    }
                }
        
        ////// filling in NA49 constrained region for pionplus ////////
        if (fluxNtuple ->tgptype==8 || fluxNtuple ->tgptype==9){ //pionplus or pionminus
            if (  Xf>0.55 ||
                ( Xf>0.45 && pT_tgp>1.5) ||
                ( Xf>0.25 && pT_tgp>1.7) ||
                ( pT_tgp>1.9) ||
                ( Xf<0.15 && Xf>0.125 && pT_tgp>1.3) ||
                ( Xf< -0.125) ||
                ( Xf< -0.125 && pT_tgp<1.2) ||
                ( Xf<0.225 && pT_tgp>1.05 && pT_tgp<1.1) ||
                ( Xf<0.125 && pT_tgp>0.1125 && pT_tgp<1.3 && pT_tgp>0.65) ||
                ( Xf<-0.0875 && pT_tgp<0.35 ) ||
                ( Xf<-0.0625 && pT_tgp<0.25 ) ||
                ( Xf<-0.055 && Xf>-0.0625 && pT_tgp<0.75 ) ||
                ( Xf<0.225 && Xf>-0.055 && pT_tgp>0.325 && pT_tgp<0.35) ||
                ( Xf<0.175 && Xf>0.17 && pT_tgp<0.65 && pT_tgp>0.35 ) ||
                ( Xf<-0.045 && pT_tgp<0.175) ||
                ( Xf<-0.025 && pT_tgp<0.075) ||
                ( pT_tgp<0.025) ||
                ( Xf<0.15 && Xf>0.1625 && pT_tgp<0.325) ||
                ( Xf>0.225 && pT_tgp< 0.045) ||
                ( pT_tgp>0.045 && Xf>0.325 && Xf<0.35 && pT_tgp<1.3) ||
                ( Xf<0.325 && Xf>0.225 && pT_tgp>0.65 && pT_tgp<0.7 ) ||
                ( Xf >0.35 && Xf<0.55 && pT_tgp>0.65 && pT_tgp<0.7)
            ){
            }
            
            else{nueFluxHisto_orig_NA49_con->Fill(enu,weight);
                nueFluxHisto_orig_NA49_pion_con->Fill(enu,weight);
            }
            
            }
    
        ////// filling in NA49 constrained region for kaonplus ////////
        if (fluxNtuple ->tgptype==11 || fluxNtuple ->tgptype==12){ //kaonplus or kaonminus
            if (  (pT_tgp<1 && pT_tgp>0.8 && Xf<0.225 && Xf>-0.0125)||
                ( pT_tgp<0.75 && pT_tgp>0.05 && Xf<0.225 && Xf>-0.0125)
            ){
                nueFluxHisto_orig_NA49_con->Fill(enu,weight);
                nueFluxHisto_orig_NA49_kaon_con->Fill(enu,weight);
            }
            
        }
      
      ////// filling in NA49 constrained region for proton ////////
      if (fluxNtuple ->tgptype==14){ //proton
          if ( Xf>-0.0125 && Xf<0.975 && pT_tgp<1.925 && pT_tgp>0){
              nueFluxHisto_orig_NA49_con->Fill(enu,weight);
              nueFluxHisto_orig_NA49_proton_con->Fill(enu,weight);
          }
      }
        
        ////// filling in NA49 constrained region for neutron ////////
        if (fluxNtuple ->tgptype==13){ //neutron
            if ( (Xf>0.475 && Xf<1 && pT_tgp<2.4 && pT_tgp>0) ||
                ( Xf>0.05 && Xf<0.475 && pT_tgp<0.6 && pT_tgp>0)){
                nueFluxHisto_orig_NA49_con->Fill(enu,weight);
                nueFluxHisto_orig_NA49_neutron_con->Fill(enu,weight);
            }
            if ( Xf>0.05 && Xf<0.475 && pT_tgp<2.4 && pT_tgp>0.6){
                double yvalue = (2.4-0.6)/(0.475-0.05) * Xf + 0.38823529;
                if ( pT_tgp < yvalue){
                    nueFluxHisto_orig_NA49_con->Fill(enu,weight);
                    nueFluxHisto_orig_NA49_neutron_con->Fill(enu,weight);
                }
            }
        }
        
    /////////////////
            
            
            
            
            if(parent_mom2>0){                          // new noDAR
                nueFluxHisto_noDAR->Fill(enu, weight);
            }
        
            // plot Ndecays 1, 5, 6, 11
            
            if (fluxNtuple->Ndecay==1)
            {
                nue_Ndecay_one->Fill(enu,weight);
            }
            
            if (fluxNtuple->Ndecay==5)
            {
                nue_Ndecay_five->Fill(enu,weight);
            }
            
            if (fluxNtuple->Ndecay==6)
            {
                nue_Ndecay_six->Fill(enu,weight);
            }
        
            if (fluxNtuple->Ndecay==11)
            {
                nue_Ndecay_eleven->Fill(enu,weight);
            }
            
        if (fluxNtuple->ptype == 5)                 // nue_parent
        {
            nue_parent_muonplus->Fill(enu,weight);
            
            if (fluxNtuple -> tgptype == 11 )// kaonplus tgptype
            {   nue_tgptype_muonplus_kaonplus->Fill(enu, weight);
            }
            
            if (fluxNtuple -> tgptype == 10 )// kaonplus tgptype
            {   nue_tgptype_muonplus_kaonlong->Fill(enu, weight);
            }
        }
        if (fluxNtuple->ptype == 11)
        {nue_parent_kaonplus->Fill(enu,weight);
        }
        if (fluxNtuple->ptype == 10)
        {nue_parent_kaonlong->Fill(enu,weight);
        }
            
        if (fluxNtuple -> tgptype == 8)
        {nue_tgptype_pionplus -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 9)
        {nue_tgptype_pionminus -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 10)
        {nue_tgptype_kaonlong -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 11)
        {nue_tgptype_kaonplus -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 12)
        {nue_tgptype_kaonminus -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 13)
        {nue_tgptype_neutron -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 14)
        {nue_tgptype_proton -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 18)
        {nue_tgptype_lambda -> Fill(enu,weight);
        }
        
        if (fluxNtuple -> tgptype == 19)
        {nue_tgptype_sigmaplus -> Fill(enu,weight);
        }
            
        break;
            
      case anue:
        anueFluxHisto->Fill(enu, weight);
            //std::cout<< "Normal: "<<"enu: "<<enu<<" weight: "<<weight<<std::endl;
            if (fluxNtuple -> tptype==11 || fluxNtuple -> tptype==12 ){ //Kplus or Kminus
                if (pZ < 20 || (pZ < 24 && pT > 1) || (pZ<31 && pT>1.2) || (pZ<42 && pT>1.55) || pZ>90 || pT>2 ){

                    anueFluxHisto_up -> Fill(enu, weight*1.75);
                    anueFluxHisto_down -> Fill(enu, weight*0.75);
                }
            }
            else {
                anueFluxHisto_up -> Fill(enu, weight);
                anueFluxHisto_down -> Fill(enu, weight);
            }
            
            
        if(parent_mom2>0){                          // new noDAR
        anueFluxHisto_noDAR->Fill(enu, weight);  
        }
            if (fluxNtuple->ptype == 11)
            {anue_parent_kaonplus->Fill(enu,weight);
            }
        break; 
    }

    // POT stuff
    if ( fluxNtuple->evtno > highest_evtno ) 
      highest_evtno = fluxNtuple->evtno;

  } // end of loop over the entries


  //***************************************
  //
  // POT scaling
  //
  //***************************************

  AccumulatedPOT += estimate_pots(highest_evtno); // To account for last tree
  double scale = NominalPOT/AccumulatedPOT;
  numuFluxHisto  -> Scale(scale);
  anumuFluxHisto -> Scale(scale);
  nueFluxHisto   -> Scale(scale);
  anueFluxHisto  -> Scale(scale);
  anueFluxHisto_noDAR -> Scale(scale);
  nueFluxHisto_noDAR -> Scale(scale);
    
    nueFluxHisto_orig->Scale(scale);
    nueFluxHisto_orig_MIPP_con->Scale(scale);
    nueFluxHisto_orig_NA49_con->Scale(scale);
    nueFluxHisto_orig_NA49_proton_con->Scale(scale);
    nueFluxHisto_orig_NA49_kaon_con->Scale(scale);
    nueFluxHisto_orig_NA49_pion_con->Scale(scale);
    nueFluxHisto_orig_NA49_neutron_con->Scale(scale);
    nueFluxHisto_tgp->Scale(scale);
    
    numuFluxHisto_up  -> Scale(scale);
    anumuFluxHisto_up -> Scale(scale);
    nueFluxHisto_up   -> Scale(scale);
    anueFluxHisto_up  -> Scale(scale);
    nueFluxHisto_up_tgp   -> Scale(scale);
    nueFluxHisto_down_tgp  -> Scale(scale);
    
    numuFluxHisto_down  -> Scale(scale);
    anumuFluxHisto_down -> Scale(scale);
    nueFluxHisto_down   -> Scale(scale);
    anueFluxHisto_down  -> Scale(scale);
    
    pionplus_MIPP-> Scale(scale);
    pionplus_MIPP_nue-> Scale(scale);
    Kplus_MIPP-> Scale(scale);
    Kplus_MIPP_nue-> Scale(scale);
    
    nue_parent_kaonlong -> Scale(scale);
    nue_parent_kaonplus -> Scale(scale);
    nue_parent_muonplus -> Scale(scale);
    
    numu_parent_kaonplus -> Scale(scale);
    numu_parent_pionplus -> Scale(scale);
    numu_parent_muonminus -> Scale(scale);
    numu_parent_kaonlong -> Scale(scale);
    
    nue_Ndecay_one -> Scale(scale);
    nue_Ndecay_six -> Scale(scale);
    nue_Ndecay_five -> Scale(scale);
    nue_Ndecay_eleven -> Scale(scale);
    
    numu_Ndecay_five -> Scale(scale);
    numu_Ndecay_seven -> Scale(scale);
    numu_Ndecay_twelve -> Scale(scale);
    numu_Ndecay_thirteen -> Scale(scale);
    
  cout << endl << ">>> TOTAL POT: " << AccumulatedPOT << endl << endl;


  //***************************************
  //
  // Apply now GENIE xsec
  //
  // source /nusoft/app/externals/setup
  // setup genie_xsec R-2_8_0   -q default
  // root -l  $GENIEXSECPATH/xsec_graphs.root
  // >  _file0->cd("nu_mu_Ar40")
  // >  tot_cc->Draw()
  //
  //***************************************
    
  const char* genieXsecPath = gSystem->ExpandPathName("$(GENIEXSECPATH)");
  if ( !genieXsecPath ) {
    std::cout << "$(GENIEXSECPATH) not defined." << std::endl;
    std::cout << "Please setup *genie_xsec*. (setup genie_xsec R-2_8_0   -q default)." << std::endl; 
  }

  if ( genieXsecPath ) {
    TString genieXsecFileName = genieXsecPath;
    genieXsecFileName += "/xsec_graphs.root";
    TFile *genieXsecFile = new TFile(genieXsecFileName,"READ");
    genieXsecFile->cd("nu_e_Ar40");
    genieXsecNueCC = (TGraph *) gDirectory->Get("tot_cc");
    genieXsecFile->Close();

    // TSpline3* genieXsecSplineNumuCC = new TSpline3("genieXsecSplineNumuCC", genieXsecNumuCC, "", 0,6);

    double value_up;
    double value_down;
    double value_up_tgp;
    double value_down_tgp;
    double value;
    double value_tgp;
      
    for(int i=1; i<histNbins+1; i++) {
        
      value = nueFluxHisto->GetBinContent(i);  // get flux value
      value_tgp= nueFluxHisto_tgp->GetBinContent(i);
        
      value_up = nueFluxHisto_up->GetBinContent(i);
      value_down = nueFluxHisto_down->GetBinContent(i);
        
      value_up_tgp = nueFluxHisto_up_tgp->GetBinContent(i);
      value_down_tgp = nueFluxHisto_down_tgp->GetBinContent(i);
        
       
      value *= genieXsecNueCC->Eval(nueFluxHisto->GetBinCenter(i)); // times flux value by xsec
      value_tgp *= genieXsecNueCC->Eval(nueFluxHisto_tgp->GetBinCenter(i));
        
      value_up *= genieXsecNueCC->Eval(nueFluxHisto_up->GetBinCenter(i));
      value_down *= genieXsecNueCC->Eval(nueFluxHisto_down->GetBinCenter(i));
        
      value_up_tgp *= genieXsecNueCC->Eval(nueFluxHisto_up_tgp->GetBinCenter(i));
      value_down_tgp *= genieXsecNueCC->Eval(nueFluxHisto_down_tgp->GetBinCenter(i));
        
       
      value *= (1e-38 * Ntarget/40.);
      value_tgp *= (1e-38 * Ntarget/40.);
        
      value_up *= (1e-38 * Ntarget/40.);
      value_down *= (1e-38 * Ntarget/40.); // 1/40 is due to I'm considering nu_mu_Ar40.
        
      value_up_tgp *= (1e-38 * Ntarget/40.);
      value_down_tgp *= (1e-38 * Ntarget/40.); // 1/40 is due to I'm considering nu_mu_Ar40.
        
      nueCCHisto->SetBinContent(i,value);
      nueCCHisto_tgp->SetBinContent(i,value_tgp);
      nueCCHisto_up->SetBinContent(i,value_up);
      nueCCHisto_down->SetBinContent(i,value_down);
      nueCCHisto_up_tgp->SetBinContent(i,value_up_tgp);
      nueCCHisto_down_tgp->SetBinContent(i,value_down_tgp);
        
    }
  } // end if ( genieXsecPath )




  //***************************************
  //
  // Writing on file
  //
  //***************************************

  f->cd();
  numuFluxHisto  -> Write();
  anumuFluxHisto -> Write();
  nueFluxHisto   -> Write();
  anueFluxHisto  -> Write();
  anueFluxHisto_noDAR -> Write();
  nueFluxHisto_noDAR -> Write();
    
    nueFluxHisto_orig_MIPP_con->Write();
    nueFluxHisto_orig_NA49_con->Write();
    nueFluxHisto_orig_NA49_proton_con->Write();
    nueFluxHisto_orig_NA49_pion_con->Write();
    nueFluxHisto_orig_NA49_kaon_con->Write();
    nueFluxHisto_orig_NA49_neutron_con->Write();
    nueFluxHisto_orig->Write();
    nueFluxHisto_tgp->Write();
    
    nue_parent_kaonplus -> Write();
    nue_parent_kaonlong -> Write();
    nue_parent_muonplus -> Write();
    
    anue_parent_kaonplus -> Write();
    anumu_parent_kaonplus -> Write();
    
    numu_parent_kaonlong ->Write();
    numu_parent_kaonplus -> Write();
    numu_parent_pionplus -> Write();
    numu_parent_muonplus -> Write();
    numu_parent_muonminus -> Write();
    
    
    nue_tgptype_muonplus_kaonlong -> Write();
    nue_tgptype_muonplus_kaonplus -> Write();
    
    nue_Ndecay_one -> Write();
    nue_Ndecay_six -> Write();
    nue_Ndecay_eleven -> Write();
    nue_Ndecay_five -> Write();
    
    numu_Ndecay_thirteen -> Write();
    numu_Ndecay_seven -> Write();
    numu_Ndecay_three -> Write();
    numu_Ndecay_twelve -> Write();
    numu_Ndecay_five -> Write();
    
    nue_tgptype_pionplus -> Write();
    nue_tgptype_pionminus -> Write();
    nue_tgptype_kaonlong -> Write();
    nue_tgptype_kaonplus -> Write();
    nue_tgptype_kaonminus -> Write();
    nue_tgptype_neutron -> Write();
    nue_tgptype_proton -> Write();
    nue_tgptype_lambda -> Write();
    nue_tgptype_sigmaplus -> Write();
    
    numu_tgptype_pionplus -> Write();
    numu_tgptype_pionminus -> Write();
    numu_tgptype_kaonlong -> Write();
    numu_tgptype_kaonplus -> Write();
    numu_tgptype_kaonminus -> Write();
    numu_tgptype_neutron -> Write();
    numu_tgptype_proton -> Write();
    numu_tgptype_lambda -> Write();
    numu_tgptype_sigmaplus -> Write();
    
    numu_tgen_two -> Write();
    numu_tgen_three -> Write();
    numu_tgen_four -> Write();
    numu_tgen_five -> Write();
    numu_tgen_six -> Write();
    numu_tgen_seven -> Write();
    numu_tgen_eight -> Write();
    numu_tgen_nine -> Write();
    numu_tgen_ten -> Write();
    numu_tgen_eleven -> Write();
    numu_tgen_twelve -> Write();
    

    pionplusMomentumHisto -> Write();
    pionplusAngleHisto -> Write();
    pionplus2DHisto -> Write();
    pionplus_MIPP -> Write();
    pionplus_MIPP_nue->Write();
    pionplus_MIPP_nue_dm1->Write();
    pionplus_MIPP_nue_dm6->Write();
    pionplus_MIPP_nue_dm5->Write();
    pionplus_MIPP_nue_dm11->Write();
    pionplus_MIPP_numu->Write();
    pionplus_MIPP_anue->Write();
    pionplus_MIPP_anumu->Write();
    pionplus_nue->Write();
    pionplus_numu->Write();
    pionplus_anue->Write();
    pionplus_anumu->Write();
    
    pionminusMomentumHisto -> Write();
    pionminusAngleHisto -> Write();
    pionminus2DHisto -> Write();
    pionminus_MIPP -> Write();
    pionminus_MIPP_nue->Write();
    pionminus_MIPP_numu->Write();
    pionminus_MIPP_anue->Write();
    pionminus_MIPP_anumu->Write();
    pionminus_nue->Write();
    pionminus_numu->Write();
    pionminus_anue->Write();
    pionminus_anumu->Write();
    
    KlongMomentumHisto -> Write();
    KlongAngleHisto -> Write();
    Klong2DHisto -> Write();
    Klong_nue->Write();
    Klong_numu->Write();
    Klong_anue->Write();
    Klong_anumu->Write();
    
    KplusMomentumHisto -> Write();
    KplusAngleHisto -> Write();
    Kplus2DHisto -> Write();
    Kplus_MIPP -> Write();
    Kplus_MIPP_nue->Write();
    Kplus_MIPP_numu->Write();
    Kplus_MIPP_anue->Write();
    Kplus_MIPP_anumu->Write();
    Kplus_nue->Write();
    Kplus_numu->Write();
    Kplus_anue->Write();
    Kplus_anumu->Write();
    
    KminusMomentumHisto -> Write();
    KminusAngleHisto -> Write();
    Kminus2DHisto -> Write();
    Kminus_MIPP -> Write();
    Kminus_MIPP_nue->Write();
    Kminus_MIPP_numu->Write();
    Kminus_MIPP_anue->Write();
    Kminus_MIPP_anumu->Write();
    Kminus_MIPP_nue_dm_eleven->Write();
    Kminus_MIPP_nue_dm_six->Write();
    Kminus_MIPP_nue_dm_one->Write();
    Kminus_nue->Write();
    Kminus_numu->Write();
    Kminus_anue->Write();
    Kminus_anumu->Write();
    
    muonminus2DHisto->Write();
    muonplus2DHisto->Write();
    
    protonMomentumHisto -> Write();
    protonAngleHisto -> Write();
    proton2DHisto->Write();
    proton_nue->Write();
    proton_numu->Write();
    proton_anue->Write();
    proton_anumu->Write();
    
    neutronMomentumHisto -> Write();
    neutronAngleHisto -> Write();
    neutron2DHisto->Write();
    neutron_nue->Write();
    neutron_numu->Write();
    neutron_anue->Write();
    neutron_anumu->Write();
    
    lambdaMomentumHisto -> Write();
    lambdaAngleHisto -> Write();
    lambda2DHisto -> Write();
    //lambda_nue->Write();
   // lambda_numu->Write();
    
    sigmaplusMomentumHisto -> Write();
    sigmaplusAngleHisto -> Write();
    sigmaplus2DHisto -> Write();
   // sigmaplus_nue->Write();
   // sigmaplus_numu->Write();
    
    ppnue_dm_one->Write();
    ppnue_dm_two->Write();
    ppnue_dm_three->Write();
    ppnue_dm_four->Write();
    ppnue_dm_five->Write();
    ppnue_dm_six->Write();
    ppnue_dm_seven->Write();
    ppnue_dm_eight -> Write();
    ppnue_dm_nine -> Write();
    ppnue_dm_ten -> Write();
    ppnue_dm_eleven -> Write();
    ppnue_dm_twelve -> Write();
    ppnue_dm_thirteen -> Write();
    ppnue_dm_fourteen -> Write();
    
    pmnue_dm_one->Write();
    pmnue_dm_two->Write();
    pmnue_dm_three->Write();
    pmnue_dm_four->Write();
    pmnue_dm_five->Write();
    pmnue_dm_six->Write();
    pmnue_dm_seven->Write();
    pmnue_dm_eight -> Write();
    pmnue_dm_nine -> Write();
    pmnue_dm_ten -> Write();
    pmnue_dm_eleven -> Write();
    pmnue_dm_twelve -> Write();
    pmnue_dm_thirteen -> Write();
    pmnue_dm_fourteen -> Write();
    
    pnue_dm_one->Write();
    pnue_dm_two->Write();
    pnue_dm_three->Write();
    pnue_dm_four->Write();
    pnue_dm_five->Write();
    pnue_dm_six->Write();
    pnue_dm_seven->Write();
    pnue_dm_eight -> Write();
    pnue_dm_nine -> Write();
    pnue_dm_ten -> Write();
    pnue_dm_eleven -> Write();
    pnue_dm_twelve -> Write();
    pnue_dm_thirteen -> Write();
    pnue_dm_fourteen -> Write();
    
    nnue_dm_one->Write();
    nnue_dm_two->Write();
    nnue_dm_three->Write();
    nnue_dm_four->Write();
    nnue_dm_five->Write();
    nnue_dm_six->Write();
    nnue_dm_seven->Write();
    nnue_dm_eight -> Write();
    nnue_dm_nine -> Write();
    nnue_dm_ten -> Write();
    nnue_dm_eleven -> Write();
    nnue_dm_twelve -> Write();
    nnue_dm_thirteen -> Write();
    nnue_dm_fourteen -> Write();
    
    klnue_dm_one->Write();
    klnue_dm_two->Write();
    klnue_dm_three->Write();
    klnue_dm_four->Write();
    klnue_dm_five->Write();
    klnue_dm_six->Write();
    klnue_dm_seven->Write();
    klnue_dm_eight -> Write();
    klnue_dm_nine -> Write();
    klnue_dm_ten -> Write();
    klnue_dm_eleven -> Write();
    klnue_dm_twelve -> Write();
    klnue_dm_thirteen -> Write();
    klnue_dm_fourteen -> Write();
    
    kpnue_dm_one->Write();
    kpnue_dm_two->Write();
    kpnue_dm_three->Write();
    kpnue_dm_four->Write();
    kpnue_dm_five->Write();
    kpnue_dm_six->Write();
    kpnue_dm_seven->Write();
    kpnue_dm_eight -> Write();
    kpnue_dm_nine -> Write();
    kpnue_dm_ten -> Write();
    kpnue_dm_eleven -> Write();
    kpnue_dm_twelve -> Write();
    kpnue_dm_thirteen -> Write();
    kpnue_dm_fourteen -> Write();
    
    kmnue_dm_one->Write();
    kmnue_dm_two->Write();
    kmnue_dm_three->Write();
    kmnue_dm_four->Write();
    kmnue_dm_five->Write();
    kmnue_dm_six->Write();
    kmnue_dm_seven->Write();
    kmnue_dm_eight -> Write();
    kmnue_dm_nine -> Write();
    kmnue_dm_ten -> Write();
    kmnue_dm_eleven -> Write();
    kmnue_dm_twelve -> Write();
    kmnue_dm_thirteen -> Write();
    kmnue_dm_fourteen -> Write();
    
    
    ppnumu_dm_one->Write();
    ppnumu_dm_two->Write();
    ppnumu_dm_three->Write();
    ppnumu_dm_four->Write();
    ppnumu_dm_five->Write();
    ppnumu_dm_six->Write();
    ppnumu_dm_seven->Write();
    ppnumu_dm_eight -> Write();
    ppnumu_dm_nine -> Write();
    ppnumu_dm_ten -> Write();
    ppnumu_dm_eleven -> Write();
    ppnumu_dm_twelve -> Write();
    ppnumu_dm_thirteen -> Write();
    ppnumu_dm_fourteen -> Write();
    
    pmnumu_dm_one->Write();
    pmnumu_dm_two->Write();
    pmnumu_dm_three->Write();
    pmnumu_dm_four->Write();
    pmnumu_dm_five->Write();
    pmnumu_dm_six->Write();
    pmnumu_dm_seven->Write();
    pmnumu_dm_eight -> Write();
    pmnumu_dm_nine -> Write();
    pmnumu_dm_ten -> Write();
    pmnumu_dm_eleven -> Write();
    pmnumu_dm_twelve -> Write();
    pmnumu_dm_thirteen -> Write();
    pmnumu_dm_fourteen -> Write();
    
    pnumu_dm_one->Write();
    pnumu_dm_two->Write();
    pnumu_dm_three->Write();
    pnumu_dm_four->Write();
    pnumu_dm_five->Write();
    pnumu_dm_six->Write();
    pnumu_dm_seven->Write();
    pnumu_dm_eight -> Write();
    pnumu_dm_nine -> Write();
    pnumu_dm_ten -> Write();
    pnumu_dm_eleven -> Write();
    pnumu_dm_twelve -> Write();
    pnumu_dm_thirteen -> Write();
    pnumu_dm_fourteen -> Write();
    
    nnumu_dm_one->Write();
    nnumu_dm_two->Write();
    nnumu_dm_three->Write();
    nnumu_dm_four->Write();
    nnumu_dm_five->Write();
    nnumu_dm_six->Write();
    nnumu_dm_seven->Write();
    nnumu_dm_eight -> Write();
    nnumu_dm_nine -> Write();
    nnumu_dm_ten -> Write();
    nnumu_dm_eleven -> Write();
    nnumu_dm_twelve -> Write();
    nnumu_dm_thirteen -> Write();
    nnumu_dm_fourteen -> Write();
    
    klnumu_dm_one->Write();
    klnumu_dm_two->Write();
    klnumu_dm_three->Write();
    klnumu_dm_four->Write();
    klnumu_dm_five->Write();
    klnumu_dm_six->Write();
    klnumu_dm_seven->Write();
    klnumu_dm_eight -> Write();
    klnumu_dm_nine -> Write();
    klnumu_dm_ten -> Write();
    klnumu_dm_eleven -> Write();
    klnumu_dm_twelve -> Write();
    klnumu_dm_thirteen -> Write();
    klnumu_dm_fourteen -> Write();
    
    kpnumu_dm_one->Write();
    kpnumu_dm_two->Write();
    kpnumu_dm_three->Write();
    kpnumu_dm_four->Write();
    kpnumu_dm_five->Write();
    kpnumu_dm_six->Write();
    kpnumu_dm_seven->Write();
    kpnumu_dm_eight -> Write();
    kpnumu_dm_nine -> Write();
    kpnumu_dm_ten -> Write();
    kpnumu_dm_eleven -> Write();
    kpnumu_dm_twelve -> Write();
    kpnumu_dm_thirteen -> Write();
    kpnumu_dm_fourteen -> Write();
    
    kmnumu_dm_one->Write();
    kmnumu_dm_two->Write();
    kmnumu_dm_three->Write();
    kmnumu_dm_four->Write();
    kmnumu_dm_five->Write();
    kmnumu_dm_six->Write();
    kmnumu_dm_seven->Write();
    kmnumu_dm_eight -> Write();
    kmnumu_dm_nine -> Write();
    kmnumu_dm_ten -> Write();
    kmnumu_dm_eleven -> Write();
    kmnumu_dm_twelve -> Write();
    kmnumu_dm_thirteen -> Write();
    kmnumu_dm_fourteen -> Write();
    
    pionplus_NA49 -> Write();
    pionplus_NA49_nue -> Write();
    
    pionminus_NA49 -> Write();
    pionminus_NA49_nue -> Write();
    
    Kplus_NA49 -> Write();
    Kplus_NA49_second -> Write();
    Kplus_NA49_nue -> Write();
    
    Kminus_NA49 -> Write();
    Kminus_NA49_nue -> Write();
    
    proton_NA49 -> Write();
    proton_NA49_nue -> Write();
    
    neutron_NA49 -> Write();
    neutron_NA49_nue -> Write();
    
    nueFluxHisto_up->Write();
    numuFluxHisto_up->Write();
    anueFluxHisto_up->Write();
    anumuFluxHisto_up->Write();
    nueFluxHisto_up_tgp->Write();
    
    nueFluxHisto_down->Write();
    nueFluxHisto_down_tgp->Write();
    numuFluxHisto_down->Write();
    anueFluxHisto_down->Write();
    anumuFluxHisto_down->Write();
    
    
  if ( genieXsecPath ) {
    nueCCHisto     -> Write();
      nueCCHisto_tgp ->Write();
    genieXsecNueCC -> Write();
    nueCCHisto_up->Write();
    nueCCHisto_down->Write();
    nueCCHisto_up_tgp->Write();
    nueCCHisto_down_tgp->Write();
  }
    
  f->Close();

}


//___________________________________________________________________________
TVector3 NuMIFlux::RandomInTPC() {
    
    TDatime *d = new TDatime;
    TRandom *r = new TRandom(d->GetTime());
    
    double xTPC = 256.35;  // cm
    double yTPC = 233.;  // cm
    double zTPC = 1036.8; // cm
    
    double x = r->Uniform(0., xTPC);
    double y = r->Uniform(-yTPC/2., yTPC/2.);
    double z = r->Uniform(0., zTPC);
    
    TVector3 det;
    det.SetXYZ(x,y,z);
   
    delete d;
    delete r;
 
    return det;
}


//___________________________________________________________________________
TVector3 NuMIFlux::FromDetToBeam( const TVector3& det ) {
    
    TVector3 beam;
    TRotation R;

    //corrected rotation matrix using the 0,0,0 position for MicroBooNE
    //Previous matrix is calculated relative to MiniBooNE, which is not in the centre of the BNB!

    TVector3 newX(0.92103853804025682, 0.0000462540012621546684, -0.38947144863934974);
    TVector3 newY(0.0227135048039241207, 0.99829162468141475, 0.0538324139386641073);
    TVector3 newZ(0.38880857519374290, -0.0584279894529063024, 0.91946400794392302);
    //old matrix
    /*
    TVector3 newX(0.921228671,   0.00136256111, -0.389019125);
    TVector3 newY(0.0226872648,  0.998103714,    0.0572211871);
    TVector3 newZ(0.388359401,  -0.061539578,    0.919450845);
    */
 
    R.RotateAxes(newX,newY,newZ);
    if (debug) {
        cout << "R_{beam to det} = " << endl;
        cout << " [ " << R.XX() << " " << R.XY() << " " << R.XZ() << " ] " << endl;
        cout << " [ " << R.YX() << " " << R.YY() << " " << R.YZ() << " ] " << endl;
        cout << " [ " << R.ZX() << " " << R.ZY() << " " << R.ZZ() << " ] " << endl;
        cout << endl;
    }
    R.Invert(); // R is now the inverse
    if (debug) {
        cout << "R_{det to beam} = " << endl;
        cout << " [ " << R.XX() << " " << R.XY() << " " << R.XZ() << " ] " << endl;
        cout << " [ " << R.YX() << " " << R.YY() << " " << R.YZ() << " ] " << endl;
        cout << " [ " << R.ZX() << " " << R.ZY() << " " << R.ZZ() << " ] " << endl;
        cout << endl;
    }
    // Now R allows to go from detector to beam coordinates.
    // NuMIDet is vector from NuMI target to uB detector (in beam coordinates)
    // Updated position - leaving old positions here (July 2018)
    //TVector3 NuMIDet (54.499, 74.461,  677.611); // m
    TVector3 NuMIDet (55.02, 72.59,  672.70); //m
    NuMIDet *= 100.; // To have NuMIDet in cm
    
    beam = R * det + NuMIDet;
    
    return beam;
}



//___________________________________________________________________________
double NuMIFlux::estimate_pots(int highest_potnum) {

  // Stolen: https://cdcvs.fnal.gov/redmine/projects/dk2nu/repository/show/trunk/dk2nu
  // looks like low counts are due to "evtno" not including                     
  // protons that miss the actual target (hit baffle, etc)
  // this will vary with beam conditions parameters
  // so we should round way up, those generating flugg files
  // aren't going to quantize less than 1000
  // though 500 would probably be okay, 100 would not.
  // also can't use _last_ potnum because muons decay->> don't
  // have theirs set

  // Marco: Trying with 10000
  const Int_t    nquant = 10000; //1000; // 500;  // 100                                 
  const Double_t rquant = nquant;

  Int_t estimate = (TMath::FloorNint((highest_potnum-1)/rquant)+1)*nquant;
  return estimate;
}



//___________________________________________________________________________
int NuMIFlux::calcEnuWgt( FluxNtuple* decay, const TVector3& xyz,
                         double& enu, double& wgt_xy)
{

    // Stolen: https://cdcvs.fnal.gov/redmine/projects/dk2nu/repository/show/trunk/dk2nu
    // Neutrino Energy and Weight at arbitrary point
    // Based on:
    //   NuMI-NOTE-BEAM-0109 (MINOS DocDB # 109)
    //   Title:   Neutrino Beam Simulation using PAW with Weighted Monte Carlos
    //   Author:  Rick Milburn
    //   Date:    1995-10-01
    
    // History:
    // jzh  3/21/96 grab R.H.Milburn's weighing routine
    // jzh  5/ 9/96 substantially modify the weighting function use dot product
    //              instead of rotation vecs to get theta get all info except
    //              det from ADAMO banks neutrino parent is in Particle.inc
    //              Add weighting factor for polarized muon decay
    // jzh  4/17/97 convert more code to double precision because of problems
    //              with Enu>30 GeV
    // rwh 10/ 9/08 transliterate function from f77 to C++
    
    // Original function description:
    //   Real function for use with PAW Ntuple To transform from destination
    //   detector geometry to the unit sphere moving with decayng hadron with
    //   velocity v, BETA=v/c, etc..  For (pseudo)scalar hadrons the decay will
    //   be isotropic in this  sphere so the fractional area (out of 4-pi) is the
    //   fraction of decay that hit the target.  For a given target point and
    //   area, and given x-y components of decay transverse location and slope,
    //   and given decay distance from target ans given decay GAMMA and
    //   rest-frame neutrino energy, the lab energy at the target and the
    //   fractional solid angle in the rest-frame are determined.
    //   For muon decay, correction for non-isotropic nature of decay is done.
    
    // Arguments:
    //    dk2nu    :: contains current decay information
    //    xyz      :: 3-vector of position to evaluate
    //                in *beam* frame coordinates  (cm units)
    //    enu      :: resulting energy
    //    wgt_xy   :: resulting weight
    // Return:
    //    (int)    :: error code
    // Assumptions:
    //    Energies given in GeV
    //    Particle codes have been translated from GEANT into PDG codes
    
    // for now ... these masses _should_ come from TDatabasePDG
    // but use these hard-coded values to "exactly" reproduce old code
    //
    const double kPIMASS = 0.13957;
    const double kKMASS  = 0.49368;
    const double kK0MASS = 0.49767;
    const double kMUMASS = 0.105658389;
    const double kOMEGAMASS = 1.67245;
    
    /*
     const int kpdg_nue       =   12;  // extended Geant 53
     const int kpdg_nuebar    =  -12;  // extended Geant 52
     const int kpdg_numu      =   14;  // extended Geant 56
     const int kpdg_numubar   =  -14;  // extended Geant 55
     
     const int kpdg_muplus     =   -13;  // Geant  5
     const int kpdg_muminus    =    13;  // Geant  6
     const int kpdg_pionplus   =   211;  // Geant  8
     const int kpdg_pionminus  =  -211;  // Geant  9
     const int kpdg_k0long     =   130;  // Geant 10  ( K0=311, K0S=310 )
     const int kpdg_k0short    =   310;  // Geant 16
     const int kpdg_k0mix      =   311;
     const int kpdg_kaonplus   =   321;  // Geant 11
     const int kpdg_kaonminus  =  -321;  // Geant 12
     const int kpdg_omegaminus =  3334;  // Geant 24
     const int kpdg_omegaplus  = -3334;  // Geant 32
     */
    
    // Marco: redefine (hopefully just for now)
    
    const int kpdg_nue       =   53;  // extended Geant 53
    const int kpdg_nuebar    =  52;  // extended Geant 52
    const int kpdg_numu      =   56;  // extended Geant 56
    const int kpdg_numubar   =  55;  // extended Geant 55
    
    const int kpdg_muplus     =   5;  // Geant  5
    const int kpdg_muminus    =    6;  // Geant  6
    const int kpdg_pionplus   =   8;  // Geant  8
    const int kpdg_pionminus  =  9;  // Geant  9
    const int kpdg_k0long     =   10;  // Geant 10  ( K0=311, K0S=310 )
    const int kpdg_k0short    =   16;  // Geant 16
    const int kpdg_k0mix      =   311;
    const int kpdg_kaonplus   =   11;  // Geant 11
    const int kpdg_kaonminus  =  12;  // Geant 12
    const int kpdg_omegaminus =  24;  // Geant 24
    const int kpdg_omegaplus  = 32;  // Geant 32
    
    // Marco: end
    
    
    
    const double kRDET = 100.0;   // set to flux per 100 cm radius
    
    double xpos = xyz.X();
    double ypos = xyz.Y();
    double zpos = xyz.Z();
    
    enu    = 0.0;  // don't know what the final value is
    wgt_xy = 0.0;  // but set these in case we return early due to error
    
    
    // in principle we should get these from the particle DB
    // but for consistency testing use the hardcoded values
    double parent_mass = kPIMASS;
    
    /*
     if ( decay->ptype == kpdg_pionminus)  parent_mass = kPIMASS;
     if ( decay->ptype == kpdg_kaonminus)  parent_mass = kKMASS;
     if ( decay->ptype == kpdg_k0mix)      parent_mass = kK0MASS;
     if ( decay->ptype == kpdg_muminus)    parent_mass = kMUMASS;
     if ( decay->ptype == kpdg_omegaplus)  parent_mass = kOMEGAMASS;
     */
    switch ( decay->ptype ) {
        case kpdg_pionplus:
        case kpdg_pionminus:
            parent_mass = kPIMASS;
            break;
        case kpdg_kaonplus:
        case kpdg_kaonminus:
            parent_mass = kKMASS;
            break;
        case kpdg_k0long:
        case kpdg_k0short:
        case kpdg_k0mix:
            parent_mass = kK0MASS;
            break;
        case kpdg_muplus:
        case kpdg_muminus:
            parent_mass = kMUMASS;
            break;
        case kpdg_omegaminus:
        case kpdg_omegaplus:
            parent_mass = kOMEGAMASS;
            break;
        default:
            std::cerr << "bsim::calcEnuWgt unknown particle type " << decay->ptype
            << std::endl << std::flush;
            assert(0);
            return 1;
    }
    
    
    
    
    double parentp2 = ( decay->pdPx*decay->pdPx +
                       decay->pdPy*decay->pdPy +
                       decay->pdPz*decay->pdPz );
    double parent_energy = TMath::Sqrt( parentp2 +
                                       parent_mass*parent_mass);
    double parentp = TMath::Sqrt( parentp2 );
    
    double gamma     = parent_energy / parent_mass;
    double gamma_sqr = gamma * gamma;
    double beta_mag  = TMath::Sqrt( ( gamma_sqr - 1.0 )/gamma_sqr );
    
    // Get the neutrino energy in the parent decay CM
    double enuzr = decay->Necm;
    // Get angle from parent line of flight to chosen point in beam frame
    double rad = TMath::Sqrt( (xpos-decay->Vx)*(xpos-decay->Vx) +
                             (ypos-decay->Vy)*(ypos-decay->Vy) +
                             (zpos-decay->Vz)*(zpos-decay->Vz) );
    
    double emrat = 1.0;
    double costh_pardet = -999. ; //, theta_pardet = -999.;
    
    // boost correction, but only if parent hasn't stopped
    if ( parentp > 0. ) {
        costh_pardet = ( decay->pdPx*(xpos-decay->Vx) +
                        decay->pdPy*(ypos-decay->Vy) +
                        decay->pdPz*(zpos-decay->Vz) )
        / ( parentp * rad);
        if ( costh_pardet >  1.0 ) costh_pardet =  1.0;
        if ( costh_pardet < -1.0 ) costh_pardet = -1.0;
      //  theta_pardet = TMath::ACos(costh_pardet);
        
        // Weighted neutrino energy in beam, approx, good for small theta
        emrat = 1.0 / ( gamma * ( 1.0 - beta_mag * costh_pardet ));
    }

    enu = emrat * enuzr;  // the energy ... normally
    
    
    
    
    // Get solid angle/4pi for detector element
    // small angle approximation, fixed by Alex Radovic
    //SAA//  double sangdet = ( kRDET*kRDET /
    //SAA//                   ( (zpos-decay->Vz)*(zpos-decay->Vz) ) ) / 4.0;
    double sanddetcomp = TMath::Sqrt( ( (xpos-decay->Vx)*(xpos-decay->Vx) ) +
                                     ( (ypos-decay->Vy)*(ypos-decay->Vy) ) +
                                     ( (zpos-decay->Vz)*(zpos-decay->Vz) )   );
    double sangdet = (1.0-TMath::Cos(TMath::ATan( kRDET / sanddetcomp )))/2.0;
    
    // Weight for solid angle and lorentz boost
    wgt_xy = sangdet * ( emrat * emrat );  // ! the weight ... normally
    
    // Done for all except polarized muon decay
    // in which case need to modify weight
    // (must be done in double precision)
    if ( decay->ptype == kpdg_muplus || decay->ptype == kpdg_muminus) {
        double beta[3], p_dcm_nu[4], p_nu[3], p_pcm_mp[3], partial;
        
        // Boost neu neutrino to mu decay CM
        beta[0] = decay->pdPx / parent_energy;
        beta[1] = decay->pdPy / parent_energy;
        beta[2] = decay->pdPz / parent_energy;
        p_nu[0] = (xpos-decay->Vx)*enu/rad;
        p_nu[1] = (ypos-decay->Vy)*enu/rad;
        p_nu[2] = (zpos-decay->Vz)*enu/rad;
        partial = gamma *
        (beta[0]*p_nu[0] + beta[1]*p_nu[1] + beta[2]*p_nu[2] );
        partial = enu - partial/(gamma+1.0);
        // the following calculation is numerically imprecise
        // especially p_dcm_nu[2] leads to taking the difference of numbers
        //  of order ~10's and getting results of order ~0.02's
        // for g3numi we're starting with floats (ie. good to ~1 part in 10^7)
        p_dcm_nu[0] = p_nu[0] - beta[0]*gamma*partial;
        p_dcm_nu[1] = p_nu[1] - beta[1]*gamma*partial;
        p_dcm_nu[2] = p_nu[2] - beta[2]*gamma*partial;
        p_dcm_nu[3] = TMath::Sqrt( p_dcm_nu[0]*p_dcm_nu[0] +
                                  p_dcm_nu[1]*p_dcm_nu[1] +
                                  p_dcm_nu[2]*p_dcm_nu[2] );
        
        
        
        
        // Boost parent of mu to mu production CM
        double particle_energy = decay->ppenergy;
        gamma = particle_energy/parent_mass;
        beta[0] = decay->ppdxdz * decay->pppz / particle_energy;
        beta[1] = decay->ppdydz * decay->pppz / particle_energy;
        beta[2] =                    decay->pppz / particle_energy;
        partial = gamma * ( beta[0]*decay->muparpx +
                           beta[1]*decay->muparpy +
                           beta[2]*decay->muparpz );
        partial = decay->mupare - partial/(gamma+1.0);
        p_pcm_mp[0] = decay->muparpx - beta[0]*gamma*partial;
        p_pcm_mp[1] = decay->muparpy - beta[1]*gamma*partial;
        p_pcm_mp[2] = decay->muparpz - beta[2]*gamma*partial;
        double p_pcm = TMath::Sqrt ( p_pcm_mp[0]*p_pcm_mp[0] +
                                    p_pcm_mp[1]*p_pcm_mp[1] +
                                    p_pcm_mp[2]*p_pcm_mp[2] );
        
        const double eps = 1.0e-30;  // ? what value to use
        if ( p_pcm < eps || p_dcm_nu[3] < eps ) {
            return 3; // mu missing parent info?
        }
        // Calc new decay angle w.r.t. (anti)spin direction
        double costh = ( p_dcm_nu[0]*p_pcm_mp[0] +
                        p_dcm_nu[1]*p_pcm_mp[1] +
                        p_dcm_nu[2]*p_pcm_mp[2] ) /
        ( p_dcm_nu[3]*p_pcm );
        if ( costh >  1.0 ) costh =  1.0;
        if ( costh < -1.0 ) costh = -1.0;
        // Calc relative weight due to angle difference
        double wgt_ratio = 0.0;
        /*
         if (decay->Ntype == kpdg_nuebar) wgt_ratio = 1.0 - costh;
         if (decay->Ntype == kpdg_numubar) {
         double xnu = 2.0 * enuzr / kMUMASS;
         wgt_ratio = ( (3.0-2.0*xnu )  - (1.0-2.0*xnu)*costh ) / (3.0-2.0*xnu);
         }
         */  
        switch ( decay->Ntype ) {
            case kpdg_nue:
            case kpdg_nuebar:
                wgt_ratio = 1.0 - costh;
                break;
            case kpdg_numu:
            case kpdg_numubar:
            {
                double xnu = 2.0 * enuzr / kMUMASS;
                wgt_ratio = ( (3.0-2.0*xnu )  - (1.0-2.0*xnu)*costh ) / (3.0-2.0*xnu);
                break;
            }
            default:
                return 2; // bad neutrino type
        }
        
        wgt_xy = wgt_xy * wgt_ratio;
        
    } // ptype is muon
    
    return 0;
}
