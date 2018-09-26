#script to plot the total nue flux histo and the constrained nue flux histo on the same graph.
import os,sys,string, time
import ROOT
from math import *
from ROOT import TTree, TObject, TFile, gDirectory, TH1D, TH2D, TH3D, TCanvas, gROOT, TGaxis, gStyle, TColor, TLegend, THStack, TChain, TLatex, TText
#from ROOT import *
from array import array
from glob import glob

# Importing rootlogon.C
ROOT.gROOT.SetMacroPath('~/');
ROOT.gROOT.Macro( os.path.expanduser( 'rootlogon.C' ) )

# Opening root file
#f = TFile("files/NuMIFlux.root")
f = TFile("../NuMIFlux.root")
f.ls()

h_nue = f.Get("nueFluxHisto_orig")

h_nue_con = f.Get("nueFluxHisto_orig_NA49_kaon_con")

### numu #####
h_nue.SetLineColor(ROOT.kBlack)
h_nue.SetLineWidth(3)
h_nue.GetXaxis().SetRangeUser(0,7)
h_nue.SetTitle("Electron Neutrino Flux")
h_nue.GetXaxis().SetTitle("Neutrino Energy [GeV]")
h_nue.GetYaxis().SetTitle("#Phi(#nu_{e}) / 50 MeV / cm^{2} / 6x10^{20} POT")

### numu_con #####
h_nue_con.SetLineColor(ROOT.kRed)
#kPink+6 - kaon
#kViolet+7 - pion
#kAzure-2 - proton
#kOrange +1 - neutron
h_nue_con.SetLineWidth(3)


c = TCanvas("c1","c1",1200,1200)
c.SetLogy()


h_nue.Draw()
h_nue_con.Draw("same")



leg = TLegend(.35, .75, .6, .9)

leg.SetFillStyle(0);
leg.AddEntry(h_nue,   "#nu_{e}",         "l");
leg.AddEntry(h_nue_con,   " NA49 covered K^{#pm} producing #nu_{e} ",         "l");

leg.Draw();


raw_input("Please press enter to exit.")
