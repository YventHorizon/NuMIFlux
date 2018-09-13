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

#h_nue_parent_totall = f.Get("nueFluxHisto")
h_nue_parent_muonplus = f.Get("nue_parent_muonplus")
h_nue_parent_kaonplus = f.Get("nue_parent_kaonplus")
h_nue_parent_kaonlong = f.Get("nue_parent_kaonlong")

ch_nue_parent = TCanvas("c1","c1",1200,1000)
ch_nue_parent.SetLogy()

h_nue_parent_kaonplus.SetLineColor(ROOT.kAzure-5)
h_nue_parent_kaonplus.SetLineWidth(3)
h_nue_parent_kaonplus.GetXaxis().SetRangeUser(0,8)
h_nue_parent_kaonplus.SetTitle("Electron Neutrino Flux by Parent Type")
h_nue_parent_kaonplus.GetXaxis().SetTitle("Neutrino Energy [GeV]")
h_nue_parent_kaonplus.GetYaxis().SetTitle("#Phi(#nu) / 50 MeV / cm^{2} / 6x10^{20} POT")


h_nue_parent_muonplus.SetLineColor(ROOT.kOrange+10)
h_nue_parent_muonplus.SetLineWidth(3)
#h_nue_parent_muonplus.SetFillColor(ROOT.kOrange+10)

#h_nue_parent_kaonplus.SetFillColor(ROOT.kOrange+8)
#h_nue_parent_kaonplus.SetFillStyle(3144)

h_nue_parent_kaonlong.SetLineColor(ROOT.kTeal+5)
h_nue_parent_kaonlong.SetLineWidth(3)
#h_nue_parent_kaonlong.SetFillColor(ROOT.kPink+5)
#h_nue_parent_kaonlong.SetFillStyle(3002)

#h_nue_parent_totall.SetLineColor(ROOT.kBlack)

h_nue_parent_kaonplus.Draw();
h_nue_parent_muonplus.Draw("same");
h_nue_parent_kaonlong.Draw("same");
#h_nue_parent_totall.Draw("Same");

leg = TLegend(.55, .5, 0.8, .70)
leg.SetFillStyle(0);
#gStyle.SetLegendTextSize(2/30.);
leg.AddEntry(h_nue_parent_muonplus,   "#mu^{+}",        "l");
leg.AddEntry(h_nue_parent_kaonplus,  "K^{+}",       "l");
leg.AddEntry(h_nue_parent_kaonlong,   "K^{0}_{L}",         "l");
#leg.AddEntry(h_nue_parent_totall,   "total #nu_{e} flux",         "l");
leg.Draw();

t = TLatex(.5, .625, "#nu_{e}");
t.SetTextColor(ROOT.kBlack);
t.SetNDC();
t.SetTextSize(2/30.);
t.SetTextAlign(32);
t.Draw();

#t2 = TLatex(.51, .48, "#splitline{Off-axis NuMI Flux}{at MicroBooNE}");
#t2.SetTextColor(ROOT.kRed+2);
#t2.SetNDC();
#t2.SetTextSize(1.4/30.);
#t2.SetTextAlign(11);
#t2.Draw();

#t3 = TLatex(.51, .40, "Anti-Neutrino Mode");
#t3.SetTextColor(ROOT.kBlack);
#t3.SetNDC();
#t3.SetTextSize(1.4/30.);
#t3.SetTextAlign(11);
#t3.Draw();


raw_input("Please press enter to exit.")
