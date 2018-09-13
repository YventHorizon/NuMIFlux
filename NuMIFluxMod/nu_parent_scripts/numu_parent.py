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

#h_numu_parent_totall = f.Get("numuFluxHisto")
h_numu_parent_pionplus = f.Get("numu_parent_pionplus")
h_numu_parent_kaonplus = f.Get("numu_parent_kaonplus")
h_numu_parent_kaonlong = f.Get("numu_parent_kaonlong")
h_numu_parent_muonminus = f.Get("numu_parent_muonminus")


ch_numu_parent = TCanvas("c1","c1",1200,1000)
ch_numu_parent.SetLogy()

h_numu_parent_pionplus.SetLineColor(ROOT.kAzure+8)
h_numu_parent_pionplus.SetLineWidth(3)
h_numu_parent_pionplus.SetTitle("Muon Neutrino Flux by Parent Type")
#h_numu_parent_pionplus.SetFillColor(ROOT.kViolet+4)
h_numu_parent_pionplus.GetXaxis().SetRangeUser(0,8)
h_numu_parent_pionplus.GetXaxis().SetTitle("Neutrino Energy [GeV]")
h_numu_parent_pionplus.GetYaxis().SetTitle("#Phi(#nu) / 50 MeV / cm^{2} / 6x10^{20} POT")

h_numu_parent_kaonplus.SetLineColor(ROOT.kOrange+6)
h_numu_parent_kaonplus.SetLineWidth(3)
#h_numu_parent_kaonplus.SetFillColor(ROOT.kViolet+6)
#h_numu_parent_kaonplus.SetFillStyle(3144)

h_numu_parent_kaonlong.SetLineColor(ROOT.kViolet+2)
h_numu_parent_kaonlong.SetLineWidth(3)
#h_numu_parent_kaonlong.SetFillColor(ROOT.kViolet+2)
#h_numu_parent_kaonlong.SetFillStyle(3006)

h_numu_parent_muonminus.SetLineColor(ROOT.kMagenta-7)
h_numu_parent_muonminus.SetLineWidth(3)
#h_numu_parent_muonminus.SetFillColor(ROOT.kMagenta-10)
#h_numu_parent_muonminus.SetFillStyle(3006)

#h_numu_parent_totall.SetLineColor(ROOT.kBlack);

h_numu_parent_pionplus.Draw();
h_numu_parent_kaonplus.Draw("same");
h_numu_parent_kaonlong.Draw("same");
h_numu_parent_muonminus.Draw("same");
#h_numu_parent_totall.Draw("same");

leg = TLegend(.5, .55, .7, .70)
leg.SetFillStyle(0);
leg.AddEntry(h_numu_parent_pionplus,   "#pi^{+}",        "l");
leg.AddEntry(h_numu_parent_kaonplus,  "K^{+}",       "l");
leg.AddEntry(h_numu_parent_kaonlong,   "K^{0}_{L}",         "l");
leg.AddEntry(h_numu_parent_muonminus,   "#mu^{-}",         "l");
#leg.AddEntry(h_numu_parent_totall,   "total #nu_{#mu} flux",         "l");
leg.Draw();

t = TLatex(.4, .625, "#nu_{#mu}");
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
