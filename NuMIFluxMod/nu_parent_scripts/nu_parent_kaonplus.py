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
h_nue_parent_kaonplus = f.Get("nue_parent_kaonplus")
h_numu_parent_kaonplus = f.Get("numu_parent_kaonplus")

ch_nue_parent_kaonplus = TCanvas("c1","c1",1200,1200)
ch_nue_parent_kaonplus.SetLogy()

h_numu_parent_kaonplus.SetLineColor(ROOT.kOrange+8)
h_numu_parent_kaonplus.SetLineWidth(3)
h_numu_parent_kaonplus.SetTitle("Neutrinos from K^{+}")
#h_nue_parent_muonplus.SetFillColor(ROOT.kOrange+10)
h_numu_parent_kaonplus.GetXaxis().SetRangeUser(0,8)
h_numu_parent_kaonplus.GetXaxis().SetTitle("Neutrino Energy [GeV]")
h_numu_parent_kaonplus.GetYaxis().SetTitle("#Phi(#nu) / 50 MeV / cm^{2} / 6x10^{20} POT")

h_nue_parent_kaonplus.SetLineColor(ROOT.kBlue - 4)
h_nue_parent_kaonplus.SetLineWidth(3)
#h_nue_parent_kaonplus.SetFillColor(ROOT.kOrange+8)
#h_nue_parent_kaonplus.SetFillStyle(3144)

h_numu_parent_kaonplus.Draw();
h_nue_parent_kaonplus.Draw("same");
#h_nue_parent_totall.Draw("Same");

leg = TLegend(.5, .55, .7, .70)
leg.SetFillStyle(0);
#gStyle.SetLegendTextSize(2/30.);
leg.AddEntry(h_nue_parent_kaonplus,   "#nu_{e}",        "l");
leg.AddEntry(h_numu_parent_kaonplus,   "#nu_{#mu}",        "l");
leg.Draw();

#t = TLatex(.4, .625, "#nu_{e}");
#t.SetTextColor(ROOT.kBlack);
#t.SetNDC();
#t.SetTextSize(2/30.);
#t.SetTextAlign(32);
#t.Draw();
raw_input("Please press enter to exit.")
