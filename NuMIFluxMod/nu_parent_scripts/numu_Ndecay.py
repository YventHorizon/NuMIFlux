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

h_numu_Ndecay_thirteen = f.Get("numu_Ndecay_thirteen")
h_numu_Ndecay_five = f.Get("numu_Ndecay_five")
h_numu_Ndecay_seven = f.Get("numu_Ndecay_seven")
h_numu_Ndecay_twelve = f.Get("numu_Ndecay_twelve")
#h_numu_Ndecay_total = f.Get("numuFluxHisto")

ch_numu_Ndecay = TCanvas("c1","c1",1200,1000)
ch_numu_Ndecay.SetLogy()

h_numu_Ndecay_thirteen.SetLineColor(ROOT.kBlue-3)
h_numu_Ndecay_thirteen.SetLineWidth(3)
#h_numu_Ndecay_thirteen.SetFillColor(ROOT.kBlue-3)
h_numu_Ndecay_thirteen.GetXaxis().SetRangeUser(0,8)
h_numu_Ndecay_thirteen.SetTitle("Muon Neutrino Flux by Neutrino Decay Type")
h_numu_Ndecay_thirteen.GetXaxis().SetTitle("Neutrino Energy [GeV]")
h_numu_Ndecay_thirteen.GetYaxis().SetTitle("#Phi(#nu) / 50 MeV / cm^{2} / 6x10^{20} POT")

h_numu_Ndecay_five.SetLineColor(ROOT.kGreen-3)
h_numu_Ndecay_five.SetLineWidth(3)
#h_numu_Ndecay_five.SetFillColor(ROOT.kGreen-3)
#h_numu_Ndecay_five.SetFillStyle(3144)

h_numu_Ndecay_seven.SetLineColor(ROOT.kMagenta-3)
h_numu_Ndecay_seven.SetLineWidth(3)
#h_numu_Ndecay_seven.SetFillColor(ROOT.kMagenta-3)
#h_numu_Ndecay_seven.SetFillStyle(3006)


h_numu_Ndecay_twelve.SetLineColor(ROOT.kRed-3)
h_numu_Ndecay_twelve.SetLineWidth(3)
#h_numu_Ndecay_twelve.SetFillColor(ROOT.kRed-3)
#h_numu_Ndecay_twelve.SetFillStyle(3006)

#h_numu_Ndecay_total.SetLineColor(ROOT.kBlack);

h_numu_Ndecay_thirteen.Draw();
h_numu_Ndecay_five.Draw("same");
h_numu_Ndecay_seven.Draw("same");
h_numu_Ndecay_twelve.Draw("same");
#h_numu_Ndecay_total.Draw("same");

leg = TLegend(.5, .55, 0.7, .70)   # x, y, x, y
leg.SetFillStyle(0);
#gStyle.SetLegendTextSize(2/30.);
leg.AddEntry(h_numu_Ndecay_thirteen,   "#pi^{+} -> #nu_{#mu} #mu^{+}",        "l");
leg.AddEntry(h_numu_Ndecay_five,  "K^{+} -> #nu_{#mu} #mu^{+}",       "l");
leg.AddEntry(h_numu_Ndecay_seven,   "K^{+} -> #nu_{#mu} #pi^{0} #mu^{+}",         "l");
leg.AddEntry(h_numu_Ndecay_twelve,   "#mu^{-} => #nu_{#mu} #bar{#nu}_{e} e^{-}",         "l");
#leg.AddEntry(h_numu_Ndecay_total,   "total flux",         "l");
leg.Draw();

t = TLatex(.45, .625, "#nu_{#mu}");
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
