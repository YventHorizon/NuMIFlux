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

h_nue_Ndecay_one = f.Get("nue_Ndecay_one")
h_nue_Ndecay_six = f.Get("nue_Ndecay_six")
h_nue_Ndecay_eleven = f.Get("nue_Ndecay_eleven")
#h_nue_Ndecay_total = f.Get("nueFluxHisto")

ch_nue_Ndecay = TCanvas("c1","c1",1200,1000)
ch_nue_Ndecay.SetLogy()

h_nue_Ndecay_six.SetLineColor(ROOT.kMagenta-3)
h_nue_Ndecay_six.SetLineWidth(3)
h_nue_Ndecay_six.GetXaxis().SetRangeUser(0,7)
h_nue_Ndecay_six.SetTitle("Electron Neutrino Flux by Neutrino Decay Type")
h_nue_Ndecay_six.GetXaxis().SetTitle("Neutrino Energy [GeV]")
h_nue_Ndecay_six.GetYaxis().SetTitle("#Phi(#nu) / 50 MeV / cm^{2} / 6x10^{20} POT")

h_nue_Ndecay_one.SetLineColor(ROOT.kBlue-3)
h_nue_Ndecay_one.SetLineWidth(3)
#h_nue_Ndecay_one.SetFillColor(ROOT.kBlue-3)


h_nue_Ndecay_eleven.SetLineColor(ROOT.kCyan-3)
h_nue_Ndecay_eleven.SetLineWidth(3)
#h_nue_Ndecay_eleven.SetFillColor(ROOT.kCyan-3)
#h_nue_Ndecay_eleven.SetFillStyle(3006)

#h_nue_Ndecay_total.SetLineColor(ROOT.kBlack);

h_nue_Ndecay_six.Draw();
h_nue_Ndecay_one.Draw("same");
h_nue_Ndecay_eleven.Draw("same");
#h_nue_Ndecay_total.Draw("same");

leg = TLegend(.5, .45, 0.8, .70)
leg.SetFillStyle(0);
#gStyle.SetLegendTextSize(2/30.);
leg.AddEntry(h_nue_Ndecay_six,   "K^{+} -> #nu_{e} #pi^{0} e^{+}",         "l");
leg.AddEntry(h_nue_Ndecay_one,   "K^{0}_{L} -> #nu_{e} #pi^{-} e^{+}",        "l");
leg.AddEntry(h_nue_Ndecay_eleven,   "#mu^{+} -> #bar{#nu}_{#mu} #nu_{e} e^{+}",         "l");
#leg.AddEntry(h_nue_Ndecay_total,   "total flux",         "l");
leg.Draw();

t = TLatex(.45, .625, "#nu_{e}");
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
