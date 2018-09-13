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

#h_nue_tgptype_total = f.Get("nueFluxHisto")
h_nue_tgptype_pionplus = f.Get("nue_tgptype_pionplus")
h_nue_tgptype_pionminus = f.Get("nue_tgptype_pionminus")
h_nue_tgptype_kaonlong = f.Get("nue_tgptype_kaonlong")
h_nue_tgptype_kaonplus = f.Get("nue_tgptype_kaonplus")
h_nue_tgptype_kaonminus = f.Get("nue_tgptype_kaonminus")
h_nue_tgptype_neutron = f.Get("nue_tgptype_neutron")
h_nue_tgptype_proton = f.Get("nue_tgptype_proton")
h_nue_tgptype_lambda = f.Get("nue_tgptype_lambda")
h_nue_tgptype_sigmaplus = f.Get("nue_tgptype_sigmaplus")

ch_nue_Ndecay = TCanvas("c1","c1",1500,1200)
ch_nue_Ndecay.SetLogy()

#h_nue_tgptype_total.SetLineColor(ROOT.kBlack)
#h_nue_Ndecay_thirteen.SetFillColor(ROOT.kBlue-3)
#h_nue_tgptype_total.GetXaxis().SetRangeUser(0,14)
#h_nue_tgptype_total.GetXaxis().SetTitle("Neutrino Energy [GeV]")
#h_nue_tgptype_total.GetYaxis().SetTitle("#Phi(#nu) / 50 MeV / cm^{2} / 6x10^{20} POT")

h_nue_tgptype_proton.SetLineColor(ROOT.kCyan - 7)
h_nue_tgptype_proton.GetXaxis().SetRangeUser(0,6)
h_nue_tgptype_proton.GetXaxis().SetTitle("Neutrino Energy [GeV]")
h_nue_tgptype_proton.GetYaxis().SetTitle("#Phi(#nu) / 50 MeV / cm^{2} / 6x10^{20} POT")
h_nue_tgptype_proton.SetLineWidth(3)
#h_nue_Ndecay_five.SetFillColor(ROOT.kGreen-3)
#h_nue_Ndecay_five.SetFillStyle(3144)

h_nue_tgptype_pionminus.SetLineColor(ROOT.kGreen - 7)
h_nue_tgptype_pionminus.SetLineWidth(3)
#h_nue_Ndecay_seven.SetFillColor(ROOT.kMagenta-3)
#h_nue_Ndecay_seven.SetFillStyle(3006)


h_nue_tgptype_kaonlong.SetLineColor(ROOT.kOrange-2)
h_nue_tgptype_kaonlong.SetLineWidth(3)
#h_nue_Ndecay_three.SetFillColor(ROOT.kCyan-3)
#h_nue_Ndecay_three.SetFillStyle(3006)

h_nue_tgptype_kaonplus.SetLineColor(ROOT.kOrange-3)
h_nue_tgptype_kaonplus.SetLineWidth(3)
#h_nue_Ndecay_twelve.SetFillColor(ROOT.kRed-3)
#h_nue_Ndecay_twelve.SetFillStyle(3006)

h_nue_tgptype_kaonminus.SetLineColor(ROOT.kRed - 7)
h_nue_tgptype_kaonminus.SetLineWidth(3)

h_nue_tgptype_neutron.SetLineColor(ROOT.kViolet - 1)
h_nue_tgptype_neutron.SetLineWidth(3)

h_nue_tgptype_pionplus.SetLineColor(ROOT.kBlue)
h_nue_tgptype_pionplus.SetLineWidth(3)

h_nue_tgptype_lambda.SetLineColor(ROOT.kTeal - 1)
h_nue_tgptype_lambda.SetLineWidth(3)

h_nue_tgptype_sigmaplus.SetLineColor(ROOT.kBlue + 7)
h_nue_tgptype_sigmaplus.SetLineWidth(3)

#h_nue_tgptype_total.Draw();
h_nue_tgptype_proton.Draw();
h_nue_tgptype_pionplus.Draw("same");
h_nue_tgptype_pionminus.Draw("same");
h_nue_tgptype_kaonlong.Draw("same");
h_nue_tgptype_kaonplus.Draw("same");
h_nue_tgptype_kaonminus.Draw("same");
h_nue_tgptype_neutron.Draw("same");
h_nue_tgptype_lambda.Draw("same");
h_nue_tgptype_sigmaplus.Draw("same");

leg = TLegend(.65, .45, .73, .75)  # x, y, x, y
leg.SetFillStyle(0);
#gStyle.SetLegendTextSize(2/30.);
leg.AddEntry(h_nue_tgptype_pionplus,   "#pi^{+}",        "l");
leg.AddEntry(h_nue_tgptype_pionminus,  "#pi^{-}",       "l");
leg.AddEntry(h_nue_tgptype_kaonlong,   "K^{0}_{L}",         "l");
leg.AddEntry(h_nue_tgptype_kaonplus,   "K^{+}",         "l");
leg.AddEntry(h_nue_tgptype_kaonminus,   "K^{-}",         "l");
leg.AddEntry(h_nue_tgptype_neutron,   "n",         "l");
leg.AddEntry(h_nue_tgptype_proton,   "p",         "l");
leg.AddEntry(h_nue_tgptype_lambda,   "#lambda",         "l");
leg.AddEntry(h_nue_tgptype_sigmaplus,   "#Sigma^{+}",         "l");
leg.Draw();

t = TLatex(.6, .60, "#nu_{e}");  # x, y
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
