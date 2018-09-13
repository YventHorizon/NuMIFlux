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

#h_numu_tgptype_total = f.Get("numuFluxHisto")
h_numu_tgptype_pionplus = f.Get("numu_tgptype_pionplus")
h_numu_tgptype_pionminus = f.Get("numu_tgptype_pionminus")
h_numu_tgptype_kaonlong = f.Get("numu_tgptype_kaonlong")
h_numu_tgptype_kaonplus = f.Get("numu_tgptype_kaonplus")
h_numu_tgptype_kaonminus = f.Get("numu_tgptype_kaonminus")
h_numu_tgptype_neutron = f.Get("numu_tgptype_neutron")
h_numu_tgptype_proton = f.Get("numu_tgptype_proton")
h_numu_tgptype_lambda = f.Get("numu_tgptype_lambda")
h_numu_tgptype_sigmaplus = f.Get("numu_tgptype_sigmaplus")

ch_numu_Ndecay = TCanvas("c1","c1",1000,1000)
ch_numu_Ndecay.SetLogy()

#h_numu_tgptype_total.SetLineColor(ROOT.kBlack)
#h_numu_Ndecay_thirteen.SetFillColor(ROOT.kBlue-3)
#h_numu_tgptype_total.GetXaxis().SetRangeUser(0,14)
#h_numu_tgptype_total.GetXaxis().SetTitle("Neutrino Energy [GeV]")
#h_numu_tgptype_total.GetYaxis().SetTitle("#Phi(#nu) / 50 MeV / cm^{2} / 6x10^{20} POT")

h_numu_tgptype_pionplus.SetLineColor(ROOT.kCyan - 7)
h_numu_tgptype_pionplus.GetXaxis().SetRangeUser(0,8)
h_numu_tgptype_pionplus.GetXaxis().SetTitle("Neutrino Energy [GeV]")
h_numu_tgptype_pionplus.GetYaxis().SetTitle("#Phi(#nu) / 50 MeV / cm^{2} / 6x10^{20} POT")
h_numu_tgptype_pionplus.SetLineWidth(3)
#h_numu_Ndecay_five.SetFillColor(ROOT.kGreen-3)
#h_numu_Ndecay_five.SetFillStyle(3144)

h_numu_tgptype_pionminus.SetLineColor(ROOT.kGreen - 7)
h_numu_tgptype_pionminus.SetLineWidth(3)
#h_numu_Ndecay_seven.SetFillColor(ROOT.kMagenta-3)
#h_numu_Ndecay_seven.SetFillStyle(3006)


h_numu_tgptype_kaonlong.SetLineColor(ROOT.kOrange-2)
h_numu_tgptype_kaonlong.SetLineWidth(3)
#h_numu_Ndecay_three.SetFillColor(ROOT.kCyan-3)
#h_numu_Ndecay_three.SetFillStyle(3006)

h_numu_tgptype_kaonplus.SetLineColor(ROOT.kOrange-3)
h_numu_tgptype_kaonplus.SetLineWidth(3)
#h_numu_Ndecay_twelve.SetFillColor(ROOT.kRed-3)
#h_numu_Ndecay_twelve.SetFillStyle(3006)

h_numu_tgptype_kaonminus.SetLineColor(ROOT.kRed - 7)
h_numu_tgptype_kaonminus.SetLineWidth(3)

h_numu_tgptype_neutron.SetLineColor(ROOT.kViolet - 1)
h_numu_tgptype_neutron.SetLineWidth(3)

h_numu_tgptype_proton.SetLineColor(ROOT.kBlue)
h_numu_tgptype_proton.SetLineWidth(3)

h_numu_tgptype_lambda.SetLineColor(ROOT.kTeal - 1)
h_numu_tgptype_lambda.SetLineWidth(3)

h_numu_tgptype_sigmaplus.SetLineColor(ROOT.kBlue + 7)
h_numu_tgptype_sigmaplus.SetLineWidth(3)

#h_numu_tgptype_total.Draw();
h_numu_tgptype_pionplus.Draw();
h_numu_tgptype_pionminus.Draw("same");
h_numu_tgptype_kaonlong.Draw("same");
h_numu_tgptype_kaonplus.Draw("same");
h_numu_tgptype_kaonminus.Draw("same");
h_numu_tgptype_neutron.Draw("same");
h_numu_tgptype_proton.Draw("same");
h_numu_tgptype_lambda.Draw("same");
h_numu_tgptype_sigmaplus.Draw("same");

leg = TLegend(.6, .55, .9, .72)
leg.SetFillStyle(0);
#gStyle.SetLegendTextSize(2/30.);
leg.AddEntry(h_numu_tgptype_pionplus,   "#pi^{+}",        "l");
leg.AddEntry(h_numu_tgptype_pionminus,  "#pi^{-}",       "l");
leg.AddEntry(h_numu_tgptype_kaonlong,   "K^{0}_{L}",         "l");
leg.AddEntry(h_numu_tgptype_kaonplus,   "K^{+}",         "l");
leg.AddEntry(h_numu_tgptype_kaonminus,   "K^{-}",         "l");
leg.AddEntry(h_numu_tgptype_neutron,   "n",         "l");
leg.AddEntry(h_numu_tgptype_proton,   "p",         "l");
leg.AddEntry(h_numu_tgptype_lambda,   "#lambda",         "l");
leg.AddEntry(h_numu_tgptype_sigmaplus,   "#Sigma^{+}",         "l");
leg.Draw();

t = TLatex(.5, .635, "#nu_{#mu}");
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
