import os,sys,string, time
import ROOT
from math import *
from ROOT import TTree, TObject, TFile, gDirectory, TH1D, TH2D, TH3D, TCanvas, gROOT, TGaxis, gStyle, TColor, TLegend, THStack, TChain, TLatex, TText, TLine
#from ROOT import *
from array import array
from glob import glob

# Importing rootlogon.C
ROOT.gROOT.SetMacroPath('~/');
ROOT.gROOT.Macro( os.path.expanduser( 'rootlogon.C' ) )

# set kBird palette
red   = [ 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764];
green = [ 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832];
blue = [ 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539];
stops = [ 0.0, 0.1, 0.2, 0.3, 0.4 ,0.5 ,0.6 ,0.7, 1.0 ];
stopsArray = array("d", stops)
redArray = array("d", red)
greenArray = array("d", green)
blueArray = array("d", blue)
ROOT.TColor.CreateGradientColorTable(9, stopsArray, redArray, greenArray, blueArray, 255)

# Opening root file
#f = TFile("files/NuMIFlux.root")
f = TFile("../NuMIFlux.root")
f.ls()

h_neutron_long_trans = f.Get("neutron_NA49_nue")

c4 = TCanvas("c1","c1",1500,1100)

h_neutron_long_trans.SetTitle("#nu_{e} s from n")
h_neutron_long_trans.GetXaxis().SetTitle("x_{F} ")
h_neutron_long_trans.GetXaxis().SetTitleFont(12)
h_neutron_long_trans.GetYaxis().SetTitle("p_{T} [GeV/c]")
h_neutron_long_trans.GetYaxis().SetTitleFont(12)
h_neutron_long_trans.GetXaxis().SetRangeUser(-0.11,1)
h_neutron_long_trans.GetYaxis().SetRangeUser(0,3.7)
h_neutron_long_trans.Draw("colz")

# Eur Phys J (2013) neutron-carbon invariant cross section measurements with pT, xF bins. (NA49, thin target, @158 GeV,  so with tgptype particle)

line1 = TLine( 0.05 ,   0,  1  ,  0)
line2 = TLine( 1  ,  0,  1   , 2.4)
line3 = TLine( 1   , 2.4,   0.475  ,  2.4)
line4 = TLine( 0.475  ,  2.4,   0.05  ,  0.6)
line5 = TLine( 0.05  ,  0.6,    0.05  ,  0)

lineVector= [ line1, line2, line3, line4, line5]
for line in lineVector:
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(3)
    line.Draw()

leg = TLegend(.54, .75 , .58, .79)
leg.SetFillStyle(0);
leg.SetLineWidth(4);
leg.Draw();


text = TLatex(.6, .75, "NA49 Coverage");
text.SetTextColor(ROOT.kBlack);
text.SetNDC();
text.SetTextSize(1.4/30.);
text.SetTextAlign(11);
#text.DrawLatex(.48, .55, "#Square");
text.Draw();

raw_input("Please press enter to exit.")
