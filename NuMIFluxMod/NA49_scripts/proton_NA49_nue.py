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

h_proton_long_trans = f.Get("proton_NA49_nue")

c4 = TCanvas("c1","c1",1500,1100)

h_proton_long_trans.SetTitle("#nu_{e} s from p")
h_proton_long_trans.GetXaxis().SetTitle("x_{F} ")
h_proton_long_trans.GetXaxis().SetTitleFont(12)
h_proton_long_trans.GetYaxis().SetTitle("p_{T} [GeV/c]")
h_proton_long_trans.GetYaxis().SetTitleFont(12)
h_proton_long_trans.GetXaxis().SetRangeUser(-1,1)
h_proton_long_trans.GetYaxis().SetRangeUser(0,4)
h_proton_long_trans.Draw("colz")

# Eur Phys J (2013) proton-carbon invariant cross section measurements with pT, xF bins. (NA49, thin target, @158 GeV,  so with tgptype particle)

#box 1
line1 = TLine(-0.0125  ,  0.05, -0.0125  ,  0.95)
line2 = TLine( -0.0125  ,  0.95,  -0.575  ,  0.95 )
line3 = TLine( -0.575  ,  0.95 , -0.575   , 0.75)
line4 = TLine( -0.575   , 0.75,  -0.625  ,  0.75)
line5 = TLine( -0.625  ,  0.75  , -0.625  ,  0.65)
line6 = TLine( -0.625  ,  0.65, -0.675  ,  0.65)
line7 = TLine( -0.675  ,  0.65, -0.675  ,  0.45)
line8 = TLine( -0.675  ,  0.45  ,-0.725  ,  0.45)
line9 = TLine( -0.725  ,  0.45, -0.725  ,  0.325)
line10 = TLine( -0.725  ,  0.325,    -0.775 ,   0.325)
line11 = TLine( -0.775 ,   0.325,    -0.775 ,   0.225)
line12 = TLine( -0.775 ,   0.225,    -0.825  ,  0.225)
line13 = TLine( -0.825  ,  0.225,    -0.825  ,  0.1)
line14 = TLine( -0.825  ,  0.1,  -0.725  ,  0.1)
line15 = TLine( -0.725  ,  0.1,  -0.725  ,  0.0125)
line16 = TLine( -0.725  ,  0.0125,   -0.475  ,  0.0125)
line17 = TLine( -0.475  ,  0.0125,   -0.475  ,  0.05)
line18 = TLine( -0.475  ,  0.05, -0.0125  ,  0.05)


#box 2
line19 = TLine( -0.025  ,  1,   -0.025 ,   2)
line20 = TLine( -0.025 ,   2,   -0.15 ,   2)
line21 = TLine( -0.15 ,   2,    -0.15  ,  1.8)
line22 = TLine(  -0.15  ,  1.8, -0.225  ,  1.8)
line23 = TLine( -0.225  ,  1.8, -0.225  ,  1.6)
line24 = TLine( -0.225  ,  1.6, -0.25  ,  1.6)
line25 = TLine( -0.25  ,  1.6, -0.25  ,  2)
line26 = TLine( -0.25  ,  2,     -0.35 ,   2)
line27 = TLine( -0.35 ,   2, -0.35  ,  1.6)
line28 = TLine( -0.35  ,  1.6,   -0.45  ,  1.6)
line29 = TLine(  -0.45  ,  1.6,  -0.45  ,  2)
line30 = TLine(  -0.45   , 2, -0.55 ,   2)
line31 = TLine(  -0.55 ,   2,    -0.55  ,  1.6)
line32 = TLine( -0.55  ,  1.6,   -0.525  ,  1.6)
line33 = TLine( -0.525  ,  1.6, -0.525   , 1)
line34 = TLine( -0.525   , 1,    -0.025 ,   1)

#box 3
line35 = TLine( -0.0125  ,  0,   0.975  ,  0)
line36 = TLine( 0.975  ,  0,    0.975   , 1.925)
line37 = TLine( 0.975   , 1.925,    -0.025  ,  1.925)
line38 = TLine( -0.025  ,  1.925,   -0.0125  ,  0)

lineVector= [ line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line13, line14, line15, line16, line17, line18, line19, line20, line21,line22,line23,line24,line25,line26,line27,line28,line29,line30,line31,line32,line33,line34, line35, line36, line37, line38]
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
