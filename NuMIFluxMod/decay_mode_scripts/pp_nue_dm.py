#!/usr/bin/env python
#
# read a root anatree using Pyroot
#
# run with python pyrootMacroTPCActive.py 14 -1
# 14 is the nu flavour (nu-mu in this case)
# -1 means loop over all the entries

###############  decay modes  #########################
# 1:    K^{0}_{L} -> #nu_{e} #pi^{-} e^{+}
# 2:    K^{0}_{L} -> #bar{#nu}_{e} #pi^{+} e^{-}
# 3:    K^{0}_{L} -> #nu_{#mu} #pi^{-} #mu^{+}
# 4:    K^{0}_{L} -> #bar{#nu}_{#mu} #pi^{+} mu^{-}

# 5:    K^{+} -> #nu_{#mu} #mu^{+}
# 6:    K^{+} -> #nu_{e} #pi^{0} e^{+}
# 7:    K^{+} -> #nu_{#mu} #pi^{0} #mu^{+}

# 8:    K^{-} -> #bar{#nu}_{#mu} #mu^{-}
# 9:    K^{-} -> #bar{#nu}_{e} #pi^{0} e^{-}
# 10:    K^{-} -> #bar{#nu}_{#mu} #pi^{0} #mu^{-}

# 11:   #mu^{+} -> #bar{#nu}_{#mu} #nu_{e} e^{+}
# 12:   #mu^{-} => #nu_{#mu} #bar{#nu}_{e} e^{-}

# 13:   #pi^{+} -> #nu_{#mu} #mu^{+}
# 14:   #pi^{-} -> #bar{#nu}_{#mu} #mu^{-}

#######################################################

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

h_pp_nue_dm = f.Get("ppnue_dm_eleven")

cpp_nue_dm = TCanvas("c1","c1",1200,1200)

h_pp_nue_dm.SetTitle("#nu_{e} from #pi^{+}, #mu^{+} -> #bar{#nu}_{#mu} #nu_{e} e^{+}")
h_pp_nue_dm.GetXaxis().SetTitle("p_{#pi^{+}} [GeV/c]")
h_pp_nue_dm.GetXaxis().SetTitleFont(12)
h_pp_nue_dm.GetYaxis().SetTitle("#theta_{#pi^{+}} [mrad]")
h_pp_nue_dm.GetYaxis().SetTitleFont(12)
h_pp_nue_dm.GetXaxis().SetRangeUser(0,110)
h_pp_nue_dm.GetYaxis().SetRangeUser(0,500)
h_pp_nue_dm.Draw("colz")

raw_input("Please press enter to exit.")

