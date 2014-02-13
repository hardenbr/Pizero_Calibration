import ROOT as rt
import rootlogon
import sys

rootlogon.style()

rt.gStyle.SetCanvasDefH(600)
rt.gStyle.SetCanvasDefW(1800)
rt.gStyle.SetOptStat(0)
rt.gStyle.SetPadBottomMargin(0.6)
rt.gStyle.SetPadTopMargin(.05)

print "USAGE: draw_s2_s1.py [file.root] [initial l1 bit] [final l1 bit]"

file = sys.argv[1]

li = int(sys.argv[2])
lf = int(sys.argv[3])

openfile = rt.TFile(file,"READ")

s2_s1_hist = openfile.Get("s2_s1_hist")

s2_s1_hist.SetFillColor(rt.kRed);
s2_s1_hist.SetFillStyle(3005);
s2_s1_hist.SetLineColor(rt.kRed);


canvas1 = rt.TCanvas("c1")
canvas1.SetLogy(1)
canvas1.SetGridy(1)
s2_s1_hist.GetXaxis().SetRangeUser(li,lf)
s2_s1_hist.Draw()

raw_input("RAW INPUT")
