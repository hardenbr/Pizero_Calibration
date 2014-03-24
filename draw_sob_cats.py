import ROOT as rt
import sys
import rootlogon

rootlogon.style()
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptTitle(0)
min_s = 0
max_s = 2.5
nbin_s = 25

min_sob = 0
max_sob = 2
nbin_sob = 50

file = rt.TFile(sys.argv[1])
cat = int(sys.argv[2])

tree = file.Get("sob_tree")


hist = rt.TH2D("hist","hist",nbin_s, min_s, max_s, nbin_sob, min_sob, max_sob)

xaxis = hist.GetXaxis()
xaxis.SetTitle("A.U. proportional to number of #pi^{0}")

yaxis = hist.GetYaxis()
yaxis.SetTitle("S/B in 2#sigma window")

tree.Draw(">>iterlist","","entrylist")
itlist = rt.gDirectory.Get("iterlist")

#loop over the tree and parse the limits                                                                                                                            
for event in range(tree.GetEntries()):
    entry = itlist.Next()
    tree.GetEntry(entry)

    bin = hist.FindBin(tree.ns/10000.,tree.sob)

    content = hist.GetBinContent(bin)

    if content == 0 and tree.cat == cat:
        hist.Fill(tree.ns/10000.,tree.sob,tree.gid)
    
hist.Draw("colztext")

raw_input("RAW INPUT")
