import ROOT as rt
import sys
import rootlogon

rootlogon.style()
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptTitle(0)
min_s = 0
max_s = 25
max_s_2 = 50
nbin_s = 25

min_sob = 0
max_sob = 2
max_sob_2 = 2

nbin_sob = 25

scale_factor = 11. * 1368. / 100000000

file = rt.TFile(sys.argv[1])

tree = file.Get("sob_tree")


hist = rt.TH2D("hist","hist",nbin_s, min_s, max_s, nbin_sob, min_sob, max_sob)
hist2 = rt.TH2D("hist2","hist",nbin_s, min_s, max_s_2, nbin_sob, min_sob, max_sob_2)

xaxis = hist.GetXaxis()
xaxis.SetTitle("Available #pi^{0} Rate (kHz)")
yaxis = hist.GetYaxis()
yaxis.SetTitle("S/B in 2#sigma window")

xaxis2 = hist2.GetXaxis()
xaxis2.SetTitle("Available #pi^{0} Rate (kHz)")
yaxis2 = hist2.GetYaxis()
yaxis2.SetTitle("S/B in 2#sigma window")

tree.Draw(">>iterlist","","entrylist")
itlist = rt.gDirectory.Get("iterlist")

gid_used = []

#loop over the tree and parse the limits                                                                                                                            
for event in range(tree.GetEntries()):
    entry = itlist.Next()
    tree.GetEntry(entry)

    bin = hist.FindBin(tree.ns * scale_factor ,tree.sob)

    content = hist.GetBinContent(bin)

    if content == 0 and tree.cat == 0:
        gid_used.append(tree.gid)
        hist.Fill(tree.ns * scale_factor ,tree.sob,tree.gid)
        print tree.ns, tree.sob

tree.Draw(">>iterlist2","","entrylist")
itlist2 = rt.gDirectory.Get("iterlist2")

for event in range(tree.GetEntries()):
    entry = itlist2.Next()
    tree.GetEntry(entry)

    if tree.cat == 2: #and tree.gid in gid_used:
       
        bin = hist2.FindBin(tree.ns * scale_factor, tree.sob)
        content = hist2.GetBinContent(bin)
        


        if content == 0:
            hist2.Fill(tree.ns * scale_factor, tree.sob, tree.gid)

canvas = rt.TCanvas("c1","c1")    
hist.Draw("colztext")

canvas2 = rt.TCanvas("c2","c2")    
hist2.Draw("colztext")

raw_input("RAW INPUT")
