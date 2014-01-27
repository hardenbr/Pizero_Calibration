from ROOT import gSystem
from  optparse  import OptionParser
import numpy as np
import itertools, math
import ROOT as rt
import rootlogon
rootlogon.style()

def addcuts(list):
    cut = ""
    for ii in list[:-1]:
        cut = ii + " && "
    cut += list[-1]
    
    return cut

def build_grids():
    #cystals scan range
    ncri1 = range(4,9)
    ncri2 = range(4,7)

    #pizero scan range
    pi_cluster_pt = np.linspace(.4,1,4)
    pi_pizero_pt = np.linspace(1.2,2.4,5)
    pi_elayer = np.linspace(0,1.5,6)
    pi_s4s9 = np.linspace(.7,.9,5)
    pi_iso = np.linspace(.1,.3,5)
    
    #eta scan range
    eta_cluster_pt = np.linspace(.4,1,4)
    eta_eta_pt = np.linspace(1.2,2.4,5)
    eta_elayer = np.linspace(0,1.5,6)
    eta_s4s9 = np.linspace(.7,.9,5)
    eta_iso = np.linspace(.1,.3,5)
    
    #list the cuts together
    pizero_cuts = (pi_cluster_pt, pi_pizero_pt, pi_elayer, pi_s4s9, pi_iso, ncri1, ncri2)
    eta_cuts = (eta_cluster_pt, eta_eta_pt, eta_elayer, eta_s4s9, eta_iso, ncri1, ncri2)

    #take all possible combinations
    pi_grid = list(itertools.product(*pizero_cuts))
    eta_grid = list(itertools.product(*eta_cuts))

    return (pi_grid, eta_grid)

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="hlt.root file name to analyze FILE",
                  action="store",type="string")

parser.add_option("-o", "--outfile", dest="outfilename",
                  help="tree.root file name to output",default="datatree.root",
                  action="store",type="string")

parser.add_option("-g", "--grid_point", dest="GRID_POINT",
                  help="grid point consisting of the cuts to use",
                  action="store",type="int",default=-1)

(options, args) = parser.parse_args()

file = rt.TFile(options.filename)
tree = file.Get("Tree_HLT")
outfile = rt.TFile(options.outfilename,"RECREATE")

grid = build_grids()[0]
cut = grid[options.GRID_POINT]




#id_cut = "STr2_ptG1_rec > %f && STr2_ptG2_rec > %f && STr2_ptPi0_rec > %f && STr2_S4S9_1 > %f && STr2_S4S9_2 > %f && STr2_IsoPi0_rec < %f && STr2_n1CrisPi0_rec > %i && STr2_n2CrisPi0_rec > %i" % (cut[0], cut[0], cut[1], cut[3], cut[3], cut[4], cut[5], cut[6])        

baseline_cut = "STr2_mPi0_rec > .05 && STr2_mPi0_rec < .3"
ncri1_cut = "STr2_n1CrisPi0_rec > %i" % cut[5]
ncri2_cut = "STr2_n2CrisPi0_rec > %i" % cut[6]
ptclus_cut = "STr2_ptG1_rec > %f && STr2_ptG2_rec > %f && STr2_ptPi0_rec > %f" % ( cut[0],cut[0],cut[1])
layer_cut = "((STr2_Es_e1_1 + STr2_Es_e2_1) > %f  && (STr2_Es_e1_2 + STr2_Es_e2_2) > %f)" % (cut[2],cut[2])    
s4s9_cut = "STr2_S4S9_1 > %f && STr2_S4S9_2 > %f" % (cut[3], cut[3])

list_cuts = [ncri1_cut, ncri2_cut, ptclus_cut, layer_cut, s4s9_cut]

nocut_hist = rt.TH1F("h0","h0",125,.05,.3)

nocut_hist.GetXaxis().SetTitle("#pi_{0} Invariant Mass [GeV]")
nocut_hist.GetYaxis().SetTitle("Events / ( 0.002 GeV )")

ncri1_hist = rt.TH1F("h1","h1",125,.05,.3)
ncri2_hist = rt.TH1F("h2","h2",125,.05,.3)
ptclus_hist = rt.TH1F("h3","h3",125,.05,.3)
layer_hist = rt.TH1F("h4","h4",125,.05,.3)
s4_s9_hist = rt.TH1F("h5","h5",125,.05,.3)

hist_list = [nocut_hist,ncri1_hist,ncri2_hist,ptclus_hist,layer_hist,s4_s9_hist]
names = ["No Cuts", "Cri1 Cut", "Cri2 Cut", "P_{t} Cluster Cut", "ES Cut","S_{4}/S_{9} Cut"]
prog_cuts = [baseline_cut]

for cut_i in range(len(list_cuts)):
    prog_cuts.append(addcuts([list_cuts[cut_i], prog_cuts[cut_i]]))

tree.Draw("STr2_mPi0_rec>>h0",prog_cuts[0])    
tree.Draw("STr2_mPi0_rec>>h1",prog_cuts[1])
tree.Draw("STr2_mPi0_rec>>h2",prog_cuts[2])
tree.Draw("STr2_mPi0_rec>>h3",prog_cuts[3])
tree.Draw("STr2_mPi0_rec>>h4",prog_cuts[4])
tree.Draw("STr2_mPi0_rec>>h5",prog_cuts[5])

#print prog_cuts
colors = [rt.kGray,rt.kOrange, rt.kGreen,rt.kViolet,rt.kPink,rt.kBlue]
legend = rt.TLegend(.415,.394,.86,.738)
legend.SetFillColor(0)

for ii in range(len(hist_list)): 
    hist_list[ii].SetFillColor(colors[ii])
    hist_list[ii].SetLineColor(rt.kBlack)
    if ii > 0:
        hist_list[ii].Draw("same")
    else:
        hist_list[ii].Draw()
        
    legend.AddEntry(hist_list[ii],names[ii],"f")


print "cut \t nevents"
print "nocut \t", int(nocut_hist.GetEntries()), "\t %.4f" % 1.0
print "ncri1\t", int(ncri1_hist.GetEntries()), "\t %.4f" % (float(ncri1_hist.GetEntries()) / nocut_hist.GetEntries())
print "ncri2\t", int(ncri2_hist.GetEntries()), "\t %.4f" % (float(ncri2_hist.GetEntries()) / nocut_hist.GetEntries())
print "ptclus\t", int(ptclus_hist.GetEntries()), "\t %.4f"% (float(ptclus_hist.GetEntries()) / nocut_hist.GetEntries())
print "elayer\t", int(layer_hist.GetEntries()), "\t %.4f" % (float(layer_hist.GetEntries()) / nocut_hist.GetEntries())
print "s4s9\t", int(s4_s9_hist.GetEntries()), "\t %.4f" % (float(s4_s9_hist.GetEntries()) / nocut_hist.GetEntries())

legend.Draw("same")


print prog_cuts[-1]

raw_input("")
