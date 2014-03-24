from  optparse  import OptionParser
import array
import ROOT as rt
from ROOT import AddressOf

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="hlt.root file name to analyze FILE",
                  action="store",type="string")

parser.add_option("-o", "--outfile", dest="outfilename",
                  help="tree.root file name to output",default="datatree.root",
                  action="store",type="string")

parser.add_option("--eta_b", dest="ETA_BEGIN",
                  help="minimum of eta range. This is only to be used for the dataset conversion i.e. the -f flag",
                  action="store",type="float",default=0)

parser.add_option("--eta_e", dest="ETA_END",
                  help="maximum of eta range.  This is only to be used for the dataset conversion i.e. the -f flag",
                  action="store",type="float",default=1000)



(options, args) = parser.parse_args()

def parse_array(vector,npizero):
    temp = []
    for ii in range(npizero): 
        temp.append(vector[ii])

    return temp

def parse_idx(eta_b,eta_e, eta_val, m_val):
    
    good_idx = []

    for ii in range(len(eta_val)):
        if abs(eta_val[ii]) > eta_b and abs(eta_val[ii]) < eta_e: #eta for the cateogry
            if m_val[ii] > .05 and m_val[ii] < .25: #we don't fit outside this region
                good_idx.append(ii)

    return good_idx
            
def parse_elem(idx, list):
    temp = []
    for ii in idx: temp.append(list[ii])

    return temp

file = rt.TFile(options.filename)
tree = file.Get("Tree_HLT")

out_file = rt.TFile(options.outfilename,"RECREATE")
out_tree = rt.TTree("Tree_HLT", "Trimmed Tree")

print "reading from file", file

  
# specify max events to keep in clouds 
nev = tree.GetEntries()

names = [ "STr2_NPi0_rec","STr2_Pi0recIsEB" ,"STr2_IsoPi0_rec", "STr2_n1CrisPi0_rec", "STr2_n2CrisPi0_rec", "STr2_mPi0_rec", "STr2_ptG1_rec", "STr2_ptG2_rec", "STr2_etaPi0_rec", "STr2_ptPi0_rec", "STr2_DeltaRG1G2", "STr2_Es_e1_1", "STr2_Es_e1_2","STr2_Es_e2_1", "STr2_Es_e2_2", "STr2_S4S9_1", "STr2_S4S9_2","L1bits"]
types = ["Int_t", "Int_t", "Float_t","Int_t","Int_t","Float_t","Float_t","Float_t","Float_t","Float_t","Float_t","Float_t","Float_t","Float_t","Float_t","Float_t","Float_t","Int_t"]
postpend = ["","[5000]","[5000]","[5000]","[5000]","[5000]","[5000]","[5000]","[5000]","[5000]","[5000]","[5000]","[5000]","[5000]","[5000]","[5000]","[5000]","[128]"]

processline = " typedef struct { "

for ii in range(len(names)):
    processline+= types[ii] + " " + names[ii] + postpend[ii]+";"
processline+= "} pizero_tree_t;"

#process the struct to build the tree object
rt.gROOT.ProcessLine(processline)

pizeros = rt.pizero_tree_t()

#build a new copy of the tree with branches in the struct
out_tree.Branch("STr2_NPi0_rec",AddressOf(pizeros,"STr2_NPi0_rec"),"STr2_NPi0_rec/I")
out_tree.Branch("STr2_Pi0recIsEB",pizeros.STr2_Pi0recIsEB,"STr2_Pi0recIsEB[STr2_NPi0_rec]/I")
out_tree.Branch("STr2_IsoPi0_rec",pizeros.STr2_IsoPi0_rec,"STr2_IsoPi0_rec[STr2_NPi0_rec]/F")
out_tree.Branch("STr2_n1CrisPi0_rec",pizeros.STr2_n1CrisPi0_rec,"STr2_n1CrisPi0_rec[STr2_NPi0_rec]/I")
out_tree.Branch("STr2_n2CrisPi0_rec",pizeros.STr2_n2CrisPi0_rec,"STr2_n2CrisPi0_rec[STr2_NPi0_rec]/I")
out_tree.Branch("STr2_mPi0_rec",pizeros.STr2_mPi0_rec,"STr2_mPi0_rec[STr2_NPi0_rec]/F")
out_tree.Branch("STr2_ptG1_rec",pizeros.STr2_ptG1_rec,"STr2_ptG1_rec[STr2_NPi0_rec]/F")
out_tree.Branch("STr2_ptG2_rec",pizeros.STr2_ptG2_rec,"STr2_ptG2_rec[STr2_NPi0_rec]/F")
out_tree.Branch("STr2_etaPi0_rec",pizeros.STr2_etaPi0_rec,"STr2_etaPi0_rec[STr2_NPi0_rec]/F")
out_tree.Branch("STr2_ptPi0_rec",pizeros.STr2_ptPi0_rec,"STr2_ptPi0_rec[STr2_NPi0_rec]/F")
out_tree.Branch("STr2_DeltaRG1G2",pizeros.STr2_DeltaRG1G2,"STr2_DeltaRG1G2[STr2_NPi0_rec]/F")
out_tree.Branch("STr2_Es_e1_1",pizeros.STr2_Es_e1_1,"STr2_Es_e1_1[STr2_NPi0_rec]/F")
out_tree.Branch("STr2_Es_e1_2",pizeros.STr2_Es_e1_2,"STr2_Es_e1_2[STr2_NPi0_rec]/F")
out_tree.Branch("STr2_Es_e2_1",pizeros.STr2_Es_e2_1,"STr2_Es_e2_1[STr2_NPi0_rec]/F")
out_tree.Branch("STr2_Es_e2_2",pizeros.STr2_Es_e2_2,"STr2_Es_e2_2[STr2_NPi0_rec]/F")
out_tree.Branch("STr2_S4S9_1",pizeros.STr2_S4S9_1,"STr2_S4S9_1[STr2_NPi0_rec]/F")
out_tree.Branch("STr2_S4S9_2",pizeros.STr2_S4S9_2,"STr2_S4S9_2[STr2_NPi0_rec]/F")
out_tree.Branch("L1bits",pizeros.L1bits,"L1bits[128]/I")

nevents = tree.Draw(">>iterlist" , "", "entrylist")

itlist = rt.gDirectory.Get("iterlist")

iev = 0

#loop over the entries in the tree
while iev < tree.GetEntries():
    if iev % 1000 == 0: print "Filling Tree...",iev
    iev += 1

    entry = itlist.Next()
    tree.GetEntry(entry)

    #produce python arrays from arrays contained inside the trees
    npiz = tree.STr2_NPi0_rec

    is_eb = parse_array(tree.STr2_Pi0recIsEB, npiz)
    iso = parse_array(tree.STr2_IsoPi0_rec, npiz)
    ncri1 = parse_array(tree.STr2_n1CrisPi0_rec, npiz)
    ncri2 = parse_array(tree.STr2_n2CrisPi0_rec, npiz)
    mpiz = parse_array(tree.STr2_mPi0_rec, npiz)
    ptg1 = parse_array( tree.STr2_ptG1_rec, npiz)
    ptg2 = parse_array( tree.STr2_ptG2_rec, npiz)
    eta =  parse_array( tree.STr2_etaPi0_rec, npiz)
    ptpi0 = parse_array( tree.STr2_ptPi0_rec, npiz)
    drg1g2 = parse_array( tree.STr2_DeltaRG1G2, npiz)
    ese1_1 = parse_array( tree.STr2_Es_e1_1, npiz)
    ese2_1 = parse_array( tree.STr2_Es_e2_1, npiz)
    ese1_2 = parse_array( tree.STr2_Es_e1_2, npiz)
    ese2_2 = parse_array( tree.STr2_Es_e2_2, npiz)
    s4s9_1 = parse_array(tree.STr2_S4S9_1, npiz)
    s4s9_2 = parse_array(tree.STr2_S4S9_2, npiz)
    l1bits = parse_array(tree.L1bits,128) 
    
    #make alist of the arrays to iterate over
    arrays = [is_eb,iso,ncri1, ncri2, mpiz, ptg1, ptg2, eta, ptpi0, drg1g2, ese1_1, ese2_1, ese1_2, ese2_2, s4s9_1, s4s9_2]

    #parse the indices of pizeros with eta in the correct regions
    idx = parse_idx(options.ETA_BEGIN, options.ETA_END, eta, mpiz)

    #keep only the elements passing eta cuts
    for ii in range(len(arrays)): arrays[ii] = parse_elem(idx, arrays[ii])

    passed_piz = len(arrays[0])
    if passed_piz > 0:
        pizeros.STr2_NPi0_rec = passed_piz
        pizeros.STr2_Pi0recIsEB = array.array("i",arrays[0])
        pizeros.STr2_IsoPi0_rec = array.array("f", arrays[1])
        pizeros.STr2_n1CrisPi0_rec = array.array("i",arrays[2])
        pizeros.STr2_n2CrisPi0_rec = array.array("i",arrays[3])
        pizeros.STr2_mPi0_rec = array.array("f",arrays[4])
        pizeros.STr2_ptG1_rec = array.array("f",arrays[5])
        pizeros.STr2_ptG2_rec =array.array("f",arrays[6])
        pizeros.STr2_etaPi0_rec =array.array("f",arrays[7])
        pizeros.STr2_ptPi0_rec = array.array("f",arrays[8])
        pizeros.STr2_DeltaRG1G2 =array.array("f",arrays[9])
        pizeros.STr2_Es_e1_1 = array.array("f",arrays[10])
        pizeros.STr2_Es_e2_1 = array.array("f",arrays[11])
        pizeros.STr2_Es_e1_2 = array.array("f",arrays[12])
        pizeros.STr2_Es_e2_2 = array.array("f",arrays[13])
        pizeros.STr2_S4S9_1 = array.array("f",arrays[14])
        pizeros.STr2_S4S9_2 = array.array("f",arrays[15])
        pizeros.L1bits = array.array("i",l1bits) #fill all of the l1 bits every time
        out_tree.Fill()
 
out_tree.Write()
out_file.Close()
