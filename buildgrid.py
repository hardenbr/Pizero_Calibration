## Author: Joshua Hardenbrook - Princeton University - joshuarh@princeton.edu
from ROOT import gSystem, RooArgSet, RooFit
gSystem.Load('libRooFit')
from  optparse  import OptionParser
import numpy as np
import itertools, math
import ROOT as rt

rt.gROOT.SetBatch(True)

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="hlt.root file name to analyze FILE",
                  action="store",type="string")

parser.add_option("-o", "--outfile", dest="outfilename",
                  help="tree.root file name to output",default="datatree.root",
                  action="store",type="string")

parser.add_option("-v", "--verbose", dest="verbose",
                  help="print more information",
                  action="store_true")

parser.add_option("-n", "--nevents", dest="nevents",
                  help="only analyze the first n events",
                  action="store",type="int",default=-1)

parser.add_option("-d", "--dataset", dest="dataset",
                  help="list of root files containing the processesed datasets. If this is not specified, the roodataset will be generated and with the name roo(NAME).root",
                  action="store",type="string",default="no_file")


(options, args) = parser.parse_args()

parser.print_help()

def apply_cut(data,cut):
    id_cut = "mpizero > %f && pt_g1 > %f && pt_g2 > %f && pi_pt > %f && pi_s4s9_1 > %f && pi_s4s9_2 > %f && pi_iso < %f && pi_ncri_1 > %i && pi_ncri_2 > %i" % (0, cut[0], cut[0], cut[1], cut[3], cut[3], cut[4], cut[5], cut[6])    

#    id_cut = "mpizero < %f && pt_g1 > %f && pt_g2 > %f && pi_pt > %f && pi_s4s9_1 > %f && pi_s4s9_2 > %f && pi_iso < %f && pi_ncri_1 > %i && pi_ncri_2 > %i" % (.25,1 , 1, 1.5, .7, .7, .3, 7, 5)    

    es_cut = "pi_iseb || ((es_e1_1 + es_e1_2) > %f  && (es_e2_1 + es_e2_2) > %f)" % (cut[2],cut[2])

    if options.verbose: print id_cut,"\n","ES cut:", es_cut,"\n"
    return (data.reduce(id_cut)).reduce(es_cut)

def set_values(set,tree,ii):
    set.setRealValue("mpizero",tree.STr2_mPi0_rec[ii])
    set.setRealValue("pi_iseb",tree.STr2_Pi0recIsEB[ii]) 
    set.setRealValue("pi_iso,",tree.STr2_IsoPi0_rec[ii])
    set.setRealValue("pi_s4s9_1",tree.STr2_S4S9_1[ii])
    set.setRealValue("pi_s4s9_2",tree.STr2_S4S9_2[ii])
    set.setRealValue("pi_ncri_1",tree.STr2_n1CrisPi0_rec[ii])
    set.setRealValue("pi_ncri_2",tree.STr2_n2CrisPi0_rec[ii])
    set.setRealValue("pt_g1",tree.STr2_ptG1_rec[ii])
    set.setRealValue("pt_g2",tree.STr2_ptG2_rec[ii])
    set.setRealValue("pi_pt",tree.STr2_ptPi0_rec[ii])
    set.setRealValue("es_e1_1",tree.STr2_Es_e1_1[ii])
    set.setRealValue("es_e1_2",tree.STr2_Es_e1_2[ii])
    set.setRealValue("es_e2_1",tree.STr2_Es_e2_1[ii])
    set.setRealValue("es_e2_2",tree.STr2_Es_e2_2[ii])
    return set

def build_workspace(input_file):
    tree = input_file.Get("Tree_HLT")

    #build a list of the tree
    tree.Draw(">>iterlist","STr2_mPi0_rec >0","entrylist")
    itlist = rt.gDirectory.Get("iterlist")
    
    #build the workspace
    workspace = rt.RooWorkspace("workspace")

    #declare our variables with ranges
    variables = ["npizero[1,0,1000]","mpizero[.1., .05., .25]","pi_iseb[0,0,1]","pi_iso[0,-10,10]","pi_s4s9_1[0,0,10]","pi_s4s9_2[0,0,10]","pi_ncri_1[0,0,10]","pi_ncri_2[0,0,10]","pt_g1[0,0,20]","pt_g2[0,0,20]","es_e1_1[0,0,10]","es_e1_2[0,0,10]","es_e2_1[0,0,10]","es_e2_2[0,0,10]","pi_pt[0,0,20]"]

    #factory all the variables
    for v in variables:
        workspace.factory(v)    

    #make the RooDataset
    args = workspace.allVars()
    data = rt.RooDataSet('PiTree','Tree of Pi0/Eta Information',args)

    iev = 0
    #loop over the tree and add it to the RooDataSet
    while iev < tree.GetEntries():
        if iev % 1000 == 0: print "Filling Tree...",iev
        iev += 1
        entry = itlist.Next()

        #break early if requested
        if iev == options.nevents: break
        tree.GetEntry(entry)
        
        #set the RooArgSet and save                                                          
        a = rt.RooArgSet(args)

        # for each reconstructed pion, add its values to the dataset
        for ii in range(tree.STr2_NPi0_rec):
            data.add(set_values(a,tree,ii))    

    getattr(workspace,'import')(data)    

    return (workspace,data)

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

def fit_dataset(dataset):
    t1 = dataset.reduce("mpizero < .2 && mpizero > .08")
    x  = rt.RooRealVar("mpizero","#pi_{0} invariant mass", .08, .2,"GeV")
    mean  = rt.RooRealVar("m","#pi_{0} peak position", .13, .12, .14,"GeV")
    sigma  = rt.RooRealVar("sigma","#pi_{0} core #sigma", .010, .005, .04,"GeV")
    gaus = rt.RooGaussian("gaus","Core Gaussian", x, mean, sigma)

    
    #build the background from a polynomial
    c0 = rt.RooRealVar("c0","c0",50,0,200)
    c1 = rt.RooRealVar("c1","c1",-1,-1000,1000)
    c2 = rt.RooRealVar("c2","c2",10,-1000,1000)
    c3 = rt.RooRealVar("c3","c3",1,-100,100)
    c4 = rt.RooRealVar("c4","c4",.1,-1,1)
    c5 = rt.RooRealVar("c5","c5",.1,-1,1)
    c6 = rt.RooRealVar("c6","c6",.3,-1,1)
    bkg_pars = rt.RooArgList(c0,c1,c2,c3)
    bkg = rt.RooPolynomial("bkg","bkg", x, bkg_pars)
    
#add the signal and the background
    n_sig = rt.RooRealVar("nsig","#pi^{0} yield",1000,500,1e5)
    n_bkg = rt.RooRealVar("nbkg","background yield",2000,1000, 1e5)
    
    model =  rt.RooAddPdf("model","sig+bkg",rt.RooArgList(gaus,bkg), rt.RooArgList(n_sig,n_bkg))
    
    
    t1.Print()
    nll = rt.RooNLLVar("nll","log likelihood var", model, t1, True)
    m = rt.RooMinuit(nll)
    m.migrad()
    
    result = m.save()
    
    #declare the sob range
    x.setRange("sobRange",mean.getVal()-2.*sigma.getVal(), mean.getVal()+2.*sigma.getVal())
    x.setRange("FULL",.08,.2)

    #calculate integrals
    integralBkg_sob = bkg.createIntegral(rt.RooArgSet(x),RooFit.NormSet(rt.RooArgSet(x)),RooFit.Range("sobRange"))
    integralBkg_full = bkg.createIntegral(rt.RooArgSet(x),RooFit.NormSet(rt.RooArgSet(x)),RooFit.Range("FULL"))
    
    normBkg_sob = integralBkg_sob.getVal()
    normBkg_full = integralBkg_full.getVal()

    bkg_scale = (normBkg_sob / normBkg_full)
    
    #make a frame
    frame = x.frame()

    #plot plints
    t1.plotOn(frame)
    #plot fit components
    model.plotOn(frame, RooFit.Components("bkg"),RooFit.LineStyle(rt.kDashed),RooFit.LineColor(rt.kBlue))
    model.plotOn(frame, RooFit.Components("gaus"),RooFit.LineColor(rt.kRed))
    #plot the full fit
    model.plotOn(frame)


    n_s_sob = n_sig.getVal() * .9546
    n_b_sob = n_bkg.getVal() * bkg_scale
    s_over_b = n_s_sob / n_b_sob 
    chi2 = frame.chiSquare()
    error_e = s_over_b * math.sqrt(math.pow(n_sig.getError() / n_sig.getVal(), 2) + math.pow(n_bkg.getError() / n_bkg.getVal(), 2))
    
    if options.verbose:
        result.Print()    
        print "normBkg_sob", normBkg_sob
        print "normBkg_full", normBkg_full
        print "signal",n_s_sob,"bkg", n_b_sob
        print "sob:", s_over_b
        print "chi^2:", chi2
        print "error_e:", error_e
 
    return (result, frame, n_s_sob, n_b_sob, s_over_b, chi2, error_e)

#################
## MAIN METHOD ##
##             ## 
#################

#build the data workspace if there is no rooDatasets Specified
if options.dataset == "no_file":
    #get the input file
    input_file = rt.TFile(options.filename)
    (workspace,data) = build_workspace(input_file)

#if the datset is specified added the roodatsets into one dataset
else:
    dataset_file = open(options.dataset) 
    dataset_file_lines = dataset_file.readlines()
    dataset_file_lines_stripped = map(lambda(x):x.rstrip("\n"),dataset_file_lines)    

    roodatasets = []
    workspace = None

    for ii in dataset_file_lines_stripped:
        workspace_file = rt.TFile(ii)
        data_temp = workspace_file.Get("tree_00")

        roodatasets.append(data_temp)
        workspace = workspace_file.Get("workspace")
    
    #add all the datasets together
    data = roodatasets[0]
    for ii in roodatasets[1:]:
        data.append(ii)

#build the cut grids
(pi_grid, eta_grid) = build_grids()

print "Number of grid points for pizero:",  len(pi_grid)

rdata = None
#output file containing fits and canvases
output = rt.TFile(options.outfilename,"RECREATE")

#output files containing cut values and fit calculations
fit_params = open("fit_result.txt","w")
fit_params.write("GRID# \t N SIGNAL \t N BKG \t \t SOB \t CHI^2 \t ERROR_E \n")
cut_values = open("cut_values.txt","w")
cut_values.write("grid# \t ga_pt \t pi_pt \t elyr \t s4s9 \t iso \t ncri1 \t ncri2 \n")

#loop over each set of cuts and reduce the dataset + fit
iev = 1
for ii in pi_grid[0:5]:    
    print "Scanning grid point", iev, "..."
    
    if options.dataset == "no_file":
        rdata = data 
        break
    else:
        #write out the data after the cut is applied and the fit result
        rdata = apply_cut(data,ii)    
        fit_result = fit_dataset(rdata)
        rdata.Write("tree_%i" % iev)
        #write the result
        fit_result[0].Write("fit_%i" % iev) 
        #write the frame
        fit_result[1].Write("frame_%i" % iev)

        #write out the values of the cuts    
        cut_values.write("%i \t" % iev)
        for cut in ii: cut_values.write("%2.3f \t" % cut)
        cut_values.write("\n")

        #write out the variables generated from the fit
        fit_params.write("%i \t" % iev)
        for jj in fit_result[2:]: fit_params.write("%10.2f \t" % jj)
        fit_params.write("\n")

    iev+=1

if options.dataset == "no_file": rdata.Write("tree_00")

workspace.Write()
output.Close()
fit_params.close()
cut_values.close()

