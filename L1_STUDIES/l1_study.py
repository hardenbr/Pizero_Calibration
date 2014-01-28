## Author: Joshua Hardenbrook - Princeton University - joshuarh@princeton.edu

##PARAMETERS TO CHANGE###
#(1) in the method build_grid() you must specify the values of the cuts to scan
#(2) in the method apply_tree_cut() you must specify how the grid cuts will be applied
#(3) in the method fit_dataset() you should put a reasonable range on the signal/bkg norm

from ROOT import gSystem, RooArgSet, RooFit
gSystem.Load('libRooFit')
from  optparse  import OptionParser
import numpy as np
import itertools, math
import ROOT as rt
#canvas style file
import rootlogon
rootlogon.style()
#dont draw anything
rt.gROOT.SetBatch(True)
#dont print  any roofit output besides errors/warnings
RooFit.PrintLevel(-1)
RooFit.Verbose(False)
rt.RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)

N_TRIGGERS = 128

parser = OptionParser()

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

parser.add_option("-w","--write", dest="do_write", help="write the fits and canvases",
                                   action="store_true", default=False)

(options, args) = parser.parse_args()

parser.print_help()

def get_hlt_hists(tree):
    print "Producting histograms for each bit"

    #determine the events in the tree and build an iteration list

    tree.Draw(">>iterlist", "", "entrylist")

    itlist = rt.gDirectory.Get("iterlist")
    
    #build the workspace
    workspace = rt.RooWorkspace("workspace")
    
    iev = 0

    l1_trig_idx = []
    l1_uniq_hists = []
    l1_raw_hists = []

    for ii in range(N_TRIGGERS):
        u_hist = rt.TH1F("u_hist_%i" % ii,"u_hist_%i" % ii, 100, .05, .25) 
        r_hist = rt.TH1F("r_hist_%i" % ii,"r_hist_%i" % ii, 100, .05, .25) 

        l1_uniq_hists.append(u_hist)
        l1_raw_hists.append(r_hist)

    #loop over the tree and add it to the RooDataSet
    while iev < tree.GetEntries():
        if iev % 5000 == 0: print "Filling L1 at event...",iev
        iev += 1
        entry = itlist.Next()
        
        #set the RooArgSet and save                                                          
        a = rt.RooArgSet(args)

        l1bits = tree.L1bits        
        sum_bits = sum(l1bits)

        #if it is uniq
        if sum_bits == 1:
            #find the uniq bit
            for il1 in range(len(l1bits)):                 
                if l1bits[il1] == 1: 
                    for ipi in range(tree.STr2_NPi0_rec): 
                        l1_uniq_hists[il1].Fill(tree.STr2_mPi0_rec)
                        l1_raw_hists[il1].Fill(tree.STr2_mPi0_rec)
                    break #we are done once we find the first one
        #if it is raw
        else:
            found_bits = 0
            for il1 in range(len(l1bits)):                 
                if l1bits[il1] == 1: 
                    sum_bits += 1
                    for ipi in range(tree.STr2_NPi0_rec): 
                        l1_raw_hists[il1].Fill(tree.STr2_mPi0_rec)
                if found_bits = sum_bits: break #we found them all

    return (l1_raw_hists, l1_uniq_hists)

def build_workspace_hist(tree,cut,l1num):
    #determine the events in the tree and build an iteration list
    hist_total = rt.TH1F("datahist_total","datahist_total",100,.05,.25)
    hist = rt.TH1F("datahist","datahist",100,.05,.25)

    cut_string = ""
    tree.Draw("STr2_mPi0_rec>>datahist",cut_string)
    tree.Draw("STr2_mPi0_rec>>datahist_total")
    
    #build the workspace
    workspace = rt.RooWorkspace("workspace")

    #declare our variables with ranges    
    variables = ["mpizero[.1., .05., .25]"]

    #factory all the variables
    for v in variables: workspace.factory(v)    

    #make the RooDataset
    args = workspace.allVars()
    data = rt.RooDataSet('PiTree','Tree of Pi0/Eta Information',args)

    #calculate the efficiency
    eff = float(hist.GetEntries())/ float(hist_total.GetEntries())

    return (workspace, hist, eff)

def fit_dataset(rdata,il1,eff):

    x  = rt.RooRealVar("mpizero","#pi_{0} invariant mass", .05, .25,"GeV")
    mean  = rt.RooRealVar("m","#pi_{0} peak position", .13, .11, .135,"GeV")
    sigma  = rt.RooRealVar("sigma","#pi_{0} core #sigma", .013, .011, .02,"GeV")
    gaus = rt.RooGaussian("gaus","Core Gaussian", x, mean, sigma)

    #t1 = rdata.reduce("mpizero < .25 && mpizero > .05")
    t1 = rt.RooDataHist("dh","#gamma#gamma invariant mass", rt.RooArgList(x), rdata)
    
    #build the background from a polynomial
    c0 = rt.RooRealVar("c0","c0",.2,-1,1)
    c1 = rt.RooRealVar("c1","c1",-.1,-1,1)
    c2 = rt.RooRealVar("c2","c2",.1,-1,1)
    c3 = rt.RooRealVar("c3","c3",-.1,-1,1)
    c4 = rt.RooRealVar("c4","c4",.1,-1,1)
    c5 = rt.RooRealVar("c5","c5",.1,-1,1)
    c6 = rt.RooRealVar("c6","c6",.3,-1,1)

    #using a polynomial background
    bkg_pars = rt.RooArgList(c0,c1,c2)
    bkg = rt.RooChebychev("bkg","bkg", x, bkg_pars)
    
    #add the signal and the background in a model
    tot = rdata.GetEntries()

    n_sig = rt.RooRealVar("nsig","#pi^{0} yield", tot*.2,tot*.05, tot*.3)
    n_bkg = rt.RooRealVar("nbkg","background yield",tot*.9, tot*.2, tot*.9)
    model =  rt.RooAddPdf("model","sig+bkg",rt.RooArgList(gaus,bkg), rt.RooArgList(n_sig,n_bkg))
    
    
    #    t1.Print()
    nll = rt.RooNLLVar("nll","log likelihood var", model, t1, True)
    m = rt.RooMinuit(nll)    
    m.setPrintLevel(-1)

    m.migrad()
    
    result = m.save()
    #    result.Print()
    
    #declare the signal over background range
    x.setRange("sobRange",mean.getVal()-2.*sigma.getVal(), mean.getVal()+2.*sigma.getVal())
    x.setRange("FULL",.05,.25)

    #calculate integrals
    integralBkg_sob = bkg.createIntegral(rt.RooArgSet(x),RooFit.NormSet(rt.RooArgSet(x)),RooFit.Range("sobRange"))
    integralBkg_full = bkg.createIntegral(rt.RooArgSet(x),RooFit.NormSet(rt.RooArgSet(x)),RooFit.Range("FULL"))
    
    normBkg_sob = integralBkg_sob.getVal()
    normBkg_full = integralBkg_full.getVal()

    bkg_scale = (normBkg_sob / normBkg_full)

    #put the frame on a canvas
    canvas = rt.TCanvas("canvas_%i" % il1)

    #make a frame
    frame = x.frame()

    #plot points
    t1.plotOn(frame)
    #plot fit components
    model.plotOn(frame, RooFit.Components("bkg"),RooFit.LineStyle(rt.kDashed),RooFit.LineColor(rt.kBlue))
    model.plotOn(frame, RooFit.Components("gaus"),RooFit.LineColor(rt.kRed))
    #plot the full fit
    model.plotOn(frame)

    #draw the frame
    frame.Draw()

    ####compute the interesting values

    #sob calculations
    n_s_sob = n_sig.getVal() * .9546
    n_b_sob = n_bkg.getVal() * bkg_scale
    s_over_b = n_s_sob / n_b_sob 

    #goodness of fit
    error_e = s_over_b * math.sqrt(math.pow(n_sig.getError() / n_sig.getVal(), 2) + math.pow(n_bkg.getError() / n_bkg.getVal(), 2))
    chi2 = frame.chiSquare()
    
    #signal calculations
    mean_val = mean.getVal()
    sigma_val = sigma.getVal()
    mu_over_err = mean.getVal() / mean.getError()
    
    if options.verbose:
        result.Print()    
        print "normBkg_sob", normBkg_sob
        print "normBkg_full", normBkg_full
        print "signal",n_s_sob,"bkg", n_b_sob
        print "sob:", s_over_b
        print "chi^2:", chi2
        print "error_e:", error_e

    #write the latex onto the canvas
    lat = rt.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(.04)
    lat.SetTextColor(1)
    
    ymin = .91
    xmin = .52
    ypass = .06
    
    line = "S/B: %.4f #pm %.5f ( %1.1f / %1.1f)" % (s_over_b, error_e, n_s_sob, n_b_sob)
    line2 = "reduced #chi^{2}: %.4f" % chi2
    line3 = "#mu: %.4f, #sigma: %.4f, #mu/err: %.1f" % (mean_val,sigma_val,mu_over_err)
    line4 = "Efficiency %.6f" % eff
    line5 = "Grid #: %i" % il1
    
    lat.DrawLatex(xmin,ymin,line)
    lat.DrawLatex(xmin,ymin-ypass,line2)
    lat.DrawLatex(xmin,ymin-2*ypass,line3)
    lat.DrawLatex(xmin,ymin-3*ypass,line4)
    lat.DrawLatex(xmin,ymin-4*ypass,line5)
 
    return (result, canvas, n_s_sob, n_b_sob, s_over_b, chi2, error_e, mean_val, sigma_val, mu_over_err, eff)

#################
## MAIN METHOD ##
##             ## 
#################

tree_set = []
sum_trees = None
output = None

#build the data workspace if there is no rooDatasets Specified

dataset_file = open(options.dataset) 
dataset_file_lines = dataset_file.readlines()
dataset_file_lines_stripped = map(lambda(x):x.rstrip("\n"),dataset_file_lines)    

#add all the trees to the tree_set
temp_file = None
#build the list of files
file_set = map(lambda(x):rt.TFile(x),dataset_file_lines_stripped)
#build the list of trees in the files
tree_set = map(lambda(x):x.Get("Tree_HLT"), file_set)
rdata = None

#output file containing fits and canvases
output = rt.TFile(options.outfilename,"RECREATE")

#output files containing cut values and fit calculations
outfile_dir = options.outfilename[:-5]

fit_params_string = "@GRID#\t\tNSIG\tNBKG\tSOB\tCHI^2\tERR_E\tMU\tSIGMA\tMU/ERR\tEFF\n"
cut_string = "@grid#\t\tga_pt\tpi_pt\telyr\ts4s9\tiso\tncri1\tncri2\n"

list = rt.TList() 
for tree in tree_set: list.Add(tree)
sum_trees = rt.TTree.MergeTrees(list)
sum_trees.SetName("Tree_HLT")


(raw_hists, uniq_hists) = get_hlt_hists(sum_trees))

#scan point by point
for il1 in range(N_TRIGGERS)
    print "Scanning L1 Bit:", il1

    uniq_data = uniq_hists[il1]
    raw_data = raw_hists[il1]

    if raw_data.GetEntries() == 0:
        fit_result = None
        print "...not fitting.."
    else:
        print "...Fitting..."

        #generate the fit result
        fit_result = fit_dataset(rdata, il1, eff)

        if options.do_write:
            #sum_trees.Write("ttree_%i" % iev)
            rdata.Write("tree_%i" % il1)
            #write the result
            fit_result[0].Write("fit_%i" % il1) 
            #write the frame
            fit_result[1].Write("frame_%i" % il1)

        #make sure tabbing is correct for output
        fit_params_string += "@@%i \t" % il1
        if iev <= 9999: fit_params_string += "\t" 

        #format the result if there was a fit result
        if fit_result != None:
            for jj in fit_result[2:]:
                if jj > 1000: fit_params_string+="%.1g\t" % jj
                else:  fit_params_string+="%.2g\t" % jj
        else:
            fit_params_string+="NO_RESULT EFFICIENCY TOO LOW eff=%f" % eff
        fit_params_string+="\n"

    print "\n"
    print fit_params_string
    print "\n"

if options.do_write:
    workspace.Write()
    output.Close()
    

