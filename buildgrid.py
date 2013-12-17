## Author: Joshua Hardenbrook - Princeton University - joshuarh@princeton.edu

##PARAMETERS TO CHANGE###
#(1) in the method build_grid() you must specify the values of the cuts to scan
#(2) in the method apply_cut() you must specify how the grid cuts will be applied
#(3) in the method fit_cut_dataset() you should put a reasonable range on the signal/bkg norm

from ROOT import gSystem, RooArgSet, RooFit
gSystem.Load('libRooFit')
from  optparse  import OptionParser
import numpy as np
import itertools, math
import ROOT as rt
import rootlogon
rootlogon.style()
rt.gROOT.SetBatch(True)
RooFit.PrintLevel(-1)
RooFit.Verbose(False)
rt.RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)


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

parser.add_option("-b", "--begin", dest="GRID_BEGIN",
                  help="starting grid point (must specify end)",
                  action="store",type="int",default=-1)

parser.add_option("-l", "--gridlist", dest="GRID_LIST",
                  help="txt file listing grid points to scan. Separated by linebreaks",
                  action="store",type="string",default="no_list")

parser.add_option("-e", "--end", dest="GRID_END",
                  help="ending grid point (-1 means to the end of the grid)",
                  action="store",type="int",default=-1)

parser.add_option("-d", "--dataset", dest="dataset",
                  help="list of root files containing the processesed datasets. If this is not specified, the roodataset will be generated and with the name roo(NAME).root",
                  action="store",type="string",default="no_file")

parser.add_option("-w","--write", dest="do_write", help="write the fits and canvases",
                                   action="store_true", default=False)

(options, args) = parser.parse_args()

#parser.print_help()

def apply_cut(data, cut):
    id_cut = "pt_g1 > %f && pt_g2 > %f && pi_pt > %f && pi_s4s9_1 > %f && pi_s4s9_2 > %f && pi_iso < %f && pi_ncri_1 > %i && pi_ncri_2 > %i" % (cut[0], cut[0], cut[1], cut[3], cut[3], cut[4], cut[5], cut[6])    
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

def fit_cut_dataset(dataset,cut,iev):
    #apply the cut
    rdata = apply_cut(dataset,cut)    

    #efficiency calculation
    eff_window_cut = "mpizero > .05 && mpizero < .3"
    
    total_events = (data.reduce(eff_window_cut)).numEntries()
    post_selection = (rdata.reduce(eff_window_cut)).numEntries()
    eff = float(post_selection) / total_events

    t1 = rdata.reduce("mpizero < .25 && mpizero > .05")
    x  = rt.RooRealVar("mpizero","#pi_{0} invariant mass", .05, .25,"GeV")
    mean  = rt.RooRealVar("m","#pi_{0} peak position", .13, .12, .14,"GeV")
    sigma  = rt.RooRealVar("sigma","#pi_{0} core #sigma", .010, .009, .016,"GeV")
    gaus = rt.RooGaussian("gaus","Core Gaussian", x, mean, sigma)

    
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
    n_sig = rt.RooRealVar("nsig","#pi^{0} yield",1001,1000,3000)
    n_bkg = rt.RooRealVar("nbkg","background yield",2000,100, 1e5)
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
    canvas = rt.TCanvas("canvas_%i" % iev)

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
    line5 = "Grid #: %i" % iev
    
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
    for ii in roodatasets[1:]: data.append(ii)

#build the cut grids
(pi_grid, eta_grid) = build_grids()

print "Number of grid points for pizero:",  len(pi_grid)

rdata = None
#output file containing fits and canvases
if options.do_write: output = rt.TFile(options.outfilename,"RECREATE")

#output files containing cut values and fit calculations
outfile_dir = options.outfilename[:-5]

fit_params_string = "@GRID#\t\tNSIG\tNBKG\tSOB\tCHI^2\tERR_E\tMU\tSIGMA\tMU/ERR\tEFF\n"
cut_string = "@grid#\t\tga_pt\tpi_pt\telyr\ts4s9\tiso\tncri1\tncri2\n"

#loop over each set of cuts and reduce the dataset + fit
p1 = options.GRID_BEGIN
p2 = options.GRID_END
grid_list = options.GRID_LIST
scan_points = None
iev_points = None

#parse out the piece of the grid to scan
#either take the begining to end or read from the input file
if p1 != -1 and p2 != -1 and grid_list == "no_list":
    iev_points = range(p1,(p2+1))
elif p1 != -1 and p2 == -1 and grid_list == "no_list":
    iev_points = range(p1,len(pi_grid)+1)
#if beginning end end are not specified and grid list is
elif grid_list != "no_list":
    #parse out the list 
    list_file = open(grid_list)
    list_file_lines = list_file.readlines()
    iev_points = map(lambda(x):int(x.rstrip("\n")),list_file_lines)        
else:
    iev_points = range(len(pi_grid)+1)
    print "------WARNING!!!----- SCANNING ALL GRID POINTS"


#scan point by point
for iev in iev_points:    
    print "Scanning grid point", iev, "..."
    
    if options.dataset == "no_file":
        rdata = data 
        break
    else:
        #write out the data after the cut is applied and the fit result
        rdata = apply_cut(data,pi_grid[iev])
        fit_result = fit_cut_dataset(data, pi_grid[iev], iev)
        
        if options.do_write:
            rdata.Write("tree_%i" % iev)
            #write the result
            fit_result[0].Write("fit_%i" % iev) 
            #write the frame
            fit_result[1].Write("frame_%i" % iev)

        #write out the values of the cuts    
        cut_string += "@@@%i \t" % iev
        for cuts in pi_grid[iev]: cut_string+="%2.3f \t" % cuts
        cut_string += "\n"        

        fit_params_string += "@@%i \t" % iev
        if iev <= 9999: fit_params_string += "\t" 

        for jj in fit_result[2:]:
            if jj > 10: fit_params_string+="%6.1f\t" % jj
            else:  fit_params_string+="%2.5f\t" % jj
        fit_params_string+="\n"

    print "\n"
    print cut_string
    print fit_params_string
    print "\n"

if options.dataset == "no_file": rdata.Write("tree_00")

if options.do_write:
    workspace.Write()
    output.Close()


