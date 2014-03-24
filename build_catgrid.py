from ROOT import gSystem, RooArgSet, RooFit, AddressOf, gROOT
gSystem.Load('libRooFit')
from  optparse  import OptionParser
import numpy as np
import itertools, math, pickle
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

class analysis:
    def __init__(self, grid, tree, grid_list):
        self.grid = grid

        self.grid_points = []  
        self.categories = []
        self.tree = tree
        self.grid_list = grid_list
        
    def add_category(self, eta_begin, eta_end):
        cat = category_def(eta_begin, eta_end)
        self.categories.append(cat)
        
    def initialize_grid(self):
        if self.categories == []:
            print " -- THERE ARE NO CATEGORIES. DEFINE BEFORE INTIALIZING GRID -- "
            exit(1)

        print "-- Initializing Grid --"
        for ii in self.grid:                        
            idx =  grid.index(ii)

            #only analyize grid points in the grid list
            if idx not in self.grid_list: continue
            
            #build the grid point
            point = grid_point(ii, idx, self.categories)

            self.grid_points.append(point)

    def build_grid_hists(self):

        print "-- Building Grid Histograms --"

        for point in self.grid_points:

            if point.grid_id % 100 == 0: print "-- Scanning Grid Point %i --" % point.grid_id
            point.get_category_hists(self.tree)

    def fit_grid_hists(self):

        for point in self.grid_points:
            point.fit_category_data()

        print "-- Fitting Complete --"
            
    def save_category_data(self, outfile):
        outfile.cd()

        for point in self.grid_points:
            for cat in point.cat_data:

                name_fit_res = "fit_result_%i_%i_%i" % (point.grid_id, cat.eta_b*100, cat.eta_e*100)
                name_hist = "hist_%i_%i_%i" % (point.grid_id, cat.eta_b*100, cat.eta_e*100)
                name_canvas = "canvas_%i_%i_%i" % (point.grid_id, cat.eta_b*100, cat.eta_e*100)

                #write the histgram for the category
                cat.rdata.Write(name_hist)
                
                #get the fit result object
                fit_result = cat.fit_result

                if fit_result == None: continue
                
                #write the canvas for tehc category
                fit_result.canvas.Write(name_canvas)
                fit_result.result.Write(name_fit_res)
                
    #store all the hard work in a histogram
    def pickle_hists(self, pickle_out):
        
        pickle.dump(self.grid_points, open(pickle_out, "wb"))
        

    #unpack multiple files of histograms and merge them                 
    def un_pickle_list(self, file_list):
        file = open(file_list,"r")
        lines = map(lambda x: x.rstrip("\n"), file.readlines())

        file_points = []

        print "-- Unpickling points from pickle list --"

        #merge all the points into a single list
        for ll in lines:
            pickled_points = pickle.load(open(ll,"rb"))
            for point in pickled_points:
                print "grid id:", point.grid_id
            
                file_points.append(point)

        #re-assign the points to the analysis object
        self.grid_points = file_points

    #build a plot of sob vs s for each category 
    def build_sob_vs_s(self, outfile):

        
        gROOT.ProcessLine("struct MyStruct{ Float_t s; Float_t b; Float_t sob; Int_t cat; Int_t gid;};")

        s = rt.MyStruct()

        outfile.cd()

        tree = rt.TTree("sob_tree", "sob info")

        tree.Branch("ns",AddressOf(s, "s"), "ns/F")
        tree.Branch("nb",AddressOf(s, "b"), "nb/F")
        tree.Branch("sob",AddressOf(s, "sob"), "sob/F")
        tree.Branch("gid",AddressOf(s, "gid"), "gid/I")
        tree.Branch("cat",AddressOf(s, "cat"), "cat/I")

        hist_list = []

        n_cats = len(self.categories)
        for cc in range(n_cats):
            hist_cc  = rt.TH2D("category_sob_vs_s_%i" % cc, "2d hist of grid performance", 1000, 0, 15000, 100, 0, 2.5)
                              
            for point in self.grid_points:
                id = point.grid_id
                cat = point.cat_data[cc]

                fit_res = cat.fit_result

                #assign values
                s.cat = cc
                s.gid = id
                s.sob = fit_res.sob
                s.ns = fit_res.ns
                s.nb = fit_res.nb
                
                tree.Fill()

        tree.Write()
        
#abstract class for category definition (eta ranges, maybe more later)
class category_def:

    def __init__(self, eta_b, eta_e):
        self.eta_b = eta_b
        self.eta_e = eta_e

#holder class for the data corresonding to a grid point
class category_data:

    def __init__(self, category_def, workspace, rdata, eff):
        self.eta_b = category_def.eta_b
        self.eta_e = category_def.eta_e
        self.cat = category_def

        self.workspace = workspace
        self.rdata = rdata
        self.eff = eff

        self.fit_result = None

#container for a fit result
class fit_result:
    def __init__(self, result, canvas, n_s, n_b, sob, chi2, error_e, mean_val, sigma, mu_over_error, eff):
        self.result, self.canvas, self.ns, self.nb, self.sob, self.chi2, self.error_e, self.mean_val, self.sigma, self.mu_over_error, self.eff =  result, canvas, n_s, n_b, sob, chi2, error_e, mean_val, sigma, mu_over_error, eff 
        
#fundamental object of the analysis, holds the categories, cut point n the grid
class grid_point:
    def __init__(self, tuple, grid_id, categories):
        self.cat_defs = categories
        self.tuple = tuple
        self.grid_id = grid_id
        self.cat_data = []

    def clear_fit(self):
        self.cat_data.fit_result = None
        self.cat_data.workspace = None            

    #for each category extract the histogram characterizing the cut
    def get_category_hists(self, tree):
        for cat in self.cat_defs:
            (workspace, rdata, eff) = build_workspace_hist(tree, self.tuple, self.grid_id, cat)

            cat_result = category_data(cat, workspace, rdata, eff)

            self.cat_data.append(cat_result)

    #fit each category for this grid point using the histogram stored previously
    def fit_category_data(self):
        for cat in self.cat_data:

            rdata = cat.rdata
            eff = cat.eff

            print "\n-- Fitting Grid Point %i Category: (%f, %f) --\n" % (self.grid_id, cat.eta_b, cat.eta_e)
            
            result = fit_dataset(rdata, self.grid_id, eff, 0, cat, self.tuple)

            cat.fit_result = result
                        
            
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

parser.add_option("-z", "--veto", dest="GRID_VETO",
                  help="txt file listing grid scan results. The points with results will not be re-fit",
                  action="store",type="string",default="no_list")

parser.add_option("-d", "--dataset", dest="dataset",
                  help="The Tree_HLT root file containing tree to be analyzed",
                  action="store",type="string",default="no_file")

parser.add_option("-p", "--pickle", dest="pickle",
                  help="name of .txt file listing pickled histogram data. If not specified, the histograms will be pickled and outputted",
                  action="store",type="string",default="no_file")

parser.add_option("-w","--write", dest="do_write", help="write the fits and canvases",
                                   action="store_true", default=False)

(options, args) = parser.parse_args()

def generate_tree_cut(cut, category):

    eta_cut = "(STr2_etaPi0_rec*STr2_etaPi0_rec) > %f && (STr2_etaPi0_rec*STr2_etaPi0_rec) < %f" % (category.eta_b*category.eta_b, category.eta_e * category.eta_e)

    id_cut = "STr2_ptG1_rec > %f && STr2_ptG2_rec > %f && STr2_ptPi0_rec > %f && STr2_S4S9_1 > %f && STr2_S4S9_2 > %f && STr2_IsoPi0_rec < %f && STr2_n1CrisPi0_rec > %i && STr2_n2CrisPi0_rec > %i" % (cut[0], cut[0], cut[1], cut[2], cut[2], cut[3], cut[4], cut[5])

    #    es_cut = "STr2_Pi0recIsEB || ((STr2_Es_e1_1 + STr2_Es_e2_2) > %f  && (STr2_Es_e2_1 + STr2_Es_e2_2) > %f)" % (cut[2],cut[2])    

    total_cut = "(" + id_cut + ") && (" + eta_cut + ")"

    print total_cut

    return total_cut


def build_workspace_hist(tree, cut, grid_id, category):

    print "-- Building workspace hist ", grid_id, "--"

    #determine the events in the tree and build an iteration list
    #hist_total = rt.TH1F("datahist_total","datahist_total",100,.05,.25)

    name = "datahist_%i_%i_%i" % (grid_id, category.eta_b*100, category.eta_e*100)
    hist = rt.TH1F(name, name, 100, .05, .25)

    cut_string = generate_tree_cut(cut, category)

    tree.Draw("STr2_mPi0_rec>>%s" % name, cut_string)
    #tree.Draw("STr2_mPi0_rec>>datahist_total")

    #print "TOTAL HIST", hist_total.Integral()
    print "SELECTED HIST", hist.Integral()
    
    #build the workspace
    workspace = rt.RooWorkspace("workspace")

    variables = ["mpizero[.1., .05., .25]"]

    #factory all the variables
    for v in variables: workspace.factory(v)    

    #make the RooDataset
    args = workspace.allVars()
    data = rt.RooDataSet('PiTree','Tree of Pi0/Eta Information',args)

    iev = 0
    nselected = 0
    ntotal = 0

    #calculate the efficiency
    eff = float(hist.GetEntries()) #/ float(hist_total.GetEntries())

    return (workspace, hist, eff)

def build_grid():
    #cystals scan range
    ncri1 = range(0, 4)
    ncri2 = range(0, 4)

    #pizero scan range
    pi_cluster_pt = np.linspace(.4, 1, 3)
    pi_pizero_pt = np.linspace(1.2, 2.4, 3)
    pi_s4s9 = np.linspace(.7, .9, 3)
    pi_iso = np.linspace(.1, .3, 3)
    
    #list the cuts together
    pizero_cuts = (pi_cluster_pt, pi_pizero_pt, pi_s4s9, pi_iso, ncri1, ncri2)

    #take all possible combinations
    pi_grid = list(itertools.product(*pizero_cuts))

    return pi_grid

def fit_dataset(rdata, il1, eff, iSamp, cat, cut):

    x  = rt.RooRealVar("mpizero","#pi^{0} invariant mass", .05, .25,"GeV")
    mean  = rt.RooRealVar("m","#pi^{0} peak position", .13, .10, .135,"GeV")
    sigma  = rt.RooRealVar("sigma","#pi^{0} core #sigma", .01, .0085, .028,"GeV")
    gaus = rt.RooGaussian("gaus","Core Gaussian", x, mean, sigma)

    #define the frame for x
    frame = x.frame()

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
    tot = rdata.Integral()#GetEntries()
    window = rdata.Integral(rdata.FindBin(.09),rdata.FindBin(.15))
    print "%i TOTAL IN WINDOW: %f" % (il1,window)

    n_sig = rt.RooRealVar("nsig","#pi^{0} yield", tot*.1,tot*.05, tot*.25)
    n_bkg = rt.RooRealVar("nbkg","background yield",tot*.7, tot*.5, tot*.95)
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

    #plot points
    t1.plotOn(frame)
    #plot fit components
    model.plotOn(frame, RooFit.Components("bkg"),RooFit.LineStyle(rt.kDashed),RooFit.LineColor(rt.kBlue))
    model.plotOn(frame, RooFit.Components("gaus"),RooFit.LineColor(rt.kRed))
    #plot the full fit
    model.plotOn(frame)


    #put the frame on a canvas
    can_name = None

    if iSamp == 0: can_name = "grid_canvas%i_cat_%2.2f_%2.2f" % (il1, cat.eta_b, cat.eta_e)

    canvas = rt.TCanvas(can_name,can_name,1000,900)
    ####compute the interesting values

    #sob calculations
    n_s_sob = n_sig.getVal() * .9546
    n_b_sob = n_bkg.getVal() * bkg_scale
    s_over_b = 0
    if n_b_sob != 0:
        s_over_b = n_s_sob / n_b_sob 

    #goodness of fit
    error_e = -99
    if n_sig.getVal() > 0 and n_bkg.getVal() > 0:
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

    #move to the higher canvas
    canvas.cd()

    #low pad for drawing the fit


    low_pad = rt.TPad("low_pad", "low_pad", 0, 0, 1, 0.8);
    low_pad.SetTopMargin(.05)
    low_pad.SetBottomMargin(.2)
    low_pad.Draw()

    #move to the low pad for drawing the fit
    low_pad.cd()

    #draw the frame
    frame.Draw()

    canvas.cd()
    #write the latex onto the canvas
    high_pad =  rt.TPad("top_pad","top_pad",0,0.8,1,1)
    high_pad.SetTopMargin(0)
    high_pad.SetBottomMargin(.4)

    high_pad.Draw()    
    high_pad.cd()

    sob_line = "S/B: %.3f #pm %.3f ( %1.1f / %1.1f)" % (s_over_b, error_e, n_s_sob, n_b_sob)
    chi_sq_line = "reduced #chi^{2}: %.2f," % chi2
    gaus_line = "#mu: %.3f, #sigma: %.3f, #mu/err: %.1f" % (mean_val,sigma_val,mu_over_err)
    cut_line = "p_{t,#gamma} > %2.2f, p_{t,#pi^{0}} > %2.2f, S_{4}/S_{9} > %2.2f, iso < %2.2f, Ncri_{1} > %i,  Ncri_{2} > %i" % cut
    #    line4 = "Efficiency %.6f" % eff

    l1_line = ""
    if iSamp == 0:
        l1_line = "Grid Point # %i  Category: %2.2f < |#eta| < %2.2f" % (il1, cat.eta_b, cat.eta_e)

    type_line = None

    if iSamp == 0: type_line = "SEL"

    lat = rt.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(.19)
    lat.SetTextColor(1)    

    ymin = .8
    xmin = .025
    ypass = .23
    
    lat.DrawLatex(xmin,ymin, l1_line)
    lat.DrawLatex(xmin,ymin-ypass, sob_line)
    lat.DrawLatex(xmin,ymin-2*ypass, chi_sq_line + " " + gaus_line)
    lat.DrawLatex(xmin,ymin-3*ypass, cut_line)

    big_text = rt.TLatex()
    big_text.SetNDC()
    big_text.SetTextFont(42)
    big_text.SetTextSize(.55)
    big_text.SetTextColor(1)    

    big_text.DrawLatex(xmin + .7, ymin-1.5*ypass, type_line)

    fit_res = fit_result(result, canvas, n_s_sob, n_b_sob, s_over_b, chi2, error_e, mean_val, sigma_val, mu_over_err, eff)
 
    return fit_res


#################
## MAIN METHOD ##
##             ## 
#################

file = rt.TFile(options.filename)
tree = file.Get("Tree_HLT")

outfile = rt.TFile(options.outfilename, "RECREATE")

grid = build_grid()
ntot = len(grid)

#trim the grid down if necessary
p1 = options.GRID_BEGIN
p2 = options.GRID_END

grid_list = None

if options.GRID_LIST != "no_list":
    #parse the lines from the file
    lines = open(options.GRID_LIST)
    lines = map( lambda x:x.rstrip("\n") , lines)
    list = []
    #add them to the array 
    for ii in lines: list.append(int(ii))

    grid_list = list
elif p1 != -1 and p2 != -1:    
    grid_list = range(p1,p2+1)
else:
    grid_list = range(0,ntot)
    print "-- Warning! Scanning Full Grid --"                                    

print "-- Analyzing # %i of total %i points -- " % (len(grid_list), ntot)

#build the analysis from the grid and tree data
analysis = analysis(grid, tree, grid_list)

#define the categories
analysis.add_category(0, 1)
analysis.add_category(1, 1.4442)

#initialize the grid points
analysis.initialize_grid()

#build the histograms for each category at each grid point
if options.pickle == "no_file":
    #build the histograms
    analysis.build_grid_hists()

    #build the name for the new pickle object
    pickle_out = options.outfilename[:-5] + "_hists_list.p"
    if p1 != -1 and p2 != -1:
        pickle_out = options.outfilename[:-5] + "_hists_%i_%i.p" % (p1, p2)

    #store the histograms
    analysis.pickle_hists(pickle_out)

else:
    analysis.un_pickle_list(options.pickle)

#fit the histograms
analysis.fit_grid_hists()    

#save the output
analysis.save_category_data(outfile)

analysis.build_sob_vs_s(outfile)
