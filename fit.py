import math
import ROOT as rt
from ROOT import gSystem, RooArgSet
gSystem.Load('libRooFit')
from ROOT import RooFit

file = rt.TFile("test.root")

t1 = file.Get("tree_1").reduce("mpizero < .2 && mpizero > .08")

file.Close()
outfile = rt.TFile("outframe.root","RECREATE")


#model the signal as a guassian
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

result.Print()

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

#frame.GetXaxis().SetTitleOffset(1.5)

frame.Draw()


n_s_sob = n_sig.getVal() * .9546
n_b_sob = n_bkg.getVal() * bkg_scale
s_over_b = n_s_sob / n_b_sob 
chi2 = frame.chiSquare()
error_e = s_over_b * math.sqrt(math.pow(n_sig.getError() / n_sig.getVal(), 2) + math.pow(n_bkg.getError() / n_bkg.getVal(), 2))

#print relavent signal calculations
#print "normSig", normSig
print "normBkg_sob", normBkg_sob
print "normBkg_full", normBkg_full
print "signal",n_s_sob,"bkg", n_b_sob
print "sob:", s_over_b
print "chi^2:", chi2
print "error_e:", error_e

frame.Write()


raw_input("enter to quit")
