import ROOT as rt
import sys
import CMS_lumi
import rootlogon

rootlogon.style()

print "draw_canvas.py [file] [grid number] [eb/ee] [category]"

(name, file_name, gid, eb_ee, cat) = sys.argv


file = rt.TFile(file_name)


eta_b = eta_e = 0

if eb_ee == "eb":
    eta_b = 0 if cat == "0" else 100
    eta_e = 100 if cat =="0" else 144
else:
    eta_b = 150 if cat == "0" else 180
    eta_e = 180 if cat =="0" else 200

canv_name = "canvas_%s_%i_%i" % (gid, eta_b, eta_e)

print "CANVAS NAME:", canv_name

canv = file.Get(canv_name)

canv.Draw()

CMS_lumi.CMS_lumi(canv,4,11)

raw_input("RAW INPUT")
