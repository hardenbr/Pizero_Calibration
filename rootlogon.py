import ROOT as rt
import array

def style():
  font = 42
  rt.gStyle.SetLegendFillColor(rt.kWhite)

  rt.gStyle.SetTitleX(0.3)
  rt.gStyle.SetTitleW(0.4)
  rt.gStyle.SetCanvasBorderMode(0)
  rt.gStyle.SetCanvasColor(rt.kWhite)
  rt.gStyle.SetCanvasDefH(600)
  rt.gStyle.SetCanvasDefW(800)
  rt.gStyle.SetCanvasDefX(0)
  rt.gStyle.SetCanvasDefY(0)
  rt.gStyle.SetPadBorderMode(0)
  rt.gStyle.SetPadColor(rt.kWhite)
  rt.gStyle.SetPadGridX(False)
  rt.gStyle.SetPadGridY(False)
  rt.gStyle.SetGridColor(0)
  rt.gStyle.SetGridStyle(3)
  rt.gStyle.SetGridWidth(1)
  rt.gStyle.SetFrameBorderMode(0)
  rt.gStyle.SetFrameBorderSize(1)
  rt.gStyle.SetFrameFillColor(0)
  rt.gStyle.SetFrameFillStyle(0)
  rt.gStyle.SetFrameLineColor(1)
  rt.gStyle.SetFrameLineStyle(1)
  rt.gStyle.SetFrameLineWidth(1)

  rt.gStyle.SetPaperSize(20,26)
  rt.gStyle.SetPadTopMargin(0.1)
  rt.gStyle.SetPadRightMargin(0.10)
  rt.gStyle.SetPadBottomMargin(0.2)
  rt.gStyle.SetPadLeftMargin(0.15)

  rt.gStyle.SetTitleBorderSize(0)
  rt.gStyle.SetTitleFont(font,"xyz")
  rt.gStyle.SetTitleFont(font," ")  
  rt.gStyle.SetTitleSize(0.06,"xyz")
  rt.gStyle.SetTitleSize(0.1," ")   
  rt.gStyle.SetLabelFont(font,"xyz")
  rt.gStyle.SetLabelSize(0.05,"xyz")
  rt.gStyle.SetLabelColor(1,"xyz")
  rt.gStyle.SetTextFont(font)
  rt.gStyle.SetTextSize(0.08)
  rt.gStyle.SetStatFont(font)
  rt.gStyle.SetTitleX(0.5)
  rt.gStyle.SetTitleW(0.4)
  rt.gStyle.SetPadTickX(1)
  rt.gStyle.SetPadTickY(1)


  NRGBs = 5
  NCont = 255
  stops = array.array("d",[0.00, 0.34, 0.61, 0.84, 1.00 ])
  red = array.array("d",[0.00, 0.00, 0.87, 1.00, 0.51 ])
  green = array.array("d",[ 0.00, 0.81, 1.00, 0.20, 0.00 ])
  blue = array.array("d",[ 0.51, 1.00, 0.12, 0.00, 0.00 ])

  rt.TColor.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont)
  rt.gStyle.SetNumberContours(NCont)
