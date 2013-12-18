import sys, array
import ROOT as rt
import rootlogon
rootlogon.style()

file = open(sys.argv[1])
lines = file.readlines()
lines_stripped = map(lambda(x):x.rstrip("\n"),lines)

grid_points = []
sob_points = []
sob_e_points = []
sob_minus_e_points = []
for line in lines_stripped:
    split = line.split("\t")
    iev = int(split[0].lstrip("@@"))
    sob = float(split[4])
    sob_e = float(split[6])

    grid_points.append(iev)
    sob_points.append(sob)
    sob_minus_e_points.append(sob-sob_e)

print grid_points, sob_points

x = array.array('f',grid_points)
y = array.array('f',sob_minus_e_points)

gr = rt.TGraph(len(grid_points)-1,x,y)
gr.Sort()
gr.Draw("A*")

print "MAX POINT SOB-SOB_E", max(sob_minus_e_points)
print "MAX POINT IEV", grid_points[sob_minus_e_points.index(max(sob_minus_e_points))]

raw_input("Press Enter")
