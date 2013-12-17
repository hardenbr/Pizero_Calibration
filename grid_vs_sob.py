import sys, array
import ROOT as rt

file = open(sys.argv[1])
lines = file.readlines()
lines_stripped = map(lambda(x):x.rstrip("\n"),lines)

grid_points = []
sob_points = []

for line in lines_stripped:
    split = line.split("\t")
    iev = int(split[0].lstrip("@@"))
    sob = float(split[4])

    grid_points.append(iev)
    sob_points.append(sob)

print grid_points, sob_points

x = array.array('f',grid_points)
y = array.array('f',sob_points)

gr = rt.TGraph(len(grid_points)-1,x,y)

gr.Draw("AC*")

raw_input("Press Enter")
