##Generates config files, and batch submission batch.src files with appropriate directories
from  optparse  import OptionParser
import sys, os
import numpy as np
import itertools, math
parser = OptionParser()


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


parser.add_option("-l", "--file_list", dest="list",
		                    help="list of files to analyze",
		                    action="store",type="string")

parser.add_option("-o", "--output", dest="output",
		                    help="output destination",
		                    action="store",type="string")

parser.add_option("-j", "--jobs", dest="jobs",
		                    help="number of jobs to run at once",
		                    action="store",type="string")

parser.add_option("-g", "--gjobs", dest="gjob",
		                    help="number of divisions of grid points ran per buildgrid.py call",
		                    action="store",type="string")


(options, args) = parser.parse_args()

parser.print_help()

#get some names and paths
pwd = os.getenv("PWD")
output_dir = options.output

#remove old instances of the output directory
os.system("rm -r " + output_dir)
#build the directories
os.system("mkdir " + output_dir)
os.system("mkdir " + output_dir+"/logs")
os.system("mkdir " + output_dir+"/src")
os.system("mkdir " + output_dir+"/res")
os.system("cp buildgrid.py %s" % output_dir)
os.system("cp ~/rootlogon.py %s" % output_dir)

commands = []

njobs = int(options.jobs)
gjobs = int(options.gjob)
pi_grid = build_grids()[0]

npoints = []

p_per_job = int(len(pi_grid)) / int(gjobs)
#group the scan based on the number of jobs we want
for ii in range(gjobs):
	npoints.append((ii*p_per_job, (ii+1)*p_per_job))

#the last job takes all the extra jobs (if they exist)
last_job_i = (gjobs) * p_per_job
last_job_f = (gjobs)  * p_per_job + len(pi_grid) % gjobs
if len(pi_grid) % gjobs != 0: npoints.append((last_job_i,last_job_f))

#make a command for each interval in the scan
for ii in npoints:
	#generate the commands to run on the raw and convert to trees
        output_name = output_dir+"res/grid_%i_%i.root" % ii
	cmd = "python %s/buildgrid.py -d %s -o %s -b %i -e %i \n" % (output_dir, pwd+"/"+options.list , output_name, ii[0], ii[1])
	commands.append(cmd)

print len(commands),njobs

job = 0
file_counter = 0
cmds_per_job = int(len(commands)) / int(njobs) #integer division
bsub_cmds = []
bsub_file = None

print "commands per job", cmds_per_job

#loop over all commands
for ii in range(len(commands)):
	#make sure the extra jobs are added to the last file
	if (file_counter == 0) and (job != njobs):		
		#write a new job
		bsub_file_name = output_dir+"/src/bsub_%i.src" % job
		bsub_file = open(bsub_file_name,"a")
		bsub_file.write('#!/bin/bash\n')
		bsub_file.write("cd /afs/cern.ch/user/h/hardenbr/CMSSW_6_2_0/\n")
		bsub_file.write("export SCRAM_ARCH=slc5_amd64_gcc462 \n")
		bsub_file.write("eval `scramv1 ru -sh`\n")

		cmd="bsub -q 1nd "
		cmd+= "-o %s/logs/job_%i.log " % (output_dir ,job)
		cmd+="source "+ bsub_file_name

		bsub_cmds.append(cmd)

		#increment the counter
		job+=1

	bsub_file.write(commands[ii])
	
	#if we are at the number of commands per job reset the counter
	file_counter+=1
	if file_counter == cmds_per_job:file_counter = 0		
	if job == njobs:
		continue

#print out the submission commands
for ii in bsub_cmds:print ii
