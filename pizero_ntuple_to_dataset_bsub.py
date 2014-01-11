##Generates config files, and batch submission batch.src files with appropriate directories
from  optparse  import OptionParser
import sys, os

parser = OptionParser()

parser.add_option("-f", "--file_list", dest="list",
		                    help="list of files to analyze",
		                    action="store",type="string")

parser.add_option("-o", "--output", dest="output",
		                    help="output destination",
		                    action="store",type="string")

parser.add_option("-j", "--jobs", dest="jobs",
		                    help="number of jobs to run at once",
		                    action="store",type="string")

parser.add_option("--eta_b", dest="ETA_BEGIN",
                  help="minimum of eta range",
                  action="store",type="float",default=0)

parser.add_option("--eta_e", dest="ETA_END",
                  help="maximum of eta range.",
                  action="store",type="float",default=1000)


(options, args) = parser.parse_args()

parser.print_help()

#get some names and paths
pwd = os.getenv("PWD")
output_dir = options.output

#add the style file
os.system("cp rootlogon.py %s" % output_dir)


#OPEN THE FILE_LIST
file_open = open(options.list,"r")

#READ THE LINES IN
files =file_open.readlines()

#Count how many config files we will need
#there better not be extra lines in the lists!
#n_reco = len(reco_lines)
n_file = len(files)

files = map(lambda(x):x.rstrip("\n"),files)

#remove old instances of the output directory
os.system("rm -r " + output_dir)
#build the directories
os.system("mkdir " + output_dir)
os.system("mkdir " + output_dir+"/logs")
os.system("mkdir " + output_dir+"/src")
os.system("mkdir " + output_dir+"/res")
os.system("cp ~/josh_scripts/buildgrid.py %s" % output_dir)

commands = []


for ii in range(n_file):

	#generate the commands to run on the raw and convert to trees
        output_name = output_dir+"/res/roo" + files[ii].split("/")[-1] 
	
	eta_b = options.ETA_BEGIN
	eta_e = options.ETA_END
	cmd = ""

	#if the eta range is specified pass this to buildgrid and change the output name of the file
	if eta_b != 0 or eta_e < 10:
		output_name != "_%2.2f_%2.2f" % (eta_b,eta_e)
		cmd = "python %s/buildgrid.py -f %s -o %s --eta_b %f --eta_e %f" % (output_dir, files[ii], output_name, eta_b, eta_e)
	else:
		cmd = "python %s/buildgrid.py -f %s -o %s" % (output_dir, files[ii], output_name)

	commands.append(cmd)

njobs = options.jobs
job = 0
file_counter = 0
cmds_per_job = int(len(commands)) / int(njobs) #integer division
bsub_cmds = []
bsub_file = None

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
