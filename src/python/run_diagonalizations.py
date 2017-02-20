#! /usr/bin/env python2
################################################################################
##
##		Author: Simon Davenport	
##
##		Last modified 06/01/2015
##
##      Purpose: Generate plots of the lowest lying interaction 
##      eigenstates of Baur and Cooper's optical flux lattice model
##      as a function of various model parameters
##
################################################################################

################################################################################
################################################################################
#   Import plot functions
import matplotlib
print ("matplotlib version:",matplotlib.__version__)
##	turn off X-forwarding [plt.show() disabled]
##matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.cm as cm
##from matplotlib import animation
from matplotlib.widgets import Slider
##  Import some 3d plotting tools
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.ticker import MultipleLocator
##  Import functions to handle math and arrays
import numpy as np
import math
import cmath
import scipy.special
import copy
##  Import functions for running scripts
import os
##  Command line argument parser
import argparse
##  SQL interface
import sqlite3 as lite
##  Import curve fitting functions
from scipy.optimize import leastsq
from scipy.optimize import curve_fit

################################################################################
##  Define plot presentation information

##  specify figure fonts

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 40}

##	set latex fonts	

matplotlib.rc('text', usetex=True)	

FILE_EXT='.pdf'		##	plot file extension

##  Colour-blindness friendly palette       

##  15 colours are defined here

COLOUR_LIST = ['#bd2309', '#bbb12d', '#1480fa', '#14fa2f', '#000000',
                  '#faf214', '#2edfea', '#ea2ec4', '#ea2e40', '#cdcdcd',
                  '#577a4d', '#2e46c0', '#f59422', '#219774', '#8086d9']

COLOUR_LIST1 = ['#762A83','#9970AB','#C2A5CF','#E7D4E8','#F7F7F7','#D9F0D3','#ACD39E','#5AAE61','#1B7837']

COLOUR_MAP='gist_stern'

##  MARKER LIST - 12 markers here

MARKER_LIST =['s','x','+','^','*','d','>','<','h','8','1','3']

################################################################################
##    Extrapolation function in a+bx^2+cx^4
##
##    Use as fitData,errors = QuarticFit(xData,yData)
##
def QuarticFit(xData,yData):    
    
    xData=np.transpose(xData)
    yData=np.transpose(yData)

    fitfunc = lambda x,a,b,c: a+b*x**(2)+c*x**(4)
    
    params,covar=curve_fit(fitfunc,xData,yData,None)
    
    fitdata=fitfunc(xData,params[0],params[1],params[2])
    
    output='\n\tFitted data:\n'
    
    for i in range(0,len(xData)):
    
        output+='\n\t'+str(xData[i])+'\t'+str(yData[i])+'\t'+str(fitdata[i])

    output+='\n\n\tQuartic fit parameter values:\n'
    output+='\n\ta\t'+str(params[0])+' +/- '+str(math.sqrt(covar[0][0]))
    output+='\n\tb\t'+str(params[1])+' +/- '+str(math.sqrt(covar[1][1]))
    output+='\n\tc\t'+str(params[2])+' +/- '+str(math.sqrt(covar[2][2]))

    print (output)

    return fitdata,math.sqrt(covar[0][0])

################################################################################
################################################################################
##  Define a discrete colour map for use in plot colour bars
##
def DiscreteColourMap(valueList):

    ##  Find the unique set of values in the list and map these to a 
    ##  discrete colour map

    colourList = []

    seen = {}
    unique = []
    
    for item in valueList:
        marker = item

        if marker in seen: continue
        seen[marker] = 1
        unique.append(item)
    
    if len(unique)>15:
        print("TOO MANY COLOURS REQUIRED!")
        quit()
        
    ##dictionary = dict(zip(unique,COLOUR_LIST[:len(unique)]))

    ##return [dictionary[k] for k in valueList]

    unique.sort()

    dictionary = dict(zip(unique,range(0,len(unique))))

    ##  Define a colour map

    colList = colour_list[:len(unique)]

    myColourMap = matplotlib.colors.ListedColormap(colList)

    # define the bins and normalize
    bounds = np.linspace(0,len(colList),len(colList)+1)
    myNorm = matplotlib.colors.BoundaryNorm(bounds, myColourMap.N)
    
    ##print(COLOUR_LIST[:len(unique)])
    ##print([dictionary[k] for k in valueList])
    
    return [dictionary[k] for k in valueList],myColourMap,bounds,myNorm

################################################################################
##  Define directories to be accessed

PROGRAM_FOLDER = '/rscratch/scd51/bin/'

##  '~/'

##'/home/scd51/bin/'

## '/usr/bin'

## 

DATA_FOLDER = '/scratch/scd51/physics_data/optical_flux_lattice_diagonalization/'

PLOTS_FOLDER   = DATA_FOLDER+'/plots/'

SQL_TABLE_NAME = 'Parameters'

################################################################################
##  Define program run information

DIAGONALIZATION_PROGRAM_NAME =  'optical_flux_lattice_diagonalization'
SQL_GEN_PROGARM_NAME =          'optical_flux_lattice_sql_gen'
SINGLE_PARTICLE_PROGRAM_NAME =  'optical_flux_lattice_single_particle_model'
MPI_CORES = 8
SUBMISSION_FILE_NAME =          'current_submission.darwin'

##  Basic program run
def RunDiag(nbrCores,options):

	command = ' nice mpirun -np '+str(nbrCores)+' '+PROGRAM_FOLDER+DIAGONALIZATION_PROGRAM_NAME+options

	print("\n\tRUNNING:   \n"+command)

	os.system(command)

##  Run sql database generator
def RunSqlGen(options):

	command = ' nice mpirun -np 1 '+PROGRAM_FOLDER+SQL_GEN_PROGARM_NAME+options

	print("\n\tRUNNING:   \n"+command)

	os.system(command)

##  Run single particle model calculation
def RunSingleParticle(nbrCores,options):

	command = 'mpirun -np '+str(nbrCores)+' '+PROGRAM_FOLDER+SINGLE_PARTICLE_PROGRAM_NAME+options

	print("\n\tRUNNING:   \n"+command)

	os.system(command)

##  Darwin cluster submission script
def ClusterRunCommand(options,jobLength):

	submissionScript='''#!/bin/bash
#!
#SBATCH -J "'''+str(jobLength)+''' job" -A COOPER --nodes=1 --ntasks=16 --time='''+str(jobLength)+''' -p sandybridge

numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\\1/')

. /etc/profile.d/modules.sh

application="'''+PROGRAM_FOLDER+DIAGONALIZATION_PROGRAM_NAME+'''"

#! Run options for the application:
options="'''+options+'''"

workdir="$SLURM_SUBMIT_DIR"

np=$[${numnodes}*${mpi_tasks_per_node}]

export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets

CMD="mpirun -ppn $mpi_tasks_per_node -np $np $application $options"

###############################################################
### You should not have to change anything below this line ####
###############################################################

cd $workdir
echo -e "Changed directory to `pwd`.\\n"

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"q
echo "Current directory: `pwd`"

if [ "$SLURM_JOB_NODELIST" ]; then
	    #! Create a machine file:
	    export NODEFILE=`generate_pbs_nodefile`
	    cat $NODEFILE | uniq > machine.file.$JOBID
	    echo -e "\\nNodes allocated:\\n================"
	    echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi

echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

echo -e "\nExecuting command:\\n==================\\n$CMD\\n"

eval $CMD
	'''

	##	Put the submission script in a file

	print("Running with options "+options)

	print("Job Length "+jobLength)

	try:
		with open(SUBMISSION_FILE_NAME,'w') as file:
			
			file.write(submissionScript)

			file.close()

			print("Submission script placed in file: "+SUBMISSION_FILE_NAME)

	except IOError as e:

		print("ERROR: Cannot generate file: "+SUBMISSION_FILE_NAME)
	
	##	Execute the submission script

	os.system("sbatch "+SUBMISSION_FILE_NAME)

################################################################################
##  Define a class to hold the command line argument data

class arguments:

    def __init__(self):

        ################################################################################
        ##  Define command line arguments to parse

        parser = argparse.ArgumentParser(description='A script to call opticalFluxLatticeDiag for various parameter sets and make some plots')
        
        ##  General options
        
        parser.add_argument('--x',default=4,type=int,help='number of k-space grid steps in x direction')
        parser.add_argument('--y',default=4,type=int,help='number of k-space grid steps in y direction')
        parser.add_argument('--x-list',default="",help='list of number of k-space grid steps in x direction to run multiple system sizes')
        parser.add_argument('--y-list',default="",help='list of number of k-space grid steps in y direction to run multiple system sizes')
        parser.add_argument('--offset-x',default='0.0',help='kx-offset of the position of the k-space grid.')
        parser.add_argument('--offset-y',default='0.0',help='ky-offset of the position of the k-space grid.')
        parser.add_argument('--nbr',default=4,type=int,help='number of particles')
        parser.add_argument('--nbr-list',default="",help='list of number of particles to run multiple system sizes')
        parser.add_argument('--nbr-eigs',default=3,type=int,help='Set number of eigenvalues to calculate/plot')
        parser.add_argument('--folder-label',help='Specify an additional string to append to the input folder name - this will append to the automatically generated folder name.')
        parser.add_argument('--add-path',help='Specify an additional string to prepend to the input folder name')
        parser.add_argument('--v0',default="1.0",help='Set lattice depth parameter V0')
        parser.add_argument('--interaction',default="1.0",help='Set interaction parameter g')
        
        ##  SQL database construction options
        
        parser.add_argument('--init-database',action='store_true',default=False,help='Generate (or append to existing) SQLite database for the given x,y, nbr parameters and for a specified cut throught parameter space.')
        parser.add_argument('--init-offset-database',action='store_true',default=False,help='Generate (or append to existing) SQLite database for the given x,y, nbr parameters and for a specified cut throught offset space.')
        parser.add_argument('--sql-v0-min',default='0.0',help='Specify minimal V0 parameter to be put in SQLite table.')
        parser.add_argument('--sql-v0-step',default='0.2',help='Specify step in V0 parameter to be put in SQLite table.')
        parser.add_argument('--sql-v0-nbr',default='6',help='Specify number of V0 parameter values to be put in SQLite table (starting from sql-v0-min and increasing in steps of sql-v0-nbr.')
        parser.add_argument('--sql-g-min',default='0.0',help='Specify minimal g parameter to be put in SQLite table.')
        parser.add_argument('--sql-g-step',default='0.2',help='Specify step in g parameter to be put in SQLite table.')
        parser.add_argument('--sql-g-nbr',default='6',help='Specify number of g parameter values to be put in SQLite table (starting from sql-g-min and increasing in steps of sql-g-nbr.')
        parser.add_argument('--sql-kx-offset-min',default='0.0',help='Specify minimal kx-offset parameter to be put in SQLite table.')
        parser.add_argument('--sql-kx-offset-step',default='0.2',help='Specify step in kx-offset parameter to be put in SQLite table.')
        parser.add_argument('--sql-kx-offset-nbr',default='6',help='Specify number of kx-offset parameter values to be put in SQLite table (starting from sql-kx-offset-min and increasing in steps of sql-kx-offset-nbr.')
        parser.add_argument('--sql-ky-offset-min',default='0.0',help='Specify minimal ky-offset parameter to be put in SQLite table.')
        parser.add_argument('--sql-ky-offset-step',default='0.2',help='Specify step in ky-offset parameter to be put in SQLite table.')
        parser.add_argument('--sql-ky-offset-nbr',default='6',help='Specify number of ky-offset parameter values to be put in SQLite table (starting from sql-ky-offset-min and increasing in steps of sql-ky-offset-nbr.')
        parser.add_argument('--sql-completed',default='0',help='Specify \'completed\' flags set in the SQL table (set to either 0 or 1)')
        
        ##  Running diagonalization program options
        
        parser.add_argument('--run-diag',action='store_true',default=False,help='Run all diagonalization rouintes for each case appearing the in the SQL database not already calcualted')
        parser.add_argument('--job-length',default="00:30:00",help='Set the Darwin cluster job length paramter in the format hh::mm::ss . Note the maximum length is 12:00:00.')
        parser.add_argument('--sector-list',default="",help='Specify a list of sectors as a list of integer pairs \'kx1 ky1 kx2 ky2\'')
        parser.add_argument('--sql-id-min',default=1,help='Set the lowest sqlid value for the set of parameters to be run.')
        parser.add_argument('--sql-id-max',default=10000,help='Set the highest sqlid value for the set of parameters to be run. If the value set exceeds the number of entries in the SQL table, then run up to the maximum sqlid value.')
        
        ##  Post-run analysis options
        
        parser.add_argument('--density-density-analysis',action='store_true',default=False,help='Set to calculate the density-density function applied to a specified map')
        parser.add_argument('--calculate-pomeranchuk',action='store_true',default=False,help='Set to apply the density-density function to a pomeranchuck map')
        parser.add_argument('--calculate-susceptibility',action='store_true',default=False,help='Set to apply the density-density function to a magnetization map')
        parser.add_argument('--calculate-r',action='store_true',default=False,help='Set to apply the density-density function to a R-60 order parameter map')
        parser.add_argument('--analyse-single-particle-levels',action='store_true',default=False,help='Perform some analysis of the single-particle energy levels. e.g. make a list of the minimal energy gaps.')
        parser.add_argument('--analyse-most-probable',action='store_true',default=False,help='Plot the configurations corresponding to the highest amplitude fock states.')
        parser.add_argument('--calculate-gap',action='store_true',default=False,help='Determine the energy gap between the highest energy low-lying state, and the lowest lying continuum state.')
        parser.add_argument('--gap-energy',default="10.0",help='Set a value to specify the approximate centre of the energy gap.')
        
        ##  Plot options
        
        parser.add_argument('--formatting',default=0,type=int,help='set option to 1,2 in order to turn on plot formatting for published version.')
        parser.add_argument('--sector',default=[0,0],nargs=2,help='A pair of integers labelling the kx_tot,ky_tot sector to be plotted')
        parser.add_argument('--cut-v0',action='store_true',default=False,help='plot energy levels for different values of V0 for fixed interaction')
        parser.add_argument('--cut-g',action='store_true',default=False,help='plot energy levels for different values of interaction for fixed V0')
        parser.add_argument('--plot-energy',action='store_true',default=False,help='make plot of energy levels for a specified cross section of parameter space')
        parser.add_argument('--plot-probability',action='store_true',default=False,help='option to generate an interactive plot showing the evolution of the probability distribution for a specified cut through parameter space')
        parser.add_argument('--plot-combined-observable',action='store_true',default=False,help='make a plot of a specified observable of the ground state, regardless of sector, over a range of parameter space.')
        parser.add_argument('--observe-susceptibility',action='store_true',default=False,help='plot the susceptibility in the combined observable plot')
        parser.add_argument('--observe-pomer',action='store_true',default=False,help='plot the Pommeramchuk parameter in the combined observable plot')
        parser.add_argument('--observe-participation',action='store_true',default=False,help='plot the participation ratio in the combined observable plot')
        parser.add_argument('--observe-r',action='store_true',default=False,help='plot the R-60 symmetry order parameter')
        parser.add_argument('--observe-rotational',action='store_true',default=False,help='plot the R-180 symmetry density-density function peak')
        parser.add_argument('--plot-magnetization-map',action='store_true',default=False,help='make a plot of the single-particle magnetization map')
        parser.add_argument('--plot-overlaps',action='store_true',default=False,help='Make plots of the overlaps with selected states as a function of model parameters.')
        parser.add_argument('--plot-spatial-wavefunction',action='store_true',default=False,help='Make plots of the spatial wave function orbitals.')
        parser.add_argument('--plot-energy-derivative',action='store_true',default=False,help='Plot the second derivative of the energy vs filling factor')
        parser.add_argument('--plot-energy-vs-offset',action='store_true',default=False,help='Plot the energy levels vs the value of the lattice offset parameter')
        parser.add_argument('--plot-energy-vs-sector',action='store_true',default=False,help='Plot the energy levels vs unfolded sector kx+nx*ky')
        parser.add_argument('--plot-combined-observable-vs-sector',action='store_true',default=False,help='Plot the an observable vs unfolded sector kx+nx*ky')
        parser.add_argument('--plot-translational-density-density',action='store_true',default=False,help='Plot the density-density function sum_{k1,k2} <c^+_{k2-G}c_k2c^+_{k1+G}c_k1> as a function of G.')
        parser.add_argument('--plot-rotational-density-density',action='store_true',default=False,help='Plot the density-density function sum_{k1,k2} <c^+_{~R60 k2}c_k2c^+_{R60 k1}c_k1> as a function of the 6-fold rotation angle.')
        parser.add_argument('--plot-theta-band-width',action='store_true',default=False,help='Make a plot of the band width/band gap in the single particle model as a function of both V0 and theta')
        parser.add_argument('--plot-epsilon-band-width',action='store_true',default=False,help='Make a plot of the band width/band gap in the single particle model as a function of both V0 and epsilon')
        parser.add_argument('--epsilon',default='0.4',help='Set the value of the epsilon parameter to use from readling from SQL database')
        parser.add_argument('--theta',default='0.3',help='Set the value of the theta parameter to use from readling from SQL database')
        
        ##  Process the arguments

        args = parser.parse_args(namespace=self)

        if self.folder_label == None:
            self.folder_label = ''
            
        if self.add_path == None:
            self.add_path = ''

        ##  Convert the #-list arguments from strings into lists of values

        if self.x_list == "" and self.y_list == "" and self.nbr_list == "":
            self.systemList = [[self.x,self.y,self.nbr]]
        elif self.x_list != "" and self.y_list != "" and self.nbr_list != "":
            
            tempX = self.x_list.split()
            tempY = self.y_list.split()
            tempN = self.nbr_list.split()

            self.systemList = []
            counter = 0
            
            for nbr in tempN:
                
                self.systemList.append([int(tempX[counter]),int(tempY[counter]),int(tempN[counter])])
                counter += 1

        else:
            print("ERROR: inconsistent x-list,y-list and nbr-list arguments!")
            quit()
            self.x_list = []

    def CppArgs(self):
        return ' -x '+str(self.x)+' -y '+str(self.y)+' -n '+str(self.nbr)+' '
        
    def FolderName(self):
        return self.add_path+'data_kx_'+str(self.x)+'_ky_'+str(self.y)+'_n_'+str(self.nbr)+self.folder_label+'/'
    
    def SqlName(self):
        return 'resultKey_kx_'+str(self.x)+'_ky_'+str(self.y)+'_n_'+str(self.nbr)
        
    def SectorName(self):
        return '_sector_'+str(self.sector[0])+'_'+str(self.sector[1])
    
    def SectorCppArgs(self):
        return ' -s '+str(self.sector[0])+' '+str(self.sector[1])+' '
        
    def GetLabel(self):
        return 'kx_'+str(self.x)+'_ky_'+str(self.y)+'_n_'+str(self.nbr)+'_'

################################################################################
##  Define a class to retrieve and manipulate data from an SQLite database
class SQLite:

    def __init__(self,dbName):
        
        self.dbName = dbName
        print(self.dbName)
        
        try:
            self.connectObject  = lite.connect(self.dbName)
        except lite.Error,e:
                print("ERROR:", e.args[0])
                quit()

    ##  Read from SQL database
    def GetData(self,fields,tableName,conditions):
        
        ## From the connection, we get the cursor object. 
        ## The cursor is used to traverse the records. 

        cursorObject = self.connectObject.cursor()  

        ##	Make the sql script to extract the relevant data
        
        sqlCommand="""
                    SELECT
                    """+str(fields)+"""
                    FROM """+str(tableName)+"""
                    WHERE
                    """+str(conditions)+"""
                    """
        try:
            cursorObject.execute(sqlCommand)
        except lite.Error,e:
                print("ERROR:", e.args[0])
                quit()
                
        colNames = [cn[0] for cn in cursorObject.description]
        
        rows = cursorObject.fetchall()
        
        return colNames,rows
    
    ##  Print out retrieved data
    
    def PrintData(self,colNames,rows):
    
        print ("\n\tSet of data retrieved:\n")

        basePrintStr = '{: <20} '
        
        printStr=''
        
        for i in range(0,len(colNames)):
            
            printStr += basePrintStr

        print (printStr.format(*colNames))

        for row in rows:
        
            print (printStr.format(*row))
    
################################################################################
##  Define a modulo function
def Modulo(a,b):
    
    if a >= 0:
        return math.fmod(a,b)
    else:
        return math.fmod(b - abs(math.fmod(a,b)),b)

################################################################################
##  Define a function to get the first binary Hamming number
def FirstBinaryHammingNumber(hammingWeight):
        
            return (1 << hammingWeight) - 1

################################################################################
##  Define a function to cycle through binary Hamming numbers     
def NextHammingNumber64(x):
    
    temp = (x & -x)  ##  Isolates the right-most bit

    return ((x ^ (x + temp)) / temp) >> 2 | (x + temp)

################################################################################
##  Define a function to convert a binary number into an array containing the
##  positions of its non-zero elements
def ConverBinaryToArray(x):

    setPos = []

    i=0
    
    while i<64:
    
        if (1 << i) & x:
        
            setPos.append(i)
    
        i+=1
    
    return setPos

################################################################################
##  Define a function to calculate binomial coefficients
def Binomial(n,r):
    return math.factorial(n) / math.factorial(r) / math.factorial(n-r)

################################################################################
##  Define a function to plot a Hexagonal lattice background within a certian
##  kx,ky parameter range, used to frame many of the plots
def PlotUnitCellBackground(ax=None):
    if ax is None:
        ax = plt.gca()

    sq3 = math.sqrt(3)
    linecol = 'black'

    ax.set_xlabel(r"$k_x$",fontsize=font_size)
    ax.xaxis.set_label_coords(1.05, -0.025)
    ax.set_xlim([-0.7*sq3,0.7*sq3])
    ax.set_xticks([-sq3/2,0,sq3/2]) 
    ax.set_xticklabels([r'$-\frac{\sqrt{3}\pi}{a}$',r'$0$',r'$\frac{\sqrt{3}\pi}{a}$'],fontsize=font_size)
    ax.tick_params(axis='x', pad=x_tick_pad)
    
    ax.set_ylabel(r"$k_y$",fontsize=font_size,rotation=0)
    ax.yaxis.set_label_coords(-0.05,0.98)
    ax.set_ylim([-2.7,0.0])
    ax.set_yticks([-2.5,-1.5,-0.5])
    ax.set_yticklabels([r'$-\frac{5\pi}{a}$',r'$-\frac{3\pi}{a}$',r'$-\frac{\pi}{a}$'],fontsize=font_size)  
    
    ## Draw on Hexagonal lattice background

    ax.plot([-sq3/2,-sq3/2],[-1,-2],'k-', lw=1,c=linecol)
    ax.plot([sq3/2,sq3/2],[-1,-2],'k-', lw=1,c=linecol)
    ax.plot([0,-sq3/2],[-1.0/2,-1.0],'k-',lw=1,c=linecol)
    ax.plot([0,sq3/2],[-1.0/2,-1.0],'k-',lw=1,c=linecol)
    ax.plot([0,-sq3/2],[-2.5,-2],'k-',lw=1,c=linecol)
    ax.plot([0,sq3/2],[-2.5,-2],'k-',lw=1,c=linecol)
    ax.plot([0,0],[-0.5,0],'k-',lw=2,c=linecol)
    ax.plot([0,0],[-2.5,-3],'k-',lw=2,c=linecol)
    #ax.plot([-sq3,-sq3/2],[-1.0/2,-1.0],'k-',lw=2,c=linecol)
    #ax.plot([sq3,sq3/2],[-1.0/2,-1.0],'k-',lw=2,c=linecol)
    #ax.plot([-sq3,-sq3/2],[-2.5,-2],'k-',lw=2,c=linecol)
    #ax.plot([sq3,sq3/2],[-2.5,-2],'k-',lw=2,c=linecol)

################################################################################
##  Define a function to plot a grid of point representing our simulation cell
def PlotSimulationGrid(ax=None):
    if ax is None:
        ax = plt.gca()
        
    ##  Draw a set of points representing the simulation cell sites
        
    xVals = []
    yVals = []
    
    sq3 = math.sqrt(3)

    for x in range(-2*args.x,2*args.x):
        for y in range(-2*args.y,2*args.y):

            xVals.append(float(x)/args.x*(-sq3/2.0)+float(y)/args.y*(sq3/2.0)+float(args.offset_x))
            yVals.append(float(x)/args.x*(-1.5)+float(y)/args.y*(-1.5)+float(args.offset_y))

    ax.scatter(xVals,yVals,marker="o",color=COLOUR_LIST[1],s=250)

################################################################################
##  Define a class to allow interactive browsing of the plot data
class PointBrowser:
    '''
    Click on a point to select and highlight it.Use the 'n'
    and 'p' keys to browse through the next and previous points
    '''
    
    ##  Constructor sets default state
    def __init__(self):
    
        self.lastind = 0

        ##  Selection text set to none
        self.text = main_axes.text(0.05, 0.95, 'Selected: none',transform=main_axes.transAxes, va='top',fontsize=20)
                            
        ##  Prepare plot highlighting   
        
        if args.cut_v0 == 1:    ##  keep g fixed       
            self.selected,  = main_axes.plot([dataV0[0]], [dataEnergy[0]], 'o', ms=30, alpha=0.4,color='red',visible=False)
        
        elif args.cut_g == 1:   ##  keep v0 fixed               
            self.selected,  = main_axes.plot([dataInteraction[0]], [dataEnergy[0]], 'o', ms=30, alpha=0.4,color='red',visible=False)
                                  
    ##  Action to perform on key press
    def onpress(self, event):
        if self.lastind is None: return
        if event.key not in ('n', 'p'): return
        if event.key=='n': inc = 1
        else:  inc = -1

        self.lastind += inc
        self.lastind = np.clip(self.lastind, 0, len(dataEnergy)-1)
        self.update()
    
    ##  Action on mouse click
    def onpick(self, event):
    
       index = event.ind
    
       for i in index:
            self.lastind = i
            self.update()
    
    ##  Plot update functions
    def update(self):
        if self.lastind is None: return

        dataind = self.lastind

        self.selected.set_visible(True)
        
        if args.cut_v0 == 1:    ##  keep g fixed
            self.selected.set_data(dataV0[dataind], dataEnergy[dataind])
        elif args.cut_g == 1:   ##  keep v0 fixed
            self.selected.set_data(dataInteraction[dataind], dataEnergy[dataind])
            
        self.text.set_text('Selected: '+spectrumData.GetText(dataind))
        main_axes.get_figure().canvas.draw()

################################################################################
##  Define a class to build a discrete slider for interactive plots
class DiscreteSlider(Slider):
    """A matplotlib slider widget with discrete steps."""
    def __init__(self, *args, **kwargs):
        """
        Identical to Slider.__init__, except for the new keyword 'allowed_vals'.
        This keyword specifies the allowed positions of the slider
        """
        self.allowed_vals = kwargs.pop('allowed_vals',None)
        self.previous_val = kwargs['valinit']
        Slider.__init__(self, *args, **kwargs)
        if self.allowed_vals==None:
            self.allowed_vals = [self.valmin,self.valmax]

    def set_val(self, val):
    
        dummy = np.zeros(len(self.allowed_vals))
                
        dummy.fill(val)
    
        discrete_val = self.allowed_vals[abs(dummy-self.allowed_vals).argmin()]
        xy = self.poly.xy
        xy[2] = discrete_val, 1
        xy[3] = discrete_val, 0
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % discrete_val)
        if self.drawon: 
            self.ax.figure.canvas.draw()
        self.val = val
        if self.previous_val!=discrete_val:
            self.previous_val = discrete_val
            if not self.eventson: 
                return
            for cid, func in self.observers.iteritems():
                func(discrete_val)

################################################################################
##  Define a class to hold data about each data point in the plot of the energy 
##  spectrum
class SpectrumDataPoint:
    
    ##  Default class parameter values

    energy       = 0.0      ##  Energy value
    interaction  = 0.0      ##  Interaction strength (g_2D) value
    v0           = 0.0      ##  Lattice depth value
    offset       = [[0.0,0.0]]  ##	Simulation grid offset along x,y direction
    degeneracy   = 1        ##  Degeneracy of this energy level
    sectors      = [[0,0]]  ##  Sector label(s)
    eigLabels    = [0]      ##  Label of energy level (s) within the given 
                            ##  sector (0 is lowest level, 1 the next lowest etc.)
    
    tol = 0.00000001        ##  Tolerance for accepting eigenvalues as degenerate
    
    ##  Define an "equals" operator such that objects with the same energy,
    ##  interaction and v) are counted as the same, regardless of the other
    ##  parameter values. Used to check for degenerate states in the energy 
    ##  spectrum
    def __eq__(self, other):
        
        return abs(self.energy-other.energy)<self.tol and self.interaction == other.interaction and self.v0 == other.v0

    ##  Define an "addition" operator such that added objects have combined degeneracy,
    ##  sectors and eigLabels variables, but other parameters are untouched. 
    ##  Used to combine degenerate states  
    def __add__(self, other):
        
        self.degeneracy += other.degeneracy
        self.sectors    += other.sectors
        self.eigLabels  += other.eigLabels
    
    ##  Define a print function   
    def Print(self):
    
        print("Energy = ",self.energy,"  g_{2D} = ",self.interaction,"  V_0 = ",self.v0," degeneracy = ",self.degeneracy,"  sectors = ",self.sectors,"  eig labels = ",self.eigLabels)

    ##  Define a function to return a string containing the class data
    def GetText(self):
    
        return r"Energy = "+str(self.energy)+r" $\,\,g_{2D} = "+str(self.interaction)+r"$ $\,\,V_0 = "+str(self.v0)+r"$"+"\n"+"$\,\,$ degeneracy = "+str(self.degeneracy)+r" $\,\,$ sectors = "+str(self.sectors)+"\n"+r" eig labels = "+str(self.eigLabels)

################################################################################
##  Define a class to hold an array of spectrum data objects, used to define 
##  a full plot of the energy spectrum
class SpectrumDataArray:

    dataArray = []

    ##  We need to define a "destructor" operator that clears the internal
    ##  list of objects explicitly. 

    def __del__(self):
        self.dataArray[:] = []

    ##  Define a function to get the number of entries in the data array
    def GetNbrPoints(self):
        return len(self.dataArray)

    ##  Define an append function such that any existing degenerate state
    ##  with the same g_2D and V_0 values is not appended, but instead
    ##  the equivalent existing point is updated. Otherwise, a new
    ##  point is appended
    
    def Append(self,newTerm):
        
        ##  Check for existing data point

        for point in self.dataArray:

            if point == newTerm:
                
                ##  Update degeneracy of existing point
                
                point += newTerm
                
                return

        ##  Otherwise append the new value to the array
    
        self.dataArray.append(newTerm)

    ##  Define a print function
    
    def Print(self):
    
        for point in self.dataArray:
        
            point.Print()

    ##  Sort the internal array in place
    
    def Sort(self):

        self.dataArray.sort(key=lambda x: x.energy)

    ##  Define a print function
    
    def GetText(self,i):
    
        return self.dataArray[i].GetText()

    ##  Get a numpy array containing all of the energy data
    
    def GetEnergyList(self):
    
        energy = np.empty([len(self.dataArray)])
    
        i = 0 
    
        for point in self.dataArray:
        
            energy[i] = point.energy
            
            i += 1
            
        return energy
    
    ##  Get a numpy array containing all of the energy data including 
    ##  repeated values for degenerate states from different sectors
       
    def GetAllEnergyList(self):
    
        energy = []

        for point in self.dataArray:
            
            for sector in point.sectors:
            
                energy.append(point.energy)

        return energy

    ##  Get a numpy array containing all of the interaction data
    
    def GetInteractionList(self):
    
        interaction = np.empty([len(self.dataArray)])
    
        i = 0 
    
        for point in self.dataArray:
        
            interaction[i] = point.interaction
            
            i += 1
            
        return interaction

    ##  Get a numpy array containing all of the V0 data

    def GetV0List(self):
    
        V0 = np.empty([len(self.dataArray)])
    
        i = 0 
    
        for point in self.dataArray:
        
            V0[i] = point.v0
            
            i += 1
            
        return V0
    
    ##  Get a numpy array containing all of the offset data
    
    def GetOffsetList(self):
    
        offset = np.empty([len(self.dataArray)])
    
        i = 0 
    
        for point in self.dataArray:
        
            offset[i] = point.offset
            
            i += 1
            
        return offset

    ##  Get a list of unfolded sectors
    
    def GetUnfoldedSectorList(self,ny):
    
        sectors = []

        for point in self.dataArray:

            for sector in point.sectors:

                sectors.append(ny*sector[0]+sector[1])
                
        return sectors
    
    ##  Get maximum values of the parameters
    
    def GetMax(self):
        
        return np.max(self.GetEnergyList()),np.max(self.GetInteractionList()),np.max(self.GetV0List())

    ##  Get minimum values of the parameters

    def GetMin(self):
        
        return np.min(self.GetEnergyList()),np.min(self.GetInteractionList()),np.min(self.GetV0List())
        
    ##  Get the energy from a given data set
    
    def GetEnergy(self,index):

        return self.dataArray[index].energy
    
    ##  Get the interaction from a given data set
    
    def GetInteraction(self,index):

        return self.dataArray[index].interaction
    
    ##  Get the v0 value from a given data set
    
    def GetV0(self,index):

        return self.dataArray[index].v0
    
    ##  Get the degeneracy from a given data set
    
    def GetDegeneracy(self,index):

        return self.dataArray[index].degeneracy
    
    ##  Get the list of sectors from a given data set
    
    def GetSectors(self,index):

        return self.dataArray[index].sectors
        
    ##  Get the list of eigenvalue labels from a given data set
    
    def GetEigLabels(self,index):

        return self.dataArray[index].eigLabels

################################################################################
##  Define a class to store single particle energy and magnetization map data 
##  associated with a given ky,ky state
class SingleParticleData:

    energy = 0.0
    magnetization = 0.0
    kx = 0.0
    ky = 0.0

################################################################################
##  Define a class to to store spatial wave function data
class SpatialWaveFunctionData:
    
    def __init__(self,kx,ky,spin):
        self.kx = kx
        self.ky = ky
        self.spin = spin
        self.realSpaceX = []
        self.realSpaceY = []
        self.values = []
    
    def Append(self,kx,ky,spin,rx,ry,value):
        
        ##  Only include a new value of the kx,ky and spin indices match 
        ##  up with the values that the class was initialized with
        
        if kx == self.kx and ky == self.ky and spin == self.spin:          
            self.realSpaceX.append(rx)
            self.realSpaceY.append(ry)
            self.values.append(value)

################################################################################
##  Define a class to store configuration amplitude data, used for reading
##  in overlap vectors and most probable states
class ConfigurationAmplitudeData:

    interaction = 0.0
    v0          = 0.0
    amplitude   = 0.0
    fockState   = 0
    dimX        = 0
    dimY        = 0
    
    ################################################################################
    ##  Class constructor given a set of kx,ky grid dimensions
    def __init__(self,dimX,dimY):
        self.dimX = dimX
        self.dimY = dimY

    ################################################################################
    ##  Define a function to get the index closet to a given kx,ky in the grid
    def GetIndex(self,kx,ky):

        sq3 = math.sqrt(3)
        
        kx -= float(args.offset_x)
        ky -= float(args.offset_y)
        
        i = -int(round(self.dimX/3.0 *(sq3*kx+ky)))
        j = int(round(self.dimY/3.0 *(sq3*kx-ky)))

        i = int(Modulo(i,self.dimX))
        
        j = int(Modulo(j,self.dimY))

        return i*self.dimY+j
    
    ################################################################################
    ##  Define a function to flip the occupation index of the kx,ky state
    def FlipState(self,kx,ky):

        ##  Update the existing fock state by flipping it's occupation bit
        ##  with an XOR operation
        
        self.fockState ^= (1 << self.GetIndex(kx,ky))
    
    ################################################################################
    ##  Define a function to get the total value of e.g. energy or
    ##  magnetization, given an array mapped to the grid of kx,ky points.
    ##  the given array will be masked by the fockState describing the
    ##  single particle occupations
    def GetTotal(self,array):
        
        accum = 0.0
        
        for i in range(0,self.dimX):
            for j in range(0,self.dimY):

                if (1 << (i*self.dimY+j)) & self.fockState:
                
                    accum += array[i*self.dimY+j]
        
        return accum
    
    ################################################################################
    ##  Define a function to convert the fock state into a list of kx,ky
    ##  co-ordinates for the set of occupied states
    def GetConfiguration(self):
        
        xData = []
        yData = []
        
        sq3 = math.sqrt(3)

        for i in range(0,self.dimX):
            for j in range(0,self.dimY):

                if (1 << (i*self.dimY+j)) & self.fockState:

                    xData.append(float(i)/self.dimX*(-sq3/2.0)+float(j)/self.dimY*(sq3/2.0)+float(args.offset_x))
                    yData.append(float(i)/self.dimX*(-1.5)+float(j)/self.dimY*(-1.5)+float(args.offset_y))
        
        dim = len(xData)
        
        ##  Convert to numpy arrays
        xData = np.array(xData)
        yData = np.array(yData)

        ##  Copy into neighbouring unit cells
        
        xData2 = xData
        yData2 = yData

        for i,j in [[2,0],[1,1],[1,-1],[-1,1],[-1,-1],[-2,0]]:

            shiftX=np.zeros(dim)
            shiftY=np.zeros(dim)

            shiftX.fill(i*sq3/2.0)
            shiftY.fill(j*1.5)

            xData2 = np.append(xData2,np.add(xData,shiftX))
            yData2 = np.append(yData2,np.add(yData,shiftY))
        
        return xData2,yData2
   
    ################################################################################
    ##  Define a function to obtain the angular momentum sector associated with the
    ##  current configuration
    def GetSector(self):
    
        iTot = 0
        jTot = 0
    
        for i in range(0,self.dimX):
            for j in range(0,self.dimY):

                if (1 << (i*self.dimY+j)) & self.fockState:
                
                    iTot += i
                    jTot += j
                    
        return Modulo(iTot,self.dimX),Modulo(jTot,self.dimY)

################################################################################
##  Define a function to get the spatial wave function distribution
def ObtainSpatialWaveFunction():

    ##  First check for the file containing the data, read it in
    ##  and arrange the data appropriately. If the file does not
    ##  exist then run the diagonalization program to calculate it
    
    gridSize = "30"
    gridSpacing = "0.25"
    basis = "1"
    
    fileName = DATA_FOLDER+'/spatial_wave_function/spatial_wave_function_kx_'+str(args.x)+'_ky_'+str(args.y)+'_grid_size_'+gridSize
    
    if basis == "1":
        fileName += '_wannier_basis'
        
    fileName += '.dat'

    try:
        with open(fileName) as file:
            pass
    except IOError as e:
    
        RunSingleParticle(MPI_CORES,' --out-path '+DATA_FOLDER+ \
        '/spatial_wave_function/'+' --in-path '+DATA_FOLDER+'/spatial_wave_function/'+\
        ' --plot-x '+str(args.x)+' --plot-y '+str(args.y)+' --x-cut 16 --y-cut 16 -v 3 '+args.SectorCppArgs()+' --calculate-spatial-wavefunctions '+gridSize+' --basis '+basis+'  --grid-spacing '+gridSpacing)
    
    ##  Initialize an array of classes to hold the spatial wave function data
    ##  for each kx,ky set of values
    
    dataArray = []

    for kx in range (0,args.x):
        for ky in range (0,args.y):
            for spin in [0,1]:
                dataArray.append(SpatialWaveFunctionData(kx,ky,spin))
    
    myRange = len(dataArray)
    
    fin = open(fileName)
    
    ##  Skip comment lines
    
    fin.readline()
    fin.readline()
    
    gridSize = int(gridSize)

    for line in fin:
        ##print(line)

        line = line.split()
        
        ##  The first column of the file contains kx*dimY+ky
        
        indexK = int(line[0])
        
        kx = math.floor(indexK/args.y)
        ky = Modulo(indexK,args.y)
        
        ##  The second column of the file contains rx*gridSize+ry
        
        indexR = int(line[1])
        
        rx = math.floor(indexR/gridSize)
        ry = Modulo(indexR,gridSize)
        
        ##  The third column contains the spin index

        indexSpin = int(line[2])
        
        ##  The 4th column of the array contains the 
        ##  wave function probability
        
        value = float(line[3])
        
        ##  Combine the kx,ky and spin indices
        
        index = indexSpin + indexK*2

        dataArray[index].Append(kx,ky,indexSpin,rx,ry,value)

    return dataArray

################################################################################
##  Define a function to read in the single particle magnetization map data
##	or generate it from scratch if unavailable
def ObtainMagnetizationMapData(sqlId):
    
    if sqlId == 0:

        outputPath = DATA_FOLDER+args.add_path

        fileName = outputPath+'magnetization_map_kx_'+str(args.x)+'_ky_'+str(args.y)+'.dat'
    
        print("searching for "+fileName)
    
        try:
            with open(fileName) as file:
                pass
        except IOError as e:

            RunSingleParticle(MPI_CORES,' --out-path '+outputPath+' --in-path '+outputPath + \
            ' --x-grid '+str(args.x)+' --y-grid '+str(args.y)+' --x-cut '+str(16)+' --y-cut '+str(16)+' -v 2 --nbr-bands 1 '+args.SectorCppArgs()+' --calculate-magnetization 1 --x-grid-shift '+args.offset_x+' --y-grid-shift '+args.offset_y)

    else:
        fileName = DATA_FOLDER+args.add_path+'magnetization_map_kx_'+str(args.x)+'_ky_'+str(args.y)+'_id_'+str(sqlId)+'.dat'

        print("searching for "+fileName)

        try:
            with open(fileName) as file:
                pass
        except IOError as e:

            RunSingleParticle(MPI_CORES,' --use-sql 1 --sql-id ' \
            +str(sqlId)+' --out-path '+DATA_FOLDER+args.FolderName() +' --sql-name '+args.SqlName() + \
            ' --sql-table-name '+SQL_TABLE_NAME+' --in-path '+DATA_FOLDER+args.FolderName()+\
            ' --x-grid '+str(args.x)+' --y-grid '+str(args.y)+' --x-cut '+str(12)+' --y-cut '+str(12)+' -v 2 --nbr-bands 1 --calculate-magnetization 1 --x-grid-shift '+args.offset_x+' --y-grid-shift '+args.offset_y)

    ##  Read in data from the file
        
    fin  = open(fileName)

    xData      = []
    yData      = []
    magnetizationData = []

    dimX = int(fin.readline()) 
    dimY = int(fin.readline())

    ##  Skip comment line
    fin.readline()

    for line in fin:
        line = line.split()
        
        xData.append(float(line[0]))
        yData.append(float(line[1]))
        magnetizationData.append(float(line[2]))

    fin.close()

    ##  Convert to np arrays

    xData = np.array(xData)
    yData = np.array(yData)
    magnetizationData = np.array(magnetizationData)
    
    ##print(magnetizationData)

    return xData,yData,magnetizationData

################################################################################
##  Define a function to read in the single particle energy data
def ObtainSingleParticleEnergyData(sqlId):

    if sqlId == 0:
		
        outputPath = DATA_FOLDER+'/single_particle_energy_data/'

        fileName = outputPath+'band_structure_kx_'+str(args.x)+'_ky_'+str(args.y)+'.dat'

        RunSingleParticle(MPI_CORES,' --out-path '+outputPath + ' --in-path '+outputPath+\
            ' --x-grid '+str(args.x)+' --y-grid '+str(args.y)+' --x-cut '+str(16)+' --y-cut '+str(16)+' -v 2 --nbr-bands 1 --calculate-bandstructure 1 --x-grid-shift '+args.offset_x+' --y-grid-shift '+args.offset_y)

    else:
        fileName = DATA_FOLDER+args.FolderName()+'band_structure_kx_'+str(args.x)+'_ky_'+str(args.y)+'_id_'+str(sqlId)+'.dat'
        
        try:
            with open(fileName) as file:
                pass
        except IOError as e:

            RunSingleParticle(MPI_CORES,' --use-sql 1 --sql-id ' \
            +str(sqlId)+' --out-path '+DATA_FOLDER+args.FolderName() +' --sql-name '+args.SqlName() + \
            ' --sql-table-name '+SQL_TABLE_NAME+' --in-path '+DATA_FOLDER+args.FolderName()+\
            ' --x-grid '+str(args.x)+' --y-grid '+str(args.y)+' --x-cut '+str(16)+' --y-cut '+str(16)+' -v 2 --nbr-bands 1 --calculate-bandstructure 1 --x-grid-shift '+args.offset_x+' --y-grid-shift '+args.offset_y)

    ##  Read in data from the file
        
    fin  = open(fileName)

    xData      = []
    yData      = []
    energyData = []
    
    dimX = int(fin.readline()) 
    dimY = int(fin.readline())

    ##  Skip comment line
    fin.readline()
    
    for line in fin:
        line = line.split()
        
        xData.append(float(line[0]))
        yData.append(float(line[1]))
        energyData.append(float(line[2]))

    fin.close()

    ##  Convert to np arrays
    
    xData = np.array(xData)
    yData = np.array(yData)
    energyData = np.array(energyData)

    return xData,yData,energyData

################################################################################
##  Define a function to turn the data from a single simulation cell into 
##  data for some additional neighbouring simulation cells also
def ExpandCellData(xData,yData,array):

    ##  Generate periodic pattern by shifting kx,ky

    sq3 = math.sqrt(3)
    
    xData2 = []
    yData2 = []
    
    array2 = array[:]
    
    for i in range(0,8):
        array2 = np.append(array2,array)

    for i in [-1,0,1]:
        for j in [-1,0,1]:
            for x in range(i*args.x,(i+1)*args.x):
                for y in range(j*args.y,(j+1)*args.y):
            
                    x0 = float(x)/args.x*(-sq3/2.0)+float(y)/args.y*(sq3/2.0)+float(args.offset_x)
                    y0 = float(x)/args.x*(-1.5)+float(y)/args.y*(-1.5)+float(args.offset_y)
                    
                    xData2.append(x0)
                    yData2.append(y0)

    xData2 = np.array(xData2)
    yData2 = np.array(yData2)
    
    return xData2,yData2,array2

################################################################################
##  Define a function to read in density-density data and filter it with a
##	given k-space map
def ApplyDensityDensityToMap(fileName,kSpaceMap):
    
    totalMappedValue = 0

    try:
        with open(fileName) as file:
            fin = open(fileName)

            ##	Read in the first two entries which are the k-space dimensions

            dimX = int(fin.readline()) 
            dimY = int(fin.readline())

            ##	Skip a comment line

            fin.readline()

            ##	Read in the density-density function and apply it to the given k-space map

            for x in range(0,dimX):
                for y in range(0,dimY):
                    for x1 in range(0,dimX):
                        for y1 in range(0,dimY):

	                        ##	Read in density-density function value from the file
	                        ##	and multiply it by the mapped values

	                        line = fin.readline()

	                        line = line.split()

	                        totalMappedValue += kSpaceMap[x*dimY+y]*kSpaceMap[x1*dimY+y1]*float(line[4])
	                        
	                        

    except IOError as e:
        print("WARNING: FILE "+fileName+" NOT FOUND")

    return totalMappedValue

################################################################################
######      MAIN        ########################################################

if __name__=='__main__':

    ############################################################################
    ##  Set up command line arguments
       
    args = arguments()

################################################################################
######      GENERATE SQL DATABASE OF PARAMETER SPACE CUT    ####################

    if args.init_database:

        print ("\n========== Initializing SQLite database "\
        +"(OR appending to existing database) ==========")

        os.system('mkdir '+DATA_FOLDER+args.FolderName())

        RunSqlGen(' --build-sql-table 1 '+ \
        args.CppArgs()+' --out-path '+DATA_FOLDER+args.FolderName()+\
        ' --sql-name '+args.SqlName()+' --sql-table-name '+SQL_TABLE_NAME+\
        ' -v 3 --sql-v0-min '+args.sql_v0_min+' --sql-v0-step '+args.sql_v0_step+\
        ' --sql-v0-nbr '+args.sql_v0_nbr+' --sql-g-min '+args.sql_g_min+\
        ' --sql-g-step '+args.sql_g_step+' --sql-g-nbr '+args.sql_g_nbr+\
        ' --sql-completed '+args.sql_completed+' --kx-shift '+args.offset_x+\
        ' --ky-shift '+args.offset_y)
    
    if args.init_offset_database:

        print ("\n========== Initializing SQLite offset database "\
        +"(OR appending to existing database) ==========")

        os.system('mkdir '+DATA_FOLDER+args.FolderName())

        RunSqlGen(' --build-sql-table-offset 1 '+ \
        args.CppArgs()+' --out-path '+DATA_FOLDER+args.FolderName()+\
        ' --sql-name '+args.SqlName()+' --sql-table-name '+SQL_TABLE_NAME+\
        ' -v 3 --sql-v0-min '+args.v0+' --sql-g-min '+args.interaction+\
        ' --sql-kx-shift-min '+args.sql_kx_offset_min+\
        ' --sql-kx-shift-step '+args.sql_kx_offset_step+\
        ' --sql-kx-shift-nbr '+args.sql_kx_offset_nbr+\
        ' --sql-ky-shift-min '+args.sql_ky_offset_min+\
        ' --sql-ky-shift-step '+args.sql_ky_offset_step+\
        ' --sql-ky-shift-nbr '+args.sql_ky_offset_nbr+\
        ' --sql-completed '+args.sql_completed)
    
################################################################################
######      CALL DIAGONALIZATION ROUTINES SEQUENTIALLY      ####################

    if args.run_diag:

        ##  Obtain a set of incomplete SQL ids from the SQL table

        sqlDb = SQLite(DATA_FOLDER+args.FolderName()+args.SqlName()+'.sql')

        colNames,rows = sqlDb.GetData("Id",SQL_TABLE_NAME,"Completed=0")    

        sqlIdList = []

        for row in rows:

	        sqlIdList.append(row[0])

        #print(sqlIdList)

        ##  Get the number of parameter sets in the SQL table by reading 
        ##  the number of lines in the text file version

        nbrId = sum(1 for line in open(DATA_FOLDER+args.FolderName()+args.SqlName()+'.txt'))

        print ("\n========== Calling diagonalization program for each"\
        +" remaining case recorded in the SQLite database from id="+\
        str(args.sql_id_min)+" to id="+str(min(nbrId,int(args.sql_id_max)))+\
        " ==========")

        sectorListArg = ""

        if args.sector_list != "":
	        sectorListArg = " -s "+args.sector_list

        for sqlId in sqlIdList:

	        if sqlId in range (int(args.sql_id_min),min(nbrId,int(args.sql_id_max))):

		        arguments = args.CppArgs()+' --use-sql 1 --sql-id '+ \
		        str(sqlId)+' --out-path '+DATA_FOLDER+args.FolderName()+' --sql-name '+args.SqlName() + \
		        ' --sql-table-name '+SQL_TABLE_NAME+' --in-path '+DATA_FOLDER+args.FolderName()+\
		        ' --method 1 --eigenvalues-file 1 --eigenvectors-file 0 -b 1 -v 3 --calculate-occupations 1 '+\
		        ' --calculate-density-density 1 -e '+str(args.nbr_eigs) +\
		        ' --diagonalize 1 --most-probable-list 1 --nbr-most-probable 12 '+\
		        ' --calculate-participation-ratio 1 --calculate-translational-density-density 1 '+\
		        ' --calculate-rotational-density-density 0 '+\
		        ' --sort-vectors 0 --kx-shift '+args.offset_x+' --ky-shift '+args.offset_y+sectorListArg

		        RunDiag(MPI_CORES,arguments)

		        ##ClusterRunCommand(arguments,args.job_length)

################################################################################
######      READ IN DENSITY-DENSITY FUNCTION AND APPLY TO A SPECIFIED     ######
######      K-SPACE MAP                   								  ######

    if args.density_density_analysis:
	
        if args.calculate_pomeranchuk:
            
            ##	Obtain the value of a "Pomeranchuk parameter", which attempts
            ##	to measure the number of Fermi surfaces in the Brillouin zone
            
            ##	Generate the Pommerahchuck map, which is a map of +1 in a region
            ##	around the centre of the Brillouin zone

            pomer_map = []

            exclusionRange = 1.0
            midX = args.x/2.0 - float(args.offset_x)
            midY = args.y/2.0 - float(args.offset_y)

            for x in range(0,args.x):
	            for y in range(0,args.y):
		            if (x-midX)**2 + (y-midY)**2 > exclusionRange:
			            pomer_map.append(0)
		            else:
			            pomer_map.append(1)
        
        if args.calculate_r:
        
            ##  Obtain the value of the R-60 parameter, which measures the
            ##  value of an order parameter that's +1 in the K point segments
            ##  of the BZ and -1 in the K' segments. It's 0 if on the boundary
            
            r_map = []
            
            ##exclusionRange = 0.2
            
            ##  Record the locations of the sublattice corners
            
            sq3 = math.sqrt(3)
            
            ##a1 = [0.0,-0.5]
            ##a2 = [-sq3/2.0,-2.0]
            ##a3 = [sq3/2.0,-2.0]
            
            ##b1 = [0,-2.5]
            ##b2 = [-sq3/2.0,-1.0]
            ##b3 = [sq3/2.0,-1.0]
            
            '''
            for x in range(0,args.x):
                for y in range(0,args.y):
	            
                    x0 = float(x)/args.x*(-sq3/2.0)+float(y)/args.y*(sq3/2.0)+float(args.offset_x)
                    y0 = float(x)/args.x*(-1.5)+float(y)/args.y*(-1.5)+float(args.offset_y)
                    
                    ##  Check for prxomity to A-sublattice corners

                    if (x0-a1[0])**2+(y0-a1[1])**2 < exclusionRange or (x0-a2[0])**2+(y0-a2[1])**2 < exclusionRange or (x0-a3[0])**2+(y0-a3[1])**2 < exclusionRange:
                    
                        r_map.append(+1)
                    
                    ##  Check for prxomity to B-sublattice corners
                       
                    elif (x0-b1[0])**2+(y0-b1[1])**2 < exclusionRange or (x0-b2[0])**2+(y0-b2[1])**2 < exclusionRange or (x0-b3[0])**2+(y0-b3[1])**2 < exclusionRange:
                    
                        r_map.append(-1)
                        
                    else:

                        r_map.append(0)
            '''
            
            ##  Define three locus lines ky = ai kx + bi divinding the hexagonal
            ##  unit cell into 6 segments, each containing a K or K' point
            
            a = [sq3,0,-sq3]
            b = [-1.5,-1.5,-1.5]

            ##  Define BZ locus
            
            c = [-1/sq3,1/sq3,-1/sq3,1/sq3]
            d = [-0.5,-2.5,-2.5,-0.5]
            
            tol = 0.000000000000001
            
            for x in range(0,args.x):
                for y in range(0,args.y):
                
                    x0 = float(x)/args.x*(-sq3/2.0)+float(y)/args.y*(sq3/2.0)+float(args.offset_x)
                    y0 = float(x)/args.x*(-1.5)+float(y)/args.y*(-1.5)+float(args.offset_y)

                    ##  Generate a grid of equivalent x0 and y0 points from 
                    ##  up to 3 adjoining equivalent BZs, then select form that
                    ##  the point which lies in the main BZ
                    
                    equivalentPoints = []

                    flag = 0
                    
                    for i in [-2,-1,0,1,2]:
                        for j in [-2,-1,0,1,2]:

                            x1 = x0+i*sq3+j*sq3/2.0
                            y1 = y0-1.5*j
                            
                            ##  Check if it's in the hexagon

                            test = []
                            for k in range(0,len(c)):
                                test.append(c[k]*x1+d[k])
                                
                            if (abs(x1-sq3/2)<tol or x1 < sq3/2) and x1 > -sq3/2 and (abs(y1-test[0])<tol or y1 < test[0]) and (abs(y1-test[1])<tol or y1 > test[1]) and y1 > test[2] and y1 < test[3]:
                                x2 = x1
                                y2 = y1
                                flag = 1
                                break
                            
                        if flag:
                            break
                    
                    if not flag:
                        print("ERROR IN MAP CALCULATION")
                        quit()
                           
                    x0 = x2
                    y0 = y2
                    
                    ##print(x0,y0)
                    
                    ##  Check for location in each segment
                    
                    test = []
                    
                    for i in range(0,len(a)):    
                        test.append(a[i]*x0+b[i])

                    ##print(test)
                    
                    if abs((y0-test[0])) < tol or abs((y0-test[1])) < tol or abs((y0-test[2])) < tol:
                        r_map.append(0)
                        ##print("SEGMENT 0")
                    elif y0 < test[0] and y0 > test[1]:
                        r_map.append(1)
                        ##print("SEGMENT 1")
                    elif y0 < test[1] and y0 > test[2]:
                        r_map.append(-1)
                        ##print("SEGMENT 2")
                    elif y0 < test[2] and y0 < test[0]:
                        r_map.append(1)
                        ##print("SEGMENT 3")
                    elif y0 > test[0] and y0 < test[1]:
                        r_map.append(-1)
                        ##print("SEGMENT 4")
                    elif y0 > test[1] and y0 < test[2]:
                        r_map.append(1)
                        ##print("SEGMENT 5")
                    elif y0 > test[2] and y0 > test[0]:
                        r_map.append(-1)
                        ##print("SEGMENT 6")
                    else:
                        print("FAIL")
                        quit()
                    
            print(r_map)       
        
        ##  Read density-density functions from SQL database

        sqlDb = SQLite(DATA_FOLDER+args.FolderName()+args.SqlName()+'.sql')

        colNames,rows = sqlDb.GetData("Id,V0,Interaction,OutFileName",SQL_TABLE_NAME,"Completed=1 AND GotDensityDensity=1")

        sqlDb.PrintData(colNames,rows)

        mappedValueData    = []
        v0Data             = []
        interactionData    = []

        for row in rows:

            if args.calculate_susceptibility:
                ##	Obtain the susceptibility from the single-particle magnetization
                ##	map and the density-density function

                xData,yData,magnetization_map = ObtainMagnetizationMapData(row[0])
            
            for sectorX in range(0,args.x):
            
                for sectorY in range(0,args.y):
                    
                    #if [sectorX,sectorY] in [[0,2],[2,0],[2,2],[0,0]]:
                    
                    ##  Set sector so that the file name is correctly generated
                    ##  by the args.FolderName() function
                    args.sector = [sectorX,sectorY]

                    ##  Read eigenvalue data from the file name indicated

                    for eig in range(0,args.nbr_eigs):

                        fileName = DATA_FOLDER+args.FolderName()+row[3]+args.SectorName()+'_eig_'+str(eig)+'_density_density.dat'

                        ##	Filter density-density function						

                        value = 0
                        outFileName = ''

                        if args.calculate_susceptibility == 1:

                            value = ApplyDensityDensityToMap(fileName,magnetization_map)
                            
                            outFileName = DATA_FOLDER+args.FolderName()+row[3]+args.SectorName()+'_eig_'+str(eig)+'_susceptibility.dat'

                            ##  Output data to a file

                            try:
                                with open(outFileName,'w') as file:
                                    fout = open(outFileName,'w')
                                    fout.write(str(value))
                                    fout.close()
                            
                            except IOError as e:
                                print("ERROR: CANNOT WRITE FILE "+outFileName)

                        if args.calculate_pomeranchuk == 1:

                            value = ApplyDensityDensityToMap(fileName,pomer_map)
                            
                            outFileName = DATA_FOLDER+args.FolderName()+row[3]+args.SectorName()+'_eig_'+str(eig)+'_pomeranchuk_parameter.dat'

                            ##  Output data to a file

                            try:
                                with open(outFileName,'w') as file:
                                    fout = open(outFileName,'w')
                                    fout.write(str(value))
                                    fout.close()
                            
                            except IOError as e:
                                print("ERROR: CANNOT WRITE FILE "+outFileName)
                
                        if args.calculate_r == 1:

                            value = ApplyDensityDensityToMap(fileName,r_map)
                            
                            print(value)
                            
                            outFileName = DATA_FOLDER+args.FolderName()+row[3]+args.SectorName()+'_eig_'+str(eig)+'_r_parameter.dat'

                            ##  Output data to a file

                            try:
                                with open(outFileName,'w') as file:
                                    fout = open(outFileName,'w')
                                    fout.write(str(value))
                                    fout.close()
                            
                            except IOError as e:
                                print("ERROR: CANNOT WRITE FILE "+outFileName)
                                
################################################################################
######      READ IN DATA AND MAKE PLOTS OF LOWEST ENERGY STATES           ######
######      FOR A CROSS SECTION THROUGH PARAMETER SPACE                   ######

    if args.plot_energy:

        spectrumData = SpectrumDataArray()

        ##  Read from SQL database

        sqlDb = SQLite(DATA_FOLDER+args.FolderName()+args.SqlName()+'.sql')

        if args.cut_v0:    ##  keep g fixed
        
            colNames,rows = sqlDb.GetData("V0,Interaction,OutFileName",SQL_TABLE_NAME,"Completed=1 AND Interaction="+args.interaction+" AND OffsetX="+str(args.offset_x)+" AND OffsetY="+str(args.offset_y))  
        
        elif args.cut_g:   ##  keep v0 fixed
        
            colNames,rows = sqlDb.GetData("V0,Interaction,OutFileName",SQL_TABLE_NAME,"Completed=1 AND V0="+args.v0+" AND OffsetX="+str(args.offset_x)+" AND OffsetY="+str(args.offset_y))
        
        else:
            print("ERROR: need to select cut-g or cut-v0 option")
            quit()
            
        sqlDb.PrintData(colNames,rows)

        for row in rows:

            ##  Read eigenvalue data from the file name indicated
            ##  Loop over all available sectors
            
            for sectorX in range(0,args.x):
            
                for sectorY in range(0,args.y):
                    
                    #if [sectorX,sectorY] in [[0,2],[2,0],[2,2],[0,0]]:
                    
                    ##  Set sector so that the file name is correctly generated
                    ##  by the args.FolderName() function
                    args.sector = [sectorX,sectorY]

                    fileName = DATA_FOLDER+args.FolderName()+row[2]+args.SectorName()+'_eigensystem.dat'
                    
                    ##print("READING DATA FROM FILE: "+fileName)
                    
                    try:
                        with open(fileName) as file:
                            
                            fin = open(fileName)
                    
                            fin.readline()  ##  skip first comment line
                            
                            nbrEigenvalues = int(fin.readline())

                            fin.readline()  ##  skip first comment line
                            
                            ##  Read in data
                            
                            for i in range(0,min(nbrEigenvalues,args.nbr_eigs)):
                                
                                dataPoint = SpectrumDataPoint()
                                
                                dataPoint.energy      = float(fin.readline())
                                dataPoint.interaction = float(row[1])
                                dataPoint.v0          = float(row[0])/0.5       ##  Divide by recoil energy
                                dataPoint.sectors     = [[sectorX,sectorY]]
                                dataPoint.eigLabels   = [i]
                                
                                spectrumData.Append(dataPoint)
                                
                            fin.close()
                            
                    except IOError as e:
                        pass

        if not spectrumData.GetNbrPoints():

            print("ERROR: Data not found in table")
            
            quit()

        ##  Get maximum and minimum values of energy, interaction and V0
        
        maxEnergy,maxInteraction,maxV0 = spectrumData.GetMax()
        minEnergy,minInteraction,minV0 = spectrumData.GetMin()
        
        energyStep      = 0.02*(maxEnergy - minEnergy)
        interactionStep = 0.02*(maxInteraction - minInteraction)
        v0Step          = 0.02*(maxV0 - minV0)

        ##  Finally we can make some plots

        ######  2D cross section of phase diagram

        main_axes = plt.axes([0.09,0.14,0.75,0.75])	# ([left, bottom, width, height])

        ##  Define a function to print the set of sectors corresponding to 
        ##  a given degenerate energy level when an in interactive event occurs
        def onpick(event):
            index = event.ind
            
            print("CURRENT SELECTION:")
            
            for i in index:
            
                spectrumData.PrintOne(i)
            
            ##print('Selected energy level ',index,' includes the following sectors: ')
        
        ##  Plot main data set
        
        print ("\n========== Generating energy level plot ==========")
        
        dataEnergy      = spectrumData.GetEnergyList()
        dataInteraction = spectrumData.GetInteractionList()
        dataV0          = spectrumData.GetV0List()
        
        if args.cut_v0:    ##  keep g fixed
            
            main_axes.plot(dataV0,dataEnergy,'_',ms=20,picker=5)
            
            main_axes.set_title(r'Interaction = '+str(float(args.interaction)),fontsize=60)
            main_axes.set_xlabel(r'$\frac{V_0}{E_R}$',fontsize=60)
            main_axes.set_xlim(minV0-v0Step,maxV0+v0Step)

        elif args.cut_g:   ##  keep v0 fixed
            
            main_axes.plot(dataInteraction,dataEnergy,'_',ms=20,picker=5)
            
            main_axes.set_title(r'$\frac{V_0}{E_R} $= '+str(float(args.v0)/0.5),fontsize=60)
            main_axes.set_xlabel(r'Interaction',fontsize=60)        
            main_axes.set_xlim(minInteraction-interactionStep,maxInteraction+interactionStep)

        main_axes.set_ylabel(r'Energy',fontsize=60)
        main_axes.set_ylim(minEnergy-energyStep,maxEnergy+energyStep)

        for t in main_axes.get_yticklabels():
            t.set_fontsize(40)
                
        for t in main_axes.get_xticklabels():
            t.set_fontsize(40)

        ##  Resize plt.show to a maximized window
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        
        ##  Set up interactive aspects of the plot
        
        browser = PointBrowser()
        
        main_axes.get_figure().canvas.mpl_connect('pick_event', browser.onpick)        
        main_axes.get_figure().canvas.mpl_connect('key_press_event', browser.onpress) 
        
        plt.show()
        
        plt.close()
        
        ##plt.savefig("energy_plot"+FILE_EXT,bbox_inches='tight')
        
        ##plt.savefig(PLOTS_FOLDER+args.FolderName()+"energy_cross_section_g_"+str(int(float(args.interaction)*1000))+FILE_EXT,bbox_inches='tight')
        ##plt.savefig(PLOTS_FOLDER+args.FolderName()+"energy_cross_section_v0_"+str(int(float(args.v0)*1000))+FILE_EXT,bbox_inches='tight')

################################################################################
######      AN INTERACTIVE PLOT SHOWING THE EVOLUTION OF THE        ############
######      PROBABILITY DISTRIBUTION OF A GIVEN SECTOR FOR A CUT    ############
######      THROUGH PARAMETER SPACE                                 ############

    if args.plot_probability:

        ##  Read from SQL database

        sqlDb = SQLite(DATA_FOLDER+args.FolderName()+args.SqlName()+'.sql')
        
        if args.cut_v0:    ##  keep g fixed
        
            colNames,rows = sqlDb.GetData("Id,V0,Interaction,OutFileName,OffsetX,OffsetY",SQL_TABLE_NAME,"Completed=1 AND GotOccupations=1 AND Interaction="+args.interaction)
        
        elif args.cut_g:   ##  keep v0 fixed
            
            colNames,rows = sqlDb.GetData("Id,V0,Interaction,OutFileName,OffsetX,OffsetY",SQL_TABLE_NAME,"Completed=1 AND GotOccupations=1 AND V0="+args.v0)

        else:
            print("ERROR: need to select cut-g or cut-v0 option")
            quit()

        sqlDb.PrintData(colNames,rows)

        probDataFull    = []
        bandFullData    = []
        v0Data          = []
        interactionData = []
        dimX = 0
        dimY = 0

        ##  Obtain plot data

        for row in rows:

            sqlId = row[0]
            v0Data.append(row[1])
            interactionData.append(row[2])
            fileName=DATA_FOLDER+args.FolderName()+row[3]+args.SectorName()+'_eig_0_occupations.dat'
            args.offset_x = str(row[4])
            args.offset_y = str(row[5])
            
            ##  Read in all of the probability distribution data

            xData = []
            yData = []
            probData = []
            bandData = []
            
            fin  = open(fileName)
            
            dimX = int(fin.readline()) 
            dimY = int(fin.readline())
            
            for line in fin:
                line = line.split()
                xData.append(float(line[0]))
                yData.append(float(line[1]))
                probData.append(float(line[2]))
                
            fin.close()

            ##  Convert to np array
            probData = np.array(probData)

            ##  Obtain the single particle band structure data

            xData,yData,bandData = ObtainSingleParticleEnergyData(sqlId)
            
            print("xData",xData)
            print("yData",yData)
            print("bandData",bandData)
            
            xData2,yData2,bandData2 = ExpandCellData(xData,yData,bandData)

            xData2,yData2,probData2 = ExpandCellData(xData,yData,probData)
            
            print("xData2",xData2)
            print("yData2",yData2)
            print("bandData2",bandData2)

            ##  Append to overall list of plot data

            probDataFull.append(probData2)
            bandFullData.append(bandData2)

        ##  Calculate triangularization
            
        triang = tri.Triangulation(xData2,yData2)
        
        triang_masked = tri.Triangulation(xData2,yData2)
        
        xmid = xData2[triang_masked.triangles].mean(axis=1)
        ymid = yData2[triang_masked.triangles].mean(axis=1)
        mask = np.where(xmid*xmid + (ymid+1.5)*(ymid+1.5) > 2.2, 1, 0)
        triang_masked.set_mask(mask)
         
        ######      Set up the interactive plot        #####################
            
        print ("\n========== Generating probability density plot ==========")

        ##  SET PLOT PARAMETERS HERE

        font_size = 30
        contour_label_size = 20
        x_tick_pad = 15
        start_pos = 0
        pdf_plot = 0
        fig_name = "probability_density_"+args.SectorName()+"_id_"+str(start_pos)+".pdf"

        contourLevels=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.8,1.0]
        contourLabels=[r'\textbf{0.0}',r'\textbf{0.05}',r'\textbf{0.1}',r'\textbf{0.15}',r'\textbf{0.2}',r'\textbf{0.25}',r'\textbf{0.3}',r'\textbf{0.4}' \
        ,r'\textbf{0.5}',r'\textbf{0.6}',r'\textbf{0.8}',r'\textbf{1.0}']
        
        if pdf_plot==1:
        
            font_size = 20
            contour_label_size = 10
            x_tick_pad = 10
        
        main_axes = plt.axes([0.0, 0.15, 1.0, 0.7]) 	# ([left, bottom, width, height])
        main_axes.set_aspect('equal')
        
        PlotUnitCellBackground(main_axes)
        
        fig = main_axes.get_figure()

        ##  Plot single particle band structure contours
        p1 = main_axes.tricontour(triang_masked,bandFullData[start_pos],cmap=cm.gray)

        ##  Write level labels on the contours
        p2 = main_axes.clabel(p1, inline=1, fontsize=contour_label_size)
        
        ##  Plot the first set of occupation probability data
        p3 = main_axes.tricontourf(triang,probDataFull[start_pos],cmap=COLOUR_MAP,levels=contourLevels)
        
        ##  Plot the occupation probability colour bar
        cb  = fig.colorbar(p3, shrink=0.8,ticks=contourLevels)
        cb.set_label(r"\textbf{Occupation}"+"\n"+r"\textbf{Probability}",fontsize=font_size-8,rotation=0)

        cb.ax.yaxis.set_label_coords(0.5,1.17)
        
        cb.set_ticklabels(contourLabels)
        
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(font_size-8)

        if pdf_plot==0:
        
            main_axes.set_title(r'$V_0 = '+str(v0Data[start_pos])+'$  $g_{2D} = '+str(interactionData[start_pos])+'$',fontsize=font_size)
        
            ##  Add a discrete slider to control the model parameter values

            axis_slider = fig.add_axes([0.48, 0.9, 0.3, 0.03])# ([left, bottom, width, height])
            
            slider = DiscreteSlider(axis_slider, '',0,len(interactionData),allowed_vals=range(0,len(interactionData)), valinit=start_pos)

            ##  Resize plt.show to a maximized window
            mng = plt.get_current_fig_manager()
            mng.resize(*mng.window.maxsize())

        ######      Define a function to update the figure 
        ######      for the ith set of data 
        def update_fig(val):

            global p1,p2,p3

            i = int(slider.val)

            ##  Update the title
            main_axes.set_title(r'$V_0 = '+str(v0Data[i])+'$  $g_{2D} = '+str(interactionData[i])+'$',fontsize=font_size)
            
            ##  Remove previous set of contour data
            for coll in p3.collections:
                main_axes.collections.remove(coll)
            
            ##  Update the probability contour data  
            p3 = main_axes.tricontourf(triang,probDataFull[i],cmap=COLOUR_MAP,levels=contourLevels)
            
            if args.cut_v0:
            
                ##  Also update the single particle band structure contours
                ##  in the case where these are not fixed
                
                ##  Remove previous set of contours
                for coll in p1.collections:
                    main_axes.collections.remove(coll)
                
                ##  Plot new contours
                p1 = main_axes.tricontour(triang_masked,bandFullData[i],cmap=cm.gray)
                
                ##  Remove previous set of contour labels
                for text in p2:
                    text.remove()
                
                ##  Plot new contour labels
                p2 = main_axes.clabel(p1, inline=1, fontsize=contour_label_size)
         
            fig.canvas.draw_idle()

        if pdf_plot==1:

            fig.savefig(fig_name,bbox_inches='tight')
            
            print ("\nDONE! find plot at "+fig_name)
        else:
        
            slider.on_changed(update_fig)

            plt.show()

################################################################################
######      PLOT A GIVEN OBSERVABLE OF THE GROUND STATE AS A        ############
######      FUNCTION OF SYSTEM PARAMETERS                           ############

    if args.plot_combined_observable:

        ##  Read in the energy level data and determine which sector has the
        ##  lowest energy at each point in parameter space. Then obtain the 
        ##  observable data for that sector
        
        if args.observe_susceptibility:
        
            gotFlag     = "GotDensityDensity"
            fileLabel   = "susceptibility"
            plotLabel   = r"$\rho_{\mbox{\small S}}$"
            norm        = args.nbr*args.nbr
        
        elif args.observe_r:
        
            gotFlag     = "GotDensityDensity"
            fileLabel   = "r_parameter"
            plotLabel   = r"$\tilde\rho_{\mbox{\small S}}$"
            norm        = args.nbr*args.nbr
            
        elif args.observe_pomer:
        
            gotFlag     = "GotDensityDensity"
            fileLabel   = "pomeranchuk_parameter"
            plotLabel   = "Centralness"
            norm        = args.nbr*args.nbr
        
        elif args.observe_participation:
        
            gotFlag     = "GotParticipationRatio"
            fileLabel   = "participation_ratio"
            plotLabel   = "$P$"
            norm        = scipy.special.binom(int(args.x)*int(args.y),int(args.nbr))/(int(args.x)*int(args.y))
            
        elif args.observe_rotational:
        
            gotFlag   = "Completed"
            fileLabel = "rotational_density_density"
            plotLabel = r"$\delta\rho_{\mbox{\small R}} (3)$"
            norm      = 1
        else:
        
            print("ERROR in plot_combined_observable: no observable selected")
            quit()
            
        observableData  = []
        degenData       = []
        v0Data          = []
        interactionData = []
        offsetData		= []
        labelList       = []

        for system in args.systemList:
            
            args.x   = system[0]
            args.y   = system[1]
            args.nbr = system[2]

            observableDataInner  = []
            degenDataInner       = []
            v0DataInner          = []
            interactionDataInner = []
            offsetDataInner		 = []
            labelListInner       = []

            sqlDb = SQLite(DATA_FOLDER+args.FolderName()+args.SqlName()+'.sql')

            if args.cut_v0:    ##  keep g fixed
                colNames,rows = sqlDb.GetData("V0,Interaction,OffsetX,OffsetY,OutFileName",SQL_TABLE_NAME,"Completed=1 AND "+gotFlag+"=1 AND Interaction="+args.interaction)
            elif args.cut_g:   ##  keep v0 fixed
                colNames,rows = sqlDb.GetData("V0,Interaction,OffsetX,OffsetY,OutFileName",SQL_TABLE_NAME,"Completed=1 AND "+gotFlag+"=1 AND V0="+args.v0)
            else:
                colNames,rows = sqlDb.GetData("V0,Interaction,OffsetX,OffsetY,OutFileName",SQL_TABLE_NAME,"Completed=1 AND "+gotFlag+"=1 AND V0 <> 0")        
                ##AND Interaction <> 0.01

            sqlDb.PrintData(colNames,rows)

            for row in rows:

                ##  Read eigenvalue data from the file name indicated
                ##  Loop over all available sectors

                spectrumData = SpectrumDataArray()

                for sectorX in range(0,args.x):
                
                    for sectorY in range(0,args.y):
                        
                        ##  Set sector so that the file name is correctly generated
                        ##  by the args.FolderName() function
                        args.sector = [sectorX,sectorY]

                        fileName = DATA_FOLDER+args.FolderName()+str(row[4])+args.SectorName()+'_eigensystem.dat'
                        
                        print("READING DATA FROM FILE: "+fileName)
                        
                        try:
                            with open(fileName) as file:
                                
                                fin = open(fileName)
                        
                                fin.readline()  ##  skip first comment line
                                
                                nbrEigenvalues = int(fin.readline())

                                fin.readline()  ##  skip first comment line
                                
                                ##  Read in data
                                
                                for i in range(0,min(nbrEigenvalues,args.nbr_eigs)):
                                    
                                    dataPoint = SpectrumDataPoint()
                                    
                                    dataPoint.energy      = float(fin.readline())
                                    dataPoint.interaction = float(row[1])##-0.01
                                    ##  Divide by recoil energy
                                    dataPoint.v0          = float(row[0])/0.5##-0.02
                                    dataPoint.offset	  = [[float(row[2]),float(row[3])]]
                                    dataPoint.sectors     = [[sectorX,sectorY]]
                                    dataPoint.eigLabels   = [i]

                                    spectrumData.Append(dataPoint)
                                    
                                fin.close()    
                                
                        except IOError as e:
                            print("ERROR READING FILE "+fileName)
                        
                            pass
                
                if not spectrumData.GetNbrPoints():

                    print("ERROR: Data not found")
                    
                    ##quit()
                
                else:
                
                    ##  Sort the spectrum data to obtain the minimum energy sector
                    ##  and eigenvalue label
                    
                    spectrumData.Sort()

                    ##  Use the first value in the spectrumData array to set the
                    ##  lowest energy sector and eigenvalue label (allowing
                    ##  for a number of degenerate states)
                    '''
                    if spectrumData.GetInteraction(0) == 0.0:
                    
                        ##  In the case of zero interaction there are multiple
                        ##  degenerate states, but for simplicity we just plot
                        ##  the 0,0 sector with the lowest eigenvalue label
                    
                        sectors = [[0,0]]
                        eigLabels = [0]
                        
                    else:
                    '''
                    sectors   = spectrumData.GetSectors(0)
                    eigLabels = spectrumData.GetEigLabels(0)

                    index = 0;

                    observableMean = 0.0

                    for s in sectors:

                        args.sector = s

                        eigLabel = eigLabels[index]

                        index += 1

                        ##  Read in the corresponding observable data from a file
                        ##  along with the associated model parameter values

                        fileName = DATA_FOLDER+args.FolderName()+row[4]+args.SectorName()+'_eig_'+str(eigLabel)+'_'+fileLabel+'.dat'
                        
                        print(fileName)
                        
                        try:
                            with open(fileName) as file:
                                
                                fin = open(fileName)
                                
                                if args.observe_rotational:
                                    
                                    ##  Skip the 0 peak
                                    fin.readline()
                                    
                                    line = fin.readline().split()
                                    peakAt1 = float(line[1])
                                    line = fin.readline().split()
                                    peakAt2 = float(line[1])
                                    line = fin.readline().split()
                                    peakAt3 = float(line[1])
                                    line = fin.readline().split()
                                    peakAt4 = float(line[1])
                                    line = fin.readline().split()
                                    peakAt5 = float(line[1])
                                    
                                    mean = (peakAt1+peakAt2+peakAt3+peakAt4+peakAt5)/5.0
                                    
                                    observableMean += (peakAt3-mean)/mean

                                else:
                                
                                    observableMean += float(fin.readline())
                                
                                fin.close()
                                
                        except IOError as e:
                            
                            print("ERROR READING FILE "+fileName)
                            quit()
                            
                            pass
                    
                    if observableMean != float('Inf'):
                        
                        ##  KLUDGE to avoid colour scale for uninteresting negative values
                        if observableMean<0:
                            observableMean=-0.01
                        
                        observableDataInner.append(observableMean/len(sectors)/(norm))
                        
                        ##/(args.nbr*args.nbr))
                        degenDataInner.append(len(sectors))
                        v0DataInner.append(dataPoint.v0)
                        interactionDataInner.append(dataPoint.interaction)
                        offsetDataInner.append(dataPoint.offset[0][0])  ##  Only get the x offset
                        labelListInner.append(r"$N="+str(args.nbr)+r"$, $\Delta_x="+str(dataPoint.offset[0][0])+"$")

            ##  Convert to np arrays
        
            observableDataInner = np.array(observableDataInner)
            degenDataInner = np.array(degenDataInner)
            v0DataInner = np.array(v0DataInner)
            interactionDataInner = np.array(interactionDataInner)
            offsetDataInner = np.array(offsetDataInner)
            labelListInner  = np.array(labelListInner)

            observableData.append(observableDataInner)
            degenData.append(degenDataInner)
            v0Data.append(v0DataInner)
            interactionData.append(interactionDataInner)
            offsetData.append(offsetDataInner)
            labelList.append(labelListInner)

        ##  Reshape the arrays to label by different offset values as well as system size
        
        uniqueOffsets = {}
    
        for offsetList in offsetData:
            for offset in offsetList:
                uniqueOffsets[offset] = 1
        
        uniqueOffsets = sorted(uniqueOffsets.keys()[:])

        print(uniqueOffsets)

        newObservableData   = []
        newDegenData        = []
        newV0Data           = []
        newInteractionData  = []
        newOffsetData       = []
        newLabelList        = []

        for testOffset in uniqueOffsets:
            
            newObservableDataInner  = []
            newDegenDataInner       = []
            newV0DataInner          = []
            newInteractionDataInner = []
            newOffsetDataInner      = []
            newLabelListInner       = []
            
            counter = 0

            for offsetList in offsetData:
                
                newObservableDataInner2     = []
                newDegenDataInner2          = []
                newV0DataInner2             = []
                newInteractionDataInner2    = []
                newOffsetDataInner2         = []
                newLabelListInner2          = []

                counter2 = 0

                for offset in offsetList:
                    if abs(offset-testOffset)<=0.000000000001:
                        newObservableDataInner2.append(observableData[counter][counter2])
                        newDegenDataInner2.append(degenData[counter][counter2])
                        newV0DataInner2.append(v0Data[counter][counter2])
                        newInteractionDataInner2.append(interactionData[counter][counter2])
                        newLabelListInner2.append(labelList[counter][counter2])
                    counter2 += 1

                ##  Convert to np arrays
            
                newObservableDataInner2 = np.array(newObservableDataInner2)
                newDegenDataInner2 = np.array(newDegenDataInner2)
                newV0DataInner2 = np.array(newV0DataInner2)
                newInteractionDataInner2 = np.array(newInteractionDataInner2)
                newOffsetDataInner2 = np.array(newOffsetDataInner2)
                newLabelListInner2 = np.array(newLabelListInner2[0])

                newObservableDataInner.append(newObservableDataInner2)
                newDegenDataInner.append(newDegenDataInner2)
                newV0DataInner.append(newV0DataInner2)
                newInteractionDataInner.append(newInteractionDataInner2)
                newOffsetDataInner.append(newOffsetDataInner2)
                newLabelListInner.append(newLabelListInner2)

                print(newObservableDataInner)

                counter += 1

            newObservableData.append(newObservableDataInner)
            newDegenData.append(newDegenDataInner)
            newV0Data.append(newV0DataInner)
            newInteractionData.append(newInteractionDataInner)
            newOffsetData.append(newOffsetDataInner)
            newLabelList.append(newLabelListInner)

        observableData  = newObservableData[:]
        degenData       = newDegenData[:]
        v0Data          = newV0Data[:]
        interactionData = newInteractionData[:]
        labelList       = newLabelList[:]
        
        print(labelList)

        ##print(observableData)
        
        ##print(degenData)

        print(v0Data[0][0])

        print(interactionData[0][0])

        ##  Make plot of the observable

        font_size = 28

        main_axes = plt.axes([0.1,0.1,0.5,0.5])	# ([left, bottom, width, height])
        
        main_axes.set_ylabel(r'$\chi / N^2$')

        if args.cut_v0 == 1:    ##  keep g fixed
        
            counter = 0
            counterTot = 0

            for data in v0Data:

                counter2 = 0                

                for val in data:

                    main_axes.scatter(v0Data[counter][counter2],susceptibilityData[counter][counter2],marker=MARKER_LIST[counterTot],color=COLOUR_LIST[counterTot])
                    main_axes.set_title(r'Interaction = '+str(float(interactionData[0][0][0])),fontsize=font_size)
                    main_axes.set_xlabel(r'$\frac{V_0}{E_R}$',fontsize=font_size)
                    ##main_axes.set_xlim(minV0-v0Step,maxV0+v0Step)
                    
                    counter2 += 1
                    counterTot += 1

                counter += 1
            
            plt.grid()

            fig_name = fileLabel+"_combined.pdf"

        elif args.cut_g == 1:   ##  keep v0 fixed
    
            counter = 0
            counterTot = 0

            for data in interactionData:
                
                print(data)

                counter2 = 0                

                for val in data:
                    
                    print(interactionData[counter][counter2])
                    print(observableData[counter][counter2])
                    print(labelList[counter][counter2])

                    main_axes.scatter(interactionData[counter][counter2],observableData[counter][counter2],marker=MARKER_LIST[counterTot],color=COLOUR_LIST[counterTot],cmap=COLOUR_MAP,label=labelList[counter][counter2])
                    
                    main_axes.set_title(r'$\frac{V_0}{E_R} $= '+str(float(v0Data[0][0][0])),fontsize=font_size)
                    main_axes.set_xlabel(r'$\frac{\tilde{g}}{2 \pi}$',fontsize=font_size)
                    
                    main_axes.legend(loc='best',scatterpoints=1,ncol=3,columnspacing=0.1,fontsize=12)
                    main_axes.set_ylim(0.0,0.2)
                    main_axes.set_xlim(0.0,1.0)

                    counter2 += 1
                    counterTot += 1

                counter += 1
            
            plt.grid()

            fig_name = fileLabel+"_combined.pdf"

        else:
            
            ##contourLevels=[np.min(observableData[0][0])]        
            ##gradation = (np.max(observableData[0][0]) - np.min(observableData))/20

            contourLevels=[0.0]
            gradation = 0.01

            while contourLevels[-1]+gradation < np.max(observableData[0][0]):
                contourLevels.append(contourLevels[-1]+gradation)

            ##contourLevels.append(contourLevels[-1]+gradation)

            if args.formatting==1:
                contourLevels = np.arange(0,0.36,0.03)
                ##contourLevels = np.arange(0,0.12,0.01)
                
            if args.formatting==2:
                contourLevels = np.arange(0,0.39,0.03)[0:13]

                contourLevels = np.concatenate((np.arange(np.min(observableData),0,np.abs(np.min(observableData))),contourLevels))

            print(contourLevels)
            
            ##  Only plot the first set of data in the list

            triang = tri.Triangulation(interactionData[0][0],v0Data[0][0])

            p = main_axes.tricontourf(triang,observableData[0][0],cmap=COLOUR_MAP,levels=contourLevels)

            

            main_axes.set_ylabel(r'$V_0/E_R$',fontsize=font_size,rotation=0)
            main_axes.yaxis.set_label_coords(0,1.03)
            main_axes.set_xlabel(r'$\tilde{g}/2\pi$',fontsize=font_size)
            main_axes.xaxis.set_label_coords(1.2,0.01)

            if args.formatting==1:
                main_axes.set_ylim([0.02,4])
                main_axes.set_yticks([0.02,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0])
                main_axes.set_yticklabels([r'$\textbf{0.02}$',r'$\textbf{0.5}$',r'$\textbf{1.0}$',r'$\textbf{1.5}$',r'$\textbf{2.0}$',r'$\textbf{2.5}$',r'$\textbf{3.0}$',r'$\textbf{3.5}$',r'$\textbf{4.0}$'],fontsize=font_size-10)
            else:
                main_axes.set_ylim([0.02,3])
                main_axes.set_yticks([0.02,0.5,1.0,1.5,2.0,2.5,3.0])
                main_axes.set_yticklabels([r'$\textbf{0.02}$',r'$\textbf{0.5}$',r'$\textbf{1.0}$',r'$\textbf{1.5}$',r'$\textbf{2.0}$',r'$\textbf{2.5}$',r'$\textbf{3.0}$'],fontsize=font_size-10)
            
            if args.formatting==3 or args.formatting==2:
                main_axes.set_xlim([0.01,1])
                main_axes.set_xticks([0.01,0.2,0.4,0.6,0.8,1.0])
                main_axes.set_xticklabels([r'$\textbf{0.01}$',r'$\textbf{0.2}$',r'$\textbf{0.4}$',r'$\textbf{0.6}$',r'$\textbf{0.8}$',r'$\textbf{1}$'],fontsize=font_size-10)
            else:
                main_axes.set_xlim([0,1])
                main_axes.set_xticks([0.0,0.2,0.4,0.6,0.8,1.0])
                main_axes.set_xticklabels([r'$\textbf{0}$',r'$\textbf{0.2}$',r'$\textbf{0.4}$',r'$\textbf{0.6}$',r'$\textbf{0.8}$',r'$\textbf{1}$'],fontsize=font_size-10)
            ##for i, txt in enumerate(degenData[0][0]):
            ##    main_axes.annotate(txt, (interactionData[0][0][i],v0Data[0][0][i]))

            main_axes.tick_params(direction='in',pad=9)

            fig = main_axes.get_figure()
            
            ##  Plot the occupation probability colour bar    
            cb  = fig.colorbar(p, shrink=0.8,ticks=contourLevels)
            cb.set_label(r"\textbf{"+plotLabel+"}",fontsize=font_size,rotation=0)

            cb.ax.yaxis.set_label_coords(1.4,1.17)

            cb.ax.set_yticks(contourLevels)
            
            if args.formatting==1:
                cb.ax.set_yticklabels([r'$\textbf{0.00}$',r'$\textbf{0.03}$',r'$\textbf{0.06}$',r'$\textbf{0.09}$',r'$\textbf{0.12}$',r'$\textbf{0.15}$',r'$\textbf{0.18}$',r'$\textbf{0.21}$',r'$\textbf{0.24}$',r'$\textbf{0.27}$',r'$\textbf{0.30}$',r'$\textbf{0.33}$',r'$\textbf{0.36}$'][0:len(contourLevels)],fontsize=font_size-10)
            elif args.formatting==2:
                cb.ax.set_yticklabels([r'$-\textbf{0.5}$',r'$\textbf{0.00}$',r'$\textbf{0.03}$',r'$\textbf{0.06}$',r'$\textbf{0.09}$',r'$\textbf{0.12}$',r'$\textbf{0.15}$',r'$\textbf{0.18}$',r'$\textbf{0.21}$',r'$\textbf{0.24}$',r'$\textbf{0.27}$',r'$\textbf{0.30}$',r'$\textbf{0.33}$',r'$\textbf{0.36}$'][0:len(contourLevels)],fontsize=font_size-10)
            else:
                cb.ax.set_yticklabels([r'$\textbf{0.00}$',r'$\textbf{0.01}$',r'$\textbf{0.02}$',r'$\textbf{0.03}$',r'$\textbf{0.04}$',r'$\textbf{0.05}$',r'$\textbf{0.06}$',r'$\textbf{0.07}$',r'$\textbf{0.08}$',r'$\textbf{0.09}$',r'$\textbf{0.10}$',r'$\textbf{0.11}$',r'$\textbf{0.12}$',r'$\textbf{0.13}$',r'$\textbf{0.14}$',r'$\textbf{0.15}$',r'$\textbf{0.16}$',r'$\textbf{0.17}$',r'$\textbf{0.18}$',r'$\textbf{0.19}$',r'$\textbf{0.20}$'][0:len(contourLevels)],fontsize=font_size-10)

            plt.grid()

            ##plt.show()

            fig_name = args.GetLabel()+fileLabel+"_combined.pdf"

        plt.savefig(fig_name,bbox_inches='tight')
        
        print ("\nDONE! find plot at "+fig_name)

################################################################################
######      PLOT THE SINGLE PARTICLE MAGNETIZATION MAP      ####################

    if args.plot_magnetization_map:

        ##  Get data from existing files, read in a plot a map of the single 
        ##  particle  magnetization

        xData,yData,magData = ObtainMagnetizationMapData(0)

        ##  Generate periodic pattern by shifting kx,ky

        xData2,yData2,magData2 = ExpandCellData(xData,yData,magData)

        ##  Set up plot as a 2d colour map 
        
        print("\n========== Generating magnetization map plot ==========")

        ##  SET PLOT PARAMETERS HERE

        font_size = 20
        x_tick_pad = 10
        start_pos = 0
        fig_name = "magnetization_map_kx_"+str(args.x)+'_ky_'+str(args.y)+".pdf"

        contourLevels = np.arange(-1.1,1.1,0.05)
        
        ##[np.min(magData)]
        
        ##gradation = (np.max(magData) - np.min(magData))/10
        
        ##print(gradation)
        
        ##for i in range(1,11):
        ##    contourLevels.append(contourLevels[i-1]+gradation)

        main_axes = plt.axes([0.0, 0.15, 1.0, 0.7]) 	# ([left, bottom, width, height])
        main_axes.set_aspect('equal')
        
        PlotUnitCellBackground(main_axes)
        
        fig = main_axes.get_figure()

        triang = tri.Triangulation(xData2,yData2)

        print(triang)

        p = main_axes.tricontourf(triang,magData2,cmap='coolwarm',levels=contourLevels)
        
        sq3 = math.sqrt(3)

        ##  Add up and down arrows on the hexagon corners

        main_axes.annotate('', xy=(0, -0.6), xytext=(0, -0.9),fontsize=font_size,
            arrowprops=dict(facecolor='black', shrink=0.05),)
        
        main_axes.annotate('', xy=(0, -2.4), xytext=(0, -2.1),fontsize=font_size,
            arrowprops=dict(facecolor='black', shrink=0.05),)
            
        main_axes.annotate('', xy=(0.6, -1.25), xytext=(0.6, -0.95),fontsize=font_size,
            arrowprops=dict(facecolor='black', shrink=0.05),)
        
        main_axes.annotate('', xy=(-0.6, -1.25), xytext=(-0.6, -0.95),fontsize=font_size,
            arrowprops=dict(facecolor='black', shrink=0.05),)
        
        main_axes.annotate('', xy=(0.6, -1.65), xytext=(0.6, -1.95),fontsize=font_size,
            arrowprops=dict(facecolor='black', shrink=0.05),)
        
        main_axes.annotate('', xy=(-0.6, -1.65), xytext=(-0.6, -1.95),fontsize=font_size,
            arrowprops=dict(facecolor='black', shrink=0.05),)
        
        ##  Mask off edges of main hexagon

        ##main_axes.axvspan(-3,-sq3/2-0.002,color='white')
        ##main_axes.axvspan(sq3/2+0.002,3,color='white')
        
        main_axes.set_ylim(-2.5-0.006,-0.5+0.006)
        main_axes.set_xlim(-sq3/2-0.006,sq3/2+0.006)
        
        lower_mask = []
        upper_mask = []
        x0 = []
        
        for x in range(-(args.x+3),args.x+3):
        
            x0.append(float(x)/args.x*(-sq3/2.0))
        
            if x0[-1] > 0:
                lower_mask.append(x0[-1]/sq3 -2.5)
                upper_mask.append(-x0[-1]/sq3 -0.5)
            else:
                lower_mask.append(-x0[-1]/sq3 -2.5)
                upper_mask.append(x0[-1]/sq3 -0.5)
      
        plt.fill_between(x0,lower_mask,-3.5,color='white')
        plt.fill_between(x0,upper_mask,3,color='white')

        main_axes.get_xaxis().set_visible(False)
        main_axes.get_yaxis().set_visible(False)
        ##main_axes.axis('off')
        main_axes.spines['bottom'].set_color('white')
        main_axes.spines['top'].set_color('white')
        main_axes.spines['left'].set_color('white')
        main_axes.spines['right'].set_color('white')
   
        print("See figure at "+fig_name)

        plt.savefig(fig_name,bbox_inches='tight')
        
        ##print ("\nDONE! find plot at "+fig_name)

################################################################################
######      ANALYSE THE SINGLE PARTICLE ENERGY LEVELS IN THE LOWEST BAND    ####

    if args.analyse_single_particle_levels:

        ##  Get single particle energy and magnetization map data from existing files

        xData,yData,magnetizationData = ObtainMagnetizationMapData(0)
        xData,yData,energyData = ObtainSingleParticleEnergyData(0)

        xData2,yData2,energyData2 = ExpandCellData(xData,yData,energyData)
        
        ##  Make an interactive plot showing the different configurations and their
        ##  associated energy and magnetizations

        print("\n========== Generating single particle data plot ==========")

        ##  SET PLOT PARAMETERS HERE

        font_size = 50
        x_tick_pad = 15
        
        contourLevels=[np.min(energyData)]        
        gradation = (np.max(energyData) - np.min(energyData))/10
        
        for i in range(1,11):
            contourLevels.append(contourLevels[i-1]+gradation)

        ##  Make background plot

        main_axes = plt.axes([0.05, 0.1, 0.9, 0.7]) 	# ([left, bottom, width, height])
        main_axes.set_aspect('equal')
        PlotUnitCellBackground(main_axes)
        PlotSimulationGrid(main_axes)

        ##main_axes.scatter(xData2,yData2,c=energyData2,s=200,cmap=COLOUR_MAP)

        ##  Plot energy level contours
        
        triang = tri.Triangulation(xData2,yData2)

        p = main_axes.tricontour(triang,energyData2,cmap=COLOUR_MAP,levels=contourLevels)

        ##  Plot the occupation probability colour bar  
        
        fig = main_axes.get_figure()
          
        cb  = fig.colorbar(p, shrink=0.8,ticks=contourLevels)
        cb.set_label(r"\textbf{Single-particle}"+'\n'+r"\textbf{Energy}",fontsize=font_size-8,rotation=0)

        cb.ax.yaxis.set_label_coords(1.1,1.2)

        cb.ax.get_children()[4].set_linewidths(45.0)

        for t in cb.ax.get_yticklabels():
            t.set_fontsize(font_size-8)

        currConfiguration = ConfigurationAmplitudeData(args.x,args.y)

        xConfigurationData,yConfigurationData = currConfiguration.GetConfiguration()

        p1 = main_axes.scatter(xConfigurationData,yConfigurationData,marker="o",color=COLOUR_LIST[0],s=250)

        ##  Set the plot title
        kx,ky = currConfiguration.GetSector()
            
        main_axes.set_title(r'Total non-interacting energy '+str(currConfiguration.GetTotal(energyData))+'\n'+r'Total magnetization '+str(currConfiguration.GetTotal(magnetizationData))+'\n'+r'Sector '+str(kx)+" "+str(ky),fontsize=font_size-10)

        ##  Define an update function to remove an occupied state
        def update(point):
            
            global p1
            
            ##  Update the current occupation
            
            currConfiguration.FlipState(point.xdata,point.ydata)
            
            ##  Update the plot

            p1.remove()

            xConfigurationData,yConfigurationData = currConfiguration.GetConfiguration()

            p1 = main_axes.scatter(xConfigurationData,yConfigurationData,marker="o",color=COLOUR_LIST[0],s=250)
         
            ##  Update the plot title
            
            kx,ky = currConfiguration.GetSector()
            
            main_axes.set_title(r'Total non-interacting energy '+str(currConfiguration.GetTotal(energyData))+'\n'+r'Total magnetization '+str(currConfiguration.GetTotal(magnetizationData))+'\n'+r'Sector '+str(kx)+" "+str(ky),fontsize=font_size-10)
         
            fig.canvas.draw_idle()

        fig.canvas.mpl_connect('button_press_event', update)
        
        ##  Resize plt.show to a maximized window
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        
        plt.show()
        
        print("\n========== Generating Lowest possible energy values ==========")
        
        ##  Obtain a list of the combined single particle energies associated
        ##  with the lowest energy configurations
        
        ##  Convert data arrays to a list of objects, to allow for sorting
        
        data = []
        
        import itertools

        for x,y,m,e in itertools.izip(xData,yData,magnetizationData,energyData):

            currData = SingleParticleData()
            
            currData.kx = x
            currData.ky = y
            currData.magnetization = m
            currData.energy = e
            
            data.append(currData)
        
        ##  First, sort the single particle data into ascending order of energy value
        
        data.sort(key=lambda x: x.energy)

        energy = []
        magnetization = []
        
        for x in data:
            energy.append(x.energy)
            magnetization.append(x.magnetization)

        ##  Construct a list of the lowest lying excitations using the ordered
        ##  binary Hamming numbers

        currConfiguration.fockState = FirstBinaryHammingNumber(args.nbr)

        i=1
        imax = min(100,Binomial(args.x*args.y,args.nbr))
        
        data2 = []
        
        while i <= imax:

            currData = SingleParticleData()
            
            currData.energy = currConfiguration.GetTotal(energy)
            currData.magnetization = currConfiguration.GetTotal(magnetization)
            
            currConfiguration.fockState = NextHammingNumber64(currConfiguration.fockState)

            data2.append(currData)
        
            i+=1

        data2.sort(key=lambda x: x.energy)
        
        i=0
        
        while i < imax:
        
            print(i," Energy: ",data2[i].energy," Magnetization: ",data2[i].magnetization)
        
            i+=1

        '''
        ##  Compare levels and determine a list of energy level differences in 
        ##  ascending order

        diffList = []
        
        for i in energyData:
            for j in energyData:
            
                if i > j:
                    
                    diffList.append(i-j)
                    
        diffList.sort()
        
        print("Lowest lying single particle states:",diffList)
        '''
        
################################################################################
######      ANALYSE THE HIGHEST AMPLITUDE FOCK STATES       ####################
        
    if args.analyse_most_probable:

        sqlDb = SQLite(DATA_FOLDER+args.FolderName()+args.SqlName()+'.sql')
    
        if args.cut_v0:    ##  keep g fixed
        
            colNames,rows = sqlDb.GetData("Id,V0,Interaction,OutFileName,OffsetX,OffsetY",SQL_TABLE_NAME,"Completed=1 AND GotMostProbable=1 AND Interaction="+args.interaction)
        
        elif args.cut_g:   ##  keep v0 fixed
            
            colNames,rows = sqlDb.GetData("Id,V0,Interaction,OutFileName,OffsetX,OffsetY",SQL_TABLE_NAME,"Completed=1 AND GotMostProbable=1 AND V0="+args.v0+" AND OffsetX="+args.offset_x+" AND OffsetY="+args.offset_y)
    
        else:
        
            print("ERROR in analyse_most_probable: need to select cut-v0 or cut-g")
            quit()
    
        sqlDb.PrintData(colNames,rows)

        data = []
        labelList = []
        v0Data = []
        interactionData = []
        ##  Obtain plot data

        for row in rows:

            fileName = DATA_FOLDER+args.FolderName()+row[3]+args.SectorName()+'_most_probable.dat'

            args.offset_x = str(row[4])
            args.offset_y = str(row[5])

            ##  Read in the list of most probable configurations
            
            fin  = open(fileName)

            breakCounter = 0

            for line in fin:

                line2 = line.split()

                breakCounter += 1
                
                if breakCounter >0:
                ##>=6*3 and breakCounter <=7*3:

                    ##  Ignore empty lines
                    if len(line2) == 0:
                        pass
                    ##  Keep comment lines in a separate list
                    elif line2[0] == "##":
                        labelList.append(' '.join(line2[1:]))
                    else:
                        currData = ConfigurationAmplitudeData(args.x,args.y)
                        
                        currData.fockState = int(line2[0])
                        currData.amplitude = complex(float(line2[1]),float(line2[2]))

                        data.append(currData)
                        
                        v0Data.append(row[1])
                        interactionData.append(row[2])
                        
                        breakCounter += 1

            fin.close()

        ##  Also read in the magnetization map
        
        ##xData,yData,magnetizationData = ObtainMagnetizationMapData()

        ##  Make an interactive plot showing the most probable configurations
        ##  as a function of system parameters
        
        print ("\n========== Generating most probable configuration plot ==========")

        ##  SET PLOT PARAMETERS HERE

        font_size = 30
        x_tick_pad = 15
        start_pos = 0

        main_axes = plt.axes([0.05, 0.2, 0.9, 0.6]) 	# ([left, bottom, width, height])
        main_axes.set_aspect('equal')
        
        PlotUnitCellBackground(main_axes)
        
        PlotSimulationGrid(main_axes)
        
        fig = main_axes.get_figure()

        xData,yData = data[start_pos].GetConfiguration()

        main_axes.set_title(r'$V_0 = '+str(v0Data[start_pos])+'$  $g_{2D} = '+str(interactionData[start_pos])+'$'+r' Sector '+str(args.sector[0])+' '+str(args.sector[1])+'\n'+labelList[start_pos]+'\n'+r"Amplitude: "+str(data[start_pos].amplitude.real)+" "+str(data[start_pos].amplitude.imag)+"I",fontsize=font_size-10)

        p1 = main_axes.scatter(xData,yData,marker="o",color=COLOUR_LIST[0],s=250)

        ##  Add a discrete slider to control the model parameter values

        axis_slider = fig.add_axes([0.35, 0.9, 0.3, 0.03])# ([left, bottom, width, height])
        
        slider = DiscreteSlider(axis_slider, '',0,len(interactionData),allowed_vals=range(0,len(interactionData)), valinit=start_pos)

        ######      Define a function to update the figure 
        ######      for the ith set of data 
        def update_fig(val):

            global p1
  
            i = int(slider.val)

            ##  Update the title
            main_axes.set_title(r'$V_0 = '+str(v0Data[i])+'$  $g_{2D} = '+str(interactionData[i])+'$'+r' Sector '+str(args.sector[0])+' '+str(args.sector[1])+'\n'+labelList[i]+'\n'+r"Amplitude: "+str(data[i].amplitude.real)+" "+str(data[i].amplitude.imag)+"I",fontsize=font_size-10)
            
            ##  Update x and y data in the plot

            p1.remove()

            xData,yData = data[i].GetConfiguration()

            p1 = main_axes.scatter(xData,yData,marker="o",color=COLOUR_LIST[0],s=250)
         
            fig.canvas.draw_idle()

        def press(event):

            if event.key=='d':

                slider.set_val(slider.previous_val+1)
                
                update_fig(slider)
                
            if event.key=='a':

                slider.set_val(slider.previous_val-1)
                
                update_fig(slider)   

        slider.on_changed(update_fig)

        fig.canvas.mpl_connect('key_press_event',press)

        ##  Resize plt.show to a maximized window
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())

        plt.show()

################################################################################
######      MAKE PLOTS OF THE OVERLAPS WITH GIVEN EIGENSTATES AS A      ########
######      FUNCTION OF MODEL PARAMETERS                                ########

    if args.plot_overlaps:
    
        ##  Read in overlap data for a given set of Fock states and plot
        ##  the overlaps as a function of system parameters
    
        pass

################################################################################       
######      MAKE PLOTS OF THE SPATIAL WAVE FUNCTION DISTRIBUTION        ########

    if args.plot_spatial_wavefunction:

        ##  Make an interactive plot of the probability distribution
        ##  on a real space grid rx,ry, with a slider as a function of
        ##  kx,ky
        
        dataArray = ObtainSpatialWaveFunction()
            
        ######  Make plots of the spatial wave function distribution
        
        ##  SET PLOT PARAMETERS HERE

        font_size  = 30
        x_tick_pad = 15
        start_pos  = 0

        main_axes = plt.axes([0.05, 0.2, 0.9, 0.6]) 	# ([left, bottom, width, height])
        main_axes.set_aspect('equal')
        
        print(dataArray[start_pos].realSpaceX)
        print(dataArray[start_pos].realSpaceY)
        print(dataArray[start_pos].values)
        
        triang = tri.Triangulation(dataArray[start_pos].realSpaceX,dataArray[start_pos].realSpaceY)

        p1 = main_axes.tricontourf(triang,dataArray[start_pos].values,cmap=COLOUR_MAP)
        
        fig = main_axes.get_figure()
        
        main_axes.set_xlabel(r"$r_x$",fontsize=font_size)
        main_axes.set_xlim([0,float(gridSize)-1])
        
        main_axes.set_ylabel(r"$r_y$",fontsize=font_size)
        main_axes.set_ylim([0,float(gridSize)-1])
        
        ##  Add a discrete slider to control the model parameter values

        axis_slider = fig.add_axes([0.35, 0.9, 0.3, 0.03])# ([left, bottom, width, height])
        
        slider = DiscreteSlider(axis_slider, '',0,myRange-1,allowed_vals=range(0,myRange),valinit=start_pos)

        ######      Define a function to update the figure 
        ######      for the ith set of data 
        def update_fig(val):

            global p1
            global dataArray
  
            i = int(slider.val)

            ##  Remove previous set of contour data
            for coll in p1.collections:
                main_axes.collections.remove(coll)
            
            triang = tri.Triangulation(dataArray[i].realSpaceX,dataArray[i].realSpaceY)
            
            ##  Update the probability contour data  
            p1 = main_axes.tricontourf(triang,dataArray[i].values,cmap=COLOUR_MAP)

        def press(event):

            if event.key=='d':

                if slider.previous_val < myRange-1:

                    slider.set_val(slider.previous_val+1)
                    
                    update_fig(slider)
                
            if event.key=='a':

                if slider.previous_val > 0:

                    slider.set_val(slider.previous_val-1)
                    
                    update_fig(slider)

        slider.on_changed(update_fig)

        fig.canvas.mpl_connect('key_press_event',press)

        ##  Resize plt.show to a maximized window
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())

        plt.show()

################################################################################
######      PLOT THE ENERGY DERIVATIVE AS A FUNCITON OF DENSITY     ############

    if args.plot_energy_derivative:
    
        ##  Obtain the energy eigenvalues at a fixed point in parameter space
        ##  but for all values of the particle number

        fillingFactors1 = []
        secondDerivs1 = []
        v0List = [1.01]
        gList = [0.01,0.26,0.51,0.76,1.01]
        
        ##[0.01,0.51,1.01,1.51]
        
        for v0 in v0List:
            for g in gList:
        
                listOfEnergies = []
            
                for nbr in range(2,args.x*args.y-1):
                    
                    ##  Update args value of nbr so that the correct file
                    ##  names are constructed
                    
                    args.nbr = nbr
                    
                    ##  Read from SQL database

                    sqlDb = SQLite(DATA_FOLDER+args.FolderName()+args.SqlName()+'.sql')

                    colNames,rows = sqlDb.GetData("V0,Interaction,OutFileName",SQL_TABLE_NAME,"Completed=1 AND Interaction="+str(g)+" AND V0="+str(v0))

                    sqlDb.PrintData(colNames,rows)
                    
                    for row in rows:

                        ##  Read eigenvalue data from the file name indicated
                        ##  Loop over all available sectors

                        spectrumData = SpectrumDataArray()

                        for sectorX in range(0,args.x):
                        
                            for sectorY in range(0,args.y):
                                
                                ##  Set sector so that the file name is correctly generated
                                ##  by the args.FolderName() function
                                args.sector = [sectorX,sectorY]

                                fileName = DATA_FOLDER+args.FolderName()+str(row[2])+args.SectorName()+'_eigensystem.dat'
                                
                                print("READING DATA FROM FILE: "+fileName)
                                
                                try:
                                    with open(fileName) as file:
                                        
                                        fin = open(fileName)
                                
                                        fin.readline()  ##  skip first comment line
                                        
                                        nbrEigenvalues = int(fin.readline())

                                        fin.readline()  ##  skip first comment line
                                        
                                        ##  Read in data
                                        
                                        for i in range(0,min(nbrEigenvalues,args.nbr_eigs)):
                                            
                                            dataPoint = SpectrumDataPoint()
                                            
                                            dataPoint.energy      = float(fin.readline())/0.5
                                            dataPoint.interaction = float(row[1])
                                            ##  Divide by recoil energy
                                            dataPoint.v0          = float(row[0])/0.5      
                                            dataPoint.sectors     = [[sectorX,sectorY]]
                                            dataPoint.eigLabels   = [i]
                                            
                                            spectrumData.Append(dataPoint)
                                            
                                        fin.close()    
                                        
                                except IOError as e:
                                    print("ERROR READING FILE "+fileName)
                                
                                    pass
                        
                        if not spectrumData.GetNbrPoints():

                            print("ERROR: Data not found in table")
                            
                            quit()
                        
                        ##  Sort the spectrum data to obtain the minimum energy sector
                        ##  and eigenvalue label
                        
                        spectrumData.Sort()

                        ##  Obtain the minimum energy value
                        
                        listOfEnergies.append(spectrumData.GetEnergy(0))
                
                print(listOfEnergies)
                        
                ##  Generate values of the second derivative and the filling
                ##  factor
                
                fillingFactors = []
                secondDerivs = []
                
                for nbr in range(1,args.x*args.y-4):

                    secondDerivs.append((listOfEnergies[nbr-1]+listOfEnergies[nbr+1]-2.0*listOfEnergies[nbr])/2.0)
                    
                    fillingFactors.append(float(nbr+2)/float(args.x*args.y))
                    
                fillingFactors = np.array(fillingFactors)
                secondDerivs   = np.array(secondDerivs)
                
                fillingFactors1.append(fillingFactors)
                secondDerivs1.append(secondDerivs)
        
        print(secondDerivs1)
        
        ##  Plot second derivative vs filling factor
        
        ##  SET PLOT PARAMETERS HERE

        font_size = 28

        fig_name = 'kx_'+str(args.x)+'_ky_'+str(args.y)+'_compressibility.pdf'
        
        main_axes = plt.axes([0.1, 0.1, 0.5, 0.9]) 	# ([left, bottom, width, height])
        ##main_axes = plt.axes([0.1, 0.1, 0.9, 0.5]) 	# ([left, bottom, width, height])
        main_axes.get_xaxis().tick_bottom()
        main_axes.get_yaxis().tick_left()
        
        ##main_axes.set_title(r'$\frac{V_0}{E_R} = '+str(float(args.v0)/0.5)+'$  $g_{2D} = '+str(args.interaction)+'$',fontsize=font_size)
        
        for i in range(0,len(fillingFactors1)):
        
            ##main_axes.plot(fillingFactors1[i],secondDerivs1[i],marker=MARKER_LIST[i+1],color=COLOUR_LIST[i],label=str(float(v0List[i])/0.5))
            
            main_axes.plot(fillingFactors1[i],secondDerivs1[i],marker=MARKER_LIST[i+3],color=COLOUR_LIST[i],ms=10,label=str(float(gList[i])))

        main_axes.tick_params(direction='in',pad=10)

        main_axes.set_xlabel(r'$\nu$',fontsize=font_size)
        main_axes.set_xlim([5.0/30,27.0/30])

        main_axes.set_xticks([0.2,0.3,0.4,0.5,0.6,0.7,0.8])
        main_axes.set_xticklabels([r'$\textbf{0.2}$',r'$\textbf{0.3}$',r'$\textbf{0.4}$',r'$\textbf{0.5}$',r'$\textbf{0.6}$',r'$\textbf{0.7}$',r'$\textbf{0.8}$'],fontsize=font_size-10)

        ##main_axes.set_ylabel(r'$\frac{\partial^2 E}{E_R\partial N^2}$ [Incompressibility] ',fontsize=font_size-4)
        
        main_axes.set_ylabel(r'$\beta^{-1}(\nu)$ ',fontsize=font_size-4)
        
        if args.formatting==1:
            main_axes.annotate(r'$\textbf{1}/\textbf{3}$', xy=(0.333,9),
                horizontalalignment='center', verticalalignment='center',fontsize=font_size-10)
        
            main_axes.annotate(r'$\textbf{1}/\textbf{2}$', xy=(0.5,18.75),
                    horizontalalignment='center', verticalalignment='center',fontsize=font_size-10)
            
            main_axes.set_ylim([-1,25])
            main_axes.set_yticks([0,5,10,15,20,25])
            main_axes.yaxis.set_minor_locator(MultipleLocator(1))
            main_axes.set_yticklabels([r'$\textbf{0}$',r'$\textbf{5}$',r'$\textbf{10}$',r'$\textbf{15}$',r'$\textbf{20}$',r'$\textbf{25}$'],fontsize=font_size-10)
            main_axes.xaxis.set_label_coords(1.0,-0.01)
            main_axes.xaxis.set_minor_locator(MultipleLocator(1/30.0))
                
        if args.formatting==2:
        
            main_axes.annotate(r'$\textbf{1}/\textbf{3}$', xy=(0.333,10),
                    horizontalalignment='center', verticalalignment='center',fontsize=font_size-10)
                    
            ##main_axes.annotate(r'$\textbf{2}/\textbf{5}$', xy=(0.43,30),
                    ##horizontalalignment='center', verticalalignment='center',fontsize=font_size-10)
                    
            main_axes.annotate(r'$\textbf{1}/\textbf{2}$', xy=(0.5,19),
                    horizontalalignment='center', verticalalignment='center',fontsize=font_size-10)
                    
            ##main_axes.annotate(r'$\textbf{3}/\textbf{5}$', xy=(0.64,30),
                    ##horizontalalignment='center', verticalalignment='center',fontsize=font_size-10)
                    
            main_axes.set_ylim([-1,22])
            main_axes.set_yticks([0,5,10,15,20])
            main_axes.set_yticklabels([r'$\textbf{0}$',r'$\textbf{5}$',r'$\textbf{10}$',r'$\textbf{15}$',r'$\textbf{20}$'],fontsize=font_size-10)
            main_axes.xaxis.set_label_coords(1.0,-0.01)

        ##main_axes.legend(loc='best',numpoints=1,scatterpoints=1,title=r'$V_0/E_R$',fontsize=font_size-10)
        main_axes.legend(loc='best',numpoints=1,scatterpoints=1,ncol=3,title=r'$\tilde{g}/2\pi$',fontsize=font_size-14)
        main_axes.get_legend().get_title().set_fontsize(font_size-8)
                
        ##plt.show()
        
        plt.savefig(fig_name,bbox_inches='tight')
        
        print ("\nDONE! find plot at "+fig_name)

################################################################################
######      PLOT THE ENERGY LEVELS VS SECTOR      ##############################

    if args.plot_energy_vs_sector:
    
        ##  Obtain energy values for a given set of system parameters
        
        spectrumData = SpectrumDataArray()

        ##  Read from SQL database

        sqlDb = SQLite(DATA_FOLDER+args.FolderName()+args.SqlName()+'.sql')

        colNames,rows = sqlDb.GetData("OutFileName",SQL_TABLE_NAME,"Completed=1 AND Interaction="+args.interaction+" AND V0="+args.v0+" AND OffsetX="+str(args.offset_x)+" AND OffsetY="+str(args.offset_y))  

        sqlDb.PrintData(colNames,rows)

        for row in rows:

            ##  Read eigenvalue data from the file name indicated
            ##  Loop over all available sectors
            
            for sectorX in range(0,args.x):
            
                for sectorY in range(0,args.y):
                    
                    #if [sectorX,sectorY] in [[0,2],[2,0],[2,2],[0,0]]:
                    
                    ##  Set sector so that the file name is correctly generated
                    ##  by the args.FolderName() function
                    args.sector = [sectorX,sectorY]

                    fileName = DATA_FOLDER+args.FolderName()+row[0]+args.SectorName()+'_eigensystem.dat'
                    
                    ##print("READING DATA FROM FILE: "+fileName)
                    
                    try:
                        with open(fileName) as file:
                            
                            fin = open(fileName)
                    
                            fin.readline()  ##  skip first comment line
                            
                            nbrEigenvalues = int(fin.readline())

                            fin.readline()  ##  skip first comment line
                            
                            ##  Read in data
                            
                            for i in range(0,min(nbrEigenvalues,args.nbr_eigs)):
                                
                                dataPoint = SpectrumDataPoint()
                                
                                dataPoint.energy      = float(fin.readline())/0.5
                                dataPoint.interaction = args.interaction
                                dataPoint.v0          = float(args.v0)/0.5       ##  Divide by recoil energy
                                dataPoint.sectors     = [[sectorX,sectorY]]
                                dataPoint.eigLabels   = [i]
                                
                                spectrumData.Append(dataPoint)
                                
                            fin.close()
                            
                    except IOError as e:
                        pass

        if not spectrumData.GetNbrPoints():

            print("ERROR: Data not found in table")
            
            quit()
            
        ##  Get maximum and minimum values of energy, interaction and V0
            
        maxEnergy,maxInteraction,maxV0 = spectrumData.GetMax()
        minEnergy,minInteraction,minV0 = spectrumData.GetMin()
        
        energyStep = 0.02*(maxEnergy - minEnergy)

        ##  Finally we can make some plots

        font_size = 28
        
        main_axes = plt.axes([0.09,0.9,0.9,0.7])	# ([left, bottom, width, height])
        
        main_axes.get_xaxis().tick_bottom()
        main_axes.get_yaxis().tick_left()

        ##  Plot main data set
        
        print ("\n========== Generating energy level plot ==========")
        
        dataEnergy = spectrumData.GetAllEnergyList()
        dataSectors = spectrumData.GetUnfoldedSectorList(args.y)

        ##print(dataEnergy)
        ##print(dataSectors)

        main_axes.plot(dataSectors,dataEnergy,'_',ms=10)
        
        ##main_axes.set_title(r'Interaction = '+str(float(args.interaction))+r',  $\frac{V_0}{E_R} $= '+str(float(args.v0)/0.5),fontsize=20)

        main_axes.set_xlabel(r'$\mathbf{k^{\mbox{\large tot.}}} = k_1^{\mbox{\large tot.}}+L_1 \times k_2^{\mbox{\large tot.}}$',fontsize=font_size,labelpad=10)
        main_axes.set_xlim(-0.5,args.x*args.y-0.5)
        main_axes.set_xticks(range(0,args.x*args.y,args.x))
        main_axes.set_xticklabels([r'$\textbf{0}$',r'$\textbf{5}$',r'$\textbf{10}$',r'$\textbf{15}$',r'$\textbf{20}$',r'$\textbf{25}$'],fontsize=font_size-10)
        main_axes.set_xticks(range(0,args.x*args.y),minor=True)

        main_axes.set_ylabel(r'$E/E_R$',fontsize=font_size)
        ##main_axes.set_ylim(minEnergy-energyStep,maxEnergy+energyStep)
        
        main_axes.set_ylim(31.5,36)
        main_axes.set_yticks([32,33,34,35,36])
        main_axes.set_yticklabels([r'$\textbf{32.0}$',r'$\textbf{33.0}$',r'$\textbf{34.0}$',r'$\textbf{35.0}$',r'$\textbf{36.0}$'],fontsize=font_size-10)
        
        main_axes.minorticks_on()

        ##for t in main_axes.get_yticklabels():
            ##t.set_fontsize(font_size-10)
                
        ##for t in main_axes.get_xticklabels():
           ## t.set_fontsize(font_size-10)
        
        fig_name = 'kx_'+str(args.x)+'_ky_'+str(args.y)+'_n_'+str(args.nbr)+'_energy_vs_sector.pdf'
        
        plt.savefig(fig_name,bbox_inches='tight')
        
        print ("\nDONE! find plot at "+fig_name)

################################################################################
######      PLOT A GIVEN OBSERVABLE OF THE GROUND STATE AS A        ############
######      FUNCTION OF SECTOR FOR FIXED V0 AND G                   ############

    if args.plot_combined_observable_vs_sector:

        ##  Read in the energy level data and determine which sector has the
        ##  lowest energy at each point in parameter space. Then obtain the 
        ##  observable data for that sector
        
        if args.observe_susceptibility:
        
            gotFlag     = "GotDensityDensity"
            fileLabel   = "susceptibility"
            plotLabel   = "Susceptibility"
        
        elif args.observe_r:
        
            gotFlag     = "GotDensityDensity"
            fileLabel   = "r_parameter"
            plotLabel   = "R-60"
            
        elif args.observe_pomer:
        
            gotFlag     = "GotDensityDensity"
            fileLabel   = "pomeranchuk_parameter"
            plotLabel   = "Centralness"
        
        elif args.observe_participation:
        
            gotFlag     = "GotParticipationRatio"
            fileLabel   = "participation_ratio"
            plotLabel   = "Participation"
        else:
        
            print("ERROR in plot_combined_observable: no observable selected")
            quit()
            
        observableData  = []
        v0Data          = []
        interactionData = []
        
        sqlDb = SQLite(DATA_FOLDER+args.FolderName()+args.SqlName()+'.sql')
    
        colNames,rows = sqlDb.GetData("OutFileName",SQL_TABLE_NAME,"Completed=1 AND "+gotFlag+"=1 AND Interaction="+args.interaction+" AND V0="+args.v0)

        sqlDb.PrintData(colNames,rows)

        for row in rows:

            ##  Read eigenvalue data from the file name indicated
            ##  Loop over all available sectors

            spectrumData = SpectrumDataArray()

            for sectorX in range(0,args.x):
            
                for sectorY in range(0,args.y):
                    
                    ##  Set sector so that the file name is correctly generated
                    ##  by the args.FolderName() function
                    args.sector = [sectorX,sectorY]

                    fileName = DATA_FOLDER+args.FolderName()+str(row[0])+args.SectorName()+'_eigensystem.dat'
                    
                    print("READING DATA FROM FILE: "+fileName)
                    
                    try:
                        with open(fileName) as file:
                            
                            fin = open(fileName)
                    
                            fin.readline()  ##  skip first comment line
                            
                            nbrEigenvalues = int(fin.readline())

                            fin.readline()  ##  skip first comment line
                            
                            ##  Read in data
                            
                            for i in range(0,min(nbrEigenvalues,args.nbr_eigs)):
                                
                                dataPoint = SpectrumDataPoint()
                                
                                dataPoint.energy      = float(fin.readline())/0.5
                                dataPoint.interaction = float(args.interaction)
                                ##  Divide by recoil energy
                                dataPoint.v0          = float(args.v0)/0.5    
                                dataPoint.sectors     = [[sectorX,sectorY]]
                                dataPoint.eigLabels   = [i]

                                spectrumData.Append(dataPoint)

                            fin.close()    
                            
                    except IOError as e:
                        print("ERROR READING FILE "+fileName)
                    
                        pass
            
            if not spectrumData.GetNbrPoints():

                print("ERROR: Data not found in table")
                
                quit()
         
            ##  Sort the spectrum data to obtain the minimum energy sector
            ##  and eigenvalue label
            
            spectrumData.Sort()

            ##  Use the first value in the spectrumData array to set the
            ##  lowest energy sector and eigenvalue label (allowing
            ##  for a number of degenerate states)
            
            for i in range(0,min(spectrumData.GetNbrPoints(),200)):
            
                sectors   = spectrumData.GetSectors(i)
                eigLabels = spectrumData.GetEigLabels(i)

                index = 0;

                observableMean = 0.0

                for s in sectors:

                    args.sector = s

                    eigLabel = eigLabels[index]

                    index += 1

                    ##  Read in the corresponding susceptibility data from a file
                    ##  along with the associated model parameter values

                    observableFileName = DATA_FOLDER+args.FolderName()+row[0]+args.SectorName()+'_eig_'+str(eigLabel)+'_'+fileLabel+'.dat'

                    print("READING DATA FROM FILE: "+observableFileName)

                    try:
                        with open(observableFileName) as file:
                            
                            fin2 = open(observableFileName)
                            
                            observableData.append(float(fin2.readline()))
                        
                            fin2.close()
                        
                    except IOError as e:
        
                        print("ERROR READING FILE "+observableFileName)
                        quit()
                        
                        pass
                        
        ##  Convert to np array
        
        observableData = np.array(observableData)
        
        print(len(observableData))
        
        dataSectors = spectrumData.GetUnfoldedSectorList(args.y)
        
        dataSectors = dataSectors[0:len(observableData)]
        
        dataEnergy = spectrumData.GetAllEnergyList()
        
        dataEnergy = dataEnergy[0:len(observableData)]
        
        print(len(dataSectors))
        
        print(len(dataEnergy))
        
        print("AVERAGE OF OBSERVABLE = "+str(np.average(observableData)))
        
        ##  Finally we can make some plots

        main_axes = plt.axes([0.09,0.14,0.75,0.75])	# ([left, bottom, width, height])
        
        p1 = main_axes.scatter(dataSectors,observableData,marker='s',c=dataEnergy,cmap=COLOUR_MAP)
        
        fig = main_axes.get_figure()
        
        ##  Plot the occupation probability colour bar
        cb  = fig.colorbar(p1, shrink=0.8)
        cb.set_label(r"\textbf{Energy}",fontsize=20,rotation=0)
        
        cb.ax.yaxis.set_label_coords(1.0,1.17)

        main_axes.set_title(r'Interaction = '+str(float(args.interaction))+r',  $\frac{V_0}{E_R} $= '+str(float(args.v0)/0.5),fontsize=20)
        main_axes.set_xlabel(r'$k_y+N_y \times k_x$',fontsize=20)
        main_axes.set_xlim(-1,args.x*args.y)
        main_axes.set_xticks(range(0,args.x*args.y+1,args.y))

        main_axes.set_ylabel(plotLabel,fontsize=20)
        ##main_axes.set_ylim(minEnergy-energyStep,maxEnergy+energyStep)

        for t in main_axes.get_yticklabels():
            t.set_fontsize(30)
                
        for t in main_axes.get_xticklabels():
            t.set_fontsize(30)
        
        fig_name = 'kx_'+str(args.x)+'_ky_'+str(args.y)+'_n_'+str(args.nbr)+'_'+fileLabel+'_vs_sector.pdf'

        plt.savefig(fig_name,bbox_inches='tight')
        
        print ("\nDONE! find plot at "+fig_name)

################################################################################
######      PLOT THE ENERGY LEVELS VS OFFSET PARAMETER      ####################

    if args.plot_energy_vs_offset:
        
        ##  Obtain the energy values as a function of the x offset parameter
        
        listOfEnergies = []
        listOfOffsets  = []
        
        energyGapList   = []
        energyWidthList = []
        
        listOfEnergies1 = {}
        listOfOffsets1 = {}
        specialSectors = []
        
        if args.formatting==1:
            specialSectors = [[0,1],[0,3],[0,5]]
            ##[[0,0],[0,2],[0,4]]
        elif args.formatting==2:
            specialSectors = [[0,0],[0,args.y/2],[args.x/2,0],[args.x/2,args.y/2]]

        for sector in specialSectors:
            listOfEnergies1[sector[0]+args.x*sector[1]] = []
            listOfOffsets1[sector[0]+args.x*sector[1]] = []
        
        sqlDb = SQLite(DATA_FOLDER+args.FolderName()+args.SqlName()+'.sql')

        colNames,rows = sqlDb.GetData("offsetX,offsetY,OutFileName",SQL_TABLE_NAME,"Completed=1 AND Interaction="+args.interaction+" AND V0="+args.v0)
        
        sqlDb.PrintData(colNames,rows)
        
        for row in rows:

            ##  Read eigenvalue data from the file name indicated
            ##  Loop over all available sectors

            energyBySector = []

            for sectorX in range(0,args.x):
            
                for sectorY in range(0,args.y):
                    
                    ##  Set sector so that the file name is correctly generated
                    ##  by the args.FolderName() function
                    args.sector = [sectorX,sectorY]

                    fileName = DATA_FOLDER+args.FolderName()+str(row[2])+args.SectorName()+'_eigensystem.dat'
                    
                    print("READING DATA FROM FILE: "+fileName)
                    
                    try:
                        with open(fileName) as file:
                            
                            fin = open(fileName)
                    
                            fin.readline()  ##  skip first comment line
                            
                            nbrEigenvalues = int(fin.readline())

                            ##args.nbr_eigs = nbrEigenvalues

                            fin.readline()  ##  skip first comment line
                            
                            ##  Read in data
                            
                            for i in range(0,min(nbrEigenvalues,args.nbr_eigs)):

                                energy = float(fin.readline())/0.5
                                offsetX = float(row[0])
                                offsetY = float(row[1])
                                
                                offset = math.sqrt(offsetX**2+offsetY**2)
                            
                                if args.sector in specialSectors:
                                    listOfEnergies1[args.sector[0]+args.x*args.sector[1]].append(energy)
                                    listOfOffsets1[args.sector[0]+args.x*args.sector[1]].append(offset)     
                                ##else:
                                listOfEnergies.append(energy)
                                listOfOffsets.append(offset)

                                energyBySector.append(energy)

                            fin.close()

                    except IOError as e:
                        print("ERROR READING FILE "+fileName)
                    
                        pass

            energyBySector.sort()
            
            ##  Determine the energy spread of the low lying states
            
            maxEnergy = 0
            minEnergy = 1000
            
            for energy in energyBySector:
                if energy < float(args.gap_energy):
                    if energy > maxEnergy:
                        maxEnergy = energy
                    if energy < minEnergy:
                        minEnergy = energy
                        
            ##  Calculate the energy gap
        
            highestLowEnergyState = 0
            lowestHighEnergyState = 1000

            for energy in energyBySector:
                if energy < float(args.gap_energy):
                    if energy > highestLowEnergyState:
                        highestLowEnergyState = energy
                elif energy < lowestHighEnergyState:
                        lowestHighEnergyState = energy

            energyGapList.append(maxEnergy-minEnergy)
            energyWidthList.append(lowestHighEnergyState-highestLowEnergyState)
            
            print("energyBySector ",energyBySector)

        print("===============================================================")
        
        ##print(listOfEnergies1[0])
        
        ##print(listOfEnergies)
        ##print(listOfOffsets)
        
        ##print("list of energies[specialSectors[0]]")
        
        ##print(listOfEnergies1[specialSectors[0][0]+args.x*specialSectors[0][1]])
        
        originalEnergy = listOfEnergies[:]

        originalEnergy1 = copy.deepcopy(listOfEnergies1)
        
        ##  Renormalize by the minimum energy value from the special sectors,
        ##  for each value of the offset
        
        uniqueOffsets = {}
    
        for offset in listOfOffsets:
            uniqueOffsets[offset] = 1
        
        for sector in specialSectors:   
            for offset in listOfOffsets1[sector[0]+args.x*sector[1]]:
                uniqueOffsets[offset] = 1
        
        uniqueOffsets = sorted(uniqueOffsets.keys()[:])

        ##print("unique offsets")

        ##print(uniqueOffsets)
        
        ##print("lowest values")
        
        counter = 0
        
        for testOffset in uniqueOffsets:
            
            minValues = []
            
            for sector in specialSectors:
                energies1 = listOfEnergies1[sector[0]+args.x*sector[1]][counter*args.nbr_eigs:(counter+1)*args.nbr_eigs]
                
                print(energies1)
                
                energies1 = sorted(energies1[:])
                
                minValues.append(energies1[0])
            
            minValues = sorted(minValues[:])
            
            ##print(minValues[0:3])
            
            ##  Mean of lowest 3
            if args.formatting==1:
                minEnergy = (minValues[0]+minValues[1]+minValues[2])/3.0
        
            ##  Mean of lowest 2
            ##minEnergy = (minValues[0]+minValues[1])/2.0
        
            ##  Lowest 1
            elif args.formatting==2:
                minEnergy = minValues[0]
            
            ##print(minEnergy)
            
            for sector in specialSectors:
                for val in range(counter*args.nbr_eigs,(counter+1)*args.nbr_eigs):
                    listOfEnergies1[sector[0]+args.x*sector[1]][val] -= minEnergy
            
            for val in range(counter*args.nbr_eigs*args.x*args.y,(counter+1)*args.nbr_eigs*args.x*args.y):
                listOfEnergies[val] -= minEnergy
            
            counter += 1
        
        ##print("===============================================================")
        
        ##print(listOfEnergies)
        
        ##print(listOfEnergies1[0])
        
        ##print(listOfOffsets1[0])  
        
        ######      Plot energy vs offset       ####################################################
        
        font_size = 28

        fig_name = 'kx_'+str(args.x)+'_ky_'+str(args.y)+'_n_'+str(args.nbr)+'_energy_vs_offset.pdf'
        
        sq3 = math.sqrt(3.0)
        
        if args.formatting==2:
            nPoints = 10.0
            nx = float(args.y)/1.0
        else:
            nPoints = 12.0
            nx = float(args.y)/3.0

        xticks = [0.0,sq3/(nPoints*nx),2*sq3/(nPoints*nx),3*sq3/(nPoints*nx),4*sq3/(nPoints*nx),5*sq3/(nPoints*nx),6*sq3/(nPoints*nx),7*sq3/(nPoints*nx),8*sq3/(nPoints*nx),9*sq3/(nPoints*nx),10*sq3/(nPoints*nx),11*sq3/(nPoints*nx),12*sq3/(nPoints*nx)][0:int(nPoints)+1]

        if args.formatting==2:
        
            xticklabels = [r'$\textbf{0}$',r'$\textbf{0.1}$',r'$\textbf{0.2}$',r'$\textbf{0.3}$',r'$\textbf{0.4}$',r'$\textbf{0.5}$',r'$\textbf{0.6}$',r'$\textbf{0.7}$',r'$\textbf{0.8}$',r'$\textbf{0.9}$',r'$\textbf{1.0}$']
        
        else:
        
            xticklabels = [r'$\textbf{0}$',r'$\textbf{0.25}$',r'$\textbf{0.5}$',r'$\textbf{0.75}$',r'$\textbf{1.0}$',r'$\textbf{1.25}$',r'$\textbf{1.5}$',r'$\textbf{1.75}$',r'$\textbf{2.0}$',r'$\textbf{2.25}$',r'$\textbf{2.5}$',r'$\textbf{2.75}$',r'$\textbf{3.0}$']
        
        inset_xticks = xticks[::2]
        inset_xticklabels = xticklabels[::2]
        
        if args.formatting==1:

            main_axes = plt.axes([0.05, 0.05, 0.9, 0.7]) 	# ([left, bottom, width, height])
            main_axes.get_xaxis().tick_bottom()
            main_axes.get_yaxis().tick_left()
            inset_axes=plt.axes([0.23,0.21,0.46,0.28])	# ([left, bottom, width, height])
            inset_axes.get_xaxis().tick_bottom()
            ##inset_axes.get_yaxis().tick_left()
            
            ##  Plot "other" sectors
            main_axes.plot(listOfOffsets,originalEnergy,'_',ms=5,color=COLOUR_LIST[7],label="Other")
            
            i=1

            for sector in specialSectors:
        
                main_axes.plot(listOfOffsets1[sector[0]+args.x*sector[1]],originalEnergy1[sector[0]+args.x*sector[1]],'_',ms=8,color=COLOUR_LIST[i])
                main_axes.scatter(listOfOffsets1[sector[0]+args.x*sector[1]],originalEnergy1[sector[0]+args.x*sector[1]],marker=MARKER_LIST[i+2],s=24,color=COLOUR_LIST[i],label=sector)
                
                inset_axes.plot(listOfOffsets1[sector[0]+args.x*sector[1]],listOfEnergies1[sector[0]+args.x*sector[1]],'_',ms=8,color=COLOUR_LIST[i])
                inset_axes.scatter(listOfOffsets1[sector[0]+args.x*sector[1]],listOfEnergies1[sector[0]+args.x*sector[1]],marker=MARKER_LIST[i+2],s=24,color=COLOUR_LIST[i],label=sector)
                
                i+=1
            
            ##main_axes.set_xlabel(r"$\Delta_{\kappa'_1}/\sqrt{3}$",fontsize=font_size)
            main_axes.set_xlabel(r'$N_{\phi}$',fontsize=font_size)
            main_axes.set_ylabel(r'$E/E_R$',fontsize=font_size)
            inset_axes.set_xlabel(r'$\Delta_{kx}/\sqrt{3}$',fontsize=font_size-10)
            inset_axes.set_xlabel(r'$N_{\phi}$',fontsize=font_size-10)
            ##inset_axes.set_ylabel(r'$[E-E_0(\Delta_{kx})]/E_R$',fontsize=font_size-10)
            inset_axes.set_ylabel(r'$[E-E_0(N_{\phi})]/E_R$',fontsize=font_size-10)
            
            ##main_axes.set_title(r'$\frac{V_0}{E_R} = '+str(float(args.v0)/0.5)+'$  $g_{2D} = '+str(args.interaction)+'$',fontsize=font_size)

            main_axes.set_xlim(-0.01,math.sqrt(3)/2.0+0.01)
            main_axes.set_xticks(xticks)
            main_axes.set_xticklabels(xticklabels,fontsize=font_size-10)
            main_axes.set_ylim([35.5,42.5])
            ##main_axes.set_ylim([28.0,33.5])
            main_axes.set_yticks([36,37,38,39,40,41,42])
            main_axes.minorticks_on()
            main_axes.set_yticklabels([r'$\textbf{36.0}$',r'$\textbf{37.0}$',r'$\textbf{38.0}$',r'$\textbf{39.0}$',r'$\textbf{40.0}$',r'$\textbf{41.0}$',r'$\textbf{42.0}$'],fontsize=font_size-10)

            inset_axes.set_xlim(-0.03,math.sqrt(3)/2.0+0.03)
            inset_axes.set_xticks(inset_xticks)
            inset_axes.set_xticklabels(inset_xticklabels,fontsize=font_size-12)
            inset_axes.set_ylim(-0.016,0.01)
            inset_axes.set_yticks([-0.015,-0.01,-0.005,0,0.005,0.01])
            inset_axes.set_yticklabels([r'$-\textbf{0.015}$',r'$-\textbf{0.01}$',r'$-\textbf{0.005}$',r'$\textbf{0}$',r'$\textbf{0.005}$',r'$\textbf{0.1}$'],fontsize=font_size-12)

            axbox = main_axes.get_position()
            main_axes.legend(loc=(axbox.x0 + 0.71, axbox.y0 + 0.14),numpoints=1,scatterpoints=1,title=r'$[k_1^{\mbox{ \small tot.}},k_2^{\mbox{\small tot.}}]$',fontsize=font_size-10)
            main_axes.get_legend().get_title().set_fontsize(font_size-10)
        
        else: 
        
            main_axes = plt.axes([0.05, 0.05, 0.9, 0.7]) 	# ([left, bottom, width, height])
            main_axes.get_xaxis().tick_bottom()
            main_axes.get_yaxis().tick_left()
            
            ##  Plot "other" sectors
            main_axes.plot(listOfOffsets,listOfEnergies,'_',ms=5,color=COLOUR_LIST[7],label="Other")
            
            i=1

            for sector in specialSectors:
        
                main_axes.plot(listOfOffsets1[sector[0]+args.x*sector[1]],listOfEnergies1[sector[0]+args.x*sector[1]],'_',ms=8,color=COLOUR_LIST[i])
                main_axes.scatter(listOfOffsets1[sector[0]+args.x*sector[1]],listOfEnergies1[sector[0]+args.x*sector[1]],marker=MARKER_LIST[i+2],s=24,color=COLOUR_LIST[i],label=sector)
                
                i+=1
        
            ##main_axes.set_xlabel(r'$\Delta_{kx}/\sqrt{3}$',fontsize=font_size)
            main_axes.set_xlabel(r'$N_{\phi}$',fontsize=font_size)
            main_axes.set_ylabel(r'$E/E_R$',fontsize=font_size)
            
            if args.formatting==2:
            
                main_axes.set_xlim(-0.003,math.sqrt(3)/6.0+0.003)
                main_axes.set_xticks(xticks)
                main_axes.set_xticklabels(xticklabels,fontsize=font_size-10)
                main_axes.set_ylim([-0.5,9])
                main_axes.set_yticks([0,2,4,6,8])
                main_axes.minorticks_on()
                main_axes.set_yticklabels([r'$\textbf{0}$',r'$\textbf{2.0}$',r'$\textbf{4.0}$',r'$\textbf{6.0}$',r'$\textbf{8.0}$'],fontsize=font_size-10)
                
                main_axes.set_xlabel(r'$N_{\phi}$',fontsize=font_size)
                main_axes.set_ylabel(r'$[E-E_0(N_{\phi})]/E_R$',fontsize=font_size)

            axbox = main_axes.get_position()
            main_axes.legend(loc=1,numpoints=1,scatterpoints=1,title=r'$[k_1^{\mbox{ \small tot.}},k_2^{\mbox{\small tot.}}]$',fontsize=font_size-10)
            main_axes.get_legend().get_title().set_fontsize(font_size-10)
        
        ##plt.grid()

        ##plt.show()
        
        plt.savefig(fig_name,bbox_inches='tight')
        
        print ("\nDONE! find plot at "+fig_name)
        
        energyGapList.sort();
        energyWidthList.sort();
        
        print("ENERGY SPREAD FOUND TO BE "+str(energyGapList[-1]))
        
        print("ENERGY GAP FOUND TO BE "+str(energyWidthList[0]))
        
        ##  Determine the energy spread of the low lying states
        
        ##sortedList = listOfEnergies[:]
        
        ##sortedList.sort()
        
        ##  Determine the maximum difference between the low lying states

        '''
        maxEnergy = 0
        minEnergy = 1000
        
        for energy in sortedList:
            if energy < float(args.gap_energy):
                if energy > maxEnergy:
                    maxEnergy = energy
                if energy < minEnergy:
                    minEnergy = energy
        
        print("ENERGY SPREAD FOUND TO BE "+str(maxEnergy-minEnergy))
        
        ##  Calculate the energy gap
        
        highestLowEnergyState = 0
        lowestHighEnergyState = 1000

        for energy in sortedList:
            if energy < float(args.gap_energy):
                if energy > highestLowEnergyState:
                    highestLowEnergyState = energy
            elif energy < lowestHighEnergyState:
                    lowestHighEnergyState = energy
                
        print("ENERGY GAP FOUND TO BE "+str(lowestHighEnergyState-highestLowEnergyState))
        '''
        
        
################################################################################
######      DETERMINE THE ENERGY GAP      ######################################

    if args.calculate_gap:
    
        ##  Obtain energy values for a given set of system parameters
        
        spectrumData = SpectrumDataArray()

        ##  Read from SQL database

        sqlDb = SQLite(DATA_FOLDER+args.FolderName()+args.SqlName()+'.sql')

        colNames,rows = sqlDb.GetData("OutFileName",SQL_TABLE_NAME,"Completed=1 AND Interaction="+args.interaction+" AND V0="+args.v0+" AND OffsetX="+str(args.offset_x)+" AND OffsetY="+str(args.offset_y))  

        sqlDb.PrintData(colNames,rows)

        for row in rows:

            ##  Read eigenvalue data from the file name indicated
            ##  Loop over all available sectors
            
            for sectorX in range(0,args.x):
            
                for sectorY in range(0,args.y):
                    
                    #if [sectorX,sectorY] in [[0,2],[2,0],[2,2],[0,0]]:
                    
                    ##  Set sector so that the file name is correctly generated
                    ##  by the args.FolderName() function
                    args.sector = [sectorX,sectorY]

                    fileName = DATA_FOLDER+args.FolderName()+row[0]+args.SectorName()+'_eigensystem.dat'
                    
                    ##print("READING DATA FROM FILE: "+fileName)
                    
                    try:
                        with open(fileName) as file:
                            
                            fin = open(fileName)
                    
                            fin.readline()  ##  skip first comment line
                            
                            nbrEigenvalues = int(fin.readline())

                            fin.readline()  ##  skip first comment line
                            
                            ##  Read in data
                            
                            for i in range(0,min(nbrEigenvalues,args.nbr_eigs)):
                                
                                dataPoint = SpectrumDataPoint()
                                
                                dataPoint.energy      = float(fin.readline())
                                dataPoint.interaction = args.interaction
                                dataPoint.v0          = float(args.v0)/0.5       ##  Divide by recoil energy
                                dataPoint.sectors     = [[sectorX,sectorY]]
                                dataPoint.eigLabels   = [i]
                                
                                spectrumData.Append(dataPoint)
                                
                            fin.close()
                            
                    except IOError as e:
                        pass

        if not spectrumData.GetNbrPoints():

            print("ERROR: Data not found in table")
            
            quit()
        
        ##  Put the energy data in ascending order
        
        spectrumData.Sort()

        ##  Obtain the highest lying state below the approximate gap energy 
        ##  and the lowest lying state above the approximate gap energy
        
        highestLowEnergyState = 0
        lowestHighEnergyState = 1000
        
        energyList = spectrumData.GetEnergyList()
        
        for energy in energyList:
            if energy < float(args.gap_energy):
                if energy > highestLowEnergyState:
                    highestLowEnergyState = energy
            elif energy < lowestHighEnergyState:
                    lowestHighEnergyState = energy
                
        print("ENERGY GAP FOUND TO BE "+str(lowestHighEnergyState-highestLowEnergyState))
        
        ##  Write the result to a file

################################################################################
######      PLOT THE TRANSLATIONAL DENSITY-DENSITY FUNCTION ####################
######      sum_{k1,k2} <c^+_{k2-G}c_k2c^+_{k1+G}c_k1>      ####################

    if args.plot_translational_density_density:
    
        ##  Real in data from SQL database
        
        sqlDb = SQLite(DATA_FOLDER+args.FolderName()+args.SqlName()+'.sql')
        
        colNames,rows = sqlDb.GetData("OutFileName",SQL_TABLE_NAME,"Completed=1 AND Interaction="+args.interaction+" AND V0="+args.v0+" AND OffsetX="+args.offset_x+" AND OffsetY="+args.offset_y)

        sqlDb.PrintData(colNames,rows)

        densitydensityData = {}

        specialSectors = [[0,0],[0,args.y/2],[args.x/2,0],[args.x/2,args.y/2]]
        
        if args.formatting==2:
            specialSectors = [[0,0]]
        
        for sector in specialSectors:
            densitydensityData[sector[0]+args.x*sector[1]] = []

        for row in rows:

            ##  Read modified density density function files for each eigenstate

            for sectorX in range(0,args.x):
                for sectorY in range(0,args.y):
                    for eigLabel in range(0,args.nbr_eigs):

                        ##  Set sector so that the file name is correctly generated
                        ##  by the args.FolderName() function
                        args.sector = [sectorX,sectorY]
                    
                        fileName = DATA_FOLDER+args.FolderName()+str(row[0])+args.SectorName()+'_eig_'+str(eigLabel)+'_translational_density_density.dat'

                        print("READING DATA FROM FILE: "+fileName)
                        
                        try:
                            with open(fileName) as file:
                                
                                fin = open(fileName)

                                densitydensityDataInner = []

                                dimX = int(fin.readline()) 
                                dimY = int(fin.readline())

                                for line in fin:
                                    line = line.split()
 
                                    densitydensityDataInner.append(float(line[2])/(args.nbr))

                                fin.close()

                                ##  Convert to np arrays

                                densitydensityDataInner = np.array(densitydensityDataInner)

                                if args.sector in specialSectors:
                                    densitydensityData[args.sector[0]+args.x*args.sector[1]].append(densitydensityDataInner)

                        except IOError as e:
                            print("ERROR READING FILE "+fileName)
                        
                            pass
        
        print(densitydensityData)

        ######      Plot energy vs offset       ####################################################
        
        ##  SET PLOT PARAMETERS HERE

        font_size = 42

        fig_name = 'kx_'+str(args.x)+'_ky_'+str(args.y)+'_n_'+str(args.nbr)+'_translational_density_density.pdf'

        main_axes = plt.axes([0.05, 0.05, 0.9, 0.9]) 	# ([left, bottom, width, height])

        main_axes.get_xaxis().tick_bottom()
        main_axes.get_yaxis().tick_left()
        
        '''
        ##main_axes.set_title(r'$\frac{V_0}{E_R} = '+str(float(args.v0)/0.5)+'$  $g_{2D} = '+str(args.interaction)+'$',fontsize=font_size)
 
        absList = [] 
        
        for i in range(0,args.x):
            for j in range(0,args.y):
                if not (i==0 and j==0):
                    
                    imin = min(i,args.x-i)
                    jmin = min(j,args.y-j)
                
                    k = math.sqrt((float(imin)/args.x)**2+(float(jmin)/args.y)**2)
                    absList.append(k)
                    
        ##main_axes.plot(vals,testCurve,'-')
        '''

        ##  Calculate average values
        
        average = np.zeros(len(specialSectors))
        
        counter = 0
        
        for sector in specialSectors:

            for densityDensity in densitydensityData[sector[0]+args.x*sector[1]]:

                average[counter] += np.average(densityDensity[1:])
            
            counter += 1

        print(average)

        counter = 1
        
        for sector in specialSectors:

            for densityDensity in densitydensityData[sector[0]+args.x*sector[1]]:

                plotY = densityDensity[1:]
                plotY -= average[counter-1]
                plotY /= average[counter-1]

                if args.formatting==2:
                    main_axes.plot(range(1,args.x*args.y),plotY,MARKER_LIST[counter+2]+'-',ms=20,lw=3,color=COLOUR_LIST[counter+2],label=sector)
                else:
                    main_axes.plot(range(1,args.x*args.y),plotY,MARKER_LIST[counter+2]+'-',ms=20,lw=3,color=COLOUR_LIST[counter],label=sector)

                '''
                ##  Sort the data points in preparation for fitting

                sortedData = zip(*sorted(zip(absList,densityDensity[1:])))

                ##main_axes.plot(absList,densityDensity[1:],MARKER_LIST[counter],ms=15,color=COLOUR_LIST[counter],label=sector)

                nSort = len(sortedData[0])

                main_axes.plot(sortedData[0][:nSort],sortedData[1][:nSort],MARKER_LIST[counter],ms=15,color=COLOUR_LIST[counter],label=sector)

                fitData,errors = QuarticFit(sortedData[0][:nSort],sortedData[1][:nSort])
                
                main_axes.plot(sortedData[0][:nSort],fitData,'-',color=COLOUR_LIST[counter])
                '''
            counter += 1
            
        
        main_axes.set_xlabel(r'$\mathbf{q} = q_1 + L_1 \times q_2$',fontsize=font_size)
        main_axes.set_xlim(0.5,args.x*args.y-0.5)
        main_axes.set_xticks(range(0,args.x*args.y,args.y))
        main_axes.set_xticks(range(0,args.x*args.y),minor=True)
        
        main_axes.set_ylabel(r'$\delta\rho_{T}({\mathbf q})$',fontsize=font_size)
        main_axes.tick_params(direction='in',pad=10)
        
        if args.formatting==1:

            main_axes.set_xticklabels([r'$\textbf{0}$',r'$\textbf{6}$',r'$\textbf{12}$',r'$\textbf{18}$'],fontsize=font_size-10)

            main_axes.set_ylim(-0.9,0.7)
            main_axes.set_yticks([-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6])
            main_axes.set_yticklabels([r'$-\textbf{0.8}$',r'$-\textbf{0.6}$',r'$-\textbf{0.4}$',r'$-\textbf{0.2}$',r'$\textbf{0}$',r'$\textbf{0.2}$',r'$\textbf{0.4}$',r'$\textbf{0.6}$'],fontsize=font_size-10)

            main_axes.annotate('$\mathbf q = \,$[0,3]', xy=(3.5,0.25),
                    horizontalalignment='center', verticalalignment='center',fontsize=font_size-10)
                    
            main_axes.annotate('$\mathbf q = \,$[2,0]', xy=(12.5,0.55),
                    horizontalalignment='right', verticalalignment='center',fontsize=font_size-10)
                    
            main_axes.annotate('$\mathbf q = \,$[2,3]', xy=(14.5,0.55),
                    horizontalalignment='left', verticalalignment='center',fontsize=font_size-10)
            
            main_axes.legend(loc=4,numpoints=1,ncol=2,fontsize=font_size-10,title=r'$[k_1^{\mbox{ tot.}},k_2^{\mbox{tot.}}]$')
        
            main_axes.get_legend().get_title().set_fontsize(font_size-10)
        
        if args.formatting==2:

            main_axes.set_xticklabels([r'$\textbf{0}$',r'$\textbf{6}$',r'$\textbf{12}$',r'$\textbf{18}$',r'$\textbf{24}$',r'$\textbf{30}$'],fontsize=font_size-10)

            main_axes.set_ylim(-0.07,0.05)
            main_axes.set_yticks([-0.06,-0.04,-0.02,0,0.02,0.04])
            main_axes.set_yticklabels([r'$-\textbf{0.06}$',r'$-\textbf{0.04}$',r'$-\textbf{0.02}$',r'$\textbf{0}$',r'$\textbf{0.02}$',r'$\textbf{0.04}$'],fontsize=font_size-10)

            main_axes.legend(loc=4,numpoints=1,ncol=2,fontsize=font_size-10,title=r'$[k_1^{\mbox{ tot.}},k_2^{\mbox{tot.}}]$')
        
            main_axes.get_legend().get_title().set_fontsize(font_size-10)
        
        ##plt.grid()

        ##plt.show()
        
        plt.savefig(fig_name,bbox_inches='tight')
        
        print ("\nDONE! find plot at "+fig_name)
        
################################################################################
######      PLOT THE ROTATIONAL DENSITY-DENSITY FUNCTION    ####################
######      sum_{k1,k2} <c^+_{~R60 k2}c_k2c^+_{R60 k1}c_k1> ####################

    if args.plot_rotational_density_density:
    
        ##  Real in data from SQL database
        
        sqlDb = SQLite(DATA_FOLDER+args.FolderName()+args.SqlName()+'.sql')
        
        colNames,rows = sqlDb.GetData("OutFileName",SQL_TABLE_NAME,"Completed=1 AND Interaction="+args.interaction+" AND V0="+args.v0+" AND OffsetX="+args.offset_x+" AND OffsetY="+args.offset_y)

        sqlDb.PrintData(colNames,rows)
        
        densitydensityData = {}

        specialSectors = [[0,0],[0,args.y/2],[args.x/2,0],[args.x/2,args.y/2]]
        
        if args.formatting==2:
            specialSectors = [[0,0]]
        
        for sector in specialSectors:
            densitydensityData[sector[0]+args.x*sector[1]] = []

        for row in rows:

            ##  Read modified density density function files for each eigenstate

            for sectorX in range(0,args.x):
                for sectorY in range(0,args.y):
                    for eigLabel in range(0,args.nbr_eigs):

                        ##  Set sector so that the file name is correctly generated
                        ##  by the args.FolderName() function
                        args.sector = [sectorX,sectorY]
                    
                        fileName = DATA_FOLDER+args.FolderName()+str(row[0])+args.SectorName()+'_eig_'+str(eigLabel)+'_rotational_density_density.dat'

                        print("READING DATA FROM FILE: "+fileName)
                        
                        try:
                            with open(fileName) as file:
                                
                                fin = open(fileName)

                                densitydensityDataInner = []

                                for line in fin:
                                    line = line.split()
 
                                    densitydensityDataInner.append(float(line[1]))

                                fin.close()

                                ##  Convert to np arrays

                                densitydensityDataInner = np.array(densitydensityDataInner)

                                if args.sector in specialSectors:
                                    densitydensityData[args.sector[0]+args.x*args.sector[1]].append(densitydensityDataInner)

                        except IOError as e:
                            print("ERROR READING FILE "+fileName)
                        
                            pass
        
        print(densitydensityData)

        ##  Calculate average values
        
        average = np.zeros(len(specialSectors))
        
        counter = 0
        
        for sector in specialSectors:

            for densityDensity in densitydensityData[sector[0]+args.x*sector[1]]:

                average[counter] += np.average(densityDensity[1:])
            
            counter += 1

        print(average)

        ##  Plot energy vs offset
        
        ##  SET PLOT PARAMETERS HERE

        font_size = 42

        fig_name = 'kx_'+str(args.x)+'_ky_'+str(args.y)+'_n_'+str(args.nbr)+'_rotational_density_density.pdf'
    
        main_axes = plt.axes([0.05, 0.05, 0.9, 0.9]) 	# ([left, bottom, width, height])
        
        main_axes.get_xaxis().tick_bottom()
        main_axes.get_yaxis().tick_left()

        rData = np.array([1,2,3,4,5])
        
        barWidth = 0.1
        barCounter = -2
        
        if args.formatting==2:
            barCounter = 0
        
        counter = 1

        for sector in specialSectors:

            for densityDensity in densitydensityData[sector[0]+args.x*sector[1]]:

                plotY = densityDensity[1:]
                plotY -= average[counter-1]
                plotY /= average[counter-1]

                if args.formatting==2:
                    main_axes.bar(rData+barCounter*barWidth,plotY,color=COLOUR_LIST[counter+2],align='center',alpha=1.0,width=barWidth)
                    main_axes.plot(rData+barCounter*barWidth,plotY,MARKER_LIST[counter+2],ms=20,color=COLOUR_LIST[counter+2],label=sector)
                else:
                    main_axes.bar(rData+barCounter*barWidth,plotY,color=COLOUR_LIST[counter],alpha=1.0,width=barWidth)   
                    main_axes.plot(rData+(barCounter+0.5)*barWidth,plotY,MARKER_LIST[counter+2],ms=20,color=COLOUR_LIST[counter],label=sector)
   
            counter += 1
            barCounter += 1
        
        ##  Draw on average line

        ##main_axes.axhline(y=average,linestyle='--',color='black')

        main_axes.axhline(y=0,linestyle='--',color='black')

        if args.formatting>=1:
        
            main_axes.set_xlabel(r'$m$',fontsize=font_size)
            main_axes.set_xlim([0.5,5.5])
            main_axes.set_xticks(rData)

            main_axes.set_xticklabels([r'$\textbf{1}$',r'$\textbf{2}$',r'$\textbf{3}$',r'$\textbf{4}$',r'$\textbf{5}$'],fontsize=font_size-10)
            ##main_axes.xaxis.set_label_coords(1,-0.05)
            
            main_axes.set_ylabel(r'$\delta\rho_{R}(m)$',fontsize=font_size)
            
            if args.formatting==2:
                main_axes.set_ylim([-0.15,0.1])
                main_axes.set_yticks([-0.15,-0.1,-0.05,0,0.05,0.1])
                main_axes.set_yticklabels([r'$-\textbf{0.15}$',r'$-\textbf{0.1}$',r'$-\textbf{0.05}$',r'$\textbf{0.00}$',r'$\textbf{0.05}$',r'$\textbf{0.1}$'],fontsize=font_size-10)
            
            else:
                main_axes.set_ylim([-0.4,0.2])
                main_axes.set_yticks([-0.4,-0.3,-0.2,-0.1,0,0.1,0.2])
                main_axes.set_yticklabels([r'$-\textbf{0.4}$',r'$-\textbf{0.3}$',r'$-\textbf{0.02}$',r'$-\textbf{0.01}$',r'$\textbf{0.00}$',r'$\textbf{0.01}$',r'$\textbf{0.02}$'],fontsize=font_size-10)
            
            main_axes.tick_params(direction='in',pad=10,axis='y')
            main_axes.tick_params(direction='in',pad=20,axis='x')
            
            ##main_axes.annotate(r'$\langle \rho_{R}(m) \rangle$', xy=(1.25,1.1*average),horizontalalignment='center', verticalalignment='center',fontsize=font_size-10)
            
        main_axes.legend(loc='best',numpoints=1,ncol=2,fontsize=font_size-10,title=r'$[k_1^{\mbox{tot.}},k_2^{\mbox{tot.}}]$')
        
        main_axes.get_legend().get_title().set_fontsize(font_size-10)

        ##plt.show()
        
        plt.savefig(fig_name,bbox_inches='tight')
        
        print ("\nDONE! find plot at "+fig_name)
        
################################################################################
######      Plot band width and band gap vs v0 for a range of   ################
######      theta or epsilon parameters                         ################

    if args.plot_theta_band_width or args.plot_epsilon_band_width:
        
        ##  Real in data from SQL database

        sqlDb = SQLite(DATA_FOLDER+args.add_path+'resultFileKey.sql')
        
        if args.plot_theta_band_width:
            
            colNames,rows = sqlDb.GetData("OutFileName,Theta",SQL_TABLE_NAME,"Completed=1 AND Epsilon="+args.epsilon+" AND OffsetX="+args.offset_x+" AND OffsetY="+args.offset_y+" AND Theta > 0.0 AND Theta <0.6")
            
        elif args.plot_epsilon_band_width:
        
            colNames,rows = sqlDb.GetData("OutFileName,Epsilon",SQL_TABLE_NAME,"Completed=1 AND Theta="+args.theta+" AND OffsetX="+args.offset_x+" AND OffsetY="+args.offset_y+" AND Epsilon > 0.1 AND Epsilon <0.7")
        
        sqlDb.PrintData(colNames,rows)
        
        v0Data      = []
        thetaData   = []
        epsilonData = []
        bandWidthData = []
        bandGapData = []
        ratioData = []

        for row in rows:
            
            if args.plot_theta_band_width:
                thetaData.append(float(row[1]))
            elif args.plot_epsilon_band_width:
                epsilonData.append(float(row[1]))
            
            fileName = DATA_FOLDER+args.add_path+str(row[0])+'.dat'

            print("READING DATA FROM FILE: "+fileName)
            
            V0InnerData = []
            bandWidthInnerData = []
            bandGapInnerData = []
            ratioInnerData = []
            
            try:
                with open(fileName) as file:
                    
                    fin = open(fileName)
                    
                    ##  Skip comment line
                    fin.readline()
                    
                    for line in fin:
                    
                        line = line.split()
                        
                        V0InnerData.append(float(line[0]))
                        bandWidthInnerData.append(float(line[1]))
                        bandGapInnerData.append(float(line[2]))
                        ratioInnerData.append(float(line[2])/float(line[1]))
                        
                    fin.close()
                    
            except IOError as e:
                    print("ERROR READING FILE "+fileName)
                
                    pass       
            
            V0InnerData = np.array(V0InnerData)
            bandWidthInnerData = np.array(bandWidthInnerData)
            bandGapInnerData = np.array(bandGapInnerData)
            ratioInnerData = np.array(ratioInnerData)
            
            v0Data.append(V0InnerData)
            bandWidthData.append(bandWidthInnerData)
            bandGapData.append(bandGapInnerData)
            ratioData.append(ratioInnerData)
            
        ##  Plot band width and band gap vs V0 for a range of theta or epsilon
        
        ##  SET PLOT PARAMETERS HERE

        font_size = 28

        if args.plot_theta_band_width:
            fig_name = 'kx_'+str(args.x)+'_ky_'+str(args.y)+'_band_width_vs_theta.pdf'
        elif args.plot_epsilon_band_width:
            fig_name = 'kx_'+str(args.x)+'_ky_'+str(args.y)+'_band_width_vs_epsilon.pdf'

        main_axes = plt.axes([0.05, 0.05, 0.9, 0.9]) 	# ([left, bottom, width, height])

        main_axes.get_xaxis().tick_bottom()
        main_axes.get_yaxis().tick_left()

        ##print(thetaData)

        if args.plot_theta_band_width:

            for i in range(0,len(thetaData)):
            
                main_axes.plot(v0Data[i],ratioData[i],marker=MARKER_LIST[i+2],ls='-',label=thetaData[i],c=COLOUR_LIST[i],ms=20)
                
                ##main_axes.plot(v0Data[i],bandGapData[i],ls='--',c=COLOUR_LIST[i])
            
            main_axes.legend(loc='best',numpoints=1,title=r'$\theta$',fontsize=font_size-10)
            main_axes.set_ylim(0,8)
            main_axes.set_xlim(0,10)
            xvals=range(0,11)
            yvals=range(0,9)
            
            main_axes.annotate('This work', xy=(2.1, 5.3), xytext=(0.8, 7),fontsize=font_size,
            arrowprops=dict(facecolor='black', shrink=0.05),)
            
        elif args.plot_epsilon_band_width:
        
            for i in range(0,len(epsilonData)):
            
                main_axes.plot(v0Data[i],ratioData[i],label=epsilonData[i],marker=MARKER_LIST[i+2],ls='-',c=COLOUR_LIST[i],ms=20)
                
                ##main_axes.plot(v0Data[i],bandGapData[i],ls='--',c=COLOUR_LIST[i])

            main_axes.legend(loc='best',numpoints=1,title=r'$\epsilon$',fontsize=font_size-10)
            main_axes.set_ylim(0,7)
            main_axes.set_xlim(0,6)
            xvals=range(0,7)
            yvals=range(0,8)
            
            main_axes.annotate('This work', xy=(2.3, 5), xytext=(3, 5.5),fontsize=font_size,
            arrowprops=dict(facecolor='black', shrink=0.05),)

        main_axes.set_xlabel('$V_0/E_R$',fontsize=font_size)
        
        main_axes.set_ylabel('Gap/Width',fontsize=font_size)
        main_axes.get_legend().get_title().set_fontsize(font_size)
        
        main_axes.set_yticks(yvals)
        main_axes.set_yticklabels([r'$\textbf{0.0}$',r'$\textbf{1.0}$',r'$\textbf{2.0}$',r'$\textbf{3.0}$',r'$\textbf{4.0}$',r'$\textbf{5.0}$',r'$\textbf{6.0}$',r'$\textbf{7.0}$',r'$\textbf{8.0}$'][0:len(yvals)],fontsize=font_size-10)
        
        main_axes.set_xticks(xvals)
        main_axes.set_xticklabels([r'$\textbf{0.0}$',r'$\textbf{1.0}$',r'$\textbf{2.0}$',r'$\textbf{3.0}$',r'$\textbf{4.0}$',r'$\textbf{5.0}$',r'$\textbf{6.0}$',r'$\textbf{7.0}$',r'$\textbf{8.0}$',r'$\textbf{9.0}$',r'$\textbf{10.0}$'][0:len(xvals)],fontsize=font_size-10)
        
        main_axes.tick_params(direction='in',pad=10,axis='y')
        main_axes.tick_params(direction='in',pad=10,axis='x')
        
        ##  Grey out the region below G/W=2
        
        main_axes.axhspan(0,2, facecolor='black', alpha=0.1)
        
        plt.savefig(fig_name,bbox_inches='tight')
        
        print ("\nDONE! find plot at "+fig_name)
         
