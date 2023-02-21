from __future__ import print_function
import os
import sys
import shutil
import fcntl
import pyTestHarness.test
from   pyTestHarness.colors import NamedColors as pthcolors
from   pyTestHarness.version import getVersion
from   pyTestHarness.utils import py23input

# mpiexec has been observed to set non-blocking I/O, which
#  has been observed to cause problems on OS X with errors like
#  "BlockingIOError: [Errno 35] write could not complete without blocking"
# We use this function to (re)set blocking I/O when launching
def setBlockingIOStdout() :
    fd = sys.stdout
    flags = fcntl.fcntl(fd, fcntl.F_GETFL)
    if flags & os.O_NONBLOCK:
        fcntl.fcntl(fd, fcntl.F_SETFL, flags & ~os.O_NONBLOCK)

class PthTestHarnessLoadException(Exception) :
  pass

def FormattedHourMin(seconds):
  m, s = divmod(int(seconds), 60)
  h, m = divmod(m, 60)
  wt = "%02d:%02d" % (h, m)
  return(wt)

def FormattedHourMinSec(seconds):
  m, s = divmod(seconds, 60)
  h, m = divmod(m, 60)
  wt = "%02d:%02d:%02d" % (h, m,s)
  return(wt)

def pthFormatMPILaunchCommand(mpiLaunch,ranks,corespernode):
  launch = mpiLaunch
  launch = launch.replace("<ranks>",str(ranks))
  launch = launch.replace("<cores>",str(ranks))
  launch = launch.replace("<tasks>",str(ranks))
  launch = launch.replace("<RANKS>",str(ranks))
  if corespernode:
    launch = launch.replace("<corespernode>",str(corespernode))
    launch = launch.replace("<cores_per_node>",str(corespernode))
    launch = launch.replace("<CORESPERNODE>",str(corespernode))
    launch = launch.replace("<CORES_PER_NODE>",str(corespernode))
    launch = launch.replace("<rankspernode>",str(corespernode))
    launch = launch.replace("<ranks_per_node>",str(corespernode))
    launch = launch.replace("<RANKSPERNODE>",str(corespernode))
    launch = launch.replace("<RANKS_PER_NODE>",str(corespernode))
  return(launch)

def generateLaunch_PBS(accountname,queuename,testname,mpiLaunch,executable,ranks,ranks_per_node,walltime,outfile):
  if not ranks:
    print("<generateLaunch_PBS>: Requires the number of MPI-ranks be specified")
  if not walltime:
    print("<generateLaunch_PBS>: Requires the walltime be specified")

  filename = testname + '-pth.pbs'
  file = open(filename,"w")
  file.write("#!/bin/bash\n")

  file.write("# pth: auto-generated pbs file\n")

  if accountname:
    file.write("#PBS -A " + accountname + "\n") # account to charge
  file.write("#PBS -N \"" + testname + "\"" + "\n") # jobname
  file.write("#PBS -o " + testname + ".stdout" + "\n")
  file.write("#PBS -e " + testname + ".stderr" + "\n")

  if queuename:
    file.write("#PBS -q " + queuename + "\n")

  wt = FormattedHourMinSec(walltime*60.0)
  file.write("#PBS -l mppwidth=1024,walltime=" + wt + "\n")

  launch = pthFormatMPILaunchCommand(mpiLaunch,ranks,ranks_per_node)
  file.write(launch + " " + executable + " > " + outfile + "\n\n") # launch command
  file.close()
  return(filename)

def generateLaunch_SLURM(accountname,queuename,testname,constraint,mpiLaunch,execute,ranks,ranks_per_node,walltime,outfile):
  if not ranks:
    print("<generateLaunch_SLURM>: Requires the number of MPI-ranks be specified")
  if not walltime:
    print("<generateLaunch_SLURM>: Requires the walltime be specified")

  filename = testname + '-pth.slurm'
  file = open(filename,"w")
  file.write("#!/bin/bash -l\n")

  file.write("# pth: auto-generated slurm file\n")
  if accountname:
    file.write("#SBATCH --account=" + accountname + "\n") # account to charge
  file.write("#SBATCH --job-name=\"" + testname + "\"" + "\n") # jobname

  file.write("#SBATCH --output=" + testname + ".stdout" + "\n") # jobname.stdout
  file.write("#SBATCH --error=" + testname + ".stderr" + "\n") # jobname.stderr

  if queuename:
    file.write("#SBATCH --partition=" + queuename + "\n")

  file.write("#SBATCH --ntasks=" + str(ranks) + "\n")
  if ranks_per_node:
    file.write("#SBATCH --ntasks-per-node=" + str(ranks_per_node) + "\n")

  wt = FormattedHourMinSec(walltime*60.0)
  file.write("#SBATCH --time=" + wt + "\n")

  if constraint :
    file.write("#SBATCH --constraint=" + constraint + "\n")

  launch = pthFormatMPILaunchCommand(mpiLaunch,ranks,ranks_per_node)
  for e in execute:
    file.write(launch + " " + e + " >> " + outfile + "\n") # launch command
  file.write("\n")
  file.close()
  return(filename)

def generateLaunch_LSF(accountname,queuename,testname,mpiLaunch,execute,ranks,rusage,walltime,outfile):
  if not ranks:
    print("<generateLaunch_LSF>: Requires the number of MPI-ranks be specified")
  if not walltime:
    print("<generateLaunch_LSF>: Requires the walltime be specified")

  filename = testname + '-pth.lsf'
  file = open(filename,"w")
  file.write("#!/bin/sh\n")

  file.write("# pth: auto-generated lsf file\n")

  file.write("#BSUB -J " + testname + "\n") # jobname

  file.write("#BSUB -o " + testname + ".stdout\n") # jobname.stdout
  file.write("#BSUB -e " + testname + ".stderr\n") # jobname.stderr

  if queuename:
    file.write("#BSUB -q " + queuename + "\n")

  file.write("#BSUB -n " + str(ranks) + "\n")

  if rusage:
    file.write("#BSUB -R \'" + rusage + "\'" + "\n")

  wt = FormattedHourMin(walltime*60.0) 
  file.write("#BSUB -W " + wt + "\n")

  launch = pthFormatMPILaunchCommand(mpiLaunch,ranks,None)
  for e in execute:
    file.write(launch + " " + e + " >> " + outfile + "\n") # launch command
  file.write("\n")
  file.close()
  return(filename)

def generateLaunch_LoadLevelerBG(accountname,queuename,testname,execute,total_ranks,machine_ranks_per_node,walltime):
  if not total_ranks:
    print("<generateLaunch_LoadLeveler>: Requires the number of MPI-ranks be specified")
  if not walltime:
    print("<generateLaunch_LoadLeveler>: Requires the walltime be specified")

  print("#!/bin/sh")
  print("# pth: auto-generated llq file")
  print("# @ job_name = " + testname)
  print("# @ job_type = bluegene")
  print("# @ error = $(job_name)_$(jobid).stderr")
  print("# @ output = $(job_name)_$(jobid).stdout")
  print("# @ environment = COPY_ALL;")
  wt = FormattedHourMinSec(walltime*60.0)
  print("# @ wall_clock_limit = " + wt)
  print("# @ notification = never")
  print("# @ class = " + queuename)

  bgsize = math.ceil(total_ranks/machine_ranks_per_node)
  print("# @ bg_size = " + str(bgsize))
  print("# @ bg_connection = TORUS")
  print("# @ queue")

  for e in execute:
    print("runjob -n " + e) # launch command


class Launcher:
  defaultConfFileName = 'pthBatchQueuingSystem.conf'

  @staticmethod
  def writeDefaultDefinition(confFileNameIn=None):
    confFileName = confFileNameIn if confFileNameIn else Launcher.defaultConfFileName
    file = open(confFileName,'w')
    major,minor,patch=getVersion()
    file.write('majorVersion=' + str(major) + '\n')
    file.write('minorVersion=' + str(minor) + '\n')
    file.write('patchVersion=' + str(patch) + '\n')
    file.write('queuingSystemType=none\n' )
    file.write('mpiLaunch=none\n' )
    file.close()

  def __init__(self,confFileName=None):
    self.accountName = []
    self.queueName = []
    self.mpiLaunch = []
    self.queuingSystemType = []
    self.batchConstraint = []
    self.jobSubmissionCommand = []
    self.useBatch = False
    self.verbosity_level = 1
    if confFileName :
      self.confFileName = confFileName
    else :
      self.confFileName = Launcher.defaultConfFileName

    self.setup()

    if self.useBatch :
      if self.mpiLaunch == 'none':
        raise RuntimeError('[pth] If using a queuing system, a valid mpi launch command must be provided')

  def setVerbosityLevel(self,value):
    self.verbosity_level = value

  def setQueueName(self,name):
    self.queueName = name

  def setBatchConstraint(self,argstring):
    self.batchConstraint=argstring

  def setHPCAccountName(self,name):
    self.accountName = name

  def setMPILaunch(self,name):
    self.mpiLaunch = name
    # check for existence of "rank" keyword in the string "name"
    if self.queuingSystemType in ['none', 'None', 'local'] and name != 'none':
      keywordlist = [ '<ranks>', '<cores>', '<tasks>', '<RANKS>' ]
      # check of any of keywordlist[i] appears in name
      valid_launcher = False
      for kword in keywordlist:
        if kword in name:
          valid_launcher = True
          break

      if valid_launcher == False:
        raise RuntimeError('[pth] Your MPI launch command must contain the keyword \"<ranks>\"')

  def setQueueSystemType(self,type):
    if type in ['PBS','pbs']:
      self.queuingSystemType = 'pbs'
      self.jobSubmissionCommand = 'qsub '
      self.useBatch = True
      #print('Recognized PBS queuing system')

    elif type in ['LSF','lsf']:
      self.queuingSystemType = 'lsf'
      self.jobSubmissionCommand = 'bsub < '
      self.useBatch = True
      #print('Recognized LSF queuing system')

    elif type in ['SLURM','slurm']:
      self.queuingSystemType = 'slurm'
      self.jobSubmissionCommand = 'sbatch '
      self.useBatch = True
      #print('Recognized Slurm queuing system')

    elif type in ['LoadLeveler','load_leveler','loadleveler','llq']:
      self.queuingSystemType = 'load_leveler'
      self.jobSubmissionCommand = 'llsubmit '
      self.useBatch = True
      #print('Recognized IBM LoadLeveler queuing system')

    elif type in ['none', 'None', 'local']:
      self.queuingSystemType = 'none'
      self.jobSubmissionCommand = ''
      #print('No queuing system being used')

    else:
      print('Value found: ' + type + ' ...')
      raise RuntimeError('[pth] Unknown or unsupported batch queuing system specified')

  def setQueueName(self,name):
    self.queueName = name

  def view(self):
    print('pth: Batch queueing system configuration [',self.confFileName,']')
    major,minor,patch=getVersion()
    print('  Version:         ',str(major)+'.'+str(minor)+'.'+str(patch))
    print('  Queue system:    ',self.queuingSystemType)
    print('  MPI launcher:    ',self.mpiLaunch)
    if self.useBatch:
      print('  Submit command:  ', self.jobSubmissionCommand)
      if self.accountName:
        print('  Account:       ',self.accountName)
      if self.queueName:
        print('  Queue:         ',self.queueName)
      if self.batchConstraint :
        print('  Constraint:    ', self.batchConstraint)

  def configure(self):
    print('----------------------------------------------------------------')
    print('Creating new configuration file ',self.confFileName)
    v = None
    while not v :
      prompt = '[1] Batch queuing system type <pbs,lsf,slurm,llq,none>: '
      v = py23input(prompt)
      if not v :
        print( 'Required.')
    self.setQueueSystemType(v)

    v = None
    while not v:
      prompt = '[2] MPI launch command with num. procs. flag (required - hit enter for examples): '
      v = py23input(prompt)
      if not v :
        print(' Required. Some example MPI launch commands:')
        print('  No MPI Required           : none')
        print('  Local Machine (mpirun)    : mpirun -np <ranks>')
        print('  Local Machine (mpiexec)   : mpiexec -np <ranks>')
        print('  SLURM w/ aprun            : aprun -B')
        print('  Native SLURM              : srun -n $SLURM_NTASKS')
        print('  LSF (Euler)               : mpirun')
        PETSC_DIR=os.getenv('PETSC_DIR')
        PETSC_ARCH=os.getenv('PETSC_ARCH')
        if PETSC_DIR and PETSC_ARCH:
          print('  Current PETSc MPI wrapper :',os.path.join(PETSC_DIR,PETSC_ARCH,'bin','mpiexec'),'-n <ranks>')
        else :
          print('  Example PETSc MPI wrapper : /users/myname/petsc/arch-xxx/bin/mpiexec -n <ranks>')
        print(' Note that the string \"<ranks>\" must be included if the number of ranks is required at launch.')
        print(' The keyword <ranks> will be replaced by the actual number of MPI ranks (defined by a given test) when the test is launched.')
    self.setMPILaunch(v)

    if self.useBatch == True:
      prompt = '[3] specify a constraint (e.g. "gpu" on Piz Daint) (optional - hit enter if not applicable):'
      v = py23input(prompt)
      self.setBatchConstraint(v)

      prompt = '[4] Account to charge (optional - hit enter if not applicable): '
      v = py23input(prompt)
      self.setHPCAccountName(v)

      prompt = '[5] Name of queue to submit tests to (optional - hit enter if not applicable): '
      v = py23input(prompt)
      self.setQueueName(v)

    self.writeDefinition()
    print('\n')
    print('** If you wish to change the config for your batch system, either')
    print('**  (i) delete the file',self.confFileName,' or')
    print('** (ii) re-run with the command line arg --configure')
    print('----------------------------------------------------------------')

  def setup(self):
    try:
      self.loadDefinition()
    except PthTestHarnessLoadException :
      self.configure()
      self.writeDefinition()

  def writeDefinition(self):
    file = open(self.confFileName,'w')
    major,minor,patch=getVersion()
    file.write('majorVersion=' + str(major) + '\n')
    file.write('minorVersion=' + str(minor) + '\n')
    file.write('patchVersion=' + str(patch) + '\n')
    file.write('queuingSystemType=' + self.queuingSystemType + '\n')
    file.write('mpiLaunch=' + self.mpiLaunch + '\n')
    if self.useBatch == True:
      file.write('accountName=' + self.accountName + '\n')
      file.write('batchConstraint=' + self.batchConstraint + '\n')
      file.write('queueName=' + self.queueName + '\n')
    file.close()

  def loadDefinition(self):
    try:
      majorFile = None
      minorFile = None
      patchFile = None
      file = open(self.confFileName,'r')
      for v in file :
        key,value = v.split('=',1)
        value = value.rstrip()
        if key == 'majorVersion' :
          majorFile = int(value)
        if key == 'minorVersion' :
          minorFile = int(value)
        if key == 'patchVersion' :
          patchFile = int(value)
        if key == 'queuingSystemType' :
          self.setQueueSystemType(value)
        if key == 'mpiLaunch' :
          self.setMPILaunch(value)
        if self.useBatch == True:
          if key == 'batchConstraint' :
            self.setBatchConstraint(value)
          if key == 'queueName' :
            self.setQueueName(value)
          if key == 'accountName' :
            self.setHPCAccountName(value)
      file.close()
    except:
      raise PthTestHarnessLoadException('[pth] You must execute configure(), and or writeDefinition() first')

    # Do not accept conf files if the major.minor version is stale, or if versions are missing
    major,minor,patch = getVersion()
    if majorFile < major or (minorFile < minor and majorFile == major) or \
         majorFile==None or minorFile==None or patchFile==None :
      message = '[pth] Incompatible, outdated configuration file ' + self.confFileName + ' detected. Please delete it and re-run to reconfigure.'
      raise RuntimeError(message)

  def createSubmissionFile(self,testname,commnd,ranks,ranks_per_node,walltime,outfile):
    filename = ''
    if not self.useBatch:
      print('Warning: no submission file creation required')
      return(filename)

    if self.batchConstraint and self.queuingSystemType != 'slurm' :
      message = '[pth] Constraints are only currently supported with SLURM'
      raise RuntimeError(message)

    if self.queuingSystemType == 'pbs':
      filename = generateLaunch_PBS(self.accountName,self.queueName,testname,self.mpiLaunch,commnd,ranks,ranks_per_node,walltime,outfile)

    elif self.queuingSystemType == 'lsf':
      filename = generateLaunch_LSF(self.accountName,self.queueName,testname,self.mpiLaunch,commnd,ranks,None,walltime,outfile)

    elif self.queuingSystemType == 'slurm':
      filename = generateLaunch_SLURM(self.accountName,self.queueName,testname,self.batchConstraint,self.mpiLaunch,commnd,ranks,ranks_per_node,walltime,outfile)

    elif self.queuingSystemType == 'load_leveler':
      raise ValueError('[pth] Unsupported: LoadLeveler needs to be updated')

    print('Created submission file:',filename)
    return(filename)

  def submitJob(self,test):
    setBlockingIOStdout()

    if test.use_sandbox:
      sandboxBack = os.getcwd()
      os.mkdir(test.sandbox_path) # error if  it already exists
      os.chdir(test.sandbox_path)
    test.setVerbosityLevel(self.verbosity_level)
    if not self.useBatch:
      mpiLaunch = self.mpiLaunch

      if self.mpiLaunch == 'none' and test.ranks != 1:
        print('[Failed to launch test \"' + test.name + '\" as test uses > 1 MPI ranks and no MPI launcher was provided]')
      else:
        if self.mpiLaunch == 'none':
          launchCmd = []
          for e in test.execute:
            launchCmd.append( e + " >> " + os.path.join(test.output_path,test.output_file) )
        else:
          launch = pthFormatMPILaunchCommand(mpiLaunch,test.ranks,None)
          launchCmd = []
          for e in test.execute:
            launchCmd.append( launch + ' ' + e + " >> " + os.path.join(test.output_path,test.output_file) )
        lc_len = len(launchCmd)
        lc_count = 0
        for lc in launchCmd:
          lc_count = lc_count + 1
          if self.verbosity_level > 0:
            launch_text = pthcolors.SUBHEADER + '[Executing ' + test.name
            if lc_len > 1 :
              launch_text = launch_text + ' (' + str(lc_count) + '/' + str(lc_len) + ')'
            launch_text = launch_text + ']' + pthcolors.ENDC
            if test.use_sandbox:
              launch_text = launch_text + ' from ' + os.getcwd()
            print(launch_text)
            print(lc)
          test.errno = os.system(lc) >> 8 # TODO: fix this clobbering of errno for multiple tests
          setBlockingIOStdout()
    else:
      outfile = os.path.join(test.output_path,test.output_file)
      launchfile = self.createSubmissionFile(test.name,test.execute,test.ranks,'',test.walltime,outfile)
      launchCmd = self.jobSubmissionCommand + launchfile
      if self.verbosity_level > 0:
        if test.use_sandbox:
          print(pthcolors.SUBHEADER + '[Executing ' + test.name + '] ' + pthcolors.ENDC + 'from ' + os.getcwd())
          print(launchCmd)
        else :
          print(pthcolors.SUBHEADER + '[Executing ' + test.name + ']' + pthcolors.ENDC)
          print(launchCmd)
      os.system(launchCmd)
      setBlockingIOStdout()

    if test.use_sandbox:
        os.chdir(sandboxBack)

  def clean(self,test):
    print('[ -- Removing output for test:',test.name,'-- ]')
    if test.use_sandbox:
      sandboxBack = os.getcwd()
      if not os.path.isdir(test.sandbox_path) :
          os.mkdir(test.sandbox_path)
      os.chdir(test.sandbox_path)
    outfile = os.path.join(test.output_path,test.output_file)
    if os.path.isfile(outfile) :
      os.remove(outfile)
    if test.comparison_file != outfile and os.path.isfile(test.comparison_file) :
      foundInLocalTree = False
      cwd = os.getcwd()
      for (root, dirs, files) in os.walk(cwd) :
        for f in files :
          if os.path.abspath(test.comparison_file) == os.path.abspath(os.path.join(root,f)) :
            foundInLocalTree = True
            break
        if foundInLocalTree :
          break
      if foundInLocalTree :
        os.remove(test.comparison_file)
      else :
        message = "Refusing to remove output file " + test.comparison_file + " since it does not live in the local subtree. If you really wanted to compare with this file, please delete it yourself to proceed"
        raise RuntimeError(message)
    if self.useBatch:
      stderrFile = test.name + '.stderr'
      if os.path.isfile(stderrFile) :
        os.remove(stderrFile)
      stdoutFile = test.name + '.stdout'
      if os.path.isfile(stdoutFile) :
        os.remove(stdoutFile)
      if self.queuingSystemType == 'pbs':
        pbsFile = test.name + '-pth.pbs'
        if os.path.isfile(pbsFile) :
          os.remove(pbsFile)
      elif self.queuingSystemType == 'lsf':
        lsfFile = test.name + '-pth.lsf'
        if os.path.isfile(lsfFile) :
          os.remove(lsfFile)
      elif self.queuingSystemType == 'slurm':
        slurmFile = test.name + '-pth.slurm'
        if os.path.isfile(slurmFile) :
          os.remove(slurmFile)
      elif self.queuingSystemType == 'load_leveler':
        llqFile = test.name + '-pth.llq'
        if os.path.isfile(llqFile) :
          os.remove(llqFile)

    if test.use_sandbox:
      os.chdir(sandboxBack)
      shutil.rmtree(test.sandbox_path) # remove entire subtree

# Deprecated alias for backwards compatibility
pthLaunch = Launcher
