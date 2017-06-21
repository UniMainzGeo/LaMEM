
import os,sys
import argparse

import pyTestHarness.unittest as pth
from   pyTestHarness.colors import pthNamedColors as bcolors

isPython2 = False
isPython3 = False
if sys.version_info[0] == 2:
  isPython2 = True;
if sys.version_info[0] == 3:
  isPython3 = True;

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
  
  file.write("# pyTH: auto-generated pbs file\n")

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

def generateLaunch_SLURM(accountname,queuename,testname,mpiLaunch,executable,ranks,ranks_per_node,walltime,outfile):
  if not ranks:
    print("<generateLaunch_SLURM>: Requires the number of MPI-ranks be specified")
  if not walltime:
    print("<generateLaunch_SLURM>: Requires the walltime be specified")
  
  filename = testname + '-pth.slurm'
  file = open(filename,"w")
  file.write("#!/bin/bash -l\n")
  
  file.write("# pyTH: auto-generated slurm file\n")
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
  
  launch = pthFormatMPILaunchCommand(mpiLaunch,ranks,ranks_per_node)
  file.write(launch + " " + executable + " > " + outfile + "\n\n") # launch command
  file.close()
  return(filename)

def generateLaunch_LSF(accountname,queuename,testname,mpiLaunch,executable,ranks,rusage,walltime,outfile):
  if not ranks:
    print("<generateLaunch_LSF>: Requires the number of MPI-ranks be specified")
  if not walltime:
    print("<generateLaunch_LSF>: Requires the walltime be specified")
  
  filename = testname + '-pth.lsf'
  file = open(filename,"w")
  file.write("#!/bin/sh\n")
  
  file.write("# pyTH: auto-generated lsf file\n")

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
  file.write(launch + " " + executable + " > " + outfile + "\n\n") # launch command
  file.close()
  return(filename)

def generateLaunch_LoadLevelerBG(accountname,queuename,testname,executable,total_ranks,machine_ranks_per_node,walltime):
  if not total_ranks:
    print("<generateLaunch_LoadLeveler>: Requires the number of MPI-ranks be specified")
  if not walltime:
    print("<generateLaunch_LoadLeveler>: Requires the walltime be specified")
  
  print("#!/bin/sh")
  print("# pyTH: auto-generated llq file")
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
  
  print("runjob -n " + executable) # launch command


def performTestSuite_local(self,registered_tests):
  launcher = self
  
  print('')
  self.view()

  for test in registered_tests:
    print('-- Executing test: ' + test.name + ' --')
    launcher.submitJob(test)
    test.verifyOutput()
  
  counter = 0
  for test in registered_tests:
    test.report('summary')
    if test.passed == False:
      counter = counter + 1
  if counter > 0:
    print('-- Unit test error report --')
    for test in registered_tests:
      test.report('log')

  print('-- Unit test report summary --')
  for test in registered_tests:
    test.report('summary')
  if counter > 0:
    print('  ' + str(counter) + ' / ' + str(len(registered_tests)) + ' tests ' + bcolors.FAIL + 'failed' + bcolors.ENDC)
  else:
    print('----------------------')
    print(bcolors.OKGREEN + ' All tests passed' + bcolors.ENDC)
  


def performTestSuite_execute(self,registered_tests):
  launcher = self

  print('')
  if self.verbosity_level > 0:
    self.view()

  for test in registered_tests:
    print('[-- Executing test: ' + test.name + ' --]')
    launcher.submitJob(test)


def performTestSuite_verify(self,registered_tests):
  launcher = self
  
  print('')
  for test in registered_tests:
    print('[-- Verifying test: ' + test.name + ' --]')
    if self.mpiLaunch == 'none' and test.ranks != 1:
      print('[Skipping verification for test \"' + test.name + '\" as test uses > 1 MPI ranks and no MPI launcher was provided]')
    else:
      test.verifyOutput()
  
  print('')
  counter = 0
  for test in registered_tests:
    if test.passed == False:
      counter = counter + 1
  if counter > 0:
    print('')
    print('[--------- Unit test error report ----------------------]')
    for test in registered_tests:
      test.report('log')

  print('[--------- Unit test summary ----------------------]')
  counter = 0
  for test in registered_tests:
    if self.mpiLaunch == 'none' and test.ranks != 1:
      print('  ['+test.name+']  skipped as ranks > 1 and no MPI launcher provided')
    else:
      test.report('summary')
    if test.passed == False:
      counter = counter + 1
  if counter > 0:
    print(bcolors.FAIL + '          ********************' + bcolors.ENDC)
    print(bcolors.FAIL +  ' [status] ' + str(counter) + ' of ' + str(len(registered_tests)) + ' tests FAILED' + bcolors.ENDC)
    print(bcolors.FAIL + '          ********************' + bcolors.ENDC)
  else:
    print(bcolors.OKGREEN + '          ****************' + bcolors.ENDC)
    print(bcolors.OKGREEN + ' [status] All tests passed' + bcolors.ENDC)
    print(bcolors.OKGREEN + '          ****************' + bcolors.ENDC)
  


class pthLaunch:

  def __init__(self):
    self.accountName = []
    self.queueName = []
    self.mpiLaunch = []
    self.queuingSystemType = []
    self.jobSubmissionCommand = []
    self.use_batch = False
    self.output_path = ''
    self.verbosity_level = 1
    
    parser = argparse.ArgumentParser(description='Python Test Harness.')
    parser.add_argument('-e', '--execute', help='Perform test execution', required=False, action='store_true')
    parser.add_argument('-v', '--verify', help='Perform test verification (and not execution)', required=False, action='store_true')
    parser.add_argument('-c', '--configure', help='Configure queuing system information', required=False, action='store_true')
    parser.add_argument('-t', '--test', help='List of test names', required=False)
    parser.add_argument('-o', '--output_path', help='Directory to write stdout into', required=False)
    parser.add_argument('-p', '--purge_output', help='Delete generated output', required=False, action='store_true')
    parser.add_argument('-f', '--error_on_test_failure', help='Return exit code of 1 if any test failed', required=False, action='store_true')
    parser.add_argument('-d', '--configure_default', help='Write default queuing system config file (no mpi, no queuing system)', required=False, action='store_true')
    self.args = parser.parse_args()

    if self.args.configure:
      self.configure()
      sys.exit(0)
    elif self.args.configure_default:
      self.writeDefaultDefinition()
      sys.exit(0)
    else:
      self.setup()

    if self.use_batch == True:
      if self.mpiLaunch == 'none':
        raise RuntimeError('[pth] If using a queuing system, a valid mpi launch command must be provided')

#  def addIgnoreKeywords(self,d):
#    self.ignoreKeywords.append(d)

  def setVerbosityLevel(self,value):
    self.verbosity_level = value

  def setQueueName(self,name):
    self.queueName = name


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
      self.use_batch = True
      #print('Recognized PBS queuing system')

    elif type in ['LSF','lsf']:
      self.queuingSystemType = 'lsf'
      self.jobSubmissionCommand = 'bsub < '
      self.use_batch = True
      #print('Recognized LSF queuing system')

    elif type in ['SLURM','slurm']:
      self.queuingSystemType = 'slurm'
      self.jobSubmissionCommand = 'sbatch '
      self.use_batch = True
      #print('Recognized Slurm queuing system')

    elif type in ['LoadLeveler','load_leveler','loadleveler','llq']:
      self.queuingSystemType = 'load_leveler'
      self.jobSubmissionCommand = 'llsubmit '
      self.use_batch = True
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
    print('pth: Batch queueing system configuration [pthBatchQueingSystem.conf]')
    print('  Queue system:    ',self.queuingSystemType)
    print('  MPI launcher:    ',self.mpiLaunch)
    if self.use_batch:
      print('  Submit command:', self.jobSubmissionCommand)
      if self.accountName:
        print('  Account:       ',self.accountName)
      if self.queueName:
        print('  Queue:         ',self.queueName)


  def configure(self):
    print('----------------------------------------------------------------')
    print('Creating new pthBatchQueuingSystem.conf file')
    if isPython2:
      v = raw_input('[1] Batch queuing system type <pbs,lsf,slurm,llq,none>: ')
    if isPython3:
      v = input('[1] Batch queuing system type <pbs,lsf,slurm,llq,none>: ')
    if not v:
      raise ValueError('[pth] You must specify the type of queuing system')
    self.setQueueSystemType(v)

    v = None
    while not v:
      if isPython2:
        v = raw_input('[2] MPI launch command with num. procs. flag (required - hit enter for examples): ')
      if isPython3:
        v = input('[2] MPI launch command with num. procs. flag (required - hit enter for examples): ')
      if not v :
        print(' Required. Some example MPI launch commands:')
        print('  none','(if your tests do not use MPI)')
        print('  mpirun -np <ranks>','(local machine)')
        print('  mpiexec -n <ranks>','(local machine)')
        print('  aprun -B','(slurm with aprun)')
        print('  srun -n $SLURM_NTASKS','(native slurm)')
        print('  /users/myname/petsc/bin/petscmpiexec -n <ranks>','(typical PETSc MPI wrapper)')
        print(' Note that the string \"<ranks>\" must be included in your launch command.')
        print(' The keyword <ranks> will be replaced by the actual number of MPI ranks (defined by a given test) when the test is launched.')
    self.setMPILaunch(v)

    if self.use_batch == True:

      if isPython2:
        v = raw_input('[3] Account to charge (optional - hit enter if not applicable): ')
      if isPython3:
        v = input('[3] Account to charge (optional - hit enter if not applicable): ')
      self.setHPCAccountName(v)

      if isPython2:
        v = raw_input('[4] Name of queue to submit tests to (optional - hit enter if not applicable): ')
      if isPython3:
        v = input('[4] Name of queue to submit tests to (optional - hit enter if not applicable): ')
      self.setQueueName(v)
    
    self.writeDefinition()
    print('\n')
    print('** If you wish to change the config for your batch system, either')
    print('**   (i) delete the file pthBatchQueuingSystem.conf, or')
    print('**  (ii) re-run pth2.configure()')
    print('** (iii) re-run with the command line arg --configure')
    print('----------------------------------------------------------------')


  def setup(self):
    try:
      self.loadDefinition()
    
    except:
      self.configure()
      self.writeDefinition()

  def writeDefaultDefinition(self):
    file = open('pthBatchQueuingSystem.conf','w')
    file.write( 'none' + '\n' )
    file.write( 'none' + '\n' )
    file.close()

  def writeDefinition(self):

    file = open('pthBatchQueuingSystem.conf','w')
    #    file.write('queuingSystemType = ' + self.queuingSystemType + '\n' )
    #    file.write('accountName = ' + self.accountName + '\n' )
    #    file.write('queueName = ' + self.queueName + '\n' )
    #    file.write('mpiLaunch = ' + self.mpiLaunch + '\n' )
    file.write( self.queuingSystemType + '\n' )
    file.write( self.mpiLaunch + '\n' )
    if self.use_batch == True:
      file.write( self.queueName + '\n' )
      file.write( self.accountName + '\n' )
    file.close()


  def loadDefinition(self):
    try:
      file = open('pthBatchQueuingSystem.conf','r')

      v = file.readline()
      self.setQueueSystemType(v.rstrip())

      v = file.readline()
      self.setMPILaunch(v.rstrip())

      if self.use_batch == True:
        v = file.readline()
        self.setQueueName(v.rstrip())

        v = file.readline()
        self.setHPCAccountName(v.rstrip())

      file.close()
    except:
      raise RuntimeError('[pth] You must execute configure(), and or writeDefinition() first')


  def createSubmissionFile(self,testname,commnd,ranks,ranks_per_node,walltime,outfile):
    filename = ''
    if not self.use_batch:
      print('Warning: no submission file creation required')
      return(filename)
    
    if self.queuingSystemType == 'pbs':
      filename = generateLaunch_PBS(self.accountName,self.queueName,testname,self.mpiLaunch,commnd,ranks,ranks_per_node,walltime,outfile)

    elif self.queuingSystemType == 'lsf':
      filename = generateLaunch_LSF(self.accountName,self.queueName,testname,self.mpiLaunch,commnd,ranks,None,walltime,outfile)

    elif self.queuingSystemType == 'slurm':
      filename = generateLaunch_SLURM(self.accountName,self.queueName,testname,self.mpiLaunch,commnd,ranks,ranks_per_node,walltime,outfile)

    elif self.queuingSystemType == 'load_leveler':
      raise ValueError('[pth] Unsupported: LoadLeveler needs to be updated')

    print('Created submission file:',filename)
    return(filename)


  def submitJob(self,unittest):
    unittest.setVerbosityLevel(self.verbosity_level)
    if not self.use_batch:
      mpiLaunch = self.mpiLaunch
    
      if self.mpiLaunch == 'none' and unittest.ranks != 1:
        print('[Failed to launch test \"' + unittest.name + '\" as test uses > 1 MPI ranks and no MPI launcher was provided]')
      else:
        if self.mpiLaunch == 'none':
          launchCmd = unittest.execute + " > " + os.path.join(unittest.output_path,unittest.output_file)
        else:
          launch = pthFormatMPILaunchCommand(mpiLaunch,unittest.ranks,None)
          launchCmd = launch + ' ' + unittest.execute + " > " + os.path.join(unittest.output_path,unittest.output_file)
        if self.verbosity_level > 0:
          print('[Executing]',launchCmd)
        unittest.errno = os.system(launchCmd) >> 8
    else:
      outfile = os.path.join(unittest.output_path,unittest.output_file)
      launchfile = self.createSubmissionFile(unittest.name,unittest.execute,unittest.ranks,'',unittest.walltime,outfile)
      launchCmd = self.jobSubmissionCommand + launchfile
      if self.verbosity_level > 0:
        print('[Executing]',launchCmd)
      os.system(launchCmd)


  def executeTestSuite(self,registered_tests):
    if self.args.output_path:
      for test in registered_tests:
        test.setOutputPath(self.args.output_path)
    
    # Don't execute if we are verifying a batch run
    if not self.use_batch and not self.args.verify:
      performTestSuite_execute(self,registered_tests)


  def verifyTestSuite(self,registered_tests):
  
    # Verify, unless we are running with a batch system and are not in verify(-only) mode
    if not self.use_batch or self.args.verify :
      performTestSuite_verify(self,registered_tests)

  def clean(self,registered_tests):
    if self.use_batch:
      for test in registered_tests:
        print('<launch> rm -f ' + test.name + '.stdout')
        print('<launch> rm -f ' + test.name + '.stderr')
        if self.queuingSystemType == 'pbs':
          print('<launch> rm -f ' + test.name + '-pth.pbs')
        elif self.queuingSystemType == 'lsf':
          print('<launch> rm -f ' + test.name + '-pth.lsf')
        elif self.queuingSystemType == 'slurm':
          print('<launch> rm -f ' + test.name + '-pth.slurm')
        elif self.queuingSystemType == 'load_leveler':
          print('<launch> rm -f ' + test.name + '-pth.llq')



# < end class >
