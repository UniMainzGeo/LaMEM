import sys
import os, fnmatch

from optparse import OptionParser
from paraview import servermanager


def locate(pattern, root=os.curdir):
	#Locate all files matching supplied filename pattern in and below supplied root directory.
	for path, dirs, files in os.walk(os.path.abspath(root)):
		for filename in fnmatch.filter(files, pattern):
			yield os.path.join(path, filename)


def convertVTSasciiToVTSbinary(inputfilename,outputfilename):
	print 'converting',inputfilename, '(ascii) ==>>',outputfilename, '(binary)'
	# create reader for legacy VTK files
	reader = servermanager.sources.XMLStructuredGridReader(FileName=inputfilename)

	# create VTU writer and connect it to the reader
	writer = servermanager.writers.XMLStructuredGridWriter(Input=reader,FileName=outputfilename)

	# Trigger execution of pipeline
	writer.UpdatePipeline()


def convert_single(inputfilename,rm_file):
	INPUTFILE = inputfilename
	# create output name
	brokenline = INPUTFILE.partition('.vts')
	OUTPUTFILE = 'tmp.vts'
	# convert
	convertVTSasciiToVTSbinary(INPUTFILE,OUTPUTFILE)

	# swap names
	renamedINPUTFILE = INPUTFILE
	brokenline = renamedINPUTFILE.partition('.vts')
	renamedINPUTFILE = brokenline[0] + '_ascii_vts'

	cmd = 'mv ' + INPUTFILE + ' ' + renamedINPUTFILE
#	print cmd
	os.system(cmd)
	cmd = 'mv ' + OUTPUTFILE + ' ' + INPUTFILE
#	print cmd
	os.system(cmd)

	if bool(rm_file) == True:
		cmd = 'rm ' + renamedINPUTFILE
#		print cmd
		os.system(cmd)
	


parser = OptionParser()

parser.add_option( "-f", "--filename", dest="singleinput",
    help="Single input file to covert. Default path is current working directory.", metavar="FILE")

parser.add_option( "-p", "--prefix", dest="multipleinput",
    help="Prefix of input files to covert. Default path is current working directory.", metavar="FILE")

parser.add_option( "-r", "--remove", dest="rm_flag", default=False,
    help="Remove ascii files after conversion. By default files are not removed.", metavar="FILE")

parser.add_option( "-d", "--dir", dest="dirinput",
    help="Directory to covert. Default path is current working directory.", metavar="FILE")


(options,args) = parser.parse_args()


#help = 'StructuredGrid-XML(Paraview,ascii) => StructuredGrid-XML(Paraview,binary) converter'

if "None" not in str(options.singleinput):
		if "None" not in str(options.multipleinput):
			print 'You cannot specify both -p PREFIX or -f FILENAME. Run with -h to see all options'
			sys.exit(0)

if "None" not in str(options.dirinput):
		if "None" in str(options.multipleinput):
			print 'You must specify a file prefix via -p PREFIX with -d DIRNAME. Run with -h to see all options'
			sys.exit(0)

if "None" in str(options.singleinput):
		if "None" in str(options.multipleinput):
			print 'You must specify either -p PREFIX or -f FILENAME. Run with -h to see all options'
			sys.exit(0)



rm_converted_files = bool(options.rm_flag)

# create a built-in connection
if not servermanager.ActiveConnection:
    connection = servermanager.Connect()


if "None" not in str(options.singleinput):
	INPUTFILE = str(options.singleinput)
	convert_single(INPUTFILE,rm_converted_files)

if "None" not in str(options.multipleinput):
	target = options.multipleinput + '*.vts'
	print target
	# search for files
	for commonfile in locate(target):
		INPUTFILE = commonfile
		convert_single(INPUTFILE,rm_converted_files)


if "None" not in str(options.dirinput):
	if "None" not in str(options.multipleinput):
		targetdir = options.dirinput + '*'

		# search for dir
		for commondir in locate(target):

			# search for files
			for commonfile in locate(target,commondir):
				INPUTFILE = commonfile
				convert_single(INPUTFILE,rm_converted_files)

	else:
		print 'You must specify a prefix to convert'
		sys.exit(0)




