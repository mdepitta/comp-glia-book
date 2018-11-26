import subprocess
import multiprocessing
import os 
import argparse

# Parameters
nbRepeats = 20
nbPerJob = 5

Fig3Seed = 600000
Fig5Seed = 300000

# Paths
dirPath = os.path.dirname(os.path.realpath(__file__))

ExecPath = os.path.join(dirPath, 'src', 'AstroSim')

Fig3Dir = 'Fig3Data'
Fig5Dir = 'Fig5Data'

Fig3ParamPath = os.path.join(dirPath, 'parameters', 'Fig3Params')
Fig5ParamPath = os.path.join(dirPath, 'parameters', 'Fig5Params')

Fig3DirPath = os.path.join(dirPath, 'data', Fig3Dir)
Fig5DirPath = os.path.join(dirPath, 'data', Fig5Dir)

# Generate the job commands associated with a given figure
def GenerateJobCommands(paramPath, dirName, addParams = [], seed = 12345):
	# Load default simulation arguments
	params = [ExecPath]
	with open(paramPath) as f:
		for line in f:
			params += line.split()

	# Modify some arguments
	params += ['-Path', dirPath, '-SubDir', dirName, '-repeat', str(nbRepeats), '-SplitGrid', str(nbPerJob)]
	params += addParams

	# Prepare simulations
	print('Preparing simulations for ' + dirName + '...')
	subprocess.run(params, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
	print('Done.')

	# Generate job commands
	launch = []
	with open(os.path.join(dirPath, 'data', dirName, 'LaunchGrids.job')) as f:
		for i, line in enumerate(f):
			launch.append([ExecPath, '-simLoad'] + line.split() + 
				['-modelLoad', '-SubDir', dirName, '-Path', dirPath, '-Sim' , '-seed', str(seed + i)])
	return launch

if __name__ == '__main__':
	# Argument parsing
	parser = argparse.ArgumentParser()
	parser.add_argument('-n', '--nbCores', type = int, default=os.cpu_count(), 
		help = 'Number of cores to be used')
	parser.add_argument('-f', '--figures', type = int, nargs = '+', 
		choices = [3, 5, 6], help = 'Simulations to be run (figure number)')
	args = parser.parse_args()

	allJobs = []
	# Figure 3 and 6
	if (args.figures is None or 3 in args.figures or 6 in args.figures):
		allJobs += GenerateJobCommands(Fig3ParamPath, Fig3Dir, seed = Fig3Seed)

	# Figure 5
	if (args.figures is None or 5 in args.figures):
		allJobs += GenerateJobCommands(Fig5ParamPath, Fig5Dir, seed = Fig5Seed)

	# Simulations
	pool = multiprocessing.Pool(args.nbCores)
	print('Running simulations ({0} jobs on {1} cores)...'.format(len(allJobs), args.nbCores))
	for i, _ in enumerate(pool.imap_unordered(subprocess.run, allJobs)):
		print('Job {0} done. ({1:.2%})'.format(i,(i+1)/len(allJobs)))

