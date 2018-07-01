import sys
import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt

from utils import *
from Fig3 import *

if __name__ == '__main__':
	if len(sys.argv) > 1:
		path = sys.argv[1]
	else:
		path = Fig5DirPath

	# Load simulation data
	gs = SimulationData(path).scalarStats

	# Query the relevant data
	NactByK = gs.getGroupedData(['MeanDegree_Mean','MeanDegree_StdDev','EXCITED_STATE_CUMUL_Mean', 'EXCITED_STATE_CUMUL_StdDev'], groups1, ['MeanDegree_Mean'], asArray=True)
	NactByL = gs.getGroupedData(['MeanShortestPath_Mean','MeanShortestPath_StdDev','EXCITED_STATE_CUMUL_Mean', 'EXCITED_STATE_CUMUL_StdDev'], groups2, ['MeanShortestPath_Mean'], asArray=True)

	# Plot the figures
	scale = 3
	#############
	# Figure 5C #
	#############
	fig = plt.figure(figsize=(scale*2.10,scale*1.65), dpi=100)
	ax = plt.gca()

	for g, dat in zip(groups1, NactByK):
		ax.errorbar(dat[:,0], dat[:,2], yerr=[dat[:,2]-np.maximum(dat[:,2]-dat[:,3], 1),dat[:,3]], xerr=dat[:,1], capthick=1, label=g.name, **g.style)

	ax.set_yscale('log')
	plt.xlabel('Mean Degree $\langle k \\rangle$')
	plt.ylabel('Propagation Extent ($N_{act}$)')
	plt.xlim([2, 18])
	plt.ylim([5, 1400])

	handles, labels = ax.get_legend_handles_labels()
	handles = [h[0] for h in handles]
	plt.legend(handles, labels, numpoints = 1, fontsize=8)

	#############
	# Figure 5D #
	#############
	fig = plt.figure(figsize=(scale*2.10,scale*1.65), dpi=100)
	ax = plt.gca()

	for g, dat in zip(groups2, NactByL):
		ax.errorbar(dat[:,0], dat[:,2], yerr=[dat[:,2]-np.maximum(dat[:,2]-dat[:,3], 1),dat[:,3]], xerr=dat[:,1], capthick=1, label=g.name, **g.style)

	ax.set_yscale('log')
	plt.xlabel('Mean Shortest Path (L)')
	plt.ylabel('Propagation Extent ($N_{act}$)')
	plt.xlim([2, 16])
	plt.ylim([5, 1400])

	handles, labels = ax.get_legend_handles_labels()
	handles = [h[0] for h in handles]
	plt.legend(handles, labels, numpoints = 1, fontsize=8, loc='lower right')

	plt.show()
