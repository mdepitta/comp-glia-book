import sys

import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt

from utils import *
from RunSimulations import *

cColor  = '#ff2a2a'
ucColor = '#333333'
iColor  = '#aa0000'

groups1 = [
	Group("Erdös-Rényi", fmt='-v', color=ucColor).addFilter(
		'NetworkConstructStrat','ErdosRenyiRandomBuildingStrategy'),
	Group('Scale Free $r_c=25 \mathrm{\mu m}$', fmt='-^', color=ucColor).addFilter(
		'NetworkConstructStrat','SpatialScaleFreeBuildingStrategy').addFilter(
		'SpatialScaleFreeInteractionRange',0.000025),
	Group('Shortcut $p_s=0.3$', fmt='-d', color=ucColor).addFilter(
		'NetworkConstructStrat','SmallWorldBuildingStrategy').addFilter(
		'SmallWorldProba',0.3),
	Group('Scale Free $r_c=2 \mathrm{\mu m}$', fmt='-^', color=cColor).addFilter(
		'NetworkConstructStrat','SpatialScaleFreeBuildingStrategy').addFilter(
		'SpatialScaleFreeInteractionRange',0.000002),
	Group('Shortcut $p_s=0$', fmt='-d', color=cColor).addFilter(
		'NetworkConstructStrat','SmallWorldBuildingStrategy').addFilter(
		'SmallWorldProba',0),
	Group('Link Radius', fmt='-o', color=cColor).addFilter(
		'NetworkConstructStrat','SpatialConnectionRadiusBuildingStrategy'),
	Group('Regular', fmt='-s', color=cColor).addFilter(
		'NetworkConstructStrat','SpatialRegularBuildingStrategy')]

groups2 = [
	Group("Erdös-Rényi", fmt='-v', color=ucColor).addFilter(
		'NetworkConstructStrat','ErdosRenyiRandomBuildingStrategy'),
	Group('Shortcut $\langle k \\rangle=6$', fmt='-d', color=ucColor).addFilter(
		'NetworkConstructStrat','SmallWorldBuildingStrategy').addFilter(
		'SmallWorldNeighbDist',1),
	Group('Link Radius', fmt='-o', color=cColor).addFilter(
		'NetworkConstructStrat','SpatialConnectionRadiusBuildingStrategy'),
	Group('Regular', fmt='-s', color=cColor).addFilter(
		'NetworkConstructStrat','SpatialRegularBuildingStrategy'),
	Group('Scale Free $\langle k \\rangle=4$', fmt='-^', color=iColor).addFilter(
		'NetworkConstructStrat','SpatialScaleFreeBuildingStrategy').addFilter(
		'SpatialScaleFreeNewLinks',2)]

if __name__ == '__main__':
	if len(sys.argv) > 1:
		path = sys.argv[1]
	else:
		path = Fig3DirPath

	# Load simulation data
	gs = SimulationData(path).scalarStats

	# Query the relevant data
	NactByK = gs.getGroupedData(['MeanDegree_Mean','MeanDegree_StdDev','TotalNbCellsActivated_Mean', 'TotalNbCellsActivated_StdDev'], groups1, ['MeanDegree_Mean'], asArray=True)
	NactByL = gs.getGroupedData(['MeanShortestPath_Mean','MeanShortestPath_StdDev','TotalNbCellsActivated_Mean', 'TotalNbCellsActivated_StdDev'], groups2, ['MeanShortestPath_Mean'], asArray=True)

	# Plot the figures
	scale = 3
	#############
	# Figure 3C #
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
	# Figure 3D #
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
