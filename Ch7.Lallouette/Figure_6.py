import sys
import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt

from utils import *
from Fig3 import *
from RunSimulations import *

rmax = 20

groups3 = [
	Group('Regular', fmt='-s', color='#a7f906').addFilter(
		'NetworkConstructStrat','SpatialRegularBuildingStrategy'),
	Group('Link Radius', fmt='-o', color='#18ca35').addFilter(
		'NetworkConstructStrat','SpatialConnectionRadiusBuildingStrategy'),
	Group('Shortcut $p_s=0$', fmt='-d', color='#cd87de').addFilter(
		'NetworkConstructStrat','SmallWorldBuildingStrategy').addFilter(
		'SmallWorldProba',0),
	Group('Shortcut $p_s=0.3$', fmt='-d', color='#660080').addFilter(
		'NetworkConstructStrat','SmallWorldBuildingStrategy').addFilter(
		'SmallWorldProba',0.3),
	Group('Scale Free $r_c=2 \mathrm{\mu m}$', fmt='-^', color='#a5b2f4').addFilter(
		'NetworkConstructStrat','SpatialScaleFreeBuildingStrategy').addFilter(
		'SpatialScaleFreeInteractionRange',0.000002),
	Group('Scale Free $r_c=25 \mathrm{\mu m}$', fmt='-^', color='#1835ca').addFilter(
		'NetworkConstructStrat','SpatialScaleFreeBuildingStrategy').addFilter(
		'SpatialScaleFreeInteractionRange',0.000025),
	Group("Erdös-Rényi", fmt='-v', color='#ff8080').addFilter(
		'NetworkConstructStrat','ErdosRenyiRandomBuildingStrategy')]

groupInds = [
	('Regular $\langle k \\rangle = 3$',0,0,{}), 
	('Regular $\langle k \\rangle = 4$',0,1,{'color':'#44aa00'}),
	('Link Radius $\langle k \\rangle = 6$',1,1,{}), 
	("Erdös-Rényi $\langle k \\rangle = 5$",6,0,{}),
	('Shortcut $p_s=0$, $\langle k \\rangle = 6$',2,0,{})]

# Iteratively compute the fraction of activated astrocytes according to the shell model.
def ShellModelCompute(N, E, W, omega, S, eta, delta, A, B):
	epsilon = 10**-5

	Nhatr = Whatr = lambda rho, Wr: 2*Wr*rho*(1-rho)
	getk = lambda Nr, Erm1, Er, Wr: (Erm1+Er+2*Wr)/max(epsilon,Nr)
	actFunct = lambda Psitot, k: 0.5*(np.tanh((Psitot-(A*k+B))/max(epsilon,delta))+1)

	PsirStim = lambda r: S*sum(N[0:r+1])**(-eta)
	PsirOut = lambda rho, Nr, Nrp1, Nrm1, Er, Erm1, Wr: rho*Nr / max(epsilon,Nrp1 + (Erm1+Whatr(rho,Wr))*(Nrm1+Nhatr(rho,Wr))*(Er+omega*Nrp1)/max(epsilon, Er*(Erm1+Whatr(rho,Wr)+omega*(Nrm1+Nhatr(rho,Wr)))))
	Psirtot = lambda r, rho, Nr, Nrp1, Nrm1, Er, Erm1, Wr: PsirOut(rho, Nr, Nrp1, Nrm1, Er, Erm1, Wr) + PsirStim(r)

	rhoVals = np.zeros(max(len(N), rmax))
	for r, Nr in enumerate(N):
		Erm1 = E[r-1] if r > 0 else 0
		Erm2 = E[r-2] if r > 1 else 0
		Nrm1 = N[r-1] if r > 0 else 0
		Nrm2 = N[r-2] if r > 1 else 0
		Wrm1 = W[r-1] if r > 0 else 0
		psitot = Psirtot(r, rhoVals[r-1] if r > 0 else 0, Nrm1, N[r], Nrm2, Erm1, Erm2, Wrm1)
		rhoVals[r] = actFunct(psitot, getk(N[r], Erm1, E[r], W[r]))

		if rhoVals[r]*N[r] < 1:
			break

	return sum([N[r]*rhoVals[r] for r in range(len(N))]), rhoVals

actDensMean = ['ActivDensityInShell_' + str(i+1).zfill(5) + '_Mean' for i in range(rmax)]
actDensStdDev = ['ActivDensityInShell_' + str(i+1).zfill(5) + '_StdDev' for i in range(rmax)]

def plotShellModelFigures(sd, omega, S, eta, delta, A, B):
	scalDat, fullDat = sd.getGroupedFullData(
		['MeanDegree_Mean','TotalNbCellsActivated_Mean'] + actDensMean + actDensStdDev, 
		['r', 'ShellNodes', 'ShellIntra', 'ShellOut'], 
		groups3, scalOrderBy=['MeanDegree_Mean'], asArray=True, replaceNanBy0 = True)

	# Compute Shell model propagation and get data
	meanK = []
	meanNact = []
	rhoExpMean = []
	rhoExpStd = []
	meanNsim = []
	rhoSim = []
	# For each group
	for gsDat, gfDat in zip(scalDat, fullDat):
		meanK.append([])
		meanNact.append([])
		rhoExpMean.append([])
		rhoExpStd.append([])
		meanNsim.append([])
		rhoSim.append([])
		# For each mean degree
		for ksDat, kfDat in zip(gsDat, gfDat):
			meanK[-1].append(ksDat[0,0])
			meanNact[-1].append(ksDat[0,1])
			rhoExpMean[-1].append([1] + list(ksDat[0, 2:rmax+1]))
			rhoExpStd[-1].append([0] + list(ksDat[0, rmax+2:2*rmax+1]))
			meanNsim[-1].append(0)
			rhoSim[-1].append(np.zeros(rmax))
			# For each run
			for rData in kfDat:
				# Adding the r=0 shell
				Nvals = [1] + list(rData[:,1])
				Wvals = [0] + list(rData[:,2])
				Evals = [Nvals[1]] + list(rData[:,3])
				# Running the simulation
				Ns, rho = ShellModelCompute(Nvals, Evals, Wvals, omega, S, eta, delta, A, B)
				meanNsim[-1][-1] += Ns
				rhoSim[-1][-1] += rho[0:rmax]
			meanNsim[-1][-1] /= len(kfDat)
			rhoSim[-1][-1] /= len(kfDat)

	# Plot Figures
	scale = 3
	#############
	# Figure 6E #
	#############
	fig = plt.figure(figsize=(scale*2.10,scale*1.65), dpi=100)
	ax = plt.gca()

	for gI, g in enumerate(groups3):
		style = copy.copy(g.style)
		del style['fmt']
		ax.plot(meanK[gI], meanNsim[gI], '--', **style)
		ax.plot(meanK[gI], meanNact[gI], g.style['fmt'].replace('-',''), label = g.name, **style)
		
	ax.set_yscale('log')
	plt.xlabel('Mean Degree $\langle k \\rangle$')
	plt.ylabel('Propagation Extent ($N_{act}$)')
	plt.xlim([2, 18])
	plt.ylim([5, 1400])

	handles, labels = ax.get_legend_handles_labels()
	plt.legend(handles, labels, numpoints = 1, fontsize=8)

	#############
	# Figure 6F #
	#############
	fig = plt.figure(figsize=(scale*2.10,scale*1.65), dpi=100)
	ax = plt.gca()

	rVals = list(range(0,rmax))
	for name, gI, kI, stl in groupInds:
		style = copy.copy({**groups3[gI].style, **stl})
		del style['fmt']
		ax.errorbar(rVals, rhoExpMean[gI][kI], yerr=[rhoExpStd[gI][kI],rhoExpStd[gI][kI]], capthick=1, label=name, fmt=groups3[gI].style['fmt'].replace('-',''), **style)
		ax.plot(rVals, rhoSim[gI][kI], '--', **style)

	plt.xlabel('Shell number r')
	plt.ylabel('Fraction of activated astrocytes $\\rho_r$')
	plt.xlim([0, 15])
	plt.ylim([0, 1.2])

	handles, labels = ax.get_legend_handles_labels()
	handles = [h[0] for h in handles]
	plt.legend(handles, labels, numpoints = 1, fontsize=8, ncol=2, loc='upper right')

	plt.show(block=True)


if __name__ == '__main__':
	if len(sys.argv) > 1:
		path = sys.argv[1]
	else:
		path = Fig3DirPath
	sd = SimulationData(path, ['NetDimFullData.dat'])

	# Parameters
	omega = 1.78527
	S     = 2169.32
	eta   = 2.62
	delta = 0.05388
	A     = 2.622e-15
	B     = 0.35061

	plotShellModelFigures(sd, omega, S, eta, delta, A, B)

