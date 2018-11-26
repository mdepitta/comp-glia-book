import re
import os
from operator import itemgetter
import numpy as np
import copy

class Group:
	def __init__(self, name, **style):
		self.name = name
		self.filters = {}
		self.style = style
	def addFilter(self, fn, v):
		self.filters[fn] = v
		return self
	def setStyle(self, **styles):
		self.style = styles

class DataFile:
	def __init__(self, path):
		self.header = {}
		self.data = []
		if path:
			with open(path) as f:
				self.header = {s:i for i,s in enumerate(f.readline().split())}
				self.data = []
				for line in f:
					self.data.append([])
					for v in line.split():
						self.data[-1].append(self.castValue(v))

	def castValue(self, v):
		p = re.compile('-?[\d\.]+')
		if v == 'NULL' or v=='NaN':
			return None
		elif p.match(v):
			return float(v)
		else:
			return v
	
	def getFilteredData(self, filters):
		return [row for row in self.data if all([row[self.header[fn]] == v for fn, v in filters.items()])]

	def getNonNullFieds(self, filters = {}):
		res = self.getFilteredData(filters)	
		return [fn for fn, i in self.header.items() if all([row[i] is not None for row in res])]

	def getData(self, fieldNames, filters = {}, orderBy=[], asArray = False, replaceNanBy0 = False):
		func = (lambda x:x) if not replaceNanBy0 else (lambda x: (0 if x is None else x))
		inds = [(self.header[fn] if fn in self.header else None) for fn in fieldNames]
		res = self.getFilteredData(filters)
		if len(orderBy) > 0:
			res = sorted(res, key = itemgetter(*[self.header[ob] for ob in orderBy]))
		if len(fieldNames) > 1:
			res = [[(func(row[i]) if i is not None else (0 if replaceNanBy0 else i)) for i in inds] for row in res]
		else:
			res = [func(row[inds[0]]) for row in res]
		return np.array(res) if asArray else res

	def getGroupedData(self, fieldNames, groups, orderBy=[], filters = {}, asArray=False, replaceNanBy0 = False):
		return [self.getData(fieldNames, {**g.filters, **filters}, orderBy, asArray = asArray, replaceNanBy0 = replaceNanBy0) for g in groups]

	def getGroups(self, fieldNames, filters={}, orderBy=[]):
		diffVals = []
		for valComb in [tuple(v for v in row) for row in self.getData(fieldNames, filters, orderBy)]:
			if valComb not in diffVals:
				diffVals.append(valComb)
		groups = []
		for line in diffVals:
			groups.append(Group('defaultName'))
			for fn, v in zip(fieldNames, line):
				groups[-1].addFilter(fn, v)
		return groups

	def merge(self, other):
		selfKeys = set(self.header.keys())
		otherKeys = set(other.header.keys())
		newKeys = otherKeys.difference(selfKeys)
		for fn in newKeys:
			self.header[fn] = len(self.header)
		tmpData = []
		for row in self.data:
			tmpData.append(row + [None]*len(newKeys))
		for row in other.data:
			tmpData.append([None]*len(self.header))
			for fn, i in other.header.items():
				tmpData[-1][self.header[fn]] = row[i]
		self.data = tmpData
		return self

	def addColumn(self, fieldName, valueFunc):
		if fieldName not in self.header:
			self.header[fieldName] = len(self.data[0]) if len(self.data) > 0 else 0
			for i in range(len(self.data)):
				self.data[i].append(self.castValue(valueFunc(self.header, self.data[i], i)))
		else:
			return False

class SimulationData:
	def __init__(self, path, perRunFname = []):
		p = re.compile('Sim_[\d]+')
		p2 = re.compile('(Run)([\d]+)')
		p3 = re.compile('(.+)_([^_]+)')
		self.scalarStats = DataFile('')
		self.fullStats = DataFile('')
		for fname in os.listdir(path):
			if p.match(fname) and os.path.isdir(os.path.join(path, fname)):
				totName = os.path.join(path, fname, 'GridScalarStatistics.dat')
				if os.path.isfile(totName):
					self.scalarStats.merge(DataFile(totName))
		for root, dirs, files in os.walk(path):
			for f in files:
				if f in perRunFname:
					df = DataFile(os.path.join(root, f))
					ok = False
					#
					rt, sf = os.path.split(root)
					while len(rt) > 0 and not p.match(sf):
						m = p2.match(sf) if p2.match(sf) else p3.match(sf)
						if m:
							df.addColumn(m.group(1), lambda h, d, i: m.group(2))
						rt, sf = os.path.split(rt)
					self.fullStats.merge(df)

	def getGroupedFullData(self, scalFNames, fullFNames, groups, scalOrderBy=[], fullOrderBy=[], filters={}, asArray=False, replaceNanBy0 = False):
		scalRes = []
		fullRes = []
		for g in groups:
			filt = {**filters, **g.filters}
			groupParams = list(set(self.scalarStats.getNonNullFieds(filt)).intersection(set(self.fullStats.getNonNullFieds(filt))))
			subgroups = self.scalarStats.getGroups(groupParams, filt, scalOrderBy)
			scalRes.append([])
			fullRes.append([])
			for sb in subgroups:
				scalRes[-1].append(self.scalarStats.getData(scalFNames, {**filt, **sb.filters}, scalOrderBy, asArray, replaceNanBy0))
				fullRes[-1].append([])
				runNbs = set(self.fullStats.getData(['Run'], {**filt, **sb.filters}))
				for r in runNbs:
					fullRes[-1][-1].append(self.fullStats.getData(fullFNames, {**filt, **sb.filters, 'Run':r}, fullOrderBy, asArray, replaceNanBy0))

		return scalRes, fullRes
