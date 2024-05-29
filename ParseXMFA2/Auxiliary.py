
from ParseXMFA2.Globals import *
import ParseXMFA2.Globals as Globals

class PseudoRange:
	start : int
	stop : int
	def __init__(self, start, stop):
		self.start = start
		self.stop = stop
	
	def __lt__(self, n):
		return self.start < n
	
	def __gt__(self, n):
		return self.stop > n

class Map:
	_ranges : list[range]
	_filepos : list[int]
	_newline : str

	def __init__(self, newline):
		self._ranges = []
		self._filepos = []
		self._newline = newline

	def __contains__(self, pos):
		for r in self._ranges:
			if pos >= r.start:
				if pos <= r.stop:
					return True
				break
		return False

	def __iadd__(self, newRange : tuple[int,range]):
		for i, r in enumerate(self._ranges):
			# Overlaps should not happen?
			# if newRange.stop > r.start:
			# 	if newRange.start < r.stop:
			# 		# Overlapping
			# 		self._ranges[i] = range(min(newRange.start, r.start), max(newRange.stop, r.stop))
			if newRange[1].start > r.start:
				self._filepos.insert(i, newRange[0])
				self._ranges.insert(i, newRange[1])
	
	def __getitem__(self, pos : int):
		'''Return where in the associated file object that this position corresponds to.'''
		for i, r in enumerate(self._ranges):
			if pos > r.start:
				if pos > r.stop:
					break
				else:
					return self._filepos[i] + pos - r.start + len(self._newline)*(pos - r.start)//LINEWIDTH
		return "-"

	def find(self, pos : int):
		'''Return where in the associated file object that this position corresponds to.'''
		return self[pos]