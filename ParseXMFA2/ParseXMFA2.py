'''
ParseXMFA is a script that parse xmfa files and extract all SNPs
	The script allows an option using a flanking score that limits
	SNPs near edges of a block to have a impact on classifying decisions.
'''
from ParseXMFA2.Globals import *
import ParseXMFA2.Globals as Globals
from ParseXMFA2.Auxiliary import Map, PseudoRange

class ParseXMFA2:
	entryHeader : re.Pattern = re.compile(b"^[>] (?P<seqID>[0-9]+):(?P<start>[0-9]+)-(?P<stop>[0-9]+) - .*", re.MULTILINE)
	reference : int
	filename : str
	_mmap : mmap.mmap
	coverage : list[Map]

	def __init__(self, filename : str=None, mode : str="r", referenceID : int=1):
		if filename is not None:
			self.open(filename, mode, referenceID)
			self.mapBlocks()
	
	def __contains__(self, item):
		return any(item in map for map in self.coverage)
	
	def __getitem__(self, key : list[int]|int):
		offset = self.index(key)
		if type(offset) is list:
			return map(self._mmap.__getitem__, offset)
		else:
			return self._mmap[offset]
				
	def open(self, filename, mode : str, referenceID : int):
		self.filename = filename
		self.referenceID = referenceID
		self.file = open(self.filename, "r+b") # Makes the file writeable, but we do not intend to edit it
		self.file.readline()
		self.newlines = self.file.newlines
		self.coverage = {}
		self.file.seek(0)

		if mode == "rw":
			self._mmap = mmap.mmap(self.file.fileno(), 0)
		elif mode == "r":
			self._mmap = mmap.mmap(self.file.fileno(), 0, prot=mmap.PROT_WRITE)
		elif mode == "w":
			self._mmap = mmap.mmap(self.file.fileno(), 0, prot=mmap.PROT_READ)
		else:
			raise ValueError(f"{mode!r} is not a valid mode for opening files.")
	
	def mapBlocks(self):
		'''Can handle multiple alignments'''
		block = {}
		for m in self.entryHeader.finditer():
			ntRange = range(int(m.group("start")), int(m.group("stop")))
			filepos = m.end()+len(self.newlines)
			if int(m.group("seqID")) in block:
				if self.referenceID not in block:
					raise ValueError("XMFA file must have blocks that match the given referenceID sequence.")
				for i in block:
					if i != self.referenceID:
						if i not in self.coverage: self.coverage[i] = Map(self.newlines)
						self.coverage[i] += block[i][0], block[self.referenceID][1]
						
				self.block = {}
			else:
				block[int(m.group("seqID"))] = (filepos, ntRange)
	
	def index(self, pos : int, *mPos : int):
		'''Takes more than 0 integer positions in the sequence covered by the XMFA and returns the byte position(s) in
		the XMFA file that correspond to the given position in the aligned sequence. If only one position is given the
		return type is int. If more than one positions are ggiven, a list of byte positions is returned.'''
		if len(mPos) == 0:
			return [map[pos] for map in self.coverage]
		else:
			return [[map[p] for map in self.coverage] for p in [pos]+mPos]
