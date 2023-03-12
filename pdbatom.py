""" 

module containing class PdbAtom 

"""


#COLUMNS        DATA  TYPE    FIELD        DEFINITION
#-------------------------------------------------------------------------------------
# 1 -  6        Record name   "ATOM  "
# 7 - 11        Integer       serial       Atom  serial number.
#13 - 16        Atom          name         Atom name.
#17             Character     altLoc       Alternate location indicator.
#18 - 20        Residue name  resName      Residue name.
#22             Character     chainID      Chain identifier.
#23 - 26        Integer       resSeq       Residue sequence number.
#27             AChar         iCode        Code for insertion of residues.
#31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
#39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
#47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
#55 - 60        Real(6.2)     occupancy    Occupancy.
#61 - 66        Real(6.2)     tempFactor   Temperature  factor.
#77 - 78        LString(2)    element      Element symbol, right-justified.
#79 - 80        LString(2)    charge       Charge  on the atom.

class PdbAtom:
	"""
	Lightweight class for reading ATOM lines of pdb files
	
	:param atomSerialNumber: pdb atomSerialNumber 
	:type pdb_in: `int`
	
	:param atomName: pdb atomName 
	:type atomName: `str`
	
	:param alternateLocationIndicator: pdb alternateLocationIndicator 
	:type alternateLocationIndicator: `str`
	
	:param residueName: pdb residueName 
	:type residueName: `str`
	
	:param chainIdentifier: pdb chainIdentifier 
	:type chainIdentifier: `str`
	
	:param residueSequenceNumber: pdb residueSequenceNumber 
	:type chainIdentifier: `int`
	
	:param insertionCode: pdb insertionCode 
	:type insertionCode: `str`
	
	:param xOrthogonalCoordinate: pdb xOrthogonalCoordinate 
	:type xOrthogonalCoordinate: `float`
	
	:param yOrthogonalCoordinate: pdb yOrthogonalCoordinate 
	:type yOrthogonalCoordinate: `float`
	
	:param zOrthogonalCoordinate: pdb zOrthogonalCoordinate 
	:type zOrthogonalCoordinate: `float`
	
	:param occupancy: pdb occupancy 
	:type occupancy: `float`
	
	:param temperatureFactor: pdb temperatureFactor 
	:type temperatureFactor: `float`
	
	:param temperatureFactor: pdb temperatureFactor 
	:type temperatureFactor: `float`
	
	:param elementSymbol: pdb elementSymbol 
	:type elementSymbol: `str`

	:param charge: pdb charge 
	:type charge: `str`
	
	:Example:
	
	.. code-block:: python3
		
		# read all atom lines from pdb file
		atoms = []
		with open('mystruct.pdb') as f:
			pdb_lines = f.readlines()
		for pdb_line in pdb_lines:
			if pdb_line.startswith('ATOM'):
				atom = pdbatom.PdbAtom()
				atom.parse(pdb_line)
				atoms.append(atom)
		
	"""
	def __init__(self):
		"""
		Constructor, creates empty PdbAtom
		"""
		self.atomSerialNumber = int()
		self.atomName		  = str()
		self.alternateLocationIndicator = str()
		self.residueName = str()
		self.chainIdentifier = str()
		self.residueSequenceNumber = int()
		self.insertionCode = str()
		self.xOrthogonalCoordinate =  float()
		self.yOrthogonalCoordinate =  float()
		self.zOrthogonalCoordinate =  float()
		self.occupancy	= float()
		self.temperatureFactor = float()
		self.elementSymbol = str()
		self.charge = str()
	def parse(self,line):
		"""
		Parse PdbAtom values from pdb ATOM line
		
		:param line: pdb ATOM line 
		:type line: `str`
	
		"""
		self.atomSerialNumber = int(line[6:11])
		self.atomName = line[12:16].strip()
		self.alternateLocationIndicator = line[16].strip()
		self.residueName = line[17:20].strip()
		self.chainIdentifier = line[21].strip()
		self.residueSequenceNumber = int(line[22:26])
		self.insertionCode = line[26].strip()
		self.xOrthogonalCoordinate =  float(line[30:38])
		self.yOrthogonalCoordinate =  float(line[38:46])
		self.zOrthogonalCoordinate =  float(line[46:54])
		self.occupancy	= float(line[54:60])
		self.temperatureFactor = float(line[60:66])
		self.elementSymbol = line[76:78].strip()
		self.charge = line[78:80].strip()
	def print(self):
		"""
		Report on PdbAtom class member values
		"""
		print("atomSerialNumber = ",self.atomSerialNumber)
		print("atomName = ",self.atomName)
		print("alternateLocationIndicator = ",self.alternateLocationIndicator)
		print("residueName = ",self.residueName)
		print("chainIdentifier = ",self.chainIdentifier)
		print("residueSequenceNumber = ",self.residueSequenceNumber)
		print("insertionCode = ",self.insertionCode)
		print("xOrthogonalCoordinate = ",self.xOrthogonalCoordinate)
		print("yOrthogonalCoordinate = ",self.yOrthogonalCoordinate)
		print("zOrthogonalCoordinate = ",self.zOrthogonalCoordinate)
		print("occupancy = ",self.occupancy)
		print("temperatureFactor = ",self.temperatureFactor)
		print("elementSymbol = ",self.elementSymbol)
		print("charge = ", self.charge)
	def format(self):
		"""
		Format pdb ATOM line from PdbAtom 
		
		:return: pdb ATOM line
		:rtype: `str`
		"""
		atomSerialNumber_str = str(self.atomSerialNumber).rjust(5)
		atomName_str = self.atomName.strip().rjust(5)
		altLoc_str = self.alternateLocationIndicator.strip()
		if altLoc_str=="": altLoc_str=" "
		residuename_str = self.residueName.strip().rjust(3)
		chainIdentifier_str = self.chainIdentifier.rjust(2)
		if chainIdentifier_str=="": chainIdentifier_str="  "
		residueSequenceNumber_str = str(self.residueSequenceNumber).rjust(4)
		xOrthogonalCoordinate_str =  "{:.3f}".format(self.xOrthogonalCoordinate).rjust(8)
		yOrthogonalCoordinate_str =  "{:.3f}".format(self.yOrthogonalCoordinate).rjust(8)
		zOrthogonalCoordinate_str =  "{:.3f}".format(self.zOrthogonalCoordinate).rjust(8)
		occupancy_str =  "{:.2f}".format(self.occupancy).rjust(6)
		temperatureFactor_str = "{:.2f}".format(self.temperatureFactor).rjust(6)
		elementSymbol_str = self.elementSymbol.strip().rjust(2)
		charge_str = self.charge.strip().rjust(2)	
		pdbline = "ATOM  "+atomSerialNumber_str+atomName_str+altLoc_str+residuename_str+chainIdentifier_str+residueSequenceNumber_str
		pdbline += "    "+xOrthogonalCoordinate_str+yOrthogonalCoordinate_str+zOrthogonalCoordinate_str
		pdbline += occupancy_str+temperatureFactor_str+"          "+elementSymbol_str+""+charge_str
		return pdbline

		
		
		
		
		

