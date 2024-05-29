
from ParseXMFA2.Globals import *
from ParseXMFA2.ParseXMFA2 import ParseXMFA2

def main(argv=sys.argv):
	import argparse
	
	aparse = argparse.ArgumentParser()
	aparse.add_argument("XMFAfile", metavar="XMFAfile", help="Alignment .XMFA file to load.")
	aparse.add_argument("SNPfile", metavar="SNPfile", help="A .vcf file of the SNPs to be called.")
	aparse.add_argument("-o", "--output", metavar="outFile", help="Name of the .vcf file to be created.")
	aparse.add_argument("-rID", "--referenceID", metavar="referenceID", help="The reference sequence used to create the .XMFA file.")

	args = aparse.parse_args(argv[1:])

	print("Creating Index.")
	parser = ParseXMFA2(filename=args.XMFAfile, referenceID=args.referenceID)

	print(f"Opening {args.SNPfile!r}.")
	snpFile = openVCF(args.SNPfile, "r")
	
	print(f"Creating {args.outFile!r}.")
	outFile = openVCF(args, "w", referenceFile="")
	for entry in snpFile:
		outFile.add(POS=entry.POS, REF=parser[entry.POS])
	outFile.close()
	print("Done!")