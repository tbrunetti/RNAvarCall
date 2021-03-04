import argparse
import sys
import os

def updated_test():
	import pandas
	import cyvcf2
	import numpy as np

	def calcAverage(vafList):
		return [sum(vafList)/len(vafList)]

	def makeHeader(samples, bioType):
		header = ['chromosome', 'position', 'ref', 'alt']
		if bioType == "rna":
			nr = [ids + '_RNA_totalReads' for ids in samples]
			header.extend(nr)
			nv = [ids + '_RNA_variantReads' for ids in samples]
			header.extend(nv)
			vaf = [ids + '_RNA_vaf' for ids in samples]
			header.extend(vaf)
			header.extend(['RNA_avgVAF'])
			return header
		elif bioType == "dna":
			header.extend([samples[0] + '_DNA_totalReads', samples[0] + '_DNA_variantReads', samples[0] + '_DNA_vaf'])
			return header


	def getSamples(sampleFile):
		samples = []
		with open(sampleFile, 'r') as sampleNames:
			for line in sampleNames:
				samples.append(line.strip())
		return samples


	vcfFile = cyvcf2.VCF(args.mergedVCF)
	dfRnaFull = []
	dfDnaFull = []

	# warning everytime you iterate through VCF() object, pointer does not reset so must specify VCF() object again every iteration
	rnaOnly = cyvcf2.VCF(args.mergedVCF, samples = getSamples(args.rnaSamples))
	for variant in rnaOnly:
		dfLine = [variant.CHROM, variant.POS,variant.REF, variant.ALT[0]]
		dfLine.extend((variant.format('NR').flatten().tolist()))
		dfLine.extend((variant.format('NV').flatten().tolist()))
		dfLine.extend(np.divide(variant.format('NV'),variant.format('NR')).flatten().tolist())
		dfLine.extend(calcAverage(np.divide(variant.format('NV'),variant.format('NR')).flatten().tolist()))
		dfRnaFull.append(dfLine)

	rnaDf = pandas.DataFrame(dfRnaFull, columns = makeHeader(rnaOnly.samples, "rna"))
	
	dnaOnly = cyvcf2.VCF(args.mergedVCF, samples = getSamples(args.dnaSamples))
	for variant in dnaOnly:
		dfLine = [variant.CHROM, variant.POS, variant.REF, variant.ALT[0]]
		dfLine.extend([variant.format('NR').flatten().tolist()[0]])
		dfLine.extend([variant.format('NV').flatten().tolist()[0]])
		dfLine.extend([np.divide(variant.format('NV'),variant.format('NR')).flatten().tolist()[0]])
		dfDnaFull.append(dfLine)

	dnaDf = pandas.DataFrame(dfDnaFull, columns = makeHeader(dnaOnly.samples, "dna"))

	mergedDf = pandas.DataFrame.merge(dnaDf, rnaDf, how = 'inner', on = ["chromosome", "position", "ref", "alt"])

	mergedDf.to_csv(os.path.join(args.outDir, args.prefix + '.tsv'), sep="\t", index=False)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="")
	parser.add_argument('--mergedVCF', required=True, type=str, help='Full path to .vcf.gz containing all RNAseq variant call information and WES variant call information for a single sample in multisample VCF format.  See docs for example.')
	parser.add_argument('--rnaSamples', required=True, type=str, help='Full path to plain text file containing VCF sample IDs pertaining to RNAseq variant calling algorithms.  One sample name per line.')
	parser.add_argument('--dnaSamples', required=True, type=str, help='Name of WES sample ID')
	parser.add_argument('--outDir', default = os.getcwd(), type=str, help='Full path to output directory to return results.')
	parser.add_argument('--prefix', default = "test", type=str, help='Prefix for output file name')

	args = parser.parse_args()

	updated_test()

