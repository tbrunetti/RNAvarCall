import argparse
import collections

def checkValues(sample, allSnps):
	import pysam

	tmpCheck = []	
	vcfStore = pysam.TabixFile(sample)

	for row in vcfStore.fetch(parser=pysam.asTuple()):
		tmpCheck.append('\t'.join([str(row[0]) ,str(row[1]), str(row[3]), str(row[4])]))
	try:
		assert collections.Counter(tmpCheck).most_common(1)[0][1] == 1
		print('Confirmed all positions and alleles are unique in your vcf file {}\n'.format(sample))
		allSnps = allSnps + tmpCheck
		return allSnps
	except AssertionError:
		print('There are duplicate positions listed in your vcf file {}.\nGetting unique values only...\n'.format(sample))
		tmp = sorted(set(tmpCheck))
		try:
			assert collections.Counter(tmp).most_common(1)[0][1] == 1
			allSnps = allSnps + tmp
			print('After duplicate removal, confirmed all positions and alleles are unique in your vcf file {} and have been added to the concordance queue\n'.format(sample))
			return allSnps
		except AssertionError:
			print("There was an error removing duplicate positions from your vcf file{}. Please double check your file.\n".format(sample))
			os.Exit(42)


def allConcordant(snps, minVal, maxVal, outPrefix):
	import os

	commonVars = []
	countAllVars = collections.Counter(snps)

	#output = open(outPrefix + "_min_" + str(minVal) + "_max_" + str(maxVal) + "_concordant.txt", 'w')
	bcftoolsOutput = open(outPrefix + "_min_" + str(minVal) + "_max_" + str(maxVal) + "_concordant_bcftools_format.txt", 'w')


	print("A peek at the top 5 most concordant snps: ", countAllVars.most_common(5))
	for key,val in countAllVars.items():
		if ((val >= minVal) and (val <= maxVal)):
			commonVars.append(key)

	with open(outPrefix + "_min_" + str(minVal) + "_max_" + str(maxVal) + "_concordant.txt", 'w') as outFile:
		for snps in commonVars:
			outFile.write('\t'.join(snps.split('\t')) + '\n')
			bcftoolsOutput.write(snps.split('\t')[0] + ":" + snps.split('\t')[1] + '\n')
		outFile.flush()
		bcftoolsOutput.flush()

	bcftoolsOutput.close()

def parseSamples(sampleList, minConc, maxConc):
	allSnps = []
	allSamples = sampleList.split(",")
	
	for sample in allSamples:
		allSnps = checkValues(sample, allSnps)

	if ((minConc == 0) and (maxConc == 0)):
		allConcordant(allSnps, len(allSamples), len(allSamples), args.prefix)
	else:
		allConcordant(allSnps, minConc, maxConc, args.prefix)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Generate files to extract SNPs called across multiple aligners.')
	parser.add_argument('--samples', required=True, type=str, help="Comma-separated list of bgzip and tabix indexed vcf files for which commonly called snps are to be reported.")
	parser.add_argument('--prefix', default="", type=str, help="String prefix for file output name")
	parser.add_argument('--minimumConcordance', default=0, type=int, help="The minimum number of samples/aligners for a snp to be considered reported as concordant.  Default is all samples (set to 0).")
	parser.add_argument('--maximumConcordance', default=0, type=int, help="The maximum number of samples/aligners for a snp to be considered reported as concordant.  Default is all samples (set to 0).")
	args = parser.parse_args()

	parseSamples(args.samples, args.minimumConcordance, args.maximumConcordance)