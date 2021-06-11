#!/usr/bin/env python
import sys
import gzip

def splitvcf(file):
	anno_info = []
	snps = []
	indels = []
	if file.endswith('.gz'):
		fh = gzip.open(file,'rb')
	else:
		fh = open(file,'r')
	for lines in fh:
		if lines.startswith('#'):
			anno_info.append(lines)
		else:
			line = lines.strip().split('\t')
			if line[7].startswith('INDEL'):
				indels.append(lines)
			else:
				snps.append(lines)
	return anno_info,snps,indels

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print("python splitSNPindelVCF.py input.vcf")
		exit()
	anno = []
	snps = []
	indels = []
	anno,snps,indels = splitvcf(sys.argv[1])
	with gzip.open(sys.argv[1]+'.snp.gz','wb') as fw_snp:
		fw_snp.write(''.join(anno+snps))
	with gzip.open(sys.argv[1]+'.indel.gz','wb') as fw_indel:	
		fw_indel.write(''.join(anno+indels))
