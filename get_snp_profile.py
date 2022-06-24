#!/usr/bin/env python

#Add citation to original manuscript

#!/bin/python

resultDir = './genotypes_bcftools/'


snpData = {}

with open(resultDir+'BARCODE.vcf') as f:
    for l in f:
        if l.startswith('#'):
            continue

        snp, pos, id, ref, alt, qual, filter, info, d, dd = l.split()

        par = {}
        for p in info.split(';'):
            pv = p.split('=')
            par[pv[0]] = pv[1]

        snpData[snp+"_"+pos] = {'pos': pos, 'ref': ref, 'alt': alt, 'qual': qual, 'filter': filter, 'info': par}

# Display count
print('Got data for {} SNPs:'.format(len(snpData)))

# Save/print results
with open(resultDir + '/BARCODE.csv', 'w') as f:
    # Table header
    f.write('snp, pos, coverage, ref_allele, ref_percent, alt_allele, alt_percent, genotype\n')

    # Table data
    for s in sorted(snpData.keys()):
        totalDepth = int(snpData[s]['info']['DP'])
        if int(totalDepth) == 0:
            continue
        depthList  = [int(d) for d in snpData[s]['info']['DP4'].split(',')]
        refDepth   = sum(depthList[0:2])
        altDepth   = sum(depthList[2:4])

        # Estimate the diploid genotype: when the minor allele is more than 10 times weaker than the major allele,
        # we should ignore it for a pure sample?
        if refDepth > altDepth and altDepth/refDepth < 0.1:
            genotype = snpData[s]['ref'] + snpData[s]['ref']
        elif altDepth > refDepth and refDepth/altDepth < 0.1:
            genotype = snpData[s]['alt'] + snpData[s]['alt']
        else:
            genotype = snpData[s]['ref'] + snpData[s]['alt']

        # Two alleles were observed
        f.write(','.join([s, snpData[s]['pos'], str(totalDepth), snpData[s]['ref'], '{:.1f}'.format(refDepth), snpData[s]['alt'], '{:.1f}'.format(altDepth), genotype]) + '\n')
        print('  {} ({})  {} ({:.1f} %)  {} ({:.1f} %)'.format(s, totalDepth, snpData[s]['ref'], refDepth, snpData[s]['alt'], altDepth))
