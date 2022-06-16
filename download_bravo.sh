#!/bin/bash
curl 'https://bravo.sph.umich.edu/freeze5/hg38/download/all' -H 'Accept-Encoding: gzip, deflate, br' -H 'omitted' --compressed > bravo-dbsnp-all.vcf.gz

bcftools view -Ob bravo-dbsnp-all.vcf.gz > bravo-dbsnp-all.bcf

bcftools index bravo-dbsnp-all.bcf
