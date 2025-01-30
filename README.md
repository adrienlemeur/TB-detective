<h1 align="center"> TB-detective </h1>

### <a name="quickstart"></a>installation
you need cyvf2
> https://brentp.github.io/cyvcf2/
that's it

### usage
usage: TB-detective [-h] -i I [I ...] [-lin LIN] [-ab AB] [-cf CF]

TB-detective is a script to identify the lineage, sublineage and antibiotic
resistance of a Mycobacterium tuberculosis sample from a VCF annotated with
snpEff. It was written in python with cyvcf2.

options:
  -h, --help    show this help message and exit
  -i I [I ...]  a single sample VCF aligned on H37Rv genome
  -lin LIN      a tab separated table with 1-based SNP position, the ALT
                nucleotide and the associated lineage
  -ab AB        a tab separated table with genes, mutations (snpEff prot. or
                nuc. mutation annotation), associated antibiotic resistance
                and confidence threshold
  -cf CF        threshold level of confidence (only AMR associated variants
                with the same or higher confidence will be considered)

Written by Adrien Le Meur, v.1.2
