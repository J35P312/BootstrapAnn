# BootstrapAnn
A tool for computing ASE P values based on ASEReadCounter output files. BootstrapAnn accepts a vcf file, aswell as a GATK-ASEReadCounter output file.
The results of the GATK-ASEReadCounter analysis is added to the output vcf, together with allele specific expression P values for each analysed snp.
The P values are computed using a two sided binomial test, as well as a non-parametric exact test.

# Running

    python BootstrapAnn.py --vcf <input.vcf> --ase <ASEReadCounter.tab > <output.vcf>

The output vcf is an annotated version of the input vcf.

# Dependencies
    

BootstrapAnn requires python 2.7, as well as these two dependencies

    scipy

    numpy

These two dependencies may be installed using pip:

    pip install numpy

    pip install scipy
