import numpy
import argparse
import scipy.stats

def non_param_test(total_count,alt_count,n,all_variants):
    p=0
    selected_vars=numpy.random.randint(low=0,high=len(all_variants),size=n)
    snp_alt_ratio=alt_count/float(total_count)
    for var in selected_vars:
        alt_ratio=all_variants[var][0]/float(all_variants[var][1])
        simulated_bases=numpy.random.choice([0, 1], size=total_count, p=[alt_ratio,1-alt_ratio])
        sim_alts=list(simulated_bases).count(0)
        sim_ref=list(simulated_bases).count(1)
        sim_ratio=sim_alts/float(total_count)

        if snp_alt_ratio < 0.5 and sim_ratio <= snp_alt_ratio:
            p+=1
        elif snp_alt_ratio > 0.5 and sim_ratio >= snp_alt_ratio:
            p+=1    
    p=p/float(n)    
    return(p)

parser = argparse.ArgumentParser("""Add bootstrap and binomial P values to a vcf file, the output vcf is printed to stdout""")
parser.add_argument('--vcf'        ,required = True, type=str, help="vcf file")
parser.add_argument('--ase'        ,required = True, type=str, help="GATK-ASE file (needs to be generated using the input vcf)")
parser.add_argument('-n'        ,type=int,default=10000, help="permutations")
args = parser.parse_args()

ase_list={}
all_variants=[]

first=True
for line in open(args.ase):
    if first:
        first=False
        continue
    content=line.strip().split()
    if not content[0] in ase_list:
        ase_list[content[0]] = {}
    if not content[1] in ase_list[content[0]]:
        ase_list[content[0]][content[1]]={}
  
    p_bin=scipy.stats.binom_test(int(content[6]), n=int(content[7]), p=0.5)
    ase_list[content[0]][content[1]][content[4]]={ "ref_count":int(content[5]),"alt_count":int(content[6]),"tot_count":int(content[7]),"p_bin":p_bin,"non_param":0}
    all_variants.append([int(content[6]),int(content[7])])

all_variants=numpy.array(all_variants)

for chromosome in ase_list:
    for pos in ase_list[chromosome]:
        for alt in ase_list[chromosome][pos]:
            total_count=ase_list[chromosome][pos][alt]["tot_count"]
            alt_count=ase_list[chromosome][pos][alt]["alt_count"]
            ase_list[chromosome][pos][alt]["non_param"]=non_param_test(total_count,alt_count,args.n,all_variants)
 

for line in open(args.vcf):
    if line[0] == "#":
        if not line[1] == "#":
            print ("##INFO=<ID=BootstrapAnn,Number=2,Type=Float,Description=\"BootstrapAnn p-values and GATK-ASEcounter stats (alt_count,total_count,binomial,nonparametric)\">")
        print line.strip()
        continue

    content=line.strip().split()
    if content[0] in ase_list:
        if content[1] in ase_list[content[0]]:
            if content[4] in ase_list[content[0]][content[1]]:
                content[7]+= ";BootstrapAnn={},{},{},{}".format(ase_list[content[0]][content[1]][content[4]]["alt_count"],ase_list[content[0]][content[1]][content[4]]["tot_count"],ase_list[content[0]][content[1]][content[4]]["p_bin"],ase_list[content[0]][content[1]][content[4]]["non_param"])
                print "\t".join(content)
                continue

    print line.strip()
