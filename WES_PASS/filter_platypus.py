############################################################
# filter_platypus.py Daniel Nichol 06/03/2017.
#
# Filters the vcf file produced by platypus when used to genotype
# the sample .bam files from the Mutect2 calls. Records those variants
# that arise in difficult/bad regions (e.g. centrosomes)
# 
# 
# The filtering criteria are:
#   1) The FILTER tag for the call is in filterList (below)
#   2) The variant is not a known GL variant, or aligned to the decoy genome.
#   3) A genotype is called for all samples (i.e. no sample assigned './.')
#   4) The genotype phred scores (GQ) are >=10 for all samples
#   5) The number of reads covering the variant site is >=10 in all samples. 
#   6) The normal sample has no reads exhibing the variant
#   7) At least one tumour samples has >=3 reads for the variant.
#
# Note: We assume that the patient sample list has the normal 
# sample name as the first line. This sample list does NOT need to
# match the ordering of the columns of the vcf file.
#
############################################################
import sys

#####################################################################
# Filter list - the list of FILTER entries that are acceptable.
# Amend this list to include/exclude specific FILTER values.
#####################################################################
filterList =   ['PASS', 'alleleBias', 'Q20', 'Q20;alleleBias','QD', 'Q20;QD', 'QD;alleleBias',
                'Q20;QD;alleleBias', 'SC', 'SC;Q20', 'SC;alleleBias', 'SC;Q20;alleleBias', 'SC;QD',
                'SC;Q20;QD', 'SC;QD;alleleBias', 'SC;Q20;QD;alleleBias', 'HapScore', 'Q20;HapScore',
                'HapScore;alleleBias', 'Q20;HapScore;alleleBias', 'QD;HapScore', 'Q20;HapScore;QD',
                'QD;HapScore;alleleBias', 'Q20;HapScore;QD;alleleBias', 'SC;HapScore', 'SC;Q20;HapScore',
                'SC;HapScore;alleleBias', 'SC;Q20;HapScore;alleleBias', 'SC;HapScore;QD', 'SC;Q20;HapScore;QD',
                'SC;HapScore;QD;alleleBias', 'SC;Q20;HapScore;QD;alleleBias']


#####################################################################
# Parse the arguments
#####################################################################

if len(sys.argv) < 4 or len(sys.argv)>=6:
    print "Improper usage is python filter_platypus.py <path>/patient.platypus.vcf <path2ref> samples.txt [optional]<path>/patient.mutect.sort.vcf"
    exit()
elif len(sys.argv) == 4:
    platypus_file = sys.argv[1]
    path2Ref = sys.argv[2]
    patient_file = sys.argv[3]
    venn = False
else: 
    platypus_file = sys.argv[1]
    path2Ref = sys.argv[2]
    patient_file = sys.argv[3]
    mutect_file = sys.argv[4]
    venn = True

#####################################################################
# Create the list of regions to ignore (from the look up table)
#####################################################################
badRegs = []
with open(str(path2Ref)+'/hg19_gap_contigs_hgTables.txt','r') as cent:
    line = cent.readline() #Skip the one header line
    line = cent.readline()
    while line:
        regInfo = line.split('\t')
        #(chrom, start, end) format:
        badRegs.append((str(regInfo[0][3:]), int(regInfo[1]), int(regInfo[2])))
        line = cent.readline()

#####################################################################
# Identify the normal sample index, and number of samples
# ASSUMPTION: The normal is the first sample in the sampleList file.
#####################################################################
with open(patient_file, 'r') as sampleList:
    lines = sampleList.readlines()
    normal = lines[0]
    lines = sorted(lines)
    normIx = lines.index(normal)
    samples = len(lines)

#####################################################################
# Filter the platypus file, record those variants that arise in "bad regions"
#####################################################################
with open(platypus_file[0:-4]+"_PASS.vcf",'w') as platypus_pass:
    with open(platypus_file[0:-4]+"_badRegions.txt",'w') as badRegions:
        with open(platypus_file,'r') as platcalls:
            #Mutation type and position info, used for the venn diagram
            mutPlat, snpsPlat, indelsPlat = [], [], []

            #Copy the header of the platypus file.
            line = ''
            while line[0:6]!="#CHROM":
                line = platcalls.readline()
                platypus_pass.write(line)
            
            line = platcalls.readline()
            while line:
                record = line.split('\t')
                #Split the (colon (:) separated) individual sample information for each sample.
                for k in range(len(record[:-samples]), len(record)):
                    record[k] = record[k].split(':')
                    # If multiple alts exist, we take NR, NV to be the maximum amongst them
                    if ',' in record[k][4]:
                        record[k][4] = max(record[k][4].split(','))
                        record[k][5] = max(record[k][5].split(','))

                allSamples = record[-samples:]
                tumSamples = record[-samples:-samples + normIx] + record[-samples+normIx+1:]
                normSample = record[-samples+normIx]
                #Filters: FILTER passed and not decoy/germline
                if (record[6] in filterList) and (record[0]!='hs37d5') and (record[0][0:2]!='GL'):
                    GQsPass = all([float(el[3]) >= 10 for el in allSamples])
                    NRsPass = all([float(el[4]) >= 10 for el in allSamples])

                    # GQs (genotype phred scores) all >= 10, and
                    # NRs (number of reads at site) all >= 10
                    if GQsPass and NRsPass:
                        allGTs = [el[0] for el in allSamples]
                        tumGTs = [el[0] for el in tumSamples] 
                        normGT = normSample[0]

                        # If the genotype for all samples exist:
                        if ("./." not in allGTs) and (normGT=="0/0") and (tumGTs.count("0/0")!=(samples-1)):
                            normNVpass = (int(normSample[5]) == 0)
                            tumNVpass = any([float(el[5]) >= 3 for el in tumSamples])
                            # No variant reads for normal and
                            # >=3 variant reads for at least one tumour sample
                            if normNVpass and tumNVpass:
                                platypus_pass.write(line)
                                
                                chrom, pos = str(record[0]), int(record[1])
                                refLen, altLen = len(record[3]), len(record[4])
                                mutPlat.append(pos)
                                if (refLen == 1) and (altLen == 1):
                                    snpsPlat.append(pos)
                                else:
                                    indelsPlat.append(pos)

                                # Check if the variant is in bad regions  
                                for loc in badRegs: 
                                    if (chrom==loc[0] or chrom==loc[0][0:2]) and (pos>=loc[1]) and (pos<=loc[0]):
                                        badRegions.write(str(pos)+'\n')
                                        break 

                line = platcalls.readline()
     
                                        

#####################################################################
#Create the venn diagram
#####################################################################
if venn: 
    from pylab import *
    import matplotlib.pyplot as plt
    import matplotlib_venn as vn
    import numpy as np 

    #####################################################################
    # Extract the sites listed as SNPs or indels in the
    # VCFs produced by Mutect2.
    #####################################################################
    pos_mutect_snps = []
    pos_mutect_indels = []
    with open(str(mutect_file), 'r') as mutectCalls:
        line = ''
        #skip the header
        while line[0:6]!="#CHROM":
            line = mutectCalls.readline()
        line= mutectCalls.readline()
        while line:
            call = line.split('\t')
            refLen, altLen = len(call[3]), len(call[4])
            if refLen == 1 and altLen == 1:
                pos_mutect_snps.append(int(call[1]))
            else:
                pos_mutect_indels.append(int(call[1]))
            line = mutectCalls.readline()


        mut_plat_snps = list(set(pos_mutect_snps).intersection(snpsPlat))
        mut_plat_indels = list(set(pos_mutect_indels).intersection(indelsPlat))
        plt.figure()
        plt.subplot(211)
        v = vn.venn2(subsets=(len(pos_mutect_snps),len(snpsPlat), len(mut_plat_snps)), set_labels = ('Mutect-Snps', 'Platypus-Snps'))
        c = vn.enn2_circles(subsets=(len(pos_mutect_snps),len(snpsPlat), len(mut_plat_snps)), linestyle='dotted')
        plt.subplot(212)
        v = vn.venn2(subsets=(len(pos_mutect_indels),len(indelsPlat), len(mut_plat_indels)), set_labels = ('Mutect-Indels', 'Platypus-Indels'))
        c = vn.venn2_circles(subsets=(len(pos_mutect_indels),len(indelsPlat), len(mut_plat_indels)), linestyle='dotted')
        plt.savefig(os.path.join('', "{}.png".format(platypus_file[0:-4])))

#####################################################################