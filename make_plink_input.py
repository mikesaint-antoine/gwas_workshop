import csv
import numpy as np
import pandas as pd






snp_info = []
data = []
with open("data/example_data.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')

    for row in reader:

        snp_info.append(row[0:3])
        data.append(row[3:])





snp_header = snp_info[0]
snp_info = snp_info[1:]

samples = data[0]
data = data[1:]




snp_info = np.array(snp_info)
data = np.array(data).astype(int)

snp_header = np.array(snp_header)
samples = np.array(samples)





print(snp_info.shape)
print(data.shape)
print()
print(snp_header.shape)
print(samples.shape)






print("snp_header:")
print(snp_header)
print()

print("samples:")
print(samples[0:3])
print()

print("snp_info:")
print(snp_info[0:5])
print()

print("data:")
print(data[0:5,0:3])
print()




labels = []

for i in range(len(samples)):
    tmp = samples[i].split("_")
    labels.append(tmp[0])
    
labels = np.array(labels)




print(labels[0:5])

print(np.unique(labels))


# After you’ve loaded snp_info (shape M×3) and data (M×N)…

# Parse CHR and POS from the “chr:pos” string in column 0
chroms = [s.split(':')[0].replace('chr','') for s in snp_info[:, 0]]
positions= [int(s.split(':')[1]) for s in snp_info[:, 0]]

# Get REF/ALT from columns 1 & 2
ref_alleles = snp_info[:, 1]
alt_alleles = snp_info[:, 2]

# Your genotype matrix
snp_ids   = snp_info[:, 0]
genotypes = data    # M×N array of 0/1/2

# Write VCF
with open('plink_input/example_data.vcf','w') as vcf:
    vcf.write("##fileformat=VCFv4.2\n")
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
              + "\t".join(samples) + "\n")
    for i in range(len(snp_ids)):
        gt_strs = []
        for g in genotypes[i]:
            if   g == 0: gt = "0/0"
            elif g == 1: gt = "0/1"
            else:        gt = "1/1"
            gt_strs.append(gt)

        vcf.write(
            f"{chroms[i]}\t"
            f"{positions[i]}\t"
            f"{snp_ids[i]}\t"
            f"{ref_alleles[i]}\t"
            f"{alt_alleles[i]}\t"
            f".\t.\t.\tGT\t"
            + "\t".join(gt_strs)
            + "\n"
        )





# Map your string labels to PLINK codes: 1=control (normal), 2=case (zombie)
phen_map = {'normal': 1, 'zombie': 2}
phenotypes = [ phen_map[lab] for lab in labels ]

# Build the .fam DataFrame
fam_df = pd.DataFrame({
    'FID':       samples,
    'IID':       samples,
    'Father':    ['0'] * len(samples),
    'Mother':    ['0'] * len(samples),
    'Sex':       ['0'] * len(samples),       # 0=unknown
    'Phenotype': phenotypes
})

# Write it out (no header, space‑delimited)
fam_df.to_csv(
    'plink_input/example_data.fam',
    sep=' ',
    header=False,
    index=False
)








########################################
# steps to run plink2


# make binary input files, with MAF filter of 0.05
plink2   --vcf plink_input/example_data.vcf   --maf 0.05 --fam plink_input/example_data.fam   --make-bed   --out plink_input/example_data_maf05

# run 