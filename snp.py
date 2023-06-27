#PBS -N snp
#PBS -l nodes=1:ppn=8
#PBS -q batch
#PBS -V
#PBS -S /bin/bash
import numpy as np


data_type = np.dtype({'names':['chr' , 'pos' , 'ref' , 'alt'] , 
                      'formats':['U8' , np.int , 'U4' , 'U4']})

data_type_1 = np.dtype({'names':['chr' , 'pos' , 'ref' , 'maternal' , 'paternal'] , 
                      'formats':['U8' , np.int , 'U4' , 'U4' , 'U4']})


maternal = np.loadtxt('/public/home/shidetong/wedata/snp/C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf' , skiprows=69 , usecols=(0 , 1 , 3 , 4) , dtype = data_type)
paternal = np.loadtxt('/public/home/shidetong/wedata/snp/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf' , skiprows=69 , usecols=(0 , 1 , 3 , 4) , dtype = data_type)



chrom = [str(x) for x in range(1 , 20)] + ['X' , 'Y']



data = []
for g in chrom[1]:
    tmp_m = maternal[maternal['chr'] == g]
    tmp_p = paternal[paternal['chr'] == g]
    for i in tmp_m:
        overlap = tmp_p[tmp_p['pos'] == i['pos']]
        if overlap.size != 0:
            chro = g
            pos = i['pos']
            ref = i['ref']
            m = i['alt']
            p = overlap[0]['alt']
        else:
            chro = g
            pos = i['pos']
            ref = i['ref']
            m = i['alt']
            p = i['ref']
        data.append((chro , pos , ref , m , p))
    for i in tmp_p:
        overlap = tmp_m[tmp_m['pos'] == i['pos']]
        if overlap.size == 0:
            chro = g
            pos = i['pos']
            ref = i['ref']
            m = i['ref']
            p = i['alt']
        data.append((chro , pos , ref , m , p))
    
     
    
data = np.array(data , dtype = data_type_1)

out = open('/public/home/shidetong/wedata/snp/SNP_test.txt' , 'w')
out.writelines('\t'.join(list(data_type_1.names)) + '\n')

for i in data:
    out.writelines('\t'.join([str(x) for x in i]) + '\n')
out.close()
