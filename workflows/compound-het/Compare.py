#!/usr/bin/python

"""
    f = compare file
    f2 = ped file
    f3 = vcf file
    f4 = output file

    data: is_CH=1
    other: is_CH=0
"""

import os, sys, getopt, subprocess

def main(argv):
    usage = "usage: python Compare.py -v [vcf file] -p [ped file] -c [compare file]"

    try:
        opts, args = getopt.getopt(argv,"hv:p:c:",["vcf=","ped=","cfile="])
    except getopt.GetoptError:
        print (usage)
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print (usage)
            sys.exit()
        elif opt in ("-v", "--vcf"):
            vcffile = arg
        elif opt in ("-p", "--ped"):
            pedfile = arg
        elif opt in ("-c", "--cfile"):
            comparefile = arg

    if len(opts) != 3:
        print ("try 'python Compare.py -h' for help")
        sys.exit(2)


    # get first sample whom has been filterd in stage 1
    f = open(comparefile, 'r')
    FL = f.readline().strip().split('\t') # First Line [list]
    FS = FL[6] # First Sample
    FATHER = FL[7] # Father
    MOTHER = FL[8] # Mother
    f.close()

    # get all other samples, no matter it's affected or unaffected, including father and mother
    f2 = open(pedfile, 'r')

    OS = {} # Other Samples [dict]
    for line in f2:
        if line.startswith('#'):
            continue

        tmp = line.strip().split('\t')
        # it could be aunt or uncle
        if tmp[1] != FS and tmp[1] != FATHER and tmp[1] != MOTHER:
            OS[tmp[1]] = tmp[5]
        elif tmp[1] == FATHER:
            f_ = tmp[1]
        elif tmp[1] == MOTHER:
            m_ = tmp[1]

    if not OS:
#        print ("[Notice] Only one affected sample was founded.")
        f = open(comparefile, 'r')
        for line in f:
            print (line.strip())
        f.close()
        sys.exit()
    f2.close()


    # put here to avoid wasting time if not other affected sample found
    f = open(comparefile, 'r')
    tmp = []
    for line in f:
        if line.startswith('#'):
            pass
        else:
            tmp.append(line)
    f.close()

    clist = [x.strip().split('\t') for x in tmp]


    # seperate other samples into two categories
    affected = []
    unaffected = []
    for k, v in OS.items():
        if v == '2':
            affected.append(k)
        else:
            unaffected.append(k)

    # parse is_model=1 data, else store to other
    data = {}
    other = {}
    for i in range(len(clist)):
        if clist[i][4] == '1':
            data[clist[i][0]+'-'+clist[i][1]] = clist[i][5]
        else:
            other[clist[i][0]+'-'+clist[i][1]] = clist[i][5]

    # get new sample indexes, including father and mother
    f3 = open(vcffile, 'r')

    FL = '\t'.join(FL)
    affected_indexes = []
    unaffected_indexes = []
    for line in f3:
        if '#CHROM' in line:
            tmp = line.strip().split('\t')

    if affected:
        for i in affected:
            FL = FL + '\t' + i
            affected_indexes.append(tmp.index(i))

    if unaffected:
        for j in unaffected:
            FL = FL + '\t' + j
            unaffected_indexes.append(tmp.index(j))

    f_index = tmp.index(f_)
    m_index = tmp.index(m_)
    FS_index = tmp.index(FS)

    f3.close()

    # get gene from clist
    gene_all = {}
    for i in range(len(clist)):
        gene_all[clist[i][0]+'-'+clist[i][1]] = clist[i][5]

    # comparison
    def is_HET(genotype):
        if len(genotype) != 3:
            return False

        g1 = genotype[0]
        g2 = genotype[2]

        if g1 != g2:
            return True
        else:
            return False

    def is_HOM(genotype):
        if len(genotype) != 3:
            return False

        g1 = genotype[0]
        g2 = genotype[2]

        if g1 == g2:
            return True
        else:
            return False

    def CompareGenotype(new, first):
        if len(new) != 3 and len(first) != 3:
            if new == first:
                return True
            else:
                return False

        elif len(new) != 3 and len(first) == 3:
            if new in first:
                return True
            else:
                return False

        elif len(new) == 3 and len(first) != 3:
            if first[0] == new and first[2] == new:
                return True
            else:
                return False

        else:
            if new == first:
                return True
            else:
                return False

    def GetAllele(g):
        if len(g) != 3:
            return g
        else:
            return g[2]


    for a in affected_indexes:
        gene_dic = {} # put inside the for loop to make sure when every new sample added in, recheck that genotype should come from both father and mother in same gene

        f3 = open(vcffile, 'r')
        for line in f3:
            if line.startswith('#'):
                continue

            tmp = line.strip().split('\t')
            chrom = tmp[0]
            pos = tmp[1]
            father = tmp[f_index].split(':')[0].replace('.', '0')
            mother = tmp[m_index].split(':')[0].replace('.', '0')
            genotype = tmp[a].split(':')[0].replace('.', '0')
            FS_genotype = tmp[FS_index].split(':')[0].replace('.', '0')

            gene = gene_all[chrom+'-'+pos]

            if chrom+'-'+pos in data:
                # change data value to new sample's genotype
                data[chrom+'-'+pos] = genotype

                if not(CompareGenotype(genotype, FS_genotype)):
                    # change is_model=1 to is_model=0
                    del data[chrom+'-'+pos]
                    other[chrom+'-'+pos] = genotype
                else:
                    if gene in gene_dic:
                        if (gene_dic[gene] == 'F' and is_HET(mother)) or (gene_dic[gene] == 'M' and is_HET(father)):
                            gene_dic[gene] = 1
                    else:
                        if is_HET(father):
                            gene_dic[gene] = 'F'
                        elif is_HET(mother):
                            gene_dic[gene] ='M'

            # store is_model=0 data's genotype
            if chrom+'-'+pos in other:
                other[chrom+'-'+pos] = genotype

        f3.close()
        # put here to make sure every new samples's genotype will be added and remove variants which genotype does not come from both father and mother in same gene
        for i in range(len(clist)):
            if clist[i][0] + '-' + clist[i][1] in data:
                clist[i].append(data[clist[i][0]+'-'+clist[i][1]])
            elif clist[i][0] + '-' + clist[i][1] in other:
                clist[i][4] = '0'
                clist[i].append(other[clist[i][0]+'-'+clist[i][1]])

            if clist[i][5] in gene_dic:
                if gene_dic[clist[i][5]] != 1:
                    clist[i][4] = '0'

    for u in unaffected_indexes:
        gene_dic = {} # put inside the for loop to make sure when every new sample added in, recheck that genotype should come from both father and mother in same gene

        f3 = open(vcffile, 'r')
        for line in f3:
            if line.startswith('#'):
                continue

            tmp = line.strip().split('\t')
            chrom = tmp[0]
            pos = tmp[1]
            father = tmp[f_index].split(':')[0].replace('.', '0')
            mother = tmp[m_index].split(':')[0].replace('.', '0')
            genotype = tmp[u].split(':')[0].replace('.', '0')
            FS_genotype = tmp[FS_index].split(':')[0].replace('.', '0')
            allele = GetAllele(FS_genotype)

            gene = gene_all[chrom+'-'+pos]

            if chrom+'-'+pos in data:
                # change data value to new sample's genotype
                data[chrom+'-'+pos] = genotype

                if allele in genotype and is_HOM(genotype):
                    del data[chrom+'-'+pos]
                    other[chrom+'-'+pos] = genotype
                else:
                    if gene in gene_dic:
                        if (gene_dic[gene] == 'F' and is_HET(mother)) or (gene_dic[gene] == 'M' and is_HET(father)):
                            gene_dic[gene] = 1
                    else:
                        if is_HET(father):
                            gene_dic[gene] = 'F'
                        elif is_HET(mother):
                            gene_dic[gene] = 'M'

            # store is_model=0 data's genotype
            if chrom+'-'+pos in other:
                other[chrom+'-'+pos] = genotype

        f3.close()
        # put here to make sure every new samples's genotype will be added and remove variants which genotype does not come from both father and mother in same gene
        for i in range(len(clist)):
            if clist[i][0] + '-' + clist[i][1] in data:
                clist[i].append(data[clist[i][0]+'-'+clist[i][1]])
            elif clist[i][0] +'-' + clist[i][1] in other:
                clist[i][4] = '0'
                clist[i].append(other[clist[i][0]+'-'+clist[i][1]])

            if clist[i][5] in gene_dic:
                if gene_dic[clist[i][5]] != 1:
                    clist[i][4] = '0'

        # if new sample is unaffected, make sure it can not comes from both father and mother
        removed_gene = {}
        for i in range(len(clist)):
            if clist[i][4] == '1' and is_HET(clist[i][-1]):
                if clist[i][5] in removed_gene:
                    if removed_gene[clist[i][5]] == 'F' and is_HET(clist[i][8]) or removed_gene[clist[i][5]] == 'M' and is_HET(clist[i][7]):
                        removed_gene[clist[i][5]] = 1
                else:
                    if is_HET(clist[i][7]):
                        removed_gene[clist[i][5]] = 'F'
                    elif is_HET(clist[i][8]):
                        removed_gene[clist[i][5]] = 'M'

        # delete variants that unaffected sample is HET and comes from both father and mother
        for i in range(len(clist)):
            if clist[i][5] in removed_gene:
                if removed_gene[clist[i][5]] == 1 and is_HET(clist[i][-1]):
                    clist[i][4] = '0'

    # recheck affected sample of remained variants got both father and mother
    recheck_gene = {}
    for i in range(len(clist)):
        if clist[i][4] == '1':
            if clist[i][5] in recheck_gene:
                if recheck_gene[clist[i][5]] == 'F' and is_HET(clist[i][8]) or recheck_gene[clist[i][5]] == 'M' and is_HET(clist[i][7]):
                    recheck_gene[clist[i][5]] = 1
            else:
                if is_HET(clist[i][7]):
                    recheck_gene[clist[i][5]] = 'F'
                elif is_HET(clist[i][8]):
                    recheck_gene[clist[i][5]] = 'M'


    # output compared data
    print (FL)
    for i in range(len(clist)):
        if clist[i][5] in recheck_gene:
            if recheck_gene[clist[i][5]] != 1:
                clist[i][4] = '0'
        print ('\t'.join(clist[i]))

if __name__ == "__main__":
    main(sys.argv[1:])
