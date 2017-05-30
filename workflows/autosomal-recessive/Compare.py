#!/usr/bin/python

"""
    f = compare file 
    f2 = ped file
    f3 = vcf file
    f4 = output file 

    data: is_AR=1
    other: is_AR=0
"""

import os, sys, getopt

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
    FS = FL[5] # First Sample
    FATHER = FL[6] # Father
    MOTHER = FL[7] # Mother
    f.close()

    # get all other samples, no matter it's affected or unaffected
    f2 = open(pedfile, 'r')

    OS = {} # Other Samples [dict]
    for line in f2:
        if line.startswith('#'):
            continue

        tmp = line.strip().split('\t')
        # it could be aunt or uncle
        if tmp[1] != FS and tmp[1] != FATHER and tmp[1] != MOTHER:
            OS[tmp[1]] = tmp[5]

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

    # get new samples indexes
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

    FS_index = tmp.index(FS)
    f3.close()


    # comparison
    for a in affected_indexes:
        f3 = open(vcffile, 'r')
        for line in f3:
            if line.startswith('#'):
                continue

            tmp = line.strip().split('\t')
            chrom = tmp[0]
            pos = tmp[1]
            genotype = tmp[a].split(':')[0].replace('.', '0')
            FS_genotype = tmp[FS_index].split(':')[0].replace('.', '0') # First sample's genotype

            if chrom+'-'+pos in data:
                # change data value to new sample's genotype
                data[chrom+'-'+pos] = genotype

                # must be the same as First Sample
                if genotype != FS_genotype:
                    # change is_model=1 to is_model=0
                    del data[chrom+'-'+pos]
                    other[chrom+'-'+pos] = genotype

            # store is_model=0 data's genotype
            if chrom+'-'+pos in other:
                other[chrom+'-'+pos] = genotype

        f3.close()
        # put here to make sure every new samples's genotype will be added
        for i in range(len(clist)):
            if clist[i][0] + '-' + clist[i][1] in data:
                clist[i].append(data[clist[i][0]+'-'+clist[i][1]])
            elif clist[i][0] +'-' + clist[i][1] in other:
                clist[i][4] = '0'
                clist[i].append(other[clist[i][0]+'-'+clist[i][1]])

    for u in unaffected_indexes:
        f3 = open(vcffile, 'r')
        for line in f3:
            if line.startswith('#'):
                continue

            tmp = line.strip().split('\t')
            chrom = tmp[0]
            pos = tmp[1]
            genotype = tmp[u].split(':')[0].replace('.', '0')
            FS_genotype = tmp[FS_index].split(':')[0].replace('.', '0') # First sample's genotype

            if chrom+'-'+pos in data:
                # change data value to new sample's genotype
                data[chrom+'-'+pos] = genotype

                # can't be the same as First Sample
                if genotype == FS_genotype:
                    # change is_model=1 to is_model=0
                    del data[chrom+'-'+pos]
                    other[chrom+'-'+pos] = genotype

            # store is_model=0 data's genotype
            if chrom+'-'+pos in other:
                other[chrom+'-'+pos] = genotype

        f3.close()
        # put here to make sure every new samples's genotype will be added
        for i in range(len(clist)):
            if clist[i][0] + '-' + clist[i][1] in data:
                clist[i].append(data[clist[i][0]+'-'+clist[i][1]])
            elif clist[i][0] + '-' + clist[i][1] in other:
                clist[i][4] = '0'
                clist[i].append(other[clist[i][0]+'-'+clist[i][1]])

    # output compared data
    print (FL)
    for i in range(len(clist)):
        print ('\t'.join(clist[i]))

if __name__ == "__main__":
    main(sys.argv[1:])
