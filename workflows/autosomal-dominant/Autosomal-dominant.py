#!/usr/bin/python

"""
    gender: True -> Male
            False -> Female
    (not used)

    is_F: Father is affected
    is_M: Mother is affected

    Region: chr 1 ~ 22
"""

import os, sys, getopt

def main(argv):
    usage = "usage: python Autosomal-dominant.py -v [vcf file] -p [ped file]"

    try:
        opts, args = getopt.getopt(argv,"hv:p:",["vcf=","ped="])
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

    if len(opts) != 2:
        print ("try 'python Autosomal-dominant.py -h' for help")
        sys.exit(2)


    # get all parents from PED
    f = open(pedfile, 'r')
    parents = []
    for line in f:
        if line.startswith('#'):
            continue

        tmp = line.strip().split('\t')
        if tmp[2] != '0':
            if tmp[2] not in parents:
                parents.append(tmp[2])
        if tmp[3] != '0':
            if tmp[3] not in parents:
                parents.append(tmp[3])
    # Re-read file
    f.seek(0)
    # get first affected sample
    gender = False
    for line in f:
        if line.startswith('#'):
            continue

        tmp = line.strip().split('\t')
        if tmp[1] not in parents and tmp[5] == '2':
            c_ = tmp[1]
            f_ = tmp[2]
            m_ = tmp[3]
            if tmp[4] == '1':
                gender = True
            # print header
            print ('#chrom\tpos\tref\talt\tis_AD\t'+ c_ + '\t' + f_ + '\t' + m_)
            break # only take one affected sample
    # Re-read file
    f.seek(0)
    # check if father or mother affected
    is_F = False
    is_M = False
    for line in f:
        if line.startswith('#'):
            continue

        tmp = line.strip().split('\t')
        if tmp[1] == f_:
            if tmp[5] == '2':
                is_F = True
        elif tmp[1] == m_:
            if tmp[5] == '2':
                is_M = True
    f.close()

    # make vlist from VCF
    f2 = open(vcffile, 'r')

    tmp = []
    for line in f2:
        if '#CHROM' in line:
            cfm = line.strip().split('\t')
            c_index = cfm.index(c_)
            f_index = cfm.index(f_)
            m_index = cfm.index(m_)
        elif line.startswith('#'):
            pass
        else:
            tmp.append(line)
    f2.close()

    vlist = [x.strip().split('\t') for x in tmp]


    # Filtering
    def is_HET(genotype):
        if len(genotype) != 3:
            return False

        g1 = genotype[0]
        g2 = genotype[2]

        if (g1 != g2):
            return True
        else:
            return False

    def is_HOM(genotype):
        if len(genotype) != 3:
            return False

        g1 = genotype[0]
        g2 = genotype[2]

        if (g1 == g2):
            return True
        else:
            return False


    for i in range(len(vlist)):
        true = vlist[i][0] + '\t' + vlist[i][1] + '\t' + vlist[i][3] + '\t' + vlist[i][4] + '\t' + '1' + '\t' + vlist[i][c_index].split(':')[0] + '\t' + vlist[i][f_index].split(':')[0] + '\t' + vlist[i][m_index].split(':')[0]
        false = vlist[i][0] + '\t' + vlist[i][1] + '\t' + vlist[i][3] + '\t' + vlist[i][4] + '\t' + '0' + '\t' + vlist[i][c_index].split(':')[0] + '\t' + vlist[i][f_index].split(':')[0] + '\t' + vlist[i][m_index].split(':')[0]

        # pass if chrom = X, Y, MT
        if 'X' in vlist[i][0] or 'Y' in vlist[i][0] or 'MT' in vlist[i][0]:
            print (false)
            continue


        tmp = '\t'.join(vlist[i])
        # consider '.' as '0'
        tmp = tmp.replace('.', '0').split('\t')

        child = tmp[c_index].split(':')[0]
        father = tmp[f_index].split(':')[0]
        mother = tmp[m_index].split(':')[0]

        # if genotype in VCF is only one '.' or '0', consider as '0/0'
        if child == '0':
            child = '0/0'
            tmp[c_index] = '0/0'
        if father == '0':
            father = '0/0'
            tmp[f_index] = '0/0'
        if mother == '0':
            mother = '0/0'
            tmp[m_index] = '0/0'

        c1 = tmp[c_index][0]
        c2 = tmp[c_index][2]


        # 1. Father is affected, mother is unaffected
        if is_F and not(is_M) and gender:
            if is_HET(child) and is_HET(father) and is_HOM(mother) and (c1 in father and c2 in mother or c1 in mother and c2 in father):
                print (true)
            else:
                print (false)

        # 2. Father is unaffected, mother is affected
        elif not(is_F) and is_M and gender:
            if is_HET(child) and is_HOM(father) and is_HET(mother) and (c1 in father and c2 in mother or c1 in mother and c2 in father):
                print (true)
            else:
                print (false)


if __name__ == "__main__":
    main(sys.argv[1:])
