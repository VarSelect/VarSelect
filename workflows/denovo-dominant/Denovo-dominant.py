#!/usr/bin/python

"""
    gender: True -> Male
            False -> Female

    is_F: Father is affected
    is_M: Mother is affected

    Region: chr 1 ~ 22
            chr X:
                only PAR region

    PAR region: p1 ~ p2, p3 ~ p4
        p1: 60001
        p2: 2699520
        p3: 154931044
        p4: 155270560
"""

import os, sys, getopt

def main(argv):
    usage = "usage: python Denovo-dominant.py -v [vcf file] -p [ped file]"

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
        print ("try 'python Denovo-dominant.py -h' for help")
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
            print ('#chrom\tpos\tref\talt\tis_DD\t'+ c_ + '\t' + f_ + '\t' + m_)
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

        if g1 != g2:
            return True
        else:
            return False


    p1 = 60001
    p2 = 2699520
    p3 = 154931044
    p4 = 155270560
    for i in range(len(vlist)):
        true = vlist[i][0] + '\t' + vlist[i][1] + '\t' + vlist[i][3] + '\t' + vlist[i][4] + '\t' + '1' + '\t' + vlist[i][c_index].split(':')[0] + '\t' + vlist[i][f_index].split(':')[0] + '\t' + vlist[i][m_index].split(':')[0]
        false = vlist[i][0] + '\t' + vlist[i][1] + '\t' + vlist[i][3] + '\t' + vlist[i][4] + '\t' + '0' + '\t' + vlist[i][c_index].split(':')[0] + '\t' + vlist[i][f_index].split(':')[0] + '\t' + vlist[i][m_index].split(':')[0]

        # pass if chrom = Y or MT
        if 'Y' in vlist[i][0] or 'MT' in vlist[i][0]:
            print (false)
            continue

        # pass if in nonPAR region
        if 'X' in vlist[i][0] and ((int(vlist[i][1]) < p1) or (p2 <= int(vlist[i][1]) <= p3) or (int(vlist[i][1]) > p4)):
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

        c1 = tmp[c_index][0]
        c2 = tmp[c_index][2]

        if mother == '0':
            mother = '0/0'
            tmp[m_index] = '0/0'


        if is_HET(child) and (((c1 in father or c1 in mother) and (c2 not in father and c2 not in mother)) or ((c1 not in father and c1 not in mother) and (c2 in father or c2 in mother))):
            print (true)
        else:
            print (false)


if __name__ == "__main__":
    main(sys.argv[1:])
