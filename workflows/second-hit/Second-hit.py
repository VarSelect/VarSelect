#!/usr/bin/python

"""
    gender: True -> Male
            False -> Female

    is_F: Father is affected
    is_M: Mother is affected

    Region: chr 1 ~ 22
            chr X:
                nonPAR -> only when child is daughter

    PAR region: p1 ~ p2, p3 ~ p4
        p1: 60001
        p2: 2699520
        p3: 154931044
        p4: 155270560
"""

import os, sys, getopt

def main(argv):
    usage = "usage: python Second-hit.py -v [vcf file] -p [ped file]"

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
        print ("try 'python Two-hits.py -h' for help")
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
            print ('#chrom\tpos\tref\talt\tis_TH\t'+ c_ + '\t' + f_ + '\t' + m_)
            # only take first affected child
            break
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
    def is_HOM(genotype):
        if len(genotype) != 3:
            return False

        g1 = genotype[0]
        g2 = genotype[2]

        if g1 == g2:
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


        nonPAR = 'X' in vlist[i][0] and (int(vlist[i][1]) < p1 or (p2 <= int(vlist[i][1]) <= p3) or int(vlist[i][1] > p4))
        if nonPAR == 1:
            nonPAR = False

        tmp = '\t'.join(vlist[i])
        # consider '.' as '0'
        tmp = tmp.replace('.', '0').split('\t')

        child = tmp[c_index].split(':')[0]
        father = tmp[f_index].split(':')[0]
        mother = tmp[m_index].split(':')[0]

        # in second-hit case, only child is daughter needs to include nonPAR region
        # only father may have one genotype
        if not(nonPAR):
            # if genotype in VCF is only one '.' or '0', consider as '0/0', excluding nonPAR region
            if father == '0':
                father = '0/0'
                tmp[f_index] = '0/0'

            f1 = tmp[f_index][0]
            f2 = tmp[f_index][2]

        if child == '0':
            child = '0/0'
            tmp[c_index] = '0/0'
        if mother == '0':
            mother = '0/0'
            tmp[m_index] = '0/0'

        c1 = tmp[c_index][0]
        c2 = tmp[c_index][2]
        m1 = tmp[m_index][0]
        m2 = tmp[m_index][2]


        # 1. both parents are unaffected:
        if not(is_F) and not(is_M):
            if not(nonPAR):
                if is_HOM(child) and ((c1 == f1 and c1 != f2 and c1 not in mother) or (c1 != f1 and c1 == f2 and c1 not in mother) or (c1 == m1 and c1 != m2 and c1 not in father) or (c1 != m1 and c1 == m2 and c1 not in father)):
                    print (true)
                else:
                    print (false)
            else:
                # in two-hits case, only include nonPAR region when child is daughter
                if gender:
                    print (false)
                    continue
                else:
                    # child is daughter, only one comes from mother, mother can not be HOM.
                    if is_HOM(child) and c1 != father and c1 in mother and not(is_HOM(mother)):
                        print (true)
                    else:
                        print (false)

        # 2. Father is affected, mother is unaffected and child is son
        elif is_F and not(is_M) and gender:
            if not(nonPAR):
                if is_HOM(child) and is_HOM(father) and child == father and c1 not in mother:
                    print (true)
                else:
                    print (false)
            else:
                # pass nonPAR region if child is son
                print (false)
                continue

        # 3. Father is unaffected, mother is affected and child is son
        elif not(is_F) and is_M and gender:
            if not(nonPAR):
                if is_HOM(child) and is_HOM(mother) and c1 not in father and child == mother:
                    print (true)
                else:
                    print (false)
            else:
                # pass nonPAR region if child is son
                print (false)
                continue

        # 4. Father is affected, mother is unaffected and child is daughter
        elif is_F and not(is_M) and not(gender):
            if not(nonPAR):
                if is_HOM(child) and is_HOM(father) and child == fahter and child not in mother:
                    print (true)
                else:
                    print (false)
            else:
                if is_HOM(child) and c1 == father and c1 not in mother:
                    print (true)
                else:
                    print (false)

        # 5. Father is unaffected, mother is affected and child is daughter
        elif not(is_F) and is_M and not(gender):
            if not(nonPAR):
                if is_HOM(child) and is_HOM(mother) and child not in father and child == mother:
                    print (true)
                else:
                    print (false)
            else:
                if is_HOM(child) and is_HOM(mother) and c1 not in father and child == mother:
                    print (true)
                else:
                    print (false)

   
if __name__ == "__main__":
    main(sys.argv[1:])
