import sys,os
import math
import numpy as np


def test():
    N = 25
    n_list = [5,8,6,2,4]
    K = 5
    S1 = 0
    S2 = 0
    S3 = 0
    for n in n_list:
        # method 1
        S1 += (n/N)*math.log(n/N)
        # method 2
        S2 += n*(n-1)
        # method 3
        S3 += (n/N)**2
    Shannon_index = -S1
    Simpson_index = S2/(N*(N-1))
    ds = S3
    Gini_Simpson_index_1 = 1 - ds
    inverse_Simpson_index_1 = 1 / ds
    evenness_index_1 = (1/ds)/K
    print("{}\t{}\t{}\t{}\t{}\t{}".format(Shannon_index, Simpson_index, ds,
            Gini_Simpson_index_1, inverse_Simpson_index_1, evenness_index_1))

def main():
    print("Sample_name\tShannon_index(pi*ln*pi)\tdiversity_index(Ds)\tGini-Simpson_index(D's=1-Ds)\tinverse_Simpson's_index(D-1s=1/Ds)\tevenness(E=D-1s/K)")
    for each_file in sys.argv[1:]:
        sample_name = os.path.basename(each_file).split(".")[0]  # RPG527P2-G164E1L1.mixcr.out.clonotypes.TRA.txt

        N = 0
        n_list = []
        K = 0
        with open(each_file, "r") as f:
            for line in f:
                if line.startswith("cloneId"):
                    continue
                line = line.rstrip("\n")
                each = line.split()
                count = float(each[1])
                N += count
                K += 1
                n_list.append(count)
        # method 1 

        S1 = 0
        S2 = 0
        S3 = 0
        for n in n_list:
            # method 1
            S1 += ((n/N) * np.log(n/N))
            # method 2
            S2 += n*(n-1)
            # method 3
            S3 += (n/N)**2
        Shannon_index = -S1
        Simpson_index = S2/(N*(N-1))
        ds = S3

        Gini_Simpson_index_1 = 1 - ds
        inverse_Simpson_index_1 = 1 / ds
        evenness_index_1 = (1/ds)/K

        #####################
        '''
        Gini_Simpson_index_2 = 1-Simpson_index
        inverse_Simpson_index_2 = 1 / Simpson_index
        evenness_index_2 = (1/Simpson_index)/K
        '''
        print("{}\t{}\t{}\t{}\t{}\t{}".format(sample_name, Shannon_index, ds, 
            Gini_Simpson_index_1, inverse_Simpson_index_1, evenness_index_1))

if __name__ == '__main__':
    main()
    # test()
