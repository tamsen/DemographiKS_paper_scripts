import os
import unittest
import random

import numpy as np
from matplotlib import pyplot as plt
from mkl import second
from scipy.stats import expon


class TestModelEffectsOfCLT(unittest.TestCase):

    def test_effects_of_CLT(self):

        output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/CLT_testing"
        png_out_0=os.path.join(output_folder,"CLT_0RC.png")

        #make a bunch of exponential distributions
        num_genes=1000
        N=1000
        xmax=10000
        bin_size = xmax/50
        bins = np.arange(0, xmax, bin_size)

        list_of_expon_distributions=[]
        for i in range(0,50):
            random.seed(i*10+10)
            data=list(theoretical_coalescent(num_genes, N))
            list_of_expon_distributions.append(data)
        plot_mrca(list_of_expon_distributions[0], bins, png_out_0)
        dist_subsampled=[]
        #build a new Tc drawing from two different distributions
        for i in range(2,10):

            segment_size=int(num_genes/i)
            data_so_far=[]
            subsampled = []
            for j in range(0,i):
                data_chunk=list(list_of_expon_distributions[j][j*segment_size:(j+1)*segment_size])
                data_so_far = data_so_far+data_chunk
                subsampled.append(j)
            png_out_i = os.path.join(output_folder, "CLT_" + str(i) +"RC.png")
            print(str(i) + ": subsampled " + str(subsampled))
            dist_subsampled = dist_subsampled + subsampled
            plot_mrca(data_so_far,bins,png_out_i)

        self.assertEqual(True, True)  # add assertion here


def plot_mrca(theoretical_mrcas_by_gene, bins, png_out):

    fig = plt.figure(figsize=(10, 10), dpi=350)
    label = "Coalescent Times For Genes From Sampled Two Ancestral Genomes"
    num_genes=len(theoretical_mrcas_by_gene)
    dgx_hist_ys, bins, patches = plt.hist(theoretical_mrcas_by_gene, bins=bins, facecolor='c', alpha=0.25,
                 label='Theoretical Tcoal by gene (total: ' + str(num_genes) + ')',
                 density=False)

    #plt.xlim([0, max_mrca * (1.1)])
    plt.title(label)
    # plt.xlim([0, max_mrca])
    # plt.ylim([0, 300])
    plt.legend()
    plt.xlabel("MRCA time")
    plt.ylabel("# genes in bin")
    plt.savefig(png_out)
    plt.clf()
    plt.close()
    return dgx_hist_ys, bins

def theoretical_coalescent(num_genes,N):
    #Co.T=(1/2N)*e^-((t-1)/2N))
    loc=0
    xscale=2.0*N
    random_draws_from_distribution = expon.rvs(loc=loc, size=num_genes, scale=xscale)
    return random_draws_from_distribution

if __name__ == '__main__':
    unittest.main()
