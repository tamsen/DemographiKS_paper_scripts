import os
import unittest
import random
from matplotlib import pyplot as plt
from scipy.stats import expon


class TestModelEffectsOfCLT(unittest.TestCase):

    def test_effects_of_CLT(self):

        output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/CLT_testing"
        png_file="CLT.png"
        png_out=os.path.join(output_folder,png_file)

        #make a bunch of exponential distributions
        num_genes=1000
        N=1000
        #pairs=[[1,5],[1,10],[1,20]]
        list_of_expon_distributions=[]
        for i in range(0,5):
            data=theoretical_coalescent(num_genes, N)
            list_of_expon_distributions.append(data)

        plot_mrca(list_of_expon_distributions[0],png_out)
        self.assertEqual(True, True)  # add assertion here


def plot_mrca(theoretical_mrcas_by_gene, png_out):

    fig = plt.figure(figsize=(10, 10), dpi=350)
    # Co.T=(1/2N)*e^-((t-1)/2N))
    label = "Coalescent Times For Genes From Sampled Two Ancestral Genomes"
    num_genes=len(theoretical_mrcas_by_gene)
    plt.hist(theoretical_mrcas_by_gene, bins=100, facecolor='c', alpha=0.25,
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

def theoretical_coalescent(num_genes,N):
    #Co.T=(1/2N)*e^-((t-1)/2N))
    loc=0
    xscale=2.0*N
    random.seed(10)
    random_draws_from_distribution = expon.rvs(loc=loc, size=num_genes, scale=xscale)
    return random_draws_from_distribution

if __name__ == '__main__':
    unittest.main()
