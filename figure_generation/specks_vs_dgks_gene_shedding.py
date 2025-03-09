import math
import os.path
import unittest

import numpy as np
import tskit
from matplotlib import pyplot as plt
from scipy.stats import expon

import config
from modules import trees_file_processor


class Testspecks_vs_dgks_gene_shedding(unittest.TestCase):

    def test_specks_vs_dgks_gene_shedding(self):

        output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Twgd"
        png_out = os.path.join(output_folder, "Num shed genes vs over time.png")
        fig = plt.figure(figsize=(4, 4), dpi=100)
        label = "Num shed genes vs over time"
        specKS_shed_genes = [0,3333-3326,3333-3296,3333-3260]
        demographiKS_shed_genes = [0, 3333-3326, 3333-3296, 3333-3260]
        t_wgd = [10000, 100000, 500000, 1000000]
        avg_WGD_gene_lifespan = 31e6 / math.log(2) #half life to mean life conversion

        t_test = [10000, 50000, 100000, 500000, 1000000]
        shed_fraction = [1.0 - avg_WGD_gene_lifespan*expon.pdf(
            t, loc=0, scale=avg_WGD_gene_lifespan) for t in t_test]
        expectations=[f*3333 for f in shed_fraction]

        for i in range(0, len(specKS_shed_genes)):
            plt.scatter(t_wgd[i],specKS_shed_genes[i], c='r', marker="x", label="SpecKS")
            plt.scatter(t_wgd[i],demographiKS_shed_genes[i], c='b', marker="+", label="DemographiKS")

        plt.plot(t_test, expectations, c='gray',  label="expectations")

        plt.title(label)
        plt.legend()
        # plt.ylim([0,1200])
        plt.xlabel("Time (gen)")
        plt.ylabel("Num shed genes")
        plt.savefig(png_out)
        plt.clf()
        plt.close()

        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
