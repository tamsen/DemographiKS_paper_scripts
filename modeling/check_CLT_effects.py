import math
import os
import unittest
import random

import numpy as np
from matplotlib import pyplot as plt
from mkl import second
from scipy.stats import expon

from figure_generation import curve_fitting
from figure_generation.colors_for_figures import lighten_color


class TestModelEffectsOfCLT(unittest.TestCase):

# "A function that can be generalized to both a Gaussian (normal distribution) "
# "and an exponential distribution is the ")generalized normal distribution"
# (also known as the "exponential power distribution"); by adjusting a single parameter
# within this function, you can smoothly transition between the characteristics
# of a Gaussian and an exponential distribution depending on the parameter value"
# https://en.wikipedia.org/wiki/Generalized_normal_distribution
# https://en.wikipedia.org/wiki/Laplace_transform
    def test_effects_of_CLT(self):

        output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/CLT_testing"

        #make a bunch of exponential distributions
        num_genes=50000
        N=1000
        xmax=10000
        bin_size = xmax/50
        bins = np.arange(0, xmax, bin_size)

        list_of_expon_distributions=[]
        for num_combined_distributions in range(0,100):
            random.seed(num_combined_distributions*10+10)
            data=list(theoretical_coalescent(num_genes, N))
            list_of_expon_distributions.append(data)


        #build a new Tc drawing from multiple different distributions
        resampled_distributions_ys=[]
        data_labels=[]


        two_Ne = 2.0 * N
        print("Tc plot bin_size_in_time: " + str(bin_size))
        kingman = [min(num_genes,
                       (bin_size * num_genes / two_Ne) * math.e ** ((-1 * i) / two_Ne))
                   for i in bins[0:-1]]
        resampled_distributions_ys.append(kingman)
        data_labels.append("Perfect Kingman distribution")
        mu=1.0/two_Ne


        bin_size = bins[1] - bins[0]
        #wgd_normal(x, amp, mu, sig):
        popt = [num_genes * bin_size, two_Ne- (bin_size/2),bin_size]
        gaussian_ys=[curve_fitting.wgd_normal(x, *popt) for x in bins[0:-1]]
        resampled_distributions_ys.append(gaussian_ys)
        data_labels.append("Perfect Gaussian distribution")

        new_data = [d for d in list_of_expon_distributions[0]]
        png_out_i = os.path.join(output_folder, "CLT_0RC.png")
        dgx_hist_ys, bins = plot_mrca(new_data, bins, png_out_i)
        resampled_distributions_ys.append(dgx_hist_ys)
        data_labels.append("0 RC events simulated")

        data_colors = ['darkblue','darkred', 'black']

        list_of_fraction_RC_events_to_test = [0.05,0.1,0.5,0.75]
        blue_lightening=[0.75,0.6,0.4,0.3]
        for fraction in list_of_fraction_RC_events_to_test:

            new_data = [d for d in list_of_expon_distributions[0]]
            for k in range(0, int(num_genes*fraction)):
                Tc1=list_of_expon_distributions[0][k]
                Tc2=list_of_expon_distributions[1][k]
                new_data[k] = 0.5 * (Tc1 + Tc2)

            png_out_i = os.path.join(output_folder, "CLT_" + str(fraction) + "%RC.png")
            dgx_hist_ys, bins = plot_mrca(new_data, bins, png_out_i)
            resampled_distributions_ys.append(dgx_hist_ys)
            data_labels.append(str(fraction) + " %RC events simulated")

        new_colors=[lighten_color('blue', amount=a) for a in blue_lightening]
        data_colors= data_colors+ new_colors

        list_of_num_RC_events_to_test=[1,2,5,10,20,50,80]
        red_lightening=[0.2,0.5,0.6,0.65,0.7,0.75]
        for num_RC_events in list_of_num_RC_events_to_test:
            num_combined_distributions=num_RC_events+1
            new_data=[0 for k in range(0, num_genes)]
            for k in range(0, num_genes):
                 Tcs_for_Kth_genes=[list_of_expon_distributions[j][k]
                                    for j in range(0, num_combined_distributions)]
                 avg_Tc_for_Kth_gene= sum(Tcs_for_Kth_genes)/num_combined_distributions
                 new_data[k]=avg_Tc_for_Kth_gene#
            png_out_i = os.path.join(output_folder, "CLT_" + str(num_RC_events) +"RC.png")
            dgx_hist_ys, bins = plot_mrca(new_data,bins,png_out_i)
            resampled_distributions_ys.append(dgx_hist_ys)
            data_labels.append(str(num_RC_events) + " RC events simulated")

        data_colors.append("gray")
        new_colors=[lighten_color('red', amount=a) for a in red_lightening]
        data_colors= data_colors+ new_colors

        png_out = os.path.join(output_folder, "composite_CLT_RC.png")
        plot_composite_tc(resampled_distributions_ys,bins,data_labels,data_colors,png_out)
        self.assertEqual(True, True)  # add assertion here

def plot_composite_tc(data_list, bins, data_labels, data_colors, png_out):

    fig = plt.figure(figsize=(10, 10), dpi=350)
    # Co.T=(1/2N)*e^-((t-1)/2N))
    label = "Coalescent Times For Genes From Sampled From Two Ancestral Genomes"
    num_genes=len(data_list[0])

    for i in range(0,len(data_list)):
        data=data_list[i]
        plt.plot(bins[0:-1],data,color=data_colors[i], alpha=0.95,label=data_labels[i])

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
