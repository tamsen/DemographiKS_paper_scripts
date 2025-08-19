import math
import os
import unittest

import numpy as np
from matplotlib.patches import Rectangle
from matplotlib import pyplot as plt

from figure_generation.histogram_plotter import read_Ks_csv
from figure_generation.ks_plot_aggregations_dgks_vs_truth import make_Tc_Ks_vs_truth_fig_with_subplots
from modeling import ks_parsers



class Final_DGKS_vs_Empirical(unittest.TestCase):

    def test_poplar_final_plots(self):
        # linux
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Poplar'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Poplar/Truth'

        species_for_plot_title = 'EMP_Pop_07_and_half_11'
        out_png = os.path.join(demographiKS_out_path, species_for_plot_title + '_final_poplar.png')
        real_full_path = os.path.join(truth_out_path, 'poplar.ks.tsv')
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        sim_full_path = os.path.join(demographiKS_out_path,
                                     species_for_plot_title,
                                     'allotetraploid_bottleneck.csv')
        demographiKS_ks_results = read_Ks_csv(sim_full_path, False)
        bin_size = 0.01
        max_Ks = 0.4
        bins = np.arange(bin_size, max_Ks + 0.1, bin_size)
        hist_ys_real, bins_real, patches = plt.hist(real_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                                    density=True)

        hist_ys_sim, bins_sim, patches = plt.hist(demographiKS_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                                  density=True)

        hist_ys_1 = [hist_ys_real, bins_real]
        hist_ys_2 = [hist_ys_sim, bins_sim]
        list_of_hist_data = [hist_ys_1, hist_ys_2]

        overlay_differences_in_curves(species_for_plot_title, list_of_hist_data, out_png)

        self.assertEqual(True, True)  # add assertion here

    def test_maize_final_plots(self):
        # linux
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize/Truth'


        species_for_plot_title = 'EMP_Mays_26_29_combined'
        out_png = os.path.join(demographiKS_out_path, species_for_plot_title+ '_final_mays.png')
        real_full_path = os.path.join(truth_out_path, 'mays.ks.tsv')
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        sim_full_path = os.path.join(demographiKS_out_path,
                                     species_for_plot_title,
                                     'allotetraploid_bottleneck.csv')
        demographiKS_ks_results = read_Ks_csv(sim_full_path, False)
        bin_size = 0.01
        max_Ks = 0.4
        bins = np.arange(bin_size, max_Ks + 0.1, bin_size)
        hist_ys_real, bins_real, patches = plt.hist(real_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                                    density=True)

        hist_ys_sim, bins_sim, patches = plt.hist(demographiKS_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                                  density=True)

        hist_ys_1 = [hist_ys_real, bins_real]
        hist_ys_2 = [hist_ys_sim, bins_sim]
        list_of_hist_data = [hist_ys_1, hist_ys_2]

        overlay_differences_in_curves(species_for_plot_title, list_of_hist_data, out_png)

        self.assertEqual(True, True)  # add assertion here


    def test_coffee_final_plots(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee/Truth'



        #species_for_plot_title = 'EMP_Coff_36_m07d09y2025_h12m22s40'
        species_for_plot_title = 'EMP_Coff_35_m07d01y2025_h08m58s51'
        out_png = os.path.join(demographiKS_out_path,species_for_plot_title + '_final_coffee.png')
        real_full_path = os.path.join(truth_out_path, 'coffea.ks.tsv')
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        sim_full_path = os.path.join(demographiKS_out_path,
                                     species_for_plot_title,
                                     'allotetraploid_bottleneck.csv')
        demographiKS_ks_results = read_Ks_csv(sim_full_path, False)
        bin_size = 0.005
        max_Ks = 0.2
        bins = np.arange(bin_size, max_Ks + 0.1, bin_size)
        hist_ys_real, bins_real, patches = plt.hist(real_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                                    density=True)

        hist_ys_sim, bins_sim, patches = plt.hist(demographiKS_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                                  density=True)

        hist_ys_1 = [hist_ys_real, bins_real]
        hist_ys_2 = [hist_ys_sim, bins_sim]
        list_of_hist_data = [hist_ys_1, hist_ys_2]

        overlay_differences_in_curves(species_for_plot_title, list_of_hist_data, out_png)

        self.assertEqual(True, True)  # add assertion here


def overlay_differences_in_curves(species_for_plot_title, list_of_hist_data, out_png):

    colors = ['green','blue','brown']
    labels = ['truth','model fit']


    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    for i in range(0,len(list_of_hist_data)):
        hist_data =list_of_hist_data[i]
        [n, bins]=hist_data
        width=(bins[2]-bins[1])
        #bar_plot_xs=[b + i*width for b in bins[0:len(bins)-1]] #to match bar-plot axes
        bar_plot_xs = [b + width for b in bins[0:len(bins) - 1]]  # to match bar-plot axes

        label=labels[i]

        plt.plot(bar_plot_xs,n,c=colors[i], alpha=1, label=label)


    diffs = [(list_of_hist_data[0][0][j]-list_of_hist_data[1][0][j]) for j in range(0,len(list_of_hist_data[1][0]))]
    rmse=math.sqrt(sum([d*d for d in diffs])/len(diffs))
    hist_data0 = list_of_hist_data[1]
    ys0=hist_data0[0]
    xs0=hist_data0[1]

    for i in range(0,len(ys0)):
        if i==0:
            ax.add_patch(Rectangle((xs0[i], ys0[i]), width, diffs[i], color=colors[2],
                                alpha=0.15, label="rmse: "+ str(round(rmse,4))))
        else:
            ax.add_patch(Rectangle((xs0[i], ys0[i]), width, diffs[i], color=colors[2],
                                alpha=0.15))

    plt.legend()
    plt.xlabel("Ks")
    ax.set(ylabel="Density")
    plt.title(species_for_plot_title, style='italic')
    plt.tight_layout()
    plt.savefig(out_png)
    plt.clf()
    plt.close()

    return n, bins



if __name__ == '__main__':
    unittest.main()
