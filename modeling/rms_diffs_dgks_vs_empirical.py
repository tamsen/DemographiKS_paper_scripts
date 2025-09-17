import math
import os
import unittest
import random
from xml.etree.ElementInclude import include

import numpy as np
from matplotlib.patches import Rectangle
from matplotlib import pyplot as plt

from figure_generation.histogram_plotter import read_Ks_csv
from figure_generation.ks_plot_aggregations_dgks_vs_truth import make_Tc_Ks_vs_truth_fig_with_subplots
from modeling import ks_parsers



class Final_DGKS_vs_Empirical(unittest.TestCase):


    def test_sugar_cane_final_plots(self):
        # linux
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Sugarcane'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Sugarcane/Truth'

        species_for_plot_title = 'EMP_Sac_10_m08d26y2025_h17m19s09'
        out_png = os.path.join(demographiKS_out_path, species_for_plot_title + '_final_saccharum.png')
        real_full_path = os.path.join(truth_out_path, 'saccharum.ks.tsv')
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

    def test_poplar_final_plots(self):
        # linux
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Poplar'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Poplar/Truth'
        include_selection = True
        species_for_plot_title = 'EMP_Pop_07_and_half_11'
        out_png = os.path.join(demographiKS_out_path, species_for_plot_title + '_final_poplar.png')
        real_full_path = os.path.join(truth_out_path, 'poplar.ks.tsv')
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        sim_full_path = os.path.join(demographiKS_out_path,
                                     species_for_plot_title,
                                     'allotetraploid_bottleneck.csv')
        demographiKS_ks_results = read_Ks_csv(sim_full_path, False)
        if include_selection:
            demographiKS_ks_results = add_selection(demographiKS_ks_results)

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

        overlay_differences_in_curves(species_for_plot_title, list_of_hist_data,
                                      include_selection, out_png)

        self.assertEqual(True, True)  # add assertion here

    def test_maize_final_plots(self):
        # linux
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize/Truth'
        include_selection = True

        species_for_plot_title = 'EMP_Mays_26_29_combined'
        out_png = os.path.join(demographiKS_out_path, species_for_plot_title+ '_final_mays.png')
        real_full_path = os.path.join(truth_out_path, 'mays.ks.tsv')
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        sim_full_path = os.path.join(demographiKS_out_path,
                                     species_for_plot_title,
                                     'allotetraploid_bottleneck.csv')
        demographiKS_ks_results = read_Ks_csv(sim_full_path, False)
        if include_selection:
            demographiKS_ks_results = add_selection(demographiKS_ks_results)

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

        overlay_differences_in_curves(species_for_plot_title, list_of_hist_data,
                                      include_selection, out_png)

        self.assertEqual(True, True)  # add assertion here


    def test_coffee_final_plots(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee/Truth'
        include_selection=True


        #species_for_plot_title = 'EMP_Coff_36_m07d09y2025_h12m22s40'
        species_for_plot_title = 'EMP_Coff_35_m07d01y2025_h08m58s51'

        out_png = os.path.join(demographiKS_out_path,species_for_plot_title + '_final_coffee.png')
        real_full_path = os.path.join(truth_out_path, 'coffea.ks.tsv')
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        sim_full_path = os.path.join(demographiKS_out_path,
                                     species_for_plot_title,
                                     'allotetraploid_bottleneck.csv')
        demographiKS_ks_results = read_Ks_csv(sim_full_path, False)

        if include_selection:
            demographiKS_ks_results = add_selection(demographiKS_ks_results)

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


        overlay_differences_in_curves(species_for_plot_title, list_of_hist_data, include_selection, out_png)

        self.assertEqual(True, True)  # add assertion here

    def test_sugar_final_plots(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Sugarcane'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Sugarcane/Truth'
        include_selection=True
        proportion_retained_genes=10
        seed = 20
        max_Ks = 0.5
        bin_size = 0.010
        species_for_plot_title =  ['EMP_Sac_36_m09d12y2025_h10m11s47',
                                   'EMP_Sac_37_m09d12y2025_h10m11s50']
        #    'EMP_Sac_37_m09d12y2025_h10m11s50'

        species_for_plot_title_string="_".join(species_for_plot_title)
        out_png = os.path.join(demographiKS_out_path,species_for_plot_title_string + '_final_sugar.png')
        real_full_path = os.path.join(truth_out_path, 'saccharum.ks.tsv')
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)



        demographiKS_ks_results=[]
        for sim_run in species_for_plot_title:
            sim_full_path = os.path.join(demographiKS_out_path,
                                         sim_run,
                                         'allotetraploid_bottleneck.csv')
            sim_run_ks_results = read_Ks_csv(sim_full_path, False)
            demographiKS_ks_results = demographiKS_ks_results + sim_run_ks_results

        if include_selection:
            demographiKS_ks_results = add_selection(demographiKS_ks_results,
                                                    proportion_retained_genes,
                                                    4.0,
                                                    seed)

        bins = np.arange(bin_size, max_Ks, bin_size)
        hist_ys_real, bins_real, patches = plt.hist(real_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                                    density=True)

        hist_ys_sim, bins_sim, patches = plt.hist(demographiKS_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                                  density=True)

        hist_ys_1 = [hist_ys_real, bins_real]
        hist_ys_2 = [hist_ys_sim, bins_sim]
        list_of_hist_data = [hist_ys_1, hist_ys_2]


        overlay_differences_in_curves(species_for_plot_title, list_of_hist_data, include_selection, out_png)

        self.assertEqual(True, True)  # add assertion here


def overlay_differences_in_curves(species_for_plot_title, list_of_hist_data,
                                  include_selection,
                                  out_png):

    colors = ['green','blue','brown']
    labels = ['truth','model fit']
    if include_selection:
        out_png = out_png.replace(".png","_with_selection.png")

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

def get_maintained_gene_Ks_values(num_Ks_values_needed, max_Ks):
    start_range = 0.0
    end_range = max_Ks
    random_decimal_list = [random.uniform(start_range, end_range) for x in range(0, num_Ks_values_needed)]
    return random_decimal_list

# 
# derive exponential decay from birth death process equation with escape from decay
# This is a standard first-order linear differential equation. With an initial population size of \(N_{0}=M(0)\),
# the solution is:\(M(t)=N_{0}e^{-kt}=N_{0}e^{-(\mu +\nu -\lambda )t}\)
# Birth/birth-death processes and their computable transition probabilities with biological applications
def add_selection(demographiKS_ks_results,fraction_genes_maintained,max_Ks, seed):
    # add pairs maintained by selection
    # pick some % of the genome to be maintained,
    # from a random distribution between Ks= 0 to 4

    random.seed(seed)
    num_Ks_values_needed = int(fraction_genes_maintained*len(demographiKS_ks_results))
    maintained_gene_Ks_values = get_maintained_gene_Ks_values(num_Ks_values_needed,max_Ks)
    demographiKS_plus_selection_results = demographiKS_ks_results + maintained_gene_Ks_values

    #return maintained_gene_Ks_values
    return demographiKS_plus_selection_results

if __name__ == '__main__':
    unittest.main()
