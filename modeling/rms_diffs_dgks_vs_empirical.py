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
from modeling import ks_parsers, birth_and_death_with_escape


class Final_DGKS_vs_Empirical(unittest.TestCase):


    def test_poplar_final_plots(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Poplar'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Poplar/Truth'
        leaf_to_leaf_Ks=2.0
        include_selection = True
        include_homEx = True

        full_results = ['EMP_Pop_25_m07d08y2026_h11m17s37','EMP_Pop_26_m07d08y2026_h11m17s39']

        if include_homEx:
            species_for_plot_title = full_results
        else:
            species_for_plot_title = [full_results[0]]

        species_for_plot_title_str = "_".join(species_for_plot_title)
        #species_for_plot_title = 'EMP_Pop_07_and_half_11'
        out_png = os.path.join(demographiKS_out_path, species_for_plot_title_str + '_final_poplar_after_fix.png')
        real_full_path = os.path.join(truth_out_path, 'poplar.ks.tsv')
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)


        demographiKS_ks_results=[]
        for sim_run in species_for_plot_title:
            sim_full_path = os.path.join(demographiKS_out_path,
                                     sim_run,
                                     'allotetraploid_bottleneck.csv')
            run_results = read_Ks_csv(sim_full_path, False)
            run_results = [k * leaf_to_leaf_Ks for k in run_results]
            demographiKS_ks_results=demographiKS_ks_results+run_results

        print("Num genes w/out maintained SSEs: " + str(len(demographiKS_ks_results)))
        if include_selection:
            seed=44
            fraction_needed=1.0
            demographiKS_ks_results = add_selection(demographiKS_ks_results,
                                                    fraction_needed, seed)

        print("Num paralogs in final genome: " + str(len(demographiKS_ks_results)))
        bin_size = 0.01
        max_Ks = 0.9
        bins = np.arange(bin_size, max_Ks + 0.1, bin_size)
        hist_ys_real, bins_real, patches = plt.hist(real_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                                    density=True)

        hist_ys_sim, bins_sim, patches = plt.hist(demographiKS_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                                  density=True)

        hist_ys_1 = [hist_ys_real, bins_real]
        hist_ys_2 = [hist_ys_sim, bins_sim]
        list_of_hist_data = [hist_ys_1, hist_ys_2]

        overlay_differences_in_curves(species_for_plot_title, list_of_hist_data,
                                      include_selection, include_homEx, out_png)

        self.assertEqual(True, True)  # add assertion here

    def test_maize_final_plots(self):
        # linux
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize/Truth'
        include_selection = True
        include_homEx =  True
        leaf_to_leaf_Ks=2.0
        fraction_needed = 4#'','EMP_Mays_38_m07d06y2026_h16m52s46'
        full_results = ['EMP_Mays_80_m07d13y2026_h10m40s15',
                        'EMP_Mays_79_m07d13y2026_h10m40s11']

        if include_homEx:
            species_for_plot_title = full_results
        else:
            species_for_plot_title= [full_results[0]]

        species_for_plot_title_str = "_".join(species_for_plot_title)

        #species_for_plot_title = 'EMP_Mays_26_29_combined'
        out_png = os.path.join(demographiKS_out_path, species_for_plot_title_str + '_final_mays_v2.png')
        real_full_path = os.path.join(truth_out_path, 'mays.ks.tsv')
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        demographiKS_ks_results=[]
        for sim_run in species_for_plot_title:
            sim_full_path = os.path.join(demographiKS_out_path,
                                     sim_run,
                                     'allotetraploid_bottleneck.csv')
            run_results = read_Ks_csv(sim_full_path, False)
            run_results = [k * leaf_to_leaf_Ks for k in run_results]
            demographiKS_ks_results=demographiKS_ks_results+run_results


        if include_selection:
            seed=42
            max_ks=2,

            demographiKS_ks_results = add_selection(demographiKS_ks_results,
                                                    fraction_needed, seed)

        print("Num paralogs in final genome: " + str(len(demographiKS_ks_results)))
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
                                      include_selection,include_homEx, out_png)

        self.assertEqual(True, True)  # add assertion here


    def test_coffee_final_plots(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee/Truth'
        include_selection=True
        include_homEx=True
        leaf_to_leaf_Ks=2.0

        #max_Ks = 0.4
        #fraction_genes_maintained = 1.0

        max_Ks = 0.2
        fraction_genes_maintained = 2.5
        rseed = 42

        species_for_plot_title =[
            'EMP_Coff_73_m07d12y2026_h09m13s58',
            'EMP_Coff_74_m07d12y2026_h09m14s21']

        species_for_plot_title_str = "_".join(species_for_plot_title)
        out_png = os.path.join(demographiKS_out_path,species_for_plot_title_str
        + "_fraction" + str(fraction_genes_maintained) + "_"
                               + '_final_coffee_vKSrates_fix.png')
        real_full_path = os.path.join(truth_out_path, 'coffea.ks.tsv')
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        sim_full_path = os.path.join(demographiKS_out_path,
                                     species_for_plot_title[0],
                                     'allotetraploid_bottleneck.csv')
        #demographiKS_ks_results = read_Ks_csv(sim_full_path, False)


        demographiKS_ks_results=[]
        for sim_run in species_for_plot_title:
            sim_full_path = os.path.join(demographiKS_out_path,
                                     sim_run,
                                     'allotetraploid_bottleneck.csv')
            run_results = read_Ks_csv(sim_full_path, False)
            run_results = [k * leaf_to_leaf_Ks for k in run_results]
            demographiKS_ks_results=demographiKS_ks_results+run_results

        if include_selection:

            demographiKS_ks_results = add_selection(demographiKS_ks_results,
                                                    fraction_genes_maintained, rseed)


        print("Num paralogs in final genome: " + str(len(demographiKS_ks_results)))
        bin_size = 0.005
        #bin_size = 0.01 too smooth
        bins = np.arange(bin_size, max_Ks + 0.1, bin_size)
        hist_ys_real, bins_real, patches = plt.hist(real_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                                    density=True)

        hist_ys_sim, bins_sim, patches = plt.hist(demographiKS_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                                  density=True)

        hist_ys_1 = [hist_ys_real, bins_real]
        hist_ys_2 = [hist_ys_sim, bins_sim]
        list_of_hist_data = [hist_ys_1, hist_ys_2]


        overlay_differences_in_curves(species_for_plot_title, list_of_hist_data,
                                      include_selection, include_homEx,out_png)

        self.assertEqual(True, True)  # add assertion here

    def test_kiwi_final_plots(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Kiwi'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Kiwi/Truth'
        include_homEx =   True
        include_selection =  True
        proportion_retained_genes = 3.0
        seed= 42
        max_Ks = 0.5#0.5#1.0
        num_bins=100
        leaf_to_leaf_Ks=2.0
        bin_size = 0.01# max_Ks/num_bins

        full_results = ['EMP_Act_56wgd_m07d09y2026_h15m27s28',
                        'EMP_Act_74hx_m07d13y2026_h10m47s09',
                        'EMP_Act_48_m07d09y2026_h08m50s20']

        if include_homEx:
            species_for_plot_title = full_results
        else:
            species_for_plot_title = [full_results[0]]

        species_for_plot_title_string = "_".join(species_for_plot_title)
        out_png = os.path.join(demographiKS_out_path, species_for_plot_title_string
                               + "MaxKs" + str(max_Ks) +  '_final_kiwi_v2.png')
        real_full_path = os.path.join(truth_out_path, 'actinidia_AvB.ks.tsv')
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        demographiKS_ks_results = []
        for sim_run in species_for_plot_title:
            sim_full_path = os.path.join(demographiKS_out_path,
                                         sim_run,
                                         'allotetraploid_bottleneck.csv')
            run_results = read_Ks_csv(sim_full_path, False)
            run_results = [k * leaf_to_leaf_Ks for k in run_results]
            demographiKS_ks_results = demographiKS_ks_results + run_results

        if include_selection:
            demographiKS_ks_results = add_selection(demographiKS_ks_results,
                                                    proportion_retained_genes,
                                                    seed)

        print("Num paralogs in final genome: " + str(len(demographiKS_ks_results)))
        bins = np.arange(bin_size, max_Ks, bin_size)
        hist_ys_real, bins_real, patches = plt.hist(real_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                                    density=True)

        hist_ys_sim, bins_sim, patches = plt.hist(demographiKS_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                                  density=True)

        hist_ys_1 = [hist_ys_real, bins_real]
        hist_ys_2 = [hist_ys_sim, bins_sim]
        list_of_hist_data = [hist_ys_1, hist_ys_2]

        overlay_differences_in_curves(species_for_plot_title, list_of_hist_data,
                                      include_selection, include_homEx, out_png)

        self.assertEqual(True, True)  # add assertion here


def overlay_differences_in_curves(species_for_plot_title, list_of_hist_data,
                                  include_selection,include_homEx,
                                  out_png):

    colors = ['green','blue','brown']
    labels = ['truth','model fit']
    if include_selection:
        out_png = out_png.replace(".png","_with_selection.png")
    if include_homEx:
        out_png = out_png.replace(".png","_and_homEx.png")

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    for i in range(0,len(list_of_hist_data)):
        hist_data =list_of_hist_data[i]
        [n, bins]=hist_data
        width=(bins[2]-bins[1])/2.0
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

    plt.legend(fontsize="15")
    plt.xlabel("Ks")
    ax.set(ylabel="Density")
    #plt.title(species_for_plot_title, style='italic')
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


def add_selection(demographiKS_ks_results,fraction_genes_maintained,rseed):

    step_size = 0.001
    decay_constant = 33  # =(1/mean life expectancy of SSD gene), in Ks space
    rate_escape = 0.3  # % of SSD are maintained

    # where decay constant comes from:
    #
    #    #(Ferris and Whitt, 1977; Nadeau and Sankoff, 1997; Lynch and Conery, 2000;
    #    # #Blanc and Wolfe, 2004b; Soltis et al.,2016; Cheng et al., 2018).
    #    # if a SSD lasts a few million years (say 3x10^6)
    #    # In Ks space, that depends mutation rate.
    #    # if Mutation rate is 10^-8, then
    #    # SSD lasts a few million years (say 3x10^-2 ) in KS
    #    # so, decay constant is 1/(3x10^-2)

    # where rate escape comes from
    #  Proportion of retained duplicates (what portion of the overall genome is retained duplicates):
    #  between 0.3, 0.4  human and mouse https://pmc.ncbi.nlm.nih.gov/articles/PMC1413713/
    #  Not the same as plants, but at least we are sure those numbers not confounded w/ WGD
    
    # Retention of duplicated genes in evolution (20-40% of the whole genome)
    # so, a genome of size 100 has 20-40 retained duplcatate SSDs
    # if they were born randomly over 574 MY (~KS=4), that's 0.035-0.070 gene retained per MY
    # from a starting genome of size 100 genes.
    # So, thats (0.035 - 0.07)/100, so 0.00035 to 0.0007 dup retained per gene million year.
    # Gene birth rate is ~0.00162 per g per million years (Tiley)
    # If 0.00162 are born, and 0.00035 to 0.0007 survive, then
    # ~0.21-0.43 survive. Wow. more than I expected.
    # refs ()

    random.seed(rseed)
    include_debugging_plots = True
    my_pdf, xs = birth_and_death_with_escape.gene_birth_death_with_escape_pdf(
        step_size, decay_constant, rate_escape, rseed, include_debugging_plots)

    num_genes_needed = int(fraction_genes_maintained*len(demographiKS_ks_results))
    SSD_ks = birth_and_death_with_escape.draw_SSDs_from_pdf(my_pdf, xs,
                                                            num_genes_needed,
                                                            rseed, include_debugging_plots)

    print("Num SSD genes simulated: " + str(len(SSD_ks)))

    demographiKS_plus_selection_results = demographiKS_ks_results + SSD_ks

    return demographiKS_plus_selection_results

if __name__ == '__main__':
    unittest.main()
