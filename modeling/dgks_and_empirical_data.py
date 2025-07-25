import math
import os
import unittest

import numpy as np
from matplotlib.patches import Rectangle
from matplotlib import pyplot as plt

from figure_generation.histogram_plotter import read_Ks_csv
from figure_generation.ks_plot_aggregations_dgks_vs_truth import make_Tc_Ks_vs_truth_fig_with_subplots
from modeling import ks_parsers


class DGKS_and_Empirical_Data_Test(unittest.TestCase):

    def test_maize_as_segmental_polyploid(self):

        # linux
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize/Truth'

        # mac
        # demographiKS_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize'
        # truth_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize/Truth'

        demographics_run_list = [False,
                                 'EMP_Mays_26_m07d09y2025_h11m39s07',
                                 'EMP_Mays_28_m07d14y2025_h17m09s16',
                                 'EMP_Mays_29_m07d14y2025_h17m35s58',
                                 'EMP_Mays_26_29_combined'
                                 ]

        truth_run_list = ['mays.ks.tsv' for f in demographics_run_list ]

        xmax_Ks = [0.4 for f in demographics_run_list ]
        bin_sizes_Ks = [xmax_KS_i/25 for xmax_KS_i in xmax_Ks]

        xmax_Tc = [10000  for f in demographics_run_list]
        bin_sizes_Tc = [xmax_Tc_i/25 for xmax_Tc_i in xmax_Tc]


        ymax_KS = [False for f in demographics_run_list]
        ymax_Tc = [False for f in demographics_run_list]
        run_list_name = "Simulated_Ks_for_Zea_mays"
        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate

        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation = False
        plots_to_show_legend = [1, 2, 3, 4]
        plot_title_lamda = lambda config: "Na:" + str(config.ancestral_Ne)
        make_Tc_Ks_vs_truth_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                              demographiKS_out_path, demographics_run_list, run_list_name,
                                              truth_run_list, truth_out_path,
                                              xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                              suptitle, show_KS_predictions,
                                              include_annotation,
                                              plots_to_show_legend,
                                              plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here

    def test_maize_as_segmental_polyploid_best(self):

        #linux
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize/Truth'
        out_png= os.path.join(demographiKS_out_path, 'final.png')

        species_for_plot_title =     'EMP_Mays_26_29_combined'
        real_full_path =  os.path.join(truth_out_path, 'mays.ks.tsv')
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        sim_full_path =  os.path.join(demographiKS_out_path,
                                              species_for_plot_title,
                                              'allotetraploid_bottleneck.csv')
        demographiKS_ks_results = read_Ks_csv(sim_full_path,False)
        bin_size=0.01
        max_Ks=0.4
        bins = np.arange(bin_size, max_Ks + 0.1, bin_size)
        hist_ys_real, bins_real, patches = plt.hist(real_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                          density=True)

        hist_ys_sim, bins_sim, patches = plt.hist(demographiKS_ks_results, bins=bins, facecolor='b', alpha=0.25,
                                          density=True)

        hist_ys_1 = [hist_ys_real, bins_real]
        hist_ys_2 = [hist_ys_sim, bins_sim]
        list_of_hist_data = [hist_ys_1, hist_ys_2 ]

        overlay_differences_in_curves(species_for_plot_title, list_of_hist_data, out_png)


        self.assertEqual(True, True)  # add assertion here

    def test_coffee_allopolyploid(self):

        # linux
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee/Truth'

        # mac
        # demographiKS_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee'
        # truth_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee/Truth'

        demographics_run_list = [False,
                                 'EMP_Coff_33_m07d01y2025_h08m59s09',
                                 'EMP_Coff_34_m07d01y2025_h08m58s56',
                                 'EMP_Coff_35_m07d01y2025_h08m58s51',
                                 'EMP_Coff_36_m07d09y2025_h12m22s40',
                                 'EMP_Coff_37_m07d16y2025_h11m12s14'
                                  ]#     'EMP_Coff_36_m07d09y2025_h12m22s40',
        #                        'EMP_Coff_35_m06d30y2025_h17m25s25',
        #                                         'EMP_Coff_30_m06d26y2025_h17m18s32',
        # 'EMP_Coff_29_m06d26y2025_h12m21s24','EMP_Coff_10_m06d12y2025_h17m06s38',  'EMP_Coff_28_m06d26y2025_h12m21s27',
        # 'EMP_Coff_24_m06d24y2025_h13m57s43', 'EMP_Coff_25_m06d25y2025_h09m29s18',                                 'EMP_Coff_26_m06d25y2025_h14m39s51',
        #                                  'EMP_Coff_27_m06d25y2025_h14m45s12',
        # 'EMP_Coff_17_m06d18y2025_h12m12s04','EMP_Coff_20_m06d23y2025_h12m33s54','EMP_Coff_21_m06d23y2025_h12m42s58',
        #'EMP_Coff_16_m06d18y2025_h09m03s22','EMP_Coff_18_m06d19y2025_h17m06s11', 'EMP_Coff_19_m06d20y2025_h10m02s34',
        truth_run_list = ['coffea.ks.tsv' for f in demographics_run_list ]


        xmax_Ks = [0.2 for f in demographics_run_list ]
        bin_sizes_Ks = [xmax_KS_i/50for xmax_KS_i in xmax_Ks]

        xmax_Tc = [400000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/50 for xmax_Tc_i in xmax_Tc]


        ymax_KS = [False for f in demographics_run_list]
        ymax_Tc = [False for f in demographics_run_list]
        run_list_name = "Simulated_Ks_for_Coffea_arabica"
        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate

        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation = False
        plots_to_show_legend = [1, 2, 3, 4]
        plot_title_lamda = lambda config: "Na:" + str(config.ancestral_Ne)
        make_Tc_Ks_vs_truth_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                              demographiKS_out_path, demographics_run_list, run_list_name,
                                              truth_run_list, truth_out_path,
                                              xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                              suptitle, show_KS_predictions,
                                              include_annotation,
                                              plots_to_show_legend,
                                              plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here

    def test_poplar_allopolyploid(self):

        # linux
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Poplar'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Poplar/Truth'

        # mac
        # demographiKS_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Poplar'
        # truth_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Poplar/Truth'

        demographics_run_list = [False,
                                 'EMP_Pop_07_m07d17y2025_h16m44s40',
                                 'EMP_Pop_09_m07d18y2025_h11m06s20',
                                 'EMP_Pop_11_m07d18y2025_h13m55s45',
                                 'EMP_Pop_07_and_half_11',
                                 'EMP_Pop_13_m07d22y2025_h13m24s05']
        #'EMP_Pop_02_m07d14y2025_h11m32s24','EMP_Pop_06_m07d17y2025_h11m35s43',
        #'EMP_Pop_04_m07d14y2025_h17m21s10','EMP_Pop_10_m07d18y2025_h11m40s57',

        truth_run_list = ['poplar.ks.tsv' for f in demographics_run_list ]


        xmax_Ks = [0.5 for f in demographics_run_list ]
        bin_sizes_Ks = [xmax_KS_i/50for xmax_KS_i in xmax_Ks]

        xmax_Tc = [80000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/50 for xmax_Tc_i in xmax_Tc]


        ymax_KS = [False for f in demographics_run_list]
        ymax_Tc = [False for f in demographics_run_list]
        run_list_name = "Simulated_Ks_for_Poplar_tricocarpa"
        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate

        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation = False
        plots_to_show_legend = [1, 2, 3, 4]
        plot_title_lamda = lambda config: "Na:" + str(config.ancestral_Ne)
        make_Tc_Ks_vs_truth_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                              demographiKS_out_path, demographics_run_list, run_list_name,
                                              truth_run_list, truth_out_path,
                                              xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                              suptitle, show_KS_predictions,
                                              include_annotation,
                                              plots_to_show_legend,
                                              plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here


    def test_rice_allopolyploid(self):

        # linux
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Rice'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Rice/Truth'

        # mac
        # demographiKS_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Rice'
        # truth_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Rice/Truth'

        demographics_run_list = [False,
                                 'EMP_Ory_06_m07d17y2025_h12m08s40',
                                 'EMP_Ory_02_m07d17y2025_h16m24s51',
                                  ]

        truth_run_list = ['sativa.ks.tsv' for f in demographics_run_list ]


        xmax_Ks = [4.0 for f in demographics_run_list ]
        bin_sizes_Ks = [xmax_KS_i/50for xmax_KS_i in xmax_Ks]

        xmax_Tc = [80000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/50 for xmax_Tc_i in xmax_Tc]


        ymax_KS = [False for f in demographics_run_list]
        ymax_Tc = [False for f in demographics_run_list]
        run_list_name = "Simulated_Ks_for_Oryza_sativa"
        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate

        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation = False
        plots_to_show_legend = [1, 2, 3, 4]
        plot_title_lamda = lambda config: "Na:" + str(config.ancestral_Ne)
        make_Tc_Ks_vs_truth_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                              demographiKS_out_path, demographics_run_list, run_list_name,
                                              truth_run_list, truth_out_path,
                                              xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                              suptitle, show_KS_predictions,
                                              include_annotation,
                                              plots_to_show_legend,
                                              plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here


    def test_sugarcane_allopolyploid(self):

        # linux
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Sugarcane'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Sugarcane/Truth'

        # mac
        # demographiKS_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Rice'
        # truth_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Rice/Truth'

        demographics_run_list = [False,
                                 'EMP_Sac_01_m07d24y2025_h14m01s54',
                                  ]

        truth_run_list = ['saccharum.ks.tsv' for f in demographics_run_list ]


        xmax_Ks = [0.5 for f in demographics_run_list ]
        bin_sizes_Ks = [xmax_KS_i/100 for xmax_KS_i in xmax_Ks]

        xmax_Tc = [80000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/100 for xmax_Tc_i in xmax_Tc]


        ymax_KS = [8 for f in demographics_run_list]
        ymax_Tc = [False for f in demographics_run_list]
        run_list_name = "Simulated_Ks_for_Saccharum_spontaneum"
        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate

        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation = False
        plots_to_show_legend = [1, 2, 3, 4]
        plot_title_lamda = lambda config: "Na:" + str(config.ancestral_Ne)
        make_Tc_Ks_vs_truth_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                              demographiKS_out_path, demographics_run_list, run_list_name,
                                              truth_run_list, truth_out_path,
                                              xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                              suptitle, show_KS_predictions,
                                              include_annotation,
                                              plots_to_show_legend,
                                              plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here


def overlay_differences_in_curves(species_for_plot_title, list_of_hist_data, out_png):

    colors = ['green','blue','brown']
    labels = ['truth','model fit']


    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    for i in range(0,len(list_of_hist_data)):
        hist_data =list_of_hist_data[i]
        [n, bins]=hist_data
        width=(bins[2]-bins[1])
        bar_plot_xs=[b + i*width for b in bins[0:len(bins)-1]] #to match bar-plot axes
        sub_bins=bins[0:len(bins) - 1]
        label=labels[i]
        #if label=='truth':

            #x_value_of_ymax = get_Ks_at_max(WGD_ks, bar_plot_xs, n, sub_bins)
            #label = label +" (peak at Ks=" +str(WGD_ks) + ")"
            #plt.axvline(WGD_ks, color=config.color_blind_friendly_color_cycle_analogs['green'],
            #            linestyle=':')
        plt.plot(bar_plot_xs,n,c=colors[i], alpha=1, label=label)
        #plt.plot(sub_bins, n, c=colors[i], alpha=1, label=label)

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
