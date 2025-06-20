import unittest

from figure_generation.ks_plot_aggregations_dgks_vs_truth import make_Tc_Ks_vs_truth_fig_with_subplots


class DGKS_and_Empirical_Data_Test(unittest.TestCase):

    def test_maize_as_segmental_polyploid(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize'
        specks_out_path = 'foo'
        demographics_run_list = [False,
                                 'EMP_Mays_07_m06d17y2025_h09m34s26',
                                 'EMP_Mays_10_m06d17y2025_h10m38s26',
                                 'EMP_Mays_13_m06d17y2025_h15m46s38',
                                 'EMP_Mays_10_07_combined'
                                  ]
        specks_TE5_run_list = [False,False,False,False,False,False,False]


        xmax_Ks = [0.4 ,0.4 ,0.4,0.4,0.4 ]
        bin_sizes_Ks = [xmax_KS_i/25 for xmax_KS_i in xmax_Ks]

        #xmax_Tc = [5000 for f in demographics_run_list ]
        xmax_Tc = [10000,10000,10000,10000,10000]
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
        plot_title_lamda = lambda config: "Ks at Tnow\n" + "Na:" + str(config.ancestral_Ne)
        make_Tc_Ks_vs_truth_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                              demographiKS_out_path, demographics_run_list, run_list_name,
                                              specks_TE5_run_list, specks_out_path,
                                              xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                              suptitle, show_KS_predictions,
                                              include_annotation,
                                              plots_to_show_legend,
                                              plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here



    def test_coffee_allopolyploid(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee/Truth'
        demographics_run_list = [False,
                                 'EMP_Coff_17_m06d18y2025_h12m12s04',
                                 'EMP_Coff_16_m06d18y2025_h09m03s22',
                                 'Truth'
                                  ]# 'EMP_Coff_10_m06d12y2025_h17m06s38'
        truth_run_list = ['coffea.ks.tsv' for f in demographics_run_list ]


        xmax_Ks = [0.1 for f in demographics_run_list ]
        bin_sizes_Ks = [xmax_KS_i/50for xmax_KS_i in xmax_Ks]

        xmax_Tc = [40000 for f in demographics_run_list ]
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
        plot_title_lamda = lambda config: "Ks at Tnow\n" + "Na:" + str(config.ancestral_Ne)
        make_Tc_Ks_vs_truth_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                              demographiKS_out_path, demographics_run_list, run_list_name,
                                              truth_run_list, truth_out_path,
                                              xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                              suptitle, show_KS_predictions,
                                              include_annotation,
                                              plots_to_show_legend,
                                              plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here




if __name__ == '__main__':
    unittest.main()
