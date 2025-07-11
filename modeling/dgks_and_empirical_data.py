import unittest

from figure_generation.ks_plot_aggregations_dgks_vs_truth import make_Tc_Ks_vs_truth_fig_with_subplots


class DGKS_and_Empirical_Data_Test(unittest.TestCase):

    def test_maize_as_segmental_polyploid(self):

        # linux
        #demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize'
        #truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize/Truth'

        # mac
        demographiKS_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize'
        truth_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize/Truth'

        demographics_run_list = [False,
                                 'EMP_Mays_09_m06d17y2025_h10m05s25',
                                 'EMP_Mays_23_m07d03y2025_h16m18s38',
                                 'EMP_Mays_24_m07d03y2025_h16m22s57',
                                 'EMP_Mays_25_m07d09y2025_h11m35s37',
                                 'EMP_Mays_26_m07d09y2025_h11m39s07',
                                 'EMP_Mays_26_25_combined'
                                 ]
        #                         'EMP_Mays_22_combined',                                 'EMP_Mays_09_m06d17y2025_h10m05s25',
        #                                  'EMP_Mays_21_m07d02y2025_h11m32s40',
        #                                  'EMP_Mays_22_m07d02y2025_h12m02s50',
        #                        ]
                                 #'EMP_Mays_20_m07d02y2025_h10m54s44',

        #'EMP_Mays_18_m07d02y2025_h10m19s38','EMP_Mays_07_m06d17y2025_h09m34s26',


        #                         'EMP_Mays_10_m06d17y2025_h10m38s26',      'EMP_Mays_19_m07d02y2025_h10m20s06',
        #                         'EMP_Mays_16_m07d01y2025_h15m13s23',                              'EMP_Mays_06_m06d16y2025_h19m17s43',
        #                         'EMP_Mays_17_m07d01y2025_h15m55s20']
                                 #'EMP_Mays_13_m06d17y2025_h15m46s38',
                                 #'EMP_Mays_13_m06d18y2025_h12m23s12',
                                 #'EMP_Mays_15_m06d20y2025_h11m24s06', 'EMP_Mays_05_m06d16y2025_h19m16s05',
                                 #]
                                 #'EMP_Mays_10_07_combined',
                                 #'EMP_Mays_15_07_combined'
                                 # ]
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

    def test_maize_as_segmental_polyploid_2(self):

        #linux
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize'
        truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Maize/Truth'


        demographics_run_list = [False,
                                 'EMP_Mays_22_combined']

        truth_run_list = ['mays.ks.tsv' for f in demographics_run_list]

        xmax_Ks = [0.4 for f in demographics_run_list]
        bin_sizes_Ks = [xmax_KS_i / 50 for xmax_KS_i in xmax_Ks]

        xmax_Tc = [10000 for f in demographics_run_list]
        bin_sizes_Tc = [xmax_Tc_i / 50 for xmax_Tc_i in xmax_Tc]

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

    def test_coffee_allopolyploid(self):

        # linux
        # demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee'
        # truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee/Truth'

        # mac
        demographiKS_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee'
        truth_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Coffee/Truth'

        demographics_run_list = [False,
                                 'EMP_Coff_31_m06d26y2025_h17m18s28',
                                 'EMP_Coff_32_m06d27y2025_h11m19s33',
                                 'EMP_Coff_33_m07d01y2025_h08m59s09',
                                 'EMP_Coff_34_m07d01y2025_h08m58s56',
                                 'EMP_Coff_35_m07d01y2025_h08m58s51'
                                  ]
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
        # demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Poplar'
        # truth_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Poplar/Truth'

        # mac
        demographiKS_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Poplar'
        truth_out_path = '/Users/tamsen/Data/DemographiKS_output_from_mesx/EmpiricalDataTesting_2/Poplar/Truth'

        demographics_run_list = [False,
                                 'EMP_Pop_01_m07d10y2025_h14m17s13',
                                  ]

        truth_run_list = ['poplar.ks.tsv' for f in demographics_run_list ]


        xmax_Ks = [0.5 for f in demographics_run_list ]
        bin_sizes_Ks = [xmax_KS_i/50for xmax_KS_i in xmax_Ks]

        xmax_Tc = [400000 for f in demographics_run_list ]
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





if __name__ == '__main__':
    unittest.main()
