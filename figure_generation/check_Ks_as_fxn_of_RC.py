import unittest

from figure_generation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots


class TestKsByNe(unittest.TestCase):

    def test_Ks_for_varying_RC_1KNe_At_like(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_RC/save_Ne_1K'
        specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_RC/save_Ne_1K'

        run_list_name = "Ks_for_varying_varying_RC_At_like"
        demographics_run_list = [False,
                                 'KSvsRC7_At2_m03d05y2025_h12m12s14',
                                  'KSvsRC8_At2_m03d05y2025_h10m52s14',
                                 'KSvsRC9_At2_m03d05y2025_h14m03s28',
                                 'KSvsRC10_At2_m03d05y2025_h13m48s18'
                                 ]

        specks_TE5_run_list = [False, False, False, False, False]

        #bin_sizes_Tc = [200, 200, 200, 200, 200, 200, 200]
        #bin_sizes_Ks = [0.002, 0.002, 0.002,0.002, 0.002, 0.002, 0.002]
        xmax_Ks = [0.02 for f in demographics_run_list] #[0.025, 0.025, 0.025, 0.025, 0.025]  # 0.001  # max(demographiKS_ks_results)
        bin_sizes_Ks = [xmax_Ks_i / 40.0 for xmax_Ks_i in xmax_Ks]
        xmax_Tc = [10000 for f in demographics_run_list]
        #ymax_KS = [160 for f in demographics_run_list]
        ymax_KS = [150 for f in demographics_run_list]
        #ymax_Tc = [100 for f in demographics_run_list]
        ymax_Tc = [150 for f in demographics_run_list]
        bin_sizes_Tc = [xmax_Tc_i / 40.0 for xmax_Tc_i in xmax_Tc]
        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate

        show_KS_predictions = [False, False, False]
        suptitle = "SLiM and SpecKS Ks histograms\n"
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions)

        self.assertEqual(True, True)  # add assertion here


    def test_Ks_for_varying_RC_10KNe_At_like(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_RC/save_Ne_10K'
        specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_RC/save_Ne_10K'

        run_list_name = "Ks_for_varying_varying_RC_At_like_Fig R-RC1"
        #  'KSvsRC10_At10K_m03d07y2025_h15m56s54',
        demographics_run_list = [False,
                                 False,
                                 'KSvsRC9_At10K_m03d08y2025_h08m27s35',
                                 'KSvsRC8_At10K_m03d07y2025_h14m10s14',
                                 'KSvsRC7_At10K_m03d08y2025_h08m32s59',
                                 'KSvsRC6_At10K_m03d08y2025_h11m05s16'
                                 ]

        specks_TE5_run_list = [False, 'SpecKS_KSvsRC0_at_m03d08y2025_h15m59s40',
                               False, False,False,False]

        #bin_sizes_Tc = [200, 200, 200, 200, 200, 200, 200]
        #bin_sizes_Ks = [0.002, 0.002, 0.002,0.002, 0.002, 0.002, 0.002]
        xmax_Ks = [0.1 for f in demographics_run_list] #[0.025, 0.025, 0.025, 0.025, 0.025]  # 0.001  # max(demographiKS_ks_results)
        bin_sizes_Ks = [xmax_Ks_i / 40.0 for xmax_Ks_i in xmax_Ks]
        xmax_Tc = [80000 for f in demographics_run_list]
        #ymax_KS = [160 for f in demographics_run_list]
        ymax_KS = [250 for f in demographics_run_list]
        #ymax_Tc = [100 for f in demographics_run_list]
        ymax_Tc = [100 for f in demographics_run_list]
        bin_sizes_Tc = [xmax_Tc_i / 40.0 for xmax_Tc_i in xmax_Tc]
        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate

        which_plot_panels_to_show_legend = [1, 2, 3,4,5,6]
        show_KS_predictions = [False, False, False]
        show_Annotations = False
        suptitle = "DemographiKS and SpecKS Ks histograms\n"
        plot_title_lamda = lambda config: "Ks at Tnow\n" + "RC rate:" + str(config.recombination_rate)

        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions,
                                     show_Annotations,which_plot_panels_to_show_legend,
                                     plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here



    def test_Ks_for_varying_RC_10K_Ne(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_RC/save_Ne_10K'
        specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'

        #full, w/Ne 10K
        demographics_run_list = [False, 'KSvsRC6_m01d25y2025_h13m19s03',
                                 'KSvsRC7_m01d25y2025_h13m14s32', 'KSvsRC8_m01d25y2025_h13m13s04',
         'KSvsRC9_m01d25y2025_h13m09s54','KSvsRC10_m01d24y2025_h10m40s24_Ne_10_000']
        specks_TE5_run_list = [False,False,False,False,False, False]

        bin_sizes_Tc = [2000 for f in demographics_run_list]
        bin_sizes_Ks = [0.01 for f in demographics_run_list]
        xmax_Ks = [0.5 for f in demographics_run_list] #[0.025, 0.025, 0.025, 0.025, 0.025]  # 0.001  # max(demographiKS_ks_results)
        xmax_Tc = [80000 for f in demographics_run_list]
        ymax_KS = [100 for f in demographics_run_list]
        ymax_Tc = [100 for f in demographics_run_list]

        run_list_name = "Ks_for_varying_varying_RC_10KNe"
        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate

        which_plot_panels_to_show_legend = [1, 2, 3, 4]
        show_KS_predictions = [False, False, False]
        show_Annotations = False
        suptitle = "DemographiKS and SpecKS Ks histograms\n"
        plot_title_lamda = lambda config: "Ks at Tnow\n" + "RC:" + str(config.recombination_rate)

        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions,
                                     show_Annotations,
                                     which_plot_panels_to_show_legend,plot_title_lamda )

        self.assertEqual(True, True)  # add assertion here

    def test_Ks_for_varying_RC_1KNe(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_RC/save_Ne_1K'
        specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_RC/save_Ne_1K'


        #full, w/Ne1000
        #KSvsRC9_m01d21y2025_h14m24s03

        #demographics_run_list = [False, 'KSvsRC6_m01d21y2025_h20m06s34', 'KSvsRC7_m01d21y2025_h12m16s26',
        #                           'KSvsRC8_m01d23y2025_h10m39s34', 'KSvsRC9_m01d21y2025_h14m24s03',
        #                           'KSvsRC10_m01d24y2025_h09m16s53', 'KSvsRC11_m01d23y2025_h11m04s33']

        #ks_hist_by_Ks_for_varying_varying_RC_test_Ne1000_long_burnin_save_for_paper.png
        demographics_run_list = [False, 'KSvsRC6_m01d21y2025_h20m06s34']
                                 #,'KSvsRC7_m01d21y2025_h12m16s26',)
                                 #'KSvsRC8_m01d23y2025_h10m39s34','KSvsRC9moreBI_m03d03y2025_h14m27s22',
                                 #'KSvsRC10_m01d24y2025_h10m40s24']
        #these ones Tc sucks and need to be re-run with more burnin
        #                         #,'KSvsRC9_m01d24y2025_h08m51s36',
        #                         #'KSvsRC10_m01d24y2025_h09m16s53','KSvsRC11_m01d23y2025_h11m04s33']
        specks_TE5_run_list = [False, False, False, False, False, 'SpecKS_KSvsRC0', False]

        run_list_name = "Ks_for_varying_varying_RC_foo"
        #demographics_run_list = [False,
        #                         'KSvsRC7_At2_m03d05y2025_h10m52s10',
        #                          'KSvsRC8_At2_m03d05y2025_h10m52s14',
        #                          'KSvsRC9_At2_m03d05y2025_h10m52s17',
        #                          'KSvsRC10_At2_m03d05y2025_h10m52s20'
        #                         ]

        specks_TE5_run_list = [False, False, False, False, False, False, False]

        #bin_sizes_Tc = [200, 200, 200, 200, 200, 200, 200]
        #bin_sizes_Ks = [0.002, 0.002, 0.002,0.002, 0.002, 0.002, 0.002]
        xmax_Ks = [0.02 for f in demographics_run_list] #[0.025, 0.025, 0.025, 0.025, 0.025]  # 0.001  # max(demographiKS_ks_results)
        bin_sizes_Ks = [xmax_Ks_i / 50.0 for xmax_Ks_i in xmax_Ks]
        xmax_Tc = [50000 for f in demographics_run_list]
        #ymax_KS = [160 for f in demographics_run_list]
        ymax_KS = [200 for f in demographics_run_list]
        #ymax_Tc = [100 for f in demographics_run_list]
        ymax_Tc = [200 for f in demographics_run_list]
        bin_sizes_Tc = [xmax_Tc_i / 50.0 for xmax_Tc_i in xmax_Tc]
        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate

        show_KS_predictions = [False, False, False]
        suptitle = "SLiM and SpecKS Ks histograms\n"
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions)

        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
