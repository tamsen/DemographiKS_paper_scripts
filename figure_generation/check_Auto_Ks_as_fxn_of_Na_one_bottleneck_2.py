import unittest

from figure_generation.ks_plot_aggregations_auto_vs_allo import make_Tc_Ks_Allo_vs_Auto_fig_with_subplots


class TestKsForAuto(unittest.TestCase):

    #row #4, youn WGD
    def test_Ks_for_varying_Na_one_bottleneck_Nb_fixed_at_10K_v2(self):

        #full, w/Ne 10K
        demographiKS_allo_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/NaVaries/NbFixedQ1000'
        demographics_auto_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/NaVaries/NbFixedQ1000'

        demographics_allo_run_list = [
            False,
            'Allo_10KNb_10Na_v4_m04d08y2026_h17m18s47',
            'Allo_10KNb_50Na_v4_m04d08y2026_h17m18s47',
            'Allo_10KNb_100Na_v4_m04d08y2026_h17m18s47',
            'Allo_10KNb_500Na_v4_m04d09y2026_h10m44s19',
            'Allo_10KNa_500Nb_old_v4_m04d09y2026_h10m44s42' ]#'Allo_10KNb_1000Na_v4_m04d09y2026_h10m44s04'
             #  ]




        demographics_auto_run_list = [
            False,
            'Auto_10KNa_10Nb_v4_m04d08y2026_h17m10s00',
            'Auto_10KNa_50Nb_v4_m04d08y2026_h17m10s00',
            'Auto_10KNa_100Nb_v4_m04d08y2026_h17m10s00',
            'Auto_10KNb_500Na_v4_m04d09y2026_h10m44s21',
            'Auto_10KNa_500Nb_old_v4_m04d09y2026_h10m44s51' ,]#'Auto_10KNb_1000Na_v4_m04d09y2026_h10m44s03'
        #]

        # option A
        #'Auto_10KNa_500Nb_old_v4_m04d08y2026_h17m10s00'
        #'Allo_10KNb_500Na_old_v4_m04d08y2026_h17m18s47'

        # option B
        #Allo_10KNb_1000Na_v4_m04d09y2026_h10m17s10
        #Auto_10KNb_1000Na_v4_m04d09y2026_h10m17s13

        #xmax_Ks = [0.8  for f in demographics_auto_run_list ]
        #xmax_Ks_array = [[0.05, 0.05,0.05,0.05,0.05 ],[0.4, 0.4,0.5,0.6,0.8 ]]
        xmax_Ks_array = [[0.05 for f in demographics_auto_run_list ],
                         [0.05 for f in demographics_auto_run_list ]]
        bin_sizes_Ks_array  = [[xmax_KS_i/25 for xmax_KS_i in xmax_Ks_array[0]],
                        [xmax_KS_i / 25 for xmax_KS_i in xmax_Ks_array[1]]]

        xmax_Tc = [4000 for f in demographics_auto_run_list]
        bin_sizes_Tc = [xmax_Tc_i / 50 for xmax_Tc_i in xmax_Tc]

        ymax_Ks_array = [[200 for f in demographics_auto_run_list],
                   [200 for f in demographics_auto_run_list]]

        ymax_Tc = [False for f in demographics_auto_run_list]
        run_list_name = "Ks_for_Allo_and_Auto_varying_varying_Na_10K_Nb_Fig_R-NaNb8_v2"

        show_KS_predictions = [False, False, False]
        suptitle = "Auto and Allo Ks histograms\n"

        plot_title_lamda = lambda config: "Ancestral pop size:" + str(config.ancestral_Ne)
        include_annotation = False
        num_plot_rows = 2

        make_Tc_Ks_Allo_vs_Auto_fig_with_subplots(num_plot_rows,bin_sizes_Ks_array, bin_sizes_Tc,
        demographiKS_allo_out_path, demographics_allo_run_list, run_list_name,
        demographics_auto_run_list, demographics_auto_out_path,
        xmax_Ks_array, xmax_Tc, ymax_Ks_array, ymax_Tc,
        suptitle, show_KS_predictions,include_annotation,plot_title_lamda)

        self.assertEqual(True, True)


    def test_Ks_for_varying_Na_one_bottleneck_Nb_fixed_at_olderWGD_v2(self):

        #full, w/Ne 10K
        demographiKS_allo_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/NaVaries/NaVariesNbfixed_OldWGD'
        demographics_auto_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/NaVaries/NaVariesNbfixed_OldWGD'

        demographics_allo_run_list = [
            False,
            'Allo_10KNb_10Na_v5_m04d09y2026_h12m04s09',
            'Allo_10KNb_50Na_v5_m04d09y2026_h12m01s33',
            'Allo_10KNb_100Na_v5_m04d09y2026_h12m01s33',
            'Allo_10KNb_500Na_v5_m04d09y2026_h12m01s33',
            'Allo_10KNb_1000Na_v5_m04d09y2026_h12m01s35',
               ]




        demographics_auto_run_list = [
            False,
            'Auto_10KNa_10Nb_v5_m04d09y2026_h12m10s10',
            'Auto_10KNa_50Nb_v5_m04d09y2026_h12m10s07',
            'Auto_10KNa_100Nb_v5_m04d09y2026_h12m09s56',
            'Auto_10KNb_500Na_v5_m04d09y2026_h12m09s59',
            'Auto_10KNb_1000Na_v5_m04d09y2026_h12m13s02',
            ]


        #xmax_Ks = [0.8  for f in demographics_auto_run_list ]
        #xmax_Ks_array = [[0.05, 0.05,0.05,0.05,0.05 ],[0.4, 0.4,0.5,0.6,0.8 ]]
        xmax_Ks_array = [[0.05 for f in demographics_auto_run_list ],
                         [0.05 for f in demographics_auto_run_list ]]
        bin_sizes_Ks_array  = [[xmax_KS_i/25 for xmax_KS_i in xmax_Ks_array[0]],
                        [xmax_KS_i / 25 for xmax_KS_i in xmax_Ks_array[1]]]

        xmax_Tc = [4000 for f in demographics_auto_run_list]
        bin_sizes_Tc = [xmax_Tc_i / 50 for xmax_Tc_i in xmax_Tc]

        ymax_Ks_array = [[200 for f in demographics_auto_run_list],
                   [200 for f in demographics_auto_run_list]]

        ymax_Tc = [False for f in demographics_auto_run_list]
        run_list_name = "Ks_for_Allo_and_Auto_varying_varying_Na_10K_Nb_Fig_R-NaNb8_olderWGD_v2"

        show_KS_predictions = [False, False, False]
        suptitle = "Auto and Allo Ks histograms\n"

        plot_title_lamda = lambda config: "Ancestral pop size:" + str(config.ancestral_Ne)
        include_annotation = False
        num_plot_rows = 2

        make_Tc_Ks_Allo_vs_Auto_fig_with_subplots(num_plot_rows,bin_sizes_Ks_array, bin_sizes_Tc,
        demographiKS_allo_out_path, demographics_allo_run_list, run_list_name,
        demographics_auto_run_list, demographics_auto_out_path,
        xmax_Ks_array, xmax_Tc, ymax_Ks_array, ymax_Tc,
        suptitle, show_KS_predictions,include_annotation,plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here




if __name__ == '__main__':
    unittest.main()
