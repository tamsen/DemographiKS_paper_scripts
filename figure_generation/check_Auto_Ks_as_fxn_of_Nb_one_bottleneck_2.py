import unittest

from figure_generation.ks_plot_aggregations_auto_vs_allo import make_Tc_Ks_Allo_vs_Auto_fig_with_subplots


class TestKsForAuto(unittest.TestCase):



    #the new row 3
    def test_Ks_for_varying_Nb_one_bottleneck_Na_fixed_at_10K_v2(self):

        #full, w/Ne 10K
        demographiKS_allo_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/NbVaries/NaFixedAt10K_Q10'
        demographics_auto_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/NbVaries/NaFixedAt10K_Q10'

        demographics_allo_run_list = [
            False,
            'Allo_10KNa_10Nb_v4_m04d08y2026_h12m26s56',
            'Allo_10KNa_50Nb_v4_m04d08y2026_h12m26s56',
            'Allo_10KNa_100Nb_v4_m04d08y2026_h12m26s56',
            'Allo_10KNa_500Nb_v4_m04d08y2026_h12m26s58',
            'Allo_10KNa_500Nb_old_v4_m04d08y2026_h16m17s04']

        demographics_auto_run_list = [
            False,
            'Auto_10KNa_10Nb_v4_m04d08y2026_h12m29s29',
            'Auto_10KNa_50Nb_v4_m04d08y2026_h12m29s29',
            'Auto_10KNa_100Nb_v4_m04d08y2026_h12m29s29',
            'Auto_10KNa_500Nb_v4_m04d08y2026_h12m29s29',
            'Auto_10KNa_500Nb_old_v4_m04d08y2026_h16m16s59']
            #'Auto_10KNa_500Nb_v4_m04d08y2026_h12m18s44']
            #'Auto_10KNa_500Nb_v4_m04d08y2026_h12m29s29' ]

        #'Auto_10KNa_20KNb_v1_m04d07y2026_h10m47s14', 'Auto_10KNa_10KNb_v1_m04d07y2026_h10m47s14',
        #xmax_Ks = [0.8  for f in demographics_auto_run_list ]
        #xmax_Ks_array = [[0.05, 0.05,0.05,0.05,0.05 ],[0.4, 0.4,0.5,0.6,0.8 ]]
        xmax_Ks_array = [[0.02 for f in demographics_auto_run_list],
                         [0.02 for f in demographics_auto_run_list]]
        bin_sizes_Ks_array  = [[xmax_KS_i/25 for xmax_KS_i in xmax_Ks_array[0]],
                        [xmax_KS_i / 25 for xmax_KS_i in xmax_Ks_array[1]]]

        xmax_Tc = [1000 for f in demographics_auto_run_list]
        bin_sizes_Tc = [xmax_Tc_i / 25 for xmax_Tc_i in xmax_Tc]

        ymax_Ks_array = [[200 for f in demographics_auto_run_list],
                   [100 for f in demographics_auto_run_list]]

        ymax_Tc = [False for f in demographics_auto_run_list]
        run_list_name = "Ks_for_Allo_and_Auto_varying_varying_Nb_10K_Na_Fig_R-NaNb6_v2"

        show_KS_predictions = [False, False, False]
        suptitle = "Auto and Allo Ks histograms\n"

        plot_title_lamda = lambda config: "Polyploid pop size:" + str(config.bottleneck_Ne)
        include_annotation = False
        num_plot_rows = 2

        make_Tc_Ks_Allo_vs_Auto_fig_with_subplots(num_plot_rows,bin_sizes_Ks_array, bin_sizes_Tc,
        demographiKS_allo_out_path, demographics_allo_run_list, run_list_name,
        demographics_auto_run_list, demographics_auto_out_path,
        xmax_Ks_array, xmax_Tc, ymax_Ks_array, ymax_Tc,
        suptitle, show_KS_predictions,include_annotation,plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here


    #the new row#1. After half a million years
    def test_Ks_for_varying_Nb_one_bottleneck_Na_fixed_at_10K_after500000_years(self):

        #full, w/Ne 10K
        demographiKS_allo_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto/Auto_vs_Nb/Na10K'
        #demographics_auto_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/NbVaries/NaFixedAt100K'
        demographics_auto_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/Auto_Only_Na'

        demographics_allo_run_list = [
           False,False,False,False,False]
        #'Auto_100KNa_50Nb_500Twgd_v7_m04d06y2026_h10m39s56','Auto_100KNa_100Nb_500Twgd_v7_m04d06y2026_h10m39s51',
        demographics_auto_run_list = [
            False,
            'Auto_100KNa_500Nb_500Twgd_v7_m04d06y2026_h10m35s11',
            'Auto_100KNa_1000Nb_500Twgd_v7_m04d06y2026_h10m35s03',
            'Auto_100KNa_5000Nb_500Twgd_v7_m04d06y2026_h10m34s57',
            'Auto_100KNa_10000Nb_500Twgd_v7_m04d06y2026_h12m21s28'

        ]
        #xmax_Ks = [0.8  for f in demographics_auto_run_list ]
        xmax_Ks_array = [[0.015 for f in demographics_auto_run_list ],[0.4, 0.4,0.5,0.6,0.8 ]]
        bin_sizes_Ks_array  = [[xmax_KS_i/25 for xmax_KS_i in xmax_Ks_array[0]],
                        [xmax_KS_i /25 for xmax_KS_i in xmax_Ks_array[1]]]

        xmax_Tc = [80000 for f in demographics_auto_run_list]
        bin_sizes_Tc = [xmax_Tc_i / 10 for xmax_Tc_i in xmax_Tc]

        ymax_Ks_array = [[2000 for f in demographics_auto_run_list],
                   [2000 for f in demographics_auto_run_list]]

        ymax_Tc = [False for f in demographics_auto_run_list]
        run_list_name = "Ks_for_Allo_and_Auto_varying_varying_Nb_10K_Na_after_5000_years_Fig R-HomEx-Nb v2"

        show_KS_predictions = [True, False, False]
        suptitle = "Auto and Allo Ks histograms\n"

        plot_title_lamda = lambda config: "Polyploid pop size:" + str(config.bottleneck_Ne)
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
