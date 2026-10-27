import math
import os
import unittest

from matplotlib import pyplot as plt

from figure_generation.colors_for_figures import lighten_color
from figure_generation.ks_with_Tc_plot_aggregations import make_Tc_Ks_fig_with_subplots
import sklearn.metrics

class TestKsMultipleChr(unittest.TestCase):


    def test_Ks_for_variation_in_dip_timing(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/MultipleChr'
        run_list_name = "Variation_in_dip_timing"
        demographics_run_list = [False,
                                 'DGKS_delta0p5_figR1_stop50K.rand10_m10d22y2026_h17m11s27',
                                 'DGKS_delta0p5_figR1_stop10K.rand10_m10d22y2026_h17m08s31',
                                 'DGKS_delta0_figR1_row1.rand10.xml_m10d21y2026_h12m49s05',
                                 'delta_Overlay'
                                 ]

        # 'DGKS_0p01Mig50Y_figR1_row1.v1_m10d19y2026_h14m54s35',
        # 'DGKS_0p01Mig50Y_figR1_row1.v3_m10d19y2026_h14m38s01',
        num_bins=25
        specks_TE5_run_list = [False, False, False, False, False, False, False]
        xmax_Ks = [0.02 for f in demographics_run_list]
        bin_sizes_Ks = [xmax_KS_i / num_bins for xmax_KS_i in xmax_Ks]
        xmax_Tc = [10000 for f in demographics_run_list]
        bin_sizes_Tc = [xmax_Tc_i / num_bins for xmax_Tc_i in xmax_Tc]
        ymax_KS = [2000 for f in demographics_run_list]
        ymax_KS[4] = 6000
        ymax_Tc = [2400 for f in demographics_run_list]

        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation = False
        plot_title_lamda = lambda config: "Ks at Tnow\n" + "delta rate:" + str(config.homoeologous_exchange_rate)
        plots_to_show_legend = [1, 2, 3, 4, 5]
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, demographiKS_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions, include_annotation,
                                     plots_to_show_legend,
                                     plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here


    def test_Ks_for_variation_in_delta(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/MultipleChr'
        run_list_name = "Variation_in_delta_more_BI_v2"
        demographics_run_list = [False,
                                 'DGKS_delta0p01_figR1_row1.rand10.xml_m10d21y2026_h12m48s59',
                                 'DGKS_delta0p01_figR1_row1.rand10.xml_m10d21y2026_h12m48s56',
                                 'DGKS_delta0p00001_figR1_row1.rand10.xml_m10d22y2026_h14m33s49',
                                 'DGKS_delta0p0000001_figR1_row1.rand10.xml_m10d22y2026_h14m33s49',
                                 'DGKS_delta0_figR1_row1.rand10.xml_m10d21y2026_h12m49s05',
                                 ]

        #'DGKS_0p01Mig50Y_figR1_row1.v1_m10d19y2026_h14m54s35',
        # 'DGKS_0p01Mig50Y_figR1_row1.v3_m10d19y2026_h14m38s01',
        specks_TE5_run_list = [False, False, False, False, False,False,False]
        xmax_Ks = [0.02 for f in demographics_run_list]
        bin_sizes_Ks = [xmax_KS_i/25 for xmax_KS_i in xmax_Ks]
        xmax_Tc = [10000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/25 for xmax_Tc_i in xmax_Tc]
        ymax_KS = [2000 for f in demographics_run_list]
        ymax_Tc = [2400 for f in demographics_run_list]

        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation=False
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ "delta rate:" + str(config.homoeologous_exchange_rate)
        plots_to_show_legend=[1,2,3,4,5]
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, demographiKS_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions, include_annotation,
                                     plots_to_show_legend,
                                     plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here



    def test_Ks_for_50_gen_of_Migration(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/MultipleChr'

        run_list_name = "DGKS_no_migration_more_BI"
        demographics_run_list = [False,
                                 'DGKS_0p0Mig50KY_figR1_row1.rand10_m10d21y2026_h12m39s15',
                                 'DGKS_0p0Mig50KY_figR1_row1.rand20_m10d21y2026_h12m40s24',
                                 'DGKS_0p0Mig50KY_figR1_row1.rand30_m10d21y2026_h12m40s59'
                                 ]

        run_list_name = "With_migration_replicates_more_BI"
        demographics_run_list = [False,
            'DGKS_2p0Mig10000Y_figR1_row1.rand10_m10d21y2026_h12m44s15',
            'DGKS_2p0Mig10000Y_figR1_row1.rand20_m10d21y2026_h12m44s11',
            'DGKS_2p0Mig10000Y_figR1_row1.rand30_m10d21y2026_h12m44s09'
                                 ]

        run_list_name = "Variation_in_delta_more_BI_v2"
        demographics_run_list = [False,
            'DGKS_delta0_figR1_row1.rand10.xml_m10d21y2026_h12m49s05',
            'DGKS_delta0p01_figR1_row1.rand10.xml_m10d21y2026_h12m48s56',
            'DGKS_delta0p01_figR1_row1.rand10.xml_m10d21y2026_h12m48s59',
                                 'DGKS_delta0p00001_figR1_row1.rand10.xml_m10d22y2026_h14m33s49',
                                 'DGKS_delta0p0000001_figR1_row1.rand10.xml_m10d22y2026_h14m33s49'
                                 ]

        #'DGKS_0p01Mig50Y_figR1_row1.v1_m10d19y2026_h14m54s35',
        # 'DGKS_0p01Mig50Y_figR1_row1.v3_m10d19y2026_h14m38s01',
        specks_TE5_run_list = [False, False, False, False, False,False,False]
        xmax_Ks = [0.02 for f in demographics_run_list]
        bin_sizes_Ks = [xmax_KS_i/25 for xmax_KS_i in xmax_Ks]
        xmax_Tc = [10000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/25 for xmax_Tc_i in xmax_Tc]
        ymax_KS = [2000 for f in demographics_run_list]
        ymax_Tc = [2400 for f in demographics_run_list]

        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation=False
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ "Mig rate:" + str(config.mig_rate)
        plots_to_show_legend=[1,2,3,4,5]
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, demographiKS_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions, include_annotation,
                                     plots_to_show_legend,
                                     plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
