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
        run_list_name = "Variation_in_dip_timing_v6"
        demographics_run_list = [False,
                                 'DGKS_delta0p5_figR1_stop0K.rand10_m10d28y2026_h09m21s21',
                                'DGKS_delta0p5_figR1_stop100K.rand10_m10d28y2026_h09m20s00',
                                 'DGKS_delta0p5_figR1_stop300K.rand10_m10d28y2026_h10m24s24',
                                 'DGKS_delta0p5_figR1_stop500K.rand10_m10d28y2026_h10m25s40',
                               'diploidization_overlay_100_300_500'
                               ,

        ]

        num_bins=25
        specks_TE5_run_list = [False, False, False, False, False, False, False]
        xmax_Ks = [0.04 for f in demographics_run_list]
        bin_sizes_Ks = [xmax_KS_i / num_bins for xmax_KS_i in xmax_Ks]
        xmax_Tc = [10000 for f in demographics_run_list]
        bin_sizes_Tc = [xmax_Tc_i / num_bins for xmax_Tc_i in xmax_Tc]
        ymax_KS = [3000 for f in demographics_run_list]
        ymax_KS[5] = 3000
        ymax_Tc = [2400 for f in demographics_run_list]

        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation = False
        plot_title_lamda = lambda config: ("Ks at Tnow\n"
                                           + ( "diff stops per chr" if isinstance(config.dij_stop, str)
                                            else ("tetras stops after "
                                          + str(int((config.dij_stop-config.dij_start)*10.0)) +
                                                  " gen")))
        plots_to_show_legend = [1]
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
        run_list_name = "Variation_in_delta_more_BI_v6"
        demographics_run_list = [False,
                                 'DGKS_delta0p0_figR1_row1.rand10.xml_m10d27y2026_h16m04s05',
                                 'DGKS_delta0p000001_figR1_row1.rand10.xml_m10d28y2026_h09m49s38',
                                 'DGKS_delta0p00001_figR1_row1.rand10.xml_m10d27y2026_h13m08s04',
                                 'DGKS_delta0p0001_figR1_row1.rand10_m10d28y2026_h09m47s56',
                                 'delta_overlay_0_e6_e4']

        specks_TE5_run_list = [False, False, False, False, False,False,False]
        xmax_Ks = [0.02 for f in demographics_run_list]
        bin_sizes_Ks = [xmax_KS_i/25 for xmax_KS_i in xmax_Ks]
        xmax_Tc = [10000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/25 for xmax_Tc_i in xmax_Tc]
        ymax_KS = [2000 for f in demographics_run_list]
        ymax_Tc = [2400 for f in demographics_run_list]
        ymax_KS[5] = 4000
        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation=False
        plot_title_lamda = lambda config: ("Ks at Tnow\n"+
                                           "delta rate:" + str(config.homoeologous_exchange_rate/10.0))
        plots_to_show_legend=[1]
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, demographiKS_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions, include_annotation,
                                     plots_to_show_legend,
                                     plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here


    def test_Ks_for_No_Migration(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/MultipleChr'

        run_list_name = "DGKS_no_migration_more_BI_v6"
        demographics_run_list = [False,
                                 'DGKS_0p0Mig50KY_figR1_row1.rand10_m10d21y2026_h12m39s15',
                                 'DGKS_0p0Mig50KY_figR1_row1.rand20_m10d21y2026_h12m40s24',
                                 'DGKS_0p0Mig50KY_figR1_row1.rand30_m10d21y2026_h12m40s59',
                                 'DGKS_0p0Mig50KY_figR1_row1.rand40_m10d28y2026_h14m39s51',
                                 '0p0Mig_overlay_v1']


        specks_TE5_run_list = [False, False, False, False, False, False, False]
        xmax_Ks = [0.02 for f in demographics_run_list]
        bin_sizes_Ks = [xmax_KS_i / 25 for xmax_KS_i in xmax_Ks]
        xmax_Tc = [10000 for f in demographics_run_list]
        bin_sizes_Tc = [xmax_Tc_i / 25 for xmax_Tc_i in xmax_Tc]
        ymax_KS = [2000 for f in demographics_run_list]
        ymax_Tc = [2400 for f in demographics_run_list]
        ymax_KS[5] = 8000
        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation = False
        plot_title_lamda = lambda config: "Ks at Tnow\n" + "Mig rate:" + str(config.mig_rate)
        plots_to_show_legend = [1]
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, demographiKS_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions, include_annotation,
                                     plots_to_show_legend,
                                     plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here

    def test_Ks_for_10K_gen_of_Migration(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/MultipleChr'


        run_list_name = "With_migration_replicates_more_BI_v6"
        demographics_run_list = [False,
            'DGKS_2p0Mig10000Y_figR1_row1.rand10_m10d21y2026_h12m44s15',
            'DGKS_2p0Mig10000Y_figR1_row1.rand20_m10d21y2026_h12m44s11',
            'DGKS_2p0Mig10000Y_figR1_row1.rand30_m10d21y2026_h12m44s09',
             'DGKS_2p0Mig10000Y_figR1_row1.rand40_m10d28y2026_h14m50s42',
                                 '2p0_mig_overlay_v1'
                                 ]


        specks_TE5_run_list = [False, False, False, False, False,False,False]
        xmax_Ks = [0.02 for f in demographics_run_list]
        bin_sizes_Ks = [xmax_KS_i/25 for xmax_KS_i in xmax_Ks]
        xmax_Tc = [10000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/25 for xmax_Tc_i in xmax_Tc]
        ymax_KS = [2000 for f in demographics_run_list]
        ymax_Tc = [2400 for f in demographics_run_list]
        ymax_KS[5] = 8000

        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation=False
        plot_title_lamda = lambda config: ("Ks at Tnow\n"+ "Mig rate:"
                                           + str(config.mig_rate/10.0))
        plots_to_show_legend=[1]
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
