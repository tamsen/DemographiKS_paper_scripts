import os
import unittest

from figure_generation.multi_row_Ks_aggregation_plots import make_multi_row_Ks_fig_with_subplots
from figure_generation.multirow_aggregations_auto_vs_allo import (make_multirow_Allo_vs_Auto_fig_with_subplots,
                                                                  AUTOvsALLOPlotData)

class MyTestCase(unittest.TestCase):

    def test_fig4_inheritance_patterns(self):

        output_root_path='/home/tamsen/Data/DemographiKS_output_from_mesx'
        output_png_path=os.path.join(output_root_path,'fig4_inheritance_patterns.png')
        qvalues_by_row=[0,100,100,10,1]

        # Fig4 row#1. After half a million years
        r1_auto_data_path = os.path.join(output_root_path,  'Auto_v2','Auto_Only_Na')
        r1_auto_run_list_name = "Ks_for_Allo_and_Auto_varying_varying_Nb_10K_Na_after_5000_years_Fig R-HomEx-Nb_v42"
        r1_auto_run_list = [
            False,
            'Auto_100KNa_500Nb_500Twgd_v7_m04d06y2026_h10m35s11',
            'Auto_100KNa_1000Nb_500Twgd_v7_m04d06y2026_h10m35s03',
            'Auto_100KNa_5000Nb_500Twgd_v7_m04d06y2026_h10m34s57',
            'Auto_100KNa_10000Nb_500Twgd_v7_m04d06y2026_h12m21s28']

        # full, w/Ne 10K
        r1_allo_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto/Auto_vs_Nb/Na10K'
        r1_allo_run_list = [False, False, False, False, False]

        r1_xmax_Ks_array = [0.015 for f in r1_auto_run_list]
        r1_bin_sizes_Ks_array  = [xmax_KS_i/25 for xmax_KS_i in r1_xmax_Ks_array]
        r1_ymax_Ks_array = [2000 for f in r1_auto_run_list]

        r1_include_annotation = False
        r1_show_Ks_predictions = [True, False, False]
        r1_which_plot_panels_to_show_legend = []
        r1_suptitle = "Auto and Allo Ks histograms\n"
        r1_plot_title_lamda = lambda config: "Polyploid pop size:" + str(config.bottleneck_Ne)

        row1_data = AUTOvsALLOPlotData(r1_auto_data_path,r1_auto_run_list, r1_auto_run_list_name,
                                       r1_allo_data_path, r1_allo_run_list,
                                       r1_bin_sizes_Ks_array, r1_xmax_Ks_array, r1_ymax_Ks_array,
                                       r1_show_Ks_predictions,
                                       r1_include_annotation, r1_which_plot_panels_to_show_legend, r1_plot_title_lamda)


        make_multirow_Allo_vs_Auto_fig_with_subplots([row1_data, row1_data], output_png_path)

        self.assertEqual(True, True)  # add assertion here

if __name__ == '__main__':
    unittest.main()
