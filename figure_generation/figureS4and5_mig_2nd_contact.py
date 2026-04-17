import os
import unittest
import math
import sklearn
from matplotlib import pyplot as plt

from figure_generation.multi_row_Ks_aggregation_plots import make_multi_row_Ks_fig_with_subplots, SPKSvsDGKSPlotData


class MyTestCase(unittest.TestCase):

    def test_figS4_migration_as_secondary_contact(self):

        output_png_path='/home/tamsen/Data/DemographiKS_output_from_mesx/figS4_migration_as_secondary_contact.png'
        qvalues_by_row=[0,100,100,10,1]
        r1_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2_50Y'
        r1_run_list_name = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_0p01percent_v2."
        r1_demographics_run_list = [
            'DGKS_0p0Mig0KY_fig2_row1.v1_m04d02y2026_h14m14s39',
            'DGKS_0p01Mig50Y_figS4_row1.v1_m04d16y2026_h17m30s04',
            'DGKS_0p1Mig50Y_figS4_row1.v1_m04d16y2026_h17m30s03',
            'DGKS_10p0Mig50Y_figS4_row1.v1_m04d16y2026_h17m27s05',
            'DGKS_1p0Mig50Y_figS4_row1.v1_m04d16y2026_h17m30s03',
                                 ]

        r1_specks_run_list = [False, False, False, False, False, False]
        r1_xmax_Ks = [0.02 for f in r1_demographics_run_list]
        r1_bins = [xmax_KS_i / 25 for xmax_KS_i in r1_xmax_Ks]
        r1_ymax_Ks = [250 for f in r1_demographics_run_list]
        r1_show_Ks_predictions = [False, False, False]
        r1_include_annotation = False
        r1_which_plot_panels_to_show_legend = []
        r1_plot_title_lamda = lambda config: "Ks at Tnow\n" + \
                                          "Mig duration:" + \
                                          str(qvalues_by_row[1]*(config.mig_stop - config.mig_start)) + " gen"

        row1_data = SPKSvsDGKSPlotData(r1_demographiKS_data_path, r1_demographics_run_list, r1_run_list_name,
                                       r1_demographiKS_data_path, r1_specks_run_list,
                                       r1_bins, r1_xmax_Ks, r1_ymax_Ks,
                                       r1_show_Ks_predictions,
                                       r1_include_annotation, r1_which_plot_panels_to_show_legend, r1_plot_title_lamda)

        make_multi_row_Ks_fig_with_subplots([row1_data,row1_data],
                                            output_png_path)



if __name__ == '__main__':
    unittest.main()
