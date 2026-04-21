import os
import unittest
import math
import sklearn
from matplotlib import pyplot as plt

from figure_generation.multi_row_Ks_aggregation_plots import make_multi_row_Ks_fig_with_subplots, SPKSvsDGKSPlotData


class MyTestCase(unittest.TestCase):

    def test_figS4_migration_as_secondary_contact(self):

        output_png_path='/home/tamsen/Data/DemographiKS_output_from_mesx/figS4_migration_as_secondary_contact_v2.png'
        qvalues_by_row=[0,1,1,1,1]
        r1_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2_50Y'
        r1_run_list_name = " "
        r1_demographics_run_list = [
            'DGKS_0p0Mig0KY_fig2_row1.v1_m04d02y2026_h14m14s39',
            'DGKS_0p01Mig50Y_figS4_row1.v1_m04d17y2026_h12m10s44',
            'DGKS_0p1Mig50Y_figS4_row1.v1_m04d17y2026_h12m11s10',
            'DGKS_1p0Mig50Y_figS4_row1.v1_m04d17y2026_h12m11s36',
            'DGKS_10p0Mig50Y_figS4_row1.v1_m04d17y2026_h12m12s10']

        r1_specks_run_list = [False, False, False, False, False, False]
        r1_xmax_Ks = [0.02 for f in r1_demographics_run_list]
        r1_bins = [xmax_KS_i / 25 for xmax_KS_i in r1_xmax_Ks]
        r1_ymax_Ks = [250 for f in r1_demographics_run_list]
        r1_show_Ks_predictions = [False, False, False]
        r1_include_annotation = False
        r1_which_plot_panels_to_show_legend = [4]
        #r1_plot_title_lamda = lambda config: "Ks at Tnow\n" + \
        #                                  "Mig duration:" + \
        #                                  str(qvalues_by_row[1]*(config.mig_stop - config.mig_start)) + " gen"

        r1_plot_title_lamda = lambda config: "Ks at Tnow\n" + \
                                          "Mig duration: " + \
                                          str(qvalues_by_row[1]*(config.mig_stop - config.mig_start)) + " gen, " +\
                                          "rate: " +  str(qvalues_by_row[1]*(config.mig_rate)*100) + "%"

        row1_data = SPKSvsDGKSPlotData(r1_demographiKS_data_path, r1_demographics_run_list, r1_run_list_name,
                                       r1_demographiKS_data_path, r1_specks_run_list,
                                       r1_bins, r1_xmax_Ks, r1_ymax_Ks,
                                       r1_show_Ks_predictions,
                                       r1_include_annotation, r1_which_plot_panels_to_show_legend, r1_plot_title_lamda)


        r2_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2_10KY'
        r2_run_list_name = " "
        r2_demographics_run_list = [
            'DGKS_0p0Mig0KY_fig2_row1.v1_m04d02y2026_h14m14s39',
            'DGKS_0p01Mig10KY_figS5_row1.v2_m04d20y2026_h17m51s37',
            'DGKS_0p1Mig10KY_figS5_row1.v2_m04d20y2026_h17m51s37',
            'DGKS_10p0Mig10KY_figS5_row1.v2_m04d20y2026_h17m51s37',
            'DGKS_1p0Mig10KY_figS5_row1.v2_m04d20y2026_h17m51s37',
            ]

        r2_plot_title_lamda = lambda config: "Ks at Tnow\n" + \
                                          "Mig duration: " + \
                                          str(qvalues_by_row[1]*(config.mig_stop - config.mig_start)) + " gen, " +\
                                          "rate: " +  str(qvalues_by_row[1]*(config.mig_rate)*100) + "%"


        row2_data = SPKSvsDGKSPlotData(r2_demographiKS_data_path, r2_demographics_run_list, r2_run_list_name,
                                       r1_demographiKS_data_path, r1_specks_run_list,
                                       r1_bins, r1_xmax_Ks, r1_ymax_Ks,
                                       r1_show_Ks_predictions,
                                       r1_include_annotation, r1_which_plot_panels_to_show_legend, r2_plot_title_lamda)


        make_multi_row_Ks_fig_with_subplots([row1_data,row2_data],
                                            output_png_path)



if __name__ == '__main__':
    unittest.main()
