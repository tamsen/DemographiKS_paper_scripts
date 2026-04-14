import os
import unittest
from figure_generation.multi_row_Ks_aggregation_plots import make_multi_row_Ks_fig_with_subplots, MulitPlotData


class MyTestCase(unittest.TestCase):

    def test_fig3_migration_as_gradual_speciation(self):

        output_png_path='/home/tamsen/Data/DemographiKS_output_from_mesx/fig3_migration_as_gradual_speciation.png'

        r1_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2/0p01percent'
        #r3_run_list_name = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_0p01percent_v2."
        r1_demographics_run_list = [
            'DGKS_0p0Mig0KY_fig2_row1.v1_m04d02y2026_h14m14s39',
                                 'DGKS_0p01Mig50KY_fig2_row1.v1_m04d02y2026_h18m14s58',
                                 'DGKS_0p01Mig100KY_fig2_row1.v1_m04d02y2026_h18m14s58',
                                 'DGKS_0p01Mig250KY_fig2_row1.v1_m04d02y2026_h18m14s58',
                                 'DGKS_0p01Mig500KY_fig2_row1.v1_m04d02y2026_h18m14s58',
                                 ]
        r1_q_value=100
        r1_specks_run_list = [False, False, False, False, False, False, False]
        r1_xmax_Ks = [0.02 for f in r1_demographics_run_list]
        r1_bins = [xmax_KS_i / 25 for xmax_KS_i in r1_xmax_Ks]
        r1_ymax_Ks = [250 for f in r1_demographics_run_list]
        r1_show_Ks_predictions = [False, False, False]
        r1_include_annotation = False
        r1_which_plot_panels_to_show_legend = []
        r1_plot_title_lamda = lambda config: "Ks at Tnow\n" + \
                                          "Mig duration:" + \
                                          str(r1_q_value*(config.mig_stop - config.mig_start)) + " gen"

        row1_data = MulitPlotData(r1_demographiKS_data_path, r1_demographics_run_list,
                                  r1_demographiKS_data_path, r1_specks_run_list,
                                  r1_bins, r1_xmax_Ks, r1_ymax_Ks,
                                  r1_show_Ks_predictions,
                                  r1_include_annotation, r1_which_plot_panels_to_show_legend, r1_plot_title_lamda)

        r2_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2/0p10percent'
        r2_run_list_name = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_0p1percent.v2."
        r2_demographics_run_list = ['DGKS_0p0Mig0KY_fig2_row1.v1_m04d02y2026_h14m14s39',
                                     'DGKS_0p1Mig50KY_fig2_row1.v1_m04d02y2026_h14m12s46',
                                     'DGKS_0p1Mig100KY_fig2_row1.v1_m04d02y2026_h14m12s46',
                                     'DGKS_0p1Mig250KY_fig2_row1.v1_m04d02y2026_h14m12s46',
                                     'DGKS_0p1Mig500KY_fig2_row1.v1_m04d02y2026_h14m12s47']
        r2_q_value=100
        r2_plot_title_lamda = lambda config: "Ks at Tnow\n" + \
                                          "Mig duration:" + \
                                          str(r2_q_value*(config.mig_stop - config.mig_start)) + " gen"

        row2_data = MulitPlotData(r2_demographiKS_data_path, r2_demographics_run_list,
                                  r1_demographiKS_data_path, r1_specks_run_list,
                                  r1_bins, r1_xmax_Ks, r1_ymax_Ks,
                                  r1_show_Ks_predictions,
                                  r1_include_annotation, r1_which_plot_panels_to_show_legend, r2_plot_title_lamda )

        r3_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2/1p00percent'
        r3_run_list_name = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_1percent_v2."
        r3_demographics_run_list = ['DGKS_0p0Mig0KY_fig2_row1.v1_m04d02y2026_h14m14s39',
                                     'DGKS_1p0Mig50KY_fig2_row1.v1_m04d02y2026_h14m44s20',
                                     'DGKS_1p0Mig100KY_fig2_row1.v1_m04d02y2026_h14m48s10',
                                     'DGKS_1p0Mig250KY_fig2_row1.v1_m04d02y2026_h14m48s14',
                                     'DGKS_1p0Mig500KY_fig2_row1.v1_m04d02y2026_h14m48s19',
                                     ]
        r3_q_value=10
        r3_plot_title_lamda = lambda config: "Ks at Tnow\n" + \
                                          "Mig duration:" + \
                                          str(r2_q_value*(config.mig_stop - config.mig_start)) + " gen"

        row3_data = MulitPlotData(r3_demographiKS_data_path, r3_demographics_run_list,
                                  r1_demographiKS_data_path, r1_specks_run_list,
                                  r1_bins, r1_xmax_Ks, r1_ymax_Ks,
                                  r1_show_Ks_predictions,
                                  r1_include_annotation, r1_which_plot_panels_to_show_legend,r3_plot_title_lamda)


        r4_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2/10percent'
        r4_run_list_name = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_10percent_v2."
        r4_demographics_run_list = ['DGKS_0p0Mig0KY_fig2_row1.v1_m04d02y2026_h14m14s39',
                                     'DGKS_10p0Mig50KY_fig2_row1.v1_m04d02y2026_h18m08s27',
                                     'DGKS_10p0Mig100KY_fig2_row1.v1_m04d02y2026_h18m08s27',
                                     'DGKS_10p0Mig250KY_fig2_row1.v1_m04d02y2026_h18m08s27',
                                     'DGKS_10p0Mig500KY_fig2_row1.v1_m04d02y2026_h18m08s27',
                                     ]

        r4_q_value=1
        r4_plot_title_lamda = lambda config: "Ks at Tnow\n" + \
                                          "Mig duration:" + \
                                          str(r2_q_value*(config.mig_stop - config.mig_start)) + " gen"
        row4_data = MulitPlotData(r4_demographiKS_data_path, r4_demographics_run_list,
                                  r1_demographiKS_data_path, r1_specks_run_list,
                                  r1_bins, r1_xmax_Ks, r1_ymax_Ks,
                                  r1_show_Ks_predictions,
                                  r1_include_annotation, r1_which_plot_panels_to_show_legend, r4_plot_title_lamda)

        make_multi_row_Ks_fig_with_subplots([row1_data,row2_data,row3_data,row4_data],
                                            output_png_path)

if __name__ == '__main__':
    unittest.main()
