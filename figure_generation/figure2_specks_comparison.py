import os
import unittest
from figure_generation.multi_row_Ks_aggregation_plots import make_multi_row_Ks_fig_with_subplots, SPKSvsDGKSPlotData


class MyTestCase(unittest.TestCase):

    def test_fig2_DGKS_vs_SPKS_Ks(self):


        q_scaling=[0,1000,100]#time dimension
        inv_q_scaling=[0.01]# 1/Time dimension

        output_png_path='/home/tamsen/Data/DemographiKS_output_from_mesx/fig1_Ne_Tdiv_RC_withQscaling.png'

        r1_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Ne_v2'
        r1_specks_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Ne_v2'

        r1_demographics_run_list = [
                                    'DGKS_Ne10_Fig1row1_v2_m03d16y2026_h17m56s18',
                                    'DGKS_Ne50_Fig1row1_v2_m03d27y2026_h10m28s27',
                                    'DGKS_Ne100_Fig1row1_v2_m03d16y2026_h17m56s20',
                                    'DGKS_Ne500_Fig1row1_v2_m03d18y2026_h10m11s04',
                                    'DGKS_Ne1000_Fig1row1_v5_m04d01y2026_h14m06s03']


        r1_specks_run_list = [
        'specks_TE10_m03d09y2026_h12m56s47','specks_TE50_m03d27y2026_h12m41s14',
                              'specks_TE100_m03d09y2026_h13m28s52',
                              'specks_TE500_m03d09y2026_h14m00s55',
                              'specks_TE1000_m03d09y2026_h13m30s53',
                              'specks_TE5000_m03d09y2026_h14m00s57']


        r1_xmax_Ks = [0.05 for f in r1_demographics_run_list]
        r1_ymax_Ks = [400 for f in r1_demographics_run_list ]
        r1_bins = [xmax_Ks_i /25 for xmax_Ks_i in r1_xmax_Ks]

        r1_show_KS_predictions=[False,False,False]
        r1_include_annotation=False
        r1_plot_title_lamda = lambda config: "Ks at Tnow\n"+ "Ne:" + str(q_scaling[1]*config.ancestral_Ne)
        r1_which_plot_panels_to_show_legend = [4]

        row1_data = SPKSvsDGKSPlotData(r1_demographiKS_data_path, r1_demographics_run_list,'',
                                       r1_specks_data_path, r1_specks_run_list,
                                       r1_bins, r1_xmax_Ks, r1_ymax_Ks,
                                       r1_show_KS_predictions,
                                       r1_include_annotation, r1_which_plot_panels_to_show_legend,
                                       r1_plot_title_lamda)


        r2_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Tdiv_v2'
        r2_specks_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Tdiv_v2'

        r2_demographics_run_list=[
            'DGKS_Tdiv100_Fig1row2_v7_m04d20y2026_h12m38s00',
            'DGKS_Tdiv1000_Fig1row2_v7_m04d20y2026_h12m38s00',
            'DGKS_Tdiv5000_Fig1row2_v7_m04d20y2026_h12m38s00',
            'DGKS_Tdiv10000_Fig1row2_v7_m04d20y2026_h12m38s00',
            'DGKS_Tdiv100000_Fig1row2_v7_m04d20y2026_h12m38s02']

        r2_specks_run_list=[ 'specks_Tdiv100_v4_m03d18y2026_h17m13s07',
                             'specks_Tdiv1000_v4_m03d18y2026_h17m15s39',
                             'specks_Tdiv5000_v4_m03d27y2026_h17m16s26',
                             'specks_Tdiv10000_v4_m03d18y2026_h18m59s33',
                             'specks_Tdiv100000_v4_m03d19y2026_h09m33s06']

        r2_xmax_Ks = [0.02,0.02,0.02,0.02,0.1,0.1,0.1]
        r2_bins  = [f /25 for f in r2_xmax_Ks]

        #r2_bins = [0.001 for f in r2_demographics_run_list]
        #r2_xmax_Ks = [0.02,0.02,0.02,0.04,0.4]
        r2_ymax_Ks = [10000 for f in r2_demographics_run_list]
        r2_show_Ks_predictions=[False,False,False]
        r2_include_annotation=False
        r2_plot_title_lamda = lambda config: "Ks at Tnow\n"+ "Twgd&Tdiv:" + str(q_scaling[2]*config.DIV_time_Ge)
        r2_which_plot_panels_to_show_legend = [4]

        row2_data = SPKSvsDGKSPlotData(r2_demographiKS_data_path, r2_demographics_run_list,'',
                                       r2_specks_data_path, r2_specks_run_list,
                                       r2_bins, r2_xmax_Ks, r2_ymax_Ks,
                                       r2_show_Ks_predictions,
                                       r2_include_annotation, r2_which_plot_panels_to_show_legend, r2_plot_title_lamda)


        r3_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_RC_v2'
        r3_specks_data_path  = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_RC_v2'

        r3_demographics_run_list = [
            'DGKS_RC10_Fig1row3_v6_m03d26y2026_h14m37s59',
                                 'DGKS_RC9_Fig1row3_v6_m03d26y2026_h14m37s59',
                                 'DGKS_RC8_Fig1row3_v6_m03d26y2026_h14m37s59',
                                 'DGKS_RC7_Fig1row3_v6_m03d26y2026_h14m37s59',
                                 'DGKS_RC6_Fig1row3_v6_m03d26y2026_h16m18s00']


        r3_specks_run_list = ['SpecKS_KSvsRC0_v1_m03d26y2026_h15m05s31',
                               False, False,False, False]

        r3_xmax_Ks = [0.05 for f in r3_demographics_run_list] #[0.025, 0.025, 0.025, 0.025, 0.025]  # 0.001  # max(demographiKS_ks_results)
        r3_ymax_Ks = [250 for f in r3_demographics_run_list]
        r3_bins = [xmax_Ks_i / 25 for xmax_Ks_i in r3_xmax_Ks]
        r3_which_plot_panels_to_show_legend = [0]
        r3_show_Ks_predictions = [False, False, False]
        r3_include_annotation = False
        #r3_plot_title_lamda = lambda config: "Ks at Tnow\n" + "RC rate:" + str(0.01*config.recombination_rate)
        r3_plot_title_lamda = lambda config: "Ks at Tnow\n" + "RC rate:"+\
                                             f"{inv_q_scaling[0]*config.recombination_rate:.1e}"

        row3_data = SPKSvsDGKSPlotData(r3_demographiKS_data_path, r3_demographics_run_list,'',
                                       r3_specks_data_path, r3_specks_run_list,
                                       r3_bins, r3_xmax_Ks, r3_ymax_Ks,
                                       r3_show_Ks_predictions,
                                       r3_include_annotation, r3_which_plot_panels_to_show_legend, r3_plot_title_lamda)

        make_multi_row_Ks_fig_with_subplots([row1_data,row2_data, row3_data],
                                            output_png_path)

        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
