import unittest

from figure_generation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots


class TestGeneLossRates(unittest.TestCase):

    def test_Ks_for_varying_varying_Tdiv_times(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Tdiv_v2'
        specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Tdiv_v2'


        demographics_TE5_run_list=[False,
                                   'DGKS_Tdiv100_Fig1row2_v6_m03d17y2026_h13m57s07',
                                   'DGKS_Tdiv1000_Fig1row2_v6_m03d17y2026_h11m21s21',
                                   'DGKS_Tdiv5000_Fig1row2_v6_m03d27y2026_h10m10s50',
                                   'DGKS_Tdiv10000_Fig1row2_v6_m03d17y2026_h11m21s32',
                                   'DGKS_Tdiv100000_Fig1row2_v6_m03d19y2026_h09m37s47']

        specks_TE5_run_list=[False,
                             'specks_Tdiv100_v4_m03d18y2026_h17m13s07',
                             'specks_Tdiv1000_v4_m03d18y2026_h17m15s39',
                             'specks_Tdiv5000_v4_m03d27y2026_h17m16s26',
                             'specks_Tdiv10000_v4_m03d18y2026_h18m59s33',
                             'specks_Tdiv100000_v4_m03d19y2026_h09m33s06']




        bin_sizes_Ks = [0.001 for f in demographics_TE5_run_list]
        #bin_sizes_Ks = [0.0004,0.0002,0.001,0.002,0.01,0.08]
        #xmax_Ks = [0.015 for f in demographics_TE5_run_list]
        xmax_Ks = [0.02,0.02,0.02,0.02,0.04,0.4,0.4,0.4]
        #xmax_Ks = [0.01, 0.001, 0.004, 0.01, 0.02, 0.2, 0.4]
        #xmax_Ks = [0.0000015*10**i for i in range(0,len(demographics_TE5_run_list))]
        #xmax_Ks = [0.15 for f in demographics_TE5_run_list]
        bin_sizes_Ks = [f/25 for f in xmax_Ks]
        bin_sizes_Ks = [f /20 for f in xmax_Ks]
        #xmax_Ks = [0.025,0.025,0.025,0.025,0.025] #0.001  # max(demographiKS_ks_results)
        xmax_Tc = [500 for f in demographics_TE5_run_list]
        bin_sizes_Tc = [f/25 for f in xmax_Tc]
        run_list_name="Ks_for_varying_varying_Tdiv_FigR-Tdiv1"
        ymax_KS = [10000 for f in demographics_TE5_run_list]
        #ymax_KS = [800,800,400,400,200,200]
        ymax_Tc = [False for f in demographics_TE5_run_list]
        show_KS_predictions=[False,False,False]
        suptitle = "SLiM and SpecKS Ks histograms"
        include_annotation=False
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ "Tdiv:" + str(config.DIV_time_Ge)
        which_plot_panels_to_show_legend = []#1, 2, 3, 4,5]
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                          demographiKS_out_path, demographics_TE5_run_list, run_list_name,
                                          specks_TE5_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,suptitle,
                                     show_KS_predictions,include_annotation,
                                     which_plot_panels_to_show_legend,plot_title_lamda
                                     )

        self.assertEqual(True, True)  # add assertion here

if __name__ == '__main__':
    unittest.main()

