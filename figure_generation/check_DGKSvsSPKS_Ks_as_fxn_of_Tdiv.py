import unittest

from figure_generation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots


class TestGeneLossRates(unittest.TestCase):

    def test_Ks_for_varying_varying_Tdiv_times(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Tdiv'
        #specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'
        specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Tdiv'
        demographics_TE5_run_list=[False,
          'TE05fix__m01d06y2025_h11m17s15','TE07_fix__m01d08y2025_h15m09s22',
           'TE08_fix_m01d14y2025_h09m17s20','TE09_fix_m01d13y2025_h14m13s12']

        specks_TE5_run_list=[False,
                             "specks_TE05_m01d14y2025_h09m50s16" ,
        "specks_TE07_m01d14y2025_h09m50s16",
        "specks_TE08_m01d14y2025_h09m50s16" ,
        "specks_TE09_m01d14y2025_h09m50s16" ]


        bin_sizes_Tc = [200,200, 200, 200,200]
        bin_sizes_Ks = [0.0002, 0.0002, 0.0002, 0.0002, 0.0002]
        xmax_Ks = [0.025,0.025,0.025,0.025,0.025] #0.001  # max(demographiKS_ks_results)
        xmax_Tc = [False,False,False,False,False]
        run_list_name="Ks_for_varying_varying_Tdiv_FigR-Tdiv1"
        ymax_KS = [800,800,800,800,800]
        ymax_Tc = [False for f in demographics_TE5_run_list]
        show_KS_predictions=[False,False,False]
        suptitle = "SLiM and SpecKS Ks histograms"
        include_annotation=False
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ "Tdiv:" + str(config.DIV_time_Ge)
        which_plot_panels_to_show_legend = [1, 2, 3, 4]
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

