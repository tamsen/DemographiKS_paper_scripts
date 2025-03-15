import unittest

from figure_generation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots


class TestKsByNe(unittest.TestCase):

    def test_Ks_for_varying_Ne_Tdiv_1000(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Ne'
        specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Ne'

        demographics_TE_run_list = [False,
                                    'DGKS_Ne10_m03d14y2025_h10m17s36',
                                    "DGKS_Ne100_m03d14y2025_h09m25s29",
                                    'DGKS_Ne500_m03d14y2025_h17m53s50',
                                    'DGKS_Ne1000_m03d14y2025_h10m17s36']
        #                            'DGKS_Ne5000_m03d14y2025_h10m17s36']
        #                            "DGKS_Ne10000_m03d13y2025_h09m37s35"]

        #demographics_TE_run_list = [False, "DGKS_10_10_m5_RC7_m01d14y2025_h09m14s47",
        #                            "DGKS_100_100_m5_RC7_m01d14y2025_h09m14s18",
        #                            "DGKS_1000_1000_m5_BI40_RC7_m01d13y2025_h15m36s22",
        #                            "DGKS_5000_5000_m5_BI_40K_RC7_m01d14y2025_h09m13s52"]
        specks_TE_run_list = [False,
        'specks_TE10_m03d14y2025_h13m33s57','specks_TE100_m03d14y2025_h13m33s57',
                              'specks_TE500_m03d14y2025_h13m33s57',
         'specks_TE1000_m03d14y2025_h13m24s06']
        #'specks_TE5000_m03d14y2025_h13m24s06']
        #    ,'specks_TE10000_m03d14y2025_h13m24s06']
        #




        #specks_TE_run_list = [False,'specks_TE10_m01d13y2025_h13m18s28',
        #                      'specks_TE100_m01d13y2025_h13m17s56',
        #                      'specks_TE1000_m01d13y2025_h13m17s53',
        #                      'specks_TE5000_m01d13y2025_h13m18s40',False]

        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate
        Ks_per_YR = 0.833*10**-5
        #bin_sizes_Tc = [80, 80, 80, 400, 800]  # looks good
        xmax_Ks = [0.15 for f in demographics_TE_run_list]
        #xmax_Ks[5] = 0.3
        #xmax_Ks = [0.05, 0.05, 0.05, 0.1, 0.2, 0.4]
        #bin_sizes_Ks = [0.001, 0.001, 0.001, 0.004, 0.008]
        bin_sizes_Ks = [xmax_Ks_i / 50 for xmax_Ks_i in xmax_Ks]
        xmax_Tc = [1000,1000,2000,10000,20000, 80000]
        #xmax_Tc = [20000 for f in demographics_TE_run_list]
        bin_sizes_Tc =[xmax_Tc_i / 50 for xmax_Tc_i in xmax_Tc]
        ymax_Tc = [False for f in demographics_TE_run_list]
        run_list_num = "DGKS_1000_gen_by_Ne_fast_mut_rate_Fig R-Ne1."
        ymax_Ks = [200 for f in demographics_TE_run_list ]

        suptitle = "SLiM vs SpecKS, Tcoal and Ks"
        show_KS_predictions=[False,False,False]
        include_annotation=False
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ "Ne:" + str(config.ancestral_Ne)
        which_plot_panels_to_show_legend = [1,2,3,4]
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_TE_run_list, run_list_num,
                                     specks_TE_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_Ks, ymax_Tc,
                                      suptitle, show_KS_predictions,
                                     include_annotation,which_plot_panels_to_show_legend,plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
