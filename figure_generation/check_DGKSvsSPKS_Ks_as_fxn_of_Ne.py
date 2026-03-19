import unittest

from figure_generation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots


class TestKsByNe(unittest.TestCase):

    def test_Ks_for_varying_Ne_1000(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Ne_v2'
        specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Ne_v2'

        demographics_TE_run_list = [False,
                                    'DGKS_Ne10_Fig1row1_v2_m03d16y2026_h17m56s18',
                                    'DGKS_Ne100_Fig1row1_v2_m03d16y2026_h17m56s20',
                                    'DGKS_Ne500_Fig1row1_v2_m03d18y2026_h10m11s04',
                                    'DGKS_Ne1000_Fig1row1_v2_m03d18y2026_h10m11s10',
                                    ]

        #DGKS_Ne1000_Fig1row1_v2_m03d18y2026_h10m11s10
        #DGKS_Ne100_Fig1row1_v2_m03d18y2026_h10m11s02
        #DGKS_Ne10_Fig1row1_v2_m03d18y2026_h10m11s01
        #DGKS_Ne500_Fig1row1_v2_m03d18y2026_h10m11s04

        #                              'DGKS_Ne10_Fig1row1_v2_m03d09y2026_h10m59s14',
      #                              'DGKS_Ne100_Fig1row1_v2_m03d09y2026_h11m09s38',
      #                              'DGKS_Ne500_Fig1row1_v2_m03d09y2026_h11m13s36',
      #                              'DGKS_Ne1000_Fig1row1_v2_m03d09y2026_h11m12s10']
        #                            'DGKS_Ne5000_Fig1row1_v2_m03d10y2026_h12m59s11',
         #                           'DGKS_Ne10000_Fig1row1_v2_m03d09y2026_h11m15s37']
        #                            False]
        #                            'DGKS_Ne10000_Fig1row1_v2_m03d09y2026_h11m15s37']
        #                            'DGKS_Ne10_m03d14y2025_h10m17s36',
        #                            "DGKS_Ne100_m03d14y2025_h09m25s29",
        #                            'DGKS_Ne500_m03d14y2025_h17m53s50',
        #                            'DGKS_Ne1000_m03d14y2025_h10m17s36']
        #                            'DGKS_Ne5000_m03d14y2025_h10m17s36']
        #                            "DGKS_Ne10000_m03d13y2025_h09m37s35"]

        specks_TE_run_list = [False,
        'specks_TE10_m03d09y2026_h12m56s47','specks_TE100_m03d09y2026_h13m28s52',
                              'specks_TE500_m03d09y2026_h14m00s55',
         'specks_TE1000_m03d09y2026_h13m30s53',
                              'specks_TE5000_m03d09y2026_h14m00s57']
        #                      'specks_TE10000_m03d09y2026_h13m34s32']
        #'specks_TE10_m03d14y2025_h13m33s57','specks_TE100_m03d14y2025_h13m33s57',
        #                      'specks_TE500_m03d14y2025_h13m33s57',
        # 'specks_TE1000_m03d14y2025_h13m24s06',False]
        #'specks_TE5000_m03d14y2025_h13m24s06']
        #    ,'specks_TE10000_m03d14y2025_h13m24s06']
        #

        #xmax_Ks = [0.01 for f in demographics_TE_run_list]
        xmax_Ks = [0.15 for f in demographics_TE_run_list]

        bin_sizes_Ks = [xmax_Ks_i / 50 for xmax_Ks_i in xmax_Ks]
        xmax_Tc = [1000,1000,2000,10000,20000,40000, 80000]
        bin_sizes_Tc =[xmax_Tc_i / 50 for xmax_Tc_i in xmax_Tc]
        ymax_Tc = [False for f in demographics_TE_run_list]
        run_list_num = "DGKS_1000_gen_by_Ne_fast_mut_rate_Fig R-Ne1."
        ymax_Ks = [200 for f in demographics_TE_run_list ]
        #ymax_Ks = [300 for f in demographics_TE_run_list]

        suptitle = "DGKS vs SpecKS, Tcoal and Ks"
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
