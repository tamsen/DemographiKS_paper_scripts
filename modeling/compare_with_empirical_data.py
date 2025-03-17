import math
import os
import unittest
import numpy as np
import scipy
import tskit
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import ks_parsers
from figure_generation import curve_fitting

#actinidia.ks.tsv
#triticum.ks.tsv
#final_ks_values_TORX.fa

class MyTestCase(unittest.TestCase):

    MBE_dpi=350
    def test_check_versions(self):
        print(scipy.__version__)
        print(tskit.__version__)

    def test_coffee_histogram(self):

        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"
        hist_comparison_out_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/Empirical_data_testing"
        ksrates_csv_file="coffea.ks.tsv"

        species_run_name='Coffea arabica'
        real_full_path=os.path.join(ksrates_out_folder,ksrates_csv_file)

        #bin_size=0.001
        bin_size=0.002
        max_Ks=0.2
        color='blue'
        wgd_ks=0.015
        density = True

        shape=0.5 #sigma
        scale = 0.02#two_Ne#going from scale 1 to 2 halves the height and makes it twice as wide
        loc = 0.01 #predic loc=0 when the mean is over the scale value.
        amp= 1.0  #0.5*num_genes*two_Ne*bin_size
        p0 = [amp,shape,loc,scale]
        lognorm_fit_range=[0.001,0.2]
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        if not os.path.exists(hist_comparison_out_folder):
            os.makedirs(hist_comparison_out_folder)

        out_png1 = os.path.join(hist_comparison_out_folder,
                                "real_" + species_run_name + "_out_R-Empirical.png")
        list_of_hist_data= curve_fit_with_histogram(real_ks_results, p0, lognorm_fit_range,
                                                    species_run_name, bin_size, color, wgd_ks,
                                                    max_Ks, density, out_png1)

        out_png2 = os.path.join(hist_comparison_out_folder,
                                "real_" + species_run_name + "_overlay_R-Empirical.png")
        overlay_differences_in_curves(species_run_name, list_of_hist_data, wgd_ks, out_png2)

    def test_poplar_histogram(self):

        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"
        hist_comparison_out_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/Empirical_data_testing"
        ksrates_csv_file = "poplar.ks.tsv"
        species_run_name ='Populus trichocarpa'
        real_full_path = os.path.join(ksrates_out_folder, ksrates_csv_file)
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        bin_size = 0.002
        max_Ks = 0.5
        color = 'blue'
        wgd_ks=0.21
        density  = True

        shape=0.5 #sigma
        scale = 0.5#two_Ne#going from scale 1 to 2 halves the height and makes it twice as wide
        loc = 0.02 #predic loc=0 when the mean is over the scale value.
        amp= 200  #0.5*num_genes*two_Ne*bin_size
        p0 = [amp,shape,loc,scale]
        lognorm_fit_range=[0.12,max_Ks]
        out_png1 = os.path.join(hist_comparison_out_folder,
                                "real_" + species_run_name + "_out_R-Empirical.png")
        list_of_hist_data= curve_fit_with_histogram(real_ks_results, p0, lognorm_fit_range,
                                                    species_run_name, bin_size, color, wgd_ks,
                                                    max_Ks, density, out_png1)

        out_png2 = os.path.join(hist_comparison_out_folder,
                                "real_" + species_run_name + "_overlay_R-Empirical.png")
        overlay_differences_in_curves(species_run_name, list_of_hist_data, wgd_ks, out_png2)

    def test_maize_histogram(self):

        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"
        hist_comparison_out_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/Empirical_data_testing"
        ksrates_csv_file="mays.ks.tsv"
        species_run_name='Zea mays'
        real_full_path=os.path.join(ksrates_out_folder,ksrates_csv_file)
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        bin_size=0.002
        max_Ks=0.3
        color='blue'
        wgd_ks=0.13
        density = True

        shape=1.8 #sigma
        scale = 0.16#two_Ne#going from scale 1 to 2 halves the height and makes it twice as wide
        loc = 0.05 #predic loc=0 when the mean is over the scale value.
        amp= 300  #0.5*num_genes*two_Ne*bin_size
        lognorm_fit_range=[0.10,max_Ks]
        p0 = [amp,shape,loc,scale]

        out_png1 = os.path.join(hist_comparison_out_folder,
                                "real_" + species_run_name + "_out_R-Empirical.png")
        list_of_hist_data= curve_fit_with_histogram(real_ks_results, p0, lognorm_fit_range,
                                                    species_run_name, bin_size, color, wgd_ks,
                                                    max_Ks, density, out_png1)

        out_png2 = os.path.join(hist_comparison_out_folder,
                                "real_" + species_run_name + "_overlay_R-Empirical.png")
        overlay_differences_in_curves(species_run_name, list_of_hist_data, wgd_ks, out_png2)


def curve_fit_with_histogram(Ks_results, p0, lognorm_fit_range,
                             species_name, bin_size, color, WGD_ks, max_Ks, density, out_png):

    # MBE says: 600 - 1200 dpi for line drawings
    # and 350 dpi for color and half-tone artwork)
    #fig = plt.figure(figsize=(10, 10), dpi=MyTestCase.MBE_dpi)
    fig = plt.figure(figsize=(5, 5), dpi=MyTestCase.MBE_dpi)
    x = Ks_results
    # print(PAML_hist_out_file)
    #[n, bins] = hist_data
    label="hist for " + species_name
    if max_Ks:
        bins = np.arange(bin_size, max_Ks + 0.1, bin_size)
        hist_ys, bins, patches = plt.hist(x, bins=bins, facecolor=color, alpha=0.25,
                                    label=label, density=density)
        plt.xlim([0, max_Ks * (1.1)])
        hist_ys_1=[hist_ys, bins]

    #fit lognorm
    #bins_centers = [0.5 * (bins_1[i] + bins_1[i + 1]) for i in range(0, len(bins_1) - 1)]
    bins_centers = [0.5 * (bins[i] + bins[i + 1]) for i in range(0, len(bins) - 1)]
    fit_wgd_xs=[]
    fit_wgd_ys=[]
    fit_exp_xs=[]
    fit_exp_ys=[]
    for i in range(0,len(bins_centers)):
        x_value=bins_centers[i]
        if x_value >= lognorm_fit_range[0]:
            if x_value < lognorm_fit_range[1]:
                fit_wgd_xs.append(x_value)
                fit_wgd_ys.append(hist_ys[i])
        else:
            fit_exp_xs.append(x_value)
            fit_exp_ys.append(hist_ys[i])

    fit_curve_ys_ln2, xs_for_wgd, popt_wgd = \
        curve_fitting.fit_curve_to_xs_and_ys(fit_wgd_xs,fit_wgd_ys, curve_fitting.wgd_lognorm2, p0=p0)

    if fit_curve_ys_ln2:
        popt_in_sci_notation = ["{:.2E}".format(p) for p in popt_wgd]
        plot_label = "fit popt:" + ",".join(popt_in_sci_notation)
        #plt.plot(xs_for_wgd, fit_curve_ys_ln2, color='g', alpha=0.95, label=plot_label,linestyle=":")

        full_exp_ys = [curve_fitting.wgd_lognorm2(x, *popt_wgd) for x in bins]
        portion_wgd_genes = bin_size * sum(full_exp_ys)
        plot_label_3 = "wgd fit genes: " + "{:.2E}".format(portion_wgd_genes)
        plt.plot(bins, full_exp_ys, color='g', alpha=0.95, label=plot_label_3, linestyle="-")

    #fit exp
    #wgd_kingman(x, bin_size, num_genes, two_Ne)
    p0=[bin_size,len(Ks_results),100]
    fit_exp, xs_for_exp, popt_exp = \
        curve_fitting.fit_curve_to_xs_and_ys(fit_exp_xs,fit_exp_ys,
                                             curve_fitting.wgd_kingman, p0=p0)

    if fit_exp:
        popt_in_sci_notation = ["{:.2E}".format(p) for p in popt_exp]
        plot_label = "fit popt:" + ",".join(popt_in_sci_notation)
        full_exp_ys = [curve_fitting.wgd_kingman(x, *popt_exp) for x in bins]
        #plt.plot(fit_exp_xs, fit_exp, color='y', alpha=0.95, label=plot_label,linestyle=":")
        portion_exp_genes=bin_size*sum(full_exp_ys)
        plot_label_2 = "exp fit genes: " + "{:.2E}".format(portion_exp_genes)
        plt.plot(bins, full_exp_ys, color='r', alpha=0.95, label=plot_label_2,linestyle="-")

    #fit exp and wgd together:
    fit_both=False
    if fit_exp:
        p2=[*popt_exp,*popt_wgd]
        #print("popt_exp:" + str(popt_exp))
        fit_both, xs_for_both, popt_both = \
        curve_fitting.fit_curve_to_xs_and_ys(bins_centers,hist_ys,
                                             curve_fitting.wgd_kingman_and_ln, p0=p2)

        if fit_both:
            popt_in_sci_notation = ["{:.2E}".format(p) for p in popt_both]
            plot_label = "fit popt:" + ",".join(popt_in_sci_notation)
            print(plot_label)
            plt.plot(xs_for_both, fit_both, color='k', alpha=0.95,
                     label="all fit genes",linestyle="-")
            hist_data2=[fit_both, bins]

    if not fit_both:
        popt_in_sci_notation = ["{:.2E}".format(p) for p in popt_wgd]
        plot_label = "fit popt:" + ",".join(popt_in_sci_notation)
        print(plot_label)
        fit_curve_ys = [curve_fitting.wgd_lognorm2(x, *popt_wgd) for x in bins_centers]
        plt.plot(bins_centers, fit_curve_ys,
                 color='k', alpha=0.95,
                 label="all fit genes", linestyle="-")
        hist_data2 = [fit_curve_ys, bins]


    if not fit_both:
        percent_homeologous_exchange = 100*(1 -portion_wgd_genes )
    else:
        percent_homeologous_exchange=portion_exp_genes*100

    phex_str =  round(percent_homeologous_exchange,2)

    #plt.axvline(x=WGD_ks, color='b', linestyle='-', label="WGD paralog start")
    num_pairs=sum(hist_ys)
    num_after_wgd=sum([hist_ys[i] for i in range(0,len(hist_ys)) if bins[i] > WGD_ks])
    plt.legend()
    plt.xlabel("Ks")
    plt.ylabel("Count in Bin")
    plt.title("Ks histogram for {0}.\n{1}% estimated homeologous exchange.".format(
        species_name,phex_str))
    plt.savefig(out_png)
    plt.clf()
    plt.close()

    return [hist_ys_1,hist_data2]

def overlay_differences_in_curves(species_for_plot_title, list_of_hist_data, WGD_ks,out_png):

    colors = ['green','blue','brown']
    labels = ['truth','model fit']


    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    for i in range(0,len(list_of_hist_data)):
        hist_data =list_of_hist_data[i]
        [n, bins]=hist_data
        width=(bins[2]-bins[1])
        bar_plot_xs=[b + i*width for b in bins[0:len(bins)-1]] #to match bar-plot axes
        sub_bins=bins[0:len(bins) - 1]
        label=labels[i]
        #if label=='truth':

            #x_value_of_ymax = get_Ks_at_max(WGD_ks, bar_plot_xs, n, sub_bins)
            #label = label +" (peak at Ks=" +str(WGD_ks) + ")"
            #plt.axvline(WGD_ks, color=config.color_blind_friendly_color_cycle_analogs['green'],
            #            linestyle=':')
        plt.plot(bar_plot_xs,n,c=colors[i], alpha=1, label=label)
        #plt.plot(sub_bins, n, c=colors[i], alpha=1, label=label)

    diffs = [(list_of_hist_data[0][0][j]-list_of_hist_data[1][0][j]) for j in range(0,len(list_of_hist_data[1][0]))]
    rmse=math.sqrt(sum([d*d for d in diffs])/len(diffs))
    hist_data0 = list_of_hist_data[1]
    ys0=hist_data0[0]
    xs0=hist_data0[1]

    for i in range(0,len(ys0)):
        if i==0:
            ax.add_patch(Rectangle((xs0[i], ys0[i]), width, diffs[i], color=colors[2],
                                alpha=0.15, label="rmse: "+ str(round(rmse,4))))
        else:
            ax.add_patch(Rectangle((xs0[i], ys0[i]), width, diffs[i], color=colors[2],
                                alpha=0.15))

    plt.legend()
    plt.xlabel("Ks")
    ax.set(ylabel="Density")
    plt.title(species_for_plot_title, style='italic')
    plt.tight_layout()
    plt.savefig(out_png)
    plt.clf()
    plt.close()

    return n, bins


if __name__ == '__main__':
    unittest.main()