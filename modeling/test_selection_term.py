import math
import os
import unittest
import numpy as np
import scipy
import tskit
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import config
import ks_parsers
from figure_generation import curve_fitting

class MyTestSelectionTerm(unittest.TestCase):

    coffee_popt=[]+[9.21E-01, 9.41E-01, -1.30E-03, 3.13E-02]
    zea_popt=[9.97E-05,4.56E+03,4.67E-02]+[6.13E-01, 6.30E-01, 4.36E-02, 1.54E-01]
    populus_popt=[8.41E-05,3.88E+03,3.06E-02]+[6.33E-01, 4.17E-01, 7.99E-02, 1.61E-01]

    def test_selection_term_in_coffee(self):

        #TODO - add in popt + 4th term in Ks pdf

        coffee_num=20#coffee_num=17
        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"
        hist_comparison_out_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/Empirical_data_testing"
        specks_csv_file = "Allo_Coffea{0}_ML_rep0_Ks_by_GeneTree.csv".format(coffee_num)
        ksrates_csv_file="coffea.ks.tsv"

        splat=specks_csv_file.split("_")
        species_run_name=splat[0]+splat[1]
        real_full_path=os.path.join(ksrates_out_folder,ksrates_csv_file)

        #bin_size=0.001
        bin_size=0.002
        max_Ks=0.2
        color='blue'
        wgd_ks=0.015
        density = True
        p0 = self.coffee_popt
        const_guess=0.5
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        if not os.path.exists(hist_comparison_out_folder):
            os.makedirs(hist_comparison_out_folder)

        out_png1 = os.path.join(hist_comparison_out_folder, "3-term_" + species_run_name + "_out.png")
        list_of_hist_data= plot_histogram_with_fit(real_ks_results, p0,const_guess,
                              species_run_name, bin_size, color, wgd_ks,
                                           max_Ks, density, out_png1)

        out_png2 = os.path.join(hist_comparison_out_folder, "3-term_" + species_run_name + "_overlay.png")
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
        p0 = self.populus_popt
        const_guess=0.5
        out_png1 = os.path.join(hist_comparison_out_folder,"3-term_" + species_run_name + "_out.png")
        list_of_hist_data= plot_histogram_with_fit(real_ks_results, p0,const_guess,
                              species_run_name, bin_size, color, wgd_ks,
                                           max_Ks, density, out_png1)

        out_png2 = os.path.join(hist_comparison_out_folder,"3-term_" + species_run_name + "_overlay.png")
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
        p0 = self.zea_popt
        const_guess=0.5
        out_png1 = os.path.join(hist_comparison_out_folder, "3-term_" + species_run_name + "_out.png")
        list_of_hist_data= plot_histogram_with_fit(real_ks_results, p0,const_guess,
                              species_run_name, bin_size, color, wgd_ks,
                                           max_Ks, density, out_png1)

        out_png2 = os.path.join(hist_comparison_out_folder, "3-term_" + species_run_name + "_overlay.png")
        overlay_differences_in_curves(species_run_name, list_of_hist_data, wgd_ks, out_png2)


def plot_histogram_with_fit(Ks_results, p0,const_guess,
                          species_name, bin_size, color,WGD_ks, max_Ks, density, out_png):

    # MBE says: 600 - 1200 dpi for line drawings
    # and 350 dpi for color and half-tone artwork)
    fig = plt.figure(figsize=(5, 5), dpi=350)

    #[n, bins] = hist_data
    label="hist for " + os.path.basename(out_png).replace("_out.png","")
    if max_Ks:
        bins = np.arange(bin_size, max_Ks + 0.1, bin_size)
        hist_ys, bins, patches = plt.hist(Ks_results, bins=bins, facecolor=color, alpha=0.25,
                                    label=label, density=density)
        plt.xlim([0, max_Ks * (1.1)])
        hist_ys_1=[hist_ys, bins]


    bar_plot_xs = [0.5*(bins[i]+bins[i+1]) for i in range(0,len(bins) - 1)]  # to match bar-plot axes
    if len(p0) >4:
        exp_params=p0[0:3]
        exp_ys = [curve_fitting.wgd_kingman(x, *exp_params) for x in bar_plot_xs]
        plt.plot(bar_plot_xs, exp_ys, color='r', alpha=0.95, label="exp fit", linestyle="-")
        ln_params=p0[3:9]
        ln_ys = [curve_fitting.wgd_lognorm2(x, *ln_params) for x in bar_plot_xs]
        plt.plot(bar_plot_xs, ln_ys, color='g', alpha=0.95, label="wgd fit", linestyle="-")
    else:
        t2_ys = [curve_fitting.wgd_lognorm2(x, *p0) for x in bar_plot_xs]
        plt.plot(bar_plot_xs, t2_ys, color='g', alpha=0.95, label="wgd fit", linestyle="-")
    #hist_ys_2=[t2_ys, list(bins)]

    p0_with_const=p0+[const_guess]
    if len(p0) >4:
        fit_curve_ys, xs_for_wgd, popt=curve_fitting.fit_curve_to_xs_and_ys(bar_plot_xs, hist_ys,
                                             curve_fitting.wgd_kingman_and_lognorm_and_constant,
                                                                            p0=p0_with_const)
        t3_ys = [curve_fitting.wgd_kingman_and_lognorm_and_constant(x, *popt) for x in bar_plot_xs]
    else:
        fit_curve_ys, xs_for_wgd, popt=curve_fitting.fit_curve_to_xs_and_ys(bar_plot_xs, hist_ys,
                                             curve_fitting.wgd_lognorm2_with_constant, p0=p0_with_const)

        t3_ys = [curve_fitting.wgd_lognorm2_with_constant(x, *popt) for x in bar_plot_xs]

    print("popt:\t" + str(popt))
    hist_ys_3=[t3_ys, list(bins)]

    cont_ys = [popt[-1] for x in bar_plot_xs]
    plt.plot(bar_plot_xs, cont_ys, color='b', alpha=0.95, label="const fit", linestyle="-")

    plot_label_4 = "3-term fit"
    plt.plot(bar_plot_xs, t3_ys, color='k', alpha=0.95, label=plot_label_4, linestyle="-")

    plt.legend()
    plt.xlabel("Ks")
    plt.ylabel("Count in Bin")
    plt.title("Ks histogram for function fit\nincluding selection")
    plt.savefig(out_png)
    plt.clf()
    plt.close()

    return [hist_ys_1,hist_ys_3]

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