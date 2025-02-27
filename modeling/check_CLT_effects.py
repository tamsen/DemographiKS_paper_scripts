import math
import os
import unittest
import random
from scipy.optimize import curve_fit
from scipy.stats import lognorm, norm, chisquare, pearsonr
import numpy as np
from matplotlib import pyplot as plt
from mkl import second
from scipy.stats import expon

from figure_generation import curve_fitting
from figure_generation.colors_for_figures import lighten_color
from figure_generation.curve_fitting import fit_curve_to_xs_and_ys


class TestModelEffectsOfCLT(unittest.TestCase):

# "A function that can be generalized to both a Gaussian (normal distribution) "
# "and an exponential distribution is the ")generalized normal distribution"
# (also known as the "exponential power distribution"); by adjusting a single parameter
# within this function, you can smoothly transition between the characteristics
# of a Gaussian and an exponential distribution depending on the parameter value"
# https://en.wikipedia.org/wiki/Generalized_normal_distribution
# https://en.wikipedia.org/wiki/Laplace_transform
    def test_effects_of_CLT(self):

        output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/CLT_testing"

        #make a bunch of exponential distributions
        num_genes=50000
        N=1000
        xmax=10000
        bin_size = xmax/50
        bins = np.arange(0, xmax, bin_size)

        list_of_expon_distributions=[]
        for num_combined_distributions in range(0,100):
            random.seed(num_combined_distributions*10+10)
            data=list(theoretical_coalescent(num_genes, N))
            list_of_expon_distributions.append(data)


        #build a new Tc drawing from multiple different distributions
        resampled_distributions_ys=[]
        fit_distributions_ys=[]
        data_labels=[]


        two_Ne = 2.0 * N
        print("Tc plot bin_size_in_time: " + str(bin_size))
        kingman = [min(num_genes,
                       (bin_size * num_genes / two_Ne) * math.e ** ((-1 * i) / two_Ne))
                   for i in bins[0:-1]]
        resampled_distributions_ys.append(kingman)
        data_labels.append("Perfect Kingman distribution")
        mu=1.0/two_Ne


        bin_size = bins[1] - bins[0]
        #wgd_normal(x, amp, mu, sig):
        popt = [num_genes * bin_size, two_Ne- (bin_size/2),bin_size]
        gaussian_ys=[curve_fitting.wgd_normal(x, *popt) for x in bins[0:-1]]
        resampled_distributions_ys.append(gaussian_ys)
        data_labels.append("Perfect Gaussian distribution")

        new_data = [d for d in list_of_expon_distributions[0]]
        png_out_i = os.path.join(output_folder, "CLT_0RC.png")
        dgx_hist_ys, bins = plot_mrca(new_data, bins, png_out_i)
        resampled_distributions_ys.append(dgx_hist_ys)
        data_labels.append("0 RC events simulated")

        data_colors = ['darkblue','darkred', 'black']

        list_of_fraction_RC_events_to_test = [0.05,0.1,0.2,0.5,0.75]
        blue_lightening=[0.75,0.6,0.4,0.3]
        for fraction in list_of_fraction_RC_events_to_test:

            new_data = [d for d in list_of_expon_distributions[0]]
            for k in range(0, int(num_genes*fraction)):
                Tc1=list_of_expon_distributions[0][k]
                Tc2=list_of_expon_distributions[1][k]
                new_data[k] = 0.5 * (Tc1 + Tc2)

            png_out_i = os.path.join(output_folder, "CLT_" + str(fraction) + "%RC.png")
            dgx_hist_ys, bins = plot_mrca(new_data, bins, png_out_i)
            resampled_distributions_ys.append(dgx_hist_ys)
            data_labels.append(str(fraction) + " %RC events simulated")

        new_colors=[lighten_color('blue', amount=a) for a in blue_lightening]
        data_colors= data_colors+ new_colors

        list_of_num_RC_events_to_test=[1,2,5,10,20,50,80]
        red_lightening=[0.2,0.5,0.6,0.65,0.7,0.75]
        for num_RC_events in list_of_num_RC_events_to_test:
            num_combined_distributions=num_RC_events+1
            new_data=[0 for k in range(0, num_genes)]
            for k in range(0, num_genes):
                 Tcs_for_Kth_genes=[list_of_expon_distributions[j][k]
                                    for j in range(0, num_combined_distributions)]
                 avg_Tc_for_Kth_gene= sum(Tcs_for_Kth_genes)/num_combined_distributions
                 new_data[k]=avg_Tc_for_Kth_gene#
            png_out_i = os.path.join(output_folder, "CLT_" + str(num_RC_events) +"RC.png")
            dgx_hist_ys, bins_xs = plot_mrca(new_data,bins,png_out_i)
            resampled_distributions_ys.append(dgx_hist_ys)

            data_labels.append(str(num_RC_events) + " RC events simulated")

        data_colors.append("gray")
        new_colors=[lighten_color('red', amount=a) for a in red_lightening]
        data_colors= data_colors+ new_colors

        png_out = os.path.join(output_folder, "composite_CLT_RC.png")
        plot_composite_tc(resampled_distributions_ys,
                          bins,data_labels,data_colors,png_out)

        csv_out = os.path.join(output_folder, "composite_CLT_RC.csv")
        write_distribution_tests(csv_out, data_labels,bins, resampled_distributions_ys)

        self.assertEqual(True, True)  # add assertion here


    def test_fits_CLT(self):

        output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/CLT_testing"
        csv_in = os.path.join(output_folder, "composite_CLT_RC.csv")
        simulated_results, bins =read_distribution_tests(csv_in)
        #fit_results = {}
        bins1=bins[0:-1]



        num_genes=50
        bin_size = bins[1] - bins[0]
        two_Ne= 2000
        shape=1  #sigma
        scale = two_Ne#going from scale 1 to 2 halves the height and makes it twice as wide
        loc = 0 #predic loc=0 when the mean is over the scale value.
        amp=0.5*num_genes*two_Ne*bin_size
        #https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lognorm.html
        #https://en.wikipedia.org/wiki/Skewness
        #https://en.wikipedia.org/wiki/Log-normal_distribution
        #https://variation.com/wp-content/distribution_analyzer_help/hs128.htm#:~:text=As%20the%20skewness%20goes%20to,multiplied%20together%20as%20shown%20below.
        #https://www.itl.nist.gov/div898/handbook/eda/section3/eda3669.htm
        # mean - median / std dev (sigma)
        # loc = expected value? s = skew? scale = std dev?
        # https://stackoverflow.com/questions/8870982/how-do-i-get-a-lognormal-distribution-in-python-with-mu-and-sigma
        
        #https://stackoverflow.com/questions/8870982/how-do-i-get-a-lognormal-distribution-in-python-with-mu-and-sigma
        # here is a solution for getting the scipy.stats.lognorm distribution if the mean mu and standard deviation
        # sigma of the lognormal distribution are known
        #a = 1 + (sigma / mu) ** 2
        #s = np.sqrt(np.log(a))
        #scale = mu / np.sqrt(a)
        # https://www.quora.com/Why-does-a-lognormal-distribution-approach-a-normal-distribution-with-a-small-coefficient-of-variation
        # https://statproofbook.github.io/P/lognorm-pdf.html#:~:text=Proof:%20A%20log%2Dnormally%20distributed%20random%20variable%20is,function%20of%20a%20normal%20random%20variable:%20Y%E2%88%BCN(%CE%BC%2C%CF%832)%E2%87%92X=exp(Y)%E2%88%BClnN(%CE%BC%2C%CF%832).
        p2 = [amp,shape,loc,scale]

        png_out = os.path.join(output_folder, "fit_composite_CLT_RC.png")
        fig = plt.figure(figsize=(10, 10), dpi=350)
        test_cases=list(simulated_results.keys())
        #'Perfect Gaussian distribution' <- will not fit'0 RC events simulated',
        #'Perfect Kingman distribution' <- will not fit
        test_cases=[
            '0.05 %RC events simulated',
            '0.1 %RC events simulated', '0.5 %RC events simulated', '0.75 %RC events simulated',
            '1 RC events simulated','2 RC events simulated', '5 RC events simulated',
            '10 RC events simulated', '20 RC events simulated', '50 RC events simulated']
        for test_case in test_cases:
            test_data=simulated_results[test_case]
            fit_curve_ys_ln2, xs_for_wgd, popt = \
                curve_fitting.fit_curve_to_xs_and_ys(bins1,test_data, curve_fitting.wgd_lognorm2, p0=p2)

            popt_in_sci_notation=["{:.2E}".format(p) for p in popt]
            plot_label= test_case + ", fit ln2, popt:"  + ",".join(popt_in_sci_notation)
            plt.plot(xs_for_wgd, fit_curve_ys_ln2, color='g', alpha=0.95, label=plot_label)
            plt.plot(bins1,test_data, color='g', alpha=0.95, label="true data", linestyle=":")

        plt.title('label')
        plt.legend()
        plt.xlabel("MRCA time")
        plt.ylabel("# genes in bin")
        plt.savefig(png_out)
        plt.clf()
        plt.close()

        print(simulated_results.keys())


    def test_read_in_of_CLT(self):

        output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/CLT_testing"
        csv_in = os.path.join(output_folder, "composite_CLT_RC.csv")
        png_out = os.path.join(output_folder, "fit_CLT_RC.png")
        simulated_results, bins =read_distribution_tests(csv_in)
        fit_results = {}
        bins1=bins[0:-1]
        fig = plt.figure(figsize=(10, 10), dpi=350)
        label = "Coalescent Times For Genes From Sampled From Two Ancestral Genomes"

        test_case_name='1 RC events simulated'
        test_data=simulated_results[test_case_name]

        num_genes=len(test_data)
        bin_size = bins[1] - bins[0]
        two_Ne= 2000
        p = [num_genes * bin_size, two_Ne - (bin_size / 2), bin_size]
        fit_curve_ys_norm, xs_for_wgd, popt = \
                curve_fitting.fit_curve_to_xs_and_ys(bins1,test_data, curve_fitting.wgd_normal, p0=p)
        plt.plot(xs_for_wgd, fit_curve_ys_norm, color='r', alpha=0.95, label="fit norm")



        #popt2=[num_genes * bin_size,bin_size,two_Ne - (bin_size / 2),0.1]
        popt2 = [num_genes * bin_size, bin_size,two_Ne, 1]
        #amp, scale, x_shift, skew
        #return amp * lognorm.pdf(scale * x + x_shift, skew)
        #pdf(x, s, loc=0, scale=1)

        #veery pointy
        #shape=1
        #loc=two_Ne
        #scale=num_genes

        #less pointy
        shape=1
        loc=0.5
        scale=0.999 #halfve the h

        scale = 2 #going from scale 1 to 2 halves the height and makes it twice as wide
        scale = two_Ne
        loc = 0 #predic loc=0 when the mean is over the scale value.
        popt2 = [shape,loc,scale]
        test_bins = np.arange(0,5, 0.1)
        #fit_curve_ys_ln = [curve_fitting.wgd_lognorm(x, *popt2) for x in test_bins]
        fit_curve_ys_ln = [0.5*num_genes*two_Ne*bin_size*lognorm.pdf(x, *popt2) for x in bins1]
        plt.plot(bins1, fit_curve_ys_ln, color='c', alpha=0.95, label="fit log_loc0")

        amp=0.5*num_genes*two_Ne*bin_size
        p2 = [amp,shape,loc,scale]
        fit_curve_ys_ln2, xs_for_wgd, popt = \
                curve_fitting.fit_curve_to_xs_and_ys(bins1,test_data, curve_fitting.wgd_lognorm2, p0=p2)

        popt_in_sci_notation=["{:.2E}".format(p) for p in popt]

        plot_label= test_case_name + ", fit ln2, popt:"  + ",".join(popt_in_sci_notation)
        fit_results[plot_label]=fit_curve_ys_ln2
        plt.plot(xs_for_wgd, fit_curve_ys_ln2, color='g', alpha=0.95, label=plot_label)


        print(str(popt))

        plt.plot(bins1, test_data, color='k', alpha=0.95, label=label, linestyle=":")
        #plt.plot( xs_for_wgd,fit_curve_ys, color='gray', alpha=0.95, label="fit")

        plt.title(label)
        # plt.xlim([0, max_mrca])
        # plt.ylim([0, 300])
        plt.legend()
        plt.xlabel("MRCA time")
        plt.ylabel("# genes in bin")
        plt.title(label)
        plt.savefig(png_out)
        plt.clf()
        plt.close()

        png_out = os.path.join(output_folder, "fit_composite_CLT_RC.png")
        plot_composite_fits(simulated_results, fit_results, bins,
                            [], [], png_out)

        print(simulated_results.keys())

def write_distribution_tests(csv_out, data_labels, bins, resampled_distributions_ys):
    with open(csv_out, 'w') as f:

        bins_str = [str(b) for b in bins]
        f.write("Bins\t"+ "\t".join(bins_str) +"\n")

        f.write("RC\tDistributionData\n")
        for i in range(0, len(resampled_distributions_ys)):
            data_as_str = [str(d) for d in resampled_distributions_ys[i]]
            f.write(data_labels[i] + "\t" + "\t".join(data_as_str) + "\n")

def read_distribution_tests(csv_in):
    results={}
    with open(csv_in, 'r') as f:
        lines=f.readlines()
        bins_splat=lines[0].split('\t')
        bins_as_floats = [float(b) for b in bins_splat[1:-1]]
        for l in lines[2:-1]:
            splat=l.split('\t')
            label=splat[0]
            data=splat[1:-1]
            data_as_float=[float(d) for d in data]
            results[label]=data_as_float
    return results,bins_as_floats



def plot_composite_fits(simulated_results_dict, fit_results_dict, bins,
                        data_labels, data_colors, png_out):

    fig, ax = plt.subplots(1, 2, figsize=(40, 20))

    label = "Coalescent Times For Genes From Sampled From Two Ancestral Genomes"

    sim_keys=list(simulated_results_dict.keys())
    for i in range(0, len(sim_keys)):
        key=sim_keys[i]
        ax[0].plot(bins[0:-1], simulated_results_dict[key],label=key)
        #, color=data_colors[i], alpha=0.95, label=data_labels[i])

    sim_keys=list(fit_results_dict.keys())
    for i in range(0, len(sim_keys)):
        key=sim_keys[i]
        ax[1].plot(bins[0:-1], fit_results_dict[key],label=key)
        #, color=data_colors[i], alpha=0.95, label=data_labels[i])


    plt.title(label)

    ax[0].legend()
    ax[1].legend()
    ax[0].set(xlabel="MRCA time")
    ax[0].set(ylabel="# genes in bin")
    ax[0].set(title=label)
    plt.savefig(png_out)
    plt.clf()
    plt.close()

def plot_composite_tc(data_list_A, bins, data_labels, data_colors, png_out):

    fig = plt.figure(figsize=(10, 10), dpi=350)
    #fig, ax = plt.subplots(1, 2, figsize=(40, 20))
    #fig.suptitle(suptitle)
    # Co.T=(1/2N)*e^-((t-1)/2N))
    label = "Coalescent Times For Genes From Sampled From Two Ancestral Genomes"
    num_genes=len(data_list_A[0])

    for i in range(0, len(data_list_A)):
        plt.plot(bins[0:-1], data_list_A[i], color=data_colors[i], alpha=0.95, label=data_labels[i])

    #for i in range(0,len(data_list_B)):
    #    data=data_list_B[i]
    #    if data:
    #        ax[1].plot(bins[0:-1],data,color=data_colors[i], alpha=0.95,label=data_labels[i])

    #ax[1].plot([1,2,3],[4,5,6], color=data_colors[i], alpha=0.95, label=data_labels[i])
    #plt.xlim([0, max_mrca * (1.1)])
    plt.title(label)
    # plt.xlim([0, max_mrca])
    # plt.ylim([0, 300])
    plt.legend()
    #ax[0].set(xlim=[0, xmax])
    #ax[0].set(xlabel="MRCA time")
    #ax[0].set(ylabel="# genes in bin")
    #ax[0].set(title=label)
    plt.xlabel("MRCA time")
    plt.ylabel("# genes in bin")
    plt.title(label)
    plt.savefig(png_out)
    plt.clf()
    plt.close()

def plot_mrca(theoretical_mrcas_by_gene, bins, png_out):

    fig = plt.figure(figsize=(10, 10), dpi=350)
    label = "Coalescent Times For Genes From Sampled Two Ancestral Genomes"
    num_genes=len(theoretical_mrcas_by_gene)
    dgx_hist_ys, bins, patches = plt.hist(theoretical_mrcas_by_gene, bins=bins, facecolor='c', alpha=0.25,
                 label='Theoretical Tcoal by gene (total: ' + str(num_genes) + ')',
                 density=False)

    #plt.xlim([0, max_mrca * (1.1)])
    plt.title(label)
    # plt.xlim([0, max_mrca])
    # plt.ylim([0, 300])
    plt.legend()
    plt.xlabel("MRCA time")
    plt.ylabel("# genes in bin")
    plt.savefig(png_out)
    plt.clf()
    plt.close()
    return dgx_hist_ys, bins

def theoretical_coalescent(num_genes,N):
    #Co.T=(1/2N)*e^-((t-1)/2N))
    loc=0
    xscale=2.0*N
    random_draws_from_distribution = expon.rvs(loc=loc, size=num_genes, scale=xscale)
    return random_draws_from_distribution

if __name__ == '__main__':
    unittest.main()
