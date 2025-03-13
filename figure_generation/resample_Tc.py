import os.path
import unittest

import numpy as np
import tskit
from matplotlib import pyplot as plt
import config
from modules import trees_file_processor


class TestResampleTc(unittest.TestCase):

    def test_ResampleTc(self):

        #data from simulations at "/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Inb"
        output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/trees_file_testing"
        trees_file="allotetraploid_bottleneck_trees_at_div_32.txt"
        png_number="32"
        input_xml_file="Inb32v1.used.xml"
        full_path_to_trees=os.path.join(output_folder,trees_file)
        full_path_to_config=os.path.join(output_folder,input_xml_file)
        config_used = config.DemographiKS_config(full_path_to_config)
        genome_length = config_used .total_num_bases
        gene_length = 3 * config_used.num_codons_in_a_gene
        num_genes = int(genome_length / gene_length)

        xmax=10000
        bin_size = xmax/50
        bins = np.arange(0, xmax, bin_size)

        pairs=[[1,2],[1,5],[1,100],[1,500]]
        ts = tskit.load(full_path_to_trees)
        data_list=[]
        data_labels=[]

        for pair in pairs:
            reduced_ts = ts.simplify(pair)  # simplify to the first 10 samples
            # get tree for every gene
            gene_starts = [i * gene_length for i in range(0, num_genes)]
            gene_centers = [g + int(gene_length / 2) for g in gene_starts]
            mrcas_by_gene = []
            for gene_center in gene_centers:
                tree_for_gene = reduced_ts.at(gene_center)
                if not tree_for_gene.has_multiple_roots:
                    mrca = reduced_ts.node(tree_for_gene.root)
                    mrcas_by_gene.append(mrca.time)

            label_for_pair = str(pair[0]) + "_" + str(pair[1])
            png_out = os.path.join(output_folder, label_for_pair + "_" + str(png_number) + ".png")
            dgx_hist_ys, bins_xs = plot_mrca(mrcas_by_gene, bins, png_out)
            data_list.append(dgx_hist_ys)
            data_labels.append(label_for_pair)

        png_out = os.path.join(output_folder,"mrca_hist_overlay_"
                               + str(png_number) + ".png")
        plot_composite_mrca(data_list, bins, data_labels, png_out)
        self.assertEqual(True, True)  # add assertion here


def plot_composite_mrca(data_list, bins, data_labels, png_out):

    fig = plt.figure(figsize=(10, 10), dpi=350)
    # Co.T=(1/2N)*e^-((t-1)/2N))
    label = "Coalescent Times For Genes From Sampled Two Ancestral Genomes"
    num_genes=len(data_list[0])

    for i in range(0,len(data_list)):
        data=data_list[i]
        plt.plot(bins[0:-1],data,alpha=0.95,label=data_labels[i])

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

def plot_mrca(mrcas_by_gene, bins, png_out):

    fig = plt.figure(figsize=(10, 10), dpi=350)
    # Co.T=(1/2N)*e^-((t-1)/2N))
    label = "Coalescent Times For Genes From Sampled Two Ancestral Genomes"
    num_genes=len(mrcas_by_gene)
    dgx_hist_ys, bins, patches = plt.hist(mrcas_by_gene, bins=bins, facecolor='c', alpha=0.5,
             label='Simulated Tcoal by gene (total: ' + str(num_genes) + ')',
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

if __name__ == '__main__':
    unittest.main()
