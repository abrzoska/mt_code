import pandas as pd

adapted = ["micOch1","criGri1","mesAur1","HLcasCan1","hetGla2","HLfukDam1","HLmarMar1","triMan1","chrAsi1","ornAna2","vicPac2","HLturTru3","HLorcOrc1","HLdelLeu1","lipVex1","phyCat1","balAcu1","HLbalMys1","bosMut1","panHod1","HLenhLut1","odoRosDiv1","lepWed1","neoSch1","conCri1"]
adapted_label = "hx"
adapted2 = ["hg38","pantro5","gorgor5","panpan2","ponabe2","hlcascan1","hetgla2","hlfukdam1","loxafr3","hlturtru3","hlorcorc1","hldelleu1","lipvex1","phycat1","balacu1","hlbalmys1","hlrhisin1","hlhiparm1","eptfus1","myodav1","myobra1","myoluc2","hlminnat1","hldesrot1"]
adapted_label2 = "ll"
adapted_dict = {
  adapted_label: adapted,
  adapted_label2: adapted2,
}
adapted_labels = list(adapted_dict.keys())
cis_reg_file = "/home/Timo/encode9_trimmed.bed"
gene_file = "/home/Timo/mart_export2.txt"
upstream_region = 50000
indel_file_name = "/moreAddSpace/tlin/all_chromosomes_0.bed"
#Todo Adjust these variables
query_species = ["nanGal1"]

folder = "/moreAddSpace/tlin"
cis_regs = pd.read_csv(cis_reg_file, delimiter='\t', names=["chr", "start", "end", "element_id", "rgb"])
df = cis_regs.groupby(by="chr")

header = ",".join(
    ["filename", "query_species_only_count_del", "query_species_only_count_in", "indel_non_adapted_included_count"])
for x in adapted_labels:
    header = header + "," + x
indel_groups = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chrX', 'chrY']
columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand,", "ThickStart", "ThickEnd", "ItemRGB"]

indel_file_name = f"{folder}/chr1_indels.bed"
indels = pd.read_csv(indel_file_name, sep="\t", skiprows=0, header=None, names=columns)
indels = indels.drop_duplicates(subset=['Name'])
indels["Batch"] = indels["Name"].str.split(".")
indels["Batch"] = indels["Batch"].apply(lambda x: x[-1] if len(x) > 0 else 0)
