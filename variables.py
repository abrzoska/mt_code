import list_util
##SPECIES DIE ANALYSIERT WERDEN SOLL
query_species="nanGal1"
##NAHE VERWANDTE DER SPEZIES (currently unused)
in_group= ["jacJac1","micOch1","criGri1","mesAur1","perManBai1","HLmusCar1","HLmusPah1","rn6","HLmerUng1","mm10"]
##LISTE ADAPTIERTER SPEZIES (INDELS DIE DIESEN TEILEN WERDEN HERVORGEHOBEN)
##hx
adapted=["micOch1","criGri1","mesAur1","HLcasCan1","hetGla2","HLfukDam1","HLmarMar1","triMan1","chrAsi1","ornAna2","vicPac2","HLturTru3","HLorcOrc1","HLdelLeu1","lipVex1","phyCat1","balAcu1","HLbalMys1","bosMut1","panHod1","HLenhLut1","odoRosDiv1","lepWed1","neoSch1","conCri1"]
adapted_label="hx"
####ll
adapted2=["hg38","pantro5","gorgor5","panpan2","ponabe2","hlcascan1","hetgla2","hlfukdam1","loxafr3","hlturtru3","hlorcorc1","hldelleu1","lipvex1","phycat1","balacu1","hlbalmys1","hlrhisin1","hlhiparm1","eptfus1","myodav1","myobra1","myoluc2","hlminnat1","hldesrot1"]
adapted_label2="ll"
adapted_dict = {
  adapted_label: adapted,
  adapted_label2: adapted2,
}
adapted_labels = list(adapted_dict.keys())

#Mapping specific:
upstream_region = 50000
analysis_header = ["cis_reg_id", "query_species_only_count_del", "query_species_only_count_in"] + list_util.flatten_list([[f"{x}_del", f"{x}_in"] for x in adapted_labels])
indel_groups =['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY']

#Pipeline
MAX_BLOCKS = 144000 # if this is changed, number of blocks must be changed accordingly (file_parameters.py)
needs_preprocessing = False #some calculation needs only to be done once not for every min indel size

#Computational:
number_of_cores = 12