#############################################################################
import sys

try:
    import pandas as pd
except ImportError as e:
    sys.exit("Module pandas not found.")
try:
    import pyranges as pr
except ImportError as e:
    sys.exit("Module pyranges not found.")
try:
    from brex import region_extractor as rex
    from brex import helper as h
except ImportError as e:
    print(e)
    sys.exit("Module brex not found. Try running setup.py")

#############################################################################
##                     HIER DATEN EINFÜGEN                                 ##
#############################################################################
##                                                                         ##
##                       AUTOMATISCHES ERKENNEN                            ##
##                             ENSEMBL ID                                  ##
##                                                                         ##
#############################################################################
bio_mart_csv="mart_export.txt"
#gene_list="mouse_genes_list_reduced.txt"
#gene_list="genelists/anti_longevity_genes_list.txt"
gene_list="genelists/cross_val_short.txt"
#gene_list="genelists/pro_longevity_genes_list.txt"
#gene_list="genelists/dna_repair_genes.txt"
#gene_list="genelists/cancer_associated_genes.txt"
#gene_list="genelists/hypoxia_responsive_genes.txt"
upstream_region=50000
#############################################################################
##                                                                         ##
##                            OUTPUT                                       ##
##                                                                         ##
#############################################################################
##NAME THIS RUN 
#run_name="anti_longevity"
run_name="cross_val_short"
#run_name="pro_longevity"
#run_name="dna_repair"
#run_name="cancer_associated"
#run_name="hypoxia_responsive"
#############################################################################
##                                                                         ##
##                          MANUELLE EINGABE                               ##
##                         DES GENOMABSCHNITTS                             ##
##                                                                         ##
#############################################################################
## FÜGE HIER DIE SCAFFOLD DES GENOMABSCHNITTES EIN
#scaffold='chr8'
## FÜGE HIER DEN START DES GENOMABSCHNITTES EIN
#start=33300000 
## FÜGE HIER DAS ENDE DES GENOMABSCHNITTES EIN
#end=33400000
#############################################################################
##                                                                         ##
##                     MULTIPLE ALIGNMENT FILE                             ##
##                                                                         ##
#############################################################################
## PFAD ZUR MAF
#input_maf = "/moreAddSpace/abrzoska/mm10/group.maf"
input_maf = "/moreAddSpace/abrzoska/mm10/multi.anno.maf"
## REFERENZ SPECIES AUF DER DIE MAF BASIERT
target_species="mm10"
#############################################################################
##                                                                         ##
##                     QUERY SPECIES ETC                                   ##
##                                                                         ##
#############################################################################
##SPECIES DIE ANALYSIERT WERDEN SOLL
query_species="nanGal1"
##NAHE VERWANDTE DER SPEZIES 
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
##ORDNER AN DEM DIE BED/MAF ABGELEGT WIRD
#output_folder="../mouse"
output_folder="../" + run_name
#############################################################################
### FIND REGIONS ###
reg_ex = rex.RegionExtractor(run_name)  
reg_ex.start_create_bed_and_maf_for_list_of_genes(output_folder, gene_list, bio_mart_csv, upstream_region, input_maf, target_species, query_species, in_group, adapted_dict)
reg_ex.run_analysis(output_folder,adapted_labels,query_species,run_name)
##########################################################################
## UNVOLLSTÄNDIG: Filtern einer Bed-Datei anhand einer (oder mehrere) anderen, hier ENCODE
## input: output_names, array of file paths .append(reg_ex.find_bed_candidates(bed_file_path, bed_file_name, scaffold, start, end, output_dir))
# bfileEncode="../mt/data/encode_small.bed"
# output_names=["../pro_longevity/ENSMUSG00000054733-indels.bed_analysis.bed"]
# reg_ex.filter_bed_with_overlap(bfileEncode,output_names)
#############################################################################
### GET LIST OF SPECIES INSIDE A MAF ###
#reg_ex.get_species_in_maf(input_maf, "out.txt")
#############################################################################
### REDUCE A MAF TO LIST OF SPECIFIED SPECIES ###
#reg_ex.reduce_maf_via_list("output/mafsInRegion/wrn_species.txt","output/mafsInRegion/wrn.maf","output/mafsInRegion/wrn_reduced.maf")
#############################################################################