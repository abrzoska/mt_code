#############################################################################
import sys 
from time import perf_counter
from pathlib import Path
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
gene_list="genelists/anti_longevity_genes_list.txt"
#gene_list="genelists/pro_longevity_genes_list.txt"
#gene_list="genelists/dna_repair_genes.txt"
#gene_list="genelists/cancer_associated_genes.txt"
#gene_list="genelists/hypoxia_responsive_genes.txt"
upstream_region=50000
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
##NAHE VERWANDTE DER SPEZIES (INDELS DIE DIESEN TEILEN WERDEN IGNORIERT)
in_group= ["jacJac1","micOch1","criGri1","mesAur1","perManBai1","HLmusCar1","HLmusPah1","rn6","HLmerUng1","mm10"]
##LISTE ADAPTIERTER SPEZIES (INDELS DIE DIESEN TEILEN WERDEN HERVORGEHOBEN)
#hx
#adapted=["micOch1","criGri1","mesAur1","HLcasCan1","hetGla2","HLfukDam1","HLmarMar1","triMan1","chrAsi1","ornAna2","vicPac2","HLturTru3","HLorcOrc1","HLdelLeu1","lipVex1","phyCat1","balAcu1","HLbalMys1","bosMut1","panHod1","HLenhLut1","odoRosDiv1","lepWed1","neoSch1","conCri1"]
#ll
adapted=["hg38","panTro5","gorGor5","panPan2","ponAbe2","HLcasCan1","hetGla2","HLfukDam1","loxAfr3","HLturTru3","HLorcOrc1","HLdelLeu1","lipVex1","phyCat1","balAcu1","HLbalMys1","HLrhiSin1","HLhipArm1","eptFus1","myoDav1","myoBra1","myoLuc2","HLminNat1","HLdesRot1"]
##LABEL DAS DIE ART DER ADAPTION BESCHREIBT
adapted_label="ll"
#############################################################################
##                                                                         ##
##                            OUTPUT                                       ##
##                                                                         ##
#############################################################################
##NAME THIS RUN 
#run_name="ll_anti_longevity"
run_name="ll_pro_longevity"
##ORDNER AN DEM DIE BED/MAF ABGELEGT WIRD
#output_folder="../mouse"
output_folder="../" + run_name
#############################################################################
### FIND REGIONS ###
def get_output_name(line, output_folder):
    name_bed = "{0}/{1}-indels.bed".format(output_folder,line)
    name_maf = "{0}/{1}-indels.maf".format(output_folder,line)
    return [name_bed, name_maf]
def run_main():
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    gene_file = open(gene_list, "r")        
    df = pd.read_csv(bio_mart_csv) 
    reg_ex = rex.RegionExtractor(run_name)  
    for line in gene_file:
        gene_line_start_counter = perf_counter()
        line = line.strip()
        print("Running {0}".format(line))
        #print("Size of csv {0}".format(len(df.index)))
        gene = df[df['Gene_stable_ID'] == line]
        #print("Size of csv {0} after assignment".format(len(df.index)))
        if not gene.empty:#len(gene.index) > 0:
            #print(gene)
            column = gene["Chromosome_scaffold_name"]
            column = gene["Transcript_start"]
            latest_transcription_start = column.max()
            earliest_transcription_start = column.min()
            upstream_cutoff = earliest_transcription_start-upstream_region
            print("Upstream cutoff {0} - Latest Transcription start {1}".format(upstream_cutoff, latest_transcription_start))
            s = gene.loc[gene.index[0], 'Chromosome_scaffold_name']
            if s.isnumeric():
                scaffold = "chr{0}".format(s)
            else:
                scaffold = s        
            print("scaffold {0}".format(scaffold))
            b_header = f'track name="Indels for {query_species} {line} ({scaffold}:{upstream_cutoff}-{latest_transcription_start})" itemRgb="On"\n'
            [output_bed, output_maf] = get_output_name(line, output_folder)
            print("output_bed {0}".format(output_bed))
            reg_ex.run(upstream_cutoff, latest_transcription_start, scaffold, input_maf, target_species, query_species, adapted, in_group, output_maf, output_bed, b_header,adapted_label)
            print("...done")
        else: 
            print("Gene {0} not found in BioMart file, skipping".format(line))   
        gene_line_stop_counter = perf_counter()
        time = gene_line_stop_counter-gene_line_start_counter
        print('Elapsed time in seconds:", %.4f' % time)

reg_ex = rex.RegionExtractor(run_name)  
reg_ex.run_analysis(output_folder,adapted_label)

   
#############################################################################
### GET LIST OF SPECIES INSIDE A MAF ###
#reg_ex.get_species_in_maf(input_maf, "out.txt")
#############################################################################
### REDUCE A MAF TO LIST OF SPECIFIED SPECIES ###
#reg_ex.reduce_maf_via_list("output/mafsInRegion/wrn_species.txt","output/mafsInRegion/wrn.maf","output/mafsInRegion/wrn_reduced.maf")
#############################################################################