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
#ensembl_id="ENSMUSG00000086493"
bio_mart_csv="mart_export.txt"
gene_list="mouse_genes_list_reduced.txt"
upstream_region=50000
#############################################################################
##                                                                         ##
##                     ORT DES GENOMABSCHNITTS                             ##
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
input_maf = "/moreAddSpace/abrzoska/mm10/group.maf"
## SPECIES AUF DER DIE MAF BASIERT
target_species="mm10"
#############################################################################
##                                                                         ##
##                     QUERY SPECIES ETC                                   ##
##                                                                         ##
#############################################################################
##SPECIES DIE ANALYSIERT WERDEN SOLL
query_species="nanGal1"
##NAHE VERWANDTE DER SPEZIES (INDELS DIE DIESEN TEILEN WERDEN IGNORIERT)
in_group= ["mm10", "rn6"]
##LISTE ADAPTIERTER SPEZIES (INDELS DIE DIESEN TEILEN WERDEN HERVORGEHOBEN)
adapted=["hetGla2"]
#############################################################################
##                                                                         ##
##                            OUTPUT                                       ##
##                                                                         ##
#############################################################################
##ORDNER AN DEM DIE BED/MAF ABGELEGT WIRD
output_folder="../mouse"
#############################################################################
### FIND REGIONS ###
def get_output_name(line, output_folder):
    name_bed = "{0}/{1}-indels.bed".format(output_folder,line)
    name_maf = "{0}/{1}-indels.maf".format(output_folder,line)
    return [name_bed, name_maf]
gene_file = open(gene_list, "r")        
df = pd.read_csv(bio_mart_csv) 
for line in gene_file:
    line = line.strip()
    reg_ex = rex.RegionExtractor()
    print("Running {0}".format(line))
    #print("Size of csv {0}".format(len(df.index)))
    gene = df[df['Gene_stable_ID'] == line] 
    #print("Size of csv {0} after assignment".format(len(df.index)))
    if not gene.empty:#len(gene.index) > 0:
        #print(gene)
        column = gene["Chromosome_scaffold_name"]
        column = gene["Transcript_start"]
        end = column.max()
        start = end-upstream_region
        print("start {0} - end {1}".format(start, end))
        s = gene.loc[gene.index[0], 'Chromosome_scaffold_name']
        if s.isnumeric():
            scaffold = "chr{0}".format(s)
        else:
            scaffold = s        
        print("scaffold {0}".format(scaffold))
        b_header = f'track name="Indels for {query_species} {line} ({scaffold}:{start}-{end})" itemRgb="On"\n'
        [output_bed, output_maf] = get_output_name(line, output_folder)
        print("output_bed {0}".format(output_bed))
        reg_ex.run(start, end, scaffold, input_maf, target_species, query_species, adapted, in_group, output_maf, output_bed, b_header)
        print("...done")
    else: 
        print("Gene {0} not found in BioMart file, skipping".format(line))   
   
#############################################################################
### GET LIST OF SPECIES INSIDE A MAF ###
#reg_ex.get_species_in_maf(reduced_maf, "out.txt")
#############################################################################
### REDUCE A MAF TO LIST OF SPECIFIED SPECIES ###
#reg_ex.reduce_maf_via_list("output/mafsInRegion/wrn_species.txt","output/mafsInRegion/wrn.maf","output/mafsInRegion/wrn_reduced.maf")
#############################################################################