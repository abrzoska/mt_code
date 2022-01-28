import sys 
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
output_dir = "../mt/data"
# #chr8	33300000 33400000
# start=33300000 
# end=33400000
# scaffold='chr8'
# reg_ex = rex.RegionExtractor()
#############################################################################
##      TODO diff DS?
#############################################################################
# bfileEncode="../mt/data/encode_small.bed"
# bfileTf="../mt/data/wrn_small_tab.bed"
# bfileCNE="../mt/data/mouseCNEs.bed"
# bfileIndels="../mt/output/mafsInRegion/wrn_nanGal1_indels.bed"
# #bfiles = {"tfmotifView": bfileTf, "CNE": bfileCNE, "indels": bfileIndels}
# bfiles = {"tfmotifView": bfileTf, "CNE": bfileCNE}
#############################################################################
##      TODO    #mafsInRegion chr/chr8.bed chr/chr8.maf multi.anno.maf
##      Assert that the right chr maf is present
#############################################################################
# maffile="/moreAddSpace/abrzoska/mm10/chr/" + scaffold + ".maf"
#############################################################################
##      TODO reduce the maf to only the species of interest
#############################################################################
# reduced_maf = "../mt/output/mafsInRegion/wrn_reduced.maf"
reduced_maf = "/moreAddSpace/abrzoska/mm10/multi.anno.maf"
#############################################################################
### FIND REGIONS ###
reg_ex = rex.RegionExtractor()
output_names = []
#indel_reduced_maf = "../mt/output/mafsInRegion/wrn_nanGal1_indels.maf"
#indel_reduced_bed = "../mt/output/mafsInRegion/wrn_nanGal1_indels.bed"
# indel_reduced_maf = "../output/wrn_nanGal1_indels.maf"
# indel_reduced_bed = "../output/wrn_nanGal1_indels.bed"
#in_group= ["mm10", "rn6"]
# in_group= ["hg38"]
# #in_group=[]
# adapted=["myoDav1", "myoLuc2", "myoBra1", "myoDav1"]
# target_species="mm10"
# query_species="nanGal1"
##      create bed files for range for main and every file
#############################################################################
#for bfile_name in bfiles:
#    output_names.append(reg_ex.find_bed_candidates(bfiles[bfile_name], bfile_name, scaffold, start, end, output_dir))
#main_name = reg_ex.find_bed_candidates(bfileEncode, "ENCODE",  scaffold, start, end, output_dir)
#print(output_names)
#print(main_name)
#############################################################################
##      Input: gene id, data files
#############################################################################
#reg_ex.annotate_bed_with_overlap(main_name,output_names)
#############################################################################
##      TODO get maf and bed representation for query
##      maf
##      mm10, nanGal1
##      in_group
##      returns name of out files. maf, bed
#############################################################################
#reg_ex.run(start, end, scaffold, reduced_maf, target_species, query_species, adapted, in_group, indel_reduced_maf, indel_reduced_bed)
#reg_ex.get_species_in_maf(reduced_maf, "out-allspecies.txt")
#reg_ex.reduce_maf_via_list("group.txt","/moreAddSpace/abrzoska/mm10/multi.anno.maf","/moreAddSpace/abrzoska/mm10/group.maf")
#input_maf = 
#############################################################################