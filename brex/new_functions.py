import uuid

def make_species_string(species):
	return str(species).replace("\'", "\"")


def make_block_id(block_number, maf_id):
	return f"maf{maf_id}_block{block_number}"


def clean_and_write_entries(json_entries, file_out, last_entries=False):
	if len(json_entries) > 0:
		json_entries[-1] = json_entries[-1].rstrip(",")
	file_out.writelines(json_entries)
	if last_entries:
		file_out.write("]")
		file_out.close()


def make_indel_json_entry(chrom, start, end, len_with_gaps, len_no_gaps, indel_type, block_number, maf_id, agreeing_species, disagreeing_species):
	agreeing_species_string = make_species_string(agreeing_species)
	disagreeing_species_string = make_species_string(disagreeing_species)
	entry = f"""
	{{
	\"ID\": \"{uuid.uuid4()}\",
	\"chr\": \"{chrom}\",
	\"start\": {start},
	\"end\": {end},
	\"len_with_gaps\": {len_with_gaps},
	\"len_no_gaps\": {len_no_gaps},
	\"type\": \"{indel_type}\",
	\"block_id\": \"{make_block_id(block_number, maf_id)}\",
	\"agree\": {agreeing_species_string},
	\"disagree\": {disagreeing_species_string},
	}},"""
	return entry


def make_block_json_entry(block_number, maf_id, maf_score, maf_block_line_start, maf_block_line_end, mouse_sequence, len_no_gaps,
						  mm10_start, mm10_chr, species_count=0):
	if species_count == 0:
		entry = f"""
	{{
	\"ID\": \"{make_block_id(block_number, maf_id)}\",
	\"block_number\": {block_number},
	\"MAF_ID\": {maf_id},
	\"score\": {maf_score},
	\"file_pos_start\": {maf_block_line_start},
	\"file_pos_end\": {maf_block_line_end},
	\"seq\": \"{mouse_sequence.upper()}\",
	\"len_with_gaps\":\"{len(mouse_sequence)}\",
	\"len_no_gaps\":\"{len_no_gaps}\",
	\"mm10_start\":\"{mm10_start}\",
	\"mm10_chr\":\"{mm10_chr}\"
	}},"""
		return entry
	entry = f"""
	{{
	\"ID\": \"{make_block_id(block_number, maf_id)}\",
	\"block_number\": {block_number},
	\"MAF_ID\": {maf_id},
	\"score\": {maf_score},
	\"file_pos_start\": {maf_block_line_start},
	\"file_pos_end\": {maf_block_line_end},
	\"seq\": \"{mouse_sequence.upper()}\",
	\"species_count\": {species_count},
	\"len_with_gaps\":\"{len(mouse_sequence)}\",
	\"len_no_gaps\":\"{len_no_gaps}\",
	\"mm10_start\":\"{mm10_start}\",
	\"mm10_chr\":\"{mm10_chr}\"
	}},"""

	return entry
