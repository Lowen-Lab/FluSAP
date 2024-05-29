import os
import sys
import numpy as np
import scipy.stats


parallel_process = True
parallel_max_cpu = 50


skip_writing_all_site_summaries = True

output_suffix = "wgs"

###################################   USER DEFINED VARIABLES   ###################################
project_dir = "./"
input_read_dir =project_dir+"raw_data/"
ref_dir = project_dir+"refs/"
trim_dir = project_dir+"trim/"
map_folder = 'rep_map'
map_dir = project_dir+map_folder+"/"
qual_dir = map_dir+"qual/"
sam_dir = map_dir+"sam/"
bam_dir = map_dir+"bam/"
sub_dir = map_dir+"substitutions/"
temp_dir = map_dir+"temp/"
consensus_dir = map_dir+"consensus/"

annotation_file_path = ref_dir + "annotation.tbl"

ref_set_list = ['IAV_H1.fa','IAV_H3.fa','IAV_M.fa','IAV_N1.fa','IAV_N2.fa','IAV_NP.fa','IAV_NS.fa','IAV_PA.fa','IAV_PB1.fa','IAV_PB2.fa','IBV_HA.fa','IBV_M.fa','IBV_NA.fa','IBV_NP.fa','IBV_NS.fa','IBV_PA.fa','IBV_PB1.fa','IBV_PB2.fa']


sigma_flank_len = 100
min_freq_sigma_score = 3.0
min_observations = 50
min_proportion = 0.001
min_qual = 35.0
min_avg_map_qual = 30.0
min_avg_read_loc = 25.0
max_avg_read_mismatch = 3.0
max_avg_read_indel = 2.0
min_major_map_qual = min_avg_map_qual

min_cov = 500
min_cov_pass_site_prop = 0.5
end_skip_len = 10

min_summary_cov = 1

############################################## SETUP ##############################################
ref_list = ref_set_list
segment_list = []
for ref_set in ref_set_list:
	refs = ref_set.split("-")
	segment = ref_set.split("-")[0].split(".f")[0]
	segment_list.append(segment)
segment_list = sorted(segment_list)


amino_dict = {'A':['GCT','GCC','GCA','GCG'],
'R':['CGT','CGC','CGA','CGG','AGA','AGG'],
'N':['AAT','AAC'],
'D':['GAT','GAC'],
'C':['TGT','TGC'],
'Q':['CAA','CAG'],
'E':['GAA','GAG'],
'G':['GGT','GGC','GGA','GGG'],
'H':['CAT','CAC'],
'I':['ATT','ATC','ATA'],
'L':['TTA','TTG','CTT','CTC','CTA','CTG'],
'K':['AAA','AAG'],
'M':['ATG'],
'F':['TTT','TTC'],
'P':['CCT','CCC','CCA','CCG'],
'S':['TCT','TCC','TCA','TCG','AGT','AGC'],
'T':['ACT','ACC','ACA','ACG'],
'W':['TGG'],
'Y':['TAT','TAC'],
'V':['GTT','GTC','GTA','GTG'],
'*':['TAG','TGA','TAA']}
codon_to_aa_dict = {}
for aa in amino_dict:
	codons = amino_dict[aa]
	for codon in codons:
		codon_to_aa_dict[codon] = aa


sites_to_exclude = {}
samples_to_exclude = {}
try:
	filename = "samples_to_exclude.txt"
	infile = open(project_dir+filename,"r")
	for line in infile:
		line = line.strip()
		samples_to_exclude[line] = ''
	infile.close
except:
	pass

if parallel_process == True:
	import multiprocessing
	from joblib import Parallel, delayed
	if parallel_max_cpu == 0:
		num_cores = multiprocessing.cpu_count()
	elif parallel_max_cpu > 0:
		num_cores = min(parallel_max_cpu,multiprocessing.cpu_count())
	elif parallel_max_cpu < 0:
		num_cores = multiprocessing.cpu_count()+parallel_max_cpu
	print("Number of cores requested: "+str(num_cores))
else:
	num_cores = 1

############################################## Pull sample list from sequencing files ##############################################

def pull_sample_list(temp_dir,input_read_dir):
	forward_read_suffix = ''
	reverse_read_suffix = ''
	suffix_info_file_path = temp_dir+"file_extensions.txt"
	try:
		file_suffix_infile = open(suffix_info_file_path)
		for line in file_suffix_infile:
			line = line.strip().split("\t")
			if line[0] == "forward":
				forward_read_suffix = line[1]
			elif line[0] == "reverse":
				reverse_read_suffix = line[1]
			continue_running = True
	except:
		sys.exit("Exiting.")

	file_list = [f for f in os.listdir(input_read_dir) if f.endswith(forward_read_suffix)]
	sampleID_list = []
	for files in file_list:
		sampleID = files.split(forward_read_suffix)[0]
		sampleID_list.append(sampleID)
	sampleID_list = list(set(sampleID_list))
	return sampleID_list


sampleID_list = pull_sample_list(temp_dir,input_read_dir)
flagged_sample_list = []
unique_sample_dict = {}
for sampleID in sampleID_list:
	try:
		samples_to_exclude[sampleID]
	except:
		flagged_sample_list.append(sampleID)
		sample_root = sampleID.split("_run")[0]
		try:
			unique_sample_dict[sample_root].append(sampleID)
		except:
			unique_sample_dict[sample_root] = []
			unique_sample_dict[sample_root].append(sampleID)
sampleID_list = []
for sample_root in unique_sample_dict:
	if len(unique_sample_dict[sample_root]) >=1:
		for sampleID in unique_sample_dict[sample_root]:
			sampleID_list.append(sampleID)
del flagged_sample_list
sample_root_list = list(unique_sample_dict.keys())

print("Found "+str(len(sampleID_list))+" files, with "+str(len(unique_sample_dict))+" unique sample root IDs.")

############################################ FUNCTIONS ############################################
def run_command(command,mode="quiet"):
	run_status = False
	if mode == "verbose":
		return_code = os.system(command)
	elif mode == "quiet":
		return_code = os.system(command+" >/dev/null 2>&1")
	if return_code == 0:
		run_status = True
	return run_status

def nucleotide_character_edit(seqin):#,ambiguous_nucleotide_list):
	ambiguous_nucleotide_list = ['B','b','D','d','H','h','K','k','M','m','R','r','S','s','V','v','W','w','Y','y']
	seqin = seqin.replace("a","A")
	seqin = seqin.replace("c","C")
	seqin = seqin.replace("g","G")
	seqin = seqin.replace("t","T")
	seqin = seqin.replace("u","T")
	seqin = seqin.replace("U","T")
	seqin = seqin.replace("n","N")
	for code in ambiguous_nucleotide_list:
		seqin = seqin.replace(code,"N")
	return seqin

def is_number(s):
	try:
		float(s)
	except ValueError:
		return False
	return True

def round_to_n_sig_figs(val,num_sig_figs):
	if is_number(val):
		val = float(val)
		if val == 0.0:
			return '0.0'
		num_sig_figs = int(num_sig_figs)
		if num_sig_figs == 0:
			num_sig_figs = 1
		sci_val = "{:.10e}".format(val)
		split_sci_val = sci_val.split("e")
		if len(split_sci_val) == 2:
			rounded_base_number = round(float(split_sci_val[0]),num_sig_figs-1)
			exponent = int(split_sci_val[1])
			if exponent == 0:
				val_out = str(rounded_base_number) + ((num_sig_figs)-1)*'0'
			elif exponent < 0:
				exponent*=-1
				val_out = '0.' + (exponent-1)*'0' + str(rounded_base_number).replace(".","")
				val_out = str(float(val_out))
			elif exponent > 0:
				val_out = str(rounded_base_number) +'e'+ (str(exponent))
			return val_out
		else:
			return val
	else:
		return val

def nt_diversity(nt_cov_list):
	N = 0.0
	D = 0.0
	if len(nt_cov_list) > 1:
		for nt in range(0,len(nt_cov_list)):
			nt_cov = int(float(nt_cov_list[nt]))
			if nt_cov > 0.0:
				D += nt_cov*(nt_cov-1.0)
				N += nt_cov
		if N > 0.0:
			D = ((N*(N-1.0)) - D)/((N*(N-1.0)))
	return D

def split_multi_space(string_in):
	string_out = string_in
	keep_running = True
	num_replaced = 0
	while keep_running == True:
		len_before = len(string_out)
		string_out = string_out.replace('  ',' ')
		len_after = len(string_out)
		num_replaced = len_before - len_after
		if num_replaced == 0:
			keep_running = False
	list_out = string_out.split(' ')
	return list_out

def find_cds_locs(sampleID,segment_list):
	cds_loc_num_dict = {}
	cds_num_info_filename = ref_dir+'cds_num_info.txt'
	cds_num_info = open(cds_num_info_filename,'r')
	for line in cds_num_info:
		line = line.strip().split("\t")
		ref_segment = line[0]
		cds_ID = line[1]
		for i in range(2,len(line)):
			loc_range = line[i].split(",")
			start_loc_num = int(loc_range[0])
			stop_loc_num = int(loc_range[1])
			try:
				cds_loc_num_dict[ref_segment][cds_ID].append((start_loc_num,stop_loc_num))
			except:
				try:
					cds_loc_num_dict[ref_segment][cds_ID] = [(start_loc_num,stop_loc_num)]
				except:
					cds_loc_num_dict[ref_segment] = {}
					cds_loc_num_dict[ref_segment][cds_ID] = [(start_loc_num,stop_loc_num)]
	cds_num_info.close()
	annot_dict = {}
	for s in range(0,len(segment_list)):
		segmentID = segment_list[s]
		segment_seqfile_path = consensus_dir+segmentID+'/'+sampleID+'.'+segmentID+'.refine.fa'
		if os.path.isfile(segment_seqfile_path) == True:
			cds_loc_file = consensus_dir+segmentID+'/'+sampleID+'.'+segmentID+'.refine.cds_loc.txt'
			if os.path.isfile(cds_loc_file) == True:
				annot_infile = open(cds_loc_file,"r")
				for line in annot_infile:
					line = line.strip().split('\t')
					ref_segment = line[0]
					cds_ID = line[1]
					for i in range(2,len(line)):
						loc_range = line[i].split(",")
						start_loc = int(loc_range[0])
						stop_loc = int(loc_range[1])
						try:
							annot_dict[segmentID][cds_ID].append((start_loc,stop_loc))
						except:
							try:
								annot_dict[segmentID][cds_ID] = [(start_loc,stop_loc)]
							except:
								annot_dict[segmentID] = {}
								annot_dict[segmentID][cds_ID] = [(start_loc,stop_loc)]
				annot_infile.close()

			else:
				command = 'hmmscan --domE 1e-5 --domtblout '+temp_dir+sampleID+'.'+segmentID+'.dtbloud.txt --notextw ./refs/cds_hmm/'+segmentID+'.cds_loc.hmm '+consensus_dir+segmentID+'/'+sampleID+'.'+segmentID+'.refine.fa'
				run_status = run_command(command)
				if run_status == False:
					sys.exit("cds hmmscan failed: "+command)

				hmm_loc_dict = {}
				hmm_dtblout_filename = temp_dir+sampleID+'.'+segmentID+'.dtbloud.txt'
				annot_infile = open(hmm_dtblout_filename,'r')
				for line in annot_infile:
					line = line.rstrip('\n')
					if len(line) >0:
						if line[0]!= '#':
							split = split_multi_space(line)
							cds_loc_ID = split[0]
							cds_loc_num = int(cds_loc_ID.split('-')[1])
							hmm_start = int(split[15])
							hmm_end = int(split[16])
							align_start = int(split[19])
							align_end = int(split[20])
							if cds_loc_num == 0:
								if cds_loc_ID in hmm_loc_dict:
									print(line)
									print(hmm_loc_dict[cds_loc_ID])
								hmm_loc_dict[cds_loc_ID] = align_start
							else:
								hmm_loc_dict[cds_loc_ID] = align_end
				annot_infile.close()
				
				outlines = ''
				cur_line = ''
				for cds_ID in cds_loc_num_dict[segmentID]:
					loc_num_list = cds_loc_num_dict[segmentID][cds_ID]
					cur_line = segmentID+'\t'+cds_ID
					for i in range(0,len(loc_num_list)):
						loc_num1 = loc_num_list[i][0]
						loc_num2 = loc_num_list[i][1]
						cds_loc_ID1 = segmentID+'-'+str(loc_num1)
						cds_loc_ID2 = segmentID+'-'+str(loc_num2)
						try:
							loc1 = hmm_loc_dict[cds_loc_ID1]
							loc2 = hmm_loc_dict[cds_loc_ID2]
						except:
							loc1 = ''
							loc2 = ''
						if loc1 != '' and loc1 != '':
							try:
								annot_dict[segmentID][cds_ID].append((loc1,loc2))
							except:
								try:
									annot_dict[segmentID][cds_ID] = [(loc1,loc2)]
								except:
									annot_dict[segmentID] = {}
									annot_dict[segmentID][cds_ID] = [(loc1,loc2)]
							cur_line += '\t'+str(loc1)+','+str(loc2)
						else:
							cur_line = ''
					if cur_line != '':
						outlines += cur_line+'\n'

				temp_cds_loc_file = temp_dir+sampleID+'.'+segmentID+'.refine.cds_loc.txt'
				cds_loc_file = consensus_dir+segmentID+'/'+sampleID+'.'+segmentID+'.refine.cds_loc.txt'
				annot_outfile = open(temp_cds_loc_file,"w")
				annot_outfile.write(outlines)
				annot_outfile.close()
				os.rename(temp_cds_loc_file,cds_loc_file)
	return annot_dict


def pull_cds_from_annot_file(annotation_filename):
	annot_dict = {}
	found_CDS = False
	cur_cds_locs = []
	cur_product_ID = ''
	cur_gene_ID = ''
	annot_infile = open(annotation_filename,"r")
	for line in annot_infile:
		line = line.strip("\n")
		if len(line.strip())>0:
			if line[0]==">":
				if line[0:8]==">Feature":
					contig = line.split(">Feature ")[1]
				else:
					contig = line[1:len(line)]
				found_CDS = False
				cur_cds_locs = []
				cur_product_ID = ''
				cur_gene_ID = ''
			elif found_CDS == False:
				if is_number(line[0]) == True:
					split = line.split("\t")
					product_type = split[2]
					if product_type == "CDS":
						if is_number(split[0]) == True and is_number(split[1]) == True:
							found_CDS = True
							start = int(split[0])
							stop = int(split[1])
							cur_cds_locs = []
							cur_product_ID = ''
							cur_gene_ID = ''
							cur_cds_locs.append((start,stop))
			elif found_CDS == True:
				if is_number(line[0]) == True:
					split = line.split("\t")
					if split[2] == '':
						if is_number(split[0]) == True and is_number(split[1]) == True:
							start = int(split[0])
							stop = int(split[1])
							cur_cds_locs.append((start,stop))
						else:
							cur_cds_locs = []
					else:
						found_CDS = False
				elif line[0] == "\t":
					if line[2] == "\t":
						split = line.split("\t")
						if split[3] == "protein_id":
							cur_product_ID = split[4]
						elif split[3] == "gene":
							cur_gene_ID = split[4]
							if cur_cds_locs != []:
								try:
									annot_dict[contig][cur_gene_ID] = cur_cds_locs
								except:
									annot_dict[contig] = {}
									annot_dict[contig][cur_gene_ID] = cur_cds_locs
							found_CDS = False
	return annot_dict

def codon_finder(site, coding_start):
	codon_number = int(np.floor((site-coding_start)/3.0))
	codon_site = int((((float((site-coding_start))/3.0)-float((np.floor((site-coding_start)/3.0))))/0.333333+0.1))
	return codon_number,codon_site

def sub_site_count(seqin, coding_start, coding_stop,codon_to_aa_dict):
	nt_list = ['A','T','C','G']
	viable_codon_count = 0.0
	ns_count = 0.0
	s_count = 0.0
	for loc in range(coding_start,coding_stop+1):
		codon_number,codon_site = codon_finder(loc, coding_start)
		codon_bases = seqin[(codon_number*3+coding_start):(codon_number*3+coding_start)+3]
		if "-" not in codon_bases and "N" not in codon_bases and len(codon_bases)==3:
			viable_codon_count += 1.0
			for nt in nt_list:
				if nt != seqin[loc]:
					sub_codon = list(codon_bases)
					sub_codon[codon_site] = nt
					sub_codon = "".join(sub_codon)
					ref_aa = codon_to_aa_dict[codon_bases]
					sub_aa = codon_to_aa_dict[sub_codon]
					if ref_aa != sub_aa:
						ns_count += 0.3333
					elif ref_aa == sub_aa:
						s_count += 0.3333
	if viable_codon_count > 0:
		ns_prop = round(float(ns_count)/float(viable_codon_count),3)
		s_prop = round(float(s_count)/float(viable_codon_count),3)
	else:
		ns_prop = 0
		s_prop = 0
	return ns_prop,s_prop,round(ns_count,2),round(s_count,2),int(viable_codon_count)

def translate_cds(seqin, coding_start, coding_stop,codon_to_aa_dict):
	aa_seq_out = ''
	for loc in range(coding_start,coding_stop+1,3):
		codon_number,codon_site = codon_finder(loc, coding_start)
		codon_bases = seqin[(codon_number*3+coding_start):(codon_number*3+coding_start)+3]
		if "-" not in codon_bases and "N" not in codon_bases and len(codon_bases)==3 and "n" not in codon_bases:
			aa = codon_to_aa_dict[codon_bases]
		elif codon_bases != '':
			aa = 'X'
		else:
			aa = ''
		aa_seq_out += aa
	return aa_seq_out

def calc_lognormal_array(freq_list):
	lognorm_freq_list = []
	for i in range(0,len(freq_list)):
		lognorm_val = -1*np.log10(freq_list[i])
		lognorm_freq_list.append(lognorm_val)
	return lognorm_freq_list

def read_in_allele_info(sampleID,segment_list,map_dir,consensus_dir):
	allele_out_dict = {}
	refseq_out_dict = {}
	obs_consensus_seq_out_dict = {}
	sub_dict = {}
	cov_pass_count_dict = {}
	bases = ['A','T','C','G']
	for s in range(0,len(segment_list)):
		segment = segment_list[s]
		sub_dict[segment] = []
		cov_pass_count_dict[segment] = 0
		substitutions_filename = map_dir+"substitutions/"+segment+"/"+sampleID+'.'+segment+".substitions.txt"
		refseq_filename = consensus_dir+segment+'/'+sampleID+'.'+segment+'.refine.fa'
		if os.path.isfile(refseq_filename) == True:
			ref_seq = ''
			infile = open(refseq_filename,'r')
			for line in infile:
				line = line.strip()
				if line[0]!=">":
					ref_seq += line
			infile.close()
			refseq_out_dict[segment] = ref_seq


		if os.path.isfile(substitutions_filename) == True:
			obs_consensus_seq_dict = {}
			sub_infile = open(substitutions_filename,"r")
			for sub_line in sub_infile:
				if len(sub_line) > 0:
					if sub_line[0] != "#" and sub_line[0] != "loc":
						sub_line = sub_line.strip().split("\t")
						loc = int(sub_line[0])
						cov = float(sub_line[1])
						A_prop,T_prop,C_prop,G_prop = float(sub_line[2]),float(sub_line[3]),float(sub_line[4]),float(sub_line[5])
						A_qual,T_qual,C_qual,G_qual = float(sub_line[6]),float(sub_line[7]),float(sub_line[8]),float(sub_line[9])
						A_map_qual,T_map_qual,C_map_qual,G_map_qual = float(sub_line[10]),float(sub_line[11]),float(sub_line[12]),float(sub_line[13])
						A_read_len,T_read_len,C_read_len,G_read_len = float(sub_line[14]),float(sub_line[15]),float(sub_line[16]),float(sub_line[17])
						A_read_loc,T_read_loc,C_read_loc,G_read_loc = float(sub_line[18]),float(sub_line[19]),float(sub_line[20]),float(sub_line[21])
						A_read_R_loc,T_read_R_loc,C_read_R_loc,G_read_R_loc = float(sub_line[22]),float(sub_line[23]),float(sub_line[24]),float(sub_line[25])
						A_read_prop_loc,T_read_prop_loc,C_read_prop_loc,G_read_prop_loc = float(sub_line[26]),float(sub_line[27]),float(sub_line[28]),float(sub_line[29])
						A_read_mismatch,T_read_mismatch,C_read_mismatch,G_read_mismatch = float(sub_line[30]),float(sub_line[31]),float(sub_line[32]),float(sub_line[33])
						A_read_indel,T_read_indel,C_read_indel,G_read_indel = float(sub_line[34]),float(sub_line[35]),float(sub_line[36]),float(sub_line[37])

						info_tup = (cov,(A_prop,T_prop,C_prop,G_prop),
							(A_qual,T_qual,C_qual,G_qual),
							(A_map_qual,T_map_qual,C_map_qual,G_map_qual),
							(A_read_len,T_read_len,C_read_len,G_read_len),
							(A_read_loc,T_read_loc,C_read_loc,G_read_loc),
							(A_read_R_loc,T_read_R_loc,C_read_R_loc,G_read_R_loc),
							(A_read_prop_loc,T_read_prop_loc,C_read_prop_loc,G_read_prop_loc),
							(A_read_mismatch,T_read_mismatch,C_read_mismatch,G_read_mismatch),
							(A_read_indel,T_read_indel,C_read_indel,G_read_indel))
						try:
							allele_out_dict[segment][loc] = info_tup
						except:
							allele_out_dict[segment] = {}
							allele_out_dict[segment][loc] = info_tup
						
						sort_prop_nt_list = sorted([(A_prop,"A"),(T_prop,"T"),(C_prop,"C"),(G_prop,"G")],reverse=True)
						major_nt = sort_prop_nt_list[0][1]
						major_prop = sort_prop_nt_list[0][0]
						if cov >= 3 and major_prop >= 0.5:
							obs_consensus_seq_dict[loc] = major_nt
						else:
							obs_consensus_seq_dict[loc] = 'N'
						if cov >= min_cov:
							cov_pass_count_dict[segment] += 1
						
							prop_list = [A_prop,T_prop,C_prop,G_prop]
							qual_list = [A_qual,T_qual,C_qual,G_qual]
							map_qual_list = [A_map_qual,T_map_qual,C_map_qual,G_map_qual]
							read_loc_list = [A_read_loc,T_read_loc,C_read_loc,G_read_loc]
							mismatch_list = [A_read_mismatch,T_read_mismatch,C_read_mismatch,G_read_mismatch]
							indel_list = [A_read_indel,T_read_indel,C_read_indel,G_read_indel]

							local_min_prop = 1.0/float(cov)
							local_null_lognormprop = -1*np.log10(1.0/float(cov))

							prop_vals = []
							for b in range(0,len(prop_list)):
								focal_nt = bases[b]
								focal_prop = prop_list[b]
								if focal_nt != major_nt and focal_prop <= 0.5:
									focal_qual = qual_list[b]
									focal_map_qual = map_qual_list[b]
									focal_read_loc = read_loc_list[b]
									focal_mismatch = mismatch_list[b]
									focal_indel = indel_list[b]
									if focal_prop>=local_min_prop:
										focal_lognorm_prop = -1*np.log10(focal_prop)
										prop_vals.append(focal_lognorm_prop)

							if len(prop_vals) == 1:
								sub_dict[segment].append(prop_vals[0])
							elif len(prop_vals) > 1:
								prop_vals = sorted(prop_vals,reverse=False)
								sub_dict[segment].append(prop_vals[0])
							else:
								sub_dict[segment].append(local_null_lognormprop)
			try:
				refseq_len = len(refseq_out_dict[segment])
			except:
				refseq_len = 0
			if refseq_len>0 and obs_consensus_seq_dict != {}:
				obs_consensus_seq = ''
				for loc in range(0,refseq_len):
					try:
						base = obs_consensus_seq_dict[loc]
					except:
						base = 'N'
					obs_consensus_seq += base
				obs_consensus_seq_out_dict[segment] = obs_consensus_seq
	return (sampleID,allele_out_dict,refseq_out_dict,sub_dict,cov_pass_count_dict,obs_consensus_seq_out_dict)



############################################   Read in allele info   #############################################
def map_sites_between_replicate_ref_seqs(sample_root,replicate_sample_list):
	min_align_prop = 0.95
	min_align_pid = 0.95
	sample_loc_dict = {}
	replicate_sample_list = sorted(replicate_sample_list,reverse=False)
	if len(replicate_sample_list)>1:
		reference_replicate = replicate_sample_list[0]
		rep_filepath = consensus_dir+"all_segments/"+reference_replicate+'.refine.fa'
		for i in range(1,len(replicate_sample_list)):
			next_sample = replicate_sample_list[i]
			next_filepath = consensus_dir+"all_segments/"+next_sample+'.refine.fa'
			if os.path.isfile(rep_filepath)==True and os.path.isfile(next_filepath)==True:
				blast_outfile_name = temp_dir+reference_replicate+"."+next_sample+'.refine.blastn.txt'
				command = "blastn -query "+rep_filepath+" -subject "+next_filepath+' -evalue 1E-10 -outfmt "6 qseqid sseqid pident evalue qlen slen length qstart qend sstart send qseq sseq" -out '+blast_outfile_name
				run_status = run_command(command)
				if run_status == False:
					sys.exit("blastn failed: "+command)
				infile = open(blast_outfile_name,"r")
				for line in infile:
					split = line.strip().split("\t")
					'''
					   0      1      2      3    4    5      6     7      8    9     10   11   12
					qseqid sseqid pident evalue qlen slen length qstart qend sstart send qseq sseq
					'''
					qseqid = split[0]
					sseqid = split[1]
					qseq = split[11]
					sseq = split[12]
					pair = qseqid+"\t"+sseqid
					pident = float(split[2])/100.0 #% identity
					qlen = int(split[4]) #length of query
					slen = int(split[5]) #length of subject 
					length = int(split[6]) #length of alignment
					qstart = int(split[7])-1 #align start site in query
					qstop = int(split[8])#-1 #align stop site in query
					sstart = int(split[9])-1 #align start site in subject
					sstop = int(split[10])#-1 #align stop site in subject

					perclength = float(length)/float(max(qlen, slen))
					if perclength >= min_align_prop and pident >= min_align_pid:
						q_loc = qstart
						s_loc = sstart
						for i in range(0,len(qseq)):
							q_base = qseq[i]
							s_base = sseq[i]
							if q_base != '-':
								q_loc += 1
							if s_base != '-':
								s_loc += 1
							if q_base != '-' and s_base != '-':
								try:
									sample_loc_dict[sseqid][s_loc] = q_loc
								except:
									sample_loc_dict[sseqid] = {}
									sample_loc_dict[sseqid][s_loc] = q_loc
					try:
						sample_loc_dict[qseqid]
					except:
						sample_loc_dict[qseqid] = {}
						for loc in range(0,qlen+1):
							sample_loc_dict[qseqid][loc] = loc
				infile.close()
				file_rm_command = 'rm '+blast_outfile_name
				run_command(file_rm_command)
	return (sample_root,sample_loc_dict)

processed_list = []
if parallel_process == True:
	processed_list = Parallel(n_jobs=num_cores)(delayed(map_sites_between_replicate_ref_seqs)(sample_root,unique_sample_dict[sample_root]) for sample_root in sample_root_list)
elif parallel_process == False:
	for sampleID in sampleID_list:
		output = map_sites_between_replicate_ref_seqs(sample_root,unique_sample_dict[sample_root])
		processed_list.append(output)

link_loc_dict = {}
for tup in processed_list:
	sample_root = tup[0]
	link_loc_dict[sample_root] = tup[1]
del processed_list

print("Finished mapping nucleotide positions between replicate samples.")



processed_list = []
if parallel_process == True:
	processed_list = Parallel(n_jobs=num_cores)(delayed(read_in_allele_info)(sampleID,segment_list,map_dir,consensus_dir) for sampleID in sampleID_list)
elif parallel_process == False:
	for sampleID in sampleID_list:
		output = read_in_allele_info(sampleID,segment_list,map_dir,consensus_dir)
		processed_list.append(output)

max_length_dict = {}
all_allele_info_dict = {}
consensus_seq_dict = {}
background_freq_dict = {}
cov_pass_dict = {}
sample_consensus_seq_dict = {}
for tup in processed_list:
	sampleID = tup[0]
	all_allele_info_dict[sampleID] = tup[1]
	consensus_seq_dict[sampleID] = tup[2]
	background_freq_dict[sampleID] = tup[3]
	cov_pass_dict[sampleID] = tup[4]
	temp_sample_consensus_seq_dict = tup[5]
	for segment in temp_sample_consensus_seq_dict:
		acc_seg = sampleID+"_"+segment
		sample_consensus_seq_dict[acc_seg] = temp_sample_consensus_seq_dict[segment]
	for segment in consensus_seq_dict[sampleID]:
		seq_len = len(consensus_seq_dict[sampleID][segment])
		try:
			max_seqlen = max_length_dict[segment]
			if seq_len > max_seqlen:
				max_length_dict[segment] = seq_len
		except:
			max_length_dict[segment] = seq_len
del processed_list
# print(sample_consensus_seq_dict)

print("Finished reading in allele info from substitution files.")


all_cds_loc_dict = pull_cds_from_annot_file(annotation_file_path)

cds_first_loc_dict = {}
cds_loc_string = ''
all_cds_site_loc_dict = {}
all_cds_refs_dict = {}
for acc_seq in all_cds_loc_dict:
	try:
		seg_seq = sample_consensus_seq_dict[acc_seq]
	except:
		seg_seq = ''
	if seg_seq != '':
		local_loc_list = []
		all_cds_site_loc_dict[acc_seq] = {}
		all_cds_refs_dict[acc_seq] = {} 
		for CDS in all_cds_loc_dict[acc_seq]:
			ref_cds = ''
			all_cds_site_loc_dict[acc_seq][CDS] = {}
			cds_loc_string+= acc_seq+'\t'+CDS
			for CDS_region in range(0,len(all_cds_loc_dict[acc_seq][CDS])):
				cur_cds_loc = len(ref_cds)
				cds_start = all_cds_loc_dict[acc_seq][CDS][CDS_region][0]
				cds_stop = all_cds_loc_dict[acc_seq][CDS][CDS_region][1]
				local_loc_list.append(cds_start)
				local_loc_list.append(cds_stop)
				ref_cds += seg_seq[cds_start-1:cds_stop]
				if CDS_region == 0:
					cds_loc_string+= '\t'+str(cds_start)+'-'+str(cds_stop)
				elif CDS_region > 0:
					cds_loc_string+= ','+str(cds_start)+'-'+str(cds_stop)
				for l in range(cds_start-1,cds_stop):
					all_cds_site_loc_dict[acc_seq][CDS][l] = cur_cds_loc
					cur_cds_loc += 1
			cds_loc_string+= '\n'
			all_cds_refs_dict[acc_seq][CDS] = ref_cds
		cds_first_loc_dict[acc_seq] = min(local_loc_list)
cds_parse_outfile = open(project_dir+"CDS_locs."+output_suffix+".txt","w")
cds_parse_outfile.write(cds_loc_string)
cds_parse_outfile.close()


##############################################################################################################################
#### Find iSNVs
# {'cov':0,'prop':1,'qual':2,'map_qual':3,'read_len':4,'read_loc':5,'read_R_loc':6,'read_prop_loc':7,'mismatch':8,'indel':9}

info_tup_dict = {'cov':0,'prop':1,'qual':2,'map_qual':3,'read_len':4,'read_loc':5,'read_R_loc':6,'read_prop_loc':7,'mismatch':8,'indel':9}
base_tup_dict = {'A':0,'T':1,'C':2,'G':3}

all_sites_info = open(project_dir+"all_sites.base_info."+output_suffix+".txt","w")
all_sites_info.write("sampleID\tsegment\tloc\tcov\tmajor_nt\tminor_nt\tmajor_prop\tmajor_qual\tmajor_map_qual\tmajor_loc\tmajor_mismatch\tmajor_indel\tminor_prop\tminor_qual\tminor_map_qual\tminor_loc\tminor_mismatch\tminor_indel\t\n")
# all_sites_info.write("sampleID\tsegment\tloc\tcov\tmajor_nt\tminor_nt\tmajor_prop\tmajor_qual\tmajor_map_qual\tmajor_read_len\tmajor_loc\tmajor_Rloc\tmajor_read_prop_loc\tmajor_mismatch\tmajor_indel\tminor_prop\tminor_qual\tminor_map_qual\tminor_read_len\tminor_loc\tminor_Rloc\tminor_read_prop_loc\tminor_mismatch\tminor_indel\tmed_read_loc\tminor_qual_pass\tmajor_qual_pass\n")
all_sites_info.close()
iSNV_info = open(project_dir+"iSNV_info."+output_suffix+".txt","w")
iSNV_info.write("sampleID\tsegment\tloc\tcov\tmajor_nt\tminor_nt\tmajor_prop\tmajor_qual\tmajor_map_qual\tmajor_read_len\tmajor_base_loc\tmajor_base_R_loc\tmajor_prop_loc\tmajor_mismatches\tmajor_indels\tminor_prop\tminor_qual\tminor_map_qual\tminor_read_len\tminor_base_loc\tminor_base_R_loc\tminor_prop_loc\tminor_mismatches\tminor_indels\n")
iSNV_info.close()

def calc_prop_stat(prop, background_prop_list):
	if prop == 0 or len(background_prop_list) < 3:
		prop_stat = 0
	else:
		background_mean = np.mean(background_prop_list)
		background_stdv = np.std(background_prop_list)
		abs_diff = np.absolute(prop - background_mean)
		if background_stdv == 0:
			prop_stat = 5
		elif abs_diff == 0:
			prop_stat = 0
		else:
			prop_stat = abs_diff/background_stdv
	return prop_stat


def identify_putative_iSNVs(sampleID,segment_list,allele_info_dict,background_prop_dict,cov_count_dict,min_proportion,min_qual,min_avg_map_qual,min_avg_read_loc,max_avg_read_mismatch,max_avg_read_indel,min_major_map_qual,min_cov,cds_loc_dict,max_length_dict,sample_consensus_seq_dict,end_skip_len,sites_to_exclude,info_tup_dict,base_tup_dict,output_suffix):
	sample_iSNV_dict = {}
	avg_cov_dict = {}
	qual_pass_dict = {}
	iSNV_info_string = ''
	all_sites_info_string = ''
	for s in range(0,len(segment_list)):
		segment = segment_list[s]
		sample_root = segment.split("_run")[0]
		background_freq_list = background_prop_dict[segment]
		avg_cov_dict[segment] = []
		acc_seg = sampleID+"_"+segment
		try:
			maxlen = max_length_dict[segment]
			cov_pass_prop = cov_count_dict[segment]/maxlen
		except:
			maxlen = 1
			cov_pass_prop = 0
		try:
			acc_seg_seq = sample_consensus_seq_dict[acc_seg]
		except:
			acc_seg_seq = "-"*maxlen

		for loc in range(end_skip_len,maxlen-1-end_skip_len):
			try:
				info_tup = allele_info_dict[segment][loc]
			except:
				info_tup = ()
			if info_tup != () and cov_pass_prop >= min_cov_pass_site_prop:
				cov = info_tup[0]
				avg_cov_dict[segment].append(cov)
				if cov >= min_cov:
					A_prop,T_prop,C_prop,G_prop = info_tup[1][0],info_tup[1][1],info_tup[1][2],info_tup[1][3]
					A_qual,T_qual,C_qual,G_qual = info_tup[2][0],info_tup[2][1],info_tup[2][2],info_tup[2][3]
					A_map_qual,T_map_qual,C_map_qual,G_map_qual = info_tup[3][0],info_tup[3][1],info_tup[3][2],info_tup[3][3]
					A_read_len,T_read_len,C_read_len,G_read_len = info_tup[4][0],info_tup[4][1],info_tup[4][2],info_tup[4][3]
					A_read_loc,T_read_loc,C_read_loc,G_read_loc = info_tup[5][0],info_tup[5][1],info_tup[5][2],info_tup[5][3]
					A_read_R_loc,T_read_R_loc,C_read_R_loc,G_read_R_loc = info_tup[6][0],info_tup[6][1],info_tup[6][2],info_tup[6][3]
					A_read_prop_loc,T_read_prop_loc,C_read_prop_loc,G_read_prop_loc = info_tup[7][0],info_tup[7][1],info_tup[7][2],info_tup[7][3]
					A_read_mismatch,T_read_mismatch,C_read_mismatch,G_read_mismatch = info_tup[8][0],info_tup[8][1],info_tup[8][2],info_tup[8][3]
					A_read_indel,T_read_indel,C_read_indel,G_read_indel = info_tup[9][0],info_tup[9][1],info_tup[9][2],info_tup[9][3]
					
					prop_nt_list = sorted([(A_prop,"A"),(T_prop,"T"),(C_prop,"C"),(G_prop,"G")],reverse=True)
					
					major_nt = prop_nt_list[0][1]
					minor_nt = prop_nt_list[1][1]
					minor_prop = info_tup[1][base_tup_dict[minor_nt]]
					minor_qual = info_tup[2][base_tup_dict[minor_nt]]
					minor_mapq = info_tup[3][base_tup_dict[minor_nt]]
					minor_readloc = info_tup[5][base_tup_dict[minor_nt]]
					minor_mismatch = info_tup[8][base_tup_dict[minor_nt]]
					minor_indel = info_tup[8][base_tup_dict[minor_nt]]

					major_mapq = info_tup[3][base_tup_dict[major_nt]]
					major_readloc = info_tup[5][base_tup_dict[major_nt]]

					major_qual_pass = False
					if major_mapq >= min_major_map_qual and major_readloc >= min_avg_read_loc:
						major_qual_pass = True
						try:
							qual_pass_dict[segment] += 1
						except:
							qual_pass_dict[segment] = 1
					
					if minor_prop >= 0.001 and minor_qual >= 25.0 and minor_mapq >= 25.0 and minor_mismatch <= 10 and minor_indel <= 5:
						all_sites_info_string += sampleID+"\t"+segment+"\t"+str(loc)+"\t"+str(int(cov))+"\t"+major_nt+"\t"+minor_nt
						info_to_write = [1,2,3,5,8,9]
						for v in range(0,len(info_to_write)):
							val = info_to_write[v]
							major_val = info_tup[val][base_tup_dict[major_nt]]
							all_sites_info_string += "\t"+str(major_val)
						for v in range(0,len(info_to_write)):
							val = info_to_write[v]
							minor_val = info_tup[val][base_tup_dict[minor_nt]]
							all_sites_info_string += "\t"+str(minor_val)
						all_sites_info_string += "\n"

					if prop_nt_list[1][0] >= min_proportion:
						minor_subID = segment+"-"+str(loc)+minor_nt
						major_subID = segment+"-"+str(loc)+major_nt
						
						iSNV_info_string += sampleID+"\t"+segment+"\t"+str(loc)+"\t"+str(int(cov))+"\t"+major_nt+"\t"+minor_nt
						for val in range(1,10):
							major_val = info_tup[val][base_tup_dict[major_nt]]
							iSNV_info_string += "\t"+str(major_val)
						for val in range(1,10):
							minor_val = info_tup[val][base_tup_dict[minor_nt]]
							iSNV_info_string += "\t"+str(minor_val)
						try:
							med_read_loc = np.median(major_read_loc_dict[segment][loc])
						except:
							med_read_loc = 'na'
						iSNV_info_string += "\n"

						nt_list = ['A','T','C','G']
						for nt_num in range(0,4):
							local_nt = nt_list[nt_num]
							local_sub_ID = segment+"-"+str(loc)+local_nt
							local_prop = allele_info_dict[segment][loc][info_tup_dict['prop']][nt_num]
							local_qual = allele_info_dict[segment][loc][info_tup_dict['qual']][nt_num]
							local_map_qual = allele_info_dict[segment][loc][info_tup_dict['map_qual']][nt_num]
							local_read_loc = allele_info_dict[segment][loc][info_tup_dict['read_loc']][nt_num]
							local_read_mismatch = allele_info_dict[segment][loc][info_tup_dict['mismatch']][nt_num]
							local_read_indel = allele_info_dict[segment][loc][info_tup_dict['indel']][nt_num]
							local_nt_count = round(local_prop*cov,0)

							local_qual_pass = False
							if local_qual >= min_qual and local_map_qual >= min_avg_map_qual and local_read_loc >= min_avg_read_loc and local_read_mismatch <= max_avg_read_mismatch and local_read_indel <= max_avg_read_indel and local_nt_count >= min_observations:
								local_qual_pass = True

							if local_qual_pass == True and local_prop >= min_proportion and major_qual_pass == True:
								local_background_freq_list = background_freq_list[max(0,loc-int(sigma_flank_len/2)):loc]
								local_background_freq_list.extend(background_freq_list[loc+1:min(len(background_freq_list),loc+int(sigma_flank_len/2))])
								local_prop_pval = calc_prop_stat(-1*np.log10(local_prop), local_background_freq_list)
								if local_prop_pval >= min_freq_sigma_score:
									local_allele_tup = (local_prop,local_qual,local_nt,local_prop_pval)
									try:
										sample_iSNV_dict[segment][loc].append(local_allele_tup)
									except:
										try:
											sample_iSNV_dict[segment][loc] = [local_allele_tup]
										except:
											sample_iSNV_dict[segment] = {}
											sample_iSNV_dict[segment][loc] = [local_allele_tup]

		all_sites_info = open(project_dir+"all_sites.base_info."+output_suffix+".txt","a")
		all_sites_info.write(all_sites_info_string)
		all_sites_info.close()
		all_sites_info_string = ''
		try:
			if len(avg_cov_dict[segment])>10:
				avg_cov = round(np.average(avg_cov_dict[segment]),2)
			else:
				avg_cov = 0.0
		except:
			avg_cov = 0.0
	iSNV_info = open(project_dir+"iSNV_info."+output_suffix+".txt","a")
	iSNV_info.write(iSNV_info_string)
	iSNV_info.close()

	out_tup = (sampleID,sample_iSNV_dict, avg_cov_dict)
	return out_tup



processed_list = []
if parallel_process == True:
	processed_list = Parallel(n_jobs=num_cores)(delayed(identify_putative_iSNVs)(sampleID,segment_list,all_allele_info_dict[sampleID],background_freq_dict[sampleID],cov_pass_dict[sampleID],min_proportion,min_qual,min_avg_map_qual,min_avg_read_loc,max_avg_read_mismatch,max_avg_read_indel,min_major_map_qual,min_cov,all_cds_refs_dict,max_length_dict,sample_consensus_seq_dict,end_skip_len,sites_to_exclude,info_tup_dict,base_tup_dict,output_suffix) for sampleID in sampleID_list)
elif parallel_process == False:
	for sampleID in sampleID_list:
		output = identify_putative_iSNVs(sampleID,segment_list,all_allele_info_dict[sampleID],background_freq_dict[sampleID],cov_pass_dict[sampleID],min_proportion,min_qual,min_avg_map_qual,min_avg_read_loc,max_avg_read_mismatch,max_avg_read_indel,min_major_map_qual,min_cov,all_cds_refs_dict,max_length_dict,sample_consensus_seq_dict,end_skip_len,sites_to_exclude,info_tup_dict,base_tup_dict,output_suffix)
		processed_list.append(output)
print("Finished identifying putative iSNVs")

iSNV_dict = {}
avg_cov_dict = {}
for tup in processed_list:
	sampleID = tup[0]
	iSNV_dict[sampleID] = tup[1]
	avg_cov_dict[sampleID] = tup[2]
del processed_list

outfile = open(project_dir+"avg_cov_table."+output_suffix+".txt","w")
for s in range(0,len(ref_list)):
	segment = ref_list[s].split(".")[0]
	outfile.write("\t"+segment)
outfile.write("\n")
for a in range(0,len(sampleID_list)):
	sampleID = sampleID_list[a]
	outfile.write(sampleID)
	for s in range(0,len(ref_list)):
		segment = ref_list[s].split(".")[0]
		try:
			if len(avg_cov_dict[sampleID][segment])>10:
				avg_cov = round(np.average(avg_cov_dict[sampleID][segment]),2)
			else:
				avg_cov = 0
		except:
			avg_cov = 0.0
		avg_cov_dict[sampleID][segment] = avg_cov
		outfile.write("\t"+str(avg_cov))
	outfile.write("\n")
outfile.close()
outfile = open(project_dir+"cov_pass_prop_table."+output_suffix+".txt","w")
for s in range(0,len(ref_list)):
	segment = ref_list[s].split(".")[0]
	outfile.write("\t"+segment)
outfile.write("\n")
for a in range(0,len(sampleID_list)):
	sampleID = sampleID_list[a]
	outfile.write(sampleID)
	for s in range(0,len(ref_list)):
		segment = ref_list[s].split(".")[0]
		try:
			maxlen = max_length_dict[segment]
		except:
			maxlen = 1
		try:
			cov_pass_count = float(cov_pass_dict[sampleID][segment])
		except:
			cov_pass_count = 0
		cov_pass_prop = round_to_n_sig_figs(cov_pass_count/maxlen,2)
		outfile.write("\t"+str(cov_pass_prop))
		cov_pass_dict[sampleID][segment] = float(cov_pass_prop)
	outfile.write("\n")
outfile.close()


num_iSNVs_found = 0
iSNV_count_dict = {}
iSNV_outlines = ''
for sample_root in unique_sample_dict:
	iSNV_count_dict[sample_root] = 0
	try:
		replicate_sample_list = unique_sample_dict[sample_root]
	except:
		replicate_sample_list = []
	if len(replicate_sample_list) >= 2:
		replicate_sample_list = sorted(replicate_sample_list,reverse=False)
		sampleID_1 = replicate_sample_list[0]
		sampleID_2 = replicate_sample_list[1]
		s1_segment_list = list(iSNV_dict[sampleID_1].keys())
		s2_segment_list = list(iSNV_dict[sampleID_2].keys())
		local_segment_list = []
		for s in range(0,len(s1_segment_list)):
			segment = s1_segment_list[s]
			try:
				iSNV_dict[sampleID_2][segment]
				local_segment_list.append(segment)
			except:
				pass
		if local_segment_list != []:
			for s in range(0,len(local_segment_list)):
				segment = local_segment_list[s]
				s1_iSNV_dict = iSNV_dict[sampleID_1][segment]
				s2_iSNV_dict = iSNV_dict[sampleID_2][segment]
				
				s1_link_locs = {}
				for s1_loc in s1_iSNV_dict:
					try:
						link_loc = link_loc_dict[sample_root][sampleID_1+'_'+segment][s1_loc]
						s1_link_locs[link_loc] = s1_loc
					except:
						pass
				s2_link_locs = {}
				shared_link_loc_dict = {}
				for s2_loc in s2_iSNV_dict:
					try:
						link_loc = link_loc_dict[sample_root][sampleID_2+'_'+segment][s2_loc]
						s2_link_locs[link_loc] = s2_loc
						try:
							s1_loc = s1_link_locs[link_loc]
							shared_link_loc_dict[link_loc] = (s1_loc,s2_loc)
						except:
							pass
					except:
						pass
				if shared_link_loc_dict != {}:
					s1_ref_seq = sample_consensus_seq_dict[sampleID_1+'_'+segment]
					s2_ref_seq = sample_consensus_seq_dict[sampleID_2+'_'+segment]
					for link_loc in shared_link_loc_dict:
						include_site = False
						s1_cov = all_allele_info_dict[sampleID_1][segment][s1_loc][0]
						s2_cov = all_allele_info_dict[sampleID_2][segment][s2_loc][0]
						total_cov = s1_cov + s2_cov
						
						loc_tup = shared_link_loc_dict[link_loc]
						s1_loc = loc_tup[0]
						s2_loc = loc_tup[1]
						s1_allele_list = sorted(s1_iSNV_dict[s1_loc],reverse=True) #local_allele_tup = (local_prop,local_qual,local_nt,local_prop_pval)
						s2_allele_list = sorted(s2_iSNV_dict[s2_loc],reverse=True) #local_allele_tup = (local_prop,local_qual,local_nt,local_prop_pval)

						s1_ref_nt = s1_ref_seq[s1_loc]
						s2_ref_nt = s2_ref_seq[s2_loc]

						if len(s1_allele_list) == 2 and len(s2_allele_list) == 2:

							s1_major_nt = s1_allele_list[0][2]
							s1_major_qual = s1_allele_list[0][1]
							s1_major_prop = s1_allele_list[0][0]
							s1_major_sigma = s1_allele_list[0][3]
							s1_minor_nt = s1_allele_list[1][2]
							s1_minor_qual = s1_allele_list[1][1]
							s1_minor_prop = s1_allele_list[1][0]
							s1_minor_sigma = s1_allele_list[1][3]

							s2_major_nt = s2_allele_list[0][2]
							s2_major_qual = s2_allele_list[0][1]
							s2_major_prop = s2_allele_list[0][0]
							s2_major_sigma = s2_allele_list[0][3]
							s2_minor_nt = s2_allele_list[1][2]
							s2_minor_qual = s2_allele_list[1][1]
							s2_minor_prop = s2_allele_list[1][0]
							s2_minor_sigma = s2_allele_list[1][3]

							if s1_minor_nt == s2_minor_nt and s1_major_nt == s2_major_nt:
								include_site = True
								minor_avg_prop = ((s1_minor_prop*s1_cov)+(s2_minor_prop*s2_cov))/total_cov
								major_avg_prop = ((s1_major_prop*s1_cov)+(s2_major_prop*s2_cov))/total_cov
								root_minor_nt = s1_minor_nt
								root_major_nt = s1_major_nt
							elif s1_minor_nt == s2_major_nt and s1_major_nt == s2_minor_nt:
								temp_minor_avg_prop = ((s1_minor_prop*s1_cov)+(s2_major_prop*s2_cov))/total_cov
								temp_major_avg_prop = ((s1_major_prop*s1_cov)+(s2_minor_prop*s2_cov))/total_cov
								if temp_minor_avg_prop < temp_major_avg_prop:
									minor_avg_prop = temp_minor_avg_prop
									major_avg_prop = temp_major_avg_prop
									root_minor_nt = s1_minor_nt
									root_major_nt = s1_major_nt
								else:
									minor_avg_prop = temp_major_avg_prop
									major_avg_prop = temp_minor_avg_prop
									root_minor_nt = s1_major_nt
									root_major_nt = s1_minor_nt
								include_site = True

						if include_site == True:
							num_iSNVs_found += 1
							s1_root_minor_prop = all_allele_info_dict[sampleID_1][segment][s1_loc][1][base_tup_dict[root_minor_nt]]
							s1_root_major_prop = all_allele_info_dict[sampleID_1][segment][s1_loc][1][base_tup_dict[root_major_nt]]
							s2_root_minor_prop = all_allele_info_dict[sampleID_2][segment][s2_loc][1][base_tup_dict[root_minor_nt]]
							s2_root_major_prop = all_allele_info_dict[sampleID_2][segment][s2_loc][1][base_tup_dict[root_major_nt]]
							
							try:
								iSNV_count_dict[sample_root] +=1
							except:
								iSNV_count_dict[sample_root] = 1

							s1_covpass = cov_pass_dict[sampleID_1][segment]
							s2_covpass = cov_pass_dict[sampleID_2][segment]
							s1_avgcov = avg_cov_dict[sampleID_1][segment]
							s2_avgcov = avg_cov_dict[sampleID_2][segment]
							
							s1_seq_qual = s1_covpass/s2_covpass + np.log10(s1_avgcov)/np.log10(s2_avgcov) + np.log10(s1_cov)/np.log10(s2_cov)
							s2_seq_qual = s2_covpass/s1_covpass + np.log10(s2_avgcov)/np.log10(s1_avgcov) + np.log10(s2_cov)/np.log10(s1_cov)
							if s1_seq_qual > s2_seq_qual:
								ref_sample = sampleID_1
								ref_loc = s1_loc
							elif s1_seq_qual < s2_seq_qual:
								ref_sample = sampleID_2
								ref_loc = s2_loc
							else:
								ref_sample = sampleID_1
								ref_loc = s1_loc

							sub_ID = segment+"-"+str(link_loc)+root_minor_nt
							sub_label = segment+"-"+root_major_nt+str(link_loc)+root_minor_nt
							ref_sub_label = segment+"-"+root_major_nt+str(link_loc)+root_minor_nt

							
							cds_sub_list = []
							ref_sample_segment = ref_sample+"_"+segment
							try:
								cds_ref_dict = all_cds_refs_dict[ref_sample_segment]
							except:
								cds_ref_dict = {}
							if cds_ref_dict != {}:
								for CDS in all_cds_refs_dict[ref_sample_segment]:
									ref_cds_seq = all_cds_refs_dict[ref_sample_segment][CDS]
									try:
										cds_loc = all_cds_site_loc_dict[ref_sample_segment][CDS][ref_loc]
									except:
										cds_loc = ''
									if cds_loc != '':
										ref_cds_nt = ref_cds_seq[cds_loc]
										cds_start = 0
										codon_number,codon_site = codon_finder(cds_loc, cds_start)
										ref_cds_codon = ref_cds_seq[(codon_number*3+cds_start):(codon_number*3+cds_start)+3]
										sub_cds_codon = list(ref_cds_codon)
										if root_minor_nt != ref_cds_nt and root_major_nt == ref_cds_nt:
											sub_cds_codon[codon_site] = root_minor_nt

										elif root_minor_nt == ref_cds_nt and root_major_nt != ref_cds_nt:
											sub_cds_codon[codon_site] = root_major_nt
										else:
											sub_cds_codon[codon_site] = root_minor_nt
										sub_cds_codon = "".join(sub_cds_codon)
										ref_cds_aa = codon_to_aa_dict[ref_cds_codon]
										sub_cds_aa = codon_to_aa_dict[sub_cds_codon].replace("*","X")
										if ref_cds_aa == sub_cds_aa:
											cds_sub_type = "S"
										elif ref_cds_aa != sub_cds_aa:
											cds_sub_type = "NS"
										sub_info_tup = (CDS,codon_number,ref_cds_aa,sub_cds_aa,cds_sub_type)
										cds_sub_list.append(sub_info_tup)
							subtype_string = ''
							if len(cds_sub_list) == 0:
								subtype_string = "na\tna"
							else:
								for t in range(0,len(cds_sub_list)):
									sub_tup = cds_sub_list[t]
									if t>0:
										subtype_string += ", "
									if sub_tup[4] == "NS":
										subtype_string += sub_tup[0]+"-"+sub_tup[2]+str(sub_tup[1])+sub_tup[3]+" ("+sub_tup[4]+")"
									elif sub_tup[4] == "S":
										subtype_string += sub_tup[0]+"-"+str(sub_tup[1])+sub_tup[3]+" ("+sub_tup[4]+")"
								subtype_string += "\t"
								for t in range(0,len(cds_sub_list)):#
									sub_tup = cds_sub_list[t]
									if t>0:
										subtype_string += ","
									subtype_string += sub_tup[2]
								subtype_string += "\t"
								for t in range(0,len(cds_sub_list)):
									sub_tup = cds_sub_list[t]
									if t>0:
										subtype_string += ","
									subtype_string += sub_tup[3]

							iSNV_outlines += sample_root +'\t'+ segment +'\t'+ str(link_loc) +'\t'+ str(total_cov) +'\t'+ root_minor_nt +'\t'+ root_major_nt +'\t'+ str(minor_avg_prop) +'\t'+ str(major_avg_prop) +'\t'+ str(s1_loc) +'\t'+ str(s1_cov) +'\t'+ str(s1_root_minor_prop) +'\t'+ str(s1_root_major_prop) +'\t'+ str(s1_minor_sigma) +'\t'+ str(s1_major_sigma) +'\t'+ str(s2_loc) +'\t'+ str(s2_cov) +'\t'+ str(s2_root_minor_prop) +'\t'+ str(s2_root_major_prop) +'\t'+ str(s2_minor_sigma) +'\t'+ str(s2_major_sigma) +'\t'+subtype_string+'\n'
						

isnv_outfile = open(project_dir+"iSNV_replicates."+output_suffix+".txt","w")
isnv_outfile.write(iSNV_outlines)
isnv_outfile.close()

print(str(num_iSNVs_found)+ " iSNVs detected in both replicates")


if skip_writing_all_site_summaries == False:
	nt_list = ['A','T','C','G']
	num_allele = 1
	for num_allele in range(0,2):
		cov_lines = ''
		prop_lines = ''
		readloc_lines = ''
		qual_lines = ''
		mapqual_lines = ''
		mismatch_lines = ''
		indel_lines = ''
		if num_allele == 1:
			cov_outfile = open(project_dir+"all_cov."+output_suffix+".txt","w")
		prop_outfile = open(project_dir+"all_sites.prop."+output_suffix+"."+str(num_allele)+".txt","w")
		readloc_outfile = open(project_dir+"all_sites.read_loc."+output_suffix+"."+str(num_allele)+".txt","w")
		qual_outfile = open(project_dir+"all_sites.qual."+output_suffix+"."+str(num_allele)+".txt","w")
		mapqual_outfile = open(project_dir+"all_sites.map_qual."+output_suffix+"."+str(num_allele)+".txt","w")
		mismatch_outfile = open(project_dir+"all_sites.mismatch."+output_suffix+"."+str(num_allele)+".txt","w")
		indel_outfile = open(project_dir+"all_sites.indel."+output_suffix+"."+str(num_allele)+".txt","w")
		if num_allele == 1:
			cov_lines += "segment\tloc"
		prop_lines += "segment\tloc"
		readloc_lines += "segment\tloc"
		qual_lines += "segment\tloc"
		mapqual_lines += "segment\tloc"
		mismatch_lines += "segment\tloc"
		indel_lines += "segment\tloc"


		for a in range(0,len(sampleID_list)):
			sampleID = sampleID_list[a]
			if num_allele == 1:
				cov_lines += "\t"+sampleID
			prop_lines += "\t"+sampleID
			readloc_lines += "\t"+sampleID
			qual_lines += "\t"+sampleID
			mapqual_lines += "\t"+sampleID
			mismatch_lines += "\t"+sampleID
			indel_lines += "\t"+sampleID
		if num_allele == 1:
			cov_lines += "\n"
		prop_lines += "\n"
		readloc_lines += "\n"
		qual_lines += "\n"
		mapqual_lines += "\n"
		mismatch_lines += "\n"
		indel_lines += "\n"
		for seg_ref in ref_set_list:
			seg = seg_ref.split("-")[0]
			segment = seg.split(".")[0]
			try:
				maxlen = max_length_dict[segment]
			except:
				maxlen = 1
			for loc in range(0,maxlen):
				if num_allele == 1:
					cov_lines += segment+"\t"+str(loc)
				prop_lines += segment+"\t"+str(loc)
				readloc_lines += segment+"\t"+str(loc)
				qual_lines += segment+"\t"+str(loc)
				mapqual_lines += segment+"\t"+str(loc)
				mismatch_lines += segment+"\t"+str(loc)
				indel_lines += segment+"\t"+str(loc)
				for a in range(0,len(sampleID_list)):
					sampleID = sampleID_list[a]
					cov_thresh = min_summary_cov

					nt_prop_list = []
					for nt_num in range(0,4):
						try:
							nt_prop = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['prop']][nt_num]
							nt_prop_list.append((nt_prop,nt_list[nt_num]))
						except:
							pass

					cov = 0
					allele_freq = ''
					read_loc = ''
					qual = ''
					map_qual = ''
					mismat = ''
					indel = ''
					if nt_prop_list != []:
						nt_prop_list = sorted(nt_prop_list,reverse=True)
						major_nt = nt_prop_list[0][1]
						minor_nt = nt_prop_list[num_allele][1]
						try:
							cov = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['cov']]
							if cov >= cov_thresh:
								allele_freq = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['prop']][base_tup_dict[minor_nt]]
								read_loc = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['read_loc']][base_tup_dict[minor_nt]]
								qual = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['qual']][base_tup_dict[minor_nt]]
								map_qual = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['map_qual']][base_tup_dict[minor_nt]]
								mismat = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['mismatch']][base_tup_dict[minor_nt]]
								indel = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['indel']][base_tup_dict[minor_nt]]
							else:
								allele_freq = ''
								read_loc = ''
								qual = ''
								map_qual = ''
								mismat = ''
								indel = ''
						except:
							pass
					if num_allele == 1:
						cov_lines += "\t"+str(cov)
					prop_lines += "\t"+str(allele_freq)
					readloc_lines += "\t"+str(read_loc)
					qual_lines += "\t"+str(qual)
					mapqual_lines += "\t"+str(map_qual)
					mismatch_lines += "\t"+str(mismat)
					indel_lines += "\t"+str(indel)
				if num_allele == 1:
					# cov_outfile.write("\n")
					cov_lines += "\n"
				prop_lines += "\n"
				readloc_lines += "\n"
				qual_lines += "\n"
				mapqual_lines += "\n"
				mismatch_lines += "\n"
				indel_lines += "\n"
				if loc%200==0:
					if num_allele == 1:
						cov_outfile.write(cov_lines)
					prop_outfile.write(prop_lines)
					readloc_outfile.write(readloc_lines)
					qual_outfile.write(qual_lines)
					mapqual_outfile.write(mapqual_lines)
					mismatch_outfile.write(mismatch_lines)
					indel_outfile.write(indel_lines)
					cov_lines = ''
					prop_lines = ''
					readloc_lines = ''
					qual_lines = ''
					mapqual_lines = ''
					mismatch_lines = ''
					indel_lines = ''
			if num_allele == 1:
				cov_outfile.write(cov_lines)
			prop_outfile.write(prop_lines)
			readloc_outfile.write(readloc_lines)
			qual_outfile.write(qual_lines)
			mapqual_outfile.write(mapqual_lines)
			mismatch_outfile.write(mismatch_lines)
			indel_outfile.write(indel_lines)
			cov_lines = ''
			prop_lines = ''
			readloc_lines = ''
			qual_lines = ''
			mapqual_lines = ''
			mismatch_lines = ''
			indel_lines = ''
		if num_allele == 1:
			cov_outfile.close()
		prop_outfile.close()
		readloc_outfile.close()
		qual_outfile.close()
		mapqual_outfile.close()
		mismatch_outfile.close()
		indel_outfile.close()

print("Complete.")
