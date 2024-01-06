import os
import sys
import numpy as np
import scipy.stats

parallel_process = True
parallel_max_cpu = -5

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

output_suffix = "IAV_WGS"

annotation_file_path = ref_dir + "annotation.tbl"


exclude_samples = {}
try:
	filename = "samples_to_exclude.txt"
	infile = open(project_dir+filename,"r")
	for line in infile:
		line = line.strip()
		exclude_samples[line] = ''
	infile.close
except:
	pass

skip_translate = True
skip_writing_all_site_summaries = True

ref_set_list = ['A_HA_H1.fa','A_HA_H3.fa','A_NA_N1.fa','A_NA_N2.fa','A_M.fa','A_NP.fa','A_NS.fa','A_PA.fa','A_PB1.fa','A_PB2.fa','B_HA.fa','B_M.fa','B_NA.fa','B_NP.fa','B_NS.fa','B_PA.fa','B_PB1.fa','B_PB2.fa']


ref_list = ref_set_list
segment_list = []
for ref_set in ref_set_list:
	refs = ref_set.split("-")
	segment = ref_set.split("-")[0].split(".f")[0]
	segment_list.append(segment)
segment_list = sorted(segment_list)

min_proportion = 0.03
min_qual = 35.0
min_avg_map_qual = 40.0
min_avg_read_loc = 30.0
max_avg_read_mismatch = 2.0
max_avg_read_indel = 2.0


min_cov = 500
end_skip_len = 20

neighbour_step_size = 2
max_neighbour_SNP_count = 9999


##################################################################################################

min_major_map_qual = 40 #not currently in use
max_major_read_mismatch = 3.0 #not currently in use
max_major_read_indel = 3.0 #not currently in use

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

samples_to_exclude = []
sites_to_exclude = {"segment-1234":''}


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


############################################ FUNCTIONS ############################################
def run_command(command):
	return_code = os.system(command)
	# return_code = os.system(command+" >/dev/null 2>&1")
	if return_code == 0:
		return command
	else:
		sys.exit("Command returned with code: "+str(return_code)+"\n"+command)


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
				contig = line.split(">Feature ")[1]
				found_CDS = False
				cur_cds_locs = []
				cur_product_ID = ''
				cur_gene_ID = ''
			elif found_CDS == False:
				if is_number(line[0]) == True:
					split = line.split("\t")
					product_type = split[2]
					if product_type == "CDS":
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
						start = int(split[0])
						stop = int(split[1])
						cur_cds_locs.append((start,stop))
					else:
						found_CDS = False
				elif line[0] == "\t":
					if line[2] == "\t":
						split = line.split("\t")
						if split[3] == "protein_id":
							cur_product_ID = split[4]
						elif split[3] == "gene":
							cur_gene_ID = split[4]
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
		# print(codon_bases)
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
		if "-" not in codon_bases and "N" not in codon_bases and len(codon_bases)==3:
			aa = codon_to_aa_dict[codon_bases]
		else:
			aa = 'X'
		aa_seq_out += aa
	return aa_seq_out


def pull_consensus_seq(sampleID,segment_list,map_dir,max_length_dict,min_cov=5,min_prop=0.5,min_qual=15.0):
	out_dict = {}
	for s in range(0,len(segment_list)):
		segment = segment_list[s]
		base_dict = {}
		substitutions_filename = map_dir+"substitutions/"+segment+"/"+sampleID+'.'+segment+".substitions.txt"
		if os.path.isfile(substitutions_filename) == True:
			sub_infile = open(substitutions_filename,"r")
			for sub_line in sub_infile:
				if len(sub_line) > 0:
					if sub_line[0] != "#" and sub_line[0] != "loc":
						sub_line = sub_line.strip().split("\t")
						loc = int(sub_line[0])
						cov = float(sub_line[1])
						A_prop,T_prop,C_prop,G_prop = float(sub_line[2]),float(sub_line[3]),float(sub_line[4]),float(sub_line[5])
						A_qual,T_qual,C_qual,G_qual = float(sub_line[6]),float(sub_line[7]),float(sub_line[8]),float(sub_line[9])
						nt_list = sorted([(A_prop,A_qual,"A"),(T_prop,T_qual,"T"),(C_prop,C_qual,"C"),(G_prop,G_qual,"G")],reverse=True)
						major_allele_prop = nt_list[0][0]
						major_allele_qual = nt_list[0][1]
						major_allele_base = nt_list[0][2]
						if major_allele_prop >= min_prop and major_allele_qual >= min_qual and cov >= min_cov:
							base_dict[loc] = major_allele_base
						else:
							base_dict[loc] = "N"
		len_seq = max_length_dict[segment]
		out_seq = ''
		for loc in range(0,len_seq):
			try:
				base = base_dict[loc]
			except:
				base = "-"
			out_seq += base
		del base_dict
		out_dict[segment] = out_seq
	return (sampleID,out_dict)

def read_in_allele_info(sampleID,segment_list,map_dir):
	out_dict = {}
	for s in range(0,len(segment_list)):
		segment = segment_list[s]
		substitutions_filename = map_dir+"substitutions/"+segment+"/"+sampleID+'.'+segment+".substitions.txt"
		if os.path.isfile(substitutions_filename) == True:
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
							out_dict[segment][loc] = info_tup
						except:
							out_dict[segment] = {}
							out_dict[segment][loc] = info_tup
	return (sampleID,out_dict)

def count_neighbouring_SNPs(sampleID,allele_info_dict,segment_list,max_length_dict,step_size,min_proportion,min_cov,max_neighbour_SNP_count):
	minor_prop_loc_dict = {}
	for s in range(0,len(segment_list)):
		segment = segment_list[s]
		minor_prop_loc_dict[segment] = {}
		for loc in range(0,max_length_dict[segment]):
			try:
				info_tup = allele_info_dict[segment][loc]
			except:
				info_tup = ()
			if info_tup != ():
				cov = info_tup[0]
				if cov >= min_cov:
					A_prop,T_prop,C_prop,G_prop = info_tup[1][0],info_tup[1][1],info_tup[1][2],info_tup[1][3]
					prop_nt_list = sorted([(A_prop,"A"),(T_prop,"T"),(C_prop,"C"),(G_prop,"G")],reverse=True)
					minor_prop = prop_nt_list[1][0]
					if minor_prop >= min_proportion:
						minor_prop_loc_dict[segment][loc] = 0
		for loc in minor_prop_loc_dict[segment]:
			for l in range(max(0,loc-step_size),min(max_length_dict[segment],loc+step_size)):
				try:
					c = minor_prop_loc_dict[segment][l]
					minor_prop_loc_dict[segment][loc] += 1
				except:
					pass
		num_over_max = 0
		for loc in minor_prop_loc_dict[segment]:
			if minor_prop_loc_dict[segment][loc] > max_neighbour_SNP_count:
				num_over_max += 1
		minor_prop_loc_dict[segment]["over"] = num_over_max
	return (sampleID,minor_prop_loc_dict)


############################################   CORE   ############################################
### Load the reference sequences for each segment 
consensus_seq_dict = {}
max_length_dict = {}
infile = open(ref_dir+"all_refs.fa","r")
for line in infile:
	line = line.strip()
	if line[0] == ">":
		header = line[1:len(line)]
		segment = header
	else:
		try:
			consensus_seq_dict[segment] += line
		except:
			consensus_seq_dict[segment] = line
for segment in consensus_seq_dict:
	max_length_dict[segment] = len(consensus_seq_dict[segment])

all_cds_loc_dict = pull_cds_from_annot_file(annotation_file_path)
all_cds_site_loc_dict = {}
all_cds_refs_dict = {}
for segment in all_cds_loc_dict:
	seg_seq = consensus_seq_dict[segment]
	all_cds_site_loc_dict[segment] = {}
	all_cds_refs_dict[segment] = {}
	for CDS in all_cds_loc_dict[segment]:
		ref_cds = ''
		all_cds_site_loc_dict[segment][CDS] = {}
		for CDS_region in range(0,len(all_cds_loc_dict[segment][CDS])):
			cur_cds_loc = len(ref_cds)
			cds_start = all_cds_loc_dict[segment][CDS][CDS_region][0]
			cds_stop = all_cds_loc_dict[segment][CDS][CDS_region][1]
			ref_cds += seg_seq[cds_start-1:cds_stop]
			for l in range(cds_start-1,cds_stop):
				all_cds_site_loc_dict[segment][CDS][l] = cur_cds_loc
				cur_cds_loc += 1
		all_cds_refs_dict[segment][CDS] = ref_cds

# load file that stores the location of CDS start and stop within the reference sequences that reads were aligned to
cds_loc_dict = {}
infile = open(ref_dir+"ref_cds_loc.txt","r")
for line in infile:
	line = line.strip().split("\t")
	segment = line[1]
	locs = line[2]
	cds_loc_dict[segment] = locs
infile.close()

sampleID_list = pull_sample_list(temp_dir,input_read_dir)
flagged_sample_list = []
unique_sample_dict = {}
for sampleID in sampleID_list:
	try:
		exclude_samples[sampleID]
	except:
		flagged_sample_list.append(sampleID)
		sample_root = sampleID.split("-run")[0]
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


processed_list = []
if parallel_process == True:
	processed_list = Parallel(n_jobs=num_cores)(delayed(read_in_allele_info)(sampleID,segment_list,map_dir) for sampleID in sampleID_list)
elif parallel_process == False:
	for sampleID in sampleID_list:
		output = read_in_allele_info(sampleID,segment_list,map_dir)
		processed_list.append(output)

all_allele_info_dict = {}
for tup in processed_list:
	sampleID = tup[0]
	all_allele_info_dict[sampleID] = tup[1]
del processed_list

processed_list = []
if parallel_process == True:
	processed_list = Parallel(n_jobs=num_cores)(delayed(pull_consensus_seq)(sampleID,segment_list,map_dir,max_length_dict) for sampleID in sampleID_list)
elif parallel_process == False:
	for sampleID in sampleID_list:
		output = read_in_allele_info(sampleID,segment_list,map_dir)
		processed_list.append(output)

sample_consensus_seq_dict = {}
for tup in processed_list:
	sampleID = tup[0]
	temp_seq_dict = tup[1]
	for segment in temp_seq_dict:
		acc_seg = sampleID+"\t"+segment
		sample_consensus_seq_dict[acc_seg] = temp_seq_dict[segment]
del processed_list



processed_list = []
if parallel_process == True:
	processed_list = Parallel(n_jobs=num_cores)(delayed(count_neighbouring_SNPs)(sampleID,all_allele_info_dict[sampleID],segment_list,max_length_dict,neighbour_step_size,min_proportion,min_cov,max_neighbour_SNP_count) for sampleID in sampleID_list)
elif parallel_process == False:
	for sampleID in sampleID_list:
		output = read_in_allele_info(sampleID,segment_list,map_dir)
		processed_list.append(output)

neighbouring_SNP_dict = {}
for tup in processed_list:
	sampleID = tup[0]
	neighbouring_SNP_dict[sampleID] = tup[1]
del processed_list


##############################################################################################################################
#### Find iSNVs
# {'cov':0,'prop':1,'qual':2,'map_qual':3,'read_len':4,'read_loc':5,'read_R_loc':6,'read_prop_loc':7,'mismatch':8,'indel':9}

info_tup_dict = {'cov':0,'prop':1,'qual':2,'map_qual':3,'read_len':4,'read_loc':5,'read_R_loc':6,'read_prop_loc':7,'mismatch':8,'indel':9}
base_tup_dict = {'A':0,'T':1,'C':2,'G':3}

all_sites_info = open(project_dir+"all_sites.base_info."+output_suffix+".txt","w")
all_sites_info.write("sampleID\tsegment\tloc\tcov\tmajor_nt\tminor_nt\tmajor_prop\tmajor_qual\tmajor_map_qual\tmajor_read_len\tmajor_loc\tmajor_Rloc\tmajor_read_prop_loc\tmajor_mismatch\tmajor_indel\tminor_prop\tminor_qual\tminor_map_qual\tminor_read_len\tminor_loc\tminor_Rloc\tminor_read_prop_loc\tminor_mismatch\tminor_indel\tmed_read_loc\tminor_qual_pass\tmajor_qual_pass\n")
all_sites_info.close()

iSNV_info = open(project_dir+"iSNV_info."+output_suffix+".txt","w")
iSNV_info.write("sampleID\tsegment\tloc\tcov\tmajor_nt\tminor_nt\tmajor_prop\tmajor_qual\tmajor_map_qual\tmajor_read_len\tmajor_base_loc\tmajor_base_R_loc\tmajor_prop_loc\tmajor_mismatches\tmajor_indels\tminor_prop\tminor_qual\tminor_map_qual\tminor_read_len\tminor_base_loc\tminor_base_R_loc\tminor_prop_loc\tminor_mismatches\tminor_indels\n")
iSNV_info.close()

nt_diversity_outfile = open(project_dir+"nt_diversity."+output_suffix+".txt","w")
nt_diversity_outfile.write("sampleID\tsegment\tgroupID\tavg_cov\tseg_site_count\tseg_S_count\tseg_NS_count\tD_sub\tD_sub_S\tD_sub_NS\tpi_sub\tpi_S\tpi_NS\n")
nt_diversity_outfile.close()


def identify_putative_iSNVs(sampleID,segment_list,allele_info_dict,neighbouring_SNP_dict,max_neighbour_SNP_count,min_proportion,min_qual,min_avg_map_qual,min_avg_read_loc,max_avg_read_mismatch,max_avg_read_indel,min_major_map_qual,max_major_read_mismatch,max_major_read_indel,min_cov,cds_loc_dict,max_length_dict,sample_consensus_seq_dict,end_skip_len,sites_to_exclude,info_tup_dict,base_tup_dict,output_suffix):
	sample_iSNV_dict = {}
	avg_cov_dict = {}
	S_NS_dict = {}
	sample_nt_diversity_dict = {}
	sub_codon_site_dict = {}
	codon_aa_dict = {}
	nt_diversity_string = ''
	iSNV_info_string = ''
	for s in range(0,len(segment_list)):
		all_sites_info_string = ''
		segment = segment_list[s]

		sample_iSNV_dict[segment] = {}
		S_NS_dict[segment] = {}
		sub_codon_site_dict[segment] = {}
		avg_cov_dict[segment] = []
		# codon_aa_dict[segment] = {}
		
		acc_seg = sampleID+"\t"+segment
		cds_start = int(cds_loc_dict[segment].split("-")[0])-1
		cds_stop = int(cds_loc_dict[segment].split("-")[1])-1
		segregating_site_count = 0
		sub_seg_site_count = 0
		pass_cov_site_count = 0
		seg_S_count = 0
		seg_NS_count = 0
		poly_S_count = 0
		poly_NS_count = 0
		D_all = 0.0
		D_sub = 0.0
		D_sub_S = 0.0
		D_sub_NS = 0.0

		try:
			acc_seg_seq = sample_consensus_seq_dict[acc_seg]
		except:
			acc_seg_seq = "-"*max_length_dict[segment]
		try:
			ref_seq = consensus_seq_dict[segment]
		except:
			ref_seq = "-"*max_length_dict[segment]

		for loc in range(end_skip_len,max_length_dict[segment]-1):
			try:
				sites_to_exclude[segment+"-"+str(loc)]
			except:
				try:
					info_tup = allele_info_dict[segment][loc]
				except:
					info_tup = ()
			if info_tup != ():
				cov = info_tup[0]
				avg_cov_dict[segment].append(cov)
				if cov >= min_cov:
					pass_cov_site_count += 1
					try:
						neighbour_SNP_count = neighbouring_SNP_dict[segment][loc]
					except:
						neighbour_SNP_count = 0
					try:
						neighbour_over_count = neighbouring_SNP_dict[segment]["over"]
					except:
						neighbour_over_count = 0
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
					
					all_sites_info_string += sampleID+"\t"+segment+"\t"+str(loc)+"\t"+str(int(cov))+"\t"+major_nt+"\t"+minor_nt
					for val in range(1,10):
						major_val = info_tup[val][base_tup_dict[major_nt]]
						all_sites_info_string += "\t"+str(major_val)
					for val in range(1,10):
						minor_val = info_tup[val][base_tup_dict[minor_nt]]
						all_sites_info_string += "\t"+str(minor_val)
					all_sites_info_string += "\t"+str(neighbour_SNP_count)+"\t"+str(neighbour_over_count)
					all_sites_info_string += "\n"

					if prop_nt_list[1][0] >= min_proportion or prop_nt_list[0][1] != ref_seq[loc]:
						
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
						# iSNV_info_string += "\t"+str(minor_qual_pass)
						iSNV_info_string += "\t"+str(neighbour_SNP_count)+"\t"+str(neighbour_over_count)
						iSNV_info_string += "\n"


						# info_tup_dict = {'cov':0,'prop':1,'qual':2,'map_qual':3,'read_len':4,'read_loc':5,'read_R_loc':6,'read_prop_loc':7,'mismatch':8,'indel':9}
						nt_list = ['A','T','C','G']
						NS_counter = 0
						S_counter = 0
						pos_nt_list = []
						nt_count_list_sub = []
						S_nt_count_list_sub = []
						NS_nt_count_list_sub = []
						for nt_num in range(0,4):
							local_nt = nt_list[nt_num]
							local_sub_ID = segment+"-"+str(loc)+local_nt
							local_prop = allele_info_dict[segment][loc][info_tup_dict['prop']][nt_num]
							local_qual = allele_info_dict[segment][loc][info_tup_dict['qual']][nt_num]
							local_map_qual = allele_info_dict[segment][loc][info_tup_dict['map_qual']][nt_num]
							local_read_loc = allele_info_dict[segment][loc][info_tup_dict['read_loc']][nt_num]
							local_read_mismatch = allele_info_dict[segment][loc][info_tup_dict['mismatch']][nt_num]
							local_read_indel = allele_info_dict[segment][loc][info_tup_dict['indel']][nt_num]

							local_qual_pass = False
							if local_qual >= min_qual and local_map_qual >= min_avg_map_qual and local_read_loc >= min_avg_read_loc and local_read_mismatch <= max_avg_read_mismatch and local_read_indel <= max_avg_read_indel:
								local_qual_pass = True

							#for counting total S and NS sites
							if loc >= cds_start and loc <= cds_stop and neighbour_SNP_count <= max_neighbour_SNP_count:# and local_rep_qual_pass == True:# and major_qual_pass == True
								codon_number,codon_site = codon_finder(loc, cds_start)
								ref_acc_seg_codon = ref_seq[(codon_number*3+cds_start):(codon_number*3+cds_start)+3]
								ref_sub_codon = list(ref_acc_seg_codon)
								ref_sub_codon[codon_site] = local_nt
								ref_sub_codon = "".join(ref_sub_codon)
								sub_codon_site_dict[segment][loc] = codon_site
								if "-" not in ref_sub_codon and "N" not in ref_sub_codon and "-" not in ref_acc_seg_codon and "N" not in ref_acc_seg_codon:
									local_ns_type = "na"
									ref_aa = codon_to_aa_dict[ref_acc_seg_codon]
									sub_aa = codon_to_aa_dict[ref_sub_codon]
									codon_aa_dict[local_sub_ID] = (ref_aa,sub_aa)
									if ref_aa != sub_aa:
										S_NS_dict[segment][local_sub_ID] = "NS"
										if local_nt != ref_acc_seg_codon[codon_site]:
											NS_counter += 0.3333334
										if local_qual_pass == True:
											NS_nt_count_list_sub.append(local_prop*cov)
									else:
										S_NS_dict[segment][local_sub_ID] = "S"
										if local_nt != ref_acc_seg_codon[codon_site]:
											S_counter += 0.3333334
										if local_qual_pass == True:
											S_nt_count_list_sub.append(local_prop*cov)
							if local_qual_pass == True and local_prop >= min_proportion and neighbour_SNP_count <= max_neighbour_SNP_count:# and local_rep_qual_pass == True:# and major_qual_pass == True
								nt_count_list_sub.append(local_prop*cov)
								pos_nt_list.append(local_nt)
								local_allele_tup = (local_prop,local_qual,local_nt)
								try:
									sample_iSNV_dict[segment][loc].append(local_allele_tup)
								except:
									sample_iSNV_dict[segment][loc] = [local_allele_tup]

						pos_nt_list = list(set(pos_nt_list))
						seg_S_count += S_counter
						seg_NS_count += NS_counter
						include_site = False
						if len(pos_nt_list)>1:# and minor_prop >= min_proportion:
							include_site = True
						if include_site == True:
							segregating_site_count += 1
							segregating_site_count += 1
							try:
								D_sub += nt_diversity(nt_count_list_sub)
							except:
								D_sub += 0
							try:
								D_sub_S += nt_diversity(S_nt_count_list_sub)
							except:
								D_sub_S += 0
							try:
								D_sub_NS += nt_diversity(NS_nt_count_list_sub)
							except:
								D_sub_NS += 0
		all_sites_info = open(project_dir+"all_sites.base_info."+output_suffix+".txt","a")
		all_sites_info.write(all_sites_info_string)
		all_sites_info.close()
		all_sites_info_string = ''
		if D_sub > 0 and segregating_site_count > 0:
			pi_subs = D_sub/float(segregating_site_count)
		else:
			pi_subs = 0
		if D_sub_S > 0 and seg_S_count > 0:
			pi_S_subs = D_sub_S/float(round(seg_S_count,2))
		else:
			pi_S_subs = 0
		if D_sub_NS > 0 and seg_NS_count > 0:
			pi_NS_subs = D_sub_NS/float(round(seg_NS_count,2))
		else:
			pi_NS_subs = 0
		try:
			if len(avg_cov_dict[segment])>10:
				avg_cov = round(np.average(avg_cov_dict[segment]),2)
			else:
				avg_cov = 0.0
		except:
			avg_cov = 0.0

		nt_diversity_string += sampleID+"\t"+segment+"\t"+str(avg_cov)+"\t"+str(pass_cov_site_count)
		nt_diversity_string += "\t"+str(segregating_site_count)+"\t"+str(round(seg_S_count,2))+"\t"+str(round(seg_NS_count,2))
		nt_diversity_string += "\t"+str(round_to_n_sig_figs(D_sub,3))+"\t"+str(round_to_n_sig_figs(D_sub_S,3))+"\t"+str(round_to_n_sig_figs(D_sub_NS,3))
		nt_diversity_string += "\t"+str(round_to_n_sig_figs(pi_subs,3))+"\t"+str(round_to_n_sig_figs(pi_S_subs,3))+"\t"+str(round_to_n_sig_figs(pi_NS_subs,3))+"\n"
		
		nt_diversity_outfile = open(project_dir+"nt_diversity."+output_suffix+".txt","a")
		nt_diversity_outfile.write(nt_diversity_string)
		nt_diversity_outfile.close()
		nt_diversity_string = ''
		try:
			div_tup_old = div_tup
			div_tup = (D_sub+div_tup_old[0],D_sub_S+div_tup_old[1],D_sub_NS+div_tup_old[2],segregating_site_count+div_tup_old[3],seg_S_count+div_tup_old[4],seg_NS_count+div_tup_old[5])
			sample_nt_diversity_dict[segment] = div_tup
		except:
			div_tup = (D_sub,D_sub_S,D_sub_NS,segregating_site_count,seg_S_count,seg_NS_count)
			sample_nt_diversity_dict[segment] = div_tup

	iSNV_info = open(project_dir+"iSNV_info."+output_suffix+".txt","a")
	iSNV_info.write(iSNV_info_string)
	iSNV_info.close()

	out_tup = (sampleID,sample_iSNV_dict, avg_cov_dict, S_NS_dict, sample_nt_diversity_dict,sub_codon_site_dict,codon_aa_dict)
	return out_tup


processed_list = []
if parallel_process == True:
	processed_list = Parallel(n_jobs=num_cores)(delayed(identify_putative_iSNVs)(sampleID,segment_list,all_allele_info_dict[sampleID],neighbouring_SNP_dict[sampleID],max_neighbour_SNP_count,min_proportion,min_qual,min_avg_map_qual,min_avg_read_loc,max_avg_read_mismatch,max_avg_read_indel,min_major_map_qual,max_major_read_mismatch,max_major_read_indel,min_cov,cds_loc_dict,max_length_dict,sample_consensus_seq_dict,end_skip_len,sites_to_exclude,info_tup_dict,base_tup_dict,output_suffix) for sampleID in sampleID_list)
elif parallel_process == False:
	for sampleID in sampleID_list:
		output = identify_putative_iSNVs(sampleID,segment_list,all_allele_info_dict[sampleID],neighbouring_SNP_dict[sampleID],max_neighbour_SNP_count,min_proportion,min_qual,min_avg_map_qual,min_avg_read_loc,max_avg_read_mismatch,max_avg_read_indel,min_major_map_qual,max_major_read_mismatch,max_major_read_indel,min_cov,cds_loc_dict,max_length_dict,sample_consensus_seq_dict,end_skip_len,sites_to_exclude,info_tup_dict,base_tup_dict,output_suffix)
		processed_list.append(output)

iSNV_dict = {}
avg_cov_dict = {}
S_NS_dict = {}
sample_nt_diversity_dict = {} #(D_sub,D_sub_S,D_sub_NS,sub_seg_site_count,seg_S_count,seg_NS_count)
iSNV_codon_site_dict = {}
codon_allele_aa_dict = {}
for tup in processed_list:
	sampleID = tup[0]
	iSNV_dict[sampleID] = tup[1]
	avg_cov_dict[sampleID] = tup[2]
	S_NS_dict[sampleID] = tup[3]
	sample_nt_diversity_dict[sampleID] = tup[4]
	if iSNV_codon_site_dict == {}:
		iSNV_codon_site_dict = tup[5]
	else:
		temp_codon_dict = tup[5]
		for segment in temp_codon_dict:
			for loc in temp_codon_dict[segment]:
				try:
					iSNV_codon_site_dict[segment][loc] = temp_codon_dict[segment][loc]
				except:
					iSNV_codon_site_dict[segment] = {}
					iSNV_codon_site_dict[segment][loc] = temp_codon_dict[segment][loc]
	codon_allele_aa_dict[sampleID] = tup[6]
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
				avg_cov = 0.0
		except:
			avg_cov = 0.0
		outfile.write("\t"+str(avg_cov))
	outfile.write("\n")
outfile.close()



num_putative_iSNVs_found = 0
outfile = open(project_dir+"sublabel_prior."+output_suffix+".txt","w")
outfile.write("ref_sub_label\tsub_label\tsampleID\tsegment\tloc\tcov\tref_nt\tminor_nt\tmajor_nt\tref_prop\tminor_prop\tmajor_prop\tsubs\n")
sub_label_dict = {}
core_sub_label_dict = {}
for sampleID in iSNV_dict:
	sub_label_dict[sampleID] = {}
	for segment in iSNV_dict[sampleID]:
		try:
			neighbour_over_count = neighbouring_SNP_dict[sampleID][segment]["over"]
		except:
			neighbour_over_count = 0
		ref_seq = consensus_seq_dict[segment]
		# try:
		# except:

		if neighbour_over_count <= max_neighbour_SNP_count:
			for loc in iSNV_dict[sampleID][segment]:
				allele_list = sorted(iSNV_dict[sampleID][segment][loc],reverse=True)
				include_site = False
				ref_nt = ref_seq[loc]
				# if len(allele_list) == 2 and allele_list[0][2] != ref_nt and allele_list[1][2] != ref_nt:
				# 	print("ref not found:\t"+sampleID+"\t"+segment+"\t"+str(loc)+"\t"+str(allele_list))
				# elif len(allele_list) == 1 and allele_list[0][2] != ref_nt:
				# 	print("ref not found:\t"+sampleID+"\t"+segment+"\t"+str(loc)+"\t"+str(allele_list))
				
				if len(allele_list) == 2:  #(local_prop,local_qual,local_nt)
					include_site = True
					major_nt = allele_list[0][2]
					major_qual = allele_list[0][1]
					major_prop = allele_list[0][0]

					minor_nt = allele_list[1][2]
					minor_qual = allele_list[1][1]
					minor_prop = allele_list[1][0]

					# major_nt = ref_nt
					ref_prop = all_allele_info_dict[sampleID][segment][loc][1][base_tup_dict[ref_nt]]
					ref_qual = all_allele_info_dict[sampleID][segment][loc][2][base_tup_dict[ref_nt]]
				elif len(allele_list) == 1:
					if allele_list[0][2] != ref_nt:
						include_site = True
						
						minor_nt = allele_list[0][2]
						minor_qual = allele_list[0][1]
						minor_prop = allele_list[0][0]

						major_nt = ref_nt
						major_prop = all_allele_info_dict[sampleID][segment][loc][1][base_tup_dict[ref_nt]]
						major_qual = all_allele_info_dict[sampleID][segment][loc][2][base_tup_dict[ref_nt]]

						ref_prop = all_allele_info_dict[sampleID][segment][loc][1][base_tup_dict[ref_nt]]
						ref_qual = all_allele_info_dict[sampleID][segment][loc][2][base_tup_dict[ref_nt]]
						
				if include_site == True:
					cov = all_allele_info_dict[sampleID][segment][loc][0]
					sub_ID = segment+"-"+str(loc)+minor_nt
					sub_label = segment+"-"+major_nt+str(loc)+minor_nt
					ref_sub_label = segment+"-"+ref_nt+str(loc)+minor_nt

					cds_sub_list = []
					for CDS in all_cds_refs_dict[segment]:
						ref_cds_seq = all_cds_refs_dict[segment][CDS]
						try:
							cds_loc = all_cds_site_loc_dict[segment][CDS][loc]
						except:
							cds_loc = ''
						if cds_loc != '':
							ref_cds_nt = ref_cds_seq[cds_loc]
							cds_start = 0
							codon_number,codon_site = codon_finder(cds_loc, cds_start)
							ref_cds_codon = ref_cds_seq[(codon_number*3+cds_start):(codon_number*3+cds_start)+3]
							sub_cds_codon = list(ref_cds_codon)
							if minor_nt != ref_cds_nt and major_nt == ref_cds_nt:
								sub_cds_codon[codon_site] = minor_nt

							elif minor_nt == ref_cds_nt and major_nt != ref_cds_nt:
								sub_cds_codon[codon_site] = major_nt
							else:
								sub_cds_codon[codon_site] = minor_nt
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
						for t in range(0,len(cds_sub_list)):
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
					outfile.write(ref_sub_label+"\t"+sub_label+"\t"+sampleID+"\t"+segment+"\t"+str(loc)+"\t"+str(cov)+"\t"+ref_nt+"\t"+minor_nt+"\t"+major_nt+"\t"+str(ref_prop)+"\t"+str(minor_prop)+"\t"+str(major_prop)+"\t"+subtype_string+"\n")
					num_putative_iSNVs_found += 1
outfile.close()

print("Number of putative iSNVs found: "+str(num_putative_iSNVs_found))

nt_diversity_outfile = open(project_dir+"nt_diversity.all_segments_prior."+output_suffix+".txt","w")
nt_diversity_outfile.write("sampleID\tseg_site_count\tseg_S_count\tseg_NS_count\tD_sub\tD_sub_S\tD_sub_NS\tpi_sub\tpi_S\tpi_NS\n")
for sampleID in sample_nt_diversity_dict:
	D_sub = 0
	D_sub_S = 0
	D_sub_NS = 0
	sub_seg_site_count = 0
	seg_S_count = 0
	seg_NS_count = 0
	for s in range(0,len(segment_list)):
		segment = segment_list[s]
		div_tup = sample_nt_diversity_dict[sampleID][segment]
		D_sub += div_tup[0]
		D_sub_S += div_tup[1]
		D_sub_NS += div_tup[2]
		sub_seg_site_count += div_tup[3]
		seg_S_count += div_tup[4]
		seg_NS_count += div_tup[5]
	if D_sub > 0 and sub_seg_site_count > 0:
		pi_subs = D_sub/float(sub_seg_site_count)
	else:
		pi_subs = 0
	if D_sub_S > 0 and seg_S_count > 0:
		pi_S_subs = D_sub_S/float(round(seg_S_count,2))
	else:
		pi_S_subs = 0
	if D_sub_NS > 0 and seg_NS_count > 0:
		pi_NS_subs = D_sub_NS/float(round(seg_NS_count,2))
	else:
		pi_NS_subs = 0

	nt_diversity_outfile.write(sampleID+"\t"+
		str(sub_seg_site_count)+"\t"+str(round(seg_S_count,2))+"\t"+str(round(seg_NS_count,2))+"\t"+
		str(round_to_n_sig_figs(D_sub,3))+"\t"+str(round_to_n_sig_figs(D_sub_S,3))+"\t"+str(round_to_n_sig_figs(D_sub_NS,3))+"\t"+
		str(round_to_n_sig_figs(pi_subs,3))+"\t"+str(round_to_n_sig_figs(pi_S_subs,3))+"\t"+str(round_to_n_sig_figs(pi_NS_subs,3))+"\n")
nt_diversity_outfile.close()


#(D_sub,D_sub_S,D_sub_NS,sub_seg_site_count,seg_S_count,seg_NS_count)
nt_diversity_outfile = open(project_dir+"nt_diversity."+output_suffix+".txt","w")
for sample_root in replicate_iSNV_dict:
	nt_div_dict = {}
	D_sub,D_sub_S,D_sub_NS = 0,0,0
	sub_seg_site_count,seg_S_count,seg_NS_count = 0,0,0
	for sub_label in replicate_iSNV_dict[sample_root]:
		tup = replicate_iSNV_dict[sample_root][sub_label]
		minor_prop = tup[0]
		sub_type = tup[1]
		segment = sub_label.split("-")[0]
		major_nt = sub_label.split("-")[1][0]
		minor_nt = sub_label.split("-")[1][-1]
		loc = int(sub_label.split("-")[1][1:-1])
		cov = allele_info_merge_dict[sample_root][segment][loc]
		nt_cov_list = [minor_prop*cov,(1-minor_prop)*cov]
		sub_seg_site_count += 1
		try:
			D_sub += nt_diversity(nt_cov_list)
		except:
			D_sub += 0
		if sub_type == "S":
			seg_S_count += 1
			try:
				D_sub_S += nt_diversity(nt_cov_list)
			except:
				D_sub_S += 0
		
		elif sub_type == "NS":
			seg_NS_count += 1
			try:
				D_sub_NS += nt_diversity(nt_cov_list)
			except:
				D_sub_NS += 0
		
	if D_sub > 0 and sub_seg_site_count > 0:
		pi_subs = D_sub/float(sub_seg_site_count)
	else:
		pi_subs = 0
	if D_sub_S > 0 and seg_S_count > 0:
		pi_S_subs = D_sub_S/float(round(seg_S_count,2))
	else:
		pi_S_subs = 0
	if D_sub_NS > 0 and seg_NS_count > 0:
		pi_NS_subs = D_sub_NS/float(round(seg_NS_count,2))
	else:
		pi_NS_subs = 0
		nt_diversity_outfile.write(sample_root+"\t"+
			str(sub_seg_site_count)+"\t"+str(round(seg_S_count,2))+"\t"+str(round(seg_NS_count,2))+"\t"+
			str(round_to_n_sig_figs(pi_subs,3))+"\t"+str(round_to_n_sig_figs(pi_S_subs,3))+"\t"+str(round_to_n_sig_figs(pi_NS_subs,3))+"\n")
nt_diversity_outfile.close()


#Process consensus sequences to translate them, also to count the expectation for NS/S ratio
if skip_translate == False:
	outfile = open(project_dir+"S_NS_subcount."+output_suffix+".txt","w")
	nucl_seq_outfile = open(project_dir+"consensus_seqs."+output_suffix+".fna","w")
	seq_outfile = open(project_dir+"translated_seqs."+output_suffix+".faa","w")
	for acc_seg in sample_consensus_seq_dict:
		sampleID = acc_seg.split("\t")[0]
		segment = acc_seg.split("\t")[1]
		seq = sample_consensus_seq_dict[acc_seg]
		cds_start = int(cds_loc_dict[segment].split("-")[0])-1
		cds_stop = int(cds_loc_dict[segment].split("-")[1])-1
		ns_prop,s_prop,ns_count,s_count,viable_codon_count = sub_site_count(seq, cds_start, cds_stop,codon_to_aa_dict)
		outfile.write(segment+"\t"+sampleID+"\t"+str(ns_prop)+"\t"+str(s_prop)+"\t"+str(ns_count)+"\t"+str(s_count)+"\t"+str(viable_codon_count)+"\n")
		translated_seq = translate_cds(seq, cds_start, cds_stop,codon_to_aa_dict)
		seq_outfile.write(">"+sampleID+"_"+segment+"\n"+translated_seq+"\n")
		nucl_seq_outfile.write(">"+sampleID+"_"+segment+"\n"+seq+"\n")


	for seg_refs in ref_set_list:
		refs = seg_refs.split("-")
		for ref in refs:
			seq = ''
			segment = ref.split(".")[0]
			ref_segment = seg_refs.split("-")[0].split(".")[0]
			temp_infile = open(ref_dir+ref,"r")
			for line in temp_infile:
				line = line.strip()
				if line[0] != ">":
					seq += line
			cds_start = int(cds_loc_dict[ref_segment].split("-")[0])-1
			cds_stop = int(cds_loc_dict[ref_segment].split("-")[1])-1
			ns_prop,s_prop,ns_count,s_count,viable_codon_count = sub_site_count(seq, cds_start, cds_stop,codon_to_aa_dict)
			outfile.write(segment+"\tref"+"\t"+str(ns_prop)+"\t"+str(s_prop)+"\t"+str(ns_count)+"\t"+str(s_count)+"\t"+str(viable_codon_count)+"\n")
			translated_seq = translate_cds(seq, cds_start, cds_stop,codon_to_aa_dict)
			seq_outfile.write(">ref_"+segment+"\n"+translated_seq+"\n")
	outfile.close()
	seq_outfile.close()
	nucl_seq_outfile.close()
del skip_translate


if skip_writing_all_site_summaries == False:
	nt_list = ['A','T','C','G']
	num_allele = 1
	for num_allele in range(0,2):
		if num_allele == 1:
			cov_outfile = open(project_dir+"all_cov."+output_suffix+".txt","w")
		prop_outfile = open(project_dir+"all_sites.prop."+output_suffix+"."+str(num_allele)+".txt","w")
		readloc_outfile = open(project_dir+"all_sites.read_loc."+output_suffix+"."+str(num_allele)+".txt","w")
		qual_outfile = open(project_dir+"all_sites.qual."+output_suffix+"."+str(num_allele)+".txt","w")
		mapqual_outfile = open(project_dir+"all_sites.map_qual."+output_suffix+"."+str(num_allele)+".txt","w")
		mismatch_outfile = open(project_dir+"all_sites.mismatch."+output_suffix+"."+str(num_allele)+".txt","w")
		indel_outfile = open(project_dir+"all_sites.indel."+output_suffix+"."+str(num_allele)+".txt","w")
		if num_allele == 1:
			cov_outfile.write("segment\tloc")
		prop_outfile.write("segment\tloc")
		readloc_outfile.write("segment\tloc")
		qual_outfile.write("segment\tloc")
		mapqual_outfile.write("segment\tloc")
		mismatch_outfile.write("segment\tloc")
		indel_outfile.write("segment\tloc")

		for a in range(0,len(sampleID_list)):
			sampleID = sampleID_list[a]
			if num_allele == 1:
				cov_outfile.write("\t"+sampleID)
			prop_outfile.write("\t"+sampleID)
			readloc_outfile.write("\t"+sampleID)
			qual_outfile.write("\t"+sampleID)
			mapqual_outfile.write("\t"+sampleID)
			mismatch_outfile.write("\t"+sampleID)
			indel_outfile.write("\t"+sampleID)
		if num_allele == 1:
			cov_outfile.write("\n")
		prop_outfile.write("\n")
		readloc_outfile.write("\n")
		qual_outfile.write("\n")
		mapqual_outfile.write("\n")
		mismatch_outfile.write("\n")
		indel_outfile.write("\n")
		for seg_ref in ref_set_list:
			seg = seg_ref.split("-")[0]
			segment = seg.split(".")[0]
			for loc in range(0,max_length_dict[segment]):
				if num_allele == 1:
					cov_outfile.write(segment+"\t"+str(loc))
				prop_outfile.write(segment+"\t"+str(loc))
				readloc_outfile.write(segment+"\t"+str(loc))
				qual_outfile.write(segment+"\t"+str(loc))
				mapqual_outfile.write(segment+"\t"+str(loc))
				mismatch_outfile.write(segment+"\t"+str(loc))
				indel_outfile.write(segment+"\t"+str(loc))
				for a in range(0,len(sampleID_list)):
					sampleID = sampleID_list[a]
					cov_thresh = min_cov

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
						cov_outfile.write("\t"+str(cov))
					prop_outfile.write("\t"+str(allele_freq))
					readloc_outfile.write("\t"+str(read_loc))
					qual_outfile.write("\t"+str(qual))
					mapqual_outfile.write("\t"+str(map_qual))
					mismatch_outfile.write("\t"+str(mismat))
					indel_outfile.write("\t"+str(indel))
				if num_allele == 1:
					cov_outfile.write("\n")
				prop_outfile.write("\n")
				readloc_outfile.write("\n")
				qual_outfile.write("\n")
				mapqual_outfile.write("\n")
				mismatch_outfile.write("\n")
				indel_outfile.write("\n")
		if num_allele == 1:
			cov_outfile.close()
		prop_outfile.close()
		readloc_outfile.close()
		qual_outfile.close()
		mapqual_outfile.close()
		mismatch_outfile.close()
		indel_outfile.close()


sub_label_list = []
for sampleID in iSNV_dict:
	for segment in iSNV_dict[sampleID]:
		try:
			neighbour_over_count = neighbouring_SNP_dict[sampleID][segment]["over"]
		except:
			neighbour_over_count = 0

		if neighbour_over_count <= 20:
			for loc in iSNV_dict[sampleID][segment]:
				allele_list = sorted(iSNV_dict[sampleID][segment][loc],reverse=True)
				if len(allele_list) == 2:  #(local_prop,local_qual,local_nt)
					major_nt = allele_list[0][2]
					major_qual = allele_list[0][1]
					major_prop = allele_list[0][0]

					minor_nt = allele_list[1][2]
					minor_qual = allele_list[1][1]
					minor_prop = allele_list[1][0]

					sub_ID = segment+"-"+str(loc)+minor_nt
					sub_label = segment+"-"+major_nt+str(loc)+minor_nt
					try:
						sub_type = S_NS_dict[sampleID][segment][sub_ID]
					except:
						sub_type = "na"
					sub_label_list.append(sub_label)

sub_label_list = sorted(list(set(sub_label_list)))


allele_pass_dict = {}
outfile = open(project_dir+"iSNV_table."+output_suffix+".txt","w")
for j in range(0,len(sub_label_list)):
	outfile.write("\t"+sub_label_list[j])
outfile.write("\n")
for i in range(0,len(sampleID_list)):
	sampleID = sampleID_list[i]
	subject = sampleID.split("-")[0]
	cov_thresh = min_cov
	outfile.write(sampleID)
	for j in range(0,len(sub_label_list)):
		sub_label = sub_label_list[j]
		segment = sub_label.split("-")[0]
		major_nt = sub_label.split("-")[1][0]
		minor_nt = sub_label.split("-")[1][-1]
		loc = int(sub_label.split("-")[1][1:-1])
		alt_sub_label = segment+"-"+minor_nt+str(loc)+major_nt
		minor_sub_prop = "#"
		try:
			site_cov = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['cov']]
		except:
			site_cov = 0
		if site_cov < cov_thresh:
			minor_sub_prop = "nan"
		else:
			try:
				minor_sub_prop = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['prop']][base_tup_dict[minor_nt]]
				minor_sub_qual = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['qual']][base_tup_dict[minor_nt]]
				minor_sub_map_qual = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['map_qual']][base_tup_dict[minor_nt]]
				minor_sub_read_loc = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['read_loc']][base_tup_dict[minor_nt]]
				minor_sub_read_mismatch = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['mismatch']][base_tup_dict[minor_nt]]
				minor_sub_read_indel = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['indel']][base_tup_dict[minor_nt]]
				
				# major_sub_prop = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['prop']][base_tup_dict[major_nt]]
				# major_sub_qual = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['qual']][base_tup_dict[major_nt]]
				# major_sub_map_qual = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['map_qual']][base_tup_dict[major_nt]]
				# major_sub_read_loc = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['read_loc']][base_tup_dict[major_nt]]
				# major_sub_read_mismatch = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['mismatch']][base_tup_dict[major_nt]]
				# major_sub_read_indel = all_allele_info_dict[sampleID][segment][loc][info_tup_dict['indel']][base_tup_dict[major_nt]]
				
				minor_qual_pass = False
				if minor_sub_qual >= min_qual and minor_sub_map_qual >= min_avg_map_qual and minor_sub_read_loc >= min_avg_read_loc and minor_sub_read_mismatch <= max_avg_read_mismatch and minor_sub_read_indel <= max_avg_read_indel:
					minor_qual_pass = True
				# if major_sub_qual >= min_qual and major_sub_map_qual >= min_major_map_qual and major_sub_read_loc >= min_avg_read_loc and major_sub_read_mismatch <= max_major_read_mismatch and major_sub_read_indel <= max_major_read_indel:
				# 	major_qual_pass = True
				# if minor_sub_prop < min_proportion:

				if minor_qual_pass == True:
					if minor_prop > 0:
						minor_sub_prop = round_to_n_sig_figs(minor_sub_prop,3)
					else:
						minor_sub_prop = 0
				else:
					minor_sub_prop = 0

			except:
				minor_sub_prop = "nan"
		outfile.write("\t"+str(minor_sub_prop))
		if minor_prop != "nan":
			try:
				allele_pass_dict[subject][sub_label] = minor_prop
			except:
				allele_pass_dict[subject] = {}
				allele_pass_dict[subject][sub_label] = minor_prop
			allele_pass_dict[subject][alt_sub_label] = (1.0-minor_prop)
	outfile.write("\n")
outfile.close()
