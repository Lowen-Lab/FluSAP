import os
import sys
import numpy as np
# import random

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

min_coverage = 100
min_major_allele_prop = 0.6
min_major_allele_qual = 30

max_N_count_dict = {'A_HA_H1':0,'A_HA_H3':0,'A_NA_N1':0,'A_NA_N2':0,'A_M':0,'A_NP':0,'A_NS':0,'A_PA':0,'A_PB1':0,'A_PB2':0,'B_HA':9,'B_M':5,'B_NA':0,'B_NP':0,'B_NS':3,'B_PA':0,'B_PB1':0,'B_PB2':0}
# max_N_count_dict = {'HA':0,'NA':0,'M':0,'NP':0,'NS':0,'PA':0,'PB1':0,'PB2':0}

ref_set_list = ['A_HA_H1.fa','A_HA_H3.fa','A_M_1.fa-A_M_2.fa','A_NA_N1.fa','A_NA_N2.fa','A_NP_1.fa-A_NP_2.fa','A_NS.fa','A_PA_1.fa-A_PA_2.fa','A_PB1_1.fa-A_PB1_2.fa','A_PB2_1.fa-A_PB2_2.fa','B_HA.fa','B_M.fa','B_NA.fa','B_NP.fa','B_NS.fa','B_PA.fa','B_PB1.fa','B_PB2.fa']
# ref_set_list = ['PB2.fa','PB1.fa','PA.fa','NP.fa','M.fa','NS.fa','HA.fa','NA.fa']


########################################################################################

def run_command(command):
	os.system(command)
	# subprocess.check_output(command,stderr=subprocess.STDOUT)
	return command

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

def trim_trailing_Ns(seq_in):
	first_non_N_site = -1
	last_non_N_site = -1
	for loc in range(0,len(seq_in)):
		if first_non_N_site == -1:
			nt = seq_in[loc]
			if nt == "A" or nt == "T" or nt == "G" or nt == "C":
				first_non_N_site = loc
	for loc in range(len(seq_in)-1,-1,-1):
		if last_non_N_site == -1:
			nt = seq_in[loc]
			if nt == "A" or nt == "T" or nt == "G" or nt == "C":
				last_non_N_site = loc
	seq_out = seq_in[first_non_N_site:last_non_N_site+1]
	return seq_out

########################################################################################


consensus_seq_dict = {}
max_length_dict = {}
infile = open(ref_dir+"all_refs.fa","r")
for line in infile:
	line = line.strip()
	if line[0] == ">":
		head = line[1:len(line)]
		segment = head#.split("-")[1]#+"_"+head.split("_")[1]
	else:
		try:
			consensus_seq_dict[segment] += line
		except:
			consensus_seq_dict[segment] = line
for segment in consensus_seq_dict:
	max_length_dict[segment] = len(consensus_seq_dict[segment])


sampleID_list = pull_sample_list(temp_dir,input_read_dir)

consensus_outfile = open("./consensus_seqs/"+output_suffix+".consensus_seqs.fasta","w")
for seg_ref in ref_set_list:
	ref_filename = seg_ref.split("-")[0]
	segment = ref_filename.split(".f")[0]#.split("_")[1]
	segment_consensus_outfile = open("./consensus_seqs/"+"consensus_seq."+segment+"."+output_suffix+".fasta","w")
	write_count = 0
	N_count_list = {}
	max_N_count = 12#max_N_count_dict[segment]
	for sampleID in sampleID_list:
		substitutions_filename = sub_dir+segment+"/"+sampleID+'.'+segment+".substitions.txt"
		if os.path.isfile(substitutions_filename) == True:
			sub_infile = open(substitutions_filename,"r")
			allele_dict = {}
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

						qual_dict = {"A":A_qual,"T":T_qual,"C":C_qual,"G":G_qual}
						allele_list = [(A_prop,"A"),(T_prop,"T"),(C_prop,"C"),(G_prop,"G")]
						allele_list = sorted(allele_list,reverse=True)
						major_allele = allele_list[0][1]
						major_prop = allele_list[0][0]
						major_qual = qual_dict[major_allele]

						if cov >= min_coverage and major_prop >= min_major_allele_prop and major_qual >= min_major_allele_qual:
							allele_dict[loc] = major_allele

			seq_string = ''
			for loc in range(0,max_length_dict[segment]):
				try:
					nt = allele_dict[loc]
				except:
					nt = "N"
				seq_string += nt
			seq_string = trim_trailing_Ns(seq_string)
			if len(seq_string) > (max_length_dict[segment]*0.9):
				N_count = len(seq_string)-len(seq_string.replace("N",""))
				if N_count > 20:
					N_count_val = '>20'
				else:
					N_count_val = N_count
				try:
					N_count_list[N_count_val] += 1
				except:
					N_count_list[N_count_val] = 1
				if N_count <= max_N_count:
					consensus_outfile.write(">"+sampleID+"_"+segment+"\n"+seq_string+"\n")
					segment_consensus_outfile.write(">"+sampleID+"\n"+seq_string+"\n")
					write_count += 1
	print(segment+"\t"+str(write_count))
	# for N_count_val in N_count_list:
	# 	print("\t"+str(N_count_val)+"\t"+str(N_count_list[N_count_val]))
	segment_consensus_outfile.close()
consensus_outfile.close()
