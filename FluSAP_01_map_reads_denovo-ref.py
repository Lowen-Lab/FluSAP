import os
import sys
import math
import random
import numpy as np
import time

parallel_process = True
parallel_max_cpu = 20
nhmmer_cores_per_job = 5

####
#Set up parallel processing if enabled
if parallel_process == True:
	import multiprocessing
	from joblib import Parallel, delayed
	max_cores = multiprocessing.cpu_count()
	if parallel_max_cpu == 0:
		num_cores = max_cores
	elif parallel_max_cpu > 0:
		num_cores = min(parallel_max_cpu,max_cores)
	elif parallel_max_cpu < 0:
		num_cores = max_cores+parallel_max_cpu
else:
	num_cores = 1
# nhmmer_cores_per_job = max(1,int((max_cores-10)/(num_cores))-1)
####

##################################  USER DEFINED VARIABLES  #######################################
project_dir = "./"
input_read_dir =project_dir+"raw_data/"
ref_dir = project_dir+"refs/"
hmm_dir = ref_dir+"hmm/"
trim_dir = project_dir+"trim/"
map_folder = 'rep_map'
map_dir = project_dir+map_folder+"/"
temp_dir = map_dir+"temp/"
qual_dir = map_dir+"qual/"
sam_dir = map_dir+"sam/"
bam_dir = map_dir+"bam/"
nhmmer_dir = map_dir+"nhmmer/"
sub_dir = map_dir+"substitutions/"
consensus_dir = map_dir+"consensus/"


nhmmer_options = '--notextw -E 1e-10 --cpu '+str(nhmmer_cores_per_job) #--popen 0.02

raw_read_input_type = "mixed" # 'paired', 'single', or 'mixed' --- When 'single' the pipeline will look for barcodes in one file, when 'paired' the pipeline will look for forward and reverse reads
forward_read_suffix = "_R1.fastq" #if empty, the script will attempt to predict these values, but it only works under specific assumptions
reverse_read_suffix = "_R2.fastq" #this value is ignored if 'reads_already_screened_for_quality' is 'True' or 'raw_read_input_type' is 'single' - NOTE: this variable cannot be removed without causing missing value errors


command_prefix = "bash "
repair_options = 'fint=t repair=t overwrite=t threads=1 -da -Xmx3000m'
bbduk_adapter_trim_options = 'ktrim=r k=23 mink=11 hdist=1 minlen=80 tpe tbo overwrite=t threads=4 -Xmx1000m'
bbmerge_options = 'qin=33 ecco mix overwrite=t threads=4 -Xmx1000m'
bbduk_options = 'qin=33 tpe=f tbo=f overwrite=t threads=4 -Xmx1000m'#maq=30  lhist=lhist.txt qhist=qhist.txt bhist=bhist.txt aqhist=aqhist.txt'
bbmap_options = 'qin=33 local=f touppercase=t overwrite=t threads=4 nodisk -Xmx2000m'#ambiguous=toss 

###################################################################################################
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

def fastq_to_fasta(filename_in,filename_out):
	line_counter = -1
	infile = open(filename_in,"r")
	outfile = open(filename_out,"w")
	for line in infile:
		line = line.strip()
		line_counter += 1
		if line_counter == 0:
			outfile.write(">"+line[1:len(line)].replace(" ","___")+"\n")
		elif line_counter == 1:
			outfile.write(nucleotide_character_edit(line)+"\n")
		elif line_counter == 3:
			line_counter = -1
	infile.close()
	outfile.close()
	return filename_out

def load_read_counts(sampleID,segment_list,trim_dir):
	read_count_dict = {}
	for segmentID in segment_list:
		try:
			read_count_file = open(trim_dir+segmentID+"/"+sampleID+'.'+segmentID+'.read_count.txt','r')
			for line in read_count_file:
				line = line.strip()
				read_count_dict[segmentID] = int(line)
			read_count_file.close()
		except:
			read_count_dict[segmentID] = 0
	return read_count_dict

def check_adapter(filename_in):
	infile = open(filename_in,'r')
	adapter_pass = False
	for line in infile:
		line = line.strip()
		if len(line) >0:
			if line[0] == ">":
				header = line
			else:
				seq = line
				if seq != "N":
					try:
						Nprop = (len(seq)-len(seq.replace("N","")))/len(seq)
					except:
						Nprop = 1.0
					if Nprop <=0.25:
						adapter_pass = True
	return adapter_pass

def check_if_input_reads_exist(sampleID,input_read_dir,local_raw_read_input_type,forward_read_suffix,reverse_read_suffix):
	input_reads_found = False
	if local_raw_read_input_type == "single":
		if os.path.isfile(input_read_dir+sampleID+forward_read_suffix) == True:
			input_reads_found = True
	elif local_raw_read_input_type == "paired":
		if os.path.isfile(input_read_dir+sampleID+forward_read_suffix) == True and os.path.isfile(input_read_dir+sampleID+reverse_read_suffix) == True:
			input_reads_found = True
	else:
		print("Unable to find input reads for sample: " +sampleID)
	return input_reads_found

def check_trimmed_reads_exist(sampleID,trim_dir):
	trimmed_reads_found = os.path.isfile(trim_dir+sampleID+".fq") == True
	return trimmed_reads_found

def check_if_reads_mapped(sampleID,segment_list,sam_dir,read_count_dict,ref_suffix=''):
	reads_mapped = True
	missing_ref_sets = []
	for segmentID in segment_list:
		if os.path.isfile(sam_dir+segmentID+"/"+sampleID+"."+segmentID+ref_suffix+'.sam') == False:
			if read_count_dict[segmentID] >= 1000:
				reads_mapped = False
				missing_ref_sets.append(segmentID)
	missing_ref_sets = list(set(missing_ref_sets))
	return reads_mapped,missing_ref_sets

def check_if_ref_refined(sampleID,consensus_dir,read_count_dict):
	ref_refined = True
	if os.path.isfile(consensus_dir+"all_segments/"+sampleID+'.refine.fa') == False:
		ref_refined = False
	return ref_refined

def determine_raw_read_input_type(sampleID,input_read_dir,forward_read_suffix,reverse_read_suffix):
	forward_input_reads_found = False
	reverse_input_reads_found = False
	if os.path.isfile(input_read_dir+sampleID+forward_read_suffix) == True:
		forward_input_reads_found = True
	if os.path.isfile(input_read_dir+sampleID+reverse_read_suffix) == True:
		reverse_input_reads_found = True

	if forward_input_reads_found == True and reverse_input_reads_found == True:
		local_raw_read_input_type = "paired"
		return local_raw_read_input_type
	elif forward_input_reads_found == True and reverse_input_reads_found == False:
		local_raw_read_input_type = "single"
		return local_raw_read_input_type
	else:
		sys.exit("Unable to find input reads for sample: " +sampleID)

def raw_read_processing_single(sampleID,input_read_dir,ref_dir,trim_dir,temp_dir,forward_read_suffix,command_prefix):
	default_adapter_file_path = ref_dir+'adapters.fa'
	input_reads_forward = input_read_dir+sampleID+forward_read_suffix
	trimadapt_forward = temp_dir+sampleID+'.trimadapt.fq'
	temp_adapter_filename = temp_dir+sampleID+".adapters.fa"
	clean_forward = temp_dir+sampleID+'.repair'+forward_read_suffix
	permanent_fq_forward = trim_dir+sampleID+".fq"
	
	command = command_prefix+'bbduk.sh in="'+input_reads_forward+'" out="'+trimadapt_forward+'" ref="'+default_adapter_file_path+'" '+bbduk_adapter_trim_options
	command_status = run_command(command)
	if command_status == False:
		sys.exit("bbduk.sh failed: "+sampleID)
	time.sleep(1)

	command = command_prefix+'bbmerge.sh in="'+trimadapt_forward+'" outadapter="'+temp_adapter_filename+'" '+bbmerge_options
	out = run_command(command)
	time.sleep(1)
	adapter_pass = check_adapter(temp_adapter_filename)
	if adapter_pass == True:
		command = command_prefix+'bbduk.sh in="'+trimadapt_forward+'" out="'+clean_forward+'" adapters="'+temp_adapter_filename+'" '+bbmerge_options
		command_status = run_command(command)
		if command_status == False:
			sys.exit("bbduk.sh failed: "+temp_adapter_filename)
		time.sleep(1)
	else:
		clean_forward = trimadapt_forward

	os.rename(clean_forward,permanent_fq_forward)
	if adapter_pass == True:
		os.remove(trimadapt_forward)
	os.remove(temp_adapter_filename)
	return permanent_fq_forward

def raw_read_processing_paired(sampleID,input_read_dir,ref_dir,trim_dir,temp_dir,forward_read_suffix,reverse_read_suffix,command_prefix):
	input_reads_forward = input_read_dir+sampleID+forward_read_suffix
	input_reads_reverse = input_read_dir+sampleID+reverse_read_suffix
	
	default_adapter_file_path = ref_dir+'adapters.fa'
	repair_forward = temp_dir+sampleID+'.repair_R1.fq'
	repair_reverse = temp_dir+sampleID+'.repair_R2.fq'
	trimadapt_forward = temp_dir+sampleID+'.trimadapt_R1.fq'
	trimadapt_reverse = temp_dir+sampleID+'.trimadapt_R2.fq'
	temp_adapter_filename = temp_dir+sampleID+".adapters.fa"
	clean_forward = temp_dir+sampleID+'.clean_R1.fq'
	clean_reverse = temp_dir+sampleID+'.clean_R2.fq'
	clean_repair = temp_dir+sampleID+'.clean.fq'
	permanent_fq = trim_dir+sampleID+".fq"

	command = command_prefix+'repair.sh in="'+input_reads_forward+'" in2="'+input_reads_reverse+'" out="'+repair_forward+'" out2="'+repair_reverse+'" '+repair_options
	command_status = run_command(command)#,"verbose")
	if command_status == False:
		sys.exit("repair.sh failed: "+command)
	time.sleep(1)

	command = command_prefix+'bbduk.sh in="'+repair_forward+'" in2="'+repair_reverse+'" out="'+trimadapt_forward+'" out2="'+trimadapt_reverse+'" ref="'+default_adapter_file_path+'" '+bbduk_adapter_trim_options
	command_status = run_command(command)
	if command_status == False:
		sys.exit("bbduk.sh failed: "+command)
	time.sleep(1)

	command = command_prefix+'bbmerge.sh in="'+trimadapt_forward+'" in2="'+trimadapt_reverse+'" outadapter="'+temp_adapter_filename+'" '+bbmerge_options
	command_status = run_command(command)
	if command_status == False:
		sys.exit("bbmerge.sh -outadapt failed: "+command)
	time.sleep(1)

	adapter_pass = check_adapter(temp_adapter_filename)
	if adapter_pass == True:
		command = command_prefix+'bbduk.sh in="'+trimadapt_forward+'" in2="'+trimadapt_reverse+'" out="'+clean_forward+'" out2="'+clean_reverse+'" adapters="'+temp_adapter_filename+'" '+bbmerge_options
		command_status = run_command(command)
		if command_status == False:
			sys.exit("bbmerge.sh failed: "+command)
		time.sleep(1)
	else:
		clean_forward = trimadapt_forward
		clean_reverse = trimadapt_reverse

	command = command_prefix+'repair.sh in="'+clean_forward+'" in2="'+clean_reverse+'" out="'+clean_repair+'" '+repair_options
	command_status = run_command(command)
	if command_status == False:
		sys.exit("repair.sh failed: "+command)
	time.sleep(1)

	os.rename(clean_repair,permanent_fq)
	os.remove(repair_forward)
	os.remove(repair_reverse)
	if adapter_pass == True:
		os.remove(trimadapt_forward)
		os.remove(trimadapt_reverse)
	os.remove(clean_forward)
	os.remove(clean_reverse)
	os.remove(temp_adapter_filename)
	return permanent_fq

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

def parse_nhmmer_alignment(sampleID,alignment_filepath,read_match_dict,min_read_align_len=80,min_depth=20):
	alignment_length = {}
	allele_dict = {}
	infile = open(alignment_filepath,"r")
	for line in infile:
		line = line.rstrip('\n')
		if len(line) >= 8:
			if line[0:7] == '#=GF ID':
				current_segment = line[8:len(line)]
				alignment_length[current_segment] = 0
				allele_dict[current_segment] = {}
		if len(line) > 300:
			if line[0] != '#':
				split_line = split_multi_space(line)
				readID = split_line[0].split("|")[-1]
				try:
					read_segment = read_match_dict[readID.split("/")[0]]
				except:
					read_segment = ''
				if read_segment == current_segment:
					align_loc_list = readID.split("/")[-1].split('-')
					align_string = split_line[1]
					read_len_aligned = max(int(align_loc_list[1]),int(align_loc_list[0]))-min(int(align_loc_list[1]),int(align_loc_list[0]))
					if alignment_length[current_segment] == 0:
						alignment_length[current_segment] = len(align_string)
					first_base_found = False
					if read_len_aligned >= min_read_align_len:
						search_sites_remaining = read_len_aligned+2
						for loc in range(0,alignment_length[current_segment]):
							if search_sites_remaining > 0:
								try:
									char = align_string[loc]
								except:
									sys.exit(alignment_filepath+" "+str(loc)+' '+str(alignment_length[current_segment])+' '+read_segment+' '+current_segment+' '+str(read_len_aligned)+'\n'+line)
								if first_base_found == False:
									if char != "-" and char != ".":
										if char == "A" or char == "T" or char == "C" or char == "G":
											first_base_found = True
								if first_base_found == True:
									if char != "-":
										if char != ".":
											search_sites_remaining -= 1
										if search_sites_remaining >= 0:
											try:
												allele_dict[current_segment][loc][char] += 1
											except:
												try:
													allele_dict[current_segment][loc][char] = 1
												except:
													allele_dict[current_segment][loc] = {}
													allele_dict[current_segment][loc][char] = 1
	char_count_string = 'segment\tloc\tA\tT\tC\tG\tgap\n'
	seq_out_dict = {}
	for current_segment in allele_dict:
		seq_out = ''
		for loc in range(0,alignment_length[current_segment]):
			try:
				A_count = allele_dict[current_segment][loc]['A']
			except:
				A_count = 0
			try:
				T_count = allele_dict[current_segment][loc]['T']
			except:
				T_count = 0
			try:
				C_count = allele_dict[current_segment][loc]['C']
			except:
				C_count = 0
			try:
				G_count = allele_dict[current_segment][loc]['G']
			except:
				G_count = 0
			try:
				dot_count = allele_dict[current_segment][loc]['.']
			except:
				dot_count = 0
			char_count_string += current_segment+"\t"+str(loc)+"\t"+str(A_count)+"\t"+str(T_count)+"\t"+str(C_count)+"\t"+str(G_count)+"\t"+str(dot_count)+"\n"
			ATCG_count = A_count+T_count+C_count+G_count
			if ATCG_count >0:
				base_prop = ATCG_count/(dot_count+ATCG_count)
			else:
				base_prop = 0.0
			if dot_count >0:
				dot_prop = dot_count/(dot_count+ATCG_count)
			else:
				dot_prop = 0
			if dot_prop < 0.01 and ATCG_count >= min_depth:
				nt_list = [(A_count,'A'),(T_count,'T'),(C_count,'C'),(G_count,'G')]
				nt_list = sorted(nt_list, reverse=True)
				major_nt = nt_list[0][1]
				seq_out += major_nt
			elif dot_prop < 0.01:
				seq_out += '-'
		if seq_out != '':
			seq_out_dict[current_segment] = seq_out

	if seq_out_dict != {}:
		outfile = open(alignment_filepath.replace(".txt",".summary.txt"),"w")
		outfile.write(char_count_string)
		outfile.close()

		temp_fasta_name = temp_dir+sampleID+'.rough_consensus.fa'
		storage_fasta_name = consensus_dir+"all_segments/"+sampleID+'.rough_consensus.fa'
		temp_fasta = open(temp_fasta_name,"w")
		for segmentID in seq_out_dict:
			temp_fasta.write(">"+segmentID+"\n"+seq_out_dict[segmentID]+"\n")
		temp_fasta.close()
		os.rename(temp_fasta_name,storage_fasta_name)
	return seq_out_dict

def parse_nhmmer_tblout(tblout_filepath):
	hit_dict = {}
	infile = open(tblout_filepath,"r")
	for line in infile:
		line = line.rstrip('\n')
		if len(line) > 1:
			if line[0] != "#":
				split_line = split_multi_space(line)
				try:
					evalue = float(split_line[12])
				except:
					evalue = 9.9
				# score = float(split_line[13])
				# bias = float(split_line[14])
				if evalue <= 1e-5:
					readID = split_line[0]
					segmentID = split_line[2]
					hit_tup = (evalue,segmentID)
					try:
						hit_dict[readID].append(hit_tup)
					except:
						hit_dict[readID] = [hit_tup]
	infile.close()
	segment_sorted_dict = {}
	segments_dict = {}
	for readID in hit_dict:
		read_hits = sorted(hit_dict[readID])
		# top_hit_evalue = read_hits[0][0]
		top_hit_segment = read_hits[0][1]
		segment_sorted_dict[readID] = top_hit_segment
		segments_dict[top_hit_segment] = ''
	segment_list = list(segments_dict.keys())
	del hit_dict
	return segment_sorted_dict,segment_list

def split_fastq_by_segment(sampleID,segment_list,readID_to_segment_dict,temp_dir,trim_dir,ref_dir,nhmmer_dir,read_input_type):
	alignment_filepath = nhmmer_dir+sampleID+'.nhmmer_align.txt'
	temp_alignment_filepath = temp_dir+sampleID+'.nhmmer_align.txt'
	tblout_filepath = nhmmer_dir+sampleID+'.nhmmer_tblout.txt'
	temp_tblout_filepath = temp_dir+sampleID+'.nhmmer_tblout.txt'
	fastq_in_filepath = trim_dir+sampleID+".fq"
	fasta_temp_filepath = temp_dir+sampleID+".fasta"

	other_lineout = ''
	lineout_dict = {}
	read_count_dict = {}
	other_count = 0
	for current_segment in segment_list:
		if not os.path.exists(trim_dir+current_segment):
			os.makedirs(trim_dir+current_segment)
		temp_fq = temp_dir+sampleID+'.'+current_segment+'.fq'
		ref_outfile = open(temp_fq,"w")
		ref_outfile.close()
		lineout_dict[current_segment] = ''
		read_count_dict[current_segment] = 0
	temp_fq = temp_dir+sampleID+'.nhmmer_remain.fq'
	ref_outfile = open(temp_fq,"w")
	ref_outfile.close()
		
	fastq_infile = open(fastq_in_filepath,"r")
	line_counter = -1
	current_segment = ''
	for line in fastq_infile:
		line = line.rstrip('\n')
		if len(line) >= 1:
			line_counter += 1
			if line_counter == 0:
				if line[0]!="@":
					sys.exit("fastq parse failed: read header expected:\n"+line)
			elif line_counter == 2: # plus
				if line[0]!="+":
					sys.exit("fastq parse failed: '+' expected:\n"+line)

			if line_counter == 0: # read header
				readID = line[1:len(line)].replace(" ","___")
				try:
					current_segment = readID_to_segment_dict[readID]
				except:
					current_segment = ''
					other_count += 1
			elif line_counter == 3: # qual string
				line_counter = -1

			if current_segment != '':
				lineout_dict[current_segment] += line+"\n"
				if line_counter == 0:
					read_count_dict[current_segment] += 1
				elif line_counter == -1:
					if len(lineout_dict[current_segment]) >= 100000:
						temp_fq = temp_dir+sampleID+'.'+current_segment+'.fq'
						ref_outfile = open(temp_fq,"a")
						ref_outfile.write(lineout_dict[current_segment])
						ref_outfile.close()
						lineout_dict[current_segment] = ''
						# print('intermediate write fq lines: '+current_segment)
			elif current_segment == '':
				other_lineout += line+"\n"
				if len(other_lineout) >= 100000:
					temp_fq = temp_dir+sampleID+'.nhmmer_remain.fq'
					ref_outfile = open(temp_fq,"a")
					ref_outfile.write(other_lineout)
					ref_outfile.close()
					other_lineout = ''

	fastq_infile.close()
	if other_lineout != '':
		temp_fq = temp_dir+sampleID+'.nhmmer_remain.fq'
		permanent_fq = trim_dir+"remain/"+sampleID+'.nhmmer_remain.fq'
		ref_outfile = open(temp_fq,"a")
		ref_outfile.write(other_lineout)
		ref_outfile.close()
		other_lineout = ''
		os.rename(temp_fq,permanent_fq)
		temp_read_count_filename = temp_dir+sampleID+'.nhmmer_remain.read_count.txt'
		permanent_read_count_filename = trim_dir+"remain/"+sampleID+'.nhmmer_remain.read_count.txt'
		read_count_file = open(temp_read_count_filename,"w")
		read_count_file.write(str(other_count)+"\n")
		read_count_file.close()
		os.rename(temp_read_count_filename,permanent_read_count_filename)
		
	for current_segment in lineout_dict:
		temp_fq = temp_dir+sampleID+'.'+current_segment+'.fq'
		repair_fq = temp_dir+sampleID+'.'+current_segment+'.repair.fq'
		permanent_fq = trim_dir+current_segment+"/"+sampleID+'.'+current_segment+'.fq'
		
		temp_read_count_filename = temp_dir+sampleID+'.'+current_segment+'.read_count.txt'
		permanent_read_count_filename = trim_dir+current_segment+"/"+sampleID+'.'+current_segment+'.read_count.txt'
		
		if lineout_dict[current_segment] != '':
			ref_outfile = open(temp_fq,"a")
			ref_outfile.write(lineout_dict[current_segment])
			ref_outfile.close()
			lineout_dict[current_segment] = ''
		
		try:
			seg_read_count = read_count_dict[current_segment]
		except:
			seg_read_count = 0
		read_count_file = open(temp_read_count_filename,"w")
		read_count_file.write(str(seg_read_count)+"\n")
		read_count_file.close()
		if os.path.isfile(permanent_read_count_filename) == True:
			os.remove(permanent_read_count_filename)
		os.rename(temp_read_count_filename,permanent_read_count_filename)
		if seg_read_count >= 100:
			if read_input_type == 'paired':
				repair_command = 'repair.sh in="'+temp_fq+'" out="'+repair_fq+'" '+repair_options
				command_status = run_command(repair_command)
				if command_status == False:
					sys.exit("Repair.sh failed on segment split file: "+sampleID+" "+current_segment)
				os.remove(temp_fq)
				if os.path.isfile(permanent_fq) == True:
					os.remove(permanent_fq)
				os.rename(repair_fq,permanent_fq)
			elif read_input_type == 'single':
				os.rename(temp_fq,permanent_fq)
		else:
			os.remove(temp_fq)
	return read_count_dict

def load_alignment(filepath):
	seq_dict = {}
	infile = open(filepath,"r")
	for line in infile:
		line = line.strip()
		if len(line) > 0:
			if line[0] == ">":
				header = line[1:len(line)]
				seq_dict[header] = ''
			else:
				seq_dict[header] += line
	for header in seq_dict:
		seqout = seq_dict[header].replace("a","A").replace("t","T").replace("c","C").replace("g","G")
		seq_dict[header] = seqout
	return seq_dict

def retrieve_seq_from_fasta(alignment_filepath,queryID):
	seqout = ''
	infile = open(alignment_filepath,"r")
	for line in infile:
		line = line.strip()
		if len(line) > 0:
			if line[0] == ">":
				header = line[1:len(line)]
			elif header == queryID:
				seqout += line
	return seqout

def hamming_dist(seq1,seq2):
	mismatches = 0
	matches = 0
	for site in range(0,len(seq1)):
		nt1 = seq1[site]
		nt2 = seq2[site]
		if nt1 != nt2 and nt1 != "-":
			mismatches += 1
		elif nt1 == nt2:
			matches += 1
	return mismatches,matches

def add_val_to_ranges(range_list_in,subject_tup):
	range_list_in = sorted(range_list_in)
	new_range_list = []
	sstart = subject_tup[0]
	sstop = subject_tup[1]
	percident = subject_tup[2]
	merge_formed = False
	#hit     |--------------|
	#A   ----|--------------|---
	#B   ----|----------    |
	#C       |   -----------|---
	#D   --- |              |
	#E       |              |  ----
	for i in range(0,len(range_list_in)):
		hit_start = range_list_in[i][0]
		hit_stop = range_list_in[i][1]
		hit_pid = range_list_in[i][2]
		range_len = hit_stop-hit_start
		if sstart <= hit_start and sstop >= hit_stop:					    #A   ----|--------------|---
			new_hit_tup = (sstart,sstop,percident)
			new_range_list.append(new_hit_tup)
			merge_formed = True
			break
		elif sstart <= hit_start and sstop < hit_stop and sstop > hit_start:    #B   ----|----------    |
			new_range_len = hit_stop-sstart
			prev_range = hit_stop-hit_start
			s_range_len = sstop-sstart
			overlap_range_len = sstop-hit_start

			overlap_pid = ((hit_pid+percident)/2)*(overlap_range_len/new_range_len)
			prev_pid = round(hit_pid*((prev_range-overlap_range_len)/new_range_len),4)
			add_pid = round(percident*((s_range_len-overlap_range_len)/new_range_len),4)
			new_pid = prev_pid+overlap_pid+add_pid
			
			new_hit_tup = (sstart,hit_stop,new_pid)
			new_range_list.append(new_hit_tup)
			merge_formed = True
			break
		elif sstart >= hit_start and sstop > hit_stop and sstart < hit_stop:    #C       |   -----------|---
			new_range_len = sstop-hit_start
			prev_range = hit_stop-hit_start
			s_range_len = sstop-sstart
			overlap_range_len = hit_stop-sstart

			overlap_pid = ((hit_pid+percident)/2)*(overlap_range_len/new_range_len)
			prev_pid = round(hit_pid*((prev_range-overlap_range_len)/new_range_len),4)
			add_pid = round(percident*((s_range_len-overlap_range_len)/new_range_len),4)
			new_pid = prev_pid+overlap_pid+add_pid
			
			new_hit_tup = (hit_start,sstop,new_pid)
			new_range_list.append(new_hit_tup)
			merge_formed = True
			break
		else:																	#D   --- |              |
			new_range_list.append(subject_tup)									#E       |              |  ----

	return new_range_list,merge_formed

def merge_ranges(range_list_in):
	if len(range_list_in) == 1:
		return range_list_in,0
	else:
		merged_ranges = []
		merges = 0
		range_list_in = sorted(range_list_in)
		for i in range(0,len(range_list_in)):
			hit = range_list_in[i]
			start = hit[0]
			stop = hit[1]
			pid = hit[2]
			if len(merged_ranges) == 0:
				merged_ranges.append((start,stop,pid))
			else:
				merged_ranges,merge_formed = add_val_to_ranges(merged_ranges,hit)
				if merge_formed == True:
					merges+=1
		return merged_ranges,merges

def process_one_blast_file(blastfilepath):
	min_prop_minlength = 0.8
	min_avg_pid = 0.95
	min_local_pid = 0.8
	match = {}
	len_dict = {}
	blastinfile = open(blastfilepath,"r")
	for line in blastinfile:
		line = line.strip()
		split = line.split("\t")
		'''
		   0      1      2      3    4    5      6     7      8    9     10
		qseqid sseqid pident evalue qlen slen length qstart qend sstart send
		'''
		qseqid = split[0]
		sseqid = split[1]
		pair = qseqid+"\t"+sseqid
		pident = float(split[2]) #% identity
		qlen = int(split[4]) #length of query
		slen = int(split[5]) #length of subject 
		length = int(split[6]) #length of alignment
		qstart = int(split[7]) #align start site in query
		qstop = int(split[8]) #align stop site in query
		if int(qstart) > int(qstop):
			x = int(qstop)
			y = int(qstart)
			qstart = str(x)
			qstop = str(y)
		sstart = int(split[9]) #align start site in subject
		sstop = int(split[10]) #align stop site in subject
		if sstart > sstop:
			x = sstop
			y = sstart
			sstart = x
			sstop = y
		
		percident = round(pident/100.0,4)
		qrange = (qstart,qstop,percident)
		srange = (sstart,sstop,percident)		
		if percident >= min_local_pid:
			merge_formed = False
			hits_found = False
			try:
				hits = match[pair]
				hits_found = True
			except:
				match[pair] = []
				match[pair].append(srange)
			if hits_found == True:
				try:
					merged_hits,merge_formed = add_val_to_ranges(hits,srange)
				except:
					print(hits)
					print(srange)
					print(srange)
					sys.exit()
				match[pair] = merged_hits
		try:
			len_dict[sseqid]
		except:
			len_dict[sseqid] = int(slen)
		try:
			len_dict[qseqid]
		except:
			len_dict[qseqid] = int(qlen)
	blastinfile.close()

	for pair in match:
		hits = match[pair]
		qseqid = pair.split("\t")[0]
		sseqid = pair.split("\t")[1]
		qlen = len_dict[qseqid]
		slen = len_dict[sseqid]
		rangesum = 0
		for i in range(0,len(hits)):
			hit = hits[i]
			hit_start = hit[0]
			hit_stop = hit[1]
			hit_pid = hit[2]
			rangesum += hit_stop-hit_start
		avg_pid = 0
		for i in range(0,len(hits)):
			hit = hits[i]
			hit_start = hit[0]
			hit_stop = hit[1]
			hit_pid = hit[2]
			range_len_prop = (hit_stop-hit_start)/rangesum
			avg_pid += range_len_prop*hit_pid
		prop_query_length = round(float(rangesum)/float(qlen),4) # %length of query seq
		# prop_minlength = round(float(rangesum)/float(min(qlen, slen)),4) # %length of min(seq)
		avg_pid = round(avg_pid,4)
		rangesum = round(rangesum,1)

		# qperclength = float(rangesum)/float(qlen) # %length of query
		# sperclength = float(rangesum)/float(slen) # %length of subject
		rank_match = []
		if prop_query_length >= min_prop_minlength and avg_pid >= min_avg_pid:
			hit_tup = (avg_pid*prop_query_length,sseqid)
			rank_match.append(hit_tup)
		rank_match = sorted(rank_match,reverse=True)
		if len(rank_match) >0:
			top_hit_subjectID = rank_match[0][1]
		else:
			top_hit_subjectID = ''
		return top_hit_subjectID

def trim_down_hmm_alignment(filepath_in,fielpath_out):
	infile = open(filepath_in,"r")
	outfile = open(fielpath_out,"w")
	outlines = ''
	for line in infile:
		line = line.rstrip('\n')
		if len(line) >= 2:
			if line[0:2] == '//':
				outlines += line+"\n\n"
				outfile.write(outlines)
				outlines = ''
			if len(line) >= 8:
				if line[0:7] == '#=GF ID':
					outlines += line+'\n'
				elif line[0:7] == '#=GC PP':
					outlines += line+'\n'
				elif line[0:7] == '#=GC RF':
					outlines += line+'\n'
				elif len(line) > 100 and line[0] != '#':
					outlines += line+'\n'
	outfile.write(outlines)
	infile.close()
	outfile.close()

def CIGAR_parse(map_position,cigar_string,quality_string,read_seq,mapping_qual): #mapping position is indexed to 1, not zero
	cigar_dict = {'M':'imperfect_match','D':'insertion_in_ref','I':'insertion_in_read','S':'soft_clip_reference','H':'soft_clip_read','=':'perfect_match','X':'mismatch'}
	max_single_deletion = 24
	read_key = ''
	ref_deletion_map = {}
	ref_insert_map = {}
	ref_site_map = []
	match_count_dict = {}
	
	cur_position = int(map_position)-2
	cur_num = ''
	cur_operator = ''
	all_valid_operators = True
	for i in range(0,len(cigar_string)):
		character = cigar_string[i]
		try:
			cigar_dict[character]
			cur_operator = character

			cur_num = int(cur_num)
			try:
				match_count_dict[cur_operator] += cur_num
			except:
				match_count_dict[cur_operator] = cur_num
			for num in range(0,cur_num):
				if cur_operator == "D": #ref base not present in read
					cur_position += 1
					ref_deletion_map[cur_position] = len(read_key)
				elif cur_operator == "I": #insertion of base in read relative to reference
					cur_position += 0
					ref_site_map.append('')
					read_key += cur_operator
					try:
						ref_insert_map[cur_position].append(len(read_key)-1)
					except:
						ref_insert_map[cur_position] = [len(read_key)-1]
				elif cur_operator == "S": #base in read but soft clipped in ref
					cur_position += 0
					ref_site_map.append('')
					read_key += cur_operator
				elif cur_operator == "H": #base in ref but soft clipped in read
					all_valid_operators = False
				elif cur_operator == "=": #perfect match
					cur_position += 1
					ref_site_map.append(cur_position)
					read_key += cur_operator
				elif cur_operator == "M": #imperfect match (N)
					cur_position += 1
					ref_site_map.append(cur_position)
					read_key += cur_operator
				elif cur_operator == "X": #mismatch
					cur_position += 1
					ref_site_map.append(cur_position)
					read_key += cur_operator
				else:
					all_valid_operators = False
			cur_num = ''
			cur_operator = ''
		except:
			cur_num += character
	try:
		deletion_count = match_count_dict['D']
	except:
		deletion_count = 0
	base_dict_fun = {}
	read_end_site = 0
	if all_valid_operators == True:
		if deletion_count <= max_single_deletion:
			for num in range(0,len(read_key)):
				operator = read_key[num]
				ref_loc = ref_site_map[num]
				if ref_loc != '':
					ref_loc = int(ref_loc)
					try:
						qual_char = quality_string[num]
					except:
						sys.exit("Index "+str(num))
					qual_score = ord(qual_char)-33
					read_nt = read_seq[num]
					read_length = len(read_seq)
					nt_prop_location = round(num/read_length,2)
					R_end_dist = read_length-num
					tup = (read_nt,qual_score,mapping_qual,read_length,num,R_end_dist,nt_prop_location)
					base_dict_fun[ref_loc] = tup
		if ref_deletion_map != {} and deletion_count <= max_single_deletion:
			for ref_loc in ref_deletion_map:
				read_nt = '-'
				qual_score = 40
				read_loc = ref_deletion_map[ref_loc]
				read_length = len(read_seq)
				R_end_dist = read_length-read_loc
				nt_prop_location = round(read_loc/read_length,2)
				tup = (read_nt,qual_score,mapping_qual,read_length,read_loc,R_end_dist,nt_prop_location)
				base_dict_fun[ref_loc] = tup
		if ref_insert_map != {} and deletion_count <= max_single_deletion:
			for ref_loc in ref_insert_map:
				insert_base_string = base_dict_fun[ref_loc][0]
				read_loc_list = sorted(ref_insert_map[ref_loc])
				read_length = len(read_seq)
				avg_qual_score = base_dict_fun[ref_loc][1]
				avg_read_loc = base_dict_fun[ref_loc][4]
				for rl in range(0,len(read_loc_list)):
					read_loc = read_loc_list[rl]
					avg_read_loc += read_loc
					try:
						qual_char = quality_string[read_loc]
					except:
						sys.exit("Error: insertion key")
					qual_score = ord(qual_char)-33
					avg_qual_score += qual_score
					read_nt = read_seq[read_loc]
					insert_base_string += read_nt
				avg_qual_score = avg_qual_score/(len(read_loc_list)+1)
				avg_read_loc = avg_read_loc/(len(read_loc_list)+1)
				avg_R_end_dist = read_length-avg_read_loc
				avg_nt_prop_location = round(avg_read_loc/read_length,2)
				tup = (insert_base_string,avg_qual_score,mapping_qual,read_length,avg_read_loc,avg_R_end_dist,avg_nt_prop_location)
				base_dict_fun[ref_loc] = tup
	return base_dict_fun,match_count_dict

def parse_map_file(sam_filepath):
	ref_map_dict = {}
	min_mapq_score = 20
	min_phred = 20
	sam_infile = open(sam_filepath,"r")
	for sam_line in sam_infile:
		if len(sam_line) > 0:
			if sam_line[0] != "@":
				sam_line = sam_line.strip().split("\t")
				map_qual = int(sam_line[4])
				read_seq = sam_line[9]
				if read_seq !="*":# and len(read_seq) >= min_read_length:#map_qual >= min_mapping_qual and 
					readID = sam_line[0]
					map_site = int(sam_line[3])
					cigar = sam_line[5] #{'M':'imperfect_match','D':'insertion_in_ref','I':'insertion_in_read','S':'soft_clip_reference','H':'soft_clip_read','=':'perfect_match','X':'mismatch'}
					read_seq_qual = sam_line[10]
					if map_qual >= min_mapq_score:
						read_base_dict,CIGAR_operator_count_dict = CIGAR_parse(map_site,cigar,read_seq_qual,read_seq,map_qual)
						if len(read_base_dict) > 0: #per base in read: (read_nt,qual_score,mapping_qual,read_length,num,R_end_dist,nt_prop_location)
							for loc in read_base_dict:
								base = read_base_dict[loc][0]
								base_qual = read_base_dict[loc][1]
								read_len = read_base_dict[loc][3]
								loc_in_read = read_base_dict[loc][4]
								loc_from_R_end = read_base_dict[loc][5]
								prop_loc = read_base_dict[loc][6]
								if base_qual >= min_phred:
									try:
										ref_map_dict[loc][base]+=1
									except:
										try:
											ref_map_dict[loc][base]=1
										except:
											ref_map_dict[loc] = {}
											ref_map_dict[loc][base]=1
	sam_infile.close()
	return ref_map_dict

def recall_consensus(allele_dict_in,ref_seq_in,sampleID,segmentID):
	sampleID_segmentID = sampleID+"_"+segmentID
	min_major_prop=0.4
	min_deletion_prop=0.7
	min_insertion_prop=min_deletion_prop
	min_deletion_count=100
	min_insertion_count=min_deletion_count
	mapinfo_string = ''

	seq_out = ''
	loc_list = sorted(list(allele_dict_in.keys()))
	max_loc = max(loc_list)
	for loc in range(0,max_loc):
		ref_allele = ref_seq_in[loc]
		ref_count = 10
		local_allele_list = []
		try:
			local_allele_dict = allele_dict_in[loc]
		except:
			local_allele_dict = {}
		allele_sum = 0
		if local_allele_dict != {}:
			try:
				local_allele_dict[ref_allele] += ref_count
			except:
				local_allele_dict[ref_allele] = ref_count
			for allele in local_allele_dict:
				count = local_allele_dict[allele]
				allele_sum += count
				tup = (count,allele)
				local_allele_list.append(tup)
		else:
			ref_tup = (ref_count,ref_allele)
			local_allele_list.append(ref_tup)
			allele_sum += ref_count

		mapinfo_string += str(loc)+'\t'+str(allele_sum)
		local_allele_list = sorted(local_allele_list,reverse=True)
		temp_allele_dict = {}
		allele_freq_string = ''
		for i in range(0,len(local_allele_list)):
			n_allele = local_allele_list[i][1]
			n_allele_count = local_allele_list[i][0]
			n_allele_freq = round(float(n_allele_count)/float(allele_sum),3)
			temp_allele_dict[n_allele] = n_allele_freq
			allele_freq_string+= '\t'+n_allele+'='+str(n_allele_count)

		major_allele = local_allele_list[0][1]
		major_freq = local_allele_list[0][0]/allele_sum
		if len(local_allele_list) == 1:
			if major_allele != '-':
				seq_out += major_allele
				mapinfo_string += '\t'+major_allele

		elif len(local_allele_list) >1:
			minor_allele = local_allele_list[1][1]
			minor_freq = local_allele_list[1][0]/allele_sum
			if major_freq >= min_major_prop:
				if major_allele != '-' and len(major_allele) == 1:
					seq_out += major_allele
					mapinfo_string += '\t'+major_allele
				elif major_allele != '-' and len(major_allele) > 1 and major_freq>=min_insertion_prop and local_allele_list[0][0] >= min_insertion_count:
					seq_out += major_allele
					mapinfo_string += '\t'+major_allele
				elif major_allele == '-':
					if major_freq >= min_deletion_prop and local_allele_list[0][0] >= min_deletion_count:
						mapinfo_string += '\t'+'DEL'
					else:
						if minor_allele != '-' and (local_allele_list[1][0]/(allele_sum-local_allele_list[0][0])) >= min_major_prop and len(minor_allele) == 1:
							seq_out += minor_allele
							mapinfo_string += '\t'+minor_allele
						elif minor_allele != '-' and (local_allele_list[1][0]/(allele_sum-local_allele_list[0][0])) >= max(min_insertion_prop,min_major_prop) and len(minor_allele) > 1 and local_allele_list[1][0] >= min_insertion_count:
							seq_out += minor_allele
							mapinfo_string += '\t'+minor_allele
						else:
							seq_out += ref_allele
							mapinfo_string += '\t'+ref_allele

			else:
				if major_allele != '-' and minor_allele != '-':
					seq_out += ref_allele
					mapinfo_string += '\t'+ref_allele
				elif major_allele != '-' and minor_allele == '-':
					if (local_allele_list[0][0]/(allele_sum-local_allele_list[1][0])) >= min_major_prop and len(major_allele)==1:
						seq_out += major_allele
						mapinfo_string += '\t'+major_allele
					elif (local_allele_list[0][0]/(allele_sum-local_allele_list[1][0])) >= max(min_insertion_prop,min_major_prop) and len(major_allele)>1 and local_allele_list[0][0] >= min_insertion_count:
						seq_out += major_allele
						mapinfo_string += '\t'+major_allele
					else:
						seq_out += ref_allele
						mapinfo_string += '\t'+ref_allele
				elif major_allele == '-' and minor_allele != '-':
					if (local_allele_list[1][0]/(allele_sum-local_allele_list[0][0])) >= min_major_prop and len(minor_allele) == 1:
						seq_out += minor_allele
						mapinfo_string += '\t'+minor_allele
					elif (local_allele_list[1][0]/(allele_sum-local_allele_list[0][0])) >= max(min_insertion_prop,min_major_prop) and len(minor_allele) > 1 and local_allele_list[1][0] >= min_insertion_count:
						seq_out += minor_allele
						mapinfo_string += '\t'+minor_allele
					else:
						seq_out += ref_allele
						mapinfo_string += '\t'+ref_allele
				else:
					seq_out += ref_allele
					mapinfo_string += '\t'+ref_allele

		mapinfo_string += allele_freq_string+'\n'
	mapinfo_outfile = open(temp_dir+sampleID_segmentID+".refine_info.txt",'w')
	mapinfo_outfile.write(mapinfo_string)
	mapinfo_outfile.close()
	os.rename(temp_dir+sampleID_segmentID+".refine_info.txt",consensus_dir+segmentID+'/'+sampleID+'.'+segmentID+".refine_info.txt")
	return seq_out

def check_if_reads_split(sampleID,trim_dir):
	read_count_dict = {}
	for subdir in os.walk(trim_dir):
		temp_dir = subdir[0]+"/"
		segmentID = subdir[0].split("/")[-1]
		readcount_filelist = [f for f in os.listdir(temp_dir) if f.endswith(".read_count.txt") and f.startswith(sampleID)]
		for files in readcount_filelist:
			infile = open(temp_dir+files,"r")
			for line in infile:
				line = line.strip()
				if len(line)>0:
					try:
						count = int(line)
					except:
						count = 0
					if count >= 1000:
						read_count_dict[segmentID] = count
	if read_count_dict == {}:
		segment_list = []
	else:
		segment_list = list(read_count_dict.keys())
	return segment_list,read_count_dict


def refine_consensus_from_SAM(sampleID,segment_list,primary_consensus_dict,sam_dir,consensus_dir,temp_dir):
	new_consensus_dict = {}
	for s in range(0,len(segment_list)):
		segmentID = segment_list[s]
		sam_filepath = sam_dir+segmentID+"/"+sampleID+"."+segmentID+".sam"
		if os.path.isfile(sam_filepath) == True:
			allele_map_dict = parse_map_file(sam_filepath)
			ref_seq = primary_consensus_dict[segmentID]
			new_consensus = recall_consensus(allele_map_dict,ref_seq,sampleID,segmentID)
			new_consensus_dict[segmentID] = new_consensus
			print(sampleID+" "+segmentID+" ref:"+str(len(ref_seq))+" refined:"+str(len(new_consensus)))
			temp_fasta_path = temp_dir+sampleID+"."+segmentID+'.consensus_refine.fa'
			permanent_fasta_path = consensus_dir+segmentID+"/"+sampleID+"."+segmentID+'.refine.fa'
			temp_fasta = open(temp_fasta_path,"w")
			temp_fasta.write(">"+sampleID+"_"+segmentID+'\n'+new_consensus+'\n')
			temp_fasta.close()
			os.rename(temp_fasta_path,permanent_fasta_path)
	if new_consensus_dict != {}:
		temp_fasta_path = temp_dir+sampleID+'.consensus_refine.fa'
		permanent_fasta_path = consensus_dir+'all_segments/'+sampleID+'.refine.fa'
		temp_fasta = open(temp_fasta_path,"w")
		for segmentID in new_consensus_dict:
			temp_fasta.write(">"+sampleID+"_"+segmentID+'\n'+new_consensus_dict[segmentID]+'\n')
		temp_fasta.close()
		os.rename(temp_fasta_path,permanent_fasta_path)
	return new_consensus_dict


def hmm_align_reads(sampleID,hmm_dir,trim_dir,consensus_dir,nhmmer_dir,ref_dir,temp_dir,nhmmer_options,read_input_type,hmm_target="influenza"):
	alignment_filepath = nhmmer_dir+sampleID+'.nhmmer_align.txt'
	temp_alignment_filepath = temp_dir+sampleID+'.nhmmer_align.txt'
	tblout_filepath = nhmmer_dir+sampleID+'.nhmmer_tblout.txt'
	temp_tblout_filepath = temp_dir+sampleID+'.nhmmer_tblout.txt'
	fastq_in_filepath = trim_dir+sampleID+".fq"
	fasta_temp_filepath = temp_dir+sampleID+".fasta"
	if os.path.isfile(tblout_filepath) == False:
		if os.path.isfile(temp_alignment_filepath) == True:
			os.remove(temp_alignment_filepath)
			try:
				os.remove(temp_tblout_filepath)
			except:
				pass
		if os.path.isfile(temp_alignment_filepath) == False:
			fastq_to_fasta(fastq_in_filepath,fasta_temp_filepath)
			command = 'nhmmer '+nhmmer_options+' -A '+temp_alignment_filepath+' --tblout '+temp_tblout_filepath+' '+hmm_dir+hmm_target+'.hmm '+fasta_temp_filepath
			run_status = run_command(command)
			if run_status == False:
				sys.exit("nhmmer run failed: "+command)
			print(sampleID+" - nhmmer done")
		
		trim_down_hmm_alignment(temp_alignment_filepath,temp_alignment_filepath.replace("nhmmer_align","nhmmer_align_short"))
		os.rename(temp_alignment_filepath.replace("nhmmer_align","nhmmer_align_short"),alignment_filepath)
		os.rename(temp_tblout_filepath,tblout_filepath)
		os.remove(temp_alignment_filepath)
		os.remove(fasta_temp_filepath)
		print(sampleID+" - trim down done")
	

	permanent_full_seqout_filename = consensus_dir+"all_segments/"+sampleID+".fa"
	seq_out_dict = {}
	segment_list = []
	if os.path.isfile(permanent_full_seqout_filename) == True:
		temp_seq_out_dict = load_alignment(permanent_full_seqout_filename)
		for key in temp_seq_out_dict:
			segmentID = key.replace(sampleID+'_','')
			segment_list.append(segmentID)
			seq_out_dict[segmentID] = temp_seq_out_dict[key]
		del temp_seq_out_dict
	elif os.path.isfile(permanent_full_seqout_filename) == False:
		rough_consensus_fasta_path = consensus_dir+"all_segments/"+sampleID+'.rough_consensus.fa'
		split_segment_list,read_count_dict = check_if_reads_split(sampleID,trim_dir)
		if read_count_dict == {} or os.path.isfile(rough_consensus_fasta_path) == False:
			readID_to_segment_dict,segment_list = parse_nhmmer_tblout(tblout_filepath)
			print(sampleID+" - tblout parse done - "+str(len(readID_to_segment_dict))+" - "+" ".join(segment_list))
			if os.path.isfile(trim_dir+"remain/"+sampleID+'.nhmmer_remain.read_count.txt') == False:
				read_count_dict = split_fastq_by_segment(sampleID,segment_list,readID_to_segment_dict,temp_dir,trim_dir,ref_dir,nhmmer_dir,read_input_type)
				print(sampleID+" - read split done - "+str(read_count_dict))
			if os.path.isfile(rough_consensus_fasta_path) == False:
				nhmmer_seq_dict = parse_nhmmer_alignment(sampleID,alignment_filepath,readID_to_segment_dict)
				print(sampleID+" - alignment parse done - "+" ".join(list(nhmmer_seq_dict.keys())))
			else:
				nhmmer_seq_dict = load_alignment(rough_consensus_fasta_path)
				print(sampleID+" - loaded rough consensus alignment - "+" ".join(list(nhmmer_seq_dict.keys())))
			del readID_to_segment_dict
		else:
			nhmmer_seq_dict = load_alignment(rough_consensus_fasta_path)
			segment_list = list(nhmmer_seq_dict.keys())
		for segmentID in nhmmer_seq_dict:
			rough_consensus_seq = nhmmer_seq_dict[segmentID]
			if rough_consensus_seq != '':
				clean_consensus_filename = consensus_dir+segmentID+"/"+sampleID+"."+segmentID+".fa"
				if not os.path.exists(consensus_dir+segmentID):
					os.makedirs(consensus_dir+segmentID)
				if os.path.isfile(clean_consensus_filename) == False:
					if len(rough_consensus_seq.replace(" ","").replace("-","")) > 400:
						rough_consensus_filename = temp_dir+sampleID+"."+segmentID+".consensus.fa"
						rough_consensus_align_filename = temp_dir+sampleID+"."+segmentID+".consensus.afa"
						seqout_file = open(rough_consensus_filename,"w")
						seqout_file.write(">"+sampleID+"_"+segmentID+"\n"+rough_consensus_seq.replace('-','')+"\n")
						seqout_file.close()

						full_align_path = hmm_dir+"full_influenza_alignments/"+segmentID+".fasta"
						full_align_blast_db_path = hmm_dir+"full_influenza_alignments/blast_db/"+segmentID
						blast_outfile_name = temp_dir+sampleID+'.'+segmentID+'.blastn.txt'
						command = "blastn -query "+rough_consensus_filename+" -db "+full_align_blast_db_path+' -evalue 1E-10 -outfmt "6 qseqid sseqid pident evalue qlen slen length qstart qend sstart send" -out '+blast_outfile_name
						run_status = run_command(command)
						if run_status == False:
							sys.exit("blastn failed: "+command)
						
						hmm_refseq_dict = load_alignment(hmm_dir+'alignments/'+segmentID+'.fasta')
						seqout_file = open(rough_consensus_filename,"a")
						for ref_seqID in hmm_refseq_dict:
							seqout_file.write(">"+ref_seqID+"\n"+hmm_refseq_dict[ref_seqID].replace('-','')+"\n")
						del hmm_refseq_dict
						top_hit_ref_seqID = process_one_blast_file(blast_outfile_name)
						if top_hit_ref_seqID != '':
							top_hit_seq = retrieve_seq_from_fasta(full_align_path,top_hit_ref_seqID)
							if top_hit_seq != '' and top_hit_seq != None:
								seqout_file.write(">"+top_hit_ref_seqID+"\n"+top_hit_seq.replace('-','')+"\n")
						seqout_file.close()
						os.remove(blast_outfile_name)

						command = 'mafft --op 1.75 --quiet --thread 4 '+rough_consensus_filename+' > '+rough_consensus_align_filename
						run_status = run_command(command,"verbose")
						if run_status == False:
							sys.exit("mafft failed: "+command)
						os.remove(rough_consensus_filename)

						hmm_align_dict = load_alignment(rough_consensus_align_filename)
						aligned_rough_seq = hmm_align_dict[sampleID+"_"+segmentID]

						dist_list = []
						for ref_seqID in hmm_align_dict:
							if ref_seqID != sampleID+"_"+segmentID:
								ref_seq_string = hmm_align_dict[ref_seqID]
								mismatches,matches = hamming_dist(aligned_rough_seq,ref_seq_string)
								dist_tup = (mismatches,1.0/matches,ref_seqID)
								dist_list.append(dist_tup)
						dist_list = sorted(dist_list,reverse=False)
						nearest_refseqID_dist = dist_list[0][0]
						nearest_refseqID_matches = int(1.0/dist_list[0][1])
						nearest_refseqID = dist_list[0][2]
						nearest_refseq = hmm_align_dict[nearest_refseqID]
						print(sampleID+" "+segmentID+" "+nearest_refseqID+" "+str(len(nearest_refseq))+" "+str(nearest_refseqID_matches)+" "+str(nearest_refseqID_dist))
						cleaned_consensus = ''
						for loc in range(0,len(nearest_refseq)):
							if aligned_rough_seq[loc] == '-' and nearest_refseq[loc] != '-':
								cleaned_consensus += nearest_refseq[loc]
							elif aligned_rough_seq[loc] != '-':
								cleaned_consensus += aligned_rough_seq[loc]

						temp_clean_consensus_filename = temp_dir+sampleID+"."+segmentID+".consensus.fa"
						
						
						seqout_file = open(temp_clean_consensus_filename,"w")
						seqout_file.write(">"+sampleID+"_"+segmentID+"\n"+cleaned_consensus+"\n")
						seqout_file.close()
						
						seq_out_dict[segmentID] = cleaned_consensus
						os.remove(rough_consensus_align_filename)
						os.rename(temp_clean_consensus_filename,clean_consensus_filename)
				else:
					cleaned_consensus = ''
					temp_infile = open(clean_consensus_filename,"r")
					for line in temp_infile:
						line = line.rstrip("\n")
						if line[0]!=">":
							cleaned_consensus += line
					temp_infile.close()
					seq_out_dict[segmentID] = cleaned_consensus
		seqout_lines = ''
		if seq_out_dict != {}:
			for segmentID in seq_out_dict:
				seqout_lines += ">"+sampleID+"_"+segmentID+"\n"+seq_out_dict[segmentID]+"\n"
		temp_full_seqout_filename = temp_dir+sampleID+".all_consensus.fa"
		
		seqout_file = open(temp_full_seqout_filename,"w")
		seqout_file.write(seqout_lines)
		seqout_file.close()
		os.rename(temp_full_seqout_filename,permanent_full_seqout_filename)
	# print(sampleID+" - complete")
	return seq_out_dict,segment_list

def map_reads_to_ref(sampleID,segment_list,trim_dir,sam_dir,consensus_dir,read_count_dict,ref_suffix='',overwrite_existing_files=False):
	'''
	After the fastq file has been split into reads that map to each segment,
	map one set of reads to the reference sequence(s) provided in a reference set
	'''
	for segmentID in segment_list:
		try:
			read_num = read_count_dict[segmentID]
		except:
			read_num = 0
		if read_num >= 1000:
			split_fastq = trim_dir+segmentID+"/"+sampleID+"."+segmentID+".fq"
			sample_consensus_ref = consensus_dir+segmentID+"/"+sampleID+"."+segmentID+ref_suffix+".fa"
			### Map reads to reference
			input_filename_1 = split_fastq
			output_filename_1 = temp_dir+sampleID+"."+segmentID+ref_suffix+".sam"
			final_output_filename = sam_dir+segmentID+"/"+sampleID+"."+segmentID+ref_suffix+".sam"
			continue_map = True
			if os.path.isfile(input_filename_1) == False or os.path.isfile(sample_consensus_ref) == False:
				continue_map = False
			if os.path.isfile(final_output_filename) == True and continue_map == True:
				if overwrite_existing_files == True:
					os.remove(final_output_filename)
				else:
					continue_map = False
			if os.path.isfile(output_filename_1) == True and continue_map == True:
				os.remove(output_filename_1)
			if os.path.isfile(input_filename_1) == True and continue_map == True:
				print("Mapping processed reads - "+segmentID+' '+sampleID)
				command = command_prefix+"bbmap.sh in="+input_filename_1+" outm="+output_filename_1+" ref="+sample_consensus_ref+" "+bbmap_options
				command_status = run_command(command)
				if command_status == False:
					sys.exit("bbmap failed: "+input_filename_1)
				
				if os.path.isfile(output_filename_1) == True:
					### Filter out unmapped reads from the .bam file
					input_filename_1 = output_filename_1
					output_filename_1 = temp_dir+sampleID+"."+segmentID+ref_suffix+".bam"
					if os.path.isfile(output_filename_1) == True:
						os.remove(output_filename_1)
					command = "samtools view -b -F 4 "+input_filename_1+" > "+output_filename_1
					command_status = os.system(command)
					if command_status != 0:
						sys.exit("samtools view failed: "+input_filename_1)

					if os.path.isfile(input_filename_1) == True:
						os.remove(input_filename_1)
					if os.path.isfile(output_filename_1) == True:
						### Convert BAM file into SAM file and sort reads according to their location in the segment
						input_filename_1 = output_filename_1
						output_filename_1 = temp_dir+sampleID+"."+segmentID+ref_suffix+".sam"
						if os.path.isfile(output_filename_1) == True:
							os.remove(output_filename_1)
						command = "samtools sort -O sam -o "+output_filename_1+" "+input_filename_1
						command_status = os.system(command)
						if command_status != 0:
							sys.exit("samtools view failed: "+input_filename_1)
						os.remove(input_filename_1)

						### Files are stored separately according to the reference sequence they were aligned to
						### If a reference-specific storage directory doesn't exist, make it
						if not os.path.exists(sam_dir+segmentID):
							os.makedirs(sam_dir+segmentID)
						
						### Move the sorted SAM file to the storage directory
						input_filename_1 = output_filename_1
						output_filename_1 = sam_dir+segmentID+"/"+sampleID+"."+segmentID+ref_suffix+".sam"
						readcount_filename = sam_dir+segmentID+"/"+sampleID+"."+segmentID+ref_suffix+".read_count.txt"
						if os.path.isfile(output_filename_1) == True:
							os.remove(output_filename_1)
						if os.path.isfile(readcount_filename) == True:
							os.remove(readcount_filename)
						try:
							os.rename(input_filename_1,output_filename_1)
						except:
							print("no SAM output for: "+sampleID+" - "+segmentID+" - "+str(read_num))



###########################################################


def main_pipeline(sampleID,input_read_dir,ref_dir,trim_dir,map_dir,qual_dir,sam_dir,bam_dir,temp_dir,sub_dir,hmm_dir,consensus_dir,nhmmer_dir,forward_read_suffix,reverse_read_suffix):
	print(sampleID+' - starting')
	if raw_read_input_type == "mixed":
		local_raw_read_input_type = determine_raw_read_input_type(sampleID,input_read_dir,forward_read_suffix,reverse_read_suffix)
	else:
		local_raw_read_input_type = raw_read_input_type

	input_reads_found = check_if_input_reads_exist(sampleID,input_read_dir,local_raw_read_input_type,forward_read_suffix,reverse_read_suffix)
	trimmed_reads_found = check_trimmed_reads_exist(sampleID,trim_dir)
	if input_reads_found == True and trimmed_reads_found == False:
		if local_raw_read_input_type == "paired":
			trimmed_fastq = raw_read_processing_paired(sampleID,input_read_dir,ref_dir,trim_dir,temp_dir,forward_read_suffix,reverse_read_suffix,command_prefix)
		elif local_raw_read_input_type == "single":
			trimmed_fastq = raw_read_processing_single(sampleID,input_read_dir,ref_dir,trim_dir,temp_dir,forward_read_suffix,command_prefix)
		else:
			sys.exit('Invaild read input type. Only valid options are: "paired" or "single" - Exiting.')
	permanent_full_seqout_filename = consensus_dir+"all_segments/"+sampleID+".fa"
	if os.path.isfile(permanent_full_seqout_filename) == True:
		consensus_dict = {}
		segment_list = []
		temp_consensus_dict = load_alignment(permanent_full_seqout_filename)
		temp_segment_list = list(temp_consensus_dict.keys())
		for item in temp_segment_list:
			segmentID = item.replace(sampleID+'_','')
			segment_list.append(segmentID)
			consensus_dict[segmentID] = temp_consensus_dict[item]
	else:
		consensus_dict,segment_list = hmm_align_reads(sampleID,hmm_dir,trim_dir,consensus_dir,nhmmer_dir,ref_dir,temp_dir,nhmmer_options,local_raw_read_input_type)
	read_count_dict = load_read_counts(sampleID,segment_list,trim_dir)
	if read_count_dict == {}:
		return None
	
	reads_mapped,missing_ref_sets = check_if_reads_mapped(sampleID,segment_list,sam_dir,read_count_dict)
	ref_refined = check_if_ref_refined(sampleID,consensus_dir,read_count_dict)
	refined_reads_mapped,refine_missing_ref_sets = check_if_reads_mapped(sampleID,segment_list,sam_dir,read_count_dict,'.refine')
	if reads_mapped == False and ref_refined == False and refined_reads_mapped == False:
		print(sampleID+" - map reads")#- "+str(reads_mapped)+' '+str(ref_refined)+' '+str(refined_reads_mapped))
		map_reads_to_ref(sampleID,missing_ref_sets,trim_dir,sam_dir,consensus_dir,read_count_dict)

	reads_mapped,missing_ref_sets = check_if_reads_mapped(sampleID,segment_list,sam_dir,read_count_dict)
	ref_refined = check_if_ref_refined(sampleID,consensus_dir,read_count_dict)
	# refined_reads_mapped,refine_missing_ref_sets = check_if_reads_mapped(sampleID,segment_list,sam_dir,read_count_dict,'.refine')
	if reads_mapped == True and ref_refined == False:
		print(sampleID+" - refine consensus")#- "+str(reads_mapped)+' '+str(ref_refined)+' '+str(refined_reads_mapped))
		print(segment_list)
		refine_consensus_dict = refine_consensus_from_SAM(sampleID,segment_list,consensus_dict,sam_dir,consensus_dir,temp_dir)

	# reads_mapped,missing_ref_sets = check_if_reads_mapped(sampleID,segment_list,sam_dir,read_count_dict)
	ref_refined = check_if_ref_refined(sampleID,consensus_dir,read_count_dict)
	refined_reads_mapped,refine_missing_ref_sets = check_if_reads_mapped(sampleID,segment_list,sam_dir,read_count_dict,'.refine')
	if ref_refined == True and refined_reads_mapped == False:
		print(sampleID+" - remap reads")# - "+str(reads_mapped)+' '+str(ref_refined)+' '+str(refined_reads_mapped))
		map_reads_to_ref(sampleID,refine_missing_ref_sets,trim_dir,sam_dir,consensus_dir,read_count_dict,'.refine')

	reads_mapped,missing_ref_sets = check_if_reads_mapped(sampleID,segment_list,sam_dir,read_count_dict)
	ref_refined = check_if_ref_refined(sampleID,consensus_dir,read_count_dict)
	refined_reads_mapped,refine_missing_ref_sets = check_if_reads_mapped(sampleID,segment_list,sam_dir,read_count_dict,'.refine')

	print(sampleID+' - complete!')
	return sampleID



####################### MAIN ##############################


if not os.path.exists(trim_dir):
	os.makedirs(trim_dir)
if not os.path.exists(trim_dir+'remain'):
	os.makedirs(trim_dir+'remain')
if not os.path.exists(map_dir):
	os.makedirs(map_dir)
if not os.path.exists(qual_dir):
	os.makedirs(qual_dir)
if not os.path.exists(sam_dir):
	os.makedirs(sam_dir)
if not os.path.exists(sub_dir):
	os.makedirs(sub_dir)
if not os.path.exists(consensus_dir):
	os.makedirs(consensus_dir)
if not os.path.exists(consensus_dir+"all_segments"):
	os.makedirs(consensus_dir+"all_segments")
if not os.path.exists(nhmmer_dir):
	os.makedirs(nhmmer_dir)
if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)


## Create list of samples to process
file_list = [f for f in os.listdir(input_read_dir) if f.endswith(forward_read_suffix)]
sampleID_list = []
print_space_offset = 0
for files in file_list:
	sampleID = files.split(forward_read_suffix)[0]
	sampleID_list.append(sampleID)
sampleID_list = list(set(sampleID_list))

if len(sampleID_list)<num_cores:
	num_cores = len(sampleID_list)

print(str(len(sampleID_list))+" total items to process, "+str(num_cores)+" cores requested, "+str(nhmmer_cores_per_job)+" cores each job has access to")

processed_list = Parallel(n_jobs=num_cores)(delayed(main_pipeline)(sampleID,input_read_dir,ref_dir,trim_dir,map_dir,qual_dir,sam_dir,bam_dir,temp_dir,sub_dir,hmm_dir,consensus_dir,nhmmer_dir,forward_read_suffix,reverse_read_suffix) for sampleID in sampleID_list)
