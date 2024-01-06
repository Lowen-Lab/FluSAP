import os
import sys
import math
import random
import numpy as np
import time


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


verbose = True
overwrite_existing_files = False #when False, pipeline will avoid re-running samples that have already been processed 
reads_already_screened_for_quality = False # 'True' or 'False' --- If reads were already merged and screened for average quality and length (such as during a demultiplexing step), the script will not attempt to filter or merge the input reads

parallel_process = True
parallel_max_cpu = -10


raw_read_input_type = "paired" # 'paired', 'single', or 'mixed' --- When 'single' the pipeline will look for barcodes in one file, when 'paired' the pipeline will look for forward and reverse reads
forward_read_suffix = "" #if empty, the script will attempt to predict these values, but it only works under specific assumptions
reverse_read_suffix = "" #this value is ignored if 'reads_already_screened_for_quality' is 'True' or 'raw_read_input_type' is 'single' - NOTE: this variable cannot be removed without causing missing value errors

ref_set_list = ['A_HA_H1.fa','A_HA_H3.fa','A_M_1.fa-A_M_2.fa','A_NA_N1.fa','A_NA_N2.fa','A_NP_1.fa-A_NP_2.fa','A_NS.fa','A_PA_1.fa-A_PA_2.fa','A_PB1_1.fa-A_PB1_2.fa','A_PB2_1.fa-A_PB2_2.fa','B_HA.fa','B_M.fa','B_NA.fa','B_NP.fa','B_NS.fa','B_PA.fa','B_PB1.fa','B_PB2.fa']
# ref_set_list = ['PB2.fa','PB1.fa','PA.fa','HA.fa','NP.fa','NA.fa','M.fa','NS.fa']

#BLAT Settings
min_evalue = 1.0E-4
min_bitscore = 50.0

command_prefix = "bash "
repair_options = 'fint=t repair=t overwrite=t -da -Xmx1000m'
bbmerge_options = 'qin=33 mix=t overwrite=t -Xmx1000m threads=4'
bbduk_options = 'qin=33 tpe=f tbo=f overwrite=t -Xmx1000m'
bbmap_options = 'qin=33 local=f touppercase=t overwrite=t threads=4 nodisk -Xmx1000m'


############################## Functions #############################################
def run_command(command):
	run_status = False
	# return_code = os.system(command)
	return_code = os.system(command+" >/dev/null 2>&1")
	if return_code == 0:
		run_status = True
	return run_status

def unzip_gz(filename,path_to_file):
	command = "gzip -d "+path_to_file+filename
	return_code = run_command(command)
	if return_code == False:
		sys.exit()
	return return_code

def unzip_tar_gz(filename,path_to_file):
	command = "tar -xvf "+path_to_file+filename
	return_code = run_command(command)
	if return_code == False:
		sys.exit()
	return return_code

def zip_gz(filename,path_to_file):
	command = "gzip -9 "+path_to_file+filename
	return_code = run_command(command)
	if return_code == False:
		sys.exit()
	return return_code

def zip_tar_gz(filename,path_to_file,zipped_file_extension):
	command = "tar -cvzf "+path_to_file+filename.split(zipped_file_extension)[0]+".tar.gz "+path_to_file+filename
	return_code = run_command(command)
	if return_code == False:
		sys.exit()
	return return_code

def check_if_input_zipped(input_read_dir,parallel_process,num_cores):
	files_found = False
	parallel_process = False
	#Unzip any files that need unzipping
	gz_inputs = [f for f in os.listdir(input_read_dir) if f.endswith(".gz")]
	if len(gz_inputs) >0:  #unzip .gz files
		if parallel_process == True:
			processed_list = Parallel(n_jobs=num_cores)(delayed(unzip_gz)(i,input_read_dir) for i in gz_inputs)
		else:
			for gzfile in gz_inputs:
				unzip_gz(gzfile,input_read_dir)
	tar_gz_inputs = [f for f in os.listdir(input_read_dir) if f.endswith(".tar.gz")]
	if len(tar_gz_inputs) > 0:  #unzip .tar.gz files
		if parallel_process == True:
			processed_list = Parallel(n_jobs=num_cores)(delayed(unzip_tar_gz)(i,input_read_dir) for i in tar_gz_inputs)
		else:
			for targzfile in tar_gz_inputs:
				unzip_tar_gz(gzfile,input_read_dir)
	#Find what the 
	fq_inputs = [f for f in os.listdir(input_read_dir) if f.endswith(".fastq") or f.endswith(".fq")]
	if len(fq_inputs) >0:
		files_found = True
	return files_found

def map_reads_to_ref(sampleID,map_ref_set_list,trim_dir,sam_dir,read_count_dict,overwrite_existing_files):
	'''
	After the fastq file has been split into reads that map to each segment,
	map one set of reads to the reference sequence(s) provided in a reference set
	'''
	# update_print_line("Mapping processed reads",sampleID)
	for ref_set in map_ref_set_list:
		refs = ref_set.split("-")
		primary_ref = refs[0].split(".f")[0]
		for ref_filename in refs:
			ref_ID = ref_filename.split(".f")[0]
			try:
				read_num = read_count_dict[ref_ID]
			except:
				read_num = 0
			if read_num >= 1000:
				split_fastq = trim_dir+primary_ref+"/"+sampleID+"."+ref_ID+".fq"
				### Map reads to reference
				input_filename_1 = split_fastq
				output_filename_1 = temp_dir+sampleID+"."+ref_ID+".sam"
				final_output_filename = sam_dir+ref_ID+"/"+sampleID+"."+ref_ID+".sam"
				# unmapped_read_filename = temp_dir+sampleID+"."+ref_ID+".outu.fq"
				continue_map = True
				if os.path.isfile(input_filename_1) == False:
					continue_map = False
				if os.path.isfile(final_output_filename) == True and continue_map == True:
					if overwrite_existing_files == True:
						os.remove(final_output_filename)
					else:
						continue_map = False
				if os.path.isfile(output_filename_1) == True and continue_map == True:
					os.remove(output_filename_1)
				if os.path.isfile(input_filename_1) == True and continue_map == True:
					update_print_line("Mapping processed reads - "+ref_ID,sampleID)
					command = command_prefix+"bbmap.sh in="+input_filename_1+" outm="+output_filename_1+" ref="+ref_dir+ref_filename+" "+bbmap_options# minid=0.90
					command_status = run_command(command)
					if command_status == False:
						sys.exit("bbmap failed: "+input_filename_1)
					
					if os.path.isfile(output_filename_1) == True:
						### Filter out unmapped reads from the .bam file
						input_filename_1 = output_filename_1
						output_filename_1 = temp_dir+sampleID+"."+ref_ID+".bam"
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
							output_filename_1 = temp_dir+sampleID+"."+ref_ID+".sam"
							if os.path.isfile(output_filename_1) == True:
								os.remove(output_filename_1)
							command = "samtools sort -O sam -o "+output_filename_1+" "+input_filename_1
							command_status = os.system(command)
							if command_status != 0:
								sys.exit("samtools view failed: "+input_filename_1)
							# if os.path.isfile(output_filename_1) == True:
							os.remove(input_filename_1)

							### Files are stored separately according to the reference sequence they were aligned to
							### If a reference-specific storage directory doesn't exist, make it
							if not os.path.exists(sam_dir+ref_ID):
								os.makedirs(sam_dir+ref_ID)
							
							### Move the sorted SAM file to the storage directory
							input_filename_1 = output_filename_1
							output_filename_1 = sam_dir+ref_ID+"/"+sampleID+"."+ref_ID+".sam"
							readcount_filename = sam_dir+ref_ID+"/"+sampleID+"."+ref_ID+".read_count.txt"
							if os.path.isfile(output_filename_1) == True:
								os.remove(output_filename_1)
							if os.path.isfile(readcount_filename) == True:
								os.remove(readcount_filename)
							try:
								os.rename(input_filename_1,output_filename_1)
							except:
								print("no SAM output for: "+sampleID+" - "+ref_filename+" - "+str(read_num))


def split_fastq_by_segment(sampleID,ref_set_list,temp_dir,trim_dir,ref_dir,min_evalue,min_bitscore):
	update_print_line("Splitting input reads by segment",sampleID)
	fastq_filename = trim_dir+sampleID+".fq"
	fasta_filename = temp_dir+sampleID+".fasta"

	#Convert fastq file into temporary fasta file
	out = fastq_to_fasta(fastq_filename,fasta_filename)

	#map reads roughly to references using BLAT
	blat_filename = temp_dir+sampleID+".blat"
	command = 'blat -t=dna -q=dna -minIdentity=70 -out=blast8 '+ref_dir+'all_refs.fa '+fasta_filename+' '+blat_filename
	command_out = run_command(command)
	if command_out == False:
		sys.exit("blat failed: "+sampleID)
	os.remove(fasta_filename)

	blat_dict = {}
	blat_infile = open(blat_filename,"r")
	for line in blat_infile:
		line = line.strip()
		if len(line) > 0:
			line = line.split("\t")
			# PID,align_len,mismatch,gap_open = float(line[2]),float(line[3]),float(line[4]),float(line[5])
			# q_start,q_end,s_start,s_end,evalue,bitscore = float(line[6]),float(line[7]),float(line[8]),float(line[9]),float(line[10]),float(line[11])
			query,subject,evalue,bitscore = line[0],line[1],float(line[10]),float(line[11])
			tup = (evalue,subject)
			if evalue <= min_evalue and bitscore >= min_bitscore:
				try:
					blat_dict[query].append(tup)
				except:
					blat_dict[query] = []
					blat_dict[query].append(tup)
	blat_infile.close()
	os.remove(blat_filename)

	#parse BLAT output
	for query in blat_dict:
		blat_list = blat_dict[query]
		if len(blat_list) == 1:
			blat_dict[query] = blat_list[0][1]
		else:
			blat_list = sorted(blat_list)
			# print(blat_list)
			if blat_list[0][0] < blat_list[1][0]:
				blat_dict[query] = blat_list[0][1]
			else:
				blat_dict[query] = ''

	#split reads into separate files according to BLAT output
	lineout_dict = {}
	read_count_dict = {}
	primary_ref_dict = {}
	for ref_set in ref_set_list:
		refs = ref_set.split("-")
		primary_ref = refs[0].split(".f")[0]
		for ref_filename in refs:
			current_segment = ref_filename.split(".f")[0]
			# primary_ref_dict[current_segment] = primary_ref
			if not os.path.exists(trim_dir+current_segment):
				os.makedirs(trim_dir+current_segment)
			temp_fq = temp_dir+sampleID+'.'+current_segment+'.fq'
			ref_outfile = open(temp_fq,"w")
			ref_outfile.close()
			lineout_dict[current_segment] = ''
			read_count_dict[current_segment] = 0

	fastq_infile = open(fastq_filename,"r")
	line_counter = -1
	current_segment = ''
	for line in fastq_infile:
		line = line.strip()
		line_counter += 1
		if line_counter == 0:
			query = line[1:len(line)].replace(" ","___")
			try:
				current_segment = blat_dict[query]
			except:
				current_segment = ''
		# elif line_counter == 1:
			# outfile.write(line+"\n")
		elif line_counter == 3:
			line_counter = -1
		# if current_segment == segment:
		if current_segment != '':
			lineout_dict[current_segment] += line+"\n"
			if line_counter == 0:
				read_count_dict[current_segment] += 1
			elif line_counter == -1:
				if len(lineout_dict[current_segment]) >= 100000:
				# if read_count_dict[current_segment] % 500 == 0 and read_count_dict[current_segment]>10:
					temp_fq = temp_dir+sampleID+'.'+current_segment+'.fq'
					ref_outfile = open(temp_fq,"a")
					ref_outfile.write(lineout_dict[current_segment])
					ref_outfile.close()
					lineout_dict[current_segment] = ''
	fastq_infile.close()

	for ref_set in ref_set_list:
		refs = ref_set.split("-")
		primary_ref = refs[0].split(".f")[0]
		for ref_filename in refs:
			current_segment = ref_filename.split(".f")[0]
			temp_fq = temp_dir+sampleID+'.'+current_segment+'.fq'
			
			if lineout_dict[current_segment] != '':
				ref_outfile = open(temp_fq,"a")
				ref_outfile.write(lineout_dict[current_segment])
				ref_outfile.close()
				lineout_dict[current_segment] = ''
			
			try:
				seg_read_count = read_count_dict[current_segment]
			except:
				seg_read_count = 0
			read_count_file = open(trim_dir+primary_ref+"/"+sampleID+'.'+current_segment+'.read_count.txt',"w")
			read_count_file.write(str(seg_read_count)+"\n")
			read_count_file.close()

			permanent_fq = trim_dir+primary_ref+"/"+sampleID+'.'+current_segment+'.fq'
			os.rename(temp_fq,permanent_fq)
	return read_count_dict


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
			outfile.write(line+"\n")
		elif line_counter == 3:
			line_counter = -1
	infile.close()
	outfile.close()
	return filename_out


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
					adapter_pass = True
	return adapter_pass


def check_refs(ref_set_list,ref_dir):
	### If there is not already a concatenated fasta of all references, make it now
	all_ref_file_exists = os.path.isfile(ref_dir+"all_refs.fa")
	if all_ref_file_exists == True:
		remake_refs = False
	else:
		remake_refs = True

	if os.path.isfile(ref_dir+"all_refs.fa") == False or remake_refs == True:
		concat_refs = open(ref_dir+"all_refs.fa","w")
		temp_seq_dict = {}
		for ref_set in ref_set_list:
			refs = ref_set.split("-")
			for ref_filename in refs:
				infile = open(ref_dir+ref_filename,"r")
				for line in infile:
					line = line.strip()
					if len(line) > 0:
						if line[0] == ">":
							head = line[1:len(line)]
							temp_seq_dict[head] = ''
						else:
							temp_seq_dict[head] += line
				infile.close()
		unqiue_ref_seqs = {}
		for head in temp_seq_dict:
			seq = temp_seq_dict[head]
			try:
				unqiue_ref_seqs[seq]
			except:
				unqiue_ref_seqs[seq] = head
				concat_refs.write(">"+head+"\n"+seq+"\n")
		concat_refs.close()


def update_print_line(string_in,sampleID):
	global print_space_offset
	print(sampleID+" "*(print_space_offset-len(sampleID))+" - "+string_in)

def determine_raw_read_input_type(sampleID,input_read_dir,forward_read_suffix,reverse_read_suffix):
	forward_input_reads_found = False
	reverse_input_reads_found = False
	if os.path.isfile(input_read_dir+sampleID+forward_read_suffix) == True:
		forward_input_reads_found = True
	else:
		sys.exit("Unable to find input reads for sample: "+sampleID)
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


def check_if_reads_split(sampleID,ref_set_list,trim_dir):
	found_splitfile = True
	for ref_set in ref_set_list:
		refs = ref_set.split("-")
		primary_ref = ref_set.split("-")[0].split(".f")[0]
		for ref_filename in refs:
			split_fastq_filename = trim_dir+primary_ref+"/"+sampleID+"."+ref_filename.split(".f")[0]+".fq"
			if os.path.isfile(split_fastq_filename) == False:
				found_splitfile = False
				# print("missing: "+split_fastq_filename)
	return found_splitfile

def check_if_reads_mapped(sampleID,ref_set_list,sam_dir,read_count_dict):
	reads_mapped = True
	missing_ref_sets = []
	for ref_set in ref_set_list:
		refs = ref_set.split("-")
		primary_ref = ref_set.split("-")[0].split(".f")[0]
		for ref_filename in refs:
			ref_ID = ref_filename.split(".f")[0]
			if os.path.isfile(sam_dir+primary_ref+"/"+sampleID+"."+ref_ID+'.sam') == False:
				if read_count_dict[ref_ID] >= 1000:
					reads_mapped = False
					missing_ref_sets.append(ref_set)
	missing_ref_sets = list(set(missing_ref_sets))
	return reads_mapped,missing_ref_sets


def load_read_counts(sampleID,ref_set_list,trim_dir):
	read_count_dict = {}
	for ref_set in ref_set_list:
		refs = ref_set.split("-")
		primary_ref = ref_set.split("-")[0].split(".f")[0]
		for ref_filename in refs:
			ref_ID = ref_filename.split(".f")[0]
			try:
				read_count_file = open(trim_dir+primary_ref+"/"+sampleID+'.'+ref_ID+'.read_count.txt','r')
				for line in read_count_file:
					line = line.strip()
					read_count_dict[ref_ID] = int(line)
				read_count_file.close()
			except:
				read_count_dict[ref_ID] = 0
	return read_count_dict


def predict_paired_file_extensions(input_read_dir):
	exten_filelist = [f for f in os.listdir(input_read_dir) if f.endswith(".fastq") or f.endswith(".fq")]
	char_dict = {}
	for i in range(-1,(-1*len(exten_filelist[0])),-1):
		char_dict[i] = {}
		for f in range(0,len(exten_filelist)):
			filename = exten_filelist[f]
			try:
				char = filename[i]
			except:
				char = "na"
			try:
				char_dict[i][char] += 1
			except:
				char_dict[i][char] = 1
	first_dichotomous = 0
	dichotomous_sites = []
	last_monomorphic = 0
	searching_for_last_monomorphic = False
	for i in range(-1,(-1*len(exten_filelist[0])),-1):
		if len(char_dict[i]) == 2:
			dichotomous_sites.append(i)
			if first_dichotomous == 0:
				first_dichotomous = i
			searching_for_last_monomorphic = True
		elif len(char_dict[i]) >1 and searching_for_last_monomorphic == True:
			searching_for_last_monomorphic = False
		elif len(char_dict[i]) == 1 and searching_for_last_monomorphic == True:
			last_monomorphic = i
	exten_count = {}
	for f in range(0,len(exten_filelist)):
		filename = exten_filelist[f]
		try:
			exten_count[filename[last_monomorphic:]] += 1
		except:
			exten_count[filename[last_monomorphic:]] = 1
	exten_pair = {1:'',2:''}
	for exten in exten_count:
		try:
			exten_pair[int(exten[first_dichotomous])] = exten
		except:
			pass
	successful_prediction = False
	if len([f for f in os.listdir(input_read_dir) if f.endswith(exten_pair[1])]) == len([f for f in os.listdir(input_read_dir) if f.endswith(exten_pair[2])]):
		if len(exten_filelist)/2 == len([f for f in os.listdir(input_read_dir) if f.endswith(exten_pair[1])]):
			successful_prediction = True
			return exten_pair[1],exten_pair[2]
	if successful_prediction == False:
		print("Unable to predict forward and reverse read suffix values. Enter manually to proceed.\n")
		return "",""

def predict_single_file_extension(input_read_dir):
	exten_filelist = [f for f in os.listdir(input_read_dir) if f.endswith(".fastq") or f.endswith(".fq")]
	char_dict = {}
	for i in range(-1,(-1*len(exten_filelist[0])),-1):
		char_dict[i] = {}
		for f in range(0,len(exten_filelist)):
			filename = exten_filelist[f]
			try:
				char = filename[i]
			except:
				char = "na"
			try:
				char_dict[i][char] += 1
			except:
				char_dict[i][char] = 1
	last_monomorphic = 0
	searching_for_last_monomorphic = True
	for i in range(-1,(-1*len(exten_filelist[0])),-1):
		if len(char_dict[i]) == 1 and searching_for_last_monomorphic == True:
			last_monomorphic = i
		elif len(char_dict[i]) > 1 and searching_for_last_monomorphic == True:
			searching_for_last_monomorphic = False
	exten_count = {}
	for f in range(0,len(exten_filelist)):
		filename = exten_filelist[f]
		try:
			exten_count[filename[last_monomorphic:]] += 1
		except:
			exten_count[filename[last_monomorphic:]] = 1
	successful_prediction = False
	if exten_count[filename[last_monomorphic:]] == len(exten_filelist):
		successful_prediction = True
		return filename[last_monomorphic:]
	if successful_prediction == False:
		print("Unable to predict read suffix value. Enter manually to proceed.\n")
		return ""


def raw_read_processing_single(sampleID,input_read_dir,trim_dir,temp_dir,forward_read_suffix,command_prefix):
	global verbose
	if verbose == True:
		update_print_line("Processing single-end input reads",sampleID)
		# print("\tProcessing input reads for: "+sampleID)
	input_reads_forward = input_read_dir+sampleID+forward_read_suffix
	temp_adapter_filename = temp_dir+sampleID+".adapters.fa"
	pretrim_fastq_filename = temp_dir+sampleID+".pretrim.fq"
	trim_fastq_filename = temp_dir+sampleID+".trim.fq"
	clean_filename_out = trim_dir+sampleID+".fq"

	command = command_prefix+'bbmerge.sh in="'+input_reads_forward+'" outadapter="'+temp_adapter_filename+'" '+bbmerge_options
	out = run_command(command)
	time.sleep(5)
	adapter_pass = check_adapter(temp_adapter_filename)
	if adapter_pass == True:
		command = command_prefix+'bbmerge.sh in="'+input_reads_forward+'" out="'+pretrim_fastq_filename+'" adapters="'+temp_adapter_filename+'" '+bbmerge_options
		command_status = run_command(command)
		if command_status == False:
			sys.exit("bbmerge failed: "+input_reads_forward)
		time.sleep(1)
		
		command = command_prefix+'bbduk.sh in="'+pretrim_fastq_filename+'" out="'+trim_fastq_filename+'" '+bbduk_options
		command_status = run_command(command)
		if command_status == False:
			sys.exit("bbduk failed: "+pretrim_fastq_filename)
		time.sleep(1)
		
		os.remove(pretrim_fastq_filename)
	else:
		command = command_prefix+'bbduk.sh in="'+input_reads_forward+'" out="'+trim_fastq_filename+'" '+bbduk_options
		command_status = run_command(command)
		if command_status == False:
			sys.exit("bbduk failed: "+input_reads_forward)
		time.sleep(1)


	os.remove(temp_adapter_filename)
	os.rename(trim_fastq_filename,clean_filename_out)
	return clean_filename_out


def raw_read_processing_paired(sampleID,input_read_dir,trim_dir,temp_dir,forward_read_suffix,reverse_read_suffix,command_prefix):
	global verbose
	if verbose == True:
		update_print_line("Processing paired-end input reads",sampleID)
	input_reads_forward = input_read_dir+sampleID+forward_read_suffix
	input_reads_reverse = input_read_dir+sampleID+reverse_read_suffix

	repair_fastq_filename = temp_dir+sampleID+".repair.fq"
	temp_adapter_filename = temp_dir+sampleID+".adapters.fa"
	pretrim_fastq_filename = temp_dir+sampleID+".merge_pretrim.fq"
	merged_trimmed_filename = temp_dir+sampleID+".merge.fq"
	clean_filename_out = trim_dir+sampleID+".fq"

	command = command_prefix+'repair.sh in="'+input_reads_forward+'" in2="'+input_reads_reverse+'" out="'+repair_fastq_filename+'" '+repair_options
	command_status = run_command(command)
	if command_status == False:
		sys.exit("repair.sh failed: "+input_reads_forward)
	time.sleep(1)
	
	command = command_prefix+'bbmerge.sh in="'+repair_fastq_filename+'" outadapter="'+temp_adapter_filename+'" '+bbmerge_options
	command_status = run_command(command)
	if command_status == False:
		sys.exit("bbmerge.sh -outadapt failed: "+repair_fastq_filename)
	time.sleep(1)

	adapter_pass = check_adapter(temp_adapter_filename)
	if adapter_pass == True:
		command = command_prefix+'bbmerge.sh in="'+repair_fastq_filename+'" out="'+pretrim_fastq_filename+'" adapters="'+temp_adapter_filename+'" '+bbmerge_options
		command_status = run_command(command)
		if command_status == False:
			sys.exit("bbmerge.sh failed: "+repair_fastq_filename)
		time.sleep(1)
	else:
		command = command_prefix+'bbmerge.sh in="'+repair_fastq_filename+'" out="'+pretrim_fastq_filename+'" '+bbmerge_options
		command_status = run_command(command)
		if command_status == False:
			sys.exit("bbmerge.sh failed: "+repair_fastq_filename)
		time.sleep(1)
	
	command = command_prefix+'bbduk.sh in="'+pretrim_fastq_filename+'" out="'+merged_trimmed_filename+'" '+bbduk_options
	command_status = run_command(command)
	if command_status == False:
		sys.exit("bbduk.sh failed: "+pretrim_fastq_filename)
	time.sleep(1)
	
	os.remove(repair_fastq_filename)
	os.remove(temp_adapter_filename)
	os.remove(pretrim_fastq_filename)
	os.rename(merged_trimmed_filename,clean_filename_out)
	
	return clean_filename_out


###########################################################

def main_pipeline(sampleID,ref_set_list,input_read_dir,ref_dir,trim_dir,map_dir,qual_dir,sam_dir,bam_dir,temp_dir,sub_dir,min_evalue,min_bitscore,raw_read_input_type,forward_read_suffix,reverse_read_suffix,command_prefix,overwrite_existing_files):
	if raw_read_input_type == "mixed":
		local_raw_read_input_type = determine_raw_read_input_type(sampleID,input_read_dir,forward_read_suffix,reverse_read_suffix)
	else:
		local_raw_read_input_type = raw_read_input_type
	# print(sampleID+"\t"+local_raw_read_input_type)
	input_reads_found = check_if_input_reads_exist(sampleID,input_read_dir,local_raw_read_input_type,forward_read_suffix,reverse_read_suffix)
	trimmed_reads_found = check_trimmed_reads_exist(sampleID,trim_dir)
	if trimmed_reads_found == False and input_reads_found == True:
		if local_raw_read_input_type == "paired":
			trimmed_filename = raw_read_processing_paired(sampleID,input_read_dir,trim_dir,temp_dir,forward_read_suffix,reverse_read_suffix,command_prefix)
		elif local_raw_read_input_type == "single":
			trimmed_filename = raw_read_processing_single(sampleID,input_read_dir,trim_dir,temp_dir,forward_read_suffix,command_prefix)
		else:
			sys.exit('Invaild read input type. Only valid options are: "paired" or "single" - Exiting.')
	elif trimmed_reads_found == False and input_reads_found == False:
		sys.exit('Unable to find trimmed reads or input reads for: '+sampleID)
	split_reads_found = check_if_reads_split(sampleID,ref_set_list,trim_dir)
	if split_reads_found == False:
		read_count_dict = split_fastq_by_segment(sampleID,ref_set_list,temp_dir,trim_dir,ref_dir,min_evalue,min_bitscore)
	else:
		read_count_dict = load_read_counts(sampleID,ref_set_list,trim_dir)
	
	reads_mapped,missing_ref_sets = check_if_reads_mapped(sampleID,ref_set_list,sam_dir,read_count_dict)
	if reads_mapped == False:
		map_reads_to_ref(sampleID,missing_ref_sets,trim_dir,sam_dir,read_count_dict,overwrite_existing_files)
	else:
		return None
	return sampleID

####################### MAIN ##############################

if not os.path.exists(trim_dir):
	os.makedirs(trim_dir)
if not os.path.exists(map_dir):
	os.makedirs(map_dir)
if not os.path.exists(qual_dir):
	os.makedirs(qual_dir)
if not os.path.exists(sam_dir):
	os.makedirs(sam_dir)
if not os.path.exists(sub_dir):
	os.makedirs(sub_dir)
if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)


#Set up parallel processing if enabled
if parallel_process == True:
	import multiprocessing
	from joblib import Parallel, delayed
	if parallel_max_cpu == 0:
		num_cores = multiprocessing.cpu_count()
	elif parallel_max_cpu > 0:
		num_cores = min(parallel_max_cpu,multiprocessing.cpu_count())
	elif parallel_max_cpu < 0:
		num_cores = multiprocessing.cpu_count()+parallel_max_cpu
else:
	num_cores = 1

# sampleID_file_name = "sampleIDs.txt"
# flagged_samples = {}
# sample_list_filename = "samples.txt"
# infile = open(project_dir+filename,"r")
# for line in infile:
# 	line = line.strip()
# 	# if line in sampleIDs_to_focus:
# 	flagged_samples[line] = ''
# infile.close


### Check if all references have been concatenated into a single file yet
check_refs(ref_set_list,ref_dir)

### Check if the files in the input directory need to be unzipped
files_found = check_if_input_zipped(input_read_dir,parallel_process,num_cores)
if files_found == False:
	sys.exit("No input files with '.fastq', '.fq', '.gz', or '.tar.gz' extensions found. Exiting.")

### Predict forward and reverse read suffix if left blank
continue_running = False

suffix_info_file_path = temp_dir+"file_extensions.txt"
if forward_read_suffix != '':
	if os.path.isfile(suffix_info_file_path) == False:
		suffix_info_file = open(suffix_info_file_path,"w")
		suffix_info_file.write("forward\t"+forward_read_suffix+"\nreverse\t"+reverse_read_suffix+"\n")
		suffix_info_file.close()
if forward_read_suffix == '':
	# if os.path.isfile(suffix_info_file_path) == True:
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
		if raw_read_input_type == "paired" and reads_already_screened_for_quality == False:
			forward_read_suffix,reverse_read_suffix = predict_paired_file_extensions(input_read_dir)
			if forward_read_suffix != "" and reverse_read_suffix != "":
				print('Predicted read suffix\n\tForward: "'+str(forward_read_suffix)+'"\n\tReverse: "'+str(reverse_read_suffix)+'"')
				print("Is this correct?")
				decision = input("(yes or no)\n")
			else:
				decision = "no"
			if decision == "Y" or decision == "y" or decision == "Yes" or decision == "yes" or decision == "YES":
				continue_running = True
			else:
				forward_read_suffix = input('Enter forward read suffix (ex "_L001_R1_001.fastq") or "exit" to cancel:\n')
				if forward_read_suffix == "exit":
					sys.exit()
				reverse_read_suffix	= input('Enter Reverse read suffix (ex "_L001_R2_001.fastq") or "exit" to cancel:\n')
				if reverse_read_suffix == "exit":
					sys.exit()
				if forward_read_suffix != "exit" and reverse_read_suffix != "exit":
					temp_file_list = [f for f in os.listdir(input_read_dir) if f.endswith(".fastq") or f.endswith(".fq")]
					example = temp_file_list[0]
					if forward_read_suffix in example:
						example = example.split(forward_read_suffix)[0]
					elif reverse_read_suffix in example:
						example = example.split(reverse_read_suffix)[0]
					else:
						sys.exit("File suffix not found in files, please re-run this script and re-enter the expected file suffix")
					
					decision = input("Example sample name after removing read suffix:\n"+example+"\nDoes this sample name look correct?\n")
					if decision == "Y" or decision == "y" or decision == "Yes" or decision == "yes" or decision == "YES":
						continue_running = True
			if continue_running == True:
				suffix_info_file = open(suffix_info_file_path,"w")
				suffix_info_file.write("forward\t"+forward_read_suffix+"\nreverse\t"+reverse_read_suffix+"\n")
				suffix_info_file.close()
		elif raw_read_input_type == "single" or reads_already_screened_for_quality == True:
			forward_read_suffix = predict_single_file_extension(input_read_dir)
			if forward_read_suffix != "":
				print('Predicted read suffix: "'+str(forward_read_suffix)+'"')
				print("Is this correct?")
				decision = input("(yes or no)\n")
			else:
				decision = "no"
			if decision == "Y" or decision == "y" or decision == "Yes" or decision == "yes" or decision == "YES":
				continue_running = True
			else:
				forward_read_suffix = input('Enter read suffix (ex "_L001_001.fastq") or "exit" to cancel:\n')
				if forward_read_suffix == "exit":
					sys.exit()
				else:
					temp_file_list = [f for f in os.listdir(input_read_dir) if f.endswith(".fastq") or f.endswith(".fq")]
					example = temp_file_list[0]
					if forward_read_suffix in example:
						example = example.split(forward_read_suffix)[0]
					else:
						sys.exit("File suffix not found in files, please re-run this script and re-enter the expected file suffix")
					
					decision = input("Example sample name after removing read suffix:\n"+example+"\nDoes this sample name look correct?\n")
					if decision == "Y" or decision == "y" or decision == "Yes" or decision == "yes" or decision == "YES":
						continue_running = True
					continue_running = True
			if continue_running == True:
				suffix_info_file = open(suffix_info_file_path,"w")
				suffix_info_file.write("forward\t"+forward_read_suffix+"\nreverse\t"+reverse_read_suffix+"\n")
				suffix_info_file.close()
		elif raw_read_input_type != "mixed":
			forward_read_suffix = input('Enter forward read suffix (ex "_L001_R1_001.fastq") or "exit" to cancel:\n')
			if forward_read_suffix == "exit":
				sys.exit()
			reverse_read_suffix	= input('Enter Reverse read suffix (ex "_L001_R2_001.fastq") or "exit" to cancel:\n')
			if reverse_read_suffix == "exit":
				sys.exit()
			if forward_read_suffix != "exit" and reverse_read_suffix != "exit":
				temp_file_list = [f for f in os.listdir(input_read_dir) if f.endswith(".fastq") or f.endswith(".fq")]
				example = temp_file_list[0]
				if forward_read_suffix in example:
					example = example.split(forward_read_suffix)[0]
				elif reverse_read_suffix in example:
					example = example.split(reverse_read_suffix)[0]
				else:
					sys.exit("File suffix not found in files, please re-run this script and re-enter the expected file suffix")
				
				decision = input("Example sample name after removing read suffix:\n"+example+"\nDoes this sample name look correct?\n")
				if decision == "Y" or decision == "y" or decision == "Yes" or decision == "yes" or decision == "YES":
					continue_running = True
		if continue_running == True:
			suffix_info_file = open(suffix_info_file_path,"w")
			suffix_info_file.write("forward\t"+forward_read_suffix+"\nreverse\t"+reverse_read_suffix+"\n")
			suffix_info_file.close()

if raw_read_input_type == "paired" or raw_read_input_type == "mixed":
	if forward_read_suffix != "" and reverse_read_suffix != "":
		continue_running = True
	else:
		continue_running = False
elif raw_read_input_type == "single":
	if forward_read_suffix != "":
		continue_running = True
	else:
		continue_running = False

if continue_running == False:
	sys.exit("No read suffix entered, please re-run and enter information when prompted")


## Create list of samples to process
file_list = [f for f in os.listdir(input_read_dir) if f.endswith(forward_read_suffix)]
sampleID_list = []
print_space_offset = 0
for files in file_list:
	sampleID = files.split(forward_read_suffix)[0]
	sampleID_list.append(sampleID)
	if len(sampleID) > print_space_offset:
		print_space_offset = len(sampleID)
sampleID_list = list(set(sampleID_list))


# sampleID_list = ['SC002-01302020-rep2']
# sampleID_list = ['CC0301-01182022-rep1']

# if len(sampleID_list) < num_cores:
# 	num_ref_sets = len(ref_set_list)
# 	num_samples = len(sampleID_list)
# 	if (num_cores/num_samples/num_ref_sets) >1:
# num_sample_process_cores = min(len(sampleID_list),num_cores)


print("A total of "+str(len(sampleID_list))+" samples were found.")
# sampleID_list = ['CC0831-802986-run1-rep1','CC0216-804417-run1-rep1']

### Process all samples
processed_list = []
if parallel_process == True:
	processed_list = Parallel(n_jobs=num_cores)(delayed(main_pipeline)(sampleID,ref_set_list,input_read_dir,ref_dir,trim_dir,map_dir,qual_dir,sam_dir,bam_dir,temp_dir,sub_dir,min_evalue,min_bitscore,raw_read_input_type,forward_read_suffix,reverse_read_suffix,command_prefix,overwrite_existing_files) for sampleID in sampleID_list)
elif parallel_process == False:
	for sampleID in sampleID_list:
		main_pipeline(sampleID,ref_set_list,input_read_dir,ref_dir,trim_dir,map_dir,qual_dir,sam_dir,bam_dir,temp_dir,sub_dir,min_evalue,min_bitscore,raw_read_input_type,forward_read_suffix,reverse_read_suffix,command_prefix,overwrite_existing_files)
		processed_list.append(sampleID)
