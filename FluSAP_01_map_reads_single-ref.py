import os
import sys
import math
import random
import numpy as np
import time

verbose = True
overwrite_existing_files = False #when False, pipeline will avoid re-running samples that have already been processed 
reads_already_screened_for_quality = False # 'True' or 'False' --- If reads were already merged and screened for average quality and length (such as during a demultiplexing step), the script will not attempt to filter or merge the input reads

command_prefix = "bash "

parallel_process = True
parallel_max_cpu = 30

#BLAT Settings
min_evalue = 1.0E-5
min_bitscore = 100.0


raw_read_input_type = "mixed" # 'paired', 'single', or 'mixed' --- When 'single' the pipeline will look for barcodes in one file, when 'paired' the pipeline will look for forward and reverse reads
forward_read_suffix = "_R1.fastq" #if empty, the script will attempt to predict these values, but it only works under specific assumptions
reverse_read_suffix = "_R2.fastq" #this value is ignored if 'reads_already_screened_for_quality' is 'True' or 'raw_read_input_type' is 'single' - NOTE: this variable cannot be removed without causing missing value errors


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
consensus_dir = map_dir+"consensus/"
temp_dir = map_dir+"temp/"

ref_set_list = ['IAV_H1.fa','IAV_H3.fa','IAV_M.fa','IAV_N1.fa','IAV_N2.fa','IAV_NP.fa','IAV_NS.fa','IAV_PA.fa','IAV_PB1.fa','IAV_PB2.fa','IBV_HA.fa','IBV_M.fa','IBV_NA.fa','IBV_NP.fa','IBV_NS.fa','IBV_PA.fa','IBV_PB1.fa','IBV_PB2.fa']

repair_options = 'fint=t repair=t overwrite=t threads=1 -da -Xmx3000m'
bbduk_adapter_trim_options = 'ktrim=r k=23 mink=11 hdist=1 minlen=80 tpe tbo overwrite=t threads=4 -Xmx1000m'
bbmerge_options = 'qin=33 ecco mix overwrite=t threads=4 -Xmx1000m'
bbduk_options = 'qin=33 tpe=f tbo=f overwrite=t threads=4 -Xmx1000m'#maq=30  lhist=lhist.txt qhist=qhist.txt bhist=bhist.txt aqhist=aqhist.txt'
bbmap_options = 'qin=33 local=f touppercase=t overwrite=t threads=4 nodisk -Xmx2000m'#ambiguous=toss 


############################## Functions #############################################
def run_command(command,mode="quiet"):
# def run_command(command,mode="verbose"):
	run_status = False
	if mode == "verbose":
		return_code = os.system(command)
	elif mode == "quiet":
		return_code = os.system(command+" >/dev/null 2>&1")
	# subprocess.check_output(command,stderr=subprocess.STDOUT)
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
	fq_inputs = [f for f in os.listdir(input_read_dir) if f.endswith(".fastq") or f.endswith(".fq")]
	if len(fq_inputs) >0:
		files_found = True
	else:
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

	return files_found

def map_reads_to_ref(sampleID,map_ref_set_list,trim_dir,sam_dir,consensus_dir,read_count_dict,overwrite_existing_files):
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
				sample_consensus_ref = consensus_dir+sampleID+"."+ref_ID+".fa"
				### Map reads to reference
				input_filename_1 = split_fastq
				output_filename_1 = temp_dir+sampleID+"."+ref_ID+".sam"
				final_output_filename = sam_dir+ref_ID+"/"+sampleID+"."+ref_ID+".sam"
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
					update_print_line("Mapping processed reads - "+ref_ID,sampleID)
					command = command_prefix+"bbmap.sh in="+input_filename_1+" outm="+output_filename_1+" ref="+sample_consensus_ref+" "+bbmap_options
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
							output_filename_1 = sam_dir+ref_ID+"/"+sampleID+"."+ref_ID+".refine.sam"
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
			if blat_list[0][0] < blat_list[1][0]:
				blat_dict[query] = blat_list[0][1]
			elif blat_list[0][0] == blat_list[1][0]:
				blat_dict[query] = blat_list[random.randint(0,1)][1]
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



###########################################################

def main_pipeline(sampleID,ref_set_list,input_read_dir,ref_dir,trim_dir,map_dir,qual_dir,sam_dir,bam_dir,temp_dir,sub_dir,consensus_dir,min_evalue,min_bitscore,raw_read_input_type,forward_read_suffix,reverse_read_suffix,command_prefix,overwrite_existing_files):
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
		map_reads_to_ref(sampleID,missing_ref_sets,trim_dir,sam_dir,consensus_dir,read_count_dict,overwrite_existing_files)
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
if not os.path.exists(consensus_dir):
	os.makedirs(consensus_dir)
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


### Check if all references have been concatenated into a single file yet
check_refs(ref_set_list,ref_dir)

### Check if the files in the input directory need to be unzipped
files_found = check_if_input_zipped(input_read_dir,parallel_process,num_cores)
if files_found == False:
	sys.exit("No input files with '.fastq', '.fq', '.gz', or '.tar.gz' extensions found. Exiting.")


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
	# if sampleID == 'SC014-242912_run2' or sampleID == 'CC0386-644596_run1' or sampleID == 'CC0386-644596_run2' or sampleID == 'CC0445-821714_run1' or sampleID == 'CC0445-821714_run2':
	sampleID_list.append(sampleID)
	if len(sampleID) > print_space_offset:
		print_space_offset = len(sampleID)
sampleID_list = list(set(sampleID_list))


print("A total of "+str(len(sampleID_list))+" samples were found.")
print("Total cores requested: "+str(num_cores))

### Process all samples
processed_list = []
if parallel_process == True:
	processed_list = Parallel(n_jobs=num_cores)(delayed(main_pipeline)(sampleID,ref_set_list,input_read_dir,ref_dir,trim_dir,map_dir,qual_dir,sam_dir,bam_dir,temp_dir,sub_dir,consensus_dir,min_evalue,min_bitscore,raw_read_input_type,forward_read_suffix,reverse_read_suffix,command_prefix,overwrite_existing_files) for sampleID in sampleID_list)
elif parallel_process == False:
	for sampleID in sampleID_list:
		main_pipeline(sampleID,ref_set_list,input_read_dir,ref_dir,trim_dir,map_dir,qual_dir,sam_dir,bam_dir,temp_dir,sub_dir,consensus_dir,min_evalue,min_bitscore,raw_read_input_type,forward_read_suffix,reverse_read_suffix,command_prefix,overwrite_existing_files)
		processed_list.append(sampleID)
