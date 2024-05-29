import os
import sys
import numpy as np

dedupe_reads = True
skip_qual = True
skip_linked_SNV_counting = True
verbose = True
min_proportion = 0.03

parallel_process = True
parallel_max_cpu = 45

''' Settings for linked iSNV analysis only '''
min_lnkSNP_qual = 30.0
min_lnkSNP_map_qual = 35.0
min_lnkSNP_cov = 100
min_lnkSNP_loc = 20
min_lnkSNP_prop = 0.03

'''  '''
###################################   Setup   ###################################

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
command_prefix = "bash "

ref_set_list = ['IAV_H1.fa','IAV_H3.fa','IAV_M.fa','IAV_N1.fa','IAV_N2.fa','IAV_NP.fa','IAV_NS.fa','IAV_PA.fa','IAV_PB1.fa','IAV_PB2.fa','IBV_HA.fa','IBV_M.fa','IBV_NA.fa','IBV_NP.fa','IBV_NS.fa','IBV_PA.fa','IBV_PB1.fa','IBV_PB2.fa']

###################################   Setup   ###################################
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

#Make output directories if they don't exist yet
if not os.path.exists(sub_dir):
	os.makedirs(sub_dir)
if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)
if skip_qual == False:
	if not os.path.exists(qual_dir):
		os.makedirs(qual_dir)
# if skip_linked_SNV_counting == False:
# 	lnk_dir = map_dir+"lnkSNV/"
# 	if not os.path.exists(lnk_dir):
# 		os.makedirs(lnk_dir)
# 	if not os.path.exists(lnk_dir+ref_set.split("-")[0].split(".f")[0]):
# 		os.makedirs(lnk_dir+ref_set.split("-")[0].split(".f")[0])


############################################ FUNCTIONS ############################################
def run_command(command):
	run_status = False
	# return_code = os.system(command)
	return_code = os.system(command+" >/dev/null 2>&1")
	# subprocess.check_output(command,stderr=subprocess.STDOUT)
	if return_code == 0:
		run_status = True
	return run_status

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


def check_if_reads_mapped(sampleID,ref_set,sam_dir,read_count_dict):
	reads_mapped = True
	# for ref_set in ref_set_list:
	refs = ref_set.split("-")
	primary_ref = ref_set.split("-")[0].split(".f")[0]
	for ref_filename in refs:
		refID = ref_filename.split('.f')[0]
		if os.path.isfile(sam_dir+refID+"/"+sampleID+"."+refID+'.refine.sam') == False:
			try:
				read_count = read_count_dict[refID]
			except:
				read_count = 0
			if read_count >= 1000:
				reads_mapped = False
	return reads_mapped


def check_if_processing_finished(sampleID,ref_set,sub_dir):
	sam_processed = True
	# for ref_set in ref_set_list:
	refs = ref_set.split("-")
	for ref_filename in refs:
		refID = ref_filename.split(".f")[0]
		out_substitutions_filename = sub_dir+refID+"/"+sampleID+"."+refID+".substitions.txt"
		if os.path.isfile(out_substitutions_filename) == False:
			sam_processed = False
	return sam_processed


def load_read_counts(sampleID,trim_dir):
	read_count_dict = {}
	for ref_set in ref_set_list:
		refs = ref_set.split("-")
		primary_ref = ref_set.split("-")[0].split(".f")[0]
		for ref_filename in refs:
			segment = ref_filename.split(".f")[0]
			read_count_filename = trim_dir+segment+"/"+sampleID+'.'+segment+'.read_count.txt'
			if os.path.isfile(read_count_filename) == True:
				read_count_file = open(read_count_filename,'r')
				for line in read_count_file:
					line = line.strip()
					read_count_dict[segment] = int(line)
				read_count_file.close()
	return read_count_dict


def dedupe_samfile(sam_filename_in,temp_sam_filename,sam_filename_out):
	sam_infile = open(sam_filename_in,"r")
	read_map_strings = {}
	for sam_line in sam_infile:
		if len(sam_line) > 0:
			if sam_line[0] != "@":
				sam_line = sam_line.strip().split("\t")
				read_seq = sam_line[9]
				if read_seq !="*":
					readID = sam_line[0].split(" ")[0]
					map_site = sam_line[3]
					map_qual = sam_line[4]
					cigar = sam_line[5]
					map_string = map_site+"_"+map_qual+"_"+cigar
					try:
						read_map_strings[readID] += " "+map_string
					except:
						read_map_strings[readID] = map_string
	sam_infile.close()
	map_string_dedupe = {}
	for readID in read_map_strings:
		map_string = read_map_strings[readID]
		try:
			one_readID = map_string_dedupe[map_string]
		except:
			map_string_dedupe[map_string] = readID
	del read_map_strings
	pass_readID_dict = {}
	for map_string in map_string_dedupe:
		readID = map_string_dedupe[map_string]
		pass_readID_dict[readID] = ''
	del map_string_dedupe
	dedupe_sam_outfile = open(temp_sam_filename,"w")
	sam_infile = open(sam_filename_in,"r")
	for sam_line in sam_infile:
		if len(sam_line) > 0:
			sam_line = sam_line.strip()#.split("\t")
			if sam_line[0] != "@":
				split_line = sam_line.split("\t")
				readID = split_line[0].split(" ")[0]
				read_seq = split_line[9]
				if read_seq !="*":
					try:
						pass_readID_dict[readID]
						dedupe_sam_outfile.write(sam_line+"\n")
					except:
						pass
			else:
				dedupe_sam_outfile.write(sam_line+"\n")
	sam_infile.close()
	dedupe_sam_outfile.close()
	os.rename(temp_sam_filename,sam_filename_out)
	return len(pass_readID_dict)


def SAM_parse(map_position,cigar_string,quality_string,read_seq,mapping_qual): #mapping position is indexed to 1, not zero
	cigar_dict = {'M':'imperfect_match','D':'insertion_in_ref','I':'insertion_in_read','S':'soft_clip_reference','H':'soft_clip_read','=':'perfect_match','X':'mismatch'}
	
	read_key = ''
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
				if cur_operator == "D":
					cur_position += 1
					# ref_site_map.append('')
					# read_key += cur_operator
				elif cur_operator == "I":
					cur_position += 0
					ref_site_map.append('')
					read_key += cur_operator
				elif cur_operator == "S":
					cur_position += 0
					ref_site_map.append('')
					read_key += cur_operator
				elif cur_operator == "H":
					all_valid_operators = False
				elif cur_operator == "=": #perfect match
					cur_position += 1
					ref_site_map.append(cur_position)
					read_key += cur_operator
				elif cur_operator == "M": #imperfect match
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
	base_dict_fun = {}
	read_end_site = 0
	if all_valid_operators == True:
		for num in range(0,len(read_key)):
			operator = read_key[num]
			ref_loc = ref_site_map[num]
			if ref_loc != '':
				ref_loc = int(ref_loc)
				try:
					qual_char = quality_string[num]
				except:
					sys.exit("Index ")
				qual_score = ord(qual_char)-33
				read_nt = read_seq[num]
				read_length = len(read_seq)
				nt_prop_location = round(num/read_length,2)
				R_end_dist = read_length-num
				# if qual_score >= 30 and num >= 30 and R_end_dist>=30 and mapping_qual>=30:
				tup = (read_nt,qual_score,mapping_qual,read_length,num,R_end_dist,nt_prop_location)
				base_dict_fun[ref_loc] = tup
				# try:
				# except:
				# 	base_dict_fun[ref_loc] = []
				# 	base_dict_fun[ref_loc].append(tup)
	# read_info_out_list = []
	# for ref_loc in base_dict_fun:
	# 	out_tup = (ref_loc,base_dict_fun[ref_loc])
	# 	read_info_out_list.append(out_tup)
	# del base_dict_fun
	return base_dict_fun,match_count_dict
	# return read_info_out_list,match_count_dict


def CIGAR_parse(map_position,cigar_string,quality_string,read_seq,mapping_qual,flank_ignore_len): #mapping position is indexed to 1, not zero
	cigar_dict = {'M':'imperfect_match','D':'insertion_in_ref','I':'insertion_in_read','S':'soft_clip_reference','H':'soft_clip_read','=':'perfect_match','X':'mismatch'}
	
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
			if len(read_key)>=flank_ignore_len and len(read_key)<=(len(read_seq)-flank_ignore_len):
				try:
					match_count_dict[cur_operator] += cur_num
				except:
					match_count_dict[cur_operator] = cur_num
			for num in range(0,cur_num):
				if cur_operator == "D": #ref base not present in read
					cur_position += 1
					ref_deletion_map[cur_position] = len(read_key)
					# ref_site_map.append(cur_position)
					# read_key += cur_operator
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
	base_dict_fun = {}
	read_end_site = 0
	if all_valid_operators == True:
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
		if ref_deletion_map != {}:
			for ref_loc in ref_deletion_map:
				read_nt = '-'
				qual_score = 40
				read_loc = ref_deletion_map[ref_loc]
				read_length = len(read_seq)
				R_end_dist = read_length-read_loc
				nt_prop_location = round(read_loc/read_length,2)
				tup = (read_nt,qual_score,mapping_qual,read_length,read_loc,R_end_dist,nt_prop_location)
				base_dict_fun[ref_loc] = tup
		if ref_insert_map != {}:
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

###########################################################

def main_pipeline(sample_ref_tup,ref_set_list,input_read_dir,ref_dir,trim_dir,map_dir,qual_dir,sam_dir,bam_dir,temp_dir,sub_dir,command_prefix,dedupe_reads):
	max_indel_count = 24
	max_mismatch_count = 12

	min_base_qual_for_counting = 15
	min_map_qual_for_counting = 15
	min_loc_for_counting = 10

	sampleID = sample_ref_tup[0]
	ref_set = sample_ref_tup[1]
	refs = ref_set.split("-")
	primary_refID = refs[0].split(".f")[0]

	read_count_dict = load_read_counts(sampleID,trim_dir)
	reads_mapped = check_if_reads_mapped(sampleID,ref_set,sam_dir,read_count_dict)
	if reads_mapped == False:
		return None
	
	map_files_processed = check_if_processing_finished(sampleID,ref_set,sub_dir)
	if map_files_processed == True:
		return None
	
	print(sampleID+" "+primary_refID)
	for ref_filename in refs:
		read_map_info_dict = {}
		map_qual_dict = {}
		indel_dict = {}
		mismatch_dict = {}
		refID = ref_filename.split(".f")[0]
		sam_filename = sam_dir+refID+"/"+sampleID+"."+refID+".refine.sam"
		cont = True
		if dedupe_reads  == False:
			sam_filename_to_open = sam_filename
			if os.path.isfile(sam_filename) == False:
				cont = False
		elif dedupe_reads == True:
			if os.path.isfile(sam_filename) == True:
				# print("\t\t"+sam_filename)
				temp_sam_filename = temp_dir+sampleID+"."+ref_filename.split(".f")[0]+".dedupe.sam"
				dedupe_samfilename = sam_dir+refID+"/"+sampleID+"."+ref_filename.split(".f")[0]+".dedupe.sam"
				num_reads_after_dedupe = dedupe_samfile(sam_filename,temp_sam_filename,dedupe_samfilename)
				sam_filename_to_open = dedupe_samfilename
			else:
				cont = False
		if cont == True:
			sam_infile = open(sam_filename_to_open,"r")
			for sam_line in sam_infile:
				if len(sam_line) > 0:
					if sam_line[0] != "@":
						sam_line = sam_line.strip().split("\t")
						map_qual = int(sam_line[4])
						read_seq = sam_line[9]
						if read_seq !="*":
							readID = sam_line[0]
							map_site = int(sam_line[3])
							cigar = sam_line[5] #{'M':'imperfect_match','D':'insertion_in_ref','I':'insertion_in_read','S':'soft_clip_reference','H':'soft_clip_read','=':'perfect_match','X':'mismatch'}
							read_seq_qual = sam_line[10]
							read_base_dict,CIGAR_operator_count_dict = CIGAR_parse(map_site,cigar,read_seq_qual,read_seq,map_qual,min_loc_for_counting)
							if len(read_base_dict) > 0: #per base in read: (read_nt,qual_score,mapping_qual,read_length,num,R_end_dist,nt_prop_location)
								try:
									in_count = CIGAR_operator_count_dict['I']
								except:
									in_count = 0
								try:
									del_count = CIGAR_operator_count_dict['D']
								except:
									del_count = 0
								indel_count = in_count + del_count
								try:
									mismatch_count = CIGAR_operator_count_dict['X']
								except:
									mismatch_count = 0
								if indel_count <= max_indel_count and mismatch_count <= max_mismatch_count:
									keep_reading = False
									try:
										prev_map_qual = read_map_info_dict[readID+"-"+read_seq[0:5]]
										if map_qual > prev_map_qual: #if higher quality mapping, keep reading to replace the entry for this read
											keep_reading = True
										else:
											keep_reading = False
									except: #if you haven't encountered the read ID-read_seq[0:6] combo, keep reading to store seq info
										keep_reading = True
									if keep_reading == True:
										read_map_info_dict[readID+"-"+read_seq[0:5]] = (map_qual,mismatch_count,indel_count,read_base_dict)
			sam_infile.close()
			if dedupe_reads == True:
				os.remove(dedupe_samfilename)
		if verbose == True:
			print(sampleID+"\t"+ref_set+"\tlen(read_map_info_dict):\t"+str(len(read_map_info_dict)))


		segment_isnv_link_count_dict = {}
		max_loc = 0
		sub_info_dict = {}
		allele_count_dict = {}
		count_dict = {}
		total_read_count = 0
		for readID_flag in read_map_info_dict:
			read_map_info_tup = read_map_info_dict[readID_flag]
			map_qual = read_map_info_dict[readID_flag][0]
			mismatch_count = read_map_info_dict[readID_flag][1]
			indel_count = read_map_info_dict[readID_flag][2]
			read_base_info_dict = read_map_info_dict[readID_flag][3]
			read_subIDs = []
			if map_qual >= min_map_qual_for_counting:
				total_read_count += 1
			for loc in read_base_info_dict:
				base = read_base_info_dict[loc][0]
				base_qual = read_base_info_dict[loc][1]
				map_qual = read_base_info_dict[loc][2]
				read_len = read_base_info_dict[loc][3]
				loc_in_read = read_base_info_dict[loc][4]
				loc_from_R_end = read_base_info_dict[loc][5]
				prop_loc = read_base_info_dict[loc][6]
				
				min_end_dist = min(loc_in_read,loc_from_R_end)
				max_end_dist = max(loc_in_read,loc_from_R_end)
				if min_end_dist >= min_loc_for_counting and base_qual >= min_base_qual_for_counting:
					loc_first = False
					try:
						count_dict[loc] += 1
					except:
						count_dict[loc] = 1
						sub_info_dict[loc] = {}
						allele_count_dict[loc] = {}
						loc_first = True
						if loc > max_loc:
							max_loc = loc
					base_first = False
					try:
						allele_count_dict[loc][base] += 1
					except:
						allele_count_dict[loc][base] = 1
						base_first = True
					if base_first == False:
						if skip_qual == True:
							sub_info_dict[loc][base][0] += base_qual
							sub_info_dict[loc][base][1] += map_qual
							sub_info_dict[loc][base][2] += read_len
							sub_info_dict[loc][base][3] += min_end_dist
							sub_info_dict[loc][base][4] += max_end_dist
							sub_info_dict[loc][base][5] += prop_loc
							sub_info_dict[loc][base][6] += mismatch_count
							sub_info_dict[loc][base][7] += indel_count
						elif skip_qual == False:
							temp_list1,temp_list2,temp_list3,temp_list4,temp_list5,temp_list6,temp_list7,temp_list8 = sub_info_dict[loc][base][0],sub_info_dict[loc][base][1],sub_info_dict[loc][base][2],sub_info_dict[loc][base][3],sub_info_dict[loc][base][4],sub_info_dict[loc][base][5],sub_info_dict[loc][base][6],sub_info_dict[loc][base][7]
							temp_list1.append(base_qual) # base_qual
							temp_list2.append(map_qual) # map_qual
							temp_list3.append(read_len) # read_len
							temp_list4.append(min_end_dist) # min_end_dist
							temp_list5.append(max_end_dist) # max_end_dist
							temp_list6.append(prop_loc) # prop_loc
							temp_list7.append(mismatch_count) # mismatch_count
							temp_list8.append(indel_count) # indel_count
							tup_out = (temp_list1,temp_list2,temp_list3,temp_list4,temp_list5,temp_list6,temp_list7,temp_list8)
							sub_info_dict[loc][base] = tup_out
					elif base_first == True:
						if skip_qual == True:
							tup_out = [base_qual,map_qual,read_len,min_end_dist,max_end_dist,prop_loc,mismatch_count,indel_count]
						elif skip_qual == False:
							tup_out = ([base_qual],[map_qual],[read_len],[loc_in_read],[loc_from_R_end],[prop_loc],[mismatch_count],[indel_count])
						sub_info_dict[loc][base] = tup_out
			if skip_linked_SNV_counting == True:
				read_map_info_dict[readID_flag] = ''
		if verbose == True:
			print(sampleID+"\t"+ref_set+"\tlen(sub_info_dict):\t"+str(len(sub_info_dict)))
			if skip_linked_SNV_counting == False:
				print("len(segment_isnv_link_count_dict): "+str(len(segment_isnv_link_count_dict)))
		if skip_linked_SNV_counting == True:
			del read_map_info_dict

		temp_substitutions_filename = temp_dir+sampleID+"."+refID+".substitions.txt"
		out_substitutions_filename = sub_dir+refID+"/"+sampleID+"."+refID+".substitions.txt"

		substitutions_outfile = open(temp_substitutions_filename,"w")
		substitutions_outfile.write("#loc\t"+str(total_read_count))
		info_list = ['base_freq','avg_qual','map_qual','read_len','min_end_dist','max_end_dist','prop_loc','mismatches','indels']
		info_string = ''
		base_string = ''
		for category in info_list:
			for i in range(0,len(bases)):
				base = bases[i]
				info_string += "\t"+category
				base_string += "\t"+base
		substitutions_outfile.write(info_string+"\n#loc\tdepth"+base_string+"\n")

		if skip_linked_SNV_counting == False:
			subID_dict = {}
		for loc in range(0,max_loc):
			try:
				loc_count = count_dict[loc]
			except:
				loc_count = 0

			substitutions_outfile.write(str(loc)+"\t"+str(loc_count))
			base_freq_string = ''
			avg_qual_string = ''
			map_qual_string = ''
			read_len_string = ''
			base_loc_string = ''
			base_R_loc_string = ''
			prop_loc_string = ''
			mismatch_count_string = ''
			indel_count_string = ''

			poly_count = 0

			for i in range(0,len(bases)):
				base = bases[i]
				try:
					base_count = allele_count_dict[loc][base]
				except:
					base_count = 0
				if base_count > 0:
					#sub_info_dict[loc][base]([base_qual],[map_qual],[read_len],[loc_in_read],[loc_from_R_end],[prop_loc],[mismatch_count],[indel_count])
					#sub_info_dict[loc][base]([    0    ],[    1   ],[    2   ],[     3     ],[      4       ],[    5   ],[      6       ],[     7     ])
					if skip_qual == True:
						prop_base = round((base_count/loc_count),4)
						avg_qual = round(sub_info_dict[loc][base][0]/base_count,1)
						avg_map_qual = round(sub_info_dict[loc][base][1]/base_count,1)
						avg_read_len = round(sub_info_dict[loc][base][2]/base_count,1)
						avg_base_loc = round(sub_info_dict[loc][base][3]/base_count,1)
						avg_base_loc_from_R = round(sub_info_dict[loc][base][4]/base_count,1)
						avg_prop_loc = round(sub_info_dict[loc][base][5]/base_count,2)
						avg_mismatch_count = float(round_to_n_sig_figs(sub_info_dict[loc][base][6]/base_count,2))
						avg_indel_count = float(round_to_n_sig_figs(sub_info_dict[loc][base][7]/base_count,2))
					elif skip_qual == False:
						prop_base = round((base_count/loc_count),4)
						avg_qual = round(np.average(sub_info_dict[loc][base][0]),1)
						avg_map_qual = round(np.average(sub_info_dict[loc][base][1]),1)
						avg_read_len = round(np.average(sub_info_dict[loc][base][2]),1)
						avg_base_loc = round(np.average(sub_info_dict[loc][base][3]),1)
						avg_base_loc_from_R = round(np.average(sub_info_dict[loc][base][4]),1)
						avg_prop_loc = round(np.average(sub_info_dict[loc][base][5]),2)
						avg_mismatch_count = float(round_to_n_sig_figs(np.average(sub_info_dict[loc][base][6]),2))
						avg_indel_count = float(round_to_n_sig_figs(np.average(sub_info_dict[loc][base][7]),2))

					if skip_linked_SNV_counting == False:
						if loc_count >= min_lnkSNP_cov and prop_base >= min_lnkSNP_prop and avg_map_qual >= min_lnkSNP_map_qual and avg_base_loc >= min_lnkSNP_loc:
							if prop_base <= (1.0-min_lnkSNP_prop):
								subID_dict[str(loc)+base] = (loc,loc_count,base,prop_base,avg_qual,avg_map_qual,avg_base_loc,avg_mismatch_count,avg_indel_count)
				else:
					prop_base = 0
					avg_qual = 0
					avg_map_qual = 0
					avg_read_len = 0
					avg_base_loc = 0
					avg_base_loc_from_R = 0
					avg_prop_loc = 0
					avg_mismatch_count = 0
					avg_indel_count = 0

				base_freq_string +="\t"+str(prop_base)
				avg_qual_string +="\t"+str(avg_qual)
				map_qual_string +="\t"+str(avg_map_qual)
				read_len_string +="\t"+str(avg_read_len)
				base_loc_string +="\t"+str(avg_base_loc)
				base_R_loc_string +="\t"+str(avg_base_loc_from_R)
				prop_loc_string +="\t"+str(avg_prop_loc)
				mismatch_count_string +="\t"+str(avg_mismatch_count)
				indel_count_string +="\t"+str(avg_indel_count)

			substitutions_outfile.write(base_freq_string+avg_qual_string+map_qual_string+read_len_string+base_loc_string+base_R_loc_string+prop_loc_string+mismatch_count_string+indel_count_string+"\n")
		substitutions_outfile.close()
		if os.path.isfile(out_substitutions_filename) == True:
			os.remove(out_substitutions_filename)
		os.rename(temp_substitutions_filename,out_substitutions_filename)
		if verbose == True:
			print(sampleID+"\t"+ref_set+"\tFinished writing substitions file")

		if skip_linked_SNV_counting == False:
			for readID_flag in read_map_info_dict:
				read_base_info_dict = read_map_info_dict[readID_flag][3]
				read_map_info_dict[readID_flag] = ''
				read_subIDs = []
				for loc in read_base_info_dict:
					base = read_base_info_dict[loc][0]
					try:
						sub_tup = subID_dict[str(loc)+base]
					except:
						sub_tup = ()
					if sub_tup != ():
						base_qual = read_base_info_dict[loc][1]
						map_qual = read_base_info_dict[loc][2]
						loc_in_read = read_base_info_dict[loc][4]
						loc_from_R_end = read_base_info_dict[loc][5]
						min_end_dist = min(loc_in_read,loc_from_R_end)
						# max_end_dist = max(loc_in_read,loc_from_R_end)
						if min_end_dist >= min_loc_for_counting and base_qual >= min_base_qual_for_counting:
							read_subIDs.append(str(loc)+base)
				read_subIDs = list(set(read_subIDs))
				if len(read_subIDs)>1:
					for num1 in range(0,len(read_subIDs)):
						for num2 in range(num1,len(read_subIDs)):
							if num1 != num2:
								sub1 = read_subIDs[num1]
								sub2 = read_subIDs[num2]
								try:
									loc1 = int(sub1[0:-1])
								except:
									loc1 = int(sub1.replace('A','').replace('T','').replace('C','').replace('G','').replace('N','').replace('-',''))
								try:
									loc2 = int(sub2[0:-1])
								except:
									loc2 = int(sub2.replace('A','').replace('T','').replace('C','').replace('G','').replace('N','').replace('-',''))
								if loc1 < loc2:
									pair = sub1+"_"+sub2
								elif loc1 > loc2:
									pair = sub2+"_"+sub1
								try:
									segment_isnv_link_count_dict[pair] += 1
								except:
									segment_isnv_link_count_dict[pair] = 1
			del read_map_info_dict
			if verbose == True:
				print(sampleID+"\t"+ref_set+"\tFinished collecting iSNV for lnkSNP analysis")

		if skip_qual == False:
			site_loc_outfilename = temp_dir+sampleID+"."+refID+".read_loc_hist.txt"
			read_len_outfilename = temp_dir+sampleID+"."+refID+".read_len_hist.txt"
			prop_loc_outfilename = temp_dir+sampleID+"."+refID+".read_prop_loc_hist.txt"
			indel_outfilename = temp_dir+sampleID+"."+refID+".indel_count_hist.txt"
			mismatch_outfilename = temp_dir+sampleID+"."+refID+".mismatch_count_hist.txt"
			map_qual_outfilename = temp_dir+sampleID+"."+refID+".base_map_qual.txt"

			out_site_loc_outfilename = qual_dir+refID+"/"+sampleID+".read_loc_hist.txt"
			out_read_len_outfilename = qual_dir+refID+"/"+sampleID+".read_len_hist.txt"
			out_prop_loc_outfilename = qual_dir+refID+"/"+sampleID+".read_prop_loc_hist.txt"
			out_indel_outfilename = qual_dir+refID+"/"+sampleID+".indel_count_hist.txt"
			out_mismatch_outfilename = qual_dir+refID+"/"+sampleID+".mismatch_count_hist.txt"
			out_map_qual_outfilename = qual_dir+refID+"/"+sampleID+".base_map_qual.txt"
			
			site_loc_outfile = open(site_loc_outfilename,"w")
			site_loc_outfile.write("loc\tbase\tfreq\tavg_qual\tavg_map_qual\tavg_base_loc\tavg_prop_loc\tavg_mismatch\tavg_indel")
			for count_val in range(0,300):
				site_loc_outfile.write("\t"+str(count_val))
			site_loc_outfile.write("\n")

			read_len_outfile = open(read_len_outfilename,"w")
			read_len_outfile.write("loc\tbase\tfreq\tavg_qual\tavg_map_qual\tavg_base_loc\tavg_prop_loc\tavg_mismatch\tavg_indel")
			for count_val in range(0,300):
				read_len_outfile.write("\t"+str(count_val))
			read_len_outfile.write("\n")

			prop_loc_outfile = open(prop_loc_outfilename,"w")
			prop_loc_outfile.write("loc\tbase\tfreq\tavg_qual\tavg_map_qual\tavg_base_loc\tavg_prop_loc\tavg_mismatch\tavg_indel")
			step_size = 0.01
			z = -1*step_size
			while z < 1.0:
				z += step_size
				prop_loc_outfile.write("\t"+str(round(z,3)))
			prop_loc_outfile.write("\n")

			indel_outfile = open(indel_outfilename,"w")
			indel_outfile.write("loc\tbase\tfreq\tavg_qual\tavg_map_qual\tavg_base_loc\tavg_prop_loc\tavg_mismatch\tavg_indel")
			for count_val in range(0,20):
				indel_outfile.write("\t"+str(count_val))
			indel_outfile.write("\t>=20\n")

			mismatch_outfile = open(mismatch_outfilename,"w")
			mismatch_outfile.write("loc\tbase\tfreq\tavg_qual\tavg_map_qual\tavg_base_loc\tavg_prop_loc\tavg_mismatch\tavg_indel")
			for count_val in range(0,40):
				mismatch_outfile.write("\t"+str(count_val))
			mismatch_outfile.write("\t>=40\n")

			map_qual_outfile = open(map_qual_outfilename,"w")
			map_qual_outfile.write("loc\tbase\tfreq\tavg_qual\tavg_map_qual\tavg_base_loc\tavg_prop_loc\tavg_mismatch\tavg_indel")
			for count_val in range(0,45):
				map_qual_outfile.write("\t"+str(count_val))
			map_qual_outfile.write("\n")
			for loc in range(0,max_loc):
				try:
					loc_count = count_dict[loc]
				except:
					loc_count = 0

				for i in range(0,len(bases)):
					base = bases[i]
					try:
						prop_base = round(float(len(sub_info_dict[loc][base][0]))/float(loc_count),4)
					except:
						prop_base = 0
					if prop_base >= min_proportion:
						try:
							avg_qual = round(np.average(sub_info_dict[loc][base][0]),1)
							avg_map_qual = round(np.average(sub_info_dict[loc][base][1]),1)
							avg_read_len = round(np.average(sub_info_dict[loc][base][2]),1)
							avg_base_loc = round(np.average(sub_info_dict[loc][base][3]),1)
							avg_base_loc_from_R = round(np.average(sub_info_dict[loc][base][4]),1)
							avg_prop_loc = round(np.average(sub_info_dict[loc][base][5]),1)
							avg_mismatch_count = float(round_to_n_sig_figs(np.average(sub_info_dict[loc][base][6]),2))
							avg_indel_count = float(round_to_n_sig_figs(np.average(sub_info_dict[loc][base][7]),2))
						except:
							avg_qual = 0
							avg_map_qual = 0
							avg_read_len = 0
							avg_base_loc = 0
							avg_base_loc_from_R = 0
							avg_prop_loc = 0
							avg_mismatch_count = 0
							avg_indel_count = 0

						try:
							loc_bincount = np.bincount(sub_info_dict[loc][base][3],minlength=300)
						except:
							loc_bincount =[]
							for z in range(0,300):
								loc_bincount.append(0) 
						site_loc_outfile.write(str(loc)+"\t"+base+"\t"+str(prop_base)+"\t"+str(avg_qual)+"\t"+str(avg_map_qual)+"\t"+str(avg_base_loc)+"\t"+str(avg_prop_loc)+"\t"+str(avg_mismatch_count)+"\t"+str(avg_indel_count))
						for count_val in range(0,300):
							site_loc_outfile.write("\t"+str(loc_bincount[count_val]))
						site_loc_outfile.write("\n")
						del loc_bincount

						try:
							len_bincount = np.bincount(sub_info_dict[loc][base][2],minlength=300)
						except:
							len_bincount =[]
							for z in range(0,300):
								len_bincount.append(0) 
						read_len_outfile.write(str(loc)+"\t"+base+"\t"+str(prop_base)+"\t"+str(avg_qual)+"\t"+str(avg_map_qual)+"\t"+str(avg_base_loc)+"\t"+str(avg_prop_loc)+"\t"+str(avg_mismatch_count)+"\t"+str(avg_indel_count))
						for count_val in range(0,300):
							read_len_outfile.write("\t"+str(len_bincount[count_val]))
						read_len_outfile.write("\n")
						del len_bincount
						
						step_size = 0.01
						prop_loc_bincount = {}
						prop_loc_list = sub_info_dict[loc][base][5]
						z = -1*step_size
						while z < 1.0:
							z += step_size
							counter = 0
							for prop_val in prop_loc_list:
								if prop_val >= z and prop_val < z+step_size:
									counter += 1
							prop_loc_bincount[z] = counter
						prop_loc_outfile.write(str(loc)+"\t"+base+"\t"+str(prop_base)+"\t"+str(avg_qual)+"\t"+str(avg_map_qual)+"\t"+str(avg_base_loc)+"\t"+str(avg_prop_loc)+"\t"+str(avg_mismatch_count)+"\t"+str(avg_indel_count))
						z = -1*step_size
						while z < 1.0:
							z += step_size
							prop_loc_outfile.write("\t"+str(prop_loc_bincount[z]))
						prop_loc_outfile.write("\n")
						del prop_loc_bincount


						try:
							qual_bincount = np.bincount(sub_info_dict[loc][base][1],minlength=45)
						except:
							qual_bincount = []
							for z in range(0,45):
								qual_bincount.append(0) 
						map_qual_outfile.write(str(loc)+"\t"+base+"\t"+str(prop_base)+"\t"+str(avg_qual)+"\t"+str(avg_map_qual)+"\t"+str(avg_base_loc)+"\t"+str(avg_prop_loc)+"\t"+str(avg_mismatch_count)+"\t"+str(avg_indel_count))
						for count_val in range(0,45):
							map_qual_outfile.write("\t"+str(qual_bincount[count_val]))
						map_qual_outfile.write("\n")
						del qual_bincount

						try:
							mismatch_bincount = np.bincount(sub_info_dict[loc][base][6],minlength=200)
						except:
							mismatch_bincount =[]
							for z in range(0,200):
								mismatch_bincount.append(0) 
						mismatch_outfile.write(str(loc)+"\t"+base+"\t"+str(prop_base)+"\t"+str(avg_qual)+"\t"+str(avg_map_qual)+"\t"+str(avg_base_loc)+"\t"+str(avg_prop_loc)+"\t"+str(avg_mismatch_count)+"\t"+str(avg_indel_count))
						for count_val in range(0,40):
							mismatch_outfile.write("\t"+str(mismatch_bincount[count_val]))
						greater_than_count = 0
						for count_val in range(40,200):
							try:
								greater_than_count += mismatch_bincount[count_val]
							except:
								pass
						mismatch_outfile.write("\t"+str(greater_than_count))
						mismatch_outfile.write("\n")
						del mismatch_bincount

						try:
							indel_bincount = np.bincount(sub_info_dict[loc][base][7],minlength=200)
						except:
							indel_bincount =[]
							for z in range(0,200):
								indel_bincount.append(0) 
						indel_outfile.write(str(loc)+"\t"+base+"\t"+str(prop_base)+"\t"+str(avg_qual)+"\t"+str(avg_map_qual)+"\t"+str(avg_base_loc)+"\t"+str(avg_prop_loc)+"\t"+str(avg_mismatch_count)+"\t"+str(avg_indel_count))
						for count_val in range(0,20):
							indel_outfile.write("\t"+str(indel_bincount[count_val]))
						greater_than_count = 0
						for count_val in range(20,200):
							try:
								greater_than_count += indel_bincount[count_val]
							except:
								pass
						indel_outfile.write("\t"+str(greater_than_count))
						indel_outfile.write("\n")
						del indel_bincount
			site_loc_outfile.close()
			map_qual_outfile.close()
			indel_outfile.close()
			mismatch_outfile.close()
			prop_loc_outfile.close()
			read_len_outfile.close()
			os.rename(site_loc_outfilename,out_site_loc_outfilename)
			os.rename(read_len_outfilename,out_read_len_outfilename)
			os.rename(prop_loc_outfilename,out_prop_loc_outfilename)
			os.rename(indel_outfilename,out_indel_outfilename)
			os.rename(mismatch_outfilename,out_mismatch_outfilename)
			os.rename(map_qual_outfilename,out_map_qual_outfilename)
			if verbose == True:
				print(sampleID+"\t"+ref_set+"\tFinished writing quality info summary files")

		if skip_linked_SNV_counting == False:
			lnk_temp_filename = temp_dir+sampleID+"."+refID+".lnSNVs.txt"
			lnk_out_filename = qual_dir+refID+"/"+sampleID+"."+refID+".lnSNVs.txt"
			if os.path.isfile(lnk_out_filename) == True:
				os.remove(lnk_out_filename)
			link_iSNV_outlines = 'refID\tloc1\tloc2\tsub1\tsub2\tlnk_count\tcov1\tprop1\tqual1\tmapq1\tmin_loc1\tmismatch1\tindel1\tcov2\tprop2\tqual2\tmapq2\tmin_loc2\tmismatch2\tindel2\n'
			link_count_outfile = open(lnk_temp_filename,"w")#lnk_dir
			link_count_outfile.write(link_iSNV_outlines)
			link_count_outfile.close()
			link_iSNV_outlines = ''
			nothing_written = True
			for pair in segment_isnv_link_count_dict:
				link_count = segment_isnv_link_count_dict[pair]
				isnv1 = pair.split("_")[0]
				isnv2 = pair.split("_")[1]
				try:
					sub1_tup = subID_dict[isnv1] #(loc,cov,base,prop_base,avg_qual,avg_map_qual,avg_base_loc,avg_mismatch_count,avg_indel_count)
				except:
					sub1_tup = ()
				try:
					sub2_tup = subID_dict[isnv2]
				except:
					sub2_tup = ()
				if sub1_tup != () and sub2_tup != ():
					cov1 = sub1_tup[1]
					cov2 = sub2_tup[1]
					loc1 = sub1_tup[0]
					nt1 = sub1_tup[2]
					sub_prop1 = sub1_tup[3]
					avg_qual_1 = sub1_tup[4]
					avg_map_qual_1 = sub1_tup[5]
					min_avg_loc_1 = sub1_tup[6]
					avg_mismatch_count_1 = sub1_tup[7]
					avg_indel_count_1 = sub1_tup[8]

					loc2 = sub2_tup[0]
					nt2 = sub2_tup[2]
					sub_prop2 = sub2_tup[3]
					avg_qual_2 = sub2_tup[4]
					avg_map_qual_2 = sub2_tup[5]
					min_avg_loc_2 = sub2_tup[6]
					avg_mismatch_count_2 = sub2_tup[7]
					avg_indel_count_2 = sub2_tup[8]

					link_iSNV_outlines += refID +"\t"+ str(loc1) +"\t"+ str(loc2) +"\t"+ pair.replace('_','\t') +"\t"+ str(link_count) +"\t"
					link_iSNV_outlines += str(cov1) +"\t"+ str(sub_prop1) +"\t"+ str(avg_qual_1) +"\t"+ str(avg_map_qual_1) +"\t"+ str(min_avg_loc_1) +"\t"+ str(avg_mismatch_count_1) +"\t"+ str(avg_indel_count_1) +"\t"
					link_iSNV_outlines += str(cov2) +"\t"+ str(sub_prop2) +"\t"+ str(avg_qual_2) +"\t"+ str(avg_map_qual_2) +"\t"+ str(min_avg_loc_2) +"\t"+ str(avg_mismatch_count_2) +"\t"+ str(avg_indel_count_2) +"\n"
					
					if len(link_iSNV_outlines) > 10000:
						link_count_outfile = open(lnk_temp_filename,"a")#lnk_dir
						link_count_outfile.write(link_iSNV_outlines)
						link_count_outfile.close()
						link_iSNV_outlines = ''
						nothing_written = False
			if link_iSNV_outlines != '':
				link_count_outfile = open(lnk_temp_filename,"a")#lnk_dir
				link_count_outfile.write(link_iSNV_outlines)
				link_count_outfile.close()
				os.rename(lnk_temp_filename,lnk_out_filename)
			elif nothing_written == False:
				os.remove(lnk_temp_filename)
		del sub_info_dict





#################################################################################################################################################################################

################################################################################## MAIN #########################################################################################

bases = ['A','T','C','G']

continue_running = False
forward_read_suffix = '_R1.fastq'
reverse_read_suffix = '_R2.fastq'
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
	forward_read_suffix = input('Read suffix info not found. Please enter the file suffix of the raw input reads to identify sampleIDs (or "exit" to stop):\n(e.g. for paired end reads that look like "GP01_day3_L001_R1.fastq", enter "_L001_R1.fastq")\n')
	if forward_read_suffix == "exit":
		sys.exit("Exiting.") 
	elif forward_read_suffix != '':
		continue_running = True
		outfile = open(suffix_info_file_path,"w")
		outfile.write("forward\t"+forward_read_suffix+'\n')
		outfile.write("reverse\t"+forward_read_suffix.replace('R1','R2')+'\n')
		outfile.close()

if continue_running == False:
	sys.exit("Exiting.")

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
sample_ref_pair_list = []
for sampleID in sampleID_list:
	for ref_set in ref_set_list:
		refs = ref_set.split("-")
		for ref_filename in refs:
			refID = ref_filename.split(".f")[0]
			if not os.path.exists(sub_dir+refID):
				os.makedirs(sub_dir+refID)
			if skip_qual == False:
				if not os.path.exists(qual_dir+refID):
					os.makedirs(qual_dir+refID)
		pair_tup = (sampleID,ref_set)
		sample_ref_pair_list.append(pair_tup)


if len(sample_ref_pair_list) < num_cores:
	num_cores = len(sample_ref_pair_list)

if parallel_process == False:
	for sample_ref_tup in sample_ref_pair_list:
		main_pipeline(sample_ref_tup,ref_set_list,input_read_dir,ref_dir,trim_dir,map_dir,qual_dir,sam_dir,bam_dir,temp_dir,sub_dir,command_prefix,dedupe_reads)
elif parallel_process == True:
	print(num_cores)
	processed_list = Parallel(n_jobs=num_cores)(delayed(main_pipeline)(sample_ref_tup,ref_set_list,input_read_dir,ref_dir,trim_dir,map_dir,qual_dir,sam_dir,bam_dir,temp_dir,sub_dir,command_prefix,dedupe_reads) for sample_ref_tup in sample_ref_pair_list)

