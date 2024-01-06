import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
import pandas
from scipy.stats import gaussian_kde


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

ref_set_list = ['A_HA_H1','A_HA_H3','A_M_1','A_NA_N1','A_NA_N2','A_NP_1','A_NS','A_PA_1','A_PB1_1','A_PB2_1','B_HA','B_M','B_NA','B_NP','B_NS','B_PA','B_PB1','B_PB2']
# ref_set_list = ['PB2.fa','PB1.fa','PA.fa','NP.fa','M.fa','NS.fa','HA.fa','NA.fa']
all_segment_list = ref_set_list

min_frequency = 0.02
min_cov = 100
min_minor_qual = 30

num_plots = 5
x = {}
y = {}
for r in range(0,num_plots+1):
	x[r] = {}
	y[r] = {}
	for i in range(0,len(all_segment_list)):
		segment = all_segment_list[i].split(".")[0]
		x[r][segment] = []
		y[r][segment] = []

infile = open("all_sites.base_info."+output_suffix+".txt","r")
first_line = True
for line in infile:
	if first_line == True:
		first_line = False
	else:
		line = line.strip().split("\t")
		sampleID = line[0]
		segment = line[1]
		
		loc = int(line[2])
		cov = int(line[3])
		
		major_prop = float(line[6])
		major_qual = float(line[7])
		major_mapq = float(line[8])
		major_loc = float(line[10])
		major_mismatch = float(line[13])
		major_indel = float(line[14])

		minor_prop = float(line[15])
		minor_qual = float(line[16])
		minor_mapq = float(line[17])
		minor_loc = float(line[19])
		minor_mismatch = float(line[22])
		minor_indel = float(line[23])
		
		# major_prop = float(line[7])
		# major_mapq = float(line[11])
		# major_loc = float(line[15])
		# major_mismatch = float(line[21])
		# major_indel = float(line[23])

		xval = major_mapq
		yval = major_loc
		r = 0
		x[r][segment].append(xval)
		y[r][segment].append(yval)

		xval = major_mapq
		yval = major_qual
		r = 1
		x[r][segment].append(xval)
		y[r][segment].append(yval)

		xval = major_mapq
		yval = major_mismatch
		r = 2
		x[r][segment].append(xval)
		y[r][segment].append(yval)

		xval = major_mapq
		yval = major_indel
		r = 3
		x[r][segment].append(xval)
		y[r][segment].append(yval)

		if minor_prop >= min_frequency and minor_qual >= min_minor_qual and cov >= min_cov:
			xval = minor_mapq
			yval = minor_loc
			r = 4
			x[r][segment].append(xval)
			y[r][segment].append(yval)
infile.close()


for r in range(0,num_plots+1):
	lineout = ''
	for i in range(0,len(all_segment_list)):
		segment = all_segment_list[i].split(".")[0]
		lineout += "\t"+str(len(x[r][segment]))
	print(lineout)


fig, axs = plt.subplots(num_plots,len(all_segment_list), sharex=False, sharey=False)
fig.set_size_inches(len(all_segment_list)*4.0,num_plots*3.5)

axis_labels = {0:("mapq","read_loc",(20,45),(0,40)),1:("mapq","qual",(20,45),(20,40)),2:("mapq","mismatch",(20,45),(0,20)),3:("mapq","indel",(20,45),(0,100)),4:("minor_mapq","minor_read_loc",(10,45),(0,50))}

for r in range(0,num_plots):
	row = r
	for i in range(0,len(all_segment_list)):
		segment = all_segment_list[i].split(".")[0]
		if len(x[r][segment]) >0:
			if len(x[r][segment]) > 10000:
				x[r][segment] = x[r][segment][0:10000]
				y[r][segment] = y[r][segment][0:10000]
			print(str(row)+"\t"+str(i))
			xy = np.vstack([x[r][segment],y[r][segment]])
			try:
				z = gaussian_kde(xy)(xy)
				axs[row, i].scatter(x[r][segment],y[r][segment], c=z, s=3)
			except:
				axs[row, i].scatter(x[r][segment],y[r][segment], s=3)

			axs[row, i].set_xlabel(axis_labels[r][0])
			axs[row, i].set_ylabel(axis_labels[r][1])
			if row == 0:
				axs[row, i].set_title(segment)

			axs[row,i].set_xlim([axis_labels[r][2][0], axis_labels[r][2][1]])
			axs[row,i].set_ylim([axis_labels[r][3][0], axis_labels[r][3][1]])
plt.tight_layout()
plt.savefig("qual_summary.major_allele."+str(output_suffix)+".png")


