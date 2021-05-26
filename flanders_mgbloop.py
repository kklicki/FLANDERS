#!/usr/bin/env python3

import subprocess
import csv
import sys
import getopt
import pandas as pd
import os
import os.path
from os import path


def get_flanders_param(argv):
	in_file='master_clusterblast.txt' ############MUST PROVIDE FULL PATH AND FIX IT SO THAT YOU DONT HAVE TO############
	#directory containing seed clusters
	out_path = 'flanders_output'
	#directory to store neighbor analysis results
	rng = 5
	#range upstream and downstream of seed cluster to analyze
	antismash = "TRUE"
	# enable antismash search
	cores = '30'
	# default cores to run antismash search on/mnt/e5514ac8-93d5-4c49-a344-0a016f364cf1/Kevin/flanders/Multigeneblast/genbank_mf
	taxon = "bacteria"
	# default number of hits per gene to display
	signalp = "TRUE"
	# toggle signalP analysis
	mgb = "TRUE"
	#toggle multigene blast of neighbor clusters
	db = '/mnt/e5514ac8-93d5-4c49-a344-0a016f364cf1/Kevin/flanders/Multigeneblast/genbank_mf' ######## NEED TO UNHARDWIRE THIS ALSO##########
	# default database for multigeneblast to search
	hitspergene = '250'
	# default number of hits per gene to display
	minseqcov = '25'
	#default minimum sequence coverage for hit
	minpercid = '30'
	#default minimum percent identity for hit
	distancekb = '2'
	#maximum kilobases between two blast hits to be considered contiguous
	syntenyweight = '0.5'
	#adjust weight given to synteny conservation in BLAST scoring
	
	

	argv= sys.argv[1:]
	try:
		opts, args = getopt.getopt(argv,"hi:o:r:a:c:t:s:m:d:g:s:p:b:w:",["help","in","out","range","antismash","cores","taxon","signalp","mgb","db","hitspergene","minseqcov","minpercid","distancekb","syntenyweight"])
	except getopt.GetoptError:
			print("usage:\n$flanders_mgbloop.py -i <input file> -o <output parent directory> -r <range of genes upstream and downstream to analyze> -a <perform antismash analysis> -c <threads for antismash computation> -t <taxon for antismash (bacteria or fungi)> -s <perform signalP analysis>")
			sys.exit()
			#looks for arguments in command line in the stated order where the directory to the input fasta is supplied first and subsequent options and arguments follow
	if len(sys.argv) == 1:
		print("usage:\n$flanders_mgbloop.py -i <input directory> -o <output parent directory> -r <range of genes upstream and downstream to analyze> -a <perform antismash analysis> -c <threads for antismash computation> -t <taxon for antismash (bacteria or fungi)> -s <perform signalP analysis>")
		sys.exit()
	else:
		pass
	for opt, arg in opts:
		if opt in ('-h',"--help",):
			print("usage:\n$flanders_mgbloop.py -i <input directory> -o <output directory> -r <range of genes upstream and downstream to analyze> -a <perform antismash analysis> -c <threads for antismash computation> -t <taxon for antismash (bacteria or fungi)> -s <perform signalP analysis>")
			sys.exit()
		elif opt in ("-i","--in"):
			in_file = arg
		elif opt in ("-o","--out"):
			out_path = arg
		elif opt in ("-r","--range"):
			rng = str(int(arg))
		elif opt in ("-a","--antismash"):
			antismash = "TRUE"
		elif opt in ("-c","--cores"):
			cores = str(int(arg))
		elif opt in ("-t","--taxon"):
			taxon = arg
		elif opt in ("-s","--signalp"):
			signalP = "TRUE"
		elif opt in ("-m","--mgb"):
			mgb = "TRUE"
		elif opt in ("-d","--db"):
			db = arg
		elif opt in ("-g","--hitspergene"):
			hitspergene = str(int(arg))
		elif opt in ("-s","--minseqcov"):
			minseqcov = str(int(arg))
		elif opt in ("-p","--minpercid"):
			minpercid = str(int(arg))
		elif opt in ("-b","--distancekb"):
			distancekb = str(int(arg))
		elif opt in ("-w","--syntenyweight"):
			if 0 <= float(arg) <= 1:
				syntenyweight = str(float(arg))
			else:
				print("Error: synteny weight must be between 0 and 1")
				sys.exit()
		
	flanders_param = [in_file, out_path, rng, antismash, cores, taxon, signalp, mgb, db, hitspergene, minseqcov, minpercid, distancekb, syntenyweight]
	#once all arguments have been collected from command line input, pass them to the rest of the functions as flanders_param
	return flanders_param



def get_neighborhood(flanders_param, line, flanders_cluster):
#given seed cluster input, generate accession numbers for genes upstream and downstream of seed cluster witin given range
	genome = ""
	fix = []
	liz = []
	line = line.rstrip("})")
	line = line.lstrip("frozenset({")
	line = line.replace("\'","")
	fix = line.split(", ")
	#fixing output from MGB_loop output so that its easier to parse
	for n in fix:
		n = n.lstrip()
		liz.append(n)
	for n in liz:
		if "_" in n:
			genome = n
			del liz[liz.index(n)]
		else:
			pass
		#Getting genome accession number from MGB_loop output
	liz.sort()
	#putting seed cluster gene accession numbers in order
	#print(liz)
	genome = genome.split("_",1)[0]
	print(genome)
	num = []
	up = []
	down = []
	for f in liz:
		f_abc = str(f[0:3])
		f_int=int(f[3:8])
		num.append(f_int)
	num.sort()
	up = []
	upstream = []
	upstream_input = genome+"_upstream"
	down = []
	downstream = []
	downstream_input = genome+"_downstream"
	dis = int(flanders_param[2])
	while dis > 0:
		q=num[0] - dis
		up.append(q)
		r=num[-1] + dis
		down.append(r)
		dis = dis - 1
	for u in up:
		upstream.append(f_abc + str(u))
	for v in down:
		downstream.append(f_abc + str(v))
		downstream.sort()
	#generate accession numbers for upstream and downstream clusters
	flanders_cluster = flanders_cluster.append({"Genome":genome,"Location":"upstream","Accession":upstream},ignore_index=True)
	flanders_cluster = flanders_cluster.append({"Genome":genome,"Location":"seed","Accession":liz},ignore_index=True)
	flanders_cluster = flanders_cluster.append({"Genome":genome,"Location":"downstream","Accession":downstream},ignore_index=True)
	#add upstream, seed, and downstream cluster gene accession numbers into dataframe
	#print(flanders_cluster)
	return(flanders_cluster, genome, upstream, downstream)

def get_genome_region(flanders_param, genome, upstream, downstream):
	#in the event that nucleotide sequence is needed for antismash/multigeneblast, retrieve nucleotide sequence from given 
	subprocess.run(["ncbi-acc-download","-m", "nucleotide","-F", "gff3", genome])
	#download gff3 of current genome
	with open(genome+".gff","r") as tsv:
		gff = csv.reader(tsv, delimiter='\t')
		for row in gff:
			if row[0].startswith("#"):
				pass
			elif upstream[0] in row[8]:
				upstream_cluster_start = int(row[3])
			elif upstream[-1] in row[8]:
				upstream_cluster_end = int(row[4])
			elif downstream[0] in row[8]:
				downstream_cluster_start = int(row[3])
			elif downstream[-1] in row[8]:
				downstream_cluster_end = int(row[4])
				break
			else:
				continue
			#find the start and end coordinates for the first and last gene in the upstream and downstream cluster respectively
	subprocess.run(["ncbi-acc-download","-m", "nucleotide","-F", "fasta", genome])
	#download the fasta file corresponding to current genome
	genome_seq = ""
	cluster_seq = ""
	with open(genome+".fa","r") as fasta:
		for line in fasta:
			if line.startswith(">"):
				title = line
				continue
			else:
				genome_seq += line.replace("\n","")
	upstream_cluster_seq=genome_seq[upstream_cluster_start-1:upstream_cluster_end-1]
	downstream_cluster_seq=genome_seq[downstream_cluster_start-1:downstream_cluster_end-1]
	#put whole genome sequence into a string and then pull out the cluster sequences based on their coordinates
	with open(genome+"_upstream_cluster.fasta","w") as w:
		w.write(title)
		w.write(upstream_cluster_seq)
	with open(genome+"_downstream_cluster.fasta","w") as w:
		w.write(title)
		w.write(downstream_cluster_seq)
	upstream_input = ""
	upstream_input = " ".join(upstream)
	downstream_input = ""
	downstream_input = " ".join(downstream)
	#prepare input for ncbi-acc-download, retrieve upstream and downstream cluster sequences
	subprocess.run(["ncbi-acc-download","-m","protein","-F","fasta","-o",genome+"_upstream.fasta",upstream_input])
	subprocess.run(["ncbi-acc-download","-m","protein","-F","fasta","-o",genome+"_downstream.fasta",downstream_input])
	return()

def clusterblast_neighbors(flanders_param, flanders_cluster, genome, upstream, downstream):
	subprocess.run(["python","/mnt/e5514ac8-93d5-4c49-a344-0a016f364cf1/Kevin/flanders/Multigeneblast/multigeneblast_py3.py","-in",genome +"_upstream.fasta", "-out",genome +"_upstream_clusterblast","-db",flanders_param[8],"-cores",flanders_param[4],"-hitspergene", flanders_param[9],'-minseqcov', flanders_param[10], '-minpercid', flanders_param[11], "-distancekb", flanders_param[12], "-syntenyweight", flanders_param[13], "-outpages", "1"])
	#NEED TO UNHARDWIRE FULL PATH TO multigeneblast.py
	os.chdir(genome+"_upstream_clusterblast")
	iteration_hits=[]
	i=0
	if path.exists("clusterblast_output.txt"):
		pass
	else:
		os.chdir("..")
		return(flanders_cluster)
	with open("clusterblast_output.txt","r") as doc:
		for line in doc:
			if "Significant hits:" in line:
				break
			else:
				i=i+1
	cluster_length=i-5
	#because there are 5 other lines at the beginning of the clusterblast output file besides the query cluster accessions and annotations
	with open("clusterblast_output.txt","r") as file: #open results file
		title = file.readlines(1) #get the file name
		head,tail=file.read().split("Details:",2)# split the file so its easy to parse the 'details'
		temp = open("details.txt", "w")#make its own details file
		temp.write(tail)#write the details portion to its own file
	with open("details.txt","r") as details:
		deets=[]
		for line in details:
			deets.append(line)
		for i in range(len(deets)):
			hits = []#initialize list of hits
			if ">>" in deets[i]:
				g=deets[i+2]
			else:
				pass
			if "Table of Blast hits (query gene, subject gene, %identity, blast score, %coverage, e-value):\n" in deets[i]:#look for each instance of this string, this string appears at 250 instances in the input file, so not sure if this is where the problem is
				n=1
				#hits.append(genome)
				while n <= int(cluster_length):# making sure its a full cluster hit (all query genes have homologs found)        
					if((i+n) < len(deets)):# checking that it doesnt exceed the number of items in deets
					# the split makes the string into a list and the [:-1] gets rid of the last item which was just '\n'
						row = deets[i+n].split('\t')[:-1] 
						hits.append(row) #use the index of the 'table of blast hits' item to get the next n items from the file list and put it on the hit list
						n += 1
				if [] not in hits:
					saturation_point = hits
					#print(hits)
					acc_list = []
					#check that its a full cluster before adding it to the full clusters list
					for h in hits:
						acc_list.append(h[1])
					space = g.find(" ")
					#return the index of the space between genome and gene accession
					g = g.replace(g[0:space],"")
					g = g.replace("\n","")
					acc_list.append(g)
					iteration_hits.append(frozenset(acc_list))
					continue
				else:
					break
			else:
				continue
		hitz = len(iteration_hits)
		# number of full clusters found by multigene blast run
		for index, r in flanders_cluster.iterrows():
		#print(index,r["Location"]) # r['Location'] points to cell in row r at column location
			if "upstream" in r["Location"]:
				r["Clusterblast_hits"] = hitz
	os.chdir("..")
	#go back to parent directory and do everything again for the downstream cluster
	subprocess.run(["python","/mnt/e5514ac8-93d5-4c49-a344-0a016f364cf1/Kevin/flanders/Multigeneblast/multigeneblast_py3.py","-in",genome +"_downstream.fasta", "-out",genome +"_downstream_clusterblast","-db",flanders_param[8],"-cores",flanders_param[4],"-hitspergene", flanders_param[9],'-minseqcov', flanders_param[10], '-minpercid', flanders_param[11], "-distancekb", flanders_param[12], "-syntenyweight", flanders_param[13], "-outpages", "1"])	#NEED TO UNHARDWIRE FULL PATH TO multigeneblast.py
	os.chdir(genome+"_downstream_clusterblast")
	iteration_hits=[]
	i=0
	if path.exists("clusterblast_output.txt"):
		pass
	else:
		os.chdir("..")
		return(flanders_cluster)
	with open("clusterblast_output.txt","r") as doc:
		for line in doc:
			if "Significant hits:" in line:
				break
			else:
				i=i+1
	cluster_length=i-5
	#because there are 5 other lines at the beginning of the clusterblast output file besides the query cluster accessions and annotations
	with open("clusterblast_output.txt","r") as file: #open results file
		title = file.readlines(1) #get the file name
		head,tail=file.read().split("Details:",2)# split the file so its easy to parse the 'details'
		temp = open("details.txt", "w")#make its own details file
		temp.write(tail)#write the details portion to its own file
	with open("details.txt","r") as details:
		deets=[]
		for line in details:
			deets.append(line)
		for i in range(len(deets)):
			hits = []#initialize list of hits
			if ">>" in deets[i]:
				g=deets[i+2]
			else:
				pass
			if "Table of Blast hits (query gene, subject gene, %identity, blast score, %coverage, e-value):\n" in deets[i]:#look for each instance of this string, this string appears at 250 instances in the input file, so not sure if this is where the problem is
				n=1
				#hits.append(genome)
				while n <= int(cluster_length):#6 because there are 6 genes in full ebo cluster, can make this into a provided argument           
					if((i+n) < len(deets)):# checking that it doesnt exceed the number of items in deets
					# I helped you split it based on tabs so that you can more easily get the stuff from the row
					# the split makes the string into a list and the [:-1] gets rid of the last item which was just '\n'
						row = deets[i+n].split('\t')[:-1] 
						hits.append(row) #use the index of the 'table of blast hits' item to get the next 6 items from the file list and put it on the hit list
						n += 1
				if [] not in hits:
					saturation_point = hits
					#print(hits)
					acc_list = []
					for h in hits:
						acc_list.append(h[1])
					space = g.find(" ")
					g = g.replace(g[0:space],"")
					g = g.replace("\n","")
					acc_list.append(g)
					iteration_hits.append(frozenset(acc_list))
					continue
				else:
					break
			else:
				continue
		hitz = len(iteration_hits)
		for index, r in flanders_cluster.iterrows():
		#print(index,r["Location"]) # r['Location'] points to cell in row r at column location
			if "downstream" in r["Location"]:
				r["Clusterblast_hits"] = hitz
	os.chdir("..")
	return(flanders_cluster)



def run_signalp(flanders_param, flanders_cluster, genome, upstream, downstream):
	upstream_input = ""
	upstream_input = " ".join(upstream)
	downstream_input = ""
	downstream_input = " ".join(downstream)
	#prepare input for ncbi-acc-download, retrieve upstream and downstream cluster sequences
	subprocess.run(["ncbi-acc-download","-m","protein","-F","fasta","-o",genome+"_upstream.fasta",upstream_input])
	subprocess.run(["ncbi-acc-download","-m","protein","-F","fasta","-o",genome+"_downstream.fasta",downstream_input])
	#convert upstream and donwstream accession lists to strings to pass to signalP input
	subprocess.run(["signalp", "-format","short", "-org", "gram-","-fasta", genome+"_upstream.fasta"])
	subprocess.run(["signalp", "-format","short", "-org", "gram-","-fasta", genome+"_downstream.fasta"])
	with open(genome + "_upstream_summary.signalp5","r") as s:
		peps = []
		s = csv.reader(s,delimiter="\t")
		#convert signalP output file to readable format
		for row in s:
			if "#" in row[0]:
				pass
			else:
				peps.append(row[0:2])
		#gather each gene and its signal peptide designation
		c = 0
		for p in peps:
			if "OTHER" in p:
				pass
			else:
				c += 1
		for index, r in flanders_cluster.iterrows():
			#print(index,r["Location"]) # r['Location'] points to cell in row r at column location
			if "upstream" in r["Location"]:
				r["Signal_peptide"] = peps
				r["Signal_peptide_count"] = str(c)
			else:
				continue
			#
	with open(genome + "_downstream_summary.signalp5","r") as s:
		peps = []
		s = csv.reader(s,delimiter="\t")
		for row in s:
			if "#" in row[0]:
				pass
			else:
				peps.append(row[0:2])
		c = 0
		for p in peps:
			if "OTHER" in p:
				pass
			else:
				c += 1
		for index, r in flanders_cluster.iterrows():
			#print(index,r["Location"]) # r['Location'] points to cell in row r at column location
			if "downstream" in r["Location"]:
				r["Signal_peptide"] = peps
				r["Signal_peptide_count"] = str(c)
			else:
				continue
	return(flanders_cluster)



def run_antismash(flanders_param, flanders_cluster, upstream, downstream, genome):
	subprocess.run(["antismash","-c",flanders_param[4], "--taxon", "bacteria","--genefinding-tool","prodigal", genome+"_upstream_cluster.fasta"])
	#run antismash on the upstream region
	os.chdir(genome+"_upstream_cluster")
	#change into antismash output directory
	smc = []#initialize lits of secondary metabolite clusters within upstream region
	with open("index.html","r") as h:
		#using the html file
		head,tail=h.read().split("<div style=\"display: flex; flex-wrap: wrap\">",2)
		mane,rest = tail.split("<!-- overview page -->",2)
		mane=mane.split("\n")
		#pull out the part delineating the buttons the html page will display
		for i in mane:
			if i.lstrip().startswith("<div class"):
				clust = i.lstrip().split()[2]
				smc.append(clust)
			else:
				continue
		#parse out smc names
	for index, r in flanders_cluster.iterrows():
	#print(index,r["Location"]) # r['Location'] points to cell in row r at column location
		if "upstream" in r["Location"]:
			r["BGC_Type"] = smc
		else:
			continue
	#update flanders cluster dataframe with smc names
	os.chdir("..")
	#go back and do it all again for the downstream cluster
	subprocess.run(["antismash","-c",flanders_param[4],"--taxon","bacteria","--genefinding-tool","prodigal", genome+"_downstream_cluster.fasta"])#figure out antismash output and stuff
	os.chdir(genome+"_downstream_cluster")
	smc = []
	with open("index.html","r") as h:
		head,tail=h.read().split("<div style=\"display: flex; flex-wrap: wrap\">",2)
		mane,rest = tail.split("<!-- overview page -->",2)
		mane=mane.split("\n")
		for i in mane:
			if i.lstrip().startswith("<div class"):
				clust = i.lstrip().split()[2]
				smc.append(clust)
			else:
				continue
		print(smc)
	for index, r in flanders_cluster.iterrows():
	#print(index,r["Location"]) # r['Location'] points to cell in row r at column location
		if "downstream" in r["Location"]:
			r["BGC_Type"] = smc
		else:
			continue
	os.chdir("..")
	return(flanders_cluster)

def main(argv):
	argv = sys.argv[1:]
	#parse provided arguments into a list
	flanders_param = get_flanders_param(argv)
	print(flanders_param)
	if path.exists(flanders_param[1]):
		os.chdir(flanders_param[1])
	else:
		os.mkdir(flanders_param[1])
		os.chdir(flanders_param[1])
	flanders_output = pd.DataFrame(columns=["Genome","Location","Accession","BGC_Type","Signal_peptide","Signal_peptide_count","Clusterblast_hits"])
	#intialize master flanders output dataframe
	with open(flanders_param[0],"r") as input:
		for l in input:
			flanders_cluster = pd.DataFrame(columns=["Genome","Location","Accession","BGC_Type","Signal_peptide","Signal_peptide_count","Clusterblast_hits"])
			flanders_cluster, genome, upstream, downstream = get_neighborhood(flanders_param,l,flanders_cluster)
			if path.exists(genome):
				continue
			else:
				pass
			os.mkdir(genome)
			os.chdir(genome)
			if flanders_param[3] == "TRUE":
				#get_genome_region(flanders_param,genome,upstream,downstream)
				print("running antismash")
				#flanders_cluster=run_antismash(flanders_param, flanders_cluster, upstream, downstream, genome)
			else:
				continue
			if flanders_param[6] == "TRUE":
				flanders_cluster = run_signalp(flanders_param, flanders_cluster, genome, upstream, downstream)
			else:
				continue
			if flanders_param[7] == "TRUE":
				flanders_cluster = clusterblast_neighbors(flanders_param,flanders_cluster,genome,upstream,downstream)
			else:
				continue
			flanders_cluster.to_csv(genome+".tsv", sep = "\t")
			flanders_output = pd.concat([flanders_output,flanders_cluster],axis = 0)
			os.chdir("..")
	flanders_output.to_csv("flanders_output.tsv", sep = '\t')



if __name__ == "__main__":
	main(sys.argv[1:])


