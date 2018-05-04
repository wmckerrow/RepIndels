"""
Copyright Wilson McKerrow, 2018
This file is part of RepIndel.
RepIndel is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
RepIndel is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with RepIndel.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy
import json
import numpy.random
import sys
from Bio import SeqIO
import copy
import datetime
from scipy.stats import norm
from sets import Set
import PEhmm2
from multiprocessing import Pool
import pickle

sys.settrace

reads_to_check = []

# !!! Set parameters

P_ERROR = 0.001 # Probability of a given sequencing error (e.g. A->C)
INDEL_OPEN = 0.001 # Probability of an indel
DEL_EXTEND = 0.999 # Probability of extending in existing deletion
INS_EXTEND = 0.9 # Probability of extending in existing insertion
DEL_END = 0.5 # If DEL_EXTEND falls below this value cut off the deletion at this position for the purpose of seeding alignments

MINP = 0.001 # Minimum value for parameters
MIN_TO_INCLUDE = 0.001 # Should be same value as above

READ_LEN = 100 # Length of a single read end

MIN_READ_PROB = P_ERROR**4*(1-P_ERROR)**96 # Add this to read likelihood so reads that align poorly are ignored

KMER_LEN = 19 # Seed length
OVERHANG = 25 # Look OVERHANG + OVERHANG_FACTOR*(# unaligned read letters on this end) position out when performing forward/backward
OVERHANG_FACTOR = 2

THREADS = 16 # Number of threads to use when parallelizing alignment. Should all be on one node.

INT_NB_SIZE = 104 # m parameter for negative binomial describing length of sequenced positions between read ends
INT_NB_P = 0.38 # p parameter for that read end
MAX_INSERT_PDF = 0.019 # maximum of that negative binomial pmf

# List the names of repeat elements to update. If not listed it is assumed that the reference is correct.
TO_UPDATE = ['source','source1','source6','source7','source13','chr3L:20361741','chr3R:22206092','chr3L:5994793','POGO_1','POGO_2','POGO_3','POGO_4','POGO_5','POGO_6','POGO_7','POGO_8','POGO_9','POGO_10','POGO_11','POGO_12','POGO_13','POGO_14','POGO_15','POGO_16','POGO_17','POGO_18','POGO_19','POGO_20','POGO_21','POGO_22','POGO_23','POGO_24','POGO_25','POGO_26','POGO_27','POGO_28']

# !!! End of parameters to set

class PEreadclass(object):
	def __init__(self, id, seq1,seq2):
		self.seq1 = seq1
		self.seq2 = seq2
		self.seqrc1 = rev_comp(seq1)
		self.seqrc2 = rev_comp(seq2)
		self.id = id.split('/')[0]
		self.alnranges_1 = list()
		self.alnranges_2 = list()
		self.alnranges_rc1 = list()
		self.alnranges_rc2 = list()
		self.prob_from_fixed = MIN_READ_PROB*MIN_READ_PROB*MAX_INSERT_PDF

class hmmclass(object):
	def __init__(self,genome_profile,del_open,del_extend,ins_open,ins_extend):
		self.genome_profile = genome_profile
		self.del_open = del_open
		self.del_extend = del_extend
		self.ins_open = ins_open
		self.ins_extend = ins_extend
		self.weight = numpy.ones(len(del_open))

# Convert sequences from string to numpy array.
def seq2array(seq):
	nuc2num = {'A':0,'C':1,'G':2,'T':3, 'N':4}
	return numpy.array([nuc2num[c] for c in seq],dtype=numpy.int8)

# Convert sequences from numpy array to string.
def array2seq(seq_array):
	seq = ''
	for letter in seq_array:
		seq += 'ACGTN'[letter]
	return seq

# Store sequence in a fasta file a dictionary of numpy arrays
def fasta_as_array(fasta_file):
	fasta_dict = SeqIO.to_dict(SeqIO.parse(open(fasta_file,'r'), "fasta"))
	for seq in fasta_dict:
		fasta_dict[seq] = seq2array(fasta_dict[seq].upper())
	return fasta_dict

# Get reverse complement of a numpy array sequence	
def rev_comp(seq_array):
	revcompseq = 3 - seq_array[::-1]
	revcompseq[revcompseq<0] = 4
	return revcompseq

# Build an initial profile guess for EM based on the reference genome. Each row is the provided probs vector translated mod 4 so that the first element of the probs vector is in the reference base position.	
def create_initial_guess(ref,probs):
	genome_profile = dict()
	for seq_name in ref:
		genome_profile[seq_name] = numpy.zeros((len(ref[seq_name]),4))
		for i in range(len(ref[seq_name])):
			if ref[seq_name][i] != 4:
				genome_profile[seq_name][i,] = probs[(numpy.arange(4)-ref[seq_name][i])%4]
			else:
				genome_profile[seq_name][i,] = numpy.array([0.25,0.25,0.25,0.25])
	return genome_profile

def forward_sum(sequence,genome_profile,del_open,del_extend,ins_open,ins_extend,fm_init,fi_init):
	
	fm = numpy.zeros((len(genome_profile),len(sequence)))
	fd = numpy.zeros((len(genome_profile),len(sequence)))
	fi = numpy.zeros((len(genome_profile),len(sequence)))
	
	fm[:,0] = fm_init*genome_profile[:,sequence[0]]
	fi[:,0] = fi_init*0.25
	
	for j in range(1,len(sequence)):
		fi[0,j] = 0.25*(fm[0,j-1]*ins_open[0] + fi[0,j-1]*ins_extend[0])
	for i in range(1,len(genome_profile)):
		fd[i,0] = fm[i-1,0]*del_open[i] + fd[i-1,0]*del_extend[i]
		for j in range(1,len(sequence)):
			fm[i,j] = fm[i-1,j-1]*(1-del_open[i]-ins_open[i-1]) + fd[i-1,j-1]*(1.0-del_extend[i]) + fi[i-1,j-1]*(1.0-ins_extend[i-1])
			fm[i,j] *= genome_profile[i,sequence[j]]
			fd[i,j] = fm[i-1,j]*del_open[i] + fd[i-1,j]*del_extend[i]
			fi[i,j] = 0.25*(fm[i,j-1]*ins_open[i] + fi[i,j-1]*ins_extend[i])
	
	return fm,fd,fi

def forward_sum_interior(fm_init,fi_init,del_open,del_extend,ins_open,ins_extend,size,p):
	
	fm = numpy.zeros((len(del_open),size+1))
	fd = numpy.zeros((len(del_open),size+1))
	fi = numpy.zeros((len(del_open),size+1))
	fm[:,0] = fm_init
	fi[:,0] = fi_init
	
	fi[0,1] = fm[0,1]*(1-p)*(ins_open[0]) + fi[0,0]
	for i in range(1,len(del_open)):
		fm[i,1] = fm[i-1,1]*(1-p)*(1-del_open[i]-ins_open[i-1]) + fd[i-1,1]*(1-del_extend[i]) + fi[i-1,1]*(1-ins_extend[i-1]) + fm[i-1,0]
		fd[i,1] = fm[i-1,1]*(1-p)*del_open[i] + fd[i-1,1]*del_extend[i]
		fi[i,1] = fm[i,1]*(1-p)*ins_open[i] + fi[i,0]
	
	for k in range(2,size+1):
		fi[0,k] = fm[0,k]*(1-p)*(ins_open[0]) + ins_extend[0]*fi[0,k-1]
		for i in range(1,len(del_open)):
			fm[i,k] = fm[i-1,k]*(1-p)*(1-del_open[i]-ins_open[i-1]) + fd[i-1,k]*(1-del_extend[i]) + fi[i-1,k]*(1-ins_extend[i-1]) + p*fm[i-1,k-1]
			fd[i,k] = fm[i-1,k]*(1-p)*del_open[i] + fd[i-1,k]*del_extend[i]
			fi[i,k] = fm[i,k]*(1-p)*ins_open[i] + ins_extend[i]*fi[i,k-1]
#			print i,k,fm[i,k],fi[i,k],fd[i,k]
	
	return fm[:,1:],fd[:,1:],fi[:,1:]

def backward_sum(sequence,genome_profile,del_open,del_extend,ins_open,ins_extend,bm_init,bi_init):
	bm = numpy.zeros((len(genome_profile),len(sequence)))
	bd = numpy.zeros((len(genome_profile),len(sequence)))
	bi = numpy.zeros((len(genome_profile),len(sequence)))
	
	bm[:,-1] = bm_init
	bi[:,-1] = bi_init
		
	for j in range(len(sequence)-1)[::-1]:
		bm[-1,j] = bi[-1,j+1]*ins_open[-1]*0.25
		bi[-1,j] = bi[-1,j+1]*ins_extend[-1]*0.25
	for i in range(len(genome_profile)-1)[::-1]:
		for j in range(len(sequence)-1)[::-1]:
			bm[i,j] = bm[i+1,j+1]*(1-del_open[i+1]-ins_open[i])*genome_profile[i+1,sequence[j+1]] + bd[i+1,j]*del_open[i+1] + bi[i,j+1]*ins_open[i]*0.25
			bd[i,j] = bm[i+1,j+1]*(1.0-del_extend[i+1])*genome_profile[i+1,sequence[j+1]] + bd[i+1,j]*del_extend[i+1]
			bi[i,j] = bm[i+1,j+1]*(1.0-ins_extend[i])*genome_profile[i+1,sequence[j+1]] + bi[i,j+1]*ins_extend[i]*0.25
			
	return bm,bd,bi
	
def backward_sum_interior(bm_init,bi_init,del_open,del_extend,ins_open,ins_extend,size,p):
	bm = numpy.zeros((len(del_open),size+1))
	bd = numpy.zeros((len(del_open),size+1))
	bi = numpy.zeros((len(del_open),size+1))
	
	bm[:,-1] = bm_init
	bi[:,-1] = bi_init
		
	for k in range(size)[::-1]:
		bm[-1,k] = bi[-1,k]*(1-p)*ins_open[-1]
		bi[-1,k] = ins_extend[-1]*bi[-1,k+1]
		for i in range(len(del_open)-1)[::-1]:
			bi[i,k] = bm[i+1,k]*(1.0-ins_extend[i]) + ins_extend[i]*bi[i,k+1]
			bm[i,k] = bm[i+1,k]*(1.0-p)*(1-del_open[i+1]-ins_open[i]) + bd[i+1,k]*(1-p)*del_open[i+1] + bi[i,k]*(1-p)*ins_open[i] + p*bm[i+1,k+1]
			bd[i,k] = bm[i+1,k]*(1.0-del_extend[i+1]) + bd[i+1,k]*del_extend[i+1]
	return bm[:,:-1],bd[:,:-1],bi[:,:-1]
	
def forward_sum_pair(seq1,ranges1,seq2,ranges2,hmm):
	
	geno_len = dict()
	fm_int_start = dict()
	fi_int_start = dict()
	for name in hmm:
		geno_len[name] = len(hmm[name].genome_profile)
		fm_int_start[name] = numpy.zeros(geno_len[name])#+MIN_READ_PROB
		fi_int_start[name] = numpy.zeros(geno_len[name])
		
	range_of_alns = dict()
	for alnrange in ranges1+ranges2:
		name,start,stop = alnrange
		if name not in range_of_alns:
			range_of_alns[name] = [start-1,stop+1]
		else:
			range_of_alns[name][0] = min(start-1,range_of_alns[name][0])
			range_of_alns[name][1] = max(stop+1,range_of_alns[name][1])
		
	fm_1 = list()
	fd_1 = list()
	fi_1 = list()
	for alnrange in ranges1:
		name,start,stop = alnrange
		fm = numpy.zeros((stop-start,len(seq1)))
		fd = numpy.zeros((stop-start,len(seq1)))
		fi = numpy.zeros((stop-start,len(seq1)))
#		fm,fd,fi = forward_sum(seq1,hmm[name].genome_profile[start:stop],hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],1.0-hmm[name].ins_open[start:stop],hmm[name].ins_open[start:stop])
		x = PEhmm2.forward_sum(seq1,hmm[name].genome_profile[start:stop],hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],(1.0-hmm[name].ins_open[start:stop])*hmm[name].weight[start:stop],hmm[name].ins_open[start:stop]*hmm[name].weight[start:stop],fm,fd,fi)
		if x:
			print 'error: forward sum seq1:', x
			print name,start,stop
			pickle.dump([seq1,hmm[name].genome_profile[start:stop],hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],hmm[name].weight,fm,fd,fi],open('error.pkl','w'))
			exit()
		fm_1.append(fm)
		fd_1.append(fd)
		fi_1.append(fi)
		fm_int_start[name][start:stop] += fm[:,-1]
		fi_int_start[name][start:stop] += fi[:,-1]
		
#		print 'seq1',name,start,stop,numpy.sum(fm[:,-1])
						
	fm_int = dict()
	fd_int = dict()
	fi_int = dict()
	for name in range_of_alns:
		fm_int[name] = numpy.zeros((geno_len[name],INT_NB_SIZE))
		fd_int[name] = numpy.zeros((geno_len[name],INT_NB_SIZE))
		fi_int[name] = numpy.zeros((geno_len[name],INT_NB_SIZE))
		start,stop = range_of_alns[name]
		x = PEhmm2.forward_sum_interior(fm_int_start[name][start:stop],fi_int_start[name][start:stop],hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],INT_NB_SIZE,INT_NB_P,fm_int[name][start:stop,:],fd_int[name][start:stop,:],fi_int[name][start:stop,:])
		if x:
			print 'error: forward sum interior:', x
			exit()
		
#		print 'int',name,start,stop, numpy.sum(fm_int_start[name][start:stop]), numpy.sum(fi_int_start[name][start:stop]), numpy.sum(fm_int[name][start:stop,-1]), numpy.sum(fd_int[name][start:stop,-1]), numpy.sum(fi_int[name][start:stop,-1])
	
#	print fm_int['POGO_20'][1459,-1]
	
	fm_2 = list()
	fd_2 = list()
	fi_2 = list()
	for alnrange in ranges2:
		name,start,stop = alnrange
		fm = numpy.zeros((stop-start,len(seq2)))
		fd = numpy.zeros((stop-start,len(seq2)))
		fi = numpy.zeros((stop-start,len(seq2)))
		fm2_init = INT_NB_P*( (1.0-hmm[name].del_open[start:stop]-hmm[name].ins_open[start-1:stop-1])*fm_int[name][start-1:stop-1,-1] + (1.0-hmm[name].del_extend[start:stop])*fd_int[name][start-1:stop-1,-1] + (1.0-hmm[name].ins_extend[start-1:stop-1])*fi_int[name][start-1:stop-1,-1])
		fi2_init = INT_NB_P*hmm[name].ins_open[start:stop]*fm_int[name][start:stop,-1] + hmm[name].ins_extend[start:stop]*fi_int[name][start:stop,-1]
#		fm,fd,fi = forward_sum(seq2,hmm[name].genome_profile[start:stop],hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],fm2_init,fi2_init)
		x = PEhmm2.forward_sum(seq2,hmm[name].genome_profile[start:stop],hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],fm2_init,fi2_init,fm,fd,fi)
		if x:
			print 'error: forward sum seq2:', x
			print name,start,stop
			pickle.dump([seq2,hmm[name].genome_profile[start:stop],hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],f_init[name][start:stop],fm,fd,fi],open('error.pkl','w'))
			exit()
		fm_2.append(fm)
		fd_2.append(fd)
		fi_2.append(fi)
	return fm_1,fd_1,fi_1,fm_int,fd_int,fi_int,range_of_alns,fm_2,fd_2,fi_2

def backward_sum_pair(seq1,ranges1,seq2,ranges2,hmm):
	geno_len = dict()
	bm_int_end = dict()
	bi_int_end = dict()
	for name in hmm:
		geno_len[name] = len(hmm[name].genome_profile)
		bm_int_end[name] = numpy.zeros(geno_len[name])#+MIN_READ_PROB
		bi_int_end[name] = numpy.zeros(geno_len[name])#+MIN_READ_PROB
	
	range_of_alns = dict()
	for alnrange in ranges1+ranges2:
		name,start,stop = alnrange
		if name not in range_of_alns:
			range_of_alns[name] = [start-1,stop+1]
		else:
			range_of_alns[name][0] = min(start-1,range_of_alns[name][0])
			range_of_alns[name][1] = max(stop+1,range_of_alns[name][1])
	bm_2 = list()
	bd_2 = list()
	bi_2 = list()
	for alnrange in ranges2:
		name,start,stop = alnrange
		bm = numpy.zeros((stop-start,len(seq2)))
		bd = numpy.zeros((stop-start,len(seq2)))
		bi = numpy.zeros((stop-start,len(seq2)))
#		bm,bd,bi = backward_sum(seq2,hmm[name].genome_profile[start:stop],hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],numpy.ones(stop-start),numpy.ones(stop-start))
		x = PEhmm2.backward_sum(seq2,hmm[name].genome_profile[start:stop],hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],numpy.ones(stop-start),numpy.ones(stop-start),bm,bd,bi)
		if x:
			print 'error: backward sum seq2:', x
			print name,start,stop
			pickle.dump([seq2,hmm[name].genome_profile[start:stop],hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],numpy.ones(stop-start),bm,bd,bi],open('error.pkl','w'))
			exit()
		bm_2.append(bm)
		bd_2.append(bd)
		bi_2.append(bi)
		bm_int_end[name][start:stop] += bm[:,0]*hmm[name].genome_profile[start:stop,seq2[0]]
		bi_int_end[name][start:stop] += bi[:,0]*0.25
	
	bm_int = dict()
	bd_int = dict()
	bi_int = dict()
	for name in range_of_alns:
		bm_int[name] = numpy.zeros((geno_len[name],INT_NB_SIZE))
		bd_int[name] = numpy.zeros((geno_len[name],INT_NB_SIZE))
		bi_int[name] = numpy.zeros((geno_len[name],INT_NB_SIZE))
		start,stop = range_of_alns[name]
		x = PEhmm2.backward_sum_interior(bm_int_end[name][start:stop],bi_int_end[name][start:stop],hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],INT_NB_SIZE,INT_NB_P,bm_int[name][start:stop,:],bd_int[name][start:stop,:], bi_int[name][start:stop,:])
		if x:
			print 'error: backward sum interior:', x
			exit()
# 		if name == 'POGO_13':
# 			print
# 			print bm_int_end[name][1855]
# 			print bi_int_end[name][1855]
# 			print bm_int[name][1387,0]
# 			print bi_int[name][1387,0]
	
#	print 'b_init',max(b_init),MIN_READ_PROB
	
	bm_1 = list()
	bd_1 = list()
	bi_1 = list()
	for alnrange in ranges1:
		name,start,stop = alnrange
#		print b_init[start:stop]
		bm = numpy.zeros((stop-start,len(seq1)))
		bd = numpy.zeros((stop-start,len(seq1)))
		bi = numpy.zeros((stop-start,len(seq1)))
		bm2_init = bm_int[name][start+1:stop+1,0]*(1.0-hmm[name].del_open[start+1:stop+1]-hmm[name].ins_open[start:stop]) + bd_int[name][start+1:stop+1,0]*hmm[name].del_open[start+1:stop+1] + bi_int[name][start:stop,0]*hmm[name].ins_open[start:stop]
		bi2_init = bm_int[name][start+1:stop+1,0]*(1.0-hmm[name].ins_extend[start:stop]) + bi_int[name][start:stop,0]*hmm[name].ins_extend[start:stop]
#		bm,bd,bi = backward_sum(seq1,hmm[name].genome_profile[start:stop],hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],bm2_init,bi2_init)
		x = PEhmm2.backward_sum(seq1,hmm[name].genome_profile[start:stop],hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],bm2_init,bi2_init,bm,bd,bi)
		if x:
			print 'error: backward sum seq2:', x
			print name,start,stop
			pickle.dump([seq1,hmm[name].genome_profile[start:stop],hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],b_init[name][start:stop],bm,bd,bi],open('error.pkl','w'))
			exit()
		bm_1.append(bm)
		bd_1.append(bd)
		bi_1.append(bi)
	return bm_1,bd_1,bi_1,bm_int,bd_int,bi_int,range_of_alns,bm_2,bd_2,bi_2

def find_seeded_alignments(sequences,del_sequences,read_kmer_dict,reads):
	
	for read in reads:
		read.alnranges_1 = list()
		read.alnranges_2 = list()
		read.alnranges_rc1 = list()
		read.alnranges_rc2 = list()
	
	for name in sequences:
		sequence = sequences[name]
		# All hits
		# Keys are read ids. Values are: [[read1,read1rc],[read2,read2rc]] where each element is a list of [extra_left,start,stop,extra_right] for potential alignments
		read_alignments = dict()
		seq_len = len(sequence)
		for i in range(seq_len-KMER_LEN):
			kmer = sequence[i:i+KMER_LEN]
			if kmer in read_kmer_dict:
				for seed_hit in read_kmer_dict[kmer]:
					read_id,read_pos,is_read2,is_reverse = seed_hit
					if not read_id in read_alignments:
						read_alignments[read_id] = [[[],[]],[[],[]]]
					added = False
					for alignment in read_alignments[read_id][is_read2][is_reverse][::-1]:
						extra_left,start,stop,extra_right = alignment
						if i+KMER_LEN > stop and READ_LEN - (read_pos+KMER_LEN) < extra_right:
							alignment[2] = i+KMER_LEN
							alignment[3] = READ_LEN - (read_pos+KMER_LEN)
							added = True
							break
					if not added:
						read_alignments[read_id][is_read2][is_reverse].append([read_pos,i,i+KMER_LEN,READ_LEN - (read_pos+KMER_LEN)])
		
		for deletion in del_sequences[name]:
			del_start,del_stop,del_seq = deletion
			for i in range(len(del_seq)-KMER_LEN):
				kmer = del_seq[i:i+KMER_LEN]
				if kmer in read_kmer_dict:
					for seed_hit in read_kmer_dict[kmer]:
						read_id,read_pos,is_read2,is_reverse = seed_hit
						if not read_id in read_alignments:
							read_alignments[read_id] = [[[],[]],[[],[]]]
						for i in range(len(read_alignments[read_id][is_read2][is_reverse])):
							extra_left,start,stop,extra_right = read_alignments[read_id][is_read2][is_reverse][i]
							seed_start = del_start-KMER_LEN+i
							seed_stop = del_stop+i
							if start < seed_start and seed_stop > stop:
								stop = seed_stop
								extra_right = READ_LEN - (read_pos+KMER_LEN)
								read_alignments[read_id][is_read2][is_reverse][i] = [extra_left,start,stop,extra_right]
							if seed_start < start and stop > seed_stop and read_pos<extra_left:
								start = seed_start
								extra_left = read_pos
								read_alignments[read_id][is_read2][is_reverse][i] = [extra_left,start,stop,extra_right]
							
		# Save hits
		for read in reads:
			
			if read.id not in read_alignments:
				continue
			if len(read_alignments[read.id][1][1]) > 0:
				for alignment in read_alignments[read.id][0][0]:
					extra_left,start,stop,extra_right = alignment
					read.alnranges_1.append((name,max(start-OVERHANG_FACTOR*extra_left-OVERHANG,1),min(stop+OVERHANG_FACTOR*extra_right+OVERHANG,seq_len-1)))
#					print read.id,'alnranges_1',name,max(start-OVERHANG_FACTOR*extra_left-OVERHANG,1),min(stop+OVERHANG_FACTOR*extra_right+OVERHANG,seq_len-1)
			
			if len(read_alignments[read.id][0][1]) > 0:
				for alignment in read_alignments[read.id][1][0]:
					extra_left,start,stop,extra_right = alignment
					read.alnranges_2.append((name,max(start-OVERHANG_FACTOR*extra_left-OVERHANG,1),min(stop+OVERHANG_FACTOR*extra_right+OVERHANG,seq_len-1)))
#					print read.id,'alnranges_2',name,max(start-OVERHANG_FACTOR*extra_left-OVERHANG,1),min(stop+OVERHANG_FACTOR*extra_right+OVERHANG,seq_len-1)
			
			if len(read_alignments[read.id][1][0]) > 0:
				for alignment in read_alignments[read.id][0][1]:
					extra_left,start,stop,extra_right = alignment
					read.alnranges_rc1.append((name,max(start-OVERHANG_FACTOR*extra_left-OVERHANG,1),min(stop+OVERHANG_FACTOR*extra_right+OVERHANG,seq_len-1)))
#					print read.id,'alnranges_rc1',name,max(start-OVERHANG_FACTOR*extra_left-OVERHANG,1),min(stop+OVERHANG_FACTOR*extra_right+OVERHANG,seq_len-1)
			
			if len(read_alignments[read.id][0][0]) > 0:
				for alignment in read_alignments[read.id][1][1]:
					extra_left,start,stop,extra_right = alignment
					read.alnranges_rc2.append((name,max(start-OVERHANG_FACTOR*extra_left-OVERHANG,1),min(stop+OVERHANG_FACTOR*extra_right+OVERHANG,seq_len-1)))
#					print read.id,'alnranges_rc2',name,max(start-OVERHANG_FACTOR*extra_left-OVERHANG,1),min(stop+OVERHANG_FACTOR*extra_right+OVERHANG,seq_len-1)
			
	return None

def save_fixed(input):
#		print 'starting save fixed'
	reads,hmm = input
	probs = list()
	for read in reads:
		
#		print 'start',read.id
		
		fm_f1,fd_f1,fi_f1,fm_int1,fd_int1,fi_int1,int_ranges1,fm_r2,fd_r2,fi_r2 = forward_sum_pair(read.seq1,read.alnranges_1,read.seqrc2,read.alnranges_rc2,hmm)
		bm_f1,bd_f1,bi_f1,bm_int1,bd_int1,bi_int1,int_ranges1,bm_r2,bd_r2,bi_r2 = backward_sum_pair(read.seq1,read.alnranges_1,read.seqrc2,read.alnranges_rc2,hmm)
		fm_f2,fd_f2,fi_f2,fm_int2,fd_int2,fi_int2,int_ranges2,fm_r1,fd_r1,fi_r1 = forward_sum_pair(read.seq2,read.alnranges_2,read.seqrc1,read.alnranges_rc1,hmm)
		bm_f2,bd_f2,bi_f2,bm_int2,bd_int2,bi_int2,int_ranges2,bm_r1,bd_r1,bi_r1 = backward_sum_pair(read.seq2,read.alnranges_2,read.seqrc1,read.alnranges_rc1,hmm)
		
		
		pr = MIN_READ_PROB*MIN_READ_PROB*MAX_INSERT_PDF
		
		for fm in fm_r2:
			pr += numpy.sum(fm[:,-1])
		for fi in fi_r2:
			pr += numpy.sum(fi[:,-1])
		for fm in fm_r1:
			pr += numpy.sum(fm[:,-1])
		for fi in fi_r1:
			pr += numpy.sum(fi[:,-1])
			
# 		pr_b = read.prob_from_fixed
# 		for i in range(len(read.alnranges_1)):
# 			name, start, stop = read.alnranges_1[i]
# 			pr_b += numpy.sum(bm_f1[i][:,0]*hmm[name].genome_profile[start:stop,read.seq1[0]]*(1.0-hmm[name].ins_open[start:stop]))
# 			pr_b += numpy.sum(bi_f1[i][:,0]*0.25*hmm[name].ins_open[start:stop])
# 		for i in range(len(read.alnranges_2)):
# 			name, start, stop = read.alnranges_2[i]
# 			pr_b += numpy.sum(bm_f2[i][:,0]*hmm[name].genome_profile[start:stop,read.seq2[0]]*(1.0-hmm[name].ins_open[start:stop]))
# 			pr_b += numpy.sum(bi_f2[i][:,0]*0.25*hmm[name].ins_open[start:stop])
# 		
# 		if max(pr,pr_b) / min(pr,pr_b) > 1.1:# and max(pr,pr_b) > 10**-10:
# 			print 'Forward/reverse sum mismatch. Using larger value.',read.id, max(pr,pr_b)/min(pr,pr_b),pr,pr_b
# 		pr = max(pr,pr_b)
		
		probs.append((read.id,pr))
		
#		print 'stop', read.id
		
#	print 'ending save fixed'
	return probs

def Estep(input):
	reads,hmm = input
	expU = dict()
	match_count = dict()
	del_count = dict()
	ins_count = dict()
	del_open_count = dict()
	del_extend_count = dict()
	ins_open_count = dict()
	ins_extend_count = dict()
	for name in hmm:
		expU[name] = numpy.zeros((len(hmm[name].genome_profile),4))
		match_count[name] = numpy.zeros(len(hmm[name].genome_profile))
		del_count[name] = numpy.zeros(len(hmm[name].genome_profile))
		ins_count[name] = numpy.zeros(len(hmm[name].genome_profile))
		del_open_count[name] = numpy.zeros(len(hmm[name].genome_profile))
		del_extend_count[name] = numpy.zeros(len(hmm[name].genome_profile))
		ins_open_count[name] = numpy.zeros(len(hmm[name].genome_profile))
		ins_extend_count[name] = numpy.zeros(len(hmm[name].genome_profile))
	
	loglik = 0
	unaligned = 0
	
	for read in reads:
		if len(reads_to_check) > 0 and read.id not in reads_to_check:
			continue
			
		if (len(read.alnranges_1)==0 or len(read.alnranges_rc2)==0) and (len(read.alnranges_2)==0 or len(read.alnranges_rc1)==0):
			unaligned += 1
				
		fm_f1,fd_f1,fi_f1,fm_int1,fd_int1,fi_int1,int_ranges1,fm_r2,fd_r2,fi_r2 = forward_sum_pair(read.seq1,read.alnranges_1,read.seqrc2,read.alnranges_rc2,hmm)
		bm_f1,bd_f1,bi_f1,bm_int1,bd_int1,bi_int1,int_ranges1,bm_r2,bd_r2,bi_r2 = backward_sum_pair(read.seq1,read.alnranges_1,read.seqrc2,read.alnranges_rc2,hmm)
		fm_f2,fd_f2,fi_f2,fm_int2,fd_int2,fi_int2,int_ranges2,fm_r1,fd_r1,fi_r1 = forward_sum_pair(read.seq2,read.alnranges_2,read.seqrc1,read.alnranges_rc1,hmm)
		bm_f2,bd_f2,bi_f2,bm_int2,bd_int2,bi_int2,int_ranges2,bm_r1,bd_r1,bi_r1 = backward_sum_pair(read.seq2,read.alnranges_2,read.seqrc1,read.alnranges_rc1,hmm)
		
		pr = read.prob_from_fixed
				
		for fm in fm_r2:
			pr += numpy.sum(fm[:,-1])
		for fi in fi_r2:
			pr += numpy.sum(fi[:,-1])
		for fm in fm_r1:
			pr += numpy.sum(fm[:,-1])
		for fi in fi_r1:
			pr += numpy.sum(fi[:,-1])
			
# 		pr_b = read.prob_from_fixed
# 		for i in range(len(read.alnranges_1)):
# 			name, start, stop = read.alnranges_1[i]
# 			pr_b += numpy.sum(bm_f1[i][:,0]*hmm[name].genome_profile[start:stop,read.seq1[0]]*(1.0-hmm[name].ins_open[start:stop]))
# 			pr_b += numpy.sum(bi_f1[i][:,0]*0.25*hmm[name].ins_open[start:stop])
# 		for i in range(len(read.alnranges_2)):
# 			name, start, stop = read.alnranges_2[i]
# 			pr_b += numpy.sum(bm_f2[i][:,0]*hmm[name].genome_profile[start:stop,read.seq2[0]]*(1.0-hmm[name].ins_open[start:stop]))
# 			pr_b += numpy.sum(bi_f2[i][:,0]*0.25*hmm[name].ins_open[start:stop])
		
# 		pr_dict = dict()
# 		prb_dict = dict()
# 		for i in range(len(read.alnranges_1)):
# 			name, start, stop = read.alnranges_1[i]
# 			if name not in pr_dict:
# 				pr_dict[name] = 0.0
# 				prb_dict[name] = 0.0
# 			prb_dict[name] += numpy.sum(bm_f1[i][:,0]*hmm[name].genome_profile[start:stop,read.seq1[0]]*(1.0-hmm[name].ins_open[start:stop]))
# 			prb_dict[name] += numpy.sum(bi_f1[i][:,0]*0.25*hmm[name].ins_open[start:stop])
# 		for i in range(len(read.alnranges_2)):
# 			name, start, stop = read.alnranges_2[i]
# 			if name not in pr_dict:
# 				pr_dict[name] = 0.0
# 				prb_dict[name] = 0.0
# 			prb_dict[name] += numpy.sum(bm_f2[i][:,0]*hmm[name].genome_profile[start:stop,read.seq2[0]]*(1.0-hmm[name].ins_open[start:stop]))
# 			prb_dict[name] += numpy.sum(bi_f2[i][:,0]*0.25*hmm[name].ins_open[start:stop])
# 		for i in range(len(read.alnranges_rc1)):
# 			name, start, stop = read.alnranges_rc1[i]
# 			if name not in pr_dict:
# 				pr_dict[name] = 0.0
# 				prb_dict[name] = 0.0
# 			pr_dict[name] += numpy.sum(fm_r1[i][:,-1])
# 			pr_dict[name] += numpy.sum(fi_r1[i][:,-1])
# 		for i in range(len(read.alnranges_rc2)):
# 			name, start, stop = read.alnranges_rc2[i]
# 			if name not in pr_dict:
# 				pr_dict[name] = 0.0
# 				prb_dict[name] = 0.0
# 			pr_dict[name] += numpy.sum(fm_r2[i][:,-1])
# 			pr_dict[name] += numpy.sum(fi_r2[i][:,-1])
# 		
# 		for name in pr_dict:
# 			if max(pr_dict[name],prb_dict[name]) / min(pr_dict[name],prb_dict[name]) > 2.0 and max(pr_dict[name],prb_dict[name]) > 10**-5:
# 				print name,read.id, max(pr_dict[name],prb_dict[name])/min(pr_dict[name],prb_dict[name]),pr_dict[name],prb_dict[name],pr,pr_b
# 				
# 				if len(hmm[name].ins_open) < 2215:
# 					continue
# 				
# 				fm_init = numpy.zeros(len(hmm[name].ins_open))
# 				fm_init[1308] = 1.0-hmm[name].ins_open[1308]
# 				fi_init = numpy.zeros(len(hmm[name].ins_open))
# 				fi_init[1308] = hmm[name].ins_open[1308]
# 				bm_init = numpy.zeros(len(hmm[name].ins_open))
# 				bm_init[2215] = 1.0-hmm[name].ins_open[2215]
# 				bi_init = numpy.zeros(len(hmm[name].ins_open))
# 				bi_init[2215] = hmm[name].ins_open[2215]
# 				fm_int = numpy.zeros((len(hmm[name].ins_open),INT_NB_SIZE))
# 				fd_int = numpy.zeros((len(hmm[name].ins_open),INT_NB_SIZE))
# 				fi_int = numpy.zeros((len(hmm[name].ins_open),INT_NB_SIZE))
# 				bm_int = numpy.zeros((len(hmm[name].ins_open),INT_NB_SIZE))
# 				bd_int = numpy.zeros((len(hmm[name].ins_open),INT_NB_SIZE))
# 				bi_int = numpy.zeros((len(hmm[name].ins_open),INT_NB_SIZE))
# 				x = PEhmm2.forward_sum_interior(fm_init,fi_init,hmm[name].del_open,hmm[name].del_extend,hmm[name].ins_open,hmm[name].ins_extend,INT_NB_SIZE,INT_NB_P,fm_int,fd_int,fi_int)
# 				if x:
# 					print 'error: forward sum interior:', x
# 					exit()
# 				x = PEhmm2.backward_sum_interior(bm_init,bi_init,hmm[name].del_open,hmm[name].del_extend,hmm[name].ins_open,hmm[name].ins_extend,INT_NB_SIZE,INT_NB_P,bm_int,bd_int,bi_int)
# 				if x:
# 					print 'error: backward sum interior:', x
# 					exit()
# 				print name,INT_NB_P*fm_int[2215,-1]+hmm[name].ins_extend[2215]*fi_int[2215,-1],bm_int[1308,0]*(1.0-hmm[name].ins_open[1308])+bi_int[1308,0]*hmm[name].ins_open[1308]
# 
# 				exit()
		
# 		if max(pr,pr_b) / min(pr,pr_b) > 1.1:# and max(pr,pr_b):
# 			print 'Forward/reverse sum mismatch. Using larger value.',read.id, max(pr,pr_b)/min(pr,pr_b),pr,pr_b
# 		pr = max(pr,pr_b)
		loglik += numpy.log(pr)
		
		for i in range(len(read.alnranges_1)):
			name,start,stop = read.alnranges_1[i]
			match_mat = fm_f1[i][:,:-1]*bm_f1[i][:,:-1]/pr
			match_count[name][start:stop] += numpy.sum(match_mat,1)
			for j in range(len(read.seq1)-1):
				expU[name][start:stop,read.seq1[j]] += match_mat[:,j]
			del_count[name][start:stop] += numpy.sum(fd_f1[i][:,:-1]*bd_f1[i][:,:-1],1)/pr
			ins_count[name][start:stop] += numpy.sum(fi_f1[i][:,:-1]*bi_f1[i][:,:-1],1)/pr
			del_open_count[name][start+1:stop] += numpy.sum(fm_f1[i][:-1,:]*numpy.transpose(numpy.tile(hmm[name].del_open[start+1:stop],(len(read.seq1),1)))*bd_f1[i][1:,:],1)/pr
			del_extend_count[name][start+1:stop] += numpy.sum(fd_f1[i][:-1,:]*numpy.transpose(numpy.tile(hmm[name].del_extend[start+1:stop],(len(read.seq1),1)))*bd_f1[i][1:,:],1)/pr
			ins_open_count[name][start:stop] += numpy.sum(fm_f1[i][:,:-1]*numpy.transpose(numpy.tile(hmm[name].ins_open[start:stop],(len(read.seq1)-1,1)))*0.25*bi_f1[i][:,1:],1)/pr
			ins_extend_count[name][start:stop] += numpy.sum(fi_f1[i][:,:-1]*numpy.transpose(numpy.tile(hmm[name].ins_extend[start:stop],(len(read.seq1)-1,1)))*0.25*bi_f1[i][:,1:],1)/pr
			
		for i in range(len(read.alnranges_rc1)):
			name,start,stop = read.alnranges_rc1[i]
			match_mat = fm_r1[i][:,:-1]*bm_r1[i][:,:-1]/pr
			match_count[name][start:stop] += numpy.sum(match_mat,1)
			for j in range(len(read.seqrc1)-1):
				expU[name][start:stop,read.seqrc1[j]] += match_mat[:,j]
			del_count[name][start:stop] += numpy.sum(fd_r1[i][:,:-1]*bd_r1[i][:,:-1],1)/pr
			ins_count[name][start:stop] += numpy.sum(fi_r1[i][:,:-1]*bi_r1[i][:,:-1],1)/pr
			del_open_count[name][start+1:stop] += numpy.sum(fm_r1[i][:-1,:]*numpy.transpose(numpy.tile(hmm[name].del_open[start+1:stop],(len(read.seqrc1),1)))*bd_r1[i][1:,:],1)/pr
			del_extend_count[name][start+1:stop] += numpy.sum(fd_r1[i][:-1,:]*numpy.transpose(numpy.tile(hmm[name].del_extend[start+1:stop],(len(read.seqrc1),1)))*bd_r1[i][1:,:],1)/pr
			ins_open_count[name][start:stop] += numpy.sum(fm_r1[i][:,:-1]*numpy.transpose(numpy.tile(hmm[name].ins_open[start:stop],(len(read.seqrc1)-1,1)))*0.25*bi_r1[i][:,1:],1)/pr
			ins_extend_count[name][start:stop] += numpy.sum(fi_r1[i][:,:-1]*numpy.transpose(numpy.tile(hmm[name].ins_extend[start:stop],(len(read.seqrc1)-1,1)))*0.25*bi_r1[i][:,1:],1)/pr
			
		for i in range(len(read.alnranges_2)):
			name,start,stop = read.alnranges_2[i]
			match_mat = fm_f2[i][:,:-1]*bm_f2[i][:,:-1]/pr
			match_count[name][start:stop] += numpy.sum(match_mat,1)
			for j in range(len(read.seq2)-1):
				expU[name][start:stop,read.seq2[j]] += match_mat[:,j]
			del_count[name][start:stop] += numpy.sum(fd_f2[i][:,:-1]*bd_f2[i][:,:-1],1)/pr
			ins_count[name][start:stop] += numpy.sum(fi_f2[i][:,:-1]*bi_f2[i][:,:-1],1)/pr
			del_open_count[name][start+1:stop] += numpy.sum(fm_f2[i][:-1,:]*numpy.transpose(numpy.tile(hmm[name].del_open[start+1:stop],(len(read.seq2),1)))*bd_f2[i][1:,:],1)/pr
			del_extend_count[name][start+1:stop] += numpy.sum(fd_f2[i][:-1,:]*numpy.transpose(numpy.tile(hmm[name].del_extend[start+1:stop],(len(read.seq2),1)))*bd_f2[i][1:,:],1)/pr
			ins_open_count[name][start:stop] += numpy.sum(fm_f2[i][:,:-1]*numpy.transpose(numpy.tile(hmm[name].ins_open[start:stop],(len(read.seq2)-1,1)))*0.25*bi_f2[i][:,1:],1)/pr
			ins_extend_count[name][start:stop] += numpy.sum(fi_f2[i][:,:-1]*numpy.transpose(numpy.tile(hmm[name].ins_extend[start:stop],(len(read.seq2)-1,1)))*0.25*bi_f2[i][:,1:],1)/pr
			
		for i in range(len(read.alnranges_rc2)):
			name,start,stop = read.alnranges_rc2[i]
			match_mat = fm_r2[i][:,:-1]*bm_r2[i][:,:-1]/pr
			match_count[name][start:stop] += numpy.sum(match_mat,1)
			for j in range(len(read.seqrc2)-1):
				expU[name][start:stop,read.seqrc2[j]] += match_mat[:,j]
			del_count[name][start:stop] += numpy.sum(fd_r2[i][:,:-1]*bd_r2[i][:,:-1],1)/pr
			ins_count[name][start:stop] += numpy.sum(fi_r2[i][:,:-1]*bi_r2[i][:,:-1],1)/pr
			del_open_count[name][start+1:stop] += numpy.sum(fm_r2[i][:-1,:]*numpy.transpose(numpy.tile(hmm[name].del_open[start+1:stop],(len(read.seqrc2),1)))*bd_r2[i][1:,:],1)/pr
			del_extend_count[name][start+1:stop] += numpy.sum(fd_r2[i][:-1,:]*numpy.transpose(numpy.tile(hmm[name].del_extend[start+1:stop],(len(read.seqrc2),1)))*bd_r2[i][1:,:],1)/pr
			ins_open_count[name][start:stop] += numpy.sum(fm_r2[i][:,:-1]*numpy.transpose(numpy.tile(hmm[name].ins_open[start:stop],(len(read.seqrc2)-1,1)))*0.25*bi_r2[i][:,1:],1)/pr
			ins_extend_count[name][start:stop] += numpy.sum(fi_r2[i][:,:-1]*numpy.transpose(numpy.tile(hmm[name].ins_extend[start:stop],(len(read.seqrc2)-1,1)))*0.25*bi_r2[i][:,1:],1)/pr
		
		for name in fm_int1:
			start,stop = int_ranges1[name]
			x = PEhmm2.add_counts_interior(fm_int1[name][start:stop],fd_int1[name][start:stop],fi_int1[name][start:stop],bm_int1[name][start:stop],bd_int1[name][start:stop],bi_int1[name][start:stop],
			                          hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],
			                          INT_NB_SIZE,INT_NB_P,pr,
			                          match_count[name][start:stop],del_count[name][start:stop],ins_count[name][start:stop],del_open_count[name][start:stop],del_extend_count[name][start:stop],ins_open_count[name][start:stop],ins_extend_count[name][start:stop])
			if x:
				print 'adding int1 error:', x
				exit()
				
			
		for name in fm_int2:
			start,stop = int_ranges2[name]
			x = PEhmm2.add_counts_interior(fm_int2[name][start:stop],fd_int2[name][start:stop],fi_int2[name][start:stop],bm_int2[name][start:stop],bd_int2[name][start:stop],bi_int2[name][start:stop],
			                          hmm[name].del_open[start:stop],hmm[name].del_extend[start:stop],hmm[name].ins_open[start:stop],hmm[name].ins_extend[start:stop],
			                          INT_NB_SIZE,INT_NB_P,pr,
			                          match_count[name][start:stop],del_count[name][start:stop],ins_count[name][start:stop],del_open_count[name][start:stop],del_extend_count[name][start:stop],ins_open_count[name][start:stop],ins_extend_count[name][start:stop])
			if x:
				print 'adding int2 error:', x
				exit()
	
	return 	loglik, unaligned, expU, match_count, del_count, ins_count, del_open_count, del_extend_count, ins_open_count, ins_extend_count
	
def EMstep(reads,hmm,sequences,master_pool):
	
	loglik = 0.0
	unaligned = 0
	expU = dict()
	match_count = dict()
	del_count = dict()
	ins_count = dict()
	del_open_count = dict()
	del_extend_count = dict()
	ins_open_count = dict()
	ins_extend_count = dict()
	for name in hmm:
		expU[name] = copy.copy(hmm[name].genome_profile)
		match_count[name] = numpy.zeros(len(hmm[name].genome_profile))+10**-5
		del_count[name] = numpy.zeros(len(hmm[name].genome_profile))+10**-5
		ins_count[name] = numpy.zeros(len(hmm[name].genome_profile))+10**-5
		del_open_count[name] = numpy.zeros(len(hmm[name].genome_profile))
		del_extend_count[name] = numpy.zeros(len(hmm[name].genome_profile))
		ins_open_count[name] = numpy.zeros(len(hmm[name].genome_profile))
		ins_extend_count[name] = numpy.zeros(len(hmm[name].genome_profile))
	
	# Divide up reads and send to separate processes for E step
	inputs = list()
	for i in range(THREADS):
		inputs.append([reads[i*len(reads)/THREADS:(i+1)*len(reads)/THREADS],hmm])
#	Estep_outputs = list(map(Estep, inputs))
	Estep_outputs = master_pool.map(Estep, inputs)
	for Estep_output in Estep_outputs:
		this_loglik, this_unaligned, this_expU, this_match_count, this_del_count, this_ins_count, this_del_open_count, this_del_extend_count, this_ins_open_count, this_ins_extend_count = Estep_output
		loglik += this_loglik
		unaligned += this_unaligned
		for name in hmm:
			expU[name] += this_expU[name]
			match_count[name] += this_match_count[name]
			del_count[name] += this_del_count[name]
			ins_count[name] += this_ins_count[name]
			del_open_count[name] += this_del_open_count[name]
			del_extend_count[name] += this_del_extend_count[name]
			ins_open_count[name] += this_ins_open_count[name]
			ins_extend_count[name] += this_ins_extend_count[name]
	
#	print ','.join([str(x) for x in del_count['POGO_21'][989:995]])
#	print ','.join([str(x) for x in del_open_count['POGO_21'][989:995]])
#	print ','.join([str(x) for x in del_extend_count['POGO_21'][989:995]])
	
	for name in hmm:
		for i in range(1000,len(expU[name])-1000):
			hmm[name].genome_profile[i,:] = numpy.maximum(expU[name][i,:]/numpy.sum(expU[name][i,:]),MINP*(1+3*MINP))
			hmm[name].genome_profile[i,:] = hmm[name].genome_profile[i,:]/numpy.sum(hmm[name].genome_profile[i,:])
		hmm[name].ins_open[1000:-1000] = numpy.minimum(numpy.maximum(ins_open_count[name][1000:-1000]/(match_count[name][1000:-1000]),MINP),1-MINP)
		hmm[name].del_open[1000:-1000] = numpy.minimum(numpy.maximum(del_open_count[name][1000:-1000]/(match_count[name][999:-1001]),MINP),1-MINP)
		hmm[name].ins_extend[1000:-1000] = numpy.minimum(numpy.maximum(ins_extend_count[name][1000:-1000]/ins_count[name][1000:-1000],MINP),1-MINP)
		hmm[name].del_extend[1000:-1000] = numpy.minimum(numpy.maximum(del_extend_count[name][1000:-1000]/del_count[name][999:-1001],MINP),1-MINP)
		hmm[name].match_count = match_count[name]
		hmm[name].del_count = del_count[name]
		hmm[name].del_open_count = del_open_count[name]
		hmm[name].del_extend_count = del_extend_count[name]
		hmm[name].ins_count = ins_count[name]
		hmm[name].ins_open_count = ins_open_count[name]
		hmm[name].ins_extend_count = ins_extend_count[name]
		
		hmm[name].weight[:] = 1.0
# 		for i in range(950,len(hmm[name].weight)-950):
# 			hmm[name].weight[i] = hmm[name].weight[i-1]*(1.0-hmm[name].del_open[i])+(1.0-hmm[name].weight[i-1])*(1.0-hmm[name].del_extend[i])
# 		
# 	totweight = 0.0
# 	totlen = 0
# 	for name in hmm:
# 		totweight += numpy.sum(hmm[name].weight)
# 		totlen += len(hmm[name].weight)
# 		
# 	for name in hmm:
# 		hmm[name].weight = hmm[name].weight/totweight*totlen
	
	return loglik,hmm,unaligned

def nEMsteps(N,hmm,sequences,read_kmer_dict,reads):
	master_pool = Pool(THREADS)
	pickle.dump(hmm,open('hmm.0.pkl','w'))
	
	hmm_fixed = dict()
	hmm_to_update = dict()
	sequences_fixed = dict()
	sequences_to_update = dict()
	del_seqeunces = dict()
	fixed_len = 0
	for name in hmm:
		if name in TO_UPDATE:
			hmm_to_update[name] = hmm[name]
			sequences_to_update[name] = sequences[name]
		else:
			hmm_fixed = hmm[name]
			sequences_fixed[name] = sequences[name]
			fixed_len += len(sequences[name])-2000
			del_seqeunces[name] = list()
	
# 	hmm_to_update = pickle.load(open('hmm.start.pkl'))
# 	for name in hmm_to_update:
# 		sequences_to_update[name] = ''
# 		for i in range(len(hmm_to_update[name].genome_profile)):
# 			sequences_to_update[name] += 'ACGT'[numpy.argmax(hmm_to_update[name].genome_profile[i,:])]
	
# 	for name in hmm_to_update:
# 		fm_init = numpy.zeros(len(hmm_to_update[name].ins_open))
# 		fm_init[0] = 1.0
# 		fi_init = numpy.zeros(len(hmm_to_update[name].ins_open))
# 		fm_int = numpy.zeros((len(hmm_to_update[name].ins_open),INT_NB_SIZE))
# 		fd_int = numpy.zeros((len(hmm_to_update[name].ins_open),INT_NB_SIZE))
# 		fi_int = numpy.zeros((len(hmm_to_update[name].ins_open),INT_NB_SIZE))
# 		bm_int = numpy.zeros((len(hmm_to_update[name].ins_open),INT_NB_SIZE))
# 		bd_int = numpy.zeros((len(hmm_to_update[name].ins_open),INT_NB_SIZE))
# 		bi_int = numpy.zeros((len(hmm_to_update[name].ins_open),INT_NB_SIZE))
# 		x = PEhmm2.forward_sum_interior(fm_init,fi_init,hmm_to_update[name].del_open,hmm_to_update[name].del_extend,hmm_to_update[name].ins_open,hmm_to_update[name].ins_extend,INT_NB_SIZE,INT_NB_P,fm_int,fd_int,fi_int)
# 		if x:
# 			print 'error: forward sum interior:', x
# 			exit()
# 		x = PEhmm2.backward_sum_interior(numpy.ones(len(hmm_to_update[name].genome_profile)),numpy.ones(len(hmm_to_update[name].genome_profile)),hmm_to_update[name].del_open,hmm_to_update[name].del_extend,hmm_to_update[name].ins_open,hmm_to_update[name].ins_extend,INT_NB_SIZE,INT_NB_P,bm_int,bd_int,bi_int)
# 		if x:
# 			print 'error: backward sum interior:', x
# 			exit()
# 		print name,numpy.sum(INT_NB_P*(1.0-hmm_to_update[name].del_open-hmm_to_update[name].ins_open)*fm_int[:,-1])+numpy.sum(hmm_to_update[name].ins_extend*fi_int[:,-1]),bm_int[0,0]*(1.0-hmm_to_update[name].ins_open[0])+bi_int[0,0]*hmm_to_update[name].ins_open[0]
#	exit()
		
	print 'calculating fixed read probs...'
	fixed_starttime = datetime.datetime.now()
	
	find_seeded_alignments(sequences_fixed,del_seqeunces,read_kmer_dict,reads)
	inputs = list()
	for i in range(THREADS):
		inputs.append([reads[i*len(reads)/THREADS:(i+1)*len(reads)/THREADS],hmm])
#	save_outputs = list(map(save_fixed, inputs))
	save_outputs = master_pool.map(save_fixed, inputs)
	
	probs_to_save = list()
	for output in save_outputs:
		probs_to_save += output
	for i in range(len(reads)):
		if probs_to_save[i][0] != reads[i].id:
			print 'wrong read', probs_to_save[i][0], reads[i].id
			exit()
		reads[i].prob_from_fixed = probs_to_save[i][1]
	
	print datetime.datetime.now()-fixed_starttime
	
	print 'running EM steps...'
	last_loglik = -numpy.inf
	last_hmm_to_update = None
	for k in range(N):
		starttime = datetime.datetime.now()
		for name in hmm_to_update:
			sequences_to_update[name] = ''
			for i in range(len(hmm_to_update[name].genome_profile)):
				sequences_to_update[name] += 'ACGT'[numpy.argmax(hmm_to_update[name].genome_profile[i,:])]
		del_sequences = dict()
		for name in hmm_to_update:
			del_sequences[name] = list()
			i = 1000
			while i < len(hmm_to_update[name].genome_profile)-1000:
				if 0.1 < hmm_to_update[name].del_open[i]:
					j=i+1
					while DEL_END < hmm_to_update[name].del_extend[j] and j<=len(hmm_to_update[name].del_open)-1000:
						j += 1
					del_sequences[name].append((i,j,sequences_to_update[name][i-KMER_LEN:i] + sequences_to_update[name][j:j+KMER_LEN]))
				i += 1
#		print del_sequences['chrX:15421973:chrUn_DS483757v1:4833-6900']
		find_seeded_alignments(sequences_to_update,del_sequences,read_kmer_dict,reads)
		loglik,hmm,unaligned = EMstep(reads,hmm_to_update,sequences_to_update,master_pool)
		total_len = fixed_len
		for name in hmm_to_update:
			total_len += len(hmm[name].genome_profile)-2000
		pickle.dump(hmm_to_update,open('hmm.'+str(k+1)+'.pkl','w'))
		loglik = loglik-len(reads)*numpy.log(total_len)
		print 'step', k,loglik, datetime.datetime.now()-starttime,len(reads)
		
		if False and k > 99 and (k+1) % 10 == 0:
			if loglik < last_loglik:
				print 'reverting'
				hmm_to_update = last_hmm_to_update
				loglik = last_loglik
				for name in hmm_to_update:
					for pos in included_dels[name]:
						hmm_to_update[name].del_open[pos] = MINP
					for pos in included_ins[name]:
						hmm_to_update[name].ins_open[pos] = MINP
			else:
				last_hmm_to_update = copy.deepcopy(hmm_to_update)
			last_loglik = loglik
			included_dels = dict()
			included_ins = dict()
			tot_indel_open = MINP
			for name in hmm_to_update:
				tot_indel_open += numpy.sum(hmm_to_update[name].ins_open-MIN_TO_INCLUDE)
				tot_indel_open += numpy.sum(hmm_to_update[name].del_open-MIN_TO_INCLUDE)
			if (k+1) % 20 != 0:
				tot_indel_open = max(tot_indel_open,1.0)
			if k == N-1:
				continue
			for name in hmm_to_update:
				i = 1000
				included_dels[name] = list()
				included_ins[name] = list()
				offset = 0
				while i < len(hmm_to_update[name].genome_profile)-1000:
					if numpy.random.random() < (hmm_to_update[name].ins_open[i]-MINP)/(1-MINP)/tot_indel_open*1:
						for insert_len in range(1,51):
							if numpy.random.random() > hmm_to_update[name].ins_extend[i]:
								break
						hmm_to_update[name].genome_profile = numpy.concatenate(( numpy.concatenate(( hmm_to_update[name].genome_profile[:i+1,:] , (numpy.zeros((insert_len,4))+0.25) )), hmm_to_update[name].genome_profile[i+1:,:] ))
						hmm_to_update[name].ins_open = numpy.append( numpy.append( hmm_to_update[name].ins_open[:i+1] , (numpy.zeros(insert_len)+INDEL_OPEN) ), hmm_to_update[name].ins_open[i+1:] )
						hmm_to_update[name].ins_extend = numpy.append( numpy.append( hmm_to_update[name].ins_extend[:i+1] , (numpy.zeros(insert_len)+INS_EXTEND) ), hmm_to_update[name].ins_extend[i+1:] )
						hmm_to_update[name].del_open = numpy.append( numpy.append( hmm_to_update[name].del_open[:i+1] , (numpy.zeros(insert_len)+INDEL_OPEN) ), hmm_to_update[name].del_open[i+1:] )
						hmm_to_update[name].del_extend = numpy.append( numpy.append( hmm_to_update[name].del_extend[:i+1] , (numpy.zeros(insert_len)+DEL_EXTEND) ), hmm_to_update[name].del_extend[i+1:] )
						hmm_to_update[name].weight = numpy.append( numpy.append( hmm_to_update[name].weight[:i+1] , (numpy.zeros(insert_len)+DEL_EXTEND) ), hmm_to_update[name].weight[i+1:] )
						print 'insert',name,i+offset,insert_len
						included_ins[name].append(i+offset)
						offset -= insert_len
						i += insert_len+1
						continue
					if numpy.random.random() < (hmm_to_update[name].del_open[i]-MINP)/(1-MINP)/tot_indel_open*1:
						j=i+1
						while numpy.random.random() < hmm_to_update[name].del_extend[j] and j<=len(hmm_to_update[name].del_open)-1000:
							j += 1
						hmm_to_update[name].genome_profile = numpy.concatenate(( hmm_to_update[name].genome_profile[:i,:],hmm_to_update[name].genome_profile[j:,:] ))
						hmm_to_update[name].ins_open = numpy.append( hmm_to_update[name].ins_open[:i],hmm_to_update[name].ins_open[j:])
						hmm_to_update[name].ins_extend = numpy.append( hmm_to_update[name].ins_extend[:i],hmm_to_update[name].ins_extend[j:])
						hmm_to_update[name].del_open = numpy.append( hmm_to_update[name].del_open[:i],hmm_to_update[name].del_open[j:])
						hmm_to_update[name].del_extend = numpy.append( hmm_to_update[name].del_extend[:i],hmm_to_update[name].del_extend[j:])
						hmm_to_update[name].weight = numpy.append( hmm_to_update[name].weight[:i],hmm_to_update[name].weight[j:])
						print 'delete',name,i+offset,j+offset
						included_dels[name].append(i+offset)
						offset += j-i
						continue
					i += 1
		
	master_pool.close()
	master_pool.join()
	
	return loglik,sequences_to_update

def main():
	
	print 'reading reference...'
	reference = fasta_as_array(sys.argv[1])
	genome_profile = create_initial_guess(reference,numpy.array([1-3*P_ERROR,P_ERROR,P_ERROR,P_ERROR]))
#	reference = reference.values()[0]
#	genome_profile = genome_profile.values()[0]
	
	print 'reading reads...'
	fastaq1_dict = SeqIO.to_dict(SeqIO.parse(open(sys.argv[2],'r'), "fastq"))
	fastaq2_dict = SeqIO.to_dict(SeqIO.parse(open(sys.argv[3],'r'), "fastq"))
	reads = list()
	for name in fastaq1_dict:
		if 'N' not in str(fastaq1_dict[name].seq).upper() and 'N' not in str(fastaq2_dict[name].seq).upper():
			reads.append(PEreadclass(name,seq2array(str(fastaq1_dict[name].seq).upper()),seq2array(str(fastaq2_dict[name].seq).upper())))
	
	print 'building hmm...'
	hmm = dict()
	for name in reference:
		del_open = numpy.zeros(len(genome_profile[name]))+INDEL_OPEN
		del_open[1000]=0.1
		ins_open = numpy.zeros(len(genome_profile[name]))+INDEL_OPEN
		del_extend = numpy.zeros(len(genome_profile[name]))+DEL_EXTEND
		del_extend[-1000]=0.1
		ins_extend = numpy.zeros(len(genome_profile[name]))+INS_EXTEND
		hmm[name] = hmmclass(genome_profile[name],del_open,del_extend,ins_open,ins_extend)

# 	hmm['chrX:15421973:chrUn_DS483757v1:4833-6900'].del_open[1627] = 0.9
# 	hmm['chrX:15421973:chrUn_DS483757v1:4833-6900'].del_extend[1628:2294] = 0.999
# 	hmm['chrX:15421973:chrUn_DS483757v1:4833-6900'].del_extend[2295] = 0.001
# 	hmm['chrX:15421973:chrUn_DS483757v1:4833-6900'].ins_open[2956] = 0.999
	
	del genome_profile
	del del_open
	del del_extend
	del ins_open
	del ins_extend
	
	print 'building kmer dictionary...'
	kmerstarttime = datetime.datetime.now()
	# Keys are kmers. Values are tuples: (read_id,read_pos,is_read2,is_reverse).
	read_kmer_dict = dict()
	for read in reads:
		readseq = array2seq(read.seq1)
		for i in range(len(readseq)-KMER_LEN):
			if readseq[i:i+KMER_LEN] in read_kmer_dict:
				read_kmer_dict[readseq[i:i+KMER_LEN]] += [(read.id,i,False,False)]
			else:
				read_kmer_dict[readseq[i:i+KMER_LEN]] = [(read.id,i,False,False)]
		readseq = array2seq(read.seqrc1)
		for i in range(len(readseq)-KMER_LEN):
			if readseq[i:i+KMER_LEN] in read_kmer_dict:
				read_kmer_dict[readseq[i:i+KMER_LEN]] += [(read.id,i,False,True)]
			else:
				read_kmer_dict[readseq[i:i+KMER_LEN]] = [(read.id,i,False,True)]
		readseq = array2seq(read.seq2)
		for i in range(len(readseq)-KMER_LEN):
			if readseq[i:i+KMER_LEN] in read_kmer_dict:
				read_kmer_dict[readseq[i:i+KMER_LEN]] += [(read.id,i,True,False)]
			else:
				read_kmer_dict[readseq[i:i+KMER_LEN]] = [(read.id,i,True,False)]
		readseq = array2seq(read.seqrc2)
		for i in range(len(readseq)-KMER_LEN):
			if readseq[i:i+KMER_LEN] in read_kmer_dict:
				read_kmer_dict[readseq[i:i+KMER_LEN]] += [(read.id,i,True,True)]
			else:
				read_kmer_dict[readseq[i:i+KMER_LEN]] = [(read.id,i,True,True)]
#	print 'store read kmers',datetime.datetime.now()-kmerstarttime
	
	sequences = dict()
	for name in reference:
		sequences[name] = array2seq(reference[name])
	
	loglik,sequences = nEMsteps(int(sys.argv[5]),hmm,sequences,read_kmer_dict,reads)
	
	outfa = open(sys.argv[4],'w')
	for name in sequences:
		outfa.write('>'+name+' reconstruction\n'+sequences[name]+'\n')
	outfa.close()
	
if __name__ == '__main__':
	main()
					
		
