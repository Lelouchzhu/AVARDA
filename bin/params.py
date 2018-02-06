# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 10:47:13 2017

@authors: Sanjay Kottapalli, Tiezheng Yuan
"""
import os
import pandas as pd
import numpy as np
from datetime import datetime as date
import flex_array

class file_IO:
	'''Class for parsing file input'''
	def __init__(self, infile, sep = '='):
		self.infile = infile
		self.sep = sep
		
	def file_to_dict(self):
		mydict = {}
		f = open(self.infile,'r')
		for line in f:
			line=line.strip()
			if 0 < line.find(self.sep) < len(line):
				key, value = line.split(self.sep)
				mydict[key] = value
		f.close()
		return mydict
	
	def dict_to_file(self, mydict):
		f = open(self.infile, 'w')
		for i in mydict:
			f.write(str(i)+self.sep+str(mydict[i])+'\n')
		f.close()
		return None
	
	def flat_file_to_df(self, df_index=[0,1,2]):
		aln_dict = {}
		#get index of rows, columns, and values
		row_index, col_index, value_index = df_index 
		#read file into nested dict
		f = open(self.infile, 'r')
		for line in f:
			line = line.strip()
			items = line.split(self.sep)
			row, col, value = items[row_index], items[col_index], items[value_index]
			if not col in aln_dict: aln_dict[col] = {}
			#assign value
			aln_dict[col][row] = float(value)
		f.close()
		#convert dict2 to df and fill NA with 0
		mydf = pd.DataFrame(aln_dict).fillna(0)
		return mydf
	
class param_dict:
	''''''
	def __init__(self, mydict):
		self.par = mydict
		
	def adjust_par(self):
		if 'organism' in self.par['file_aln']:
			self.par['organism'] = True
		else:
			self.par['organism'] = False
		if not self.par['dir_home'][-1] == '/':
			self.par['dir_home'] += '/'
		self.par['dir_input'] = self.par['dir_home'] + 'input/'
		self.par['dir_result'] = self.par['dir_home'] + 'results/'
		try:
			os.makedirs(self.par['dir_result'], exist_ok=False)
		except:
			pass
		self.par['dir_ref_seq'] = self.par['dir_home'] + 'ref_seq/'
		print("Home directory: " + self.par['dir_home'])
		print("Result directory: " + self.par['dir_result'])
		
		d = date.today()
		d = d.strftime("%d-%B-%Y_%H-%M")
		file_head = self.par['zscore_file'].split('.')[0]
		self.par['sub_dir'] = self.par['dir_result'] + file_head + ' ' + d + '/'
		os.makedirs(self.par['sub_dir'], exist_ok=True)
		
		self.par['graph_dir'] = self.par['sub_dir'] + 'sample_networks/'
		os.makedirs(self.par['graph_dir'], exist_ok=True)
		
		self.par['file_annotation'] = self.par['dir_ref_seq'] + self.par['file_annotation']
		self.par['file_aln'] = self.par['dir_ref_seq'] + self.par['file_aln']
		self.par['zscore_file'] = self.par['dir_input'] + self.par['zscore_file']
		
		if self.par['use_filter'].lower() == 'yes':
			self.par['use_filter'] = True
		elif self.par['use_filter'].lower() == 'no':
			self.par['use_filter'] = False
		
		#default parameters
		self.par['Z_threshold']=int(self.par['Z_threshold']) if 'Z_threshold' in self.par else 10
		self.par['p_threshold']=float(self.par['p_threshold']) if 'p_threshold' in self.par else 0.01
		self.par['x_threshold']=int(self.par['x_threshold']) if 'x_threshold' in self.par else 2
		self.par['bh_threshold']=float(self.par['bh_threshold']) if 'bh_threshold' in self.par else 0.05
		#
		return(self.par)
	
	def probability_ref(self):
		# Writes to file probability tables later used in binomial assessments 
		viral_peptidome = open(self.par['file_annotation'], 'r')
		peptide_lib = []
		next(viral_peptidome)
		for line in viral_peptidome:
			items = line.split('\t')
			peptide_lib.append(str(items[0]))
		viral_peptidome.close()
		peptide_lib = peptide_lib[:-1]
		peptide_lib.sort(key=int)
		
		binary_b = flex_array.binary_aln_df(self.par['file_aln'])
		binary_b = flex_array.array(binary_b).filter_aln()
		binary_b = binary_b.reindex(peptide_lib).fillna(0)
		
		virus_sums = binary_b.apply(np.count_nonzero, axis=0)
		first_round_prob = virus_sums/len(peptide_lib)
		first_round_prob.to_csv(self.par['dir_ref_seq']+"total_probabilities.csv", header=False, index=True)
		print("First probability file generated.")
		
		viruses = list(binary_b.columns)
		virus_intersections = pd.DataFrame(index=viruses, columns=viruses)
		for i in viruses:
			for j in viruses:
				a = binary_b[i]; b = binary_b[j]
				virus_intersections.loc[i,j] = np.dot(a,b)
		
		virus_unique = pd.DataFrame(index=viruses, columns=viruses)
		for i in virus_intersections.columns:
			virus_unique[i] = virus_sums[i] - virus_intersections[i]
		
		second_round_prob = pd.DataFrame(index=viruses, columns=viruses)
		for i in virus_intersections.index:
			for j in virus_intersections.columns:
				second_round_prob.loc[i,j] = virus_unique.loc[i,j]/(len(peptide_lib)-virus_intersections.loc[i,j])
		second_round_prob.to_csv(self.par['dir_ref_seq']+"unique_probabilities.csv", header=True, index=True)
		print("Second probability file generated.")
		
		third_round_prob = pd.DataFrame(index=viruses, columns=viruses)
		for i in virus_intersections.index:
			for j in virus_intersections.columns:
				third_round_prob.loc[i,j] = virus_intersections.loc[i,j]/(len(peptide_lib)-virus_unique.loc[i,j]-virus_unique.loc[j,i])
		third_round_prob.to_csv(self.par['dir_ref_seq']+"shared_probabilities.csv", header=True, index=True)
		print("Third (and last) probability file generated.")
		
		return None
		
# End