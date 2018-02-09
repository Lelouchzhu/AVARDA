# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 15:39:38 2017

@authors: Sanjay Kottapalli, Tiezheng Yuan
"""
import flex_array
import params
import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests

class phip:
	def __init__(self, par):
		self.par = par
		self.dependent_peptides()
		self.run_analysis()
	
	def dependent_peptides(self):
		self.dependent_pep = {}
		file_dependent = self.par['dir_ref_seq']+'virus_dependent_peptides_trunc.csv'
		f = open(file_dependent, 'r')
		for line in f:
			line = line.strip().split(',')
			pep1 = str(line[0])
			pep2 = str(line[1])
			if pep1 in self.dependent_pep:
				self.dependent_pep[pep1].append(pep2)
			else:
				self.dependent_pep[pep1] = [pep2]
		f.close()
		return self.dependent_pep
	
	def run_analysis(self):
		zdf = flex_array.standard_df(self.par['zscore_file'])
		
		binary_b = flex_array.sparse_aln_df(self.par['file_aln'])
		binary_b = flex_array.array(binary_b).filter_aln()
		binary_b=binary_b.reindex(zdf.index).fillna(0)
	
		sum_df = pd.DataFrame(0, index=list(binary_b), columns=list(zdf))
		glob_unique = pd.DataFrame(0, index=list(binary_b), columns=list(zdf))
		pep_df = pd.DataFrame(np.nan, index=list(binary_b), columns=list(zdf))
		p_df = pd.DataFrame(index=list(binary_b), columns=list(zdf))
		padjust_df = pd.DataFrame(index=list(binary_b), columns=list(zdf))
		orig_p = pd.DataFrame(index=list(binary_b), columns=list(zdf))
		
		hits_series = pd.Series(index=list(zdf))
		nonoverlap_hits_series = pd.Series(index=list(zdf))
		samples = list(zdf.columns)
		
		nonoverlap_dict = {}
		
		for sample_name, column in zdf.iteritems():
			hits = column[column>=self.par['Z_threshold']].copy()
			if self.par['use_filter']:
				nonoverlap_hits = flex_array.array(hits).gen_ind_hits(self.dependent_pep, 
											  self.par['graph_dir'], samples.index(sample_name))
			elif not self.par['use_filter']:
				nonoverlap_hits = hits.copy()
			input_num = len(nonoverlap_hits)
			hits_series[sample_name] = len(hits)
			nonoverlap_hits_series[sample_name] = input_num
			nonoverlap_dict[sample_name] = list(nonoverlap_hits.index)
			print("%s:\thits=%s, nonoverlapped=%s" %(sample_name, len(hits), input_num))
			
			if input_num > 0:
				zb_df=binary_b.loc[nonoverlap_hits.index]
				collapse_zb, glob_array, sim_tag, p_series, orig_pseries = flex_array.array(zb_df).binom_reassign(
						input_num, self.par['dir_ref_seq'], self.par['p_threshold'], self.par['x_threshold'], self.par['organism'])
				sum_df[sample_name]=collapse_zb.apply(sum, axis=0) + sim_tag
				glob_unique[sample_name] = glob_array.apply(sum, axis=0) + sim_tag
				pep_df[sample_name]=collapse_zb.apply(lambda x: flex_array.array(x).names_string(0.001),axis=0)
				p_df[sample_name]=p_series
				orig_p[sample_name]=orig_pseries
			
		file_head = self.par['sub_dir'] + self.par['zscore_file'].split('/')[-1].split('.')[0] #Removes file path and extension
		if self.par['organism']:
			file_head += '_organism_'
		else:
			file_head += '_species_'
			
		#Write log file
		params.file_IO(self.par['sub_dir']+'parameters.log', sep='=').dict_to_file(self.par)
		
		#Write analysis files
		sum_df.to_csv(file_head+'total-counts.txt', sep='\t', header=True, index_label='Specie')
		glob_unique.to_csv(file_head+'unique-counts.txt', sep='\t', header=True, index_label='Specie')
		pep_df.to_csv(file_head+'peptides.txt', sep='\t', header=True, index_label='Specie')
		p_df.to_csv(file_head+'p-values.txt', sep='\t', header=True, index_label='Specie')
		orig_p.to_csv(file_head+'orig-p-values.txt', sep='\t', header=True, index_label='Specie')
		
		for i in p_df.columns:
			pvals=np.array(p_df[i].values)
			if not pd.isnull(pvals).all():
				mask = [j for j in np.where(np.isfinite(pvals))[0]]
				pval_corrected = np.empty(pvals.shape)
				pval_corrected.fill(np.nan)
				pval_corrected[mask] = multipletests(pvals[mask], method='fdr_bh')[1]
				padjust_df[i] = pval_corrected
		padjust_df.to_csv(file_head+'p-adjusted.txt', sep='\t', header=True, index_label='Specie')
		
		#Write independent peptides file
		f = open(self.par['sub_dir']+'independent_peptides.txt', 'w')
		for i in samples:
			f.write(i)
			for j in nonoverlap_dict[i]:
				f.write('\t' + str(j))
			f.write('\n')
		f.close()

		#Write summary file
		f = open(self.par['sub_dir']+'results_summary.txt', 'w')
		f.write("Sample name\tVirus\tBH p-value\tRaw p-value\tOrig p-value\tAssigned counts\t")
		f.write("Assigned peptides\tTotal sample hits\tTotal filtered sample hits\n")
		for i in samples:
			BH = padjust_df[i]
			BH = BH[BH < self.par['bh_threshold']]
			p_value = p_df[i]
			p_value = p_value[BH.index]
			orig_pvalue = orig_p[i]
			orig_pvalue = orig_pvalue[BH.index]
			counts = sum_df[i]
			counts = counts[BH.index]
			peptides = pep_df[i]
			peptides = peptides[BH.index]
			
			for j in BH.index:
					if counts[j] > self.par['x_threshold']:
						f.write(i+'\t')
						f.write(j+'\t'+str(BH[j])+'\t')
						f.write(str(p_value[j])+'\t'+str(orig_pvalue[j])+'\t')
						f.write(str(counts[j])+'\t'+str(peptides[j])+'\t')
						f.write(str(hits_series[i])+'\t'+str(nonoverlap_hits_series[i])+'\n')
		f.close()
		print("End of run.")
		#End
		