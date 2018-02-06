# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 15:41:05 2017

@authors: Sanjay Kottapalli, Tiezheng Yuan
"""
import params
import pandas as pd
import numpy as np
import timeit
import networkx as nx
from scipy.stats import binom_test

def standard_df(infile):
	sep = ',' if infile.endswith('csv') else '\t'
	stand_df = pd.read_csv(infile, header=0, index_col=0, sep=sep, low_memory=False)
	stand_df.index=[str(i) for i in stand_df.index]
	stand_df.columns=[str(i) for i in stand_df.columns]
	stand_df.fillna(0, inplace=True)
	return stand_df

def binary_aln_df(infile):
	print('Reading alignment file: ' + infile)
	file_sep = '\t' if infile.endswith('txt') else ','
	aln_df = params.file_IO(infile, file_sep).flat_file_to_df([0,1,11])
	#convert to binary matrix
	binary_b=pd.DataFrame(np.where(aln_df > 0, 1, 0))
	binary_b.index=[str(i).split('_')[1] for i in aln_df.index]
	binary_b.columns=[str(i).replace(',', ';') for i in aln_df.columns]
	binary_b.fillna(0, inplace=True)
	return binary_b

def full_aln_df(infile):
	print('Reading alignment file: ' + infile)
	file_sep = '\t' if infile.endswith('txt') else ','
	aln_df = params.file_IO(infile, file_sep).flat_file_to_df([0,1,11])
	aln_df.index=[str(i).split('_')[1] for i in aln_df.index]
	aln_df.columns=[str(i).replace(',', ';') for i in aln_df.columns]
	aln_df.fillna(0, inplace=True)
	return aln_df

def sparse_aln_df(infile):
    print('Reading alignment file: ' + infile)
    aln_df = pd.read_pickle(infile)
    #aln_df = aln_df.to_dense()
    binary_b = pd.DataFrame(np.where(aln_df > 0, 1, 0), index=[str(i) for i in aln_df.index], columns=aln_df.columns)
    binary_b.fillna(0, inplace=True)
    return binary_b

class array:
	def __init__(self, array):
		self.array = array
	
	def filter_aln(self, min_alignments=10):
		time1 = timeit.default_timer()
		virus_sums = self.array.apply(np.sum, axis=0)
		viruses = list(self.array.columns)
		virus_intersections = pd.DataFrame(index=viruses, columns=viruses)
		for i in viruses:
			a = self.array[i]
			for j in viruses:
				b = self.array[j]
				virus_intersections.loc[i,j] = np.dot(a,b)

		shared_prob = virus_intersections.divide(virus_sums, axis='columns')		
		np.fill_diagonal(shared_prob.values, 0)
		
		#Create directed graph representation of virus dependencies
		shared_df = pd.DataFrame(columns=['child', 'parent'])
		child = []; parent = [];
		for i in shared_prob.columns:
			if 1.0 in shared_prob[i].values:
				parents = np.array(shared_prob.index)[np.where(shared_prob[i].values == 1.0)[0]]
				for j in parents:
					child.append(i)
					parent.append(j)
		shared_df['child'] = child
		shared_df['parent'] = parent
		G = nx.from_pandas_dataframe(shared_df, source='child', target='parent', create_using=nx.DiGraph())
		virus_lengths = {i:len(i) for i in G.nodes()}
		nx.set_node_attributes(G, name='Length', values=virus_lengths)
		#Remove complete subset viruses
		in_degree = pd.Series(dict(G.in_degree()))
		remove_viruses = list(in_degree.index[np.where(in_degree == 0)[0]])
		G.remove_nodes_from(remove_viruses)
		
		#Find cycles in G, remove all but shortest string virus in each cycle
		len_cycles = 1
		while len_cycles != 0:
			cycles = []
			for i in list(G.nodes()):
				try:
					cycles.append(nx.find_cycle(G, source=i))
				except:
					pass
			len_cycles = len(cycles)
			if len_cycles != 0:
				cycles = [np.array(i) for i in cycles]
				cycles = [np.ndarray.flatten(i) for i in cycles]
				cycles = [np.unique(i) for i in cycles]
				cycles = pd.DataFrame(cycles).drop_duplicates().values.tolist()
				lengths = nx.get_node_attributes(G, name='Length')
				for i in cycles:
					cycle_lengths = {k:lengths[k] for k in i}
					min_length = min(cycle_lengths)
					for j in i:
						if cycle_lengths[j] != min_length:
							remove_viruses.append(j)
							G.remove_node(j)
		
		#Remove viruses with most outgoing edges (denoting subset) iteratively
		if len(G.edges()) != 0:
			out_degree = pd.Series(dict(G.out_degree()))
			while max(out_degree) != 0:
				out_degree = pd.Series(dict(G.out_degree()))
				if max(out_degree) != 0:
					max_virus = out_degree.index[list(out_degree).index(max(out_degree))]
					#print(max_virus)
					remove_viruses.append(max_virus)
					G.remove_node(max_virus)
		'''
		f = open('removed_viruses.txt', 'w')
		for i in remove_viruses:
			f.write(i+'\n')
		f.close()
		'''
		self.array.drop(remove_viruses, axis=1, inplace=True)
		
		#Drop the appropriate viruses from the alignment matrix
		virus_sums = self.array.apply(np.sum, axis=0, raw=True, reduce=True)
		virus_sums_index = list(virus_sums.index)
		for i in range(len(virus_sums)):
			if virus_sums.iloc[i] <= min_alignments:
				self.array.drop(virus_sums_index[i], axis=1, inplace=True)
		
		time2 = timeit.default_timer()
		print("Time to filter the alignment matrix: " + str(time2-time1))
		
		return self.array
	
	def binom_reassign(self, input_num, ref_seq, p_threshold=0.01, hits_threshold=2, organism=False):
		#Read probability files
		if organism:
			org=ref_seq+'org_'
		else:
			org=ref_seq
		
		first_round_prob = pd.read_csv(org+"total_probabilities.csv", index_col=0, header=None, squeeze=True)
		second_round_prob = pd.read_csv(org+"unique_probabilities.csv", header=0, index_col=0)
		third_round_prob = pd.read_csv(org+"shared_probabilities.csv", header=0, index_col=0)
		
		#Function for defining significance
		def is_sig(p, n, x):
			p_value = binom_test(p=p, n=n, x=x, alternative='greater')
			if p_value < p_threshold and x > hits_threshold:
				return True
			else:
				return False
		
		#Series of p-values for each virus
		p_series = pd.Series(index=list(self.array.columns))
		
		#Number of hits to each virus
		ranked_hits = self.array.apply(sum, axis=0)
		orig_viruses = self.array.columns
		self.array = self.array[ranked_hits[ranked_hits>0].index]
		#Calculate p-values for each virus based on the number of total hits
		virus_pvalues_1 = pd.Series(index=list(self.array.columns))
		for i in virus_pvalues_1.index:
			virus_pvalues_1[i] = binom_test(p=first_round_prob[i], n=input_num, x=ranked_hits[i], alternative='greater')
		
		#Pre-greedy p-values for output
		orig_pseries = virus_pvalues_1.copy()
		
		#Sort viruses by initial p-value
		virus_pvalues_1.sort_values(inplace=True, ascending=True)
		#Sort virus hits series by p-value
		ranked_hits = ranked_hits[virus_pvalues_1.index]   
		specie = list(ranked_hits.index)
		n_rank = input_num
		
		while len(ranked_hits)>0 and ranked_hits.iloc[0]>0:
			#Comparing top hit to everything else (reassignments only)
			for i in specie[1:]:
				#Top hit virus (binary vector)
				highrank = self.array[specie[0]]
				high_num = sum(highrank)
				i_rank = self.array[i]
				i_num = sum(i_rank)
				#element-wise multiplication to get the vector of shared peptides
				shared_peps = np.multiply(i_rank, highrank)
				shared_num = sum(shared_peps)
				#Only do reassignment/sim tags if there is overlap between the two viruses
				if shared_num > 0:
					#If only one passes the threshold, reassign peptides accordingly
					if is_sig(second_round_prob.loc[i,specie[0]], n_rank-shared_num, high_num-shared_num) and not is_sig(second_round_prob.loc[specie[0],i], n_rank-shared_num, i_num-shared_num):#highrank_pvalue < p_threshold and i_pvalue >= p_threshold and high_num-shared_num > hits_threshold and i_num-shared_num > hits_threshold:
						#Subtract shared peptides from the insignificant virus
						self.array.loc[:,i] = i_rank-shared_peps
					elif not is_sig(second_round_prob.loc[i,specie[0]], n_rank-shared_num, high_num-shared_num) and is_sig(second_round_prob.loc[specie[0],i], n_rank-shared_num, i_num-shared_num):#highrank_pvalue >= p_threshold and i_pvalue < p_threshold and high_num-shared_num > hits_threshold and i_num-shared_num > hits_threshold:
						#Same as above
						self.array.loc[:,specie[0]] = highrank-shared_peps
			
			#Add p-value to series after any potential reassignments
			top_hit = specie[0]
			p_series[top_hit] = binom_test(p=first_round_prob[top_hit], n=n_rank, x=sum(self.array[top_hit]), alternative='greater')
			
			#Figure out how many peptides were globally unique to highrank
			high_peptides = np.where(self.array[specie[0]] == 1)[0]
			high_unique = len(np.where(self.array.iloc[high_peptides,:].apply(sum, axis=1) == 1)[0])
			
			#Now remove the highest hit since it will not be involved in subsequent comparisons
			specie.remove(specie[0])
			#Re-rank virus hits by total binomial
			ranked_hits = self.array.apply(sum, axis=0)
			ranked_hits = ranked_hits[specie]
			virus_pvalues_1 = pd.Series(index=specie)
			#Adjust the n used for binomial tests
			n_rank -= high_unique
			
			#Sort viruses by p-value
			for i in specie:
				virus_pvalues_1[i] = binom_test(p=first_round_prob[i], n=n_rank, x=ranked_hits[i], alternative='greater')
			virus_pvalues_1.sort_values(inplace=True, ascending=True)
			specie = list(virus_pvalues_1.index)
			ranked_hits = ranked_hits[specie]
		
		#Calculate sim tags after performing all peptide reassignments
		ranked_hits = self.array.apply(sum, axis=0)
		virus_pvalues_1 = pd.Series(index=list(self.array.columns))
		for i in virus_pvalues_1.index:
			virus_pvalues_1[i] = binom_test(p=first_round_prob[i], n=input_num, x=ranked_hits[i], alternative='greater')
		virus_pvalues_1.sort_values(inplace=True, ascending=True)
		ranked_hits = ranked_hits[virus_pvalues_1.index]
		
		#Sim tags for each virus, list of species to be examined
		specie=list(ranked_hits.index)
		sim_tag = pd.Series(float(0), index=specie)
		tag = 0
		
		n_rank = input_num
		
		while len(ranked_hits)>0 and ranked_hits.iloc[0]>0:
			for i in specie[1:]:
				highrank = self.array[specie[0]]
				high_num = sum(highrank)
				i_rank = self.array[i]
				i_num = sum(i_rank)
				shared_peps = np.multiply(i_rank, highrank)
				shared_num = sum(shared_peps)
				#Only do sim tags if there is overlap between the two viruses
				if shared_num > 0:
					#If neither passes (using shared test, symmetric probability table) and neither can stand on their own:
					if is_sig(third_round_prob.loc[i,specie[0]], n_rank-(high_num-shared_num)-(i_num-shared_num), shared_num) and not is_sig(second_round_prob.loc[i,specie[0]], n_rank-shared_num, high_num-shared_num) and not is_sig(second_round_prob.loc[specie[0],i], n_rank-shared_num, i_num-shared_num):
						#If neither of them already has a sim tag
						if np.array_equal(sim_tag[[specie[0],i]], [0,0]):
							tag += 1
							#Add a tag to the list for both viruses
							sim_tag.loc[[specie[0],i]] = 0.001*tag
						#If either already has a tag, assign that tag to the other
						elif sim_tag[specie[0]] != 0 and sim_tag[i] == 0:
							sim_tag.loc[i] = sim_tag[specie[0]]
						elif sim_tag[i] != 0 and sim_tag[specie[0]] == 0:
							sim_tag.loc[specie[0]] = sim_tag[i]
			
			#Figure out how many peptides were globally unique to highrank
			high_peptides = np.where(self.array[specie[0]] == 1)[0]
			high_unique = len(np.where(self.array.iloc[high_peptides,:].apply(sum, axis=1) == 1)[0])
			#Now remove the highest hit since it will not be involved in subsequent comparisons
			specie.remove(specie[0])
			#Re-rank virus hits by total binomial
			ranked_hits = self.array.apply(sum, axis=0)
			ranked_hits = ranked_hits[specie]
			virus_pvalues_1 = pd.Series(index=specie)
			n_rank -= high_unique
			for i in specie:
				virus_pvalues_1[i] = binom_test(p=first_round_prob[i], n=n_rank, x=ranked_hits[i], alternative='greater')
			#Sort viruses by p-value
			virus_pvalues_1.sort_values(inplace=True, ascending=True)
			#Sort virus hits series by p-value
			specie = list(virus_pvalues_1.index)
			ranked_hits = ranked_hits[specie]
		
		#Generate unique values matrix
		glob_unique = self.array.copy()
		for i in glob_unique.index:
			if sum(glob_unique.loc[i]) > 1:
				glob_unique.loc[i] -= glob_unique.loc[i]
		
		#Generate adjusted p-values for the p-value output Series (using the R package)
		#stats = importr('stats')
		#p_adjusted = stats.p_adjust(FloatVector(p_series.values), method='BH')​
		self.array = self.array.reindex(columns=orig_viruses, fill_value=0)
		glob_unique = glob_unique.reindex(columns=orig_viruses, fill_value=0)
		sim_tag = sim_tag.reindex(orig_viruses, fill_value=0)
		p_series = p_series.reindex(orig_viruses, fill_value=0)
		orig_pseries = orig_pseries.reindex(orig_viruses, fill_value=0)
		
		print("Done with reassignments!")
		return self.array, glob_unique, sim_tag, p_series, orig_pseries#, p_adjust_series
		
	def names_string(self, cutoff=0.001):
		hits = pd.Series(self.array)
		hits = hits[hits>=cutoff]
		names_str = ';'.join(list(hits.index))
		return names_str
	
	def gen_ind_hits(self, overlap_dict, graph_dir, sample_number):
		hits_series = pd.Series(self.array).astype(float)
		#Start time, for timing purposes
		start_time = timeit.default_timer()
		#Take only key-value pairs from overlap_dict for the sample hits
		hits_dict = hits_series.to_dict()
		peptide_hits = [str(i) for i in hits_dict.keys()]
		sub_dict = {i: overlap_dict[i] for i in peptide_hits if i in overlap_dict}
		#Generated a subgraph of relationships between all sample hits
		G = nx.Graph(sub_dict)
		for i in list(G.nodes()):
			if i not in peptide_hits:
				G.remove_node(i)
		#Add peptides that have no connections to others
		for i in peptide_hits:
			if i not in list(G.nodes()):
				G.add_node(i)
		#Add z-scores as attributes to the graph
		num_hits = len(hits_series)
		zscore = hits_series.to_dict()
		zscore = {i:float(zscore[i]) for i in zscore}
		nx.set_node_attributes(G, name='Z-Score', values=zscore)
		zscore = nx.get_node_attributes(G, name='Z-Score')
		
		edge = nx.get_edge_attributes(G, name='weight')
		edge = {i:float(edge[i]) for i in edge}
		nx.set_edge_attributes(G, name='weight', values=edge)
		
		#Write graphml file to appropriate directory
		#nx.write_graphml(G, graph_dir+sample_name.replace('.','-')+'.graphml')
		nx.write_graphml(G, graph_dir+'sample_' + str(sample_number) + '.graphml')
		
		#Reducing graph to max_degree 2
		if len(G.edges()) != 0:
			degree = dict(G.degree(nbunch=G.nodes()))
			degrees = [i for i in degree.values()]
			max_degree = np.max(degrees)
			while max_degree > 2:
				degree = dict(G.degree(nbunch=G.nodes()))
				vertices = [i for i in degree.keys()]
				degrees = [i for i in degree.values()]
				vertices = np.array(vertices); degrees = np.array(degrees);
				max_degree = np.max(degrees)
				#Remove nodes of highest degree (look for lowest z-score if tie)
				if max_degree > 2:
					max_vertex_indices = np.where(degrees==max_degree)
					max_vertices = vertices[max_vertex_indices]
					max_degree_scores = [zscore[i] for i in max_vertices]
					min_z = min(max_degree_scores)
					for i in max_vertices:
						if zscore[i] == min_z:
							G.remove_node(i)
							break #so that multiple nodes are not removed
				
		#Eliminates one vertex from each cycle (lowest z-score) to convert them into paths
		len_cycles = 1
		#While loop to make sure that there are no cycles left
		while len_cycles != 0:
			cycles = []
			for i in list(G.nodes()):
				try:
					cycles.append(nx.find_cycle(G, source=i))
				except:
					pass
			len_cycles = len(cycles)
			if len_cycles != 0:
				cycles = [np.array(i) for i in cycles]
				cycles = [np.ndarray.flatten(i) for i in cycles]
				cycles = [np.unique(i) for i in cycles]
				cycles = pd.DataFrame(cycles).drop_duplicates().values.tolist()
				for i in range(len(cycles)):
					cycles[i] = [str(j) for j in cycles[i] if str(j) in G.nodes()]
				node_zscores = nx.get_node_attributes(G, name='Z-Score')
				for i in cycles:
					cycle_scores=[node_zscores[k] for k in i]
					min_score = min(cycle_scores)
					for j in i:
						if node_zscores[j] == min_score:
							G.remove_node(j)
							break #otherwise it will eliminate two nodes in one cycle which both have the same z-score
		
		#Code for deleting vertices from paths based on even or odd length
		node_zscores = nx.get_node_attributes(G, name='Z-Score')
		degree = dict(G.degree(nbunch=G.nodes()))
		if len(G.edges()) != 0:
			degrees = [i for i in degree.values()]
			max_degree = max(degrees)
			while(max_degree > 0):
				components = list(nx.connected_component_subgraphs(G))
				for i in components:
					#even paths
					if len(i.nodes())%2 == 0:
						path_scores = [node_zscores[k] for k in i]
						min_score = min(path_scores)
						for j in i.nodes():
							if node_zscores[j] == min_score:
								G.remove_node(j)
								break #same as above
					#odd paths
					elif len(i.nodes())%2 == 1 and len(i.nodes()) != 1:
						endpoints=[]
						for j in i.nodes():
							if degree[j] == 1:
								endpoints.append(j)
						if len(endpoints)==2:
							path = nx.shortest_path(G, source=endpoints[0], target=endpoints[1])
							middle_indices = np.arange(1,len(path)-1,2)
							for i in middle_indices:
								G.remove_node(path[i])
				degree = dict(G.degree(nbunch=G.nodes()))
				degrees = [i for i in degree.values()]
				max_degree = max(degrees)
		
		num_nodes = len(G.nodes())
		proportion = np.divide(float(num_nodes), float(num_hits))
		print("Proportion of hits kept: " + str(proportion))
		
		#Creating pd.Series object from the nodes of the graph
		ind_hits_dict = {i: node_zscores[i] for i in G.nodes()}
		ind_hits_series = pd.Series(ind_hits_dict)
		#print("Number of peptides kept: " + str(len(ind_hits_series)))  
		
		end_time = timeit.default_timer()
		sample_time = end_time - start_time
		print("Time it took to remove overlaps: " + str(sample_time))
		
		#print("Number of edges in G: " + str(len(G.edges()))) (should be 0, always)
		
		return ind_hits_series
	
#End