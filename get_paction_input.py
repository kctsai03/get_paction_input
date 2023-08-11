#input a state tree, the mutation and the number of clusters and get the vaf array for all samples
#the vaf array has the average vaf for cluster that this mutation belongs to
from cluster import get_df_cluster, cluster, elbow
import pandas as pd

#input the state tree, the mutation, the number of clusters, and the decifer output file
#get the vaf array for that mutation (the CLUSTERED vaf array)
def get_vaf(state_tree, mutation, num_clusters, k14_output):
  df_cluster = get_df_cluster(state_tree, k14_output)
  df_cluster, cluster_centers = cluster(num_clusters, df_cluster)
  index = df_cluster[df_cluster['mutations'] == mutation].index[0]
  mut_cluster = df_cluster['cluster_label'][index]
  vaf = list(cluster_centers[mut_cluster, :])
  return vaf




##parse state tree data from file
#input:
#(1,1,0)->(1,1,1);(1,1,0)->(2,2,0);(1,1,0)->(4,0,0)
#output:
#[[['1', '1'], ['1', '1']], [['1', '1'], ['2', '2']], [['1', '1'], ['4', '0']]]
def parse_state_tree(state_tree_raw):
  state_tree_list = state_tree_raw.split(";")
  state_tree = []
  for i in state_tree_list:
    edge = i.split('->')

    edge_start = []
    edge_start.append(edge[0].split(',')[0][1])
    edge_start.append(edge[0].split(',')[1])
    edge_start.append(edge[0].split(',')[2][0])

    edge_end = []
    edge_end.append(edge[1].split(',')[0][1])
    edge_end.append(edge[1].split(',')[1])
    edge_end.append(edge[1].split(',')[2][0])

    edge_update = []
    edge_update.append(edge_start)
    edge_update.append(edge_end)

    state_tree.append(edge_update)
  return state_tree

#function that takes in
#state_tree (1,1,0)->(1,1,1);(1,1,0)->(2,1,0);(1,1,1)->(2,2,2);(2,2,2)->(3,2,2)
#mutation "10.52573509.G.GA"
#sample_to_pt {0: "MPAM06PT3", 1: "MPAM06PT4", etc.}
#vaf [0.1837795220329837,0.4155523820036714,0.12402229821802381, etc.]
#best_input (files with copy number clones and proportions for every bin)
#
#and returns a dataframe that gives us the genotype proportions
#  genotypes        s0        s1        s2        s3        s4        s5  \
# 0     1|1|0  0.735879  0.158843  0.832850  0.795854  0.349248  0.637473
# 1     1|1|1  0.033574  0.203305  0.038323  0.016625  0.183018  0.109941
# 2     2|1|0  0.000000  0.074337  0.000000  0.000000  0.044884  0.030000
# 3     2|2|2  0.000000  0.563515  0.000000  0.000000  0.422850  0.176182
# 4     3|2|2  0.230547  0.000000  0.128827  0.187521  0.000000  0.046404

#          s6        s7        s8
# 0  0.783001  0.475768  0.715027
# 1 -0.016127  0.145705  0.005792
# 2  0.000000  0.065034  0.000000
# 3  0.000000  0.000000  0.000000
# 4  0.233126  0.313493  0.279182
def get_genotype_prop(mutation, sample_to_pt, state_tree, best_input, num_clusters, op, k14_output):
  #get the number of samples we have for each mutation
  num_samples = 0
  for col in list(k14_output.columns):
    if "VAR" in col:
      num_samples += 1

  #get the vaf array (clustered)
  vaf = get_vaf(state_tree, mutation, num_clusters, k14_output)
  parsed_tree = parse_state_tree(state_tree)
  genotypes = []

  #create a list of all the genotypes seen in the genotype tree
  for edge in parsed_tree:
    for i in range(2):
      if edge[i] not in genotypes:
        genotypes.append(edge[i])

  #figure out what x*, y*, and m* are for (x*, y*, and m*)
  for i in range(len(genotypes)):
    for j in range(i, len(genotypes)):
      if genotypes[i][0] == genotypes[j][0] and genotypes[i][1] == genotypes[j][1] and i != j:
        x_star = genotypes[i][0]
        y_star = genotypes[i][1]
        m_star = genotypes[j][2]

  #locate the mutation in the copy number file
  chr_number = int(mutation.split('.')[0])
  chr_position = int(mutation.split('.')[1])
  #chr_x is the tsv file where chromosome number = x
  chr_x = best_input[best_input['#CHR'] == chr_number]
  #row_numbers will be the row index of the mutation bin for 52573509 (837-845)
  row_numbers = list(chr_x[(chr_x['START'] < chr_position) & (chr_x['END'] > chr_position)].index)

  #for each of the genotypes, add the inferred proportions
  for g in genotypes:
    geno = g[0] + "|" + g[1] + "|" + g[2]
    new_insert = [geno]
    #for each sample, we figure out the proportions of each copy number and calculate the genotype proportions according to the formula
    for index in range(num_samples):
      sample_name = sample_to_pt[index]  #sample name for i = 0 is PT3
      rows = best_input['SAMPLE'][row_numbers]
      row_number = rows[rows==sample_name].index[0] #840

      clone_proportions = {}
      cn_states = []
      cn_prop = []
      for column in best_input.columns:
        if "cn" in column:
          cn_states.append(best_input[column][row_number])
        if "u_" in column:
          cn_prop.append(best_input[column][row_number])
      for x in range(len(cn_states)):
        clone_proportions[cn_states[x]] = cn_prop[x]
      f = 0
      for key in clone_proportions.keys():
        xy = int(key.split('|')[0]) + int(key.split('|')[1])
        f += xy*clone_proportions[key]
      miu_star = clone_proportions[x_star + "|" + y_star]
      sigma = 0
      for genotype in genotypes:
        if genotype[0] == x_star and genotype[1] == y_star:
          continue
        sigma += int(genotype[2]) * clone_proportions[genotype[0] + "|" + genotype[1]]
      if miu_star != 0:
        lm = 1/(int(m_star)*miu_star) * (vaf[index]*f-sigma)


      if g[0] == x_star and g[1] == y_star:
        #if the proportion of copy number is 0, then we append 0 (avoid divide by 0 error)
        if miu_star == 0:
          new_insert.append(0)
        elif g[2] == '0':
          new_insert.append((1-lm)*miu_star)
        else:
          new_insert.append(lm*miu_star)
      else:
        new_insert.append(clone_proportions[g[0] + "|" + g[1]])
    op.loc[len(op)] = new_insert
  return op


#input the decifer output file, the copy number input file (from hatchet), the mutation 
#returns the dataframe of genotypes and the proportions for each sample
def get_paction_inputs(k14_output, best_input, mutation, sample_to_pt):
  num_samples = 0
  for col in list(k14_output.columns):
    if "VAR" in col:
      num_samples += 1
  col = ["s" for i in range(num_samples)]
  for i in range(num_samples):
    col[i] = col[i] + str(i)
  col.insert(0, "genotypes")
  row = ["g" for i in range(5)]
  for i in range(5):
    row[i] = row[i] + str(i)
  op = pd.DataFrame(columns=col)

  state_tree = k14_output[k14_output['mut_index'] == mutation]['state_tree'].values[0]

  df_cluster = get_df_cluster(state_tree, k14_output)
  num_clusters = elbow(df_cluster)
  df_genotype_prop = get_genotype_prop(mutation, sample_to_pt, state_tree, best_input, num_clusters, op, k14_output)
  return df_genotype_prop

#input is k14_output and the mutation
#returns the dataframe with the edges of the genotype tree for that mutation 
def get_mut_tree(k14_output, mutation):
  state_tree = parse_state_tree(k14_output[k14_output['mut_index'] == mutation]['state_tree'].values[0])
  col = ['edge_start', 'edge_end']
  output = pd.DataFrame(columns=col)
  for edge in state_tree:
    row = []
    for i in range(2): 
      row.append(edge[i][0] + "|" + edge[i][1] + "|" + edge[i][2])
    output.loc[len(output)] = row
  return output


k14_output = pd.read_csv('/Users/kyletsai/Desktop/SUMMER_2023/decifer/RA17_22_output/MPAM06_output_K14.tsv', sep='\t')
best_input = pd.read_csv('/Users/kyletsai/Desktop/SUMMER_2023/decifer/RA17_22_input/best.seg.ucn.tsv', sep='\t')
sample_to_pt = {}
#sample_to_pt is a dictionary where the key is the sample number and the value is the name of the sample (such as MPAM06PT3). 
# note for palash: i wasn't sure how to generalize this since there wasn't a very nice way of getting this from the decifer_input file
# so i just hard coded it - but this is so that the sample we are looking at in the copy number input file matches with the sample 
# number that we're looking at from the decifer output file
sample_to_pt[0] = "MPAM06PT3"
sample_to_pt[1] = "MPAM06PT4"
sample_to_pt[2] = "MPAM06PT5"
sample_to_pt[3] = "MPAM06PT7"
sample_to_pt[4] = "MPAM06PT8"
sample_to_pt[5] = "MPAM06PT9"
sample_to_pt[6] = "MPAM06PT1"
sample_to_pt[7] = "MPAM06PT10"
sample_to_pt[8] = "MPAM06PT2"

test = get_paction_inputs(k14_output, best_input,'10.52573509.G.GA', sample_to_pt)
x = get_mut_tree(k14_output, '10.52573509.G.GA')
print(test)
print(x)

# x.to_csv("mut1_tree.csv", sep=',', index=False, header=False)
# test.to_csv("mut1.csv", sep=',', index=False) 
