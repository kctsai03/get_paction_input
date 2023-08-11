#inputs a state tree and the k14 file and returns a dataframe with all the mutations for that state tree and the vaf's
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

import numpy as np

#this function uses elbow method to determine k (number of clusters) and returns the optimal number of cluters
def elbow(df_cluster):
  matrix_col = list(df_cluster.columns)
  matrix_col.remove('mutations')
  cluster_matrix = df_cluster[matrix_col].to_numpy()
  cluster_matrix.shape
  data = cluster_matrix
  num_clusters_to_test = min(6, len(df_cluster))

  # Generate random data points (assuming you already have the 'data' matrix)

  # Calculate cost (inertia) for different numbers of clusters
  costs = []
  for k in range(1, num_clusters_to_test + 1):
      kmeans = KMeans(n_clusters=k, random_state=42, n_init = 'auto')
      kmeans.fit(data)
      costs.append(kmeans.inertia_)

  # Calculate the second derivatives of the cost curve
  second_derivatives = np.gradient(np.gradient(costs))

  # Find the index of the point with the highest second derivative
  optimal_num_clusters = np.argmax(second_derivatives) + 1

  # Print the optimal number of clusters
  # print(f"Optimal number of clusters: {optimal_num_clusters}")

  # # Plot the elbow curve
  # plt.plot(range(1, num_clusters_to_test + 1), costs, marker='o')
  # plt.xlabel('Number of Clusters')
  # plt.ylabel('Cost (Inertia)')
  # plt.title('Elbow Curve')
  # plt.xticks(range(1, num_clusters_to_test + 1))
  # plt.show()
  return optimal_num_clusters

#input is the state tree and decifer output file
#returns the dataframe of all the mutations that have that state tree and the vaf's directly calculated
def get_df_cluster(state_tree, k14_output):
    col = ["sample" for i in range(9)]
    for i in range(9):
        col[i] = col[i] + str(i)
    col.insert(0, "mutations")

    df_cluster = pd.DataFrame(columns=col)

    for i in range(len(k14_output)):
        if k14_output['state_tree'][i] == state_tree:
            vaf = []
            vars = []
            tots = []
            for column in list(k14_output.columns):
                if "VAR" in column:
                    var = k14_output[column][i]
                    vars.append(var)
                if "TOT" in column:
                    tot = k14_output[column][i]
                    tots.append(tot)
            for j in range(len(vars)):
                vaf.append(vars[j]/tots[j])
            vaf.insert(0, k14_output['mut_index'][i])
            df_cluster.loc[len(df_cluster)] = vaf
    return df_cluster


#clusters the mutations and adds a column to df_cluster giving the cluster label 
#returns df_cluster and the cluster center vector
def cluster(num_clusters, df_cluster):
    matrix_col = list(df_cluster.columns)
    matrix_col.remove('mutations')
    cluster_matrix = df_cluster[matrix_col].to_numpy()
    # run k means
    kmeans = KMeans(n_clusters=num_clusters, random_state=0, n_init = "auto").fit(cluster_matrix)
    labels = kmeans.labels_
    cluster_centers = kmeans.cluster_centers_
    df_cluster['cluster_label'] = labels
    df_cluster.sort_values('cluster_label')
    # df_cluster['cluster_center'] = cluster_centers
    #this is the label for the first state tree)
    return df_cluster, cluster_centers

# def graph(df_cluster):
#   bins = [x/50 for x in range(0, 50, 1)]
#   for sample in df_cluster.columns:
#     if sample == 'mutations' or sample == 'cluster_label':
#       continue
#     fig, ax = plt.subplots(figsize =(5, 4))
#     ls_sample = list(df_cluster[sample])
#     ls_rm_0 = [i for i in ls_sample if i!= 0.0]
#     ax.hist(ls_rm_0, bins = bins)
#     plt.xlabel('vaf')
#     plt.ylabel('number of mutations')
#     plt.title(sample)



