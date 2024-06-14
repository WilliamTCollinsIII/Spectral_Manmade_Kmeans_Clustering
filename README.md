# Spectral_Manmade_Kmeans_Clustering
This Repository is different methods of Kmeans clustering to try and group manmade materials with similar spectral components together. For this research, the samples of materials come from two data sets (ASTER & KLUM), and all material with less than 3 total samples for a specific material have been removed (i.e. all materials have at least 3 samples with similar descriptions). Below are all of the main files that perform a specific kmeans clustering technique.

Kmeans - This script performs kmeans clustering on 3 components from PCA. 

KmeansAllComponents - This script performs kmeans clustering on all features (wavelengths) available. You can change the wavelengths to try by simply changing the wavelengths variable in the script.

KmeansAllComponentsSMOTE - This script performs kmeans clustering on all features (wavelengths) available, but also creates synthetic data via SMOTE. 

ConstrainedKmeansClustering - This script performs a constrained kmeans clustering on 3 components from PCA. The constraint is groups must be clustered in the same cluster. The groups are created by samples that have the same file name. Some samples belong to any group and kmeans clsutering with cluster them to the closest cluster/centroid.

ConstrainedKmeansClustering2 - This script performs a constrained kmeans clustering on 3 components from PCA. The constraint is groups must be clustered in the same cluster. The groups are created by samples that have a similar keyword. The 12 keywords cover all of the samples therefore every sample belongs to a group and will stay with its group when clustering.

WardsClustering - This scipt performs Wards hierarchical clustering method on 3 components from PCA. In this method, clustering is achieved by first places each sample in its own cluster. Then compute the distance between all clusters, merge clusters that would create the smallest increase in variance, and continue this process until the desired number of clusters is reached.
