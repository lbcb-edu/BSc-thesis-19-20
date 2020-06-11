import sys
import re
from sklearn.feature_extraction import DictVectorizer
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from itertools import cycle
from sklearn.metrics import silhouette_score

sam_file_name = sys.argv[1]
ref_file_name = sys.argv[2]

ref_file = open(ref_file_name, 'r', encoding='utf8')
line = ref_file.readline()
reference = ref_file.readline()
ref_file.close()

reads_list = []
measurements = []
explained_modif_per_read = []

modification = {}
no_of_reads = 0

sam_file = open(sam_file_name, 'r', encoding='utf8')
line = sam_file.readline()
while line.startswith('@'):
    line = sam_file.readline()
while line != '':
    parts = re.split('\s+', line)
    read_name = parts[0]
    mapped = parts[1]
    ref_name = parts[2]
    position = int(parts[3])
    cigar = parts[5]
    read = parts[9]
    ref_copy = reference[position - 1:]

    if mapped == '4' or ref_name != 'tetra':
        line = sam_file.readline()
        continue

    if read_name in reads_list:
        line = sam_file.readline()
        continue

    if position > 50:
        line = sam_file.readline()
        continue

    reads_list.append(read_name)
    map_for_measurements = {}
    map_for_explained_mod = {}
    start_position = position
    no_of_reads += 1

    m = re.match('^([0-9]+S).*', cigar)
    if m is not None:
        remove = len(m.group(1))
        remove_read = int(m.group(1)[:remove - 1])
        read = read[remove_read:]
        cigar = cigar[remove:]
    m = re.match('.*[MID]([0-9]+S)$', cigar)
    if m is not None:
        remove = len(m.group(1))
        remove_read = int(m.group(1)[:remove - 1])
        read = read[:len(read) - remove_read]
        cigar = cigar[:len(cigar) - remove]
    while len(cigar) != 0:
        m = re.match('^([0-9]+[MID]).*', cigar)
        if m is None:
            print('Error!')
            exit()
        if m.group(1).endswith('I'):
            remove_read = int(m.group(1)[:len(m.group(1)) - 1])
            
            if position >= 50 and position <= 408:
                map_for_measurements[position] = True
                map_for_explained_mod[position] = ('I', read[0:remove_read])
                modification[position] = modification.get(position, 0) + 1

            read = read[remove_read:]
            cigar = cigar[len(m.group(1)):]
        elif m.group(1).endswith('D'):
            remove_ref = int(m.group(1)[:len(m.group(1)) - 1])

            for i in range(remove_ref):

                if position >= 50 and position <= 408:
                    map_for_measurements[position] = True
                    map_for_explained_mod[position] = ('D') 
                    modification[position] = modification.get(position, 0) + 1

                position += 1
            ref_copy = ref_copy[remove_ref:]
            cigar = cigar[len(m.group(1)):]
        elif m.group(1).endswith('M'):
            number = int(m.group(1)[:len(m.group(1)) - 1])
            for i in range(number):
                if ref_copy[i] != read[i]:

                    if position >= 50 and position <= 408:
                        map_for_measurements[position] = True
                        map_for_explained_mod[position] = ('M', read[i])
                        modification[position] = modification.get(position, 0) + 1

                position += 1
            ref_copy = ref_copy[number:]
            read = read[number:]
            cigar = cigar[len(m.group(1)):]

    if position < 408:
        reads_list.remove(read_name)
        line = sam_file.readline()
        continue

    measurements.append(map_for_measurements)
    explained_modif_per_read.append(map_for_explained_mod)

    line = sam_file.readline()

positions = []
for p in range(50, 409):
    if p in modification:
        if modification[p]/no_of_reads > 0.1:
            positions.append(p)
positions = sorted(positions)

for i in range(len(measurements)):
    for p in range(50, 409):
        if p in positions and p not in measurements[i]:
            measurements[i][p] = False
        if p not in positions and p in measurements[i]:
            del measurements[i][p]

#for read_map in measurements:
#    print(read_map)

vec = DictVectorizer()
X = vec.fit_transform(measurements).toarray()


n_clusters = 4
kmeans_model = KMeans(n_clusters=n_clusters)
kmeans = kmeans_model.fit(X)

### PCA!!!!
def plot_figure_pca():
    pca = PCA(n_components=2)
    new_X = pca.fit_transform(X)

    k_means_cluster_centres = kmeans.cluster_centers_
    new_cluster_centres = pca.transform(k_means_cluster_centres)
    
    fig = plt.figure(1)
    k_means_labels = kmeans.labels_
    colors =  cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
    for k, col in zip(range(n_clusters), colors):
        my_members = k_means_labels == k
        cluster_center = new_cluster_centres[k]
        plt.plot(new_X[my_members, 0], new_X[my_members, 1], col+'.')
        plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
                markeredgecolor='k', markersize=10)

    plt.title('PCA - Number of clusters: %d' % n_clusters)
    plt.show()

plot_figure_pca()


### ELBOW METHOD!!!!
def elbow_method():
    distortions = []
    K = range(1,10)
    for n_clusters in K: 
        kmeans_model = KMeans(n_clusters=n_clusters)
        kmeans = kmeans_model.fit(X)
        distortions.append(kmeans_model.inertia_)

    plt.figure(figsize=(16,8))
    plt.plot(K, distortions, 'bx-')
    plt.xlabel('k')
    plt.ylabel('Distortion')
    plt.title('The Elbow Method showing the optimal k')
    plt.show()


### PRINT DATA ABOUT MODIFICATIONS
def print_data_about_clusters():
    for i in range(len(reads_list)):
        print(reads_list[i], kmeans.labels_[i])

    number_of_elems_in_cluster = {}
    clusters_performance = []
    for k in range(n_clusters):
        list_of_modif = {}
        for i in range(len(reads_list)):
            if kmeans.labels_[i] == k:
                number_of_elems_in_cluster[k] = number_of_elems_in_cluster.get(k, 0) + 1
                for p in explained_modif_per_read[i]:
                    if p not in list_of_modif:
                        list_of_modif[p] = {}
                    modif = explained_modif_per_read[i][p]
                    list_of_modif[p][modif] = list_of_modif[p].get(modif, 0) + 1
        clusters_performance.append(list_of_modif)

    cnt = 0
    for list_of_modif in clusters_performance:
        print('cluster', cnt)
        print(number_of_elems_in_cluster[cnt])
        map_of_significant = {}
        sorted_pos = sorted(list_of_modif.keys())
        for p in sorted_pos:
            for mod in list_of_modif[p]:
                if list_of_modif[p][mod]/number_of_elems_in_cluster[cnt] > 0.5:
                    print('{}\t{}\t{:.1f}'.format(p, mod, list_of_modif[p][mod]/number_of_elems_in_cluster[cnt]*100))
            print(p, list_of_modif[p])
        print()
        cnt += 1


### TSNE!!!!!!
def plot_figure_tsne():
    tsne = TSNE(n_components=2)
    X_embedded = tsne.fit_transform(X)

    plt.figure(1)
    plt.clf()
    k_means_labels = kmeans.labels_ 
    colors =  cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
    for k, col in zip(range(n_clusters), colors):
        my_members = k_means_labels == k
        plt.plot(X_embedded[my_members, 0], X_embedded[my_members, 1], col+'.')

    plt.title('TSNE: k = ' + str(n_clusters))
    plt.show()

