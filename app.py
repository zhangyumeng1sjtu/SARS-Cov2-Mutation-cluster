import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
from sklearn.cluster import KMeans, SpectralClustering, Birch
import matplotlib.cm as cm
from sklearn.preprocessing import StandardScaler,MinMaxScaler
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from PIL import Image
import geocoder
import random

st.title('Mutation clustering for SARS-CoV2')
image = Image.open("fig1.jpg")
st.image(image, width=1000, use_column_width=True)
data = pd.read_csv("data.csv", header=None)
data.columns = ["seq_id","type","start","end","ref","alt","mutation","pro_id"]
default = True

if st.checkbox("Show Raw Data"):
    st.dataframe(data,1500,300)

SNP_data = data[data["type"] == "SNP"]
pro_id = np.unique(SNP_data['pro_id'].astype(str))
option_1 = st.multiselect('Select SNPs on Your Interested Proteins',pro_id[:-1])
if option_1:
    filtered_data = SNP_data[SNP_data['pro_id'].isin(option_1)]
    st.write('You have selected {} SNPs on : {}'.format(len(filtered_data),str(option_1)))
    st.dataframe(filtered_data,1500,300)

    if len(filtered_data) > 200:
        new_data = filtered_data[filtered_data['mutation'].isin(['missense_variant','synonymous_variant'])]
        st.write("Focus on Two Key Variants.")
        st.dataframe(new_data['mutation'].value_counts())
    
        if st.checkbox("Show Mutation matrix"):
            matrix = pd.pivot_table(new_data,index=['ref'],columns=['alt'],values=['mutation'], aggfunc= ['count']).fillna(0)
            matrix.columns = matrix.columns.droplevel((0,1))
            st.table(matrix)
        df1 = pd.pivot_table(new_data,index=['seq_id'],columns=['mutation'],values=['start','end'],aggfunc= ['mean','std']).fillna(0)
        df2 = pd.pivot_table(new_data,index=['seq_id'],columns=['mutation'],values=['type'],aggfunc= ['count']).fillna(0)
        df = pd.concat([df1,df2],axis=1)
        default = False
        if st.checkbox("Show Correlation Heatmap"):
            plt.figure(figsize=(10,10))
            sns.heatmap(df.corr(), square=True,annot=True,cmap='viridis')
            plt.title("Correlation Heatmap")
            st.pyplot()

st.sidebar.header("Cluster Options")
if default == True:
    value = pd.read_csv("default.csv")[3:].values
    X = value[:,1:]
    default_idx = value[:,0]
else:
    X = df.values

option_2 = st.sidebar.selectbox("Preprocessing Method",('None','Nomalize','Standardize'))
if option_2 == 'Nomalize':
    X = MinMaxScaler().fit_transform(X)
elif option_2 == 'Standardize':
    X = StandardScaler().fit_transform(X)
else:
    pass

option_3 = st.sidebar.selectbox("Dimension Reduction Method",('t-SNE','PCA'))
if option_3 == 't-SNE':
    tsne=TSNE(n_components=2,init='random',random_state=10)
    X_new = tsne.fit_transform(X)
elif option_3 == 'PCA':
    pca = PCA(n_components=2, random_state=10)
    X_new = pca.fit_transform(X)
if st.checkbox("Visualize the Dimension Reduction Result"):
    plt.figure(figsize=(8, 6))
    plt.scatter(X_new[:,0], X_new[:,1])
    plt.title("{} Visualization Result".format(option_3))
    st.pyplot()

@st.cache
def train(X,method,num):
    cluster = method(n_clusters=num).fit_predict(X)
    return cluster 

def get_loc(loc_name):
    g = geocoder.arcgis(loc_name)
    return g.latlng
    
option_4 = st.sidebar.selectbox("Cluster Algorithm",('K-means','Birch','SpectralClustering'))
dic = {'K-means': KMeans,'Birch':Birch, 'SpectralClustering':SpectralClustering}
option_5 = st.sidebar.slider("Numbers of Clusters",2,8,3)


if st.checkbox('Let\'s Start the Clustering'): 
    cluster_ = train(X,dic[option_4],option_5)
    plt.figure(figsize=(8, 6))
    plt.scatter(X_new[:,0], X_new[:,1], c=cluster_)
    plt.title("{} Result for {} Clusters".format(option_4,option_5))
    st.pyplot()

    index =  default_idx if default else df.index
    label_info = pd.DataFrame(cluster_).set_index(index)

    # annotation = pd.read_csv("annotations.csv").set_index('Accession ID')
    # lat_list = []
    # lon_list = []
    # for loc in annotation['Location']:
    #     lat_list.append(get_loc(loc)[0])
    #     lon_list.append(get_loc(loc)[1])
    # annotation['lat'] = lat_list
    # annotation['lon'] = lon_list
    # annotation.to_csv("annotations.csv")
    annotation = pd.read_csv("annotations.csv").set_index('Accession ID')
        
    result = pd.concat([annotation, label_info],axis=1).reset_index().dropna(axis=0,how='any')
    result.columns = ['Accession ID','Sample Time','Location','lat','lon','Cluster Label']
    st.subheader("Cluster Info Table")
    st.dataframe(result)


    option_6 = st.sidebar.slider("Select one Cluster: ", 0, option_5-1,0)
    select_result = result[result['Cluster Label']==option_6]

    lat = np.array(select_result['lat']) + np.random.randn(len(select_result))/5-0.1
    lon = np.array(select_result['lon']) + np.random.randn(len(select_result))/5-0.1
    map_data = pd.DataFrame({'lat': lat, 'lon': lon})
    st.subheader("Mutation Location Map")
    st.map(map_data)

    

 


