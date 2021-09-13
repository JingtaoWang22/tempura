import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import networkx as nx
import pandas as pd
from collections import defaultdict
import os
import scanpy
from time import time
from sys import getsizeof

################################################################################
#check if it possible to just load relevant parts of database
if 'adata' not in globals():
    print('We are loading the data now')
    duration = time()
    adata = scanpy.read_h5ad("database-2300.h5ad")
    observations = adata.obs
    print('Done loading data. It took', round(time()-duration,2),'seconds\n')
    print('The size of adata is:',round(getsizeof(adata)/1000000,2),'mb')

def extractPath(path):
    filenames = os.listdir(path)
    filenames.sort()
    clusters = {}
    for each in filenames:
        clusters.update({each:[[],[],[],[]]})
        temp = each.split('.')[0].split('-')
        for i,item in enumerate(temp):
            temp1 = item.split('n')
            clusters[each][i] = temp1
    return clusters

idremPaths = extractPath('./assets/')

#makes it easier to create the dropdowns if the cells are saved in arrays:
controlCells = []
IPF1Cells = []
IPF2Cells = []
IPF3Cells = []
for j in range(len(adata.uns['clusterType'])):
    for i in range(len(adata.uns['clusterType'][str(j)])):
        cell = '-'.join([str(j), str(i), adata.uns['clusterType'][str(j)][i]])
        if j == 0:
            controlCells.append(cell)
        elif j == 1:
            IPF1Cells.append(cell)
        elif j == 2:
            IPF2Cells.append(cell)
        else:
            IPF3Cells.append(cell)


#the following saves all the top genes in each stage to an array
print('Finding top Genes', end=' ')
duration = time()
topGenes = {}
for stage in ['1','2','3','0']:
    for index,cellsTopGenes in enumerate(adata.uns['topGene'][stage]):
        for gene in cellsTopGenes:
            coordinates = [int(stage),index]
            if gene not in topGenes.keys():
                topGenes[gene]=[coordinates]
            elif coordinates not in topGenes[gene]:
                topGenes[gene].append(coordinates)
print('Finished finding top Genes, it took', round(time()-duration,2),'seconds\n')

def find_cell_name(coordinates):
    stage, index = coordinates
    return [controlCells,IPF1Cells,IPF2Cells,IPF3Cells][stage][index]

def create_graph():
    G = nx.Graph()
    for j in range(len(adata.uns['clusterType'])):
        for i in range(len(adata.uns['clusterType'][str(j)])):
            node_index = 100*j
            node_index += i
            cell = '-'.join([str(j),str(i),adata.uns['clusterType'][str(j)][i]])
            if j>0:
                G.add_nodes_from([(node_index, {"pos": [j, 1.7*i],'stage': node_index//100,'index': i,'cellType':cell}), ])
            else:
                G.add_nodes_from([(node_index, {"pos": [j, i], 'stage': node_index // 100, 'index': i,'cellType':cell}), ])
    for key in sorted(adata.uns['edges']):
        for edge in adata.uns['edges'][key]:
            G.add_edge(int(key) * 100 + edge[0], (int(key) + 1) * 100 + edge[1])
    return G

#the following dictionary contains the cells as keys and an array of cells that they have edges to. This makes the path selection easier.
paths = defaultdict(list)

def add_path(startingNode,endNode, paths):
    #check if both are in, and add them.
    if type(startingNode) != str or type(endNode) != str:
        raise TypeError
    if endNode not in paths[startingNode]:
        paths[startingNode].append(endNode)
    if startingNode not in paths[endNode]:
        paths[endNode].append(startingNode)

def make_edge(x, y, text, width,color='cornflowerblue'):
    return  go.Scatter(x         = x,
                       y         = y,
                       line      = dict(width = width,
                                   color = color),
                       hoverinfo = 'skip',
                       text      = text,
                       mode      = 'lines')

def get_paths():
    """
    the stucture of the returned list is as follows: list[index] returns a list containing information about the path that has
    the cluster (0,index) as root. list[index][2] is a list that contains the indices of all clusters in stage 2 that belong to
    the path that has (0,index) as a root.
    :return: list
    """
    pathList = []
    for index in range(len(adata.uns['clusterType']['0'])):
        pathList.append([])
        for stage in range(len(adata.uns['clusterType'])):
            pathList[index].append([index]) if stage == 0 else pathList[index].append([])
    for stage in sorted(adata.uns['edges']):
        for edge in adata.uns['edges'][stage]:
            for pathIndex in range(len(pathList)):
                if edge[0] in pathList[pathIndex][int(stage)]:
                    pathList[pathIndex][int(stage) + 1].append(edge[1])
    return pathList

pathList = get_paths()

def create_cluster_coordinates(path):
    coordinates = []
    for stage in range(len(path)):
        coordinates.extend([[stage,index] for index in path[stage]])
    return coordinates

def get_list_clusters(pos):
    stage = pos[0]
    index = pos[1]
    for path in pathList:
        if index in path[stage]:
            return create_cluster_coordinates(path)

def create_cluster_path_dict():
    """
    :return: dict with cluster (stage,index) tuple as key and list of tuples in the same path as value
    """
    paths = get_paths()

def create_graph_fig(G,selected_cells=0):
    #This function greates the no.fig object of the graph

    #creates scatter plot containing the edges
    edge_x = []
    edge_y = []
    edge_trace = []

    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        #updating the dictionary
        add_path(G.nodes[edge[0]]['cellType'],G.nodes[edge[1]]['cellType'],paths)
        cellPos1 = [G.nodes[edge[0]]['stage'],G.nodes[edge[0]]['index']]
        cellPos2 = [G.nodes[edge[1]]['stage'],G.nodes[edge[1]]['index']]
        text = G.nodes[edge[0]]['cellType'] + '--' + G.nodes[edge[1]]['cellType']
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)
        if selected_cells and cellPos1 in selected_cells and cellPos2 in selected_cells:
            edge_trace.append(make_edge([x0,x1,None],[y0,y1,None],text,3,'#000000'))
        else:
            edge_trace.append(make_edge([x0,x1,None],[y0,y1,None],text,0.5,'#808080'))

    #creates scatter plot containing the nodes
    node_x = []
    node_y = []
    node_text = []
    node_color=[]
    node_opacity=[]
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_text.append(G.nodes[node]['cellType'])
        nodePos = [G.nodes[node]['stage'],G.nodes[node]['index']]
        if selected_cells and nodePos not in selected_cells:
            node_opacity.append(0.3)
        else:
            node_opacity.append(1)
        node_color.append(G.nodes[node]['stage'])
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            # colorscale options
            #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale='Blackbody',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Stage',
                xanchor='left',
                titleside='right',
                dtick=1,
                tick0=0,
            ),
            line_width=2))

    node_trace.text = node_text
    node_trace.marker.color = node_color
    node_trace.marker.opacity = node_opacity

    #creates the figure object using the two scatter plots of nodes and edges from above
    fig = go.Figure(data=node_trace,
                    layout=go.Layout(
                        title='<br>Progression of cell behaviour in IPF patients',
                        titlefont_size=16,
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20,l=5,r=5,t=40),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        height=600,
                        #width=1200,
                    )
                    )
    for trace in edge_trace:
        fig.add_trace(trace)

    fig.add_annotation(
        x=3,  # arrows' head
        y=-1,  # arrows' head
        ax=0,  # arrows' tail
        ay=-1,  # arrows' tail
        xref='x',
        yref='y',
        axref='x',
        ayref='y',
        text='',
        showarrow=True,
        arrowhead=5,
        arrowsize=1,
        arrowwidth=2,
        arrowcolor='#578018'
    )
    fig.add_annotation(
        x=1.5,
        y=-2,
        text='Progression of IPF',
        showarrow=False,
    )

    return fig

def create_graph_fig_cellGenes(G,selected_cells=0):
    #This function greates the no.fig object of the graph

    #creates scatter plot containing the edges
    edge_x = []
    edge_y = []
    edge_trace = []

    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        #updating the dictionary
        add_path(G.nodes[edge[0]]['cellType'],G.nodes[edge[1]]['cellType'],paths)
        cell1 = G.nodes[edge[0]]['cellType']
        cell2 = G.nodes[edge[1]]['cellType']
        text = cell1 + '--' + cell2
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)
        edge_trace.append(make_edge([x0, x1, None], [y0, y1, None],text, 0.5, '#808080'))

    #creates scatter plot containing the nodes
    node_x = []
    node_y = []
    node_text = []
    node_color=[]
    node_opacity=[]
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_text.append(G.nodes[node]['cellType'])
        nodePos = [G.nodes[node]['stage'],G.nodes[node]['index']]
        if selected_cells and nodePos not in selected_cells:
            node_opacity.append(0.3)
        else:
            node_opacity.append(1)
        node_color.append(G.nodes[node]['stage'])
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            colorscale='Blackbody',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Stage',
                xanchor='left',
                titleside='right',
                dtick=1,
                tick0=0,
            ),
            line_width=2))

    node_trace.text = node_text
    node_trace.marker.color = node_color
    node_trace.marker.opacity = node_opacity

    #creates the figure object using the two scatter plots of nodes and edges from above
    fig = go.Figure(data=node_trace,
                    layout=go.Layout(
                        title='<br>Progression of cell behaviour in IPF patients',
                        titlefont_size=16,
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20,l=5,r=5,t=40),
                        annotations=[ dict(
                            showarrow=False,
                            xref="paper", yref="paper",
                            x=0.005, y=-0.002 ) ],
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        height=600,
                        width=1200,
                    )
                    )
    for trace in edge_trace:
        fig.add_trace(trace)

    fig.add_annotation(
        x=3,  # arrows' head
        y=-1,  # arrows' head
        ax=0,  # arrows' tail
        ay=-1,  # arrows' tail
        xref='x',
        yref='y',
        axref='x',
        ayref='y',
        text='',  # if you want only the arrow
        showarrow=True,
        arrowhead=5,
        arrowsize=1,
        arrowwidth=2,
        arrowcolor='#578018'
    )
    fig.add_annotation(
        x=1.5,
        y=-2,
        text='Progression of IPF',
        showarrow=False,
    )
    return fig

def load_cluster_gene_expr_matrix(stage,index):
    """
    :param stage: char representing the stage ('c',"1",'2' or '3')
    :param index: int representing the index of the cluster in a stage
    :return: gene expression matrix of the specified cluster
    """
    index = int(index)
    stage = 0 if stage == 'c' else int(stage)

    print('We are doing filtering',(stage,index), 'now',end='    ')
    duration = time()
    g = observations[observations.leiden.eq(index) & observations.stage.eq(stage)].index.tolist()
    array = adata[g].X.toarray()
    print('It took',round(time()-duration,2),'seconds')
    return array

def gene_expression_mean(X):
    """
    :param X: array of matrixes. Each matrix contains gene Expression values of a cluster.
    :return: np.matrix (gene expression matrix)
    Note X is a 3D array, but we want it to be 2D. Hence we create meanMatrix
    """
    meanMatrix = X[0]
    for i in X[1:]:
        meanMatrix = np.concatenate([meanMatrix,i])
    mean =np.mean(meanMatrix,axis=0)
    return mean

def get_gene_index(gene):
    """
    :param gene: str representing the gene of which then index in adata.var['features'] is found
    :return: int representing index in adata.var['features']
    """
    for i in range(len(adata.var['features'])):
        if gene == adata.var['features'][i]:
            return i
    raise Exception('Gene not found')

def get_gene_path_df(path,geneList):
    """
    :param stage: int representing stageof cluster
    :param index: int representing index of cluster
    :param geneList: list genes of which the expression should be plotted
    :return: list,list,list,list representing the data of each gene
    """
    #add a way to differentiate for paths that end before stage 3
    stages=[]
    index=[]
    means=[]
    g=[]
    pathIndicator = []
    clusterMemory=[]
    for cluster in get_end_clusters(path):
        singularPath = get_cells_leading_to_cluster(cluster)
        for c in singularPath[::-1]:
            #if and else statement stops the code from redundant work. ie. the mean gene expr. is not calculated twice or more often, instead taken from previous paths
            if c not in clusterMemory:
                clusterMemory.append(c)
                X = load_cluster_gene_expr_matrix(*c)
                for gene in geneList:
                    s,i = c
                    pathIndicator.append(str(cluster[1])+gene)
                    stages.append(s)
                    index.append(i)
                    means.append(np.mean(X, axis=0)[get_gene_index(gene)])
                    g.append(gene)
            else:
                indicator = 0
                while(c[0] !=stages[indicator] or c[1]!=index[indicator]): indicator+=1
                for gene in geneList:
                    pathIndicator.append(str(cluster[1])+gene)
                    stages.append(stages[indicator])
                    index.append(index[indicator])
                    means.append(means[indicator])
                    g.append(gene)
                    indicator+=1 #trusts in how the lists are made, the means have to be in order in the existing list items

        X = load_cluster_gene_expr_matrix(*cluster)
        for gene in geneList:
            s, i = cluster
            pathIndicator.append(str(cluster[1])+gene)
            stages.append(s)
            index.append(i)
            means.append(np.mean(X, axis=0)[get_gene_index(gene)])
            g.append(gene)
    df = pd.DataFrame(data={'stage': stages, 'index': index, 'means': means, 'path': pathIndicator,'gene':g})
    return df

def create_genexpr_mean_df(stage,index):
    geneMatrix = load_cluster_gene_expr_matrix(stage,index)
    genes = adata.var['features'].tolist()
    mean = gene_expression_mean([geneMatrix,])
    return pd.DataFrame(list(zip(genes, mean)), columns=['Gene', 'Expression'])

def get_top_gene_expression(stage,index):
    """
    :param stage: int, representing the stage of node
    :param index: int, representing the index of node
    :return: dict with the gene expression values of the nodes top genes
    """
    stage = '0' if stage == 'c' else stage
    df = create_genexpr_mean_df(stage,index)
    #filter such that only top genes are displayed
    tGenes= adata.uns['topGene'][stage][int(index)]
    df = df[df.Gene.isin(tGenes)].sort_values(by=['Gene'])
    return df.to_dict('records')

def create_homemade_graph(G,edgeList=0):
    #This function greates the no.fig object of the graph

    for edge in edgeList:
        s1,i1 = edge[0]
        s2, i2 = edge[1]
        vertice1=int(s1) * 100 + int(i1)
        vertice2=int(s2)  * 100 + int(i2)
        G.add_edge(vertice1, vertice2)
    #creates scatter plot containing the edges
    edge_x = []
    edge_y = []
    edge_trace = []

    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        #updating the dictionary
        add_path(G.nodes[edge[0]]['cellType'],G.nodes[edge[1]]['cellType'],paths)
        cell1 = G.nodes[edge[0]]['cellType']
        cell2 = G.nodes[edge[1]]['cellType']
        text = cell1 + '--' + cell2
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)
        edge_trace.append(make_edge([x0, x1, None], [y0, y1, None],text, 3, '#000000'))

    #creates scatter plot containing the nodes
    node_x = []
    node_y = []
    node_text = []
    node_color=[]
    node_opacity=[]
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_text.append(G.nodes[node]['cellType'])
        nodePos = [G.nodes[node]['stage'],G.nodes[node]['index']]
        node_color.append(G.nodes[node]['stage'])
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            colorscale='Blackbody',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Stage',
                xanchor='left',
                titleside='right',
                dtick=1,
                tick0=0,
            ),
            line_width=2))

    node_trace.text = node_text
    node_trace.marker.color = node_color
    node_trace.marker.opacity = node_opacity

    #creates the figure object using the two scatter plots of nodes and edges from above
    fig = go.Figure(data=node_trace,
                    layout=go.Layout(
                        title='<br>Progression of cell behaviour in IPF patients',
                        titlefont_size=16,
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20,l=5,r=5,t=40),
                        annotations=[ dict(
                            showarrow=False,
                            xref="paper", yref="paper",
                            x=0.005, y=-0.002 ) ],
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        height=600,
                        width=1200,
                    )
                    )
    for trace in edge_trace:
        fig.add_trace(trace)
    return fig

def get_number_cells(stage,index):
    stage = 0 if stage == 'c' else stage
    num =len(adata.obs[(adata.obs['stage'] == int(stage)) & (adata.obs['leiden'] == int(index))].index.tolist())
    return num

def get_number_cells_stage(stage):
    totalCells = 0
    stage = [controlCells,IPF1Cells,IPF2Cells,IPF3Cells][stage]
    for c in stage:
        s,i = c.split('-')[:2]
        totalCells+=get_number_cells(s,i)
    return totalCells

def get_cells_leading_to_cluster(cluster):
    """
    :param cluster: list of length 2, indicating the clusters coordinates. it is assumed the cluster is a leaf
    :return: list of clusters leading to input cluster
    """
    #####################################################check if it works when a cluster in control is chosen that does not have children
    previous = cluster
    s, i = cluster
    s = s - 1
    clusterList = []
    for curStage in range(s, -1, -1):
        edgeTo = adata.uns['edges'][str(curStage)][previous[1]]
        previous = [curStage, edgeTo[0]]
        clusterList.append(previous)
    return clusterList

def has_child(cluster):
    s,i = cluster
    if s == 3:
        return False
    for edge in adata.uns['edges'][str(s)]:
        if edge[0]==i:
            return True
    return False

def get_end_clusters(path):
    endClusters=[]
    for cluster in path:
        if not has_child(cluster):
            endClusters.append(cluster)
    return endClusters

def get_number_cells_df(path):
    """
    :param path: list containing coordinates of clusters in selected path
    :return: pd.DataFrame,pd.DataFrame containing number of cells in each stage of path and the percentage of cells that are in the selected cluster relative to its stage. df2 is the sum of the cells in a path
    """
    stages = []
    index = []
    percentages=[]
    pathIndicator = []
    for cluster in get_end_clusters(path):
        singularPath = get_cells_leading_to_cluster(cluster)
        for c in singularPath[::-1]:
            s,i = c
            pathIndicator.append(str(cluster))
            stages.append(s)
            index.append(i)
            percentages.append(get_number_cells(s,i)/get_number_cells_stage(s))

        s,i = cluster
        pathIndicator.append(str(cluster))
        stages.append(s)
        index.append(i)
        percentages.append(get_number_cells(s,i)/get_number_cells_stage(s))
    df1 = pd.DataFrame(data={'stage': stages, 'index': index, 'percentage':percentages, 'path':pathIndicator})
    df2 = df1.drop_duplicates(['stage','index'])
    df2 = df2.groupby("stage").sum()
    df2['stage'] = [i for i in range(len(df2))]
    path=str(path[0])
    df2['path'] = [path for i in range(len(df2))]
    df2.index.name = None
    return df1,df2

colorscale=['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6',
            '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3',
            '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff', '#000000']

def getIndexesForUmap(stage,leiden):
    """
    :param stage: int gives stage
    :param leiden: int gives cluster
    :return: list of indexes. the indexes are the indexes of the cells in the cluster passed to the function
    """
    g = observations[observations.stage.eq(stage)]
    g.reset_index(inplace=True)
    return g[g.leiden.eq(leiden)].index.tolist()

def makeListOfLeidenForDF(stage):
    """
    :param stage: int representing stage
    :param leiden: int representing cluster index
    :return: list with values of cluster indices, whose index corresponds to the coordinates in the umap
    """
    stageObservations = observations[observations.stage.eq(stage)]
    stageIndexes = sorted(stageObservations.leiden.unique())
    leidenList = [-1 for i in range(len(stageObservations))]
    for leiden in stageIndexes:
        indexes = getIndexesForUmap(stage, leiden)
        for i in indexes:
            leidenList[int(i)] = str(leiden)
    if -1 in leidenList:
        raise Exception('Shouldnt be in the list anymore')
    return leidenList

def getUmapDf(stage):
    umap =adata.uns['X_umap'][str(stage)]
    df = pd.DataFrame(umap, columns=['x', 'y'])
    df['leiden'] = makeListOfLeidenForDF(stage)
    return df

def creat_umap_fig(stage):
    return px.scatter(getUmapDf(stage),x='x',y='y',color='leiden',color_discrete_sequence =colorscale)

G = create_graph()
baseFig = create_graph_fig(G)
baseFig.layout.plot_bgcolor = '#FFFFFF'

import dash
from dash.exceptions import PreventUpdate
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_table

app = dash.Dash(__name__,external_stylesheets=[dbc.themes.BOOTSTRAP])

controls = dbc.Card(
    [
        dbc.FormGroup(
            [
                dbc.Label('Selected clusters:'),
                dbc.Row(
                    [
                        dbc.Col(html.Button(' Reset ', id='reset', n_clicks=0)),
                        dbc.Col(html.Button(' Draw ', id='drawPath', n_clicks=0)),
                        dbc.Col(
                            dcc.Checklist(
                                id='hover-path',
                                options=[
                                    {'label': 'HoverPath', 'value': 'active'},
                                ],
                                value=['active'],
                                labelStyle={'display': 'inline-block'}
                            )
                        )
                    ]
                )
            ]
        ),
        dbc.FormGroup(
            [
                dbc.Label('Control clusters:'),
                dcc.Dropdown(
                    id='select-path-control',
                    options=[{'label': i, 'value': i} for i in controlCells],
                    value=[],
                    multi=True,
                ),
            ]
        ),
        dbc.FormGroup(
            [
                dbc.Label('IPF1 clusters:'),
                dcc.Dropdown(
                    id='select-path-1',
                    options=[{'label': i, 'value': i} for i in IPF1Cells],
                    value=[],
                    multi=True,
                ),
            ]
        ),
        dbc.FormGroup(
            [
                dbc.Label('IPF2 clusters:'),
                dcc.Dropdown(
                    id='select-path-2',
                    options=[{'label': i, 'value': i} for i in IPF2Cells],
                    value=[],
                    multi=True,
                ),
            ]
        ),
        dbc.FormGroup(
            [
                dbc.Label('IPF3 clusters:'),
                dcc.Dropdown(
                    id='select-path-3',
                    options=[{'label': i, 'value': i} for i in IPF3Cells],
                    value=[],
                    multi=True,
                ),
            ]
        ),
        html.Button('Confirm selected Nodes', id='submit_nodes', n_clicks=0),
        html.A(id='link', children='0-9n11-9n11-12n13n19',
               href='/assets/0-9n11-9n11-12n13n19.txt_viz/idrem_result.html', hidden=True),

    ],
    body=True,
)

geneTable=dash_table.DataTable(
    id='gene-table',
    columns=[{"name": i, "id": i} for i in ['Gene','Expression']],
    sort_action="native",
)

displayGeneTable = dbc.Modal(
    [
        dbc.ModalHeader(id='modalHeader',children="Data of"),
        dbc.ModalBody(
            geneTable,
        ),
        dbc.ModalFooter([
            dbc.Button(
                "Close",
                id="close-body-scroll",
                className="ml-auto",
                n_clicks=0,
            ),

            dbc.Button(
                "Download",
                id="download-data",
                className="ml-auto",
                n_clicks=0,
            )]
        ),
    ],
    id="modal-body-scroll",
    size="lg",
    scrollable=True,
    is_open=False,
)

choseGraphCheckbox = dcc.Checklist(
    id='ccp',
    options=[
        {'label': 'Customized path','value':'ccp'}
    ],
)

########################################################################

graphTab=dcc.Tab(label='graph-tab',children=
    dbc.Row([
        dbc.Col(dcc.Graph(id="yumin-graph",animate=True, clear_on_unhover=True,figure =baseFig),width =7),
        dbc.Col([
            dcc.Dropdown(
                id='select-umap',
                options=[{'label': i, 'value': i} for i in ['Stage 0','Stage 1','Stage 2','Stage 3']],
                value='Stage 0',
                multi=False,
            ),
            dcc.Graph(id='umap-plot')
        ],width = 5)
    ]),
)

########################################################################
geneTab=dcc.Tab(label='gene-tab',children=[
    dcc.Graph(id='gene-graph'),
])

percentageCellTab=dcc.Tab(label='percentage-tab',children=[
    dcc.Graph(id='percentage-graph',style={'height': '87vh'}),
 ])

tabs= dcc.Tabs(
    [
        graphTab,
        geneTab,
        percentageCellTab,
    ]
)

searchGene = dbc.Card(
    [
        dcc.Dropdown(
            id='search-top-gene',
            options=[{'label': i, 'value': i} for i in adata.var['features'].to_dict()],
            value=[],
            multi=True,
            placeholder="Select a gene",
        ),
        html.Button('Highlight relevant clusters', id='drawSelectedGenes', n_clicks=0),
        dcc.Checklist(
            id='onlyTopGenes',
            options=[
                {'label': 'only TopGenes','value':''}
            ],
        )
    ]
)

chooseGraph = dbc.Card(
    [
        dbc.Row([
            dbc.Col(choseGraphCheckbox),
            dbc.Col(html.Div(id='click-data',children='Your last clicked cluster:',hidden=True)),
        ]),
    ]
)

chooseClickEvent = dcc.Dropdown(
    id='select-click-event',
    clearable=False,
    options=[{'label': i, 'value': i} for i in ["iDrem", "Top_differnetial_genes"]],
    value='Top_differnetial_genes',
)

invisibleStuff = [
    dcc.Store(id='selectedCells'),
    dcc.Store(id='hoverCells'),
    dcc.Store(id='selectedGenes'),
    dcc.Store(id='selectedPath'),
    dcc.Store(id='previousFig'),
    html.Div(id="idrem",children="init",hidden = True),
    dcc.Download(id="download-topGene-expression"),
    html.Div(id="hidden-content-global2",children='test',hidden=True),
]

app.layout = dbc.Container(
    [
        html.H4('Progression of cell behaviour in IPF patients'),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(
                    [
                        controls,
                        searchGene,
                        chooseGraph,
                        chooseClickEvent
                    ],md=2),
                dbc.Col(tabs, md=10),
            ]
        ),
        displayGeneTable,
        *invisibleStuff,
    ],
    fluid=True,
)


#This callback updates the graph when to dropdown value changes
@app.callback(
    Output('yumin-graph','figure'),
    Output('click-data','hidden'),
    Output('previousFig','data'),
    Input('drawSelectedGenes', 'n_clicks'),
    Input('drawPath', 'n_clicks'),
    Input('reset', 'n_clicks'),
    Input('yumin-graph','hoverData'),
    Input('ccp','value'),
    State('selectedCells','data'),
    State('selectedGenes','data'),
    State('selectedPath','data'),
    State('hover-path','value'),
    State('previousFig','data'),
)
def show_dropdown_selection(c1,c2,c3,hoverData,customPath,dataPath,dataGene,edgeList,hoverpath,previousFig):
    #depending on which button was clicked, the top gene clusters are highlighted or the path is drawn
    ctx =dash.callback_context
    trigger = ctx.triggered[0]['prop_id'].split('.')[0]
    showLastClick = True

    if not trigger or trigger == 'reset':
        fig = baseFig
        previousFig = fig
    elif customPath:
        showLastClick = False
        if not edgeList:
            fig= baseFig
        elif edgeList and len(edgeList[-1]) < 2:
            raise PreventUpdate
        else:
            localGraph = nx.create_empty_copy(G)
            fig = create_homemade_graph(localGraph,edgeList)
    elif trigger == 'yumin-graph' and hoverpath:
        if hoverData:
            pos = [int(i) for i in hoverData['points'][0]['text'].split('-')[:2]]
            fig = create_graph_fig(G,get_list_clusters(pos))
        else:
            fig = go.Figure(previousFig)
    elif trigger == 'drawSelectedGenes':
        fig = create_graph_fig_cellGenes(G, dataGene)
        previousFig = fig
    else:
        #the list in the following are the dictionary values flattened
        fig = create_graph_fig(G, [cluster for sublist in list(dataPath.values()) for cluster in sublist] if dataPath else None)
        previousFig = fig

    fig.layout.plot_bgcolor = '#FFFFFF'

    return fig,showLastClick,previousFig

#This callback updates the Dropdown from the nodes clicked on and selected in select-path portion
@app.callback(
    Output('selectedCells','data'),
    Input('select-path-control','value'),
    Input('select-path-1','value'),
    Input('select-path-2','value'),
    Input('select-path-3','value'),
)
def display_selected_cell(valueC,value1,value2,value3):
    updatedValue = {}
    dropdownValues = [valueC,value1,value2,value3]
    #add values that were there before
    for i,valuei in enumerate(dropdownValues):
        if valuei == None:
            continue
        updatedValue[i] = []
        for cluster in valuei:
            cluster = [int(i) for i in cluster.split('-')[:2]]
            updatedValue[i].append(cluster)
    return updatedValue

#update IF1 value after control value changed
@app.callback(
    Output('select-path-1','value'),
    Input('reset','n_clicks'),
    Input('select-path-control','value'),
    State('select-path-1','value'),
)
def update1(n_clicks,value,value1):
    ctx = dash.callback_context
    trigger = ctx.triggered[0]['prop_id'].split('.')[0]
    if not ctx.triggered or trigger == 'reset':
        return []
    elif not value:
        return value1
    value1=[]
    for i in value:
        for cell in paths[i]:
            if cell in IPF1Cells and cell not in value1:
                value1.append(cell)
    return value1

#update IF2 value after IPF1 value changed
@app.callback(
    Output('select-path-2','value'),
    Input('reset','n_clicks'),
    Input('select-path-1','value'),
    State('select-path-2','value'),
)
def update2(n_clicks,value,value2):
    ctx = dash.callback_context
    trigger = ctx.triggered[0]['prop_id'].split('.')[0]
    if not ctx.triggered or trigger == 'reset':
        return []
    elif not value:
        return value2
    value2=[]
    for i in value:
        for cell in paths[i]:
            if cell in IPF2Cells and cell not in value2:
                value2.append(cell)
    return value2

#update IF3 value after IPF3 value changed
@app.callback(
    Output('select-path-3','value'),
    Input('reset','n_clicks'),
    Input('select-path-2','value'),
    State('select-path-3','value'),
)
def update3(n_clicks,value,value3):
    ctx = dash.callback_context
    trigger = ctx.triggered[0]['prop_id'].split('.')[0]
    if not ctx.triggered or trigger == 'reset':
        return []
    elif not value:
        return value3
    value3=[]
    for i in value:
        for cell in paths[i]:
            if cell in IPF3Cells and cell not in value3:
                value3.append(cell)
    return value3

#updates control value according to selected cells in following stages
@app.callback(
    Output('select-path-control','value'),
    Input('submit_nodes','n_clicks'),
    Input('reset','n_clicks'),
    State('select-path-1','value'),
    State('select-path-2', 'value'),
    State('select-path-3','value'),
)
def update_control(n_clicks1,n_clicksR,value1,value2,value3):
    ctx = dash.callback_context
    trigger = ctx.triggered[0]['prop_id'].split('.')[0]
    if not ctx.triggered or trigger == 'reset':
        return []
    valueC = []
    for i in value3:
        for cell in paths[i]:
            if cell in IPF2Cells and cell not in value2:
                value2.append(cell)

    for i in value2:
        for cell in paths[i]:
            if cell in IPF1Cells and cell not in value1:
                value1.append(cell)

    for i in value1:
        for cell in paths[i]:
            if cell in controlCells and cell not in valueC:
                valueC.append(cell)

    return valueC

@app.callback(
    Output('selectedGenes','data'),
    Input('search-top-gene','value')
)
def update_gene_Store(selectedGenes):
    correctGenes=[]
    pos = []
    for gene in selectedGenes:
        if gene in topGenes:
            for p in topGenes[gene]:
                if p not in pos:
                    pos.append(p)
    return pos

@app.callback(
    Output('search-top-gene','value'),
    Input('reset','n_clicks')
)
def empty_top_gene_dropdown(clicks):
    return []

@app.callback(
    Output('gene-graph','figure'),
    Input('search-top-gene','value'),
    State('selectedCells','data'),
)
def draw_gene_expression(geneList,data):
    if not geneList or not data['0']:
        #might not work, check later
        return go.Figure()

    df = pd.DataFrame()
    for cluster in data['0']:
        tmp = get_gene_path_df(get_list_clusters(cluster),geneList)
        df = pd.concat([df,tmp],ignore_index=True)
    fig = px.line(df, x='stage', y='means', color='path')
    fig.update_xaxes(nticks=6)
    print('fig should be there')
    return fig

@app.callback(
    Output('modalHeader','children'),
    Output('yumin-graph','clickData'),
    Output('gene-table','data'),
    Input('yumin-graph', 'clickData'),
    State('ccp','value'),
    State('select-click-event','value'),
)
def toggle_modal(clickData,valueCCP,clickEvent):
    if not clickData or valueCCP or clickEvent == 'iDrem':
        raise PreventUpdate
    table= get_top_gene_expression(*clickData['points'][0]['text'].split('-')[:2])
    return 'Gene expression of top-genes in cluster {}'.format(clickData['points'][0]['text']),None,table

@app.callback(
    Output("modal-body-scroll", "is_open"),
    Input('modalHeader','children'),
    Input('close-body-scroll','n_clicks'),
    State('modal-body-scroll','is_open'),
    State('ccp','value'),
)
def show_gene_expr_table(header,n1,is_open,value):
    ctx = dash.callback_context
    if ctx.triggered and (n1 or header) and not value:
        return not is_open
    return is_open

@app.callback(
    Output('link','children'),
    Output('link','href'),
    Output('link','hidden'),
    Input('selectedCells','data'),
)
def update_link(data):
    ctx =dash.callback_context
    hidden = True
    if data['0'] == []:
        return 'nothing to see','/assets/0-9n11-9n11-12n13n19.txt_viz/idrem_result.html', hidden
    href = '/assets/'
    for i in idremPaths:
        cluster = int(idremPaths[i][0][0])
        if data['0'][0][1] == cluster:
            href +=i
            hidden = False
            break
    href += '/idrem_result.html'
    return f'Show path from root {cluster}',href,hidden

@app.callback(
    Output('selectedPath','data'),
    Input("yumin-graph",'clickData'),
    Input('reset', 'n_clicks'),
    State('selectedPath','data'),
    State('ccp','value')
)
def update_path_data(clickData,reset,data,value):
    ctx = dash.callback_context
    trigger = ctx.triggered[0]['prop_id'].split('.')[0]
    if not value or trigger=='reset':
        return []
    stage,index = clickData['points'][0]['text'].split('-')[:2]
    if data and len(data[-1]) <2:
        data[-1].append([stage,index])
    else:
        data.append([[stage,index]])
    return data

@app.callback(
    Output('click-data','children'),
    Input("yumin-graph",'clickData'),
)
def update_click_info(clickData):
    if not clickData:
        raise PreventUpdate
    cluster = clickData['points'][0]['text']
    return f'Your last clicked cluster: {cluster}'

@app.callback(
    Output('hover-path','value'),
    Input('drawPath', 'n_clicks'),
    Input('reset', 'n_clicks'),
)
def change_hoverPath(c1,c2):
    ctx = dash.callback_context
    trigger = ctx.triggered[0]['prop_id'].split('.')[0]
    if trigger == 'reset' or not ctx.triggered:
        return ['active']
    else:
        return []

@app.callback(
    Output('percentage-graph','figure'),
    Input('submit_nodes','n_clicks'),
    State('selectedCells','data')
)
def make_percentage_graph(c1,data):
    """
    edge case that data is empty not yet considered!!! Also only works if merely one path is chosen
    :param c1: int representing number of clicks( not used)
    :param controlCluster: list representing coordinates of selected cell(s) in control stage
    :return: go.figure   showing the percentage of the cells in each stage in relation to the whole path
    """
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    df1 = pd.DataFrame()
    df2 = pd.DataFrame()
    for cluster in data['0']:
        tmp1,tmp2 = get_number_cells_df(get_list_clusters(cluster))
        df1 = pd.concat([df1,tmp1],ignore_index=True)
        df2 = pd.concat([df2,tmp2])
    fig = make_subplots(rows=2, cols=1,
                    shared_xaxes=True,)
    for path in df1['path'].unique():
        tmp = df1[df1['path'] == path]
        fig.add_trace(
            go.Scatter(x=tmp['stage'], y=tmp['percentage'],name=path),
            row=1, col=1
        )

    for path in df2['path'].unique():
        tmp = df2[df2['path'] == path]
        fig.add_trace(
            go.Scatter(x=tmp['stage'], y=tmp['percentage'], name='Path ' + path),
            row=2, col=1,
        )

    fig.update_layout(
        xaxis = dict(
        tickmode = 'array',
        tickvals = [0, 1, 2, 3],
        ticktext = [0, 1, 2, 3],
    ))

    return fig

@app.callback(
    Output("idrem","children"),
    Input('yumin-graph','clickData'),
    State('select-click-event','value'),
)
def update_link(clickData,clickEvent):
    ctx =dash.callback_context
    if not ctx.triggered or clickEvent!='iDrem':
        raise PreventUpdate
    stage,index = clickData['points'][0]['text'].split('-')[:2]
    href = "/assets/"
    for i in idremPaths:
        clusters = idremPaths[i][int(stage)]
        if index in clusters:
            href +=i
            break
    href += "/idrem_result.html"
    return href

app.clientside_callback(
"""
function(href,value) {
if (href !="init" && value == "iDrem") {
window.open(href)
} 
return href
}
""",
Output("hidden-content-global2", "children"),
Input("idrem", "children"),
State("select-click-event","value")
)

@app.callback(
    Output('download-topGene-expression','data'),
    Input('download-data','n_clicks'),
    [State('modalHeader','children'),
    State('gene-table','data')],
    prevent_initial_call=True,
)
def download_data(n_clicks,header,df):
    df = pd.DataFrame(df)
    name = header.split(' ')[-1]
    return dcc.send_data_frame(df.to_csv,name+'.csv')

@app.callback(
    Output('search-top-gene','options'),
    Input('onlyTopGenes','value'),
    prevent_initial_call=True,
)
def change_geneselection_options(value):
    if value:
        return [{'label': i, 'value': i} for i in sorted(topGenes)]
    else:
        return [{'label': i, 'value': i} for i in adata.var['features'].to_dict()]

@app.callback(
    Output('umap-plot','figure'),
    Input('select-umap','value')
)
def change_umap(stage):
    if stage:
        stage= stage.split(' ')[1]
    else:
        print('is enterd')
        stage = 0
    fig=creat_umap_fig(int(stage))
    fig.update_yaxes(visible=False, showticklabels=False)

    return fig


#app.run_server(host="0.0.0.0",debug=True,port=3003)


"""
Dash requires the knowledge of the path used to access the app.
ShinyProxy makes this path available as an environment variable,
which we expose to dash below.
"""
if 'SHINYPROXY_PUBLIC_PATH' in os.environ:
    app.config.update({
        'routes_pathname_prefix': os.environ['SHINYPROXY_PUBLIC_PATH'],
        'requests_pathname_prefix': os.environ['SHINYPROXY_PUBLIC_PATH']
    })

if __name__ == "__main__":
    dev = False  # Set to True if in development
    app.run_server(debug=dev, use_reloader=dev, host='0.0.0.0',port=3003)


#http://127.0.0.1:3003/
#local link
