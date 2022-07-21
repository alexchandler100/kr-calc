import os
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
pd.set_option('display.max_rows', 1000)
os.chdir('/Users/alexchandler/kr-calc/csv_files/') #change to csv file directory
knot_info = pd.read_csv('Knotinfo_data.csv', dtype='object')

import sklearn
import scipy
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split, KFold, cross_val_score
from sklearn import preprocessing
from sklearn.preprocessing import PolynomialFeatures
from sklearn.feature_selection import VarianceThreshold
import itertools


# import required module for iterating through directory
import os
import json
# assign directory
os.chdir('/Users/alexchandler/kr-calc') #change to csv file directory
directory = 'data'
homfly_data={}
#we create a dictionary with keys: names of knots
                          #values: dictionaries (keys: nonzero gradings, values: homology dimensions)homfly_data={}
# sort directory and then iterate over files in that directory, creating a dictionary for easier use
listing = os.listdir(directory)
listing.sort()
for filename in listing:
    f = os.path.join(directory, filename)
    # checking if it is a file
    if os.path.isfile(f):
        with open(f) as json_file:
            data = json.load(json_file)
            data_dict = {}
            for i in range(int(len(data)/2)):
                data_dict[tuple(data[2*i])]=data[2*i+1]
            homfly_data[filename[:-5]]=data_dict

#defining some functions
def convert(tup):
    return (tup[0],tup[1],int((tup[2]-tup[1])/2))

def convert_homology(dictionary):
    result={}
    for key in list(dictionary.keys()):
        result[convert(key)]=dictionary[key]
    return result

def convert_NS_DGR(tup):
    return (tup[0],tup[1],int((tup[1]-tup[2])/2))

def convert_homology_NS_DGR(dictionary):
    result={}
    for key in list(dictionary.keys()):
        result[convert_NS_DGR(key)]=dictionary[key]
    return result

#still not sure if this is computing the right thing...
def check_parity(knot):
    homology = homfly_data[knot]
    dictionary = convert_homology(homology)
    obstructions = 0
    agradings = []
    for key in list(dictionary.keys()):
        agradings.append(key[1])
    agradings.sort()
    tgradings = []
    for key in list(dictionary.keys()):
        tgradings.append(key[2])
    tgradings.sort()
    tgradings = list(set(tgradings))
    a_and_t_dict={}
    for a in agradings:
        parities = []
        a_and_t_dict[a]=[key[2] for key in list(dictionary.keys()) if key[1]==a]
        for i in range(len(a_and_t_dict[a])):
            for j in range(i+1,len(a_and_t_dict[a])):
                parities.append((a_and_t_dict[a][j]-a_and_t_dict[a][i])%2)
        if len(list(set(parities)))<=1:
            obstructions+=0
        else:
            obstructions+=1
    #print(a_and_t_dict)
    if obstructions==0:
        return 1.0
    else:
        return 0.0


#determines if there are any nonzero homology groups in any odd gradings
#input: homology as dictionary, e.g. check_odd(homfly_data['10_1'])
#wait... what are we supposed to be checking here?
def check_odd(dictionary):
    gradings = (0,)
    for key in list(dictionary.keys()):
        gradings+=key
    for grad in gradings:
        if grad%2==1:
            #print(dictionary)
            return True
            break;
    return False


#returns a vector from point1 to point2
def vector(point1,point2):
    return (point2[0]-point1[0],point2[1]-point1[1],point2[2]-point1[2]);

def most_frequent(List):
    return max(set(List), key = List.count)

def normalize(vect):
    if vect[0]!=0:
        return (vect[0]/vect[0],vect[1]/vect[0], vect[2]/vect[0])
    else:
        return vect

#wrote this function before I knew the planes were always (1,1,1)... better now to use planar_support_alt
#takes as input the homfly homology and outputs the number of planes the homology is supported on
def planar_support(homology):
    planes = 0
    points = list(homology.keys())
    while len(points)>0:
        base_point = points[0]
        #remove the base point
        points=points[1:]
        #def: given a point, its "maximal plane" is the plane containing it which contains the maximal number of other points possible
        #lets find the maximal plane containing the base point
        normals=[]
        for i in range(len(points)):
            for j in range(i+1,len(points)):
                vect1 = vector(base_point,points[i])
                vect2 = vector(base_point,points[j])
                cross_prod = tuple(np.cross(vect1,vect2))
                if (cross_prod!=(0,0,0)):
                    normals.append(normalize(tuple(cross_prod)))
        if (len(points)>1):
            base_normal = most_frequent(normals)
            #get rid of all points on the maximal plane of the base point
            points = [point for point in points if np.dot(base_normal,(point[0]-base_point[0],point[1]-base_point[1],point[2]-base_point[2]))!=0]
            planes+=1
        else:
            planes+=1
            points=[]
    print('normal vector: '+str(base_normal))
    return planes

#takes as input the homfly homology and outputs the number of planes the homology is supported on
#knowing the planes have to have normal (1,1,1) makes it a lot easier to compute this
def planar_support_alt(homology):
    planes = 0
    points = list(homology.keys())
    while len(points)>0:
        base_point=points[0]
        delta = sum(list(base_point))
        points=[point for point in points if sum(list(point))!=delta]
        planes+=1
    return planes


#0 just indicates that the knot is not in the homfly_data dataset
def compute_planar_support(knot):
    if knot in homfly_data:
        return planar_support_alt(homfly_data[knot])
    else:
        return 0

#this makes an interactive 3d plot for homology
import plotly
import plotly.graph_objs as go
import random

def plot_homology(knot):
    homology = homfly_data[knot]
    gradings = homology.keys()
    # Configure Plotly to be rendered inline in the notebook.
    plotly.offline.init_notebook_mode()
    #making colors
    colors=[]
    color_dict={}
    for i in range(max(list(homology.values()))):
        color_dict[i+1]=random.choices(range(256), k=3)
    for grading in gradings:
        colors.append(color_dict[homology[grading]])
        #if homology[grading]>2:
            #colors.append('red')
        #elif homology[grading]==2:
            #colors.append('green')
        #else:
            #colors.append('black')
    #making sizes
    sizes=[]
    for grading in gradings:
        sizes.append(3+3*homology[grading])
    # Configure the trace.
    trace = go.Scatter3d(
        x=[grading[0] for grading in gradings],
        y=[grading[1] for grading in gradings],
        z=[grading[2] for grading in gradings],
        mode='markers',
        marker={'size': sizes,'opacity': 0.8,'color': colors})
    deltas=list(set([grading[0]+grading[1]+grading[2] for grading in gradings]))
    planes=[]
    for delta in deltas:
        xs=list(set([grading[0] for grading in gradings if grading[0]+grading[1]+grading[2]==delta]))
        ys=list(set([grading[1] for grading in gradings if grading[0]+grading[1]+grading[2]==delta]))
        x = np.outer(np.linspace(min(xs), max(xs), 30), np.ones(30))
        y = np.outer(np.linspace(min(ys), max(ys), 30), np.ones(30)).copy().T
        z = delta-x-y
        planes.append(go.Surface(x=x, y=y, z=z, showscale=False, opacity=0.5))
    # Configure the layout.
    layout = go.Layout(margin={'l': 0, 'r': 0, 'b': 0, 't': 0})
    data = [trace]
    for plane in planes:
        data.append(plane)
    x_lines = list()
    y_lines = list()
    z_lines = list()

    if len(deltas)==2:
        deltas.sort()
        delta1=deltas[0]
        delta2=deltas[1]
        d1gradings=[grading for grading in gradings if grading[0]+grading[1]+grading[2]==delta1]
        d2gradings=[grading for grading in gradings if grading[0]+grading[1]+grading[2]==delta2]
        pairs=[]
        for grading1 in d1gradings:
            for grading2 in d2gradings:
                if grading2[0]==grading1[0]+4 and grading2[1]==grading1[1]-2 and grading2[2]==grading1[2]:
                    x_lines.append(grading1[0])
                    x_lines.append(grading2[0])
                    x_lines.append(None)
                    y_lines.append(grading1[1])
                    y_lines.append(grading2[1])
                    y_lines.append(None)
                    z_lines.append(grading1[2])
                    z_lines.append(grading2[2])
                    z_lines.append(None)
        trace_lines = go.Scatter3d(
            x=x_lines,
            y=y_lines,
            z=z_lines,
            line = dict(width = 2, color = 'rgb(200,0,100)'))
        data.append(trace_lines)

    annot_list=[]
    for point in homology.keys():
        annot_list.append(dict(
            showarrow=False,
            x=point[0],
            y=point[1],
            z=point[2],
            text=str(homology[point]),
            xanchor="left",
            xshift=10,
            opacity=0.7))
    plot_figure = go.Figure(data=data, layout=layout)
    plot_figure.update_layout(
        scene=dict(
            annotations=annot_list
        ),
        )
    # Render the plot.
    print('Reduced HOMFLY-PT Homology of '+knot)
    plotly.offline.iplot(plot_figure)


import matplotlib.pyplot as plt
import matplotlib.image as mpimg




compute_planar_support_vect = np.vectorize(compute_planar_support)


# n = dimension of full chain complex (must have n odd for the answer to be nonempty)
# outputs all possible chain complexes (just the list of dimensions) which are acyclic
def reduce_once(set_of_complexes):
    possibilities=[]
    for complx in set_of_complexes:
        possible_indices = [i for i in range(len(complx)-1) if complx[i]>0 and complx[i+1]>0]
        for i in possible_indices:
            temp_complx=[]
            for j in range(len(complx)):
                if j!=i and j!=i+1:
                    temp_complx.append(complx[j])
                else:
                    temp_complx.append(complx[j]-1)
            possibilities.append(tuple(temp_complx))
    return possibilities

import math

def all_reductions(cmplx):
    n=sum(cmplx)
    k=math.floor(n/2)
    cc=[cmplx]
    for i in range(k):
        cc=reduce_once(cc)
    return list(set(cc))

def possible_homology_indices(cmplx):
    if len(cmplx)==1:
        return [0]
    elif len(cmplx)==0:
        return []
    else:
        reductions=all_reductions(cmplx)
        indices=[]
        for cmplx in reductions:
            try:
                i=[j for j in range(len(cmplx)) if cmplx[j]!=0][0]
                indices.append(i)
            except:
                continue
        return indices


#these just give the possible keys
def lee_grading_d1(knot):
    if (compute_planar_support(knot)==1):
        homology = homfly_data[knot]
        gradings = homology.keys()

        #computing chains for d1 differential
        d1_heights=list(set([g[0]+g[1] for g in gradings]))
        d1_heights.sort()
        d1_chains = {}
        d1_lee_height=0
        for i in d1_heights:
            d1_chains[i]={}
            gradings_at_height_i=[g for g in gradings if g[0]+g[1]==i]
            gradings_at_height_i.sort()
            for g in gradings_at_height_i:
                d1_chains[i][g]=homfly_data[knot][g]
        #for key in d1_chains.keys():
            #print(key,d1_chains[key])
        #computing euler characteristics for d1 differential
        for i in d1_heights:
            euler_char=0
            k=0
            for g in d1_chains[i].keys():
                euler_char+=((-1)**k)*d1_chains[i][g]
                k+=1
            if euler_char**2==1:
                #print(i)
                d1_lee_height=i

        #computing x,y location of unique position of homology with respect to d1,dm1
        chains = []
        for key in d1_chains[d1_lee_height]:
            chains.append(d1_chains[d1_lee_height][key])
        #print(chains)     #prints chain complex before reduction
        #reduce chain complex
        possibilities = possible_homology_indices(chains)
        keys=[]
        for i in possibilities:
            keys.append(list(d1_chains[d1_lee_height].keys())[i])
        return keys
    elif (compute_planar_support(knot)==0):
        print(str(knot)+' is not planar')
    else:
        print('homfly homology is not known for '+str(knot))

def lee_grading_dm1(knot):
    if (compute_planar_support(knot)==1):
        homology = homfly_data[knot]
        gradings = homology.keys()

        #computing chains for dm1 differential
        dm1_heights=list(set([g[0]-g[1] for g in gradings]))
        dm1_heights.sort()
        dm1_chains = {}
        dm1_lee_height=0
        for i in dm1_heights:
            dm1_chains[i]={}
            gradings_at_height_i=[g for g in gradings if g[0]-g[1]==i]
            gradings_at_height_i.sort()
            for g in gradings_at_height_i:
                dm1_chains[i][g]=homfly_data[knot][g]
        #for key in dm1_chains.keys():
            #print(key,dm1_chains[key])
        #computing euler characteristics for dm1 differential
        for i in dm1_heights:
            euler_char=0
            k=0
            for g in dm1_chains[i].keys():
                euler_char+=((-1)**k)*dm1_chains[i][g]
                k+=1
            if euler_char**2==1:
                dm1_lee_height=i
                #print(i)
        #computing x,y location of unique position of homology with respect to d1,dm1 (no... there are two locations... this is wrong)
        chains = []
        for key in dm1_chains[dm1_lee_height]:
            chains.append(dm1_chains[dm1_lee_height][key])
        #print(chains)
        #reduce chain complex
        possibilities = possible_homology_indices(chains)
        keys=[]
        for i in possibilities:
            keys.append(list(dm1_chains[dm1_lee_height].keys())[i])
        return keys
    elif (compute_planar_support(knot)==0):
        print(str(knot)+' is not planar')
    else:
        print('homfly homology is not known for '+str(knot))



from matplotlib.patches import Ellipse
plt.rcParams.update({'figure.max_open_warning': 0})


#for planar knots we drop the z coordinate and plot just the xy and the value of the homology
def plot_homology_planar(knot,convention='NS',gridlines=True,s_inv_circles_on=True,vars='qa'):
    if (compute_planar_support(knot)==1):
        homology = homfly_data[knot]
        if convention=='DGR':
            homology=convert_homology_NS_DGR(homology)
        gradings = homology.keys()
        first = list(gradings)[0]
        if convention=='NS':
            delta = str(first[0]+first[1]+first[2])
            deltastr = 'q + a + t'
        elif convention=='DGR':
            delta = str(int(-first[0]/2-first[1]+first[2]))
            deltastr = '-q/2 - a + t_DGR'
        #making markers as ranks of homology groups
        markers=[]
        for grading in gradings:
            markers.append(homology[grading])
        xs = [grading[0] for grading in gradings]
        ys = [grading[1] for grading in gradings]
        #configuring the plot
        fig, ax = plt.subplots(1, 2, figsize=(10,10))
        #ax[0].scatter(xs, ys, s=1, c="white")
        # change default range
        ax[0].set_xlim((min(xs)-1, max(xs)+1))
        ax[0].set_ylim((min(ys)-1, max(ys)+1))
        # Minor ticks
        ax[0].set_xticks(np.arange(min(xs)-.5, max(xs)+1, 1), minor=True)
        ax[0].set_yticks(np.arange(min(ys)-.5, max(ys)+1, 1), minor=True)
        # Gridlines based on minor ticks
        if gridlines:
            ax[0].grid(which='minor', color='gray', linestyle='dotted', linewidth=1)
        #computing the signature
        base = list(gradings)[0]
        signature = str(-(base[0]+base[1]+base[2]))
        #setting the title
        ax[0].set(title='Reduced HOMFLY-PT Homology of '+knot+'\nDelta = '+delta+' = '+deltastr,
                  aspect=1, xticks=range(min(xs)-1, max(xs)+2), yticks=range(min(ys)-1, max(ys)+2))
        ax[0].set_xlabel('q')
        ax[0].set_ylabel('a')
        #setting xscale and yscale to be equal
        plt.gca().set_aspect('equal', adjustable='box')
        #circling gradings containing d1 and dm1 homology
        if s_inv_circles_on and vars=='qa' and convention == 'NS':
            for grading in lee_grading_d1(knot):
                ellipse = Ellipse(xy=(grading[0], grading[1]),
                               width=1.4, height=1.4, edgecolor='r', fc='None', lw=2)
                ax[0].add_patch(ellipse)
            for grading in lee_grading_dm1(knot):
                ellipse = Ellipse(xy=(grading[0], grading[1]),
                               width=1.4, height=1.4, edgecolor='r', fc='None', lw=2)
                ax[0].add_patch(ellipse)
        elif s_inv_circles_on:
            print("d_1 and d_-1 homology only implemented for convention=='NS' and vars='qa'")
        #making colors depending on size of homology
        colors=[]
        for txt in markers:
            if txt%6==1:
                colors.append("black")
            elif txt%6==2:
                colors.append("blue")
            elif txt%6==3:
                colors.append("green")
            elif txt%6==4:
                colors.append("red")
            elif txt%6==5:
                colors.append("purple")
            elif txt%6==0:
                colors.append("orange")
        #drawing the ranks of homology groups over top of invisible dots
        for i, txt in enumerate(markers):
            text_kwargs = dict(ha='center', va='center', fontsize=20, color=colors[i])
            ax[0].text(xs[i], ys[i], txt, **text_kwargs)
            #ax[0].annotate(txt, xy=(xs[i], ys[i]),xytext=(xs[i]-0.3, ys[i]-0.4),fontsize=20, c=colors[i])
        #showing image of knot
        ax[1].set(title=knot, aspect=1, xticks=[], yticks=[])
        os.chdir('/Users/alexchandler/kr-calc') #change to python script directory
        img = mpimg.imread('diagrams/'+knot+'.png')
        imgplot = ax[1].imshow(img)
        #plt.show()
    elif (compute_planar_support(knot)==0 or compute_planar_support(knot)>=2):
        print(str(knot)+' is not planar')
    else:
        print('homfly homology is not known for '+str(knot))


    #plt.plot(xs, ys, markers)
    #for x, y in zip(xs, ys):
        #plt.text(x, y, markers, color="black", fontsize=12)





def compute_deltas(knot):
    if knot in homfly_data:
        homology = homfly_data[knot]
        return planar_support_alt(homology)
    else:
        return 0

compute_planar_support_vect = np.vectorize(compute_planar_support)

#computes list of all delta gradings
def delta(knot,convention='NS'):
    if knot in homfly_data:
        homology = homfly_data[knot]
        deltas = []
        if convention=='NS':
            for g in list(homology.keys()):
                deltas.append(g[0]+g[1]+g[2])
            return set(deltas)
        elif convention=='DGR':
            homology=convert_homology_NS_DGR(homology)
            for g in list(homology.keys()):
                deltas.append(int(-g[0]/2-g[1]+g[2]))
            return set(deltas)
    else:
        return 0

delta_vect = np.vectorize(delta)


#these just give the possible keys
def lee_grading_d1_nonplanar(knot):
    if (compute_planar_support(knot)==2):
        homology = homfly_data[knot]
        gradings = homology.keys()
        deltas=list(delta(knot))
        deltas.sort()
        gradings1=[grad for grad in gradings if grad[0]+grad[1]+grad[2]==deltas[0]]
        gradings2=[grad for grad in gradings if grad[0]+grad[1]+grad[2]==deltas[1]]
        #compute possible gradings for first plane
        d1_heights1=list(set([g[0]+g[1] for g in gradings1]))
        d1_heights1.sort()
        d1_chains1 = {}
        d1_lee_height1=0
        for i in d1_heights1:
            d1_chains1[i]={}
            gradings_at_height_i1=[g for g in gradings1 if g[0]+g[1]==i]
            gradings_at_height_i1.sort()
            for g in gradings_at_height_i1:
                d1_chains1[i][g]=homfly_data[knot][g]
        for i in d1_heights1:
            euler_char=0
            k=0
            for g in d1_chains1[i].keys():
                euler_char+=((-1)**k)*d1_chains1[i][g]
                k+=1
            if euler_char**2==1:
                d1_lee_height1=i
        chains1 = []
        try:
            for key in d1_chains1[d1_lee_height1]:
                chains1.append(d1_chains1[d1_lee_height1][key])
        except:
            pass
        possibilities1 = possible_homology_indices(chains1)
        keys1=[]
        for i in possibilities1:
            keys1.append(list(d1_chains1[d1_lee_height1].keys())[i])
        #compute possible gradings for second plane
        d1_heights2=list(set([g[0]+g[1] for g in gradings2]))
        d1_heights2.sort()
        d1_chains2 = {}
        d1_lee_height2=0
        for i in d1_heights2:
            d1_chains2[i]={}
            gradings_at_height_i2=[g for g in gradings2 if g[0]+g[1]==i]
            gradings_at_height_i2.sort()
            for g in gradings_at_height_i2:
                d1_chains2[i][g]=homfly_data[knot][g]
        for i in d1_heights2:
            euler_char=0
            k=0
            for g in d1_chains2[i].keys():
                euler_char+=((-1)**k)*d1_chains2[i][g]
                k+=1
            if euler_char**2==1:
                d1_lee_height2=i
        chains2 = []
        try:
            for key in d1_chains2[d1_lee_height2]:
                chains2.append(d1_chains2[d1_lee_height2][key])
        except:
            pass
        possibilities2 = possible_homology_indices(chains2)
        keys2=[]
        for i in possibilities2:
            keys2.append(list(d1_chains2[d1_lee_height2].keys())[i])
        return [keys1,keys2]
    elif (compute_planar_support(knot)!=2):
        print(str(knot)+' does not have planar support 2')
    else:
        print('homfly homology is not known for '+str(knot))

def lee_grading_dm1_nonplanar(knot):
    if (compute_planar_support(knot)==2):
        homology = homfly_data[knot]
        gradings = homology.keys()
        deltas=list(delta(knot))
        deltas.sort()
        gradings1=[grad for grad in gradings if grad[0]+grad[1]+grad[2]==deltas[0]]
        gradings2=[grad for grad in gradings if grad[0]+grad[1]+grad[2]==deltas[1]]

        #computing gradings for first plane
        dm1_heights1=list(set([g[0]-g[1] for g in gradings1]))
        dm1_heights1.sort()
        dm1_chains1 = {}
        dm1_lee_height1=0
        for i in dm1_heights1:
            dm1_chains1[i]={}
            gradings_at_height_i1=[g for g in gradings1 if g[0]-g[1]==i]
            gradings_at_height_i1.sort()
            for g in gradings_at_height_i1:
                dm1_chains1[i][g]=homfly_data[knot][g]
        for i in dm1_heights1:
            euler_char=0
            k=0
            for g in dm1_chains1[i].keys():
                euler_char+=((-1)**k)*dm1_chains1[i][g]
                k+=1
            if euler_char**2==1:
                dm1_lee_height1=i
        chains1 = []
        try:
            for key in dm1_chains1[dm1_lee_height1]:
                chains1.append(dm1_chains1[dm1_lee_height1][key])
        except:
            pass
        possibilities1 = possible_homology_indices(chains1)
        keys1=[]
        for i in possibilities1:
            keys1.append(list(dm1_chains1[dm1_lee_height1].keys())[i])
        #computing gradings for second plane
        dm1_heights2=list(set([g[0]-g[1] for g in gradings2]))
        dm1_heights2.sort()
        dm1_chains2 = {}
        dm1_lee_height2=0
        for i in dm1_heights2:
            dm1_chains2[i]={}
            gradings_at_height_i2=[g for g in gradings2 if g[0]-g[1]==i]
            gradings_at_height_i2.sort()
            for g in gradings_at_height_i2:
                dm1_chains2[i][g]=homfly_data[knot][g]
        for i in dm1_heights2:
            euler_char=0
            k=0
            for g in dm1_chains2[i].keys():
                euler_char+=((-1)**k)*dm1_chains2[i][g]
                k+=1
            if euler_char**2==1:
                dm1_lee_height2=i
        chains2 = []
        try:
            for key in dm1_chains2[dm1_lee_height2]:
                chains2.append(dm1_chains2[dm1_lee_height2][key])
        except:
            pass
        possibilities2 = possible_homology_indices(chains2)
        keys2=[]
        for i in possibilities2:
            keys2.append(list(dm1_chains2[dm1_lee_height2].keys())[i])
        return [keys1,keys2]
    elif (compute_planar_support(knot)!=2):
        print(str(knot)+' does not have planar support 2')
    else:
        print('homfly homology is not known for '+str(knot))

def plot_homology_nonplanar(knot,convention='NS',gridlines=True,s_inv_circles_on=True, vars='qa'):
    if (compute_planar_support(knot)==2):
        homology = homfly_data[knot]
        if convention=='DGR':
            homology=convert_homology_NS_DGR(homology)
        gradings = homology.keys()
        deltas=list(delta(knot,convention))
        deltas.sort()
        if convention=='NS':
            gradings1=[grad for grad in gradings if grad[0]+grad[1]+grad[2]==deltas[0]]
            gradings2=[grad for grad in gradings if grad[0]+grad[1]+grad[2]==deltas[1]]
            deltastr = 'q + a + t'
        elif convention=='DGR':
            gradings1=[grad for grad in gradings if int(-grad[0]/2-grad[1]+grad[2])==deltas[0]]
            gradings2=[grad for grad in gradings if int(-grad[0]/2-grad[1]+grad[2])==deltas[1]]
            deltastr = '-q/2 - a + t_DGR'
        fig, ax = plt.subplots(1, 3, figsize=(10,10))
        #creating plot for first delta grading
        markers1=[]
        for grading in gradings1:
            markers1.append(homology[grading])
        xs1 = [grading[0] for grading in gradings1]
        ys1 = [grading[1] for grading in gradings1]
        if vars =='qt':
            xs1 = [grading[0] for grading in gradings1]
            ys1 = [grading[2] for grading in gradings1]
        elif vars =='at':
            xs1 = [grading[1] for grading in gradings1]
            ys1 = [grading[2] for grading in gradings1]
        ax[0].set_xlim((min(xs1)-1, max(xs1)+1))
        ax[0].set_ylim((min(ys1)-1, max(ys1)+1))
        # Minor ticks
        ax[0].set_xticks(np.arange(min(xs1)-.5, max(xs1)+1, 1), minor=True)
        ax[0].set_yticks(np.arange(min(ys1)-.5, max(ys1)+1, 1), minor=True)
        ax[0].set_xlabel('q')
        ax[0].set_ylabel('a')
        # Gridlines based on minor ticks
        if gridlines:
            ax[0].grid(which='minor', color='gray', linestyle='dotted', linewidth=1)
        base1 = list(gradings1)[0]
        signature1 = str((base1[0]+base1[1]+base1[2]))
        ax[0].set(title='Reduced HOMFLY-PT Homology of '+knot+'\nDelta = '+str(deltas[0])+' = '+deltastr,
                  aspect=1, xticks=range(min(xs1)-1, max(xs1)+2), yticks=range(min(ys1)-1, max(ys1)+2))
        plt.gca().set_aspect('equal', adjustable='box')
        if s_inv_circles_on and vars=='qa' and convention == 'NS':
            for grading in lee_grading_d1_nonplanar(knot)[0]:
                ellipse = Ellipse(xy=(grading[0], grading[1]),
                               width=1.4, height=1.4, edgecolor='r', fc='None', lw=2)
                ax[0].add_patch(ellipse)
            for grading in lee_grading_dm1_nonplanar(knot)[0]:
                ellipse = Ellipse(xy=(grading[0], grading[1]),
                               width=1.4, height=1.4, edgecolor='r', fc='None', lw=2)
                ax[0].add_patch(ellipse)
        colors1=[]
        for txt in markers1:
            if txt%6==1:
                colors1.append("black")
            elif txt%6==2:
                colors1.append("blue")
            elif txt%6==3:
                colors1.append("green")
            elif txt%6==4:
                colors1.append("red")
            elif txt%6==5:
                colors1.append("purple")
            elif txt%6==0:
                colors1.append("orange")
        for i, txt in enumerate(markers1):
            text_kwargs = dict(ha='center', va='center', fontsize=20, color=colors1[i])
            ax[0].text(xs1[i], ys1[i], txt, **text_kwargs)

        #creating plot for second delta grading
        markers2=[]
        for grading in gradings2:
            markers2.append(homology[grading])
        xs2 = [grading[0] for grading in gradings2]
        ys2 = [grading[1] for grading in gradings2]
        if vars =='qt':
            xs2 = [grading[0] for grading in gradings2]
            ys2 = [grading[2] for grading in gradings2]
        elif vars =='at':
            xs2 = [grading[1] for grading in gradings2]
            ys2 = [grading[2] for grading in gradings2]
        ax[1].set_xlim((min(xs2)-1, max(xs2)+1))
        ax[1].set_ylim((min(ys2)-1, max(ys2)+1))
        # Minor ticks
        ax[1].set_xticks(np.arange(min(xs2)-.5, max(xs2)+1, 1), minor=True)
        ax[1].set_yticks(np.arange(min(ys2)-.5, max(ys2)+1, 1), minor=True)
        ax[1].set_xlabel('q')
        ax[1].set_ylabel('a')
        # Gridlines based on minor ticks
        if gridlines:
            ax[1].grid(which='minor', color='gray', linestyle='dotted', linewidth=1)
        base2 = list(gradings2)[0]
        signature2 = str((base2[0]+base2[1]+base2[2]))
        ax[1].set(title='Reduced HOMFLY-PT Homology of '+knot+'\nDelta = '+str(deltas[1])+' = '+deltastr,
                  aspect=1, xticks=range(min(xs2)-1, max(xs2)+2), yticks=range(min(ys2)-1, max(ys2)+2))
        plt.gca().set_aspect('equal', adjustable='box')
        if s_inv_circles_on and vars=='qa' and convention == 'NS':
            for grading in lee_grading_d1_nonplanar(knot)[1]:
                ellipse = Ellipse(xy=(grading[0], grading[1]),
                               width=1.4, height=1.4, edgecolor='r', fc='None', lw=2)
                ax[1].add_patch(ellipse)
            for grading in lee_grading_dm1_nonplanar(knot)[1]:
                ellipse = Ellipse(xy=(grading[0], grading[1]),
                               width=1.4, height=1.4, edgecolor='r', fc='None', lw=2)
                ax[1].add_patch(ellipse)
        elif s_inv_circles_on:
            print('d_1 and d_-1 homology only implemented for vars="qa" and convention="NS"')
        colors2=[]
        for txt in markers2:
            if txt%6==1:
                colors2.append("black")
            elif txt%6==2:
                colors2.append("blue")
            elif txt%6==3:
                colors2.append("green")
            elif txt%6==4:
                colors2.append("red")
            elif txt%6==5:
                colors2.append("purple")
            elif txt%6==0:
                colors2.append("orange")
        for i, txt in enumerate(markers2):
            text_kwargs = dict(ha='center', va='center', fontsize=20, color=colors2[i])
            ax[1].text(xs2[i], ys2[i], txt, **text_kwargs)
        ax[2].set(title=knot, aspect=1, xticks=[], yticks=[])
        os.chdir('/Users/alexchandler/kr-calc') #change to python script directory
        img = mpimg.imread('diagrams/'+knot+'.png')
        imgplot = ax[2].imshow(img)
    elif (compute_planar_support(knot)==0 or compute_planar_support(knot)>=3):
        print(str(knot)+' is not on two planes... this is not yet implemented')
    else:
        print('homfly homology is not known for '+str(knot))

def plot_homology_2D(knot,convention='NS',gridlines=True,s_inv_circles_on=True, vars='qa'):
    PS=compute_planar_support(knot)
    if PS==1:
        plot_homology_planar(knot,convention,gridlines,s_inv_circles_on,vars)
    elif PS==2:
        plot_homology_nonplanar(knot,convention,gridlines,s_inv_circles_on,vars)
    else:
        print('knot data is not in the database')

#computes list of all delta gradings
def size_of_homology(knot):
    if knot in homfly_data:
        size = 0
        homology = homfly_data[knot]
        for key in list(homology.keys()):
            size+=homology[key]
        return size
    else:
        return 0

size_of_homology_vect = np.vectorize(size_of_homology)

#computes if the knot is parity (i.e. all homological gradings have same parity)
def parity(knot):
    if knot in homfly_data:
        homology = homfly_data[knot]
        return check_parity(convert_homology(homology))
    else:
        return 9999

parity_vect = np.vectorize(parity)

def detect_int(col):
    if col.dtype=='int64':
        return True
    else:
        return False

def detect_float(col):
    if col.dtype=='float64':
        return True
    else:
        return False

def floatize(number):
    return float(number)

float_vect= np.vectorize(floatize)

#turns alternating column from "Y", "N" to 1 and 0
def convert_YN(boo):
    if boo=='Y' or boo==1 or boo==1.0:
        return 1.0
    elif boo=="N" or boo==0 or boo==0.0:
        return 0.0
    else:
        return float('NAN')

convert_YN_vect = np.vectorize(convert_YN)

def convert_str_to_float(val):
    boo = True
    while boo:
        try:
            return float(val)
            boo=False
            break
        except ValueError:
            return float('NaN')
            boo=False
            break

def convert_small_or_Large(val):
    boo = True
    while boo:
        try:
            if val=='Small':
                return 0.0
            elif val=='Large':
                return 1.0
            boo=False
            break
        except ValueError:
            return float('NaN')
            boo=False
            break

def convert_L_space(val):
    boo = True
    while boo:
        try:
            if val=='No':
                return 0.0
            elif val=='Yes':
                return 1.0
            boo=False
            break
        except ValueError:
            return float('NaN')
            boo=False
            break

def float_to_boolean(flo):
    if flo==1.0:
        return True
    else:
        return False

#scores a model (inpute e.g. (all_knot_info_int, features, 'Size of Homology'))
def model_score(dataframe,features, target, algo = 'linear', scoring = 'accuracy'):
    yyy=dataframe[target]
    XXX=dataframe[features]
    XXXscaler = preprocessing.StandardScaler().fit(XXX)
    XXX_scaled = XXXscaler.transform(XXX)
    if algo=='linear':
        reg = LinearRegression()
    elif algo=='logistic':
        reg = LogisticRegression()
    if scoring=='accuracy':
        return cross_val_score(reg,XXX_scaled,yyy,cv=4).mean()
    elif scoring=='f1':
        precision = cross_val_score(reg,XXX_scaled,yyy,cv=3, scoring='precision_micro').mean()
        recall = cross_val_score(reg, XXX_scaled, yyy, cv=3, scoring='recall_macro').mean()
        return precision*recall/(precision+recall)

#CODE FOR THE GREEDY ALGORITHM FOR FEATURE SELECTION
def greedy(dataframe, features, subsetsize, target, algo = 'linear', scoring = 'accuracy'):
    i=0
    best_stats=[]
    s=set(features)
    subsets=list(map(set, itertools.combinations(s, subsetsize))) #subsets of size (subsetsize)
    possible_stat_dict = {}
    scores={0:0}
    for stat_pair in subsets:
        possible_stat_dict[tuple(stat_pair)]=0
    while (i==0) or (scores[i]>scores[i-1]):
        i+=1
        for stat_pair in list(possible_stat_dict.keys()):
            stats_temp = best_stats+list(stat_pair)
            possible_stat_dict[tuple(stat_pair)]=model_score(dataframe,stats_temp,target,algo, scoring)
        max_key = max(possible_stat_dict, key=possible_stat_dict.get)
        best_stats.extend(list(max_key))
        scores[i]=possible_stat_dict[max_key]
        possible_stat_dict.pop(max_key)
        print(best_stats,scores[i])
    return (best_stats[:-subsetsize], scores[i-1])

def is_unimodal(seq):
    M=seq.index(max(seq))
    return all(seq[i]<=seq[i+1] for i in range(M)) and all(seq[i]>=seq[i+1] for i in range(M+1,len(seq)-1))


def check_unimodal(knot):
    homology = homfly_data[knot]
    grad_dict={}
    seq_dict={}
    for x in set([grad[0] for grad in homology.keys()]):
        grad_dict[x]=[grad for grad in homology.keys() if grad[0]==x]
        grad_dict[x].sort()
        seq_dict[x]=[homology[grad] for grad in grad_dict[x]]
    return all(is_unimodal(seq) for seq in seq_dict.values())
