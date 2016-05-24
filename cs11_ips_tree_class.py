
# coding: utf-8

# In[ ]:

###Draw cell lineage tree, 2016-05-24, by Anatole.King
###It is specified to lineage tracing data, just for test GitHub.

import matplotlib
matplotlib.use('SVG')
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

path = '/Users/bioinformatics/Projects/live_cell/The_tree/draw_tree'
os.chdir(path)

#data = pd.read_table('v4_osk.tsv', index_col='Name') #bad to use pandas, should transform data to dic

Cells = []
oh = open('v4_osk.tsv', 'r')
buffer = oh.readline() # header line

while 1:
    line1 = oh.readline().strip('\n').split('\t')
    line2 = oh.readline().strip('\n').split('\t')
    
    if len(line2) == 1:
        break
    NewCell = {
        'name': line1[0],
        'p1': int(line1[1]),
        'p2-1': int(line1[2]),
        'p2-2': int(line1[3]),
        'p3-1': int(line1[4]),
        'p3-2': int(line1[5]),
        'p3-3': int(line1[6]),
        'p3-4': int(line1[7]),
        'p4-1': int(line1[8]),
        'p4-2': int(line1[9]),
        'p4-3': int(line1[10]),
        'p4-4': int(line1[11]),
        'p4-5': int(line1[12]),
        'p4-6': int(line1[13]),
        'p4-7': int(line1[14]),
        'p4-8': int(line1[15]),
        'note_p1': line2[1],
        'note_p2-1': line2[2],
        'note_p2-2': line2[3],
        'note_p3-1': line2[4],
        'note_p3-2': line2[5],
        'note_p3-3': line2[6],
        'note_p3-4': line2[7],
        'note_p4-1': line2[8],
        'note_p4-2': line2[9], 
        'note_p4-3': line2[10],
        'note_p4-4': line2[11],
        'note_p4-5': line2[12],
        'note_p4-6': line2[13],
        'note_p4-7': line2[14],
        'note_p4-8': line2[15]  
    }
    Cells.append(NewCell)
    
#draw lineage map

def ipstree(data, name, lower, upper):
    def NodeLoc(Ang):
        ym = {}
        ym['p1'] = 0.0
        ym['p2-1'] = ym['p1'] - Ang/4
        ym['p2-2'] = ym['p1'] + Ang/4
        ym['p3-1'] = ym['p2-1'] - Ang/8
        ym['p3-2'] = ym['p2-1'] + Ang/8
        ym['p3-3'] = ym['p2-2'] - Ang/8
        ym['p3-4'] = ym['p2-2'] + Ang/8
        ym['p4-1'] = ym['p3-1'] - Ang/16
        ym['p4-2'] = ym['p3-1'] + Ang/16
        ym['p4-3'] = ym['p3-2'] - Ang/16
        ym['p4-4'] = ym['p3-2'] + Ang/16
        ym['p4-5'] = ym['p3-3'] - Ang/16
        ym['p4-6'] = ym['p3-3'] + Ang/16
        ym['p4-7'] = ym['p3-4'] - Ang/16
        ym['p4-8'] = ym['p3-4'] + Ang/16
        return ym
    all_cell_labels = ['p1', 
            'p2-1', 'p2-2',
            'p3-1', 'p3-2', 'p3-3','p3-4',
            'p4-1','p4-2','p4-3','p4-4','p4-5','p4-6','p4-7','p4-8']

    all_connections = [('p1', 'p2-1'), ('p1', 'p2-2'),
        ('p2-1', 'p3-1'), ('p2-1', 'p3-2'),
        ('p2-2', 'p3-3'), ('p2-2', 'p3-4'),
        ('p3-1', 'p4-1'), ('p3-1', 'p4-2'),
        ('p3-2', 'p4-3'), ('p3-2', 'p4-4'),
        ('p3-3', 'p4-5'), ('p3-3', 'p4-6'),
        ('p3-4', 'p4-7'), ('p3-4', 'p4-8'),]
    Angle = 2 * np.pi / 8  # draw 8 trees 
    y_modifer = NodeLoc(Angle)

    red_scatters_o = [[], []] # gfp-
    gre_scatters_o = [[], []] # gfp+
    red_scatters_x = [[], []] # dead
    purple_scatters_s = [[], []] # arrest
    early_scatters_o = [[], []] # p1 -- p3
    y_labs = []

    Arrow = {}
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111, projection='polar')

    for actual_yoff, c in enumerate(data[lower:upper]):
        actual_yoff = np.pi/8 + np.pi * actual_yoff/4
        y_labs.append(c['name'])
        for node in all_cell_labels:
            if c['note_%s'% node] == 'none':
                pass
            else:
                if c[node]  > 144:
                    c[node] = 140  ## the large nummber is argly to show
                Xloc = c[node]
                Yloc = actual_yoff + y_modifer[node]
                c[node] = (Xloc, Yloc)

                # work out the appropriate color:
                if c['note_%s' % node] in ('gp', 'gpA', 'gpB'):
                    gre_scatters_o[0].append(Xloc)
                    gre_scatters_o[1].append(Yloc)
                elif c['note_%s' % node] in ('dd','ds'): # cell dead
                    red_scatters_x[0].append(Xloc)
                    red_scatters_x[1].append(Yloc)
                elif c['note_%s' % node] in ('stop','ar'):
                    purple_scatters_s[0].append(Xloc)
                    purple_scatters_s[1].append(Yloc)
                elif c['note_%s' % node] == 'gn':
                    red_scatters_o[0].append(Xloc)
                    red_scatters_o[1].append(Yloc)
                else:
                    early_scatters_o[0].append(Xloc)
                    early_scatters_o[1].append(Yloc)
        for edge in all_connections:
            if c['note_%s'% edge[0]] != 'dd':
                if c['note_%s' % edge[1]] not in ('stop', 'none'):
                    ax.annotate("", xytext=(c[edge[0]][1], c[edge[0]][0]), 
                                xy=(c[edge[1]][1], c[edge[1]][0]), alpha=0.1, color='grey', 
                                arrowprops=dict(arrowstyle="->"), zorder=0)  

    ax.scatter(red_scatters_o[1], red_scatters_o[0], c='red',  s=30, edgecolor='none')
    ax.scatter(gre_scatters_o[1], gre_scatters_o[0], c='green',  s=30, edgecolor='none')
    ax.scatter(red_scatters_x[1], red_scatters_x[0], c='Fuchsia', marker='x', s=30, edgecolor='none')
    ax.scatter(purple_scatters_s[1], purple_scatters_s[0], c='purple', marker='*', s=60, edgecolor='none')
    ax.scatter(early_scatters_o[1], early_scatters_o[0], c='DeepSkyBlue', marker='o', s=30, edgecolor='none')


    ax.set_rmax(144)
    ax.grid(True)
    ax.set_rgrids(radii=[48, 96, 144] , labels = ['48h', '96h', '144h'], angle=0)
    ax.set_xticklabels([]*len(data))
    plt.legend(labels=['GFP negative colony', 'GFP positive colony', 'Dead', 'No division later', 
                       'Intermediate cell'])
    plt.title('Reprogramming Tree')
    plt.savefig('ips_tree_%s'%name)
    #plt.show()
    plt.close()
    
## Draw 8 tress per fig
for i in np.arange(int(len(Cells)/8) +1):
    if (i+1)*8 <= len(Cells):
        lower = i*8
        upper = (i+1)*8
        ipstree(Cells, 'v4_osk_round%s'%i, lower, upper)
    else:
        lower = i*8
        upper = len(Cells)
        ipstree(Cells, 'v4_osk_round%s'%i, lower, upper)


