import os
import sys

protein_seq_map = {}

wild = "VPVNPEPDATSVENVALKTGSGDSQSDPIKADLEVKGQSALPFDVDCWAILCKGAPNVLQRVNEKTKNSNRDRSGANKGPFKDPQKWGIKALPPKNPSWSAQDFKSPEEYAFASSLQGGTNAILAPVNLASQNSQGGVLNGFYSANKVAQFDPSKPQQTKGTWFQITKFTGAAGPYCKALGSNDKSVCDKNKNIAGDWGFDPAKWAYQYDEKNNKFNYVGK"

with open("../test.csv","r") as fin:
    fin.readline()
    for aline in fin:
        linelist = aline.strip().split(',')
        seq_id.append(int(linelist[0]))
        protein_seq_map[linelist[1]] = int(linelist[0])

p = os.popen("find %s -type %s -name %s"%('.','d','\"*\"'))
folder = p.read().split('\n')
folder.pop()
del folder[0]
folder = sorted(folder)
for i in range(len(folder)): folder[i] = folder[i].replace('./', "")

fout = open("single-mutation-map.txt","w")

for i in range(len(folder)):
    if folder[i] == 'wild':
        continue
    
    if folder[i] == 'WT':
        mutation = wild
    elif ord('Z') >= ord(folder[i][-1]) >= ord('A'):
        pos = int(folder[i][1:-1])
        mutation = wild[:pos-1] + folder[i][-1] + wild[pos:]
    else:
        pos = int(folder[i][1:])
        mutation = wild[:pos-1] + wild[pos:]
    fout.write("%s %d\n"%(folder[i], protein_seq_map[mutation]))
fout.close()
