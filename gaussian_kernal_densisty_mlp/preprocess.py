import torch
import pickle
import sys
from model import EzDataset, Descriptor

torch.set_default_tensor_type(torch.DoubleTensor)

Ds_step = float(sys.argv[1])
eta = float(sys.argv[2])

seq_id = []
D = []
seq_len = []
ph = []
tm = []
dspr = Descriptor(Ds_step = Ds_step, eta = eta)
print(dspr.Ds)
with open("train.csv","r") as fin:
    keys = fin.readline()
    for aline in fin:
        linelist = aline.strip().split(',')
        if linelist[1] != '':
            if linelist[2] == '': linelist[2] = '7.0'
            seq_id_instance = int(linelist[0])
            protein_sequence = linelist[1]
            D_instance, seq_len_instance = dspr.encode(protein_sequence)
            ph_instance = float(linelist[2])-7.0
            tm_instance = float(linelist[4])

            D.append(D_instance)
            seq_len.append(seq_len_instance)
            ph.append(ph_instance)
            tm.append(tm_instance)
            seq_id.append(seq_id_instance)
            print(seq_id_instance, end='\r')
            # if len(D) > 100: break

seq_id = torch.tensor(seq_id)
D = torch.stack(D)
seq_len = torch.tensor(seq_len, dtype=torch.float64)
ph = torch.tensor(ph)
tm = torch.tensor(tm)
print("Sample size = %d"%(len(seq_id)))
ezdata = EzDataset(seq_id, D, seq_len, ph, tm)
fout = open("preprocessed_train_%s_%s.pickle"%(Ds_step, eta),"wb")
pickle.dump((dspr, ezdata), fout)