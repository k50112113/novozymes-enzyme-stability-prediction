import torch
import torch.nn as nn

torch.set_default_tensor_type(torch.DoubleTensor)

# allseq = ""
# for val in data.values():
#     allseq += val[0]
# for i in range(ord('A'),ord('Z')+1):
#     cnt = allseq.count(chr(i))
#     if cnt > 0: print(chr(i), end="")
# amino_list = 'ACDEFGHIKLMNPQRSTVWY'
# oh = torch.nn.functional.one_hot(torch.arange(0, 20), num_classes = 20)
# amino_oh = {}
# for i, a in enumerate(amino_list):
#     amino_oh[a] = oh[i]
# for keys in amino_oh:
#     print(keys, amino_oh[keys])

class EzDataset(torch.utils.data.Dataset):
    def __init__(self, seq_id, D, seq_len, ph, tm):
        self.seq_id = seq_id
        self.D = D
        self.seq_len = seq_len
        self.ph = ph
        self.tm = tm

    def to(self, device = 'cpu'):
        self.seq_id = seq_id.to(device)
        self.D = D.to(device)
        self.seq_len = seq_len.to(device)
        self.ph = ph.to(device)
        self.tm = tm.to(device)

    def __len__(self):
        return len(self.D)

    def __getitem__(self, index):
        return self.seq_id[index], self.D[index], self.seq_len[index], self.ph[index], self.tm[index]

class Descriptor:
    def __init__(self, Ds_step = 0.001, eta = 0.001):
        self.Ds = torch.linspace(0., 1., steps=int(1./Ds_step+1))
        self.eta = eta
        self.amino_list = 'ACDEFGHIKLMNPQRSTVWY'
        self.amino_id = {}
        for index, amino in enumerate(self.amino_list):
            self.amino_id[amino] = index
        pass
    
    def encode(self, protein_sequence):
        seq_len = len(protein_sequence)
        inv_seq_len = 1.0/seq_len
        D = torch.zeros(len(self.amino_list), len(self.Ds))
        for index, amino in enumerate(protein_sequence):
            D[self.amino_id[amino]] += (-seq_len*self.eta*(self.Ds - index*inv_seq_len)**2).exp()
        D *= inv_seq_len
        # with open("d.txt","w") as fout:
        #     for i in range(len(self.Ds)):
        #         fout.write("%s\t"%self.Ds[i].item())
        #         for j in range(len(self.amino_list)):
        #             fout.write("%s\t"%D[j][i].item())
        #         fout.write('\n')
        return D, seq_len
    
    def get_descriptor_size(self):
        return len(self.amino_list)*len(self.Ds)
        

class EzStabilityMLP(nn.Module):
    def __init__(self, input_size, hidden_layer_size, activation_type = 'tanh'):
        super().__init__()
        self.input_size = input_size
        self.assign_activation(activation_type)
        self.assign_layer_size(hidden_layer_size)
        self.linear_layer = nn.ModuleList()
        for layer_index in range(len(self.layer_size)-1):
            self.linear_layer.append(nn.Linear(self.layer_size[layer_index],self.layer_size[layer_index+1]))

    def assign_activation(self, activation_type):
        if activation_type == 'tanh': self.activate = nn.Tanh()
        elif activation_type == 'relu': self.activate = nn.ReLU()
        elif activation_type == 'silu': self.activate = nn.SiLU()

    def assign_layer_size(self, hidden_layer_size):
        self.layer_size = [self.input_size] + hidden_layer_size + [1]
        
    def forward(self, x):
        y = x
        for layer_index in range(len(self.layer_size)-2):
            y = self.linear_layer[layer_index](y)
            y = self.activate(y)
        y = self.linear_layer[len(self.layer_size)-2](y)
        return y

    def get_layer_size(self):
        return self.layer_size

import time

class Clock:
    def __init__(self):
        self.reset_timer()

    def get_dt(self):
        self.time2 = time.perf_counter()
        dt = self.time2-self.time1
        self.time1 = self.time2
        return dt

    def reset_timer(self):
        self.time1 = time.perf_counter()
        self.time2 = None
    

