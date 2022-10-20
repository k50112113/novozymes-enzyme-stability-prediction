import torch
import torch.nn as nn
import sys
import pickle
from model import EzDataset, Descriptor, EzStabilityMLP, Clock


torch.set_default_tensor_type(torch.DoubleTensor)

Ds_step = float(sys.argv[1])
eta = float(sys.argv[2])
test_ratio          = 0.1
batch_size          = 50
random_split_seed   = 42
number_of_epochs    = 1000
lr_start            = 1e-4
weight_decay        = 0.0
lr_decay_rate       = 0.9
lr_decay_step       = 20

fin = open("preprocessed_train_%s_%s.pickle"%(Ds_step, eta), "rb")
dspr, ezdata = pickle.load(fin)
fout = open("output_%s_%s.txt"%(Ds_step, eta),"w")

# ezdata_sorted = sorted(ezdata, key=lambda x: x[2])
# print(ezdata_sorted[0][0],ezdata_sorted[0][2],ezdata_sorted[0][-1],ezdata_sorted[-1][0],ezdata_sorted[-1][2],ezdata_sorted[-1][-1])
# exit()
total_n_test_sample  = int(test_ratio*len(ezdata))
total_n_train_sample = len(ezdata) - total_n_test_sample
print("number of training samples = %s"%(total_n_train_sample))
print("number of testing samples  = %s"%(total_n_test_sample))
fout.write("number of training samples = %s\n"%(total_n_train_sample))
fout.write("number of testing samples  = %s\n"%(total_n_test_sample))

train_dataset, test_dataset = torch.utils.data.random_split(ezdata, [total_n_train_sample, total_n_test_sample], generator=torch.Generator().manual_seed(random_split_seed))
train_loader = torch.utils.data.DataLoader(dataset=train_dataset, batch_size=batch_size,shuffle=True)
test_loader  = torch.utils.data.DataLoader(dataset=test_dataset, batch_size=batch_size,shuffle=False)

activation_type = 'relu'
model = EzStabilityMLP(input_size = dspr.get_descriptor_size() + 2, \
                       hidden_layer_size = [1200, 560, 480, 320, 240, 120], \
                       activation_type = activation_type)
print("model layer: %s"%(model.get_layer_size()))             
print("activation: %s"%(activation_type))      
fout.write("model layer: %s\n"%(model.get_layer_size()))      
fout.write("activation: %s\n"%(activation_type))

optimization = torch.optim.Adam(model.parameters(), lr=lr_start, weight_decay = weight_decay)
learning_rate_schedule = torch.optim.lr_scheduler.StepLR(optimization, step_size = lr_decay_step, gamma = lr_decay_rate)
loss_function = nn.MSELoss(reduction='mean')
number_of_batches_in_train_loader = len(train_loader)
number_of_batches_in_test_loader  = len(test_loader)
#print(number_of_batches_in_train_loader, number_of_batches_in_test_loader)

timer = Clock()
timer.get_dt()
print("epoch train_loss test_loss lr time(s)")
fout.write("epoch train_loss test_loss lr time(s)\n")
for ith_epoch in range(number_of_epochs):
    #Train an epoch
    train_loss = 0.
    for jth_batch, train_data in enumerate(train_loader):
        seq_id, D, seq_len, ph, tm = train_data
        optimization.zero_grad()
        x = torch.cat((D.flatten(start_dim = -2), seq_len.unsqueeze(-1), ph.unsqueeze(-1)), dim = 1)
        tm_predict = model(x)
        loss = loss_function(tm.unsqueeze(-1), tm_predict)
        loss.backward()
        optimization.step()
        train_loss += loss.item()*tm.shape[0]
        print("training %d/%d batches"%(jth_batch, number_of_batches_in_train_loader), end='\r')
    train_loss /= number_of_batches_in_train_loader
    #Test
    test_loss = 0.
    with torch.no_grad():
        for jth_batch, test_data in enumerate(test_loader):
            seq_id, D, seq_len, ph, tm = test_data
            x = torch.cat((D.flatten(start_dim = -2), seq_len.unsqueeze(-1), ph.unsqueeze(-1)), dim = 1)
            tm_predict = model(x)
            loss = loss_function(tm.unsqueeze(-1), tm_predict)
            test_loss += loss.item()*tm.shape[0]
            print("testing %d/%d batches"%(jth_batch, number_of_batches_in_test_loader), end='\r')
    test_loss /= number_of_batches_in_test_loader
    dt = timer.get_dt()
    print("%d %.4e %.4e %.3e %.3f"%(ith_epoch, train_loss, test_loss, optimization.param_groups[0]["lr"], dt))
    fout.write("%d %.4e %.4e %.3e %.3f\n"%(ith_epoch, train_loss, test_loss, optimization.param_groups[0]["lr"], dt))
    learning_rate_schedule.step()
    
fout.close()