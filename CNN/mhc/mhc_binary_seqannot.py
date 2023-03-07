import torch
import argparse
import os, sys
import torchvision
import torchvision.transforms as transforms
import matplotlib.pyplot as plt
import itertools
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
import pandas as pd
from sklearn.metrics import f1_score, roc_auc_score
import shap


parser = argparse.ArgumentParser(description='gRNA prediction model')
parser.add_argument('--date', help='date used in output name to versonize the results')
#parser.add_argument('--savedir', default='./', help='path to save results')
#parser.add_argument('--ckptdir', default='./ckpt', help='path to save checkpoints')
#parser.add_argument('--batch-size', type=int, default=128,
#                    help='input batch size for training (default: 128)')
#parser.add_argument('--epochs', type=int, default=100,
#                    help='number of epochs to train (default: 100)')
#parser.add_argument('--lr', type=float, default=0.001,
#                    help='learning rate (default: 0.001)')
parser.add_argument('--fold', type=int, default=1, help='which fold of data to use')
parser.add_argument('--celltype', default="k562", help='which celltype to train model')
args = parser.parse_args()
#
#savedir = args.savedir
#ckptdir = args.ckptdir


#batch_size = args.batch_size
#epochs = args.epochs
#lr = args.lr
#ngpu=1

celltype = args.celltype
datadir = '/proj/yunligrp/users/tianyou/gRNA/data/mhc/'
resultdir = os.path.join('/proj/yunligrp/users/tianyou/gRNA/result/mhc/', celltype)
batch_size = 512
epochs = 100
lr = 0.00005
ngpu=1
fold = args.fold

device = torch.device("cuda:0" if (torch.cuda.is_available() and ngpu > 0) else "cpu")
print(device)

#torch.manual_seed(762)

def preprocess_seq(data):
    print("Start preprocessing the sequence")
    length = len(data[0])
    DATA_X = np.zeros((len(data),4,length), dtype=int)
    print(DATA_X.shape)
    for l in range(len(data)):
        if l % 10000 == 0:
            print(str(l) + " sequences processed")
        for i in range(length):
            try: data[l][i]
            except: print(data[l], i, length, len(data))
            if data[l][i]in "Aa":
                DATA_X[l, 0, i] = 1
            elif data[l][i] in "Cc":
                DATA_X[l, 1, i] = 1
            elif data[l][i] in "Gg":
                DATA_X[l, 2, i] = 1
            elif data[l][i] in "Tt":
                DATA_X[l, 3, i] = 1
            else:
                print("Non-ATGC character " + data[i])
                sys.exit()
    print("Preprocessing the sequence done")
    return DATA_X



dat = pd.read_csv(os.path.join(datadir, celltype, celltype+'-pfdr0.05-pfdr0.2-binary-fold'+str(fold)+'-train.csv'), index_col = False)
sequence = dat['grna']
sequence_onehot = preprocess_seq(sequence)
sequence_sum = sequence_onehot.sum(axis=2)
dat['Acount'] = sequence_sum[:,0]
dat['Ccount'] = sequence_sum[:,1]
dat['Gcount'] = sequence_sum[:,2]
dat['Tcount'] = sequence_sum[:,3]
dat['GCcount'] = sequence_sum[:,1] + sequence_sum[:,2]
dat['GCprop'] = dat['GCcount'] / sequence_onehot.shape[2]
label = dat['significant'].to_numpy(dtype = np.float32)
class_count = dat['significant'].value_counts()
w = class_count[0] / class_count[1]
feas_sel = ["deltagb", "deltagh", "GCcount", "GCprop", "Acount", "Ccount", "Tcount", "Gcount", "OGEE_prop_Essential", "H3k27ac_CPM_1Kb_new", 
            "ATAC_CPM_1Kb_new", "H3K4me3_CPM_1Kb_new"]
annotation = dat.loc[:,feas_sel].to_numpy(dtype = np.float32)


X1 = torch.tensor(sequence_onehot, dtype=torch.float32)
#Xloader = torch.utils.data.DataLoader(X, batch_size=batch_size, shuffle=True)
X2 = torch.tensor(annotation, dtype=torch.float32)
Y = torch.tensor(label, dtype=torch.float32)
Y = Y.view(-1, 1)
#Yloader = torch.utils.data.DataLoader(Y, batch_size=batch_size, shuffle=True)
input_dat = TensorDataset(X1,X2,Y)
datloader = DataLoader(input_dat, batch_size=batch_size, shuffle=True)


## test set
def get_test_data(cell):
    test = pd.read_csv(os.path.join(datadir, cell, cell+'-pfdr0.05-pfdr0.2-binary-fold'+str(fold)+'-test.csv'), index_col = False)
    test_sequence = test['grna']
    test_sequence_onehot = preprocess_seq(test_sequence)
    test_sequence_sum = test_sequence_onehot.sum(axis=2)
    test['Acount'] = test_sequence_sum[:,0]
    test['Ccount'] = test_sequence_sum[:,1]
    test['Gcount'] = test_sequence_sum[:,2]
    test['Tcount'] = test_sequence_sum[:,3]
    test['GCcount'] = test_sequence_sum[:,1] + test_sequence_sum[:,2]
    test['GCprop'] = test['GCcount'] / test_sequence_onehot.shape[2]
    test_label = test['significant'].to_numpy(dtype = np.float32)
    test_annotation = test.loc[:,feas_sel].to_numpy(dtype = np.float32)

    test_X1 = torch.tensor(test_sequence_onehot, dtype=torch.float32).to(device)
    test_X2 = torch.tensor(test_annotation, dtype=torch.float32).to(device)
    return [test_X1, test_X2, test_label, test['grna']]

k562_X1, k562_X2, k562_label, k562_grna = get_test_data("k562")
npc_X1, npc_X2, npc_label, npc_grna = get_test_data("npc")
ipsc_X1, ipsc_X2, ipsc_label, ipsc_grna = get_test_data("ipsc")
testc_X1, testc_X2, testc_label, testc_grna = get_test_data(celltype)

dim_fc = 112

class DeepSeqCNN(nn.Module):
    def __init__(self):
        super(DeepSeqCNN, self).__init__()
        self.conv0 = nn.Sequential(
            nn.Conv1d(4, 50, 2, stride=2),
            nn.ReLU(),
            nn.Dropout(0.4),   ## for ReLU, it is interchangeable with max pooling and dropout
        )
        self.conv1 = nn.Sequential(
            nn.Conv1d(4, 100, 3, padding=1),
            nn.MaxPool1d(2),
            nn.ReLU(),
            nn.Dropout(0.4),   ## for ReLU, it is interchangeable with max pooling and dropout
        )
        self.conv2 = nn.Sequential(
            nn.Conv1d(4, 70, 5, padding=2),
            nn.MaxPool1d(2),
            nn.ReLU(),
            nn.Dropout(0.4),   ## for ReLU, it is interchangeable with max pooling and dropout
        )
        self.conv3 = nn.Sequential(
            nn.Conv1d(4, 40, 7, padding=3),
            nn.MaxPool1d(2),  ## Their codes seemed to use average pooling but texts suggest max pooling?
            nn.ReLU(),
            nn.Dropout(0.4),   ## for ReLU, it is interchangeable with max pooling and dropout
        )
        self.conv4 = nn.Sequential(
            nn.Conv1d(4, 1, 1),
            nn.ReLU(),
        )
        self.fc1 = nn.Sequential(
            nn.Linear(260*10, 80),
            nn.ReLU(),
            nn.Dropout(0.3))
        self.fc2 = nn.Sequential(
            nn.Linear(dim_fc, 80),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(80, 60),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(60, 40),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(40, 1),
            # nn.Sigmoid()  ## BCEWithLogitsLoss takes in logits directly without sigmoid
        )
        
    def forward(self, x, y):
        x0 = self.conv0(x)
        x1 = self.conv1(x)
        x2 = self.conv2(x)
        x3 = self.conv3(x)
        x4 = self.conv4(x).view(-1, 20)
        x_concat = torch.cat((x0, x1, x2, x3), dim=1) # size: [:,260,10]
        x_concat = x_concat.view(-1, 260*10)
        x_concat = self.fc1(x_concat)
        xy_concat = torch.cat((x_concat, x4, y), dim = 1)
        xy_concat = self.fc2(xy_concat)
        #for layer in self.fc:
        #    x_concat = layer(x_concat)
        #    print(x_concat.size())
        #x_concat = self.fc(x_concat)
        return xy_concat


CNN = DeepSeqCNN().to(device)
optimizer = optim.Adam(CNN.parameters(), lr=lr)
#lossfunc = nn.L1Loss().to(device)
lossfunc = nn.BCEWithLogitsLoss(pos_weight=torch.Tensor([w])).to(device)
sigmoid = nn.Sigmoid()

def train_model(model, num_epochs):
    for epoch in range(num_epochs):
        # Training
        if epoch % 2 == 0:
            model.eval()
            test_predict = sigmoid(model(testc_X1, testc_X2))
            test_predict_np = test_predict.detach().to('cpu').numpy()
            auc = roc_auc_score(testc_label, test_predict_np)
            print('Epoch [%d] AUC: %.3f' %
                    (epoch + 1, auc))
            model.train()
        running_loss = 0.0
        for i, batch in enumerate(datloader, 0):
            # Transfer to GPU
            local_x1, local_x2, local_y = batch
            local_x1, local_x2, local_y = local_x1.to(device), local_x2.to(device), local_y.to(device)
            optimizer.zero_grad()
            FC_pred = model(local_x1, local_x2)
            loss = lossfunc(FC_pred, local_y)
            loss.backward()
            optimizer.step()
            running_loss += loss.item()
            if i % 20 == 19:    # print every 200 mini-batches
                print('[%d, %5d] loss: %.3f' %
                    (epoch + 1, i + 1, running_loss / 20))
                running_loss = 0.0
    
    return model

CNN = train_model(CNN, num_epochs=epochs)

ckptPATH = os.path.join(resultdir,'models',celltype+'-binary-BCE-seqannot-'+args.date+'-fold'+str(fold)+'.pth')
torch.save(CNN.state_dict(), ckptPATH)

del testc_X1, testc_X2, testc_label, testc_grna
CNN.eval()

k562_predict = sigmoid(CNN(k562_X1, k562_X2))
k562_predict_np = k562_predict.detach().to('cpu').numpy()
roc_auc_score(k562_label, k562_predict_np)
k562_PD = pd.DataFrame(np.stack((k562_grna, k562_label, k562_predict_np[:,0]), axis=1), columns = ['grna', 'true', 'predict'])
k562_PD.to_csv(os.path.join(resultdir,celltype+'-binary-BCE-seqannot-fold'+str(fold)+'-test-k562.csv'))

npc_predict = sigmoid(CNN(npc_X1, npc_X2))
npc_predict_np = npc_predict.detach().to('cpu').numpy()
roc_auc_score(npc_label, npc_predict_np)
npc_PD = pd.DataFrame(np.stack((npc_grna, npc_label, npc_predict_np[:,0]), axis=1), columns = ['grna', 'true', 'predict'])
npc_PD.to_csv(os.path.join(resultdir,celltype+'-binary-BCE-seqannot-fold'+str(fold)+'-test-npc.csv'))

ipsc_predict = sigmoid(CNN(ipsc_X1, ipsc_X2))
ipsc_predict_np = ipsc_predict.detach().to('cpu').numpy()
roc_auc_score(ipsc_label, ipsc_predict_np)
ipsc_PD = pd.DataFrame(np.stack((ipsc_grna, ipsc_label, ipsc_predict_np[:,0]), axis=1), columns = ['grna', 'true', 'predict'])
ipsc_PD.to_csv(os.path.join(resultdir,celltype+'-binary-BCE-seqannot-fold'+str(fold)+'-test-ipsc.csv'))
