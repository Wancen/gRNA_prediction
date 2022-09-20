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
from scipy.stats import spearmanr
from sklearn.metrics import fbeta_score
import shap


datadir = '/proj/yunligrp/users/tianyou/gRNA/data/data_April_resplit/'
resultdir = '/proj/yunligrp/users/tianyou/gRNA/result_April_resplit/logfc_3cat/'
batch_size = 256
epochs = 4000 ## 4000 for promoter winsorized
lr = 0.0001
ngpu=1
grp = 'enh'
thres = 0.1

device = torch.device("cuda:0" if (torch.cuda.is_available() and ngpu > 0) else "cpu")
print(device)



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



dat = pd.read_csv(datadir+'wgCERES-gRNAs-k562-discovery-screen-'+grp+'_logFC'+str(thres)+'-train-clean.csv', index_col = False)
sequence = dat['protospacer'].reset_index(drop=True)
sequence_onehot = preprocess_seq(sequence)
log2FC_cat = dat['class']
if grp == 'pro':
    annotation = dat.iloc[:,12:54].to_numpy(dtype = np.float32) # for promoters
elif grp == 'enh':
    annotation = dat.iloc[:,12:56].to_numpy(dtype = np.float32) # for enhancers
else:
    print("Invalid group: " + grp)

class_count = np.unique(log2FC_cat, return_counts=True)[1]
#w = class_count / class_count[1]
w = class_count[1] / class_count

X1 = torch.tensor(sequence_onehot, dtype=torch.float32)
#Xloader = torch.utils.data.DataLoader(X, batch_size=batch_size, shuffle=True)
X2 = torch.tensor(annotation, dtype=torch.float32)
Y = torch.tensor(log2FC_cat, dtype=torch.long)
#Y = Y.view(-1, 1)
#Yloader = torch.utils.data.DataLoader(Y, batch_size=batch_size, shuffle=True)
input_dat = TensorDataset(X1,X2,Y)
datloader = DataLoader(input_dat, batch_size=batch_size, shuffle=True)


## test set
test = pd.read_csv(datadir+'wgCERES-gRNAs-k562-discovery-screen-'+grp+'_logFC'+str(thres)+'-test-clean.csv', index_col = False)
test_sequence = test['protospacer'].reset_index(drop=True)
test_sequence_onehot = preprocess_seq(test_sequence)
test_log2FC_cat = test['class']
if grp == 'pro':
    test_annotation = test.iloc[:,12:54].to_numpy(dtype = np.float32) #promoters
elif grp == "enh":
    test_annotation = test.iloc[:,12:56].to_numpy(dtype = np.float32) #enhancers
else:
    print("Invalid group:" + grp)

#subsample = np.random.choice(len(test_sequence), size = 3000, replace = False)
#test_X_sub = torch.tensor(test_sequence_onehot[subsample,:], dtype=torch.float32).to(device)
test_X1_sub = torch.tensor(test_sequence_onehot, dtype=torch.float32).to(device)
test_X2_sub = torch.tensor(test_annotation, dtype=torch.float32).to(device)

if grp == 'pro':
    dim_fc = 142
elif grp == 'enh':
    dim_fc = 144

class DeepSeqCNN(nn.Module):
    def __init__(self):
        super(DeepSeqCNN, self).__init__()
        #self.input_size = input_size
        self.conv0 = nn.Sequential(
            nn.Conv1d(4, 50, 1),
            nn.MaxPool1d(2),
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
        self.fc1 = nn.Sequential(
            nn.Linear(260*10, 100),
            nn.ReLU(),
            nn.Dropout(0.3))
        self.fc2 = nn.Sequential(
            nn.Linear(dim_fc, 100),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(100, 80),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(80, 60),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(60, 3),
            nn.Softmax(dim=1)
        )
        
    def forward(self, x, y):
        x0 = self.conv0(x)
        x1 = self.conv1(x)
        x2 = self.conv2(x)
        x3 = self.conv3(x)
        x_concat = torch.cat((x0, x1, x2, x3), dim=1) # size: [:,260,10]
        x_concat = x_concat.view(-1, 260*10)
        x_concat = self.fc1(x_concat)
        xy_concat = torch.cat((x_concat, y), dim = 1)
        xy_concat = self.fc2(xy_concat)
        #for layer in self.fc:
        #    x_concat = layer(x_concat)
        #    print(x_concat.size())
        #x_concat = self.fc(x_concat)
        return xy_concat


CNN = DeepSeqCNN().to(device)
optimizer = optim.Adam(CNN.parameters(), lr=lr)
#lossfunc = nn.L1Loss().to(device)
lossfunc = nn.CrossEntropyLoss(weight = torch.Tensor(w)).to(device)

def train_model(model, num_epochs):
    for epoch in range(num_epochs):
        # Training
        if epoch % 5 == 0:
            model.eval()
            test_predict = model(test_X1_sub, test_X2_sub)
            test_predict_np = test_predict.max(1).indices.detach().to('cpu').numpy()
            f2 = fbeta_score(test_log2FC_cat, test_predict_np, average='micro', beta=2)
            print('Epoch [%d] F2 score: %.3f' %
                    (epoch + 1, f2))
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
            if i % 200 == 199:    # print every 200 mini-batches
                print('[%d, %5d] loss: %.3f' %
                    (epoch + 1, i + 1, running_loss / 200))
                running_loss = 0.0
    
    return model

CNN = train_model(CNN, num_epochs=epochs)

ckptPATH = resultdir + '/models/gRNA_logFC_cat3-'+grp+str(thres)+'-CE-seqannot.pth'
torch.save(CNN.state_dict(), ckptPATH)

del test_X1_sub, test_X2_sub
test_X1 = torch.tensor(test_sequence_onehot, dtype=torch.float32).to(device)
test_X2 = torch.tensor(test_annotation, dtype=torch.float32).to(device)
CNN.eval()
test_predict = CNN(test_X1, test_X2)
test_predict_np = test_predict.detach().to('cpu').numpy()
test_predict_np_cat = test_predict.max(1).indices.detach().to('cpu').numpy()
fbeta_score(test_log2FC_cat, test_predict_np_cat, average='micro', beta=2)
PD = pd.DataFrame(np.concatenate([np.stack((test_log2FC_cat, test_predict_np_cat), axis=1), test_predict_np], axis=1), columns = ['true', 'predict', 'score_0','score_1','score_2'])
PD.to_csv(resultdir + '/gRNA_logFC_cat3-'+grp+str(thres)+'-CE-seqannot.csv')

