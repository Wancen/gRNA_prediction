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


#parser = argparse.ArgumentParser(description='gRNA prediction model')
#parser.add_argument('--savedir', default='./', help='path to save results')
#parser.add_argument('--ckptdir', default='./ckpt', help='path to save checkpoints')
#parser.add_argument('--batch-size', type=int, default=128,
#                    help='input batch size for training (default: 128)')
#parser.add_argument('--epochs', type=int, default=100,
#                    help='number of epochs to train (default: 100)')
#parser.add_argument('--lr', type=float, default=0.001,
#                    help='learning rate (default: 0.001)')
#args = parser.parse_args()
#
#savedir = args.savedir
#ckptdir = args.ckptdir


#batch_size = args.batch_size
#epochs = args.epochs
#lr = args.lr
#ngpu=1

datadir = '/proj/yunligrp/users/tianyou/gRNA/data/data_April_resplit/'
resultdir = '/proj/yunligrp/users/tianyou/gRNA/result_April_resplit/binary/'
batch_size = 256
epochs = 15 ## 50 for promoter, 150 for enhancer
lr = 0.0001
ngpu=1
grp = 'enh'

device = torch.device("cuda:0" if (torch.cuda.is_available() and ngpu > 0) else "cpu")
print(device)


dat = pd.read_csv(datadir+'wgCERES-gRNAs-k562-discovery-screen-'+grp+'_baseMean125-binary-train-clean.csv', index_col = False)
label = dat['significance'].to_numpy(dtype = np.float32)
class_count = dat['significance'].value_counts()
w = class_count[0] / class_count[1]
if grp == 'pro':
    annotation = dat.iloc[:,12:54].to_numpy(dtype = np.float32) # for promoters
elif grp == 'enh':
    annotation = dat.iloc[:,12:56].to_numpy(dtype = np.float32) # for enhancers
else:
    print("Invalid group: " + grp)


X2 = torch.tensor(annotation, dtype=torch.float32)
Y = torch.tensor(label, dtype=torch.float32)
Y = Y.view(-1, 1)
#Yloader = torch.utils.data.DataLoader(Y, batch_size=batch_size, shuffle=True)
input_dat = TensorDataset(X2,Y)
datloader = DataLoader(input_dat, batch_size=batch_size, shuffle=True)


## test set
test = pd.read_csv(datadir+'/wgCERES-gRNAs-k562-discovery-screen-'+grp+'_baseMean125-binary-test-clean.csv', index_col = False)
test_label = test['significance'].to_numpy(dtype = np.float32)
if grp == 'pro':
    test_annotation = test.iloc[:,12:54].to_numpy(dtype = np.float32) #promoters
elif grp == "enh":
    test_annotation = test.iloc[:,12:56].to_numpy(dtype = np.float32) #enhancers
else:
    print("Invalid group:" + grp)

subsample = np.random.choice(len(test_annotation), size = 4000, replace = False)
#test_X_sub = torch.tensor(test_sequence_onehot[subsample,:], dtype=torch.float32).to(device)
test_X2_sub = torch.tensor(test_annotation[subsample,:], dtype=torch.float32).to(device)

if grp == 'pro':
    dim_fc = 42
elif grp == 'enh':
    dim_fc = 44

class DeepSeqCNN(nn.Module):
    def __init__(self):
        super(DeepSeqCNN, self).__init__()
        self.fc2 = nn.Sequential(
            nn.Linear(dim_fc, 64),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(32, 8),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(8, 1),
            #nn.Sigmoid()
        )
        
    def forward(self, y):
        y_concat = self.fc2(y)
        #for layer in self.fc:
        #    x_concat = layer(x_concat)
        #    print(x_concat.size())
        #x_concat = self.fc(x_concat)
        return y_concat


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
            test_predict = sigmoid(model(test_X2_sub))
            test_predict_np = test_predict.detach().to('cpu').numpy()
            auc = roc_auc_score(test_label[subsample], test_predict_np)
            print('Epoch [%d] AUC: %.3f' %
                    (epoch + 1, auc))
            model.train()
        running_loss = 0.0
        for i, batch in enumerate(datloader, 0):
            # Transfer to GPU
            local_x2, local_y = batch
            local_x2, local_y = local_x2.to(device), local_y.to(device)
            optimizer.zero_grad()
            FC_pred = model(local_x2)
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

ckptPATH = resultdir + '/models/gRNA_binary-'+grp+'-BCE-annot-Oct04.pth'
torch.save(CNN.state_dict(), ckptPATH)

del test_X2_sub
test_X2 = torch.tensor(test_annotation, dtype=torch.float32).to(device)
CNN.eval()
test_predict = sigmoid(CNN(test_X2))
test_predict_np = test_predict.detach().to('cpu').numpy()
roc_auc_score(test_label, test_predict_np)
PD = pd.DataFrame(np.stack((test_label, test_predict_np[:,0]), axis=1), columns = ['true', 'predict'])
PD.to_csv(resultdir + '/gRNA_binary-'+grp+'-BCE-annot.csv')

