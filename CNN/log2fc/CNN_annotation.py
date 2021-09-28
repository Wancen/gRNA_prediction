import torch
import argparse
import os
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

parser = argparse.ArgumentParser(description='gRNA prediction model')
parser.add_argument('--savedir', default='./', help='path to save results')
parser.add_argument('--ckptdir', default='./ckpt', help='path to save checkpoints')
parser.add_argument('--batch-size', type=int, default=128,
                    help='input batch size for training (default: 128)')
parser.add_argument('--epochs', type=int, default=100,
                    help='number of epochs to train (default: 100)')
parser.add_argument('--lr', type=float, default=0.001,
                    help='learning rate (default: 0.001)')
args = parser.parse_args()

savedir = args.savedir
ckptdir = args.ckptdir


#batch_size = args.batch_size
#epochs = args.epochs
#lr = args.lr
#ngpu=1


batch_size = 256
epochs = 130
lr = 0.0001
ngpu=1

device = torch.device("cuda:0" if (torch.cuda.is_available() and ngpu > 0) else "cpu")
print(device)


def preprocess_seq(data):
    print("Start preprocessing the sequence")
    length = len(data[0])
    DATA_X = np.zeros((len(data),4,length), dtype=int)
    print(DATA_X.shape)
    for l in range(len(data)):
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


dat = pd.read_csv('/pine/scr/t/i/tianyou/Patrick/data/wgCERES-gRNAs-k562-discovery-screen-enh_rawp0.05-padj0.2-train-clean.csv', index_col = False)
sequence = dat['protospacer']
sequence_onehot = preprocess_seq(sequence)
log2FC = dat['log2FoldChange'].to_numpy(dtype = np.float32)
#log2FC = np.abs(log2FC)
#annotation = dat.iloc[:,13:46].to_numpy(dtype = np.float32) # for promoters
annotation = dat.iloc[:,13:48].to_numpy(dtype = np.float32) # for enhancers

X1 = torch.tensor(sequence_onehot, dtype=torch.float32)
#Xloader = torch.utils.data.DataLoader(X, batch_size=batch_size, shuffle=True)
X2 = torch.tensor(annotation, dtype=torch.float32)
Y = torch.tensor(log2FC, dtype=torch.float32)
Y = Y.view(-1, 1)
#Yloader = torch.utils.data.DataLoader(Y, batch_size=batch_size, shuffle=True)
input_dat = TensorDataset(X1,X2,Y)
datloader = DataLoader(input_dat, batch_size=batch_size, shuffle=True)


## test set
test = pd.read_csv('/pine/scr/t/i/tianyou/Patrick/data/wgCERES-gRNAs-k562-discovery-screen-enh_rawp0.05-padj0.2-test-clean.csv', index_col = False)
test_sequence = test['protospacer']
test_sequence_onehot = preprocess_seq(test_sequence)
test_log2FC = test['log2FoldChange'].to_numpy(dtype = np.float32)
#test_log2FC = np.abs(test_log2FC)
#test_annotation = test.iloc[:,13:46].to_numpy(dtype = np.float32) #promoters
test_annotation = test.iloc[:,13:48].to_numpy(dtype = np.float32) #enhancers
#subsample = np.random.choice(len(test_sequence), size = 3000, replace = False)
#test_X_sub = torch.tensor(test_sequence_onehot[subsample,:], dtype=torch.float32).to(device)
test_X1_sub = torch.tensor(test_sequence_onehot, dtype=torch.float32).to(device)
test_X2_sub = torch.tensor(test_annotation, dtype=torch.float32).to(device)


class DeepSeqCNN(nn.Module):
    def __init__(self):
        super(DeepSeqCNN, self).__init__()
        #self.input_size = input_size
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
            nn.Linear(210*10, 100),
            nn.ReLU(),
            nn.Dropout(0.3))
        self.fc2 = nn.Sequential(
            #nn.Linear(133, 100), #promoter
            nn.Linear(135, 100),  #enhancer
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(100, 80),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(80, 60),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(60, 1)
        )
        
    def forward(self, x, y):
        x1 = self.conv1(x)
        x2 = self.conv2(x)
        x3 = self.conv3(x)
        x_concat = torch.cat((x1, x2, x3), dim=1) # size: [:,210,10]
        x_concat = x_concat.view(-1, 210*10)
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
lossfunc = nn.MSELoss().to(device)

for epoch in range(epochs):
    # Training
    if epoch % 5 == 0:
        CNN.eval()
        test_predict = CNN(test_X1_sub, test_X2_sub)
        test_predict_np = test_predict.detach().to('cpu').numpy()
        #r2 = spearmanr(test_log2FC[subsample], test_predict_np)
        r2 = spearmanr(test_log2FC, test_predict_np)
        print('[%d] spearman correlation: %.3f' %
                  (epoch + 1, r2[0]))
        CNN.train()
    
    running_loss = 0.0
    for i, batch in enumerate(datloader, 0):
        # Transfer to GPU
        local_x1, local_x2, local_y = batch
        local_x1, local_x2, local_y = local_x1.to(device), local_x2.to(device), local_y.to(device)
        optimizer.zero_grad()
        FC_pred = CNN(local_x1, local_x2)
        loss = lossfunc(FC_pred, local_y)
        loss.backward()
        optimizer.step()
        
        running_loss += loss.item()
        if i % 25 == 24:    # print every 20 mini-batches
            print('[%d, %5d] loss: %.3f' %
                  (epoch + 1, i + 1, running_loss / 25))
            running_loss = 0.0


del test_X1_sub, test_X2_sub
test_X1 = torch.tensor(test_sequence_onehot, dtype=torch.float32).to(device)
test_X2 = torch.tensor(test_annotation, dtype=torch.float32).to(device)
CNN.eval()
test_predict = CNN(test_X1, test_X2)
test_predict_np = test_predict.detach().to('cpu').numpy()
spearmanr(test_log2FC, test_predict_np)
PD = pd.DataFrame(np.stack((test_log2FC, test_predict_np[:,0]), axis=1), columns = ['true', 'predict'])
PD.to_csv("./result_deltaGB/gRNA_predict-pro-padj0.2-MSE.csv")


## evaluate on the other test set: promoter on enhancer
test_oth1 = pd.read_csv('/pine/scr/t/i/tianyou/Patrick/data/wgCERES-gRNAs-k562-discovery-screen-enh_rawp0.05-padj0.2-test-clean.csv', index_col = False)
test_oth2 = pd.read_csv('/pine/scr/t/i/tianyou/Patrick/data/wgCERES-gRNAs-k562-discovery-screen-enh_rawp0.05-padj0.2-train-clean.csv', index_col = False)
test_oth = pd.concat([test_oth1,test_oth2], ignore_index=True)
test_oth_sequence = test_oth['protospacer']
test_oth_sequence_onehot = preprocess_seq(test_oth_sequence)
test_oth_log2FC = test_oth['log2FoldChange'].to_numpy(dtype = np.float32)
#test_log2FC = np.abs(test_log2FC)
test_oth_annotation = test_oth.iloc[:,np.r_[13:42,46,47,42,43]].to_numpy(dtype = np.float32) #for promoters, make the enhancer test file the same format
test_oth_X1 = torch.tensor(test_oth_sequence_onehot, dtype=torch.float32).to(device)
test_oth_X2 = torch.tensor(test_oth_annotation, dtype=torch.float32).to(device)
CNN.eval()
test_oth_predict = CNN(test_oth_X1, test_oth_X2)
test_oth_predict_np = test_oth_predict.detach().to('cpu').numpy()
spearmanr(test_oth_log2FC, test_oth_predict_np)
PD_oth = pd.DataFrame(np.stack((test_oth_log2FC, test_oth_predict_np[:,0]), axis=1), columns = ['true', 'predict'])
PD_oth.to_csv("./result_deltaGB/gRNA_predict-pro-on-enh-padj0.2-MSE.csv")


## evaluate on the other test set: enhancer on promoter
test_oth1 = pd.read_csv('/pine/scr/t/i/tianyou/Patrick/data/wgCERES-gRNAs-k562-discovery-screen-pro_rawp0.05-padj0.2-test-clean.csv', index_col = False)
test_oth2 = pd.read_csv('/pine/scr/t/i/tianyou/Patrick/data/wgCERES-gRNAs-k562-discovery-screen-pro_rawp0.05-padj0.2-train-clean.csv', index_col = False)
test_oth = pd.concat([test_oth1,test_oth2], ignore_index=True)
test_oth['promnumber']=np.mean(test['promnumber'])
test_oth['promlog10fdr']=np.mean(test['promlog10fdr'])
test_oth_sequence = test_oth['protospacer']
test_oth_sequence_onehot = preprocess_seq(test_oth_sequence)
test_oth_log2FC = test_oth['log2FoldChange'].to_numpy(dtype = np.float32)
#test_log2FC = np.abs(test_log2FC)
test_oth_annotation = test_oth.iloc[:,np.r_[13:42,44,45,47,48,42,43]].to_numpy(dtype = np.float32) #for enhancers, make the promoter test file the same format
test_oth_X1 = torch.tensor(test_oth_sequence_onehot, dtype=torch.float32).to(device)
test_oth_X2 = torch.tensor(test_oth_annotation, dtype=torch.float32).to(device)
CNN.eval()
test_oth_predict = CNN(test_oth_X1, test_oth_X2)
test_oth_predict_np = test_oth_predict.detach().to('cpu').numpy()
spearmanr(test_oth_log2FC, test_oth_predict_np)
PD_oth = pd.DataFrame(np.stack((test_oth_log2FC, test_oth_predict_np[:,0]), axis=1), columns = ['true', 'predict'])
PD_oth.to_csv("./result_deltaGB/gRNA_predict-enh-on-pro-padj0.2-MSE.csv")