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

datadir = '/proj/yunligrp/users/tianyou/gRNA/data/data_Dec_resplit/'
resultdir = '/proj/yunligrp/users/tianyou/gRNA/result_resplit/models/'
batch_size = 256
epochs = 150 ## 200 for promoter, 150 for enhancer
lr = 0.0001
ngpu=1
grp = 'enh'
take_log = False

n_shap_sel = 1000
n_shap_ref = 500

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


if take_log:
    dat = pd.read_csv(datadir+'/wgCERES-gRNAs-k562-discovery-screen-'+grp+'_rawp0.05-binary-train-clean-log.csv', index_col = False)
else:
    dat = pd.read_csv(datadir+'/wgCERES-gRNAs-k562-discovery-screen-'+grp+'_rawp0.05-binary-train-clean.csv', index_col = False)

sequence = dat['protospacer']
sequence_onehot = preprocess_seq(sequence)
label = dat['significant'].to_numpy(dtype = np.float32)
class_count = dat['significant'].value_counts()
w = class_count[0] / class_count[1]
if grp == 'pro':
    annotation = dat.iloc[:,12:53].to_numpy(dtype = np.float32) # for promoters
elif grp == 'enh':
    annotation = dat.iloc[:,12:55].to_numpy(dtype = np.float32) # for enhancers
else:
    print("Invalid group: " + grp)


X1 = torch.tensor(sequence_onehot, dtype=torch.float32)
#Xloader = torch.utils.data.DataLoader(X, batch_size=batch_size, shuffle=True)
X2 = torch.tensor(annotation, dtype=torch.float32)
Y = torch.tensor(label, dtype=torch.float32)
Y = Y.view(-1, 1)
#Yloader = torch.utils.data.DataLoader(Y, batch_size=batch_size, shuffle=True)
input_dat = TensorDataset(X1,X2,Y)
datloader = DataLoader(input_dat, batch_size=batch_size, shuffle=True)


## test set
if take_log:
    test = pd.read_csv(datadir+'/wgCERES-gRNAs-k562-discovery-screen-'+grp+'_rawp0.05-binary-test-clean-log.csv', index_col = False)
else:
    test = pd.read_csv(datadir+'/wgCERES-gRNAs-k562-discovery-screen-'+grp+'_rawp0.05-binary-test-clean.csv', index_col = False)

test_sequence = test['protospacer']
test_sequence_onehot = preprocess_seq(test_sequence)
test_label = test['significant'].to_numpy(dtype = np.float32)
if grp == 'pro':
    test_annotation = test.iloc[:,12:53].to_numpy(dtype = np.float32) #promoters
elif grp == "enh":
    test_annotation = test.iloc[:,12:55].to_numpy(dtype = np.float32) #enhancers
else:
    print("Invalid group:" + grp)

subsample = np.random.choice(len(test_sequence), size = 3000, replace = False)
#test_X_sub = torch.tensor(test_sequence_onehot[subsample,:], dtype=torch.float32).to(device)
test_X1_sub = torch.tensor(test_sequence_onehot[subsample,:], dtype=torch.float32).to(device)
test_X2_sub = torch.tensor(test_annotation[subsample,:], dtype=torch.float32).to(device)

if grp == 'pro':
    dim_fc = 141
elif grp == 'enh':
    dim_fc = 143

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
            nn.Linear(dim_fc, 100),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(100, 80),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(80, 60),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(60, 1),
            nn.Sigmoid()
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
lossfunc = nn.BCEWithLogitsLoss(pos_weight=torch.Tensor([w])).to(device)

## load in the trained model
if take_log:
    ckptPATH=resultdir + '/gRNA_binary-log-'+grp+'-BCE-seqannot.pth'
else:
    ckptPATH=resultdir + '/gRNA_binary-'+grp+'-BCE-seqannot.pth'

CNN.load_state_dict(torch.load(ckptPATH, map_location = device))



rdm = np.random.choice(len(test_sequence), size = n_shap_sel, replace = False)
test_shap_rdm1 = torch.tensor(test_sequence_onehot[rdm,:], dtype=torch.float32)
test_shap_rdm2 = torch.tensor(test_annotation[rdm,:], dtype=torch.float32)

background_1 = test_shap_rdm1[:n_shap_ref]
background_2 = test_shap_rdm2[:n_shap_ref]

e = shap.DeepExplainer(CNN, [background_1, background_2])
shap_values = e.shap_values([test_shap_rdm1, test_shap_rdm2])

shap_seq = shap_values[0]
shap_annot = shap_values[1]

#shap.summary_plot(shap_annot, feature_names = dat.columns[12:55], show = False)
#
#if take_log:
#    plt.savefig(resultdir+'shap-log-seqannot'+grp+'.png', bbox_inches='tight')
#else:
#    plt.savefig(resultdir+'shap-seqannot'+grp+'.png', bbox_inches='tight')



#dinuc_shuff_explanations = torch.tensor(np.sum(shap_seq, axis=1)[:,None,:]) * test_images1

#from deeplift.visualization import viz_sequence
#for dinuc_shuff_explanation in dinuc_shuff_explanations:
#    viz_sequence.plot_weights(dinuc_shuff_explanation)



