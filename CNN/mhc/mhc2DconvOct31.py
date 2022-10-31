import torch
import argparse
import os, sys
import torchvision
import torchvision.transforms as transforms
import matplotlib.pyplot as plt
# import itertools
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
import pandas as pd
from sklearn.metrics import f1_score, roc_auc_score
import shap

datadir = '/Users/wancen/Downloads/data/'
resultdir = '/Users/wancen/Downloads/data/res'
batch_size = 256
epochs = 300 ## 60 for promoter, 15 for enhancer
lr = 0.0001
ngpu=1

device = torch.device('mps' if torch.has_mps else 'cpu')
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

np.random.seed(1)
torch.manual_seed(1)
torch.cuda.manual_seed_all(1)
torch.backends.cudnn.deterministic = True  # 保证每次结果一样

dat = pd.read_csv(datadir+'ipsc-pfdr0.05-pfdr0.2-binary-train.csv', index_col = False)
sequence = dat['grna']
sequence_onehot = preprocess_seq(sequence)
label = dat['significant'].to_numpy(dtype = np.float32)
class_count = dat['significant'].value_counts()
w = class_count[0] / class_count[1]
column_names = dat.columns
annotation = dat.iloc[:,np.r_[23,27:34]].to_numpy(dtype = np.float32)

X1 = torch.tensor(sequence_onehot, dtype=torch.float32)
#Xloader = torch.utils.data.DataLoader(X, batch_size=batch_size, shuffle=True)
X2 = torch.tensor(annotation, dtype=torch.float32)
Y = torch.tensor(label, dtype=torch.float32)
Y = Y.view(-1, 1)
#Yloader = torch.utils.data.DataLoader(Y, batch_size=batch_size, shuffle=True)
input_dat = TensorDataset(X1,X2,Y)
datloader = DataLoader(input_dat, batch_size=batch_size, shuffle=True, num_workers=0)

## return a random batch
# dataiter = iter(datloader)
# seq, anno, labels = next(dataiter)
# print(seq.shape) # 100, 4, 20
# seq = seq.unsqueeze(1)
# print(seq.shape)
# conv0 = nn.Conv2d(1,50,(4,2),padding=(0,1)) ## kernel size with 2
# conv1 = nn.Conv2d(1,100,(4,3),padding=(0,1))
# conv2 = nn.Conv2d(1,70,(4,5),padding=(0,2))
# conv3 = nn.Conv2d(1,40,(4,7),padding=(0,3))
# def conv_and_pool(x, conv):
#     x = F.relu(conv(x)).squeeze(2) # 512, 50, 21
#     x = F.max_pool1d(x, x.size(2)).squeeze(2)
#     x = F.dropout(x,0.4)
#     return x
# x0 = conv_and_pool(seq,conv0)
# print(x0.shape)
# x1 = conv_and_pool(seq,conv1)
# print(x1.shape)
# x2 = conv_and_pool(seq,conv2)
# print(x2.shape) # 100, 70, 10
# x3 = conv_and_pool(seq,conv3)
# print(x3.shape) # 100, 40, 10
# # x_concat = torch.cat((x0, x1, x2, x3), dim=1) #100, 260,10
# x_concat = torch.cat((x0,x1,x2,x3), dim=1) # 512 260
# print(x_concat.shape)
# fc1 = nn.Sequential(
#             nn.Linear(260, 32),
#             nn.ReLU(),
#             nn.Dropout(0.3))
# # x_concat = x_concat.view(-1, 260 * 10)
# x_concat = fc1(x_concat)
# print(x_concat.shape) # 512 32
# conv4 = nn.Sequential(
#             nn.Conv2d(1, 1, (4,1)),
#             nn.ReLU(),
#         )
# x4 = conv4(seq).view(-1, 20) # 512, 20
# xy_concat = torch.cat((x_concat, x4, anno), dim=1)
# print(xy_concat.shape) # 512, 60
## test set
test = pd.read_csv(datadir+'ipsc-pfdr0.05-pfdr0.2-binary-test.csv', index_col = False)
test_sequence = test['grna']
test_sequence_onehot = preprocess_seq(test_sequence)
test_label = test['significant'].to_numpy(dtype = np.float32)
test_annotation = test.iloc[:,np.r_[23,27:34]].to_numpy(dtype = np.float32)

test_X1_sub = torch.tensor(test_sequence_onehot, dtype=torch.float32).to(device)
test_X2_sub = torch.tensor(test_annotation, dtype=torch.float32).to(device)

dim_fc = 60
class DeepSeqCNN(nn.Module):
    def __init__(self):
        super(DeepSeqCNN, self).__init__()
        self.conv0 = nn.Conv2d(1,50,(4,2),padding=(0,1))
        self.conv1 = nn.Conv2d(1,100,(4,3),padding=(0,1))
        self.conv2 = nn.Conv2d(1,70,(4,5),padding=(0,1))
        self.conv3 = nn.Conv2d(1,40,(4,7),padding=(0,1))
        self.conv4 = nn.Sequential(
            nn.Conv2d(1, 1, (4,1)),
            nn.ReLU(),
        )
        self.fc1 = nn.Sequential(
            nn.Linear(260, 32),
            nn.ReLU(),
            nn.Dropout(0.3))
        self.fc2 = nn.Sequential(
            nn.Linear(dim_fc, 16),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(16, 1),
            # nn.Sigmoid()  ## BCEWithLogitsLoss takes in logits directly without sigmoid
        )

    def conv_and_pool(self, x, conv):
        x = F.relu(conv(x)).squeeze(2)  # 512, 50, 21
        x = F.max_pool1d(x, x.size(2)).squeeze(2)
        x = F.dropout(x, 0.4)
        return x
    def forward(self, x, y):
        x = x.unsqueeze(1)
        x0 = self.conv_and_pool(x,self.conv0)
        x1 = self.conv_and_pool(x,self.conv1)
        x2 = self.conv_and_pool(x,self.conv2)
        x3 = self.conv_and_pool(x,self.conv3)
        x4 = self.conv4(x).view(-1, 20)
        # x_concat = torch.cat(
        #     (x0.view(batch_size, -1), x1.view(batch_size, -1), x2.view(batch_size, -1), x3.view(batch_size, -1)), dim=1)  # 100, 2350
        x_concat = torch.cat((x0, x1, x2, x3), dim=1)  # size: [:,260,10]
        # x_concat = x_concat.view(-1, 260 * 10)
        x_concat = self.fc1(x_concat)
        xy_concat = torch.cat((x_concat, x4, y), dim=1)
        xy_concat = self.fc2(xy_concat)
        # for layer in self.fc:
        #    x_concat = layer(x_concat)
        #    print(x_concat.size())
        # x_concat = self.fc(x_concat)
        return xy_concat

CNN = DeepSeqCNN().to(device)
optimizer = optim.Adam(CNN.parameters(), lr=lr)
scheduler = optim.lr_scheduler.StepLR(optimizer, step_size = 3, gamma =0.1)
# lossfunc = nn.L1Loss().to(device)
lossfunc = nn.BCEWithLogitsLoss(pos_weight=torch.Tensor([w])).to(device)
sigmoid = nn.Sigmoid()


def train_model(model, num_epochs):
    for epoch in range(num_epochs):
        # Training
        if epoch % 2 == 0:
            model.eval()
            test_predict = sigmoid(model(test_X1_sub, test_X2_sub))
            test_predict_np = test_predict.detach().to('cpu').numpy()
            auc = roc_auc_score(test_label, test_predict_np)
            print('Epoch [%d] AUC: %.3f' %
                  (epoch + 1, auc))
            model.train()
        # running_loss = 0.0
        for i, batch in enumerate(datloader, 0):
            # Transfer to GPU
            local_x1, local_x2, local_y = batch
            local_x1, local_x2, local_y = local_x1.to(device), local_x2.to(device), local_y.to(device)
            optimizer.zero_grad()
            FC_pred = model(local_x1, local_x2)
            loss = lossfunc(FC_pred, local_y)
            loss.backward()
            optimizer.step()
            # running_loss += loss.item()
            if i % 10 == 0:  # print every 200 mini-batches
                print('[%d, %5d] loss: %.3f' %
                      (epoch + 1, i + 1, loss.item()))
                # running_loss = 0.0
        # scheduler.step()
    return model

CNN = train_model(CNN, num_epochs=epochs)

ckptPATH = resultdir + '/models/iPSC-binary-BCE-seqannot-Oct30.pth'
torch.save(CNN.state_dict(), ckptPATH)
CNN = torch.load(ckptPATH)
del test_X1_sub, test_X2_sub
test_X1 = torch.tensor(test_sequence_onehot, dtype=torch.float32).to(device)
test_X2 = torch.tensor(test_annotation, dtype=torch.float32).to(device)
CNN.eval()
test_predict = sigmoid(CNN(test_X1, test_X2))
test_predict_np = test_predict.detach().to('cpu').numpy()
roc_auc_score(test_label, test_predict_np)
PD = pd.DataFrame(np.stack((test_label, test_predict_np[:,0]), axis=1), columns = ['true', 'predict'])
PD.to_csv(resultdir + '/iPSC-binary-BCE-seqannot.csv')