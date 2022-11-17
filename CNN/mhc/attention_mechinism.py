from numpy import array
from numpy import random
from numpy import dot
from scipy.special import softmax
 
# encoder representations of four different words
word_1 = array([1, 0, 0])
word_2 = array([0, 1, 0])
word_3 = array([1, 1, 0])
word_4 = array([0, 0, 1])
 
# stacking the word embeddings into a single array
words = array([word_1, word_2, word_3, word_4])
 
# generating the weight matrices
random.seed(42)
W_Q = random.randint(3, size=(3, 3))
W_K = random.randint(3, size=(3, 3))
W_V = random.randint(3, size=(3, 3))
 
# generating the queries, keys and values
Q = words @ W_Q
K = words @ W_K
V = words @ W_V
 
# scoring the query vectors against all key vectors
scores = Q @ K.transpose()
 
# computing the weights by a softmax operation
weights = softmax(scores / K.shape[1] ** 0.5, axis=1)
 
# computing the attention by a weighted sum of the value vectors
attention = weights @ V
 
print(attention)

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
batch_size = 64
epochs = 100 ## 60 for promoter, 15 for enhancer
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

def preprocess_seq(data):
    print("Start preprocessing the sequence")
    length = len(data[0])
    DATA_X = np.zeros((len(data),length), dtype=int)
    print(DATA_X.shape)
    for l in range(len(data)):
        if l % 10000 == 0:
            print(str(l) + " sequences processed")
        for i in range(length):
            try: data[l][i]
            except: print(data[l], i, length, len(data))
            if data[l][i]in "Aa":
                DATA_X[l, i] = 0
            elif data[l][i] in "Cc":
                DATA_X[l, i] = 1
            elif data[l][i] in "Gg":
                DATA_X[l, i] = 2
            elif data[l][i] in "Tt":
                DATA_X[l, i] = 3
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
from sklearn import preprocessing
min_max_scaler = preprocessing.MinMaxScaler()
annotation_scaled = min_max_scaler.fit_transform(annotation)

X1 = torch.LongTensor(sequence_onehot)
#Xloader = torch.utils.data.DataLoader(X, batch_size=batch_size, shuffle=True)
X2 = torch.tensor(annotation, dtype=torch.float32)
Y = torch.tensor(label, dtype=torch.float32)
Y = Y.view(-1, 1)

test = pd.read_csv(datadir+'ipsc-pfdr0.05-pfdr0.2-binary-test.csv', index_col = False)
test_sequence = test['grna']
test_sequence_onehot = preprocess_seq(test_sequence)
test_label = test['significant'].to_numpy(dtype = np.float32)
test_annotation = test.iloc[:,np.r_[23,27:34]].to_numpy(dtype = np.float32)
test_annotation_scaled = min_max_scaler.fit_transform(test_annotation)

test_X1 = torch.LongTensor(test_sequence_onehot).to(device)
test_X2 = torch.tensor(test_annotation, dtype=torch.float32).to(device)

## return a random batch
dataiter = iter(datloader)
seq, anno, labels = next(dataiter)
print(seq.shape) # 100, 4, 20
seq = seq.unsqueeze(1)


import numpy as np
import torch
from torch import nn

class SH_SelfAttention(nn.Module):
    """ single head self-attention module
    """
    def __init__(self, input_size):
        
        super().__init__()
        # define query, key and value transformation matrices
        # usually input_size is equal to embed_size
        self.embed_size = input_size
        self.Wq = nn.Linear(input_size, self.embed_size, bias=False)
        self.Wk = nn.Linear(input_size, self.embed_size, bias=False)
        self.Wv = nn.Linear(input_size, self.embed_size, bias=False)
        self.softmax = nn.Softmax(dim=2) # normalized across feature dimension
    
    def forward(self, X):
        """
        Args:
            X: tensor, (batch, sequence length, input_size)
        """
        X_q = self.Wq(X) # queries
        X_k = self.Wk(X) # keys
        X_v = self.Wv(X) # values
        
        # scaled queries and keys by forth root 
        X_q_scaled = X_q / (self.embed_size ** (1/4))
        X_k_scaled = X_k / (self.embed_size ** (1/4))
        
        attn_w = torch.bmm(X_q_scaled, X_k_scaled.transpose(1,2))
        # (batch, sequence length, sequence length)
        attn_w_normalized = self.softmax(attn_w)
        # print('attn_w_normalized.shape', attn_w_normalized.shape)
        
        # reweighted value vectors
        z = torch.bmm(attn_w_normalized, X_v)
        
        return z, attn_w_normalized
    

class MH_SelfAttention(nn.Module):
    """ multi head self-attention module
    """
    def __init__(self, input_size, num_attn_heads):
        
        super().__init__()
        
        layers = [SH_SelfAttention(input_size) for i in range(num_attn_heads)]
        
        self.multihead_pipeline = nn.ModuleList(layers)
        embed_size = input_size
        self.Wz = nn.Linear(num_attn_heads*embed_size, embed_size)
        
    
    
    def forward(self, X):
        """
        Args:
            X: tensor, (batch, sequence length, input_size)
        """
        
        out = []
        bsize, num_positions, inp_dim = X.shape
        attn_tensor = X.new_zeros((bsize, num_positions, num_positions))
        for SH_layer in self.multihead_pipeline:
            z, attn_w = SH_layer(X)
            out.append(z)
            attn_tensor += attn_w
        # concat on the feature dimension
        out = torch.cat(out, -1) 
        attn_tensor = attn_tensor/len(self.multihead_pipeline)

        # return a unified vector mapping of the different self-attention blocks
        return self.Wz(out), attn_tensor
        

class TransformerUnit(nn.Module):
    
    def __init__(self, input_size, num_attn_heads, mlp_embed_factor, nonlin_func, pdropout):
        
        super().__init__()
        
        embed_size = input_size
        self.multihead_attn = MH_SelfAttention(input_size, num_attn_heads)
        
        self.layernorm_1 = nn.LayerNorm(embed_size)

        # also known as position wise feed forward neural network
        self.MLP = nn.Sequential(
            nn.Linear(embed_size, embed_size*mlp_embed_factor),
            nonlin_func,
            nn.Linear(embed_size*mlp_embed_factor, embed_size)
        )
        
        self.layernorm_2 = nn.LayerNorm(embed_size)
        
        self.dropout = nn.Dropout(p=pdropout)
                
    
    def forward(self, X):
        """
        Args:
            X: tensor, (batch, sequence length, input_size)
        """
        # z is tensor of size (batch, sequence length, input_size)
        z, attn_mhead_tensor = self.multihead_attn(X)
        # layer norm with residual connection
        z = self.layernorm_1(z + X)
        z = self.dropout(z)
        z_ff= self.MLP(z)
        z = self.layernorm_2(z_ff + z)
        z = self.dropout(z)
        
        return z, attn_mhead_tensor

"""
      implement position encoder based on cosine and sine approach proposed 
      by original Transformers paper ('Attention is all what you need')
"""
class NucleoPosEncoding(nn.Module):
    def __init__(self, num_nucleotides, seq_len, embed_dim, pdropout=0.1):
        super().__init__()
        self.nucleo_emb = nn.Embedding(num_nucleotides, embed_dim)
        self.dropout = nn.Dropout(p=pdropout)
        # positional encoding matrix
        base_pow = 10000
        PE_matrix = torch.zeros((1, seq_len, embed_dim))
        i_num = torch.arange(0., seq_len).reshape(-1, 1) # i iterates over sequence length (i.e. sequence items)
        j_denom = torch.pow(base_pow, torch.arange(0., embed_dim, 2.) / embed_dim) # j iterates over embedding dimension
        PE_matrix[:, :, 0::2] = torch.sin(i_num/j_denom)
        PE_matrix[:, :, 1::2] = torch.cos(i_num/j_denom)
        self.register_buffer('PE_matrix', PE_matrix)
        
        
    def forward(self, X):
        """
        Args:
            X: tensor, int64,  (batch, sequence length)
        """
        X_emb = self.nucleo_emb(X)
        # (batch, sequence length, embedding dim)
        X_embpos = X_emb + self.PE_matrix
        return self.dropout(X_embpos)

class NucleoPosEmbedder(nn.Module):
    def __init__(self, num_nucleotides, seq_length, embedding_dim):
        super().__init__()
        self.nucleo_emb = nn.Embedding(num_nucleotides, embedding_dim)
        self.pos_emb = nn.Embedding(seq_length, embedding_dim)

    def forward(self, X):
        """
        Args:
            X: tensor, int64,  (batch, sequence length)
        """
        X_emb = self.nucleo_emb(X)
        bsize, seqlen, featdim = X_emb.size()
        device = X_emb.device
        positions = torch.arange(seqlen).to(device)
        positions_emb = self.pos_emb(positions)[None, :, :].expand(bsize, seqlen, featdim)
        # (batch, sequence length, embedding dim)
        X_embpos = X_emb + positions_emb
        return X_embpos

class PerBaseFeatureEmbAttention(nn.Module):
    """ Per base feature attention module
    """
    def __init__(self, input_dim, seq_len):
        
        super().__init__()
        # define query, key and value transformation matrices
        # usually input_size is equal to embed_size
        self.embed_size = input_dim
        self.Q = nn.Parameter(torch.randn((seq_len, input_dim), dtype=torch.float32), requires_grad=True)
        self.softmax = nn.Softmax(dim=-1) # normalized across feature dimension
    
    def forward(self, X):
        """
        Args:
            X: tensor, (batch, sequence length, input_size)
        """
        bsize, seqlen, featdim = X.shape
        X_q = self.Q[None, :, :].expand(bsize, seqlen, featdim) # queries
        X_k = X
        X_v = X
        # scaled queries and keys by forth root 
        X_q_scaled = X_q / (self.embed_size ** (1/4))
        X_k_scaled = X_k / (self.embed_size ** (1/4))
        # print(X_q_scaled.shape)
        # print(X_k_scaled.shape)
        
        attn_w = torch.bmm(X_q_scaled, X_k_scaled.transpose(1,2))
        # attn_w = X_q_scaled.matmul(X_k_scaled.transpose(1,0))
        # (batch, sequence length, sequence length)
        attn_w_normalized = self.softmax(attn_w)
        # print('attn_w_normalized.shape', attn_w_normalized.shape)
        
        # reweighted value vectors
        z = torch.bmm(attn_w_normalized, X_v)
        # print('z.shape', z.shape)
        
        return z, attn_w_normalized

class FeatureEmbAttention(nn.Module):
    def __init__(self, input_dim):
        '''
        Args:
            input_dim: int, size of the input vector (i.e. feature vector)
        '''

        super().__init__()
        self.input_dim = input_dim
        # use this as query vector against the transformer outputs
        self.queryv = nn.Parameter(torch.randn(input_dim, dtype=torch.float32), requires_grad=True)
        self.softmax = nn.Softmax(dim=1) # normalized across seqlen

    def forward(self, X):
        '''Performs forward computation
        Args:
            X: torch.Tensor, (bsize, seqlen, feature_dim), dtype=torch.float32
        '''

        X_scaled = X / (self.input_dim ** (1/4))
        queryv_scaled = self.queryv / (self.input_dim ** (1/4))
        # using  matmul to compute tensor vector multiplication
        
        # (bsize, seqlen)
        attn_weights = X_scaled.matmul(queryv_scaled)

        # softmax
        attn_weights_norm = self.softmax(attn_weights)

        # reweighted value vectors (in this case reweighting the original input X)
        # unsqueeze attn_weights_norm to get (bsize, 1, seqlen)
        # perform batch multiplication with X that has shape (bsize, seqlen, feat_dim)
        # result will be (bsize, 1, feat_dim)
        # squeeze the result to obtain (bsize, feat_dim)
        z = attn_weights_norm.unsqueeze(1).bmm(X).squeeze(1)
        
        # returns (bsize, feat_dim), (bsize, seqlen)
        return z, attn_weights_norm


class Categ_CrisCasTransformer(nn.Module):

    def __init__(self, input_size=64, num_nucleotides=4, 
                 seq_length=20, num_attn_heads=8, 
                 mlp_embed_factor=2, nonlin_func=nn.ReLU(), 
                 pdropout=0.3, num_transformer_units=12, 
                 pooling_mode='attn', num_classes=1, per_base=False):
        
        super().__init__()
        
        embed_size = input_size

        self.nucleopos_embedder = NucleoPosEmbedder(num_nucleotides, seq_length, embed_size)
        # self.nucleopos_embedder = NucleoPosEncoding(num_nucleotides, seq_length, embed_size)
        
        trfunit_layers = [TransformerUnit(input_size, num_attn_heads, mlp_embed_factor, nonlin_func, pdropout) 
                          for i in range(num_transformer_units)]
        # self.trfunit_layers = trfunit_layers
        self.trfunit_pipeline = nn.ModuleList(trfunit_layers)
        # self.trfunit_pipeline = nn.Sequential(*trfunit_layers)
        self.per_base = per_base
        
        if not per_base:
            self.pooling_mode = pooling_mode
            if pooling_mode == 'attn':
                self.pooling = FeatureEmbAttention(input_size)
            elif pooling_mode == 'mean':
                self.pooling = torch.mean
            self.Wy = nn.Linear(embed_size, num_classes, bias=True)

        else:
            self.pooling_mode = pooling_mode
            if pooling_mode == 'attn':
                self.pooling = PerBaseFeatureEmbAttention(input_size, seq_length)
            self.bias = nn.Parameter(torch.randn((seq_length, num_classes), dtype=torch.float32), requires_grad=True)
            self.Wy = nn.Linear(embed_size, num_classes, bias=False)

        # perform log softmax on the feature dimension
        self.log_softmax = nn.LogSoftmax(dim=-1)
        self._init_params_()
        
    def _init_params_(self):
        for p_name, p in self.named_parameters():
            param_dim = p.dim()
            if param_dim > 1: # weight matrices
                nn.init.xavier_uniform_(p)
            elif param_dim == 1: # bias parameters
                if p_name.endswith('bias'):
                    nn.init.uniform_(p, a=-1.0, b=1.0)
                    # nn.init.xavier_uniform_(p)

    def forward(self, X):
        """
        Args:
            X: tensor, int64,  (batch, sequence length)
        """
        # (batch, seqlen, embedding dim)
        X_embpos = self.nucleopos_embedder(X)
        # z is tensor of size (batch,  seqlen, embedding dim)
        bsize, num_positions, inp_dim = X_embpos.shape
        attn_tensor = X_embpos.new_zeros((bsize, num_positions, num_positions))
        xinput = X_embpos
        for trfunit in self.trfunit_pipeline:
            z, attn_mhead_tensor = trfunit(xinput)
            xinput = z
            attn_tensor += attn_mhead_tensor
        attn_tensor = attn_tensor/len(self.trfunit_pipeline)

         # pool across seqlen vectors
        if not self.per_base:
            if self.pooling_mode == 'attn':
                z, fattn_w_norm = self.pooling(z)
            # Note: z.mean(dim=1) or self.pooling(z, dim=1) will change shape of z to become (batch, embedding dim)
            # we can keep dimension by running z.mean(dim=1, keepdim=True) to have (batch, 1, embedding dim)
            elif self.pooling_mode == 'mean':
                z = self.pooling(z, dim=1)
                fattn_w_norm = None
            y = self.Wy(z)
        else:
            if self.pooling_mode == 'attn':
                z, fattn_w_norm = self.pooling(z)
            y = self.Wy(z) + self.bias
        
        return self.log_softmax(y), fattn_w_norm, attn_tensor

X = torch.LongTensor(random.randint(4, size=(10000, 20)))
Y = torch.FloatTensor(random.randint(2, size = (10000,1)))
input_dat = TensorDataset(X1,Y)
datloader = DataLoader(input_dat, batch_size=64, shuffle=True, num_workers=0)

def run_prediction(model, dloader, num_epochs):

    for epoch in range(num_epochs):
        if epoch % 2 == 0:
            model.eval()
            logsoftmax_scores_test, fattn_w_scores, hattn_w_scores = model(test_X1)
            test_predict = (torch.exp(logsoftmax_scores_test))[:,0].view(-1,1)
            test_predict_np = test_predict.detach().to('cpu').numpy()
            auc = roc_auc_score(test_label, test_predict_np)
            print('Epoch [%d] AUC: %.3f' %
                  (epoch + 1, auc))
            model.train()

        for i_batch, samples_batch in enumerate(dloader):

            X_batch, mask = samples_batch

            X_batch = X_batch.to(device)
            mask = mask.to(device)
            optimizer.zero_grad()
            logsoftmax_scores, fattn_w_scores, hattn_w_scores = model(X_batch)

            # seqid_fattnw_map.update({seqid:fattn_w_scores[c].detach().cpu() for c, seqid in enumerate(b_seqs_id)})
            # seqid_hattnw_map.update({seqid:hattn_w_scores[c].detach().cpu() for c, seqid in enumerate(b_seqs_id)})

            # __, y_pred_clss = torch.max(logsoftmax_scores, -1)

            # print('y_pred_clss.shape', y_pred_clss.shape)
            # use mask to retrieve relevant entries
            # tindx= torch.where(mask.type(torch.bool))

            # pred_class.extend(y_pred_clss[tindx].view(-1).tolist())
            # prob_scores.append((torch.exp(logsoftmax_scores.detach().cpu())).numpy())
            
            # seqs_ids_lst.extend([b_seqs_id[i] for i in tindx[0].tolist()])
            # base_pos_lst.extend(tindx[1].tolist()) # positions of target base



            # torch.cuda.ipc_collect()
            # torch.cuda.empty_cache()
            loss = lossfunc((torch.exp(logsoftmax_scores))[:,0].view(-1,1), mask)

            loss.backward()
            optimizer.step()

            if i_batch % 10 == 0:  # print every 200 mini-batches
                print('[%d, %5d] loss: %.3f' %
                      (epoch + 1, i_batch + 1, loss.item()))
        # end of epoch
        # print("+"*35)
        prob_scores_arr = np.concatenate(prob_scores, axis=0)
        # predictions_df = build_probscores_df(seqs_ids_lst, prob_scores_arr, base_pos_lst)

        
    return model

        # dump attention weights
        # if wrk_dir:
        #     ReaderWriter.dump_data(seqid_fattnw_map, os.path.join(wrk_dir, 'seqid_fattnw_map.pkl'))
        #     predictions_df.to_csv(os.path.join(wrk_dir, f'predictions.csv'))

CNN = Categ_CrisCasTransformer().to(device)
CNN = run_prediction(CNN, datloader, num_epochs= 50)

dataiter = iter(datloader)
X_batch, mask = next(dataiter)
lr = 0.01
optimizer = optim.Adam(CNN.parameters(), lr=lr)
scheduler = optim.lr_scheduler.StepLR(optimizer, step_size = 3, gamma =0.1)
# lossfunc = nn.L1Loss().to(device)
lossfunc = nn.BCEWithLogitsLoss(pos_weight=torch.Tensor([w])).to(device)
sigmoid = nn.Sigmoid()