#!/usr/bin/env python

import os
import re
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch


from botorch import fit_gpytorch_mll
from gpytorch.likelihoods import FixedNoiseGaussianLikelihood
from gpytorch.mlls import ExactMarginalLogLikelihood
from gpytorch import settings

from gskgpr import GaussianStringKernelGP
from seq2ascii import Seq2Ascii



def initialize_model(train_x, train_y, err_y, translator,num_outputs):
    model = GaussianStringKernelGP(train_x=train_x, train_y=train_y, 
                            likelihood=FixedNoiseGaussianLikelihood(noise=err_y), 
                            translator=translator, num_outputs=num_outputs)
    mll = ExactMarginalLogLikelihood(model.likelihood, model).to(device)
    return model, mll


def fit_gpytorch_model(L_values):
    best_loss = float('inf')
    best_L = None
    best_model = None
    best_mll = None

    for L in L_values:
        model, mll = initialize_model(encoded_x, train_y, err_y, translator, num_outputs=1)
        mll.model.covar_module.L = L
        
        model.train()
        mll.train()
        
        fit_gpytorch_mll(mll)
        
        model.eval()
        mll.eval()
        
        with torch.no_grad(), settings.fast_pred_var():
            output = model(*model.train_inputs)
            loss = -mll(output, model.train_targets)        
        
        print('L: %d - Loss: %.3f' % (L, loss))
        
        if loss.item() < best_loss:
            best_loss = loss.item()
            best_L = L
            best_model = model
            best_mll = mll

    print('Best L: %d - Best Loss: %.3f' % (best_L, best_loss))
    return best_model, best_mll



ELP_all=pd.read_csv('../FE_all.csv')
computed=pd.read_csv('../FE.csv')
dataset = ELP_all.merge(computed[['ELP', 'dG', 'dG_err']], on='ELP', how='left')


def append_star(name):
    if len(name) == 45:
        return name + '*****'
    else:
        return name

dataset['Sequence'] = dataset['Sequence'].apply(append_star)

translator = Seq2Ascii("./AA.blosum62.pckl")
fspace = dataset.Sequence.to_list()
translator.fit(fspace)

dataset_train=dataset.dropna(subset=['dG'], inplace=False)

dataset_train
print(len(dataset_train))

dataset_train.dG = (dataset_train.dG - dataset_train.dG.mean())/dataset_train.dG.std()
dataset_train.dG_err = dataset_train.dG_err/dataset.dG.std()
dataset_train['dG_var'] = dataset_train.dG_err**2 # pass variance instead of standard deviations

device='cpu'
encoded_x = translator.encode_to_int(dataset_train.Sequence.to_list()).to(device)
train_y = torch.tensor(dataset_train.dG.to_numpy()).float().to(device)
err_y = torch.tensor(dataset_train.dG_var.to_numpy()).float().to(device)

model, mll = fit_gpytorch_model(L_values=np.arange(25,51))

print(f'Actual sigma1: {model.covar_module.sigma1.item()}')
print(f'Actual sigma2: {model.covar_module.sigma2.item()}')
print(f'Actual L: {model.covar_module.L.item()}')

torch.save(model.state_dict(), './grid-search.pth')
