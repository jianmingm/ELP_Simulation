import matplotlib.pyplot as plt
import torch
import gpytorch
from gskernel import GenericStringKernel
from botorch.models.gpytorch import GPyTorchModel
from tqdm import tqdm


class GaussianStringKernelGP(gpytorch.models.ExactGP, GPyTorchModel):
    def __init__(self, train_x, train_y, likelihood, translator, num_outputs=1, **kwargs):
        super(GaussianStringKernelGP, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = GenericStringKernel(translator, **kwargs)
        self._num_outputs = num_outputs  # Use a private variable
    
    @property
    def num_outputs(self):
        return self._num_outputs
    
    @num_outputs.setter
    def num_outputs(self, value):
        self._num_outputs = value
    
    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)