import xgboost as xgb
from matplotlib import pyplot as plt
from xgboost import plot_importance
import shap
import numpy as np


################### Load Training Model #######################
# param = {'max_depth': 6, 'eta': 1.0, 'gamma': 5.0, 'min_child_weight': 300,'save_period': 0}
param = {'colsample_bytree': 1.0, 'eta': 0.01, 'gamma': 0.001, 'max_depth': 36, 'min_child_weight': 43.23278323565677, 'subsample': 0.4}

bst = xgb.Booster(param)
bst.load_model('../train/100000.model')

#################### Input data ##############################
dtest = xgb.DMatrix('../SOAP_testing.txt')
ypred = bst.predict(dtest)

print(ypred.shape)
################### Output data ##################################
filename = './predicted.txt'
with open(filename,'x') as output:
    for i in range(len(ypred)):
        output.write("%lf\n" % ypred[i])
