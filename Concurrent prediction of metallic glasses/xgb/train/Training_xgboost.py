import xgboost as xgb
from matplotlib import pyplot as plt
from xgboost import plot_importance
import shap
import numpy as np

dtrain = xgb.DMatrix('../SOAP_training.txt')
# dtest = xgb.DMatrix('SOAP_testing.txt')
dvalidate = xgb.DMatrix('../SOAP_validation.txt')

#param = {'max_depth': 20, 'eta': 0.0271589, 'gamma': 0.01878, 'min_child_weight': 82.011109,'subsample':0.767066, 'colsample_bytree': 0.623959}

param = {'colsample_bytree': 1.0, 'eta': 0.01, 'gamma': 0.001, 'max_depth': 36, 'min_child_weight': 43.23278323565677, 'subsample': 0.4}

evallist = [(dvalidate, 'eval'), (dtrain, 'train')]

num_round = 100000
bst = xgb.train(param, dtrain, num_round, evallist, early_stopping_rounds=100)
print(bst.best_iteration, bst.best_ntree_limit)

bst.save_model(str(num_round) + '.model')
bst.dump_model('dump.raw.txt')

######################n SHAP explainer ######################
shap.initjs()
TrainingFeature = np.loadtxt('TrainingFeatures.txt')

explainer = shap.TreeExplainer(bst,data=TrainingFeature,feature_dependence='independent')
shap_values = explainer.shap_values(TrainingFeature)
print('shap_values = ', shap_values)
shap.summary_plot(shap_values, TrainingFeature, plot_type="bar", max_display=40)
plt.savefig('SHAP.png',bbox_inches='tight')
