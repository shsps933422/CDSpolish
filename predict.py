import sys
import time
import joblib
import numpy as np
import pandas as pd
import preprocessing as p

from multiprocessing import Pool
from sklearn import preprocessing
from sklearn.metrics import confusion_matrix,classification_report

start = time.time()
predict_alignment = sys.argv[1]
model_file = sys.argv[2]

X = pd.read_feather(predict_alignment)
Y = X['label'].values
position = X['position']
X = p.double_ins_chaos(X)
X = p.preprocessing(X)
X = pd.DataFrame(X.drop(['position','label'], axis=1))

print(X.columns)
print('training shape: ', X.shape)
size = 100000
list_of_X = [X.loc[i:i+size-1,:] for i in range(0, len(X),size)]

def predict_class(input_data):
    #model = svm_load_model(model_file)
    model = joblib.load(model_file).set_params(n_jobs=1)
    p_label=model.predict(input_data.values.tolist())
    #p_label, p_acc, p_val = svm_predict([], input_data.values.tolist(), model)
    return p_label

pool = Pool(32)
results = pool.map(predict_class, list_of_X)

pool.close()
pool.join()
result = []
for i in results:
    result.extend(i)
print(len(result))
ACC, MSE, SCC = evaluations(Y.tolist(), result)
print(ACC)
print(confusion_matrix(Y.tolist(),result))

sub = X
sub['position'] = position
sub['label'] = Y
sub['predict'] = result
sub.to_feather('result.feather')

debug = sub[sub['predict'] != sub['label']]
debug.to_csv('debug.csv', index=False)
right = sub[sub['predict'] == sub['label']]
right.to_csv('right.csv', index=False)
end = time.time()
print(end-start)
