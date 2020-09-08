
import sys
import time
import joblib
import numpy as np
import pandas as pd
import preprocessing as p
import matplotlib.pyplot as plt

from itertools import cycle
from multiprocessing import Pool
from sklearn.preprocessing import label_binarize
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import auc, roc_curve, roc_auc_score
from sklearn.metrics import precision_recall_curve, average_precision_score
from sklearn.metrics import confusion_matrix,classification_report, accuracy_score
from sklearn.model_selection import cross_val_score, GridSearchCV, train_test_split


start_time = time.time()

#def read_hdf(filename):
 #   'converts a filename to a pandas dataframe'
  #  return pd.read_hdf(filename)
    
def read_feather(filename):
    'converts a filename to a pandas dataframe'
    return pd.read_feather(filename)

TRAIN = ['sawa','vgh117','goku','Pseudomonas','Ecoli','Staphylococcus','Enterococcus']
file_list = []
for strain in TRAIN:
    file_list.append('/big6_disk/shiuanrung107/cdspolish/0712/alignment/weight4/{strain}/alignment.feather'.format(strain=strain))


with Pool(processes=8) as pool:
    df_list = pool.map(read_feather, file_list)
    train = pd.concat(df_list, ignore_index=True)
    pool.close()
    pool.join()

train = p.double_ins_chaos(train)

Y_train = train['label']
print('Label count: \n',Y_train.value_counts())
Y_train.astype(object)

X_train = p.preprocessing(train)
#X_train = pd.get_dummies(X_train, columns=['draft','homopolymer'])
#X_train = pd.get_dummies(X_train, columns=['indel_3n'])
X_train = pd.DataFrame(X_train.drop(['label', 'position'], axis=1))

print(X_train.columns)
print('X_train shape: ', X_train.shape)
print('Y_train shape: ', Y_train.shape)


X_train, X_val, Y_train, Y_val = train_test_split(X_train, Y_train, test_size = 0.1, random_state=1)

 
clf = RandomForestClassifier(n_estimators=20, random_state=42,n_jobs=-1)
clf.fit(X_train, Y_train) 

joblib.dump(clf, 'RF.model')
y_pred = clf.predict(X_val) 
y_score= clf.predict_proba(X_val)
print(clf.get_params())
print("Confusion Matrix: ")
print( confusion_matrix(Y_val, y_pred))
print ("Accuracy : ", accuracy_score(Y_val,y_pred)*100)
print(classification_report(Y_val, y_pred))


import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import metrics
cm =confusion_matrix(Y_val,y_pred)
plt.figure(figsize=(15,15))
ax= plt.subplot()
sns.heatmap(cm, annot=True, linewidths=.5, square = True, cmap = 'Blues_r', fmt='g', annot_kws={"size":15}, ax=ax);
ax.set_xticklabels( ['A','T','C','G','deletion','keep'])
ax.set_yticklabels( ['A','T','C','G','deletion','keep'])

plt.ylabel('Actual label');
plt.xlabel('Predicted label');
all_sample_title = 'Accuracy Score: {0}'.format(accuracy_score(Y_val,y_pred)*100)
plt.title(all_sample_title, size = 15);
#ax.xaxis.set_ticklabels(targets); ax.yaxis.set_ticklabels(targets);
plt.savefig('confusion.png')

end_time = time.time()
print(end_time-start_time)

'''
print("======")
importances = clf.feature_importances_
feat_labels = X_train.columns
print("importances: ",importances)
indices = np.argsort(importances)[::-1]
for f in range(X_train.shape[1]):
    print("%2d) %-*s %f" % (f + 1, 30, feat_labels[indices[f]], importances[indices[f]]))
plt.figure(figsize=(20,12))
plt.title("Feature importances",fontsize = 18)
plt.ylabel("import level",fontsize = 15,rotation=90)
#plt.rcParams['font.sans-serif'] = ["SimHei"]
plt.rcParams['axes.unicode_minus'] = False
for i in range(X_train.shape[1]):
    plt.bar(i,importances[i],color='orange',align='center')
    plt.xticks(np.arange(feat_labels.shape[0]),feat_labels,rotation=90,fontsize=15)
plt.savefig("importances.png")
print("======")





# Compute ROC curve and ROC area for each class
fpr = dict()
tpr = dict()
roc_auc = dict()
# Binarize the output
Y_val = label_binarize(Y_val, classes=[0, 1, 2, 3, 4, 5])
# one vs rest
n_classes=Y_val.shape[1]
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(Y_val[:, i], y_score[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])
# Compute micro-average ROC curve and ROC area

fpr["micro"], tpr["micro"], _ = roc_curve(Y_val.ravel(), y_score.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

micro_auc = roc_auc_score(Y_val, y_score, average='micro')
lw = 2
# First aggregate all false positive rates
all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))
# Then interpolate all ROC curves at this points
mean_tpr = np.zeros_like(all_fpr)
for i in range(n_classes):
    mean_tpr += np.interp(all_fpr, fpr[i], tpr[i])
# Finally average it and compute AUC
mean_tpr /= n_classes
fpr["macro"] = all_fpr
tpr["macro"] = mean_tpr
roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

macro_auc = roc_auc_score(Y_val, y_score, average='macro')
print(roc_auc)
print('micro auc:', micro_auc)
print('macro auc:', macro_auc)
# Plot all ROC curves
plt.figure(figsize=(20,12))
plt.plot(fpr["micro"], tpr["micro"], label='micro-average ROC curve (area = {0:0.2f})'.format(roc_auc["micro"]),color='deeppink', linestyle=':', linewidth=4)
plt.plot(fpr["macro"], tpr["macro"],label='macro-average ROC curve (area = {0:0.2f})'.format(roc_auc["macro"]),color='navy', linestyle=':', linewidth=4)
colors = ['aqua', 'darkorange', 'cornflowerblue','purple','red','green']
for i, color in zip(range(n_classes), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=lw,label='ROC curve of class {0} (area = {1:0.2f})'.format(i, roc_auc[i]))
plt.plot([0, 1], [0, 1], 'k--', lw=lw)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate',size=20)
plt.ylabel('True Positive Rate',size=20)
plt.title('Extension of ROC curve to multi-class',size=20)
plt.legend(loc="lower right",fontsize='x-large')
plt.savefig("ROC.png")


#PR-curve
# For each class
precision = dict()
recall = dict()
average_precision = dict()
for i in range(n_classes):
    precision[i], recall[i], _ = precision_recall_curve(Y_val[:, i],y_score[:, i])
    average_precision[i] = average_precision_score(Y_val[:, i], y_score[:, i])
# A "micro-average": quantifying score on all classes jointly
precision["micro"], recall["micro"], _ = precision_recall_curve(Y_val.ravel(),y_score.ravel())
average_precision["micro"] = average_precision_score(Y_val, y_score,average="micro")
#print('Average precision score, micro-averaged over all classes: {0:0.2f}'.format(average_precision["micro"]))
# setup plot details
colors = cycle(['navy', 'turquoise', 'darkorange', 'cornflowerblue', 'teal','red'])
plt.figure(figsize=(20, 12))
f_scores = np.linspace(0.2, 0.8, num=4)
lines = []
labels = []
for f_score in f_scores:
    x = np.linspace(0.01, 1)
    y = f_score * x / (2 * x - f_score)
    l, = plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.2)
    plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.9, y[45] + 0.02))
lines.append(l)
labels.append('iso-f1 curves')
l, = plt.plot(recall["micro"], precision["micro"], color='purple', linestyle=':', lw=2)
lines.append(l)
labels.append('micro-average Precision-recall (area = {0:0.2f})'.format(average_precision["micro"]))
for i, color in zip(range(n_classes), colors):
    l, = plt.plot(recall[i], precision[i], color=color, lw=2)
    lines.append(l)
    labels.append('Precision-recall for class {0} (area = {1:0.2f})'.format(i, average_precision[i]))
fig = plt.gcf()
fig.subplots_adjust(bottom=0.25)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Extension of Precision-Recall curve to multi-class')
plt.legend(lines, labels, prop=dict(size=14),loc="lower left")
plt.savefig("PR_curve.png")

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

#start = time.time()
# Standardizing the features
x_trans = StandardScaler().fit_transform(X_train)
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x_trans)
#principalComponents = pca.fit_transform(X_train)
principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2'])

finalDf = pd.concat([principalDf, pd.DataFrame(Y_train)], axis = 1)
#finalDf

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
targets = [0,1,2,3,4,5] 
colors = ['black', 'blue', 'purple', 'yellow', 'lime', 'red'] 
for target, color in zip(targets,colors):
    indicesToKeep = finalDf.iloc[:,-1] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()
plt.savefig('pca_5.png')

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
targets = [0,1,2,3,4]
colors = ['black', 'blue', 'purple', 'yellow', 'lime'] 
for target, color in zip(targets,colors):
    indicesToKeep = finalDf.iloc[:,-1] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()
plt.savefig('pca_4.png')

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
targets = [0,1,2,3]
colors = ['black', 'blue', 'purple', 'yellow'] 
for target, color in zip(targets,colors):
    indicesToKeep = finalDf.iloc[:,-1] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()
plt.savefig('pca_3.png')

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
targets = [0,1,2]
colors = ['black', 'blue', 'purple'] 
for target, color in zip(targets,colors):
    indicesToKeep = finalDf.iloc[:,-1] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()
plt.savefig('pca_2.png')

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
targets = [0,1]
colors = ['black', 'blue'] 
for target, color in zip(targets,colors):
    indicesToKeep = finalDf.iloc[:,-1] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()
plt.savefig('pca_1.png')

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
targets = [0]
colors = ['black'] 
for target, color in zip(targets,colors):
    indicesToKeep = finalDf.iloc[:,-1] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()
plt.savefig('pca_0.png')
'''


