import numpy as np 
import sklearn 
import sklearn.metrics
import sklearn.cross_validation 
import sklearn.linear_model 
import sklearn.ensemble 


def _eval_classifier_cv(classifier_name, clf, X, Y, cv=5):
  acc = np.mean(sklearn.cross_validation.cross_val_score(clf, X, Y, cv = cv))
  print classifier_name, "cross-validation accuracy", acc 
  auc = np.mean(sklearn.cross_validation.cross_val_score(clf, X, Y, cv = cv, scoring='roc_auc'))
  print classifier_name, "cross-validation AUC", auc
  return acc, auc
  
def eval_cv(X, Y, logistic_regression = True, n_trees=50, cv = 5):
  lr = sklearn.linear_model.LogisticRegression()
  _eval_classifier_cv("Logistic Regression", lr, X, Y, cv)
  rf = sklearn.ensemble.RandomForestClassifier(n_trees)
  _eval_classifier_cv("Random Forest", rf, X, Y, cv)
  
def _eval_classifier(classifier_name, clf, X_train, Y_train, X_test, Y_test):
  clf.fit(X_train, Y_train)
  Y_pred = clf.predict(X_test)
  acc = np.mean(Y_test == Y_pred)
  print classifier_name, "Accuracy", acc
  Y_prob = clf.predict_proba(X_test)
  auc = sklearn.metrics.roc_auc_score(Y_test, Y_prob)
  print classifier_name, "AUC", auc
  return acc, auc 
  
def eval_split(X_train, Y_train, X_test, Y_test, n_trees = 50):
  """
  Given a dataset split into training and testing parts, 
  evaluate several classifiers for accuracy and area-under-ROC
  """
  lr = sklearn.linear_model.LogisticRegression()
  _eval_classifier("Logistic Regression", lr, X_train, Y_train, X_test, Y_test)  
  rf = sklearn.ensemble.RandomForestClassifier(n_trees)
  _eval_classifier("Random Forest", rf,    X_train, Y_train, X_test, Y_test) 
  
  
  