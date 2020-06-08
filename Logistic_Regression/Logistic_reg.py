# -*- coding: utf-8 -*-
"""
Logistic regression (yes, pretty basic)

Created on Sun Jun  7 21:05:59 2020

@author: Gurpinder
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin_tnc


X = []
y = []

#Reading the data
with open("ex2data1.txt") as f:
    lines = f.readlines()
    for line in lines:
        line = line.replace("\n", "")
        X.append([float(line.split(",")[0]) , float(line.split(",")[1]) ])
        y.append([float(line.split(",")[2])])

# converting to vectors
X = np.array(X)
y = np.array(y)
        

#sigmoid function
def sigmoid(a):
    return 1/(1+ np.exp(-a))
## Alternative way
#sigVec = np.vectorize(sigmoid)


ones = np.ones((len(X), 1))
X = np.append(ones, X, axis = 1)


initial_theta = np.zeros(3)

def costfunction(theta, X, y):
    m, n = X.shape
    first_term = -y.transpose().dot(np.log(sigmoid(X.dot(theta))))
    second_term = (1-y.transpose()).dot(np.log(1 - sigmoid(X.dot(theta))))
    return (1/m)*np.sum(first_term - second_term)

def gradient(theta, X, y):
    m, n = X.shape
    theta = theta.reshape((3,1))
    term = (sigmoid(X.dot(theta)) - y).transpose()
    return (1/m)*(term.dot(X))

#Result = fmin_ncg(costfunction, initial_theta, gradient, args=(X, y), full_output=True)
    
Result = fmin_tnc(func= costfunction, x0=initial_theta, fprime=gradient, args=(X, y))
best_theta = Result[0]

#plotting the decision boundary
def plotDecision(theta, X, y):
    X_point = np.array([min(X[:,1])-2, max(X[:,1])+2])
    plot_y = (-1/theta[2])*(theta[1]*X_point + theta[0])
    
    positives = y[:,0] == 1
    negatives = y[:,0] == 0
    
    plt.figure(figsize = (8,5))        
    plt.plot(X[(positives),1], X[(positives),2], "r+", label="Admitted")
    plt.plot(X[(negatives),1], X[(negatives),2], "bo", label="Not Admitted")
    plt.xlabel("Exam 1 Score")
    plt.ylabel("Exam 2 Score")
    plt.legend(loc = "best")
    plt.title("Student Admissions based on 2 Exam Scores.")
    plt.plot(X_point, plot_y, "g-")
    plt.show()
        



