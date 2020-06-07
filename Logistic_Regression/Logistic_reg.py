# -*- coding: utf-8 -*-
"""
Logistic regression

Created on Sun Jun  7 21:05:59 2020

@author: Gurpinder
"""
import numpy as np
import matplotlib.pyplot as plt

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
        
#data visualization
positives = y[:,0] == 1
negatives = y[:,0] == 0

plt.figure(figsize = (8,5))        
plt.plot(X[(positives),0], X[(positives),1], "r+", label="Admitted")
plt.plot(X[(negatives),0], X[(negatives),1], "bo", label="Not Admitted")
plt.xlabel("Exam 1 Score")
plt.ylabel("Exam 2 Score")
plt.legend(loc = "best")
plt.title("Student Admissions based on 2 Exam Scores.")
plt.show()

#sigmoid function
def sigmoid(a):
    return 1/(1+ np.exp(-a))
## Alternative way
#sigVec = np.vectorize(sigmoid)

ones = np.ones((len(X), 1))
scaled = sigmoid(X)  
scaledX = np.append(ones, scaled, axis = 1)  

initial_theta = np.zeros((len(scaledX[0]), 1))

def costfunction(theta, X, y):
    m, n = X.shape
    first_term = -y.transpose().dot(np.log(sigmoid(X.dot(theta))))
    second_term = (1-y.transpose()).dot(np.log(1 - sigmoid(X.dot(theta))))
    return (1/m)*np.sum(first_term - second_term)

def gradient(theta, X, y):
    m, n = X.shape
    term = (sigmoid(X.dot(theta)) - y).transpose()
    return (1/m)*(term.dot(X)).transpose()

