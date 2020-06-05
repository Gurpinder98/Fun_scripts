# -*- coding: utf-8 -*-
"""
Implementing Linear Regression from scratch
(okay, not exactly scatch, using Numpy) 

Date: 2 June 2020 12:52 PM
Author: Gurpinder Singh
"""
import numpy as np
import matplotlib.pyplot as plt

#Read the data
x = [] 
y = []
with open("ex1data1.txt") as data:
    lines = data.readlines()
for line in lines:
    line = line.replace("\n", "")
    x.append([float(line.split(",")[0])])
    y.append([float(line.split(",")[1])])    
  

#Make X vector, Y vector out of them
ones = [[1] for i in range(len(x))]
X = np.array(x) 
Y = np.array(y)
X = np.append(ones, X, axis=1)
theta = np.array([[0], [0]])


#implement cost function
def cost_function(X, Y, theta):
    """
    Cost function for particular theta values for X and Y
    """
    h_theta = X.dot(theta)
    SSE = np.square(h_theta-Y)
    MSE = sum(SSE)/(2*len(Y))
    return MSE

      
#implement gradient descent
def gradient_descent(theta, X, Y, alpha):
    """
    Implementing one step of gradient descent algorithm.
    """
    Err = X.dot(theta) - y
    inside_term = Err.transpose().dot(X)
    theta = theta - ((alpha/len(Y))*inside_term.transpose())    
    return theta

alpha = 0.01
J = []
Lines = []
for j in [10, 100, 500, 1000, 1500, 3000]:
    for i in range(j):
        J.append(cost_function(X, Y, theta))
        theta = gradient_descent(theta, X, Y, alpha)
        
    Lines.append([theta[0] + theta[1]*j for j in np.linspace(1, 30, num= 50)])    
    
    #data visualization
plt.figure()
plt.plot(X[:,1].transpose(), Y[:, 0].transpose(), "rx", label = "Data points")
for Line, j in zip(Lines, [10, 100, 500, 1000, 1500, 3000]):
    plt.plot(Line, "-", label = "{} interations".format(j))

plt.title("Visualizing Gradient Descent Fit over multiple iterations")
plt.xlabel("Population Size (in 10,000)")
plt.ylabel("Profits (in $10,000)")
plt.legend(loc = "best")
plt.savefig("linear_regression.jpg")

