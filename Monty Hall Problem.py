# -*- coding: utf-8 -*-
"""
The Monty Hall Problem Simulated

Created on Mon May 25 20:30:43 2020
@author: Gurpinder
"""
import random
import matplotlib.pyplot as plt
import numpy as np


def generate_doors(n):
    """
    Returns: a list of n random doors, one with a car behind it all others with Goats.
    """
    Doors = ["G" for i in range(n)]
    # choose a random place to put car behind
    Doors[random.randint(0, n-1)] = "C"
    return Doors

        
def play(Doors, switcher):
    """
    Plays one turn of the Monty Hall Game.
    Doors = list of Doors, one with car, all other with goats. (rhymes lol.)
    switcher = True/Fasle, tells if the player will make a switch or not.
    """
    
    choice = random.randint(0, len(Doors)-1) #player chooses 1 door
    #player chooses a door with no car behind
    if Doors[choice] != "C":
        #Monty keeps the car door closed, opens the rest
        closed_door = Doors.index("C")
    
    #if player chooses the door that contains car
    if Doors[choice] == "C":
        #Monty keeps one door closed, opens the rest
        while True:
            closed_door = random.randint(0, len(Doors)-1)        
            if closed_door != choice:
                break
    
    if switcher: #if player wants to switch it can do so here.
            choice = closed_door #Player switches to the door that monty kept CLOSED.
    
    return Doors[choice] == "C" #Return if the player choice door had a car behind it
        
    
def simulation(numTrials, numDoors = 3, toprint = False):
    """
    numTrials: number of times the game will be played.
    numDoors = number of doors in the game.
    
    """
    results_switch =  [] #storing the results of player that makes a switch
    results_not_switch = [] #same for the player that doesnot make a switch
    
    for i in range(numTrials):
        #generate new door set every Trial
        Doors = generate_doors(numDoors) 
        
        results_switch.append(play(Doors, True)) #Switcher plays the game
        results_not_switch.append(play(Doors, False)) #non switcher playes the game
    
    #printing the results
    if toprint:
        print("\nPrinting the simulation results:")
        print("For {} Doors, Iterations {}".format(numDoors, numTrials))
        print("The person who switched choice won {}% of times.".format((results_switch.count(True)/numTrials)*100))
        print("The person who didn't switch choice won {}% of times.".format((results_not_switch.count(True)/numTrials)*100))
        
    return (results_switch.count(True)/numTrials, results_not_switch.count(True)/numTrials)    

    
def door_number_exp(numDoors, numTrials = 100, numRepitions = 20):
    """
    Does the number of doors affect the wining ratio - does it faour either player?
    Answer: Yes, the player who switches gets the benefit.
    
    Inputs: numDoors - list of number of doors to experiment on
            numTrials - number of times to play the game
            numRepitions - number of times each door number is played 100 times and average is plotted.
    """
    #storing the averages for each door. 
    r_s = []
    r_ns = []
    for num in numDoors:
        #tempory lists to store numRepitions number of numTrials Trials.
        sc = []
        nsc = []
        for i in range(numRepitions):
            Switcher, nonSwitcher = simulation(numTrials, num)
            sc.append(Switcher)
            nsc.append(nonSwitcher)
        
        #append the means into the lists for plotting
        r_s.append(np.array(sc).mean())
        r_ns.append(np.array(nsc).mean())
    
    #plotting code    
    plt.figure(figsize = (7,7))
    plt.plot(numDoors, r_s, label = "Switcher")
    plt.plot(numDoors, r_ns, label= "Non-Switcher")
    plt.title("Simulating Monty hall problem over {} doors.".format(len(numDoors)))
    plt.xlabel("Number of Doors")
    plt.ylabel("Proportion of games won")
    plt.legend(loc="best")
    plt.show()
    


#Running the simulation
door_number_exp(np.arange(3, 103, dtype=int))         
    

