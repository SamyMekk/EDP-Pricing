# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 22:31:39 2023

@author: samym
"""

import numpy as np
import numpy.linalg as lng  # linear algebra
import matplotlib.pyplot as plt
import scipy.stats as stats  # pour la fonction de repartition de la loi normale
import pandas as pd  # pour un meilleur affichage des tables


class OptionPricer():
    def __init__(self,r=0.1,sigma=0.2,K=100,T=1):
        self.r=r
        self.sigma=sigma
        self.K=K
        self.T=T

    def Setting(self,Smin,Smax,N,I):
        self.Smax=Smax
        self.Smin=Smin
        self.I=I
        self.N=N
        self.h=(self.Smax-self.Smin)/(self.I)
        self.dt=self.T/(self.N)
        self.S=np.arange(self.Smin,self.Smax+self.h,self.h)
        self.S_squared = np.array([x**2 for x in self.S])
        self.grille=np.arange(0,self.T+self.dt,self.dt)
        alpha=pow(self.sigma,2)/(2*pow(self.h,2))*self.S_squared
        beta=self.r/(2*self.h)*self.S
        self.U = np.zeros((self.I,1))
        for i in range(self.I):
            #remplissage de U0
            self.U[i] = self.u0(self.S[i])

        self.A=np.zeros((self.I,self.I))
        for i in range(1,len(self.A)+1):
            self.A[i-1][i-1]=2*alpha[i]+self.r
            if i>1:
                self.A[i-1][i-2]=-alpha[i]+beta[i]

            if i<self.I:
                self.A[i-1][i]=-alpha[i]-beta[i]

        return self.A



    def u0(self,s):
        return max(self.K-s,0)

    def uleft(self,t):
        return (self.K*np.exp(-self.r*t)-self.Smin)

    def uright(self,t):
        return 0

    def q(self,t):
        y=np.zeros((self.I,1))
        alpha1=pow(self.sigma,2)/(2*pow(self.h,2))*self.S_squared[1]
        beta1=(self.r*self.S[1])/(2*self.h)
        y[0]=(-alpha1+beta1)*self.uleft(t)
        alphaI=pow(self.sigma,2)/(2*pow(self.h,2))*self.S_squared[self.I]
        betaI=(self.r*self.S[self.I])/(2*self.h)
        y[self.I-1]=(-alphaI-betaI)*self.uright(t)
        return y



    def computeScheme(self, SCHEME="EE"):
        self.SCHEME = SCHEME

        if self.SCHEME == "EE":
            for n in range(0, self.N):
                self.U = self.U - self.dt * (self.A @ self.U + self.q(self.grille[n]))

        if self.SCHEME == "IE":
            for n in range(0, self.N):
                self.U = lng.solve(np.identity(self.I ) + self.dt * self.A,
                                   -self.dt * self.q(self.grille[n]) + self.U)

        if self.SCHEME == "CN":
            theta = 0.5
            for n in range(0, self.N):
                F = -(theta * self.q(self.grille[n]) + (1 - theta) * self.q(self.grille[n]))
                part1 = self.U - (1 - theta) * self.dt * (self.A @ self.U)
                self.U = lng.solve(np.identity(self.I) + theta * self.dt * self.A, self.dt * F + part1)

    def plot_u(self):
        #print("Pricing par le schema ", self.SCHEME, " avec I= {} ,  N= {}".format(self.I, self.N))
        # Utilisation de self.I pour s'assurer que les indices sont valides
        plt.xlabel('Prix du sous-jacent (S)')
        plt.ylabel('Prix du Put')
        plt.title("Pricing du Put Européen par le schema {} avec T ={}, I= {},  N= {}".format(self.SCHEME,self.T,self.I,self.N))
        plt.grid(True)
        plt.legend()
        plt.plot(self.S[:self.I], self.U[:self.I])
        plt.plot(self.S[:self.I + 1], [self.u0(self.S[i]) for i in range(self.I + 1)])
        plt.show()


    def amplification(self):
        B = np.identity(self.I + 1) - self.dt * self.A
        return {"matrice": B, "norme": lng.norm(B, np.inf)}

    def CFLnumber(self):
        mu = self.dt * ((self.sigma * self.Smax / self.h) ** 2)
        return mu

    def P_Interpolate(self, sbar):
        i = np.argmax(self.S >= sbar)
        ubar = ((self.S[i] - sbar) * self.U[i - 1] + (sbar - self.S[i - 1]) * self.U[i]) / self.h
        return ubar

   


class PutPrice(OptionPricer):
    def __init__(self,r=0.1,sigma=0.2,K=100,T=1):
        super().__init__(r, sigma,K,T)
        
    def u0(self,s):
         return max(self.K-s,0)

    def uleft(self,t):
         return (self.K*np.exp(-self.r*t)-self.Smin)

    def uright(self,t):
         return 0
     
        
    def plot_u(self):
        #print("Pricing par le schema ", self.SCHEME, " avec I= {} ,  N= {}".format(self.I, self.N))
        # Utilisation de self.I pour s'assurer que les indices sont valides
        plt.xlabel('Prix du sous-jacent (S)')
        plt.ylabel('Prix du Put')
        plt.title("Pricing du Put Européen par le schema {} avec T ={}, I= {},  N= {}".format(self.SCHEME,self.T,self.I,self.N))
        plt.grid(True)
        plt.legend()
        plt.plot(self.S[:self.I], self.U[:self.I])
        plt.plot(self.S[:self.I + 1], [self.u0(self.S[i]) for i in range(self.I + 1)])
        plt.show()
    def BlackScholes(self, sbar):
        dplus = (np.log(sbar / self.K) + self.r + 0.5 * self.T * (self.sigma ** 2)) / (self.sigma * np.sqrt(self.T))
        dmoins = (np.log(sbar / self.K) + self.r - 0.5 * self.T * (self.sigma ** 2)) / (self.sigma * np.sqrt(self.T))
        vbar = np.exp(-self.r * self.T) * self.K * stats.norm.cdf(-dmoins) - sbar * stats.norm.cdf(-dplus)
        return vbar


class CallPrice(OptionPricer):
    def __init__(self,r=0.1,sigma=0.2,K=100,T=1):
        super().__init__(r, sigma,K,T)
      
    def u0(self,s):
        return max(s-self.K,0)
    def uleft(self,t):
        return 0
    def uright(self,t):
        return (self.Smax-self.K*np.exp(-self.r*t))
      
    def plot_u(self):
        #print("Pricing par le schema ", self.SCHEME, " avec I= {} ,  N= {}".format(self.I, self.N))
        # Utilisation de self.I pour s'assurer que les indices sont valides
        plt.xlabel('Prix du sous-jacent (S)')
        plt.ylabel('Prix du Put')
        plt.title("Pricing du Call Européen par le schema {} avec T ={}, I= {},  N= {}".format(self.SCHEME,self.T,self.I,self.N))
        plt.grid(True)
        plt.legend()
        plt.plot(self.S[:self.I], self.U[:self.I])
        plt.plot(self.S[:self.I + 1], [self.u0(self.S[i]) for i in range(self.I + 1)])
        plt.show()

    def BlackScholes(self, sbar):
        dplus = (np.log(sbar / self.K) + self.r + 0.5 * self.T * (self.sigma ** 2)) / (self.sigma * np.sqrt(self.T))
        dmoins = (np.log(sbar / self.K) + self.r - 0.5 * self.T * (self.sigma ** 2)) / (self.sigma * np.sqrt(self.T))
        vbar = -np.exp(-self.r * self.T) * self.K * stats.norm.cdf(dmoins) +sbar * stats.norm.cdf(dplus)
        return vbar

        
    