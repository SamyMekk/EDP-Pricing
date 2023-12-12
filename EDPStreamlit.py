# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 22:31:22 2023

@author: samym
"""

from mainEDP import *
import streamlit as st





st.title("Application pour le calcul d'options par méthode des différences finies")

def user_input():
    option=st.selectbox("Quelle Option voulez-vous pricer par la méthode de Monte-Carlo ",('Put Européen','Call Européen'))
    return option





df=user_input()

st.set_option('deprecation.showPyplotGlobalUse', False)
st.header("Vous avez choisi de prendre  le " + str(df))


SpotPrice=st.sidebar.number_input("Choississez le spot",value= 100)
Strike=st.sidebar.number_input("Choississez le strike K",value= 100)
Maturité=st.sidebar.number_input("Choississez la maturité T",value= 1)
Volatilité=st.sidebar.number_input("Choississez σ",value= 0.1)
TauxSansRisque=st.sidebar.number_input("Choississez r",value= 0.02)
Schéma=st.selectbox("Choissisez le schéma de discrétisation : ",("EE","IE","CN"))
data={'''SpotPrice''': SpotPrice,
         '''Strike''':Strike,
         '''Maturité''':Maturité,
         '''Volatilité''':Volatilité,
         '''Taux Sans Risque''':TauxSansRisque
         }
Parametres=pd.DataFrame(data,index=["Caractéristique de l'option"])
st.subheader("Voici les caractéristiques de l'option : ")
st.write(Parametres)   
    


if df=="Put Européen":
    Option1=PutPrice(TauxSansRisque,Volatilité,Strike,Maturité)
    N=st.number_input("Choissisez le nombre N",value=100)
    I=st.number_input("Choississez le nombre I",value=20)
    Smin=st.number_input("Choississez Smin",value=0)
    Smax=st.number_input("Choississez Smax",value=200)
    
    Option1.Setting(Smin=Smin,Smax=Smax,N=N,I=I)
    Option1.computeScheme(SCHEME=Schéma)
    
    
    st.pyplot(Option1.plot_u())
    # st.write(Option1.BlackScholes(SpotPrice))
    
    
if df=="Call Européen":
    Option1=CallPrice(TauxSansRisque,Volatilité,Strike,Maturité)
    N=st.number_input("Choissisez le nombre N",value=100)
    I=st.number_input("Choississez le nombre I",value=20)
    Smin=st.number_input("Choississez Smin",value=0)
    Smax=st.number_input("Choississez Smax",value=200)
    
    Option1.Setting(Smin=Smin,Smax=Smax,N=N,I=I)
    Option1.computeScheme(SCHEME=Schéma)
    
    
    st.pyplot(Option1.plot_u())
    # st.write(Option1.BlackScholes(SpotPrice))
    
    



