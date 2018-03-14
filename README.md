# Randomness Steering Inequality

Data files associated to  

[Maximal randomness expansion from steering inequality violations using qudits](https://arxiv.org/list/quant-ph/new)  
Paul Skrzypczyk and Daniel Cavalcanti

All code is written in MATLAB and requires:
- [CVX](http://cvxr.com/) - a Matlab-based convex modeling framework
- [QETLAB](http://www.qetlab.com/) - A MATLAB Toolbox for Quantum Entanglement

It has been tested on Matlab R2015a, and CVX 2.1 

The code comprises the following:

  - [RandomnessSteeringInequality](https://github.com/paulskrzypczyk/randomnesssteeringinequality/blob/master/RandomnessSteeringInequality.m): 
  returns the min-entropy that can be certified in a 1SDI manner from a given steering inequality violation
  - [genExpGraph](https://github.com/paulskrzypczyk/randomnesssteeringinequality/blob/master/genExpGraph.m):
  script file which reproduces Fig. 1 from the paper. 
