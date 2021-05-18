## Introduction to Q-Factor sideband subtraction at GlueX
[The Original Paper introducing the Q-Factor method - Multivariate side-band subtraction using probabilistic event weights](https://arxiv.org/pdf/0809.2548.pdf)

Sideband subtraction is a common technique used in particle physics to separate signal and background. For this technique to work the kinematics of the sideband region should be the same as the kinematics of the signal region. Also, the yield in higher dimensional sidebands shrink exponentially which can enlarge statistical uncertainties. This technique is relatively easy to implement. The Q-Factors technique is an alternative method to separate signal and background. 

An implementation of this technique is develoepd to calculating Q-factors in the reaction <img src="https://render.githubusercontent.com/render/math?math=\gamma p\rightarrow\pi^0\eta p \rightarrow 4\gamma p"> at GlueX. Q-factors is an event-by-event multivariate sideband subtraction technique. The only requirement is the knowledge of the signal and background distribution of some discriminating variable. 
1. First the nearest neighbors, under some set of phase space variables, is found for a given event
2. Distribution of the discriminating variable is filled with the nearest neighbors
3. Fit the above distribution with the known/assumed signal and bkg distribution
4. Calculate Q-factor as the signal fraction
5. Do this for all entries

See documentation for more details
