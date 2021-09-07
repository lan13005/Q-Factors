## Updates:
9/7/2021: Includes program to aggregate results from run.py's runAllPhaseCombos variable convertImagesToPDF.py

8/29/2021: run.py outputs progress by checking the output in processLog.txt files in the log directory

08/27/2021: Updated method to set fit range. We now implement it as a restriction on the potential neighbors which can be set using neighborReqs variable in run.py. neighborReqs uses a semicolon separated string of conditions (currently only accepts less than or greater signs). This has the benefit of allowing one to select the same set of neighbors even if the {discrim var, fit range} combination was different. The previous approach would select a different set of neighbors since the fitRange would drop different entries. This also required the loading of extra variables which is defined in run.py as extraVars. Also updated how the data is stored. A single vector is now used to hold all the data in RAM as a contiguous chunk as opposed to having separate arrays to hold the different variables. This allows for easier querying of data based on the name of branch/variable and more efficient for memory if we had duplicated requests for variables. 

08/24/2021: Weights allow semicolon separated branch names

08/11/2021: Reogranization of the folders to allow for looping over all potential subsets of phase space variables. Useful in determining final phase space. makePlots also draws 2D matched thrown distributions if it can

## Introduction to Q-Factor sideband subtraction at GlueX
[The Original Paper introducing the Q-Factor method - Multivariate side-band subtraction using probabilistic event weights](https://arxiv.org/pdf/0809.2548.pdf)

Sideband subtraction is a common technique used in particle physics to separate signal and background. The technique is relatively easy to implement but has some drawbacks. For this technique to work the kinematics of the sideband region should be the same as the kinematics of the signal region. Also, the yield in higher dimensional sidebands shrink exponentially which can enlarge statistical uncertainties. The Q-Factors technique is an alternative method to separate signal and background. 

An implementation of this technique is developed to calculat Q-factors in the reaction <img src="https://render.githubusercontent.com/render/math?math=\gamma p\rightarrow\pi^0\eta p \rightarrow 4\gamma p"> at GlueX. Q-factors is an event-by-event multivariate sideband subtraction technique. The only requirement is the knowledge of the signal and background distribution of some discriminating variable. 
1. First the nearest neighbors, under some set of phase space variables, is found for a given event
    - All entries in the tree will be potential neighbors. Make sure to feed only the events in the region of interest ( i.e. determined by your PDF fit range)
2. Distribution of the discriminating variable is filled with the nearest neighbors
3. Fit the above distribution with the known/assumed signal and bkg distribution
4. Calculate Q-factor as the signal fraction
5. Do this for all entries


Example usage:
1. Clone repository
2. Update run.py to match your data. Most important variables are under STANDARD comment
3. Copy one of the config files in auxilliary/pdfTemplates/ if you are looking for an example of 1D or 2D pdf fits and replace configPDFs.h with it.
    - You have to update the template to match your reaction
4. Finally type: run.py 11
    - This will run the program to do the fitting (first digit) and plotting (second digit)
    - If plotting you also need to update histsToMake variable in makePlots.C to draw w.e. variable you want

Extra...
5. convertROOTtoImage.C will aggregate the ROOT files in the histograms folder into a single pdf file
6. flat_to_amptools.C can be used to convert the flat tree outputs of the Q-factor program into amptools format. Not made general yet


See documentation for more details
