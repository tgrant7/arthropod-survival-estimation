# arthropod-survival-estimation
UPDATE FALL 2021: the .bug file has been updated with an improved, more accurate degree day model based on Asclepias syriaca development times instead of development times from Zalucki (1982). The new bug file is called MonarchModelSyriaca.bug. This file should be used instead of MonarchModel.bug (unless you are in Australia and your monarchs are feeding on Gomphocarpus fruticosus). A new R file is also provided (MonarchFieldDataAnalysiswithSyriacaModel.R), though you will notice it only has minor changes. 

The new model and testing on the new model will be described in a monarch butterfly review paper to be submitted soon. At this time I do not recommend using stage-specific survival probabilities from the Syriaca model as they appear to now be much less "identifiable." The overall cumulative survival probability is still solid, of course. 

I am leaving the older code here because it can still be useful. 


--------------------------------------------------------------------------------------------

Code for Grant et al. (2020) Estimating arthropod survival probability from field counts:  A case study with monarch butterflies.

MonarchModel.bug - JAGS model specification for survival estimation model

SecondStage2.bug - model specification for meta-analysis regression

SimulateandAnalyzeMonarchCounts.R - R code for simulation and analysis of monarch butterfly field counts.  Begin here for a new simulation or analysis.  

MonarchFieldDataAnalysis.R - R code for analysis of monarch field data reported in the paper.  Data is included.  

SimulateandAnalyzeMonarchCountswithScen1-3.R - R code for Scenario 1-3 from the manuscript

SimulateandAnalyzeScen4-5.R - R code for Scenarios 4-5 from the manuscript

Flowchart of Model Objects Simple.pptx - to aid in comprehension of the model code, a flowchart of the model code is provided.  Originally, I intended to include ADD and SDC in the model specification (the .bug file) so that only the temp data and counts would have to be provided to JAGS.  But somehow the ADD and SDC matrices threw an error in JAGS, and after discussing it with Martyn Plummer the easiest solution was to move ADD and SDC outside.  
