# arthropod-survival-estimation

Code for Grant et al. (2020) Estimating arthropod survival probability from field counts:  A case study with monarch butterflies.

MonarchModel.bug - JAGS model specification for survival estimation model

SecondStage2.bug - model specification for meta-analysis regression

SimulateandAnalyzeMonarchCounts.R - R code for simulation and analysis of monarch butterfly field counts.  Begin here for a new simulation or analysis.  

MonarchFieldDataAnalysis.R - R code for analysis of monarch field data reported in the paper.  Data is included.  

SimulateandAnalyzeMonarchCountswithScen1-3.R - R code for Scenario 1-3 from the manuscript

SimulateandAnalyzeScen4-5.R - R code for Scenarios 4-5 from the manuscript

Flowchart of Model Objects Simple.pptx - to aid in comprehension of the model code, a flowchart of the model code is provided.  Originally, I intended to include ADD and SDC in the model specification (the .bug file) so that only the temp data and counts would have to be provided to JAGS.  But somehow the ADD and SDC matrices threw an error in JAGS, and after discussing it with Martyn Plummer the easiest solution was to move ADD and SDC outside.  
