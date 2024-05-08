# reservoirdogs

This repository contains data and code to analyze the epidemiology of environmentally transmitted parasites of wild and domestic species around tropical forests. The data include information about parasitology of domestic dogs and wild felines (ocelots and pumas) in the Osa peninsula, Costa Rica. The goal of the project was to determine the likelihood of ongoing spillover between dogs and cats. We collected data about host distribution and abundance, and about prevalence of parasites in both types of hosts. The results of this project are set to be published in Ecosphere (Vargas Soto et al. 2024). 

Data files include:
  - prevalences.csv: a data file with information about estimated prevalence of hookworms in wild felines and domestic dogs. It includes information about the person collecting the sample, and an aggregated number of positive and negative samples per collection. For collectors KG and JJ, samples were collected in 2011/2012 and were analysed with visual parasitology, while for JV samples were collected in 2017-2019 and were analyzed using both visual parasitology and molecular methods. 
  - 2019stationData.csv: a data file with the information about camera-trap stations used to estimate the frequencies with which ocelots visit human-made trails and roads.
  - dog_record_tbl.csv: a data file with information about detections of individual dogs at different camera-trap stations
  - ocelot_record_tbl.csv: same but for individual ocelots.
  
Code includes:
  - A Wolfram mathematica notebook where all the ODE models are solved numerically, simulated, and the visualizations are created.
  - An R markdown file that estimates the prevalence based on parasitological data for each species, and calculates the frequencies of visitation to sites on human-modified landscapes. Th Rmd file also includes information about parameter estimates, as described in the text. 

