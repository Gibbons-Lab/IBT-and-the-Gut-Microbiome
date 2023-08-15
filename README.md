# IBT-and-the-Gut-Microbiome
Repository containing the scripts, the intermediate files, and the data used in the manuscript, "Island biogeography theory and the gut microbiome: why taller people tend to harbor more diverse gut microbiomes"

# Steps
Before we can visualize and run regressions on our data, we need to generate the diversity metrics for each cohort, which is completed in the notebooks: HumanDiversity.ipynb and VertebrateDiversity.ipynb. These notebooks use the packages: pandas and qiime2. Each notebook: 
  - Converts our ASV or OTU table into a qiime feature table
  - Summarizes the feature table
  - Rarefies the feature table
  - Quantifies the Simpson's Diversity for each sample in our cohort
  - Converts the Simpson's Diversity returned by qiime from Simpson (1-D) to Simpson (1/D)

Now that we have Simpson's Diversity (1/D) calculated for our datasets, we can next visualize this data and run regressions. The American Gut cohort and the Arivale cohort are both analyzed in All_Human_Data.ipynb. The vertebrate datasets: Godon et al., Song et al., and Groussin et al. are analyzed in All_Vertebrate_Data.ipynb. Each notebook:
  - log transforms height or mass and Simpson's Diversity
  - For our human datasets:
    - we add in relevant covariates to our datasets for inclusion in regressions
    - Standardize the data
  - Plot using the Seaborn package
  - Run OLS regressions with the Statsmodels package
