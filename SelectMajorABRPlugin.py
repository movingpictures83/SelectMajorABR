#!/usr/bin/env python
# coding: utf-8

# # Supplementary Notebook 4: Select Major ABR genes
# ## Paper: Novel Approach for Microbiome Analysis Using Bacterial Replication Rates and Causal Inference to Determine Resistome Potential
# ### Vitalii Stebliankin, Musfiqur Sazal, Camilo Valdes, Kalai Mathee, and GiriNarasimhan
# 
# #### Dataset: Gibson et al. (BioProject ID: PRJNA301903)
# 
# In this notebook we will select ABR genes that present in at least 5% of the study samples

# In[1]:


import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import os

import PyPluMA
import PyIO
class SelectMajorABRPlugin:
 def input(self, inputfile):
  self.parameters = PyIO.readParameters(inputfile)
 def run(self):
     pass
 def output(self, outputfile):
  amr = PyPluMA.prefix()+"/"+self.parameters["amr"]#"A-out/ptr_amr.csv"
  major_file=PyPluMA.prefix()+"/"+self.parameters["major"]#"analysis-out/1-FilteringPTR/PTR_species_filtered_metadata_major.csv"

  out_dir=outputfile#"analysis-out/4-Select_major_ABR"
    
  out_file = out_dir+"/PTR_species_filtered_metadata_major_AMR.csv"


  clinical_vars = ["Day_of_Life", "PostMenst_Age", "Gestational_Age",
                       "Birthweight", "Gentamicin", "Cefazolin","Ampicillin", "Trimethoprim.Sulfamathoxazole", "Meropenem",
                       "Vancomycin", "Ticarcillin.Clavulanate", "Clindamycin", "Cefotaxime", "Total_abx", "r_Gentamicin",
                       "r_Meropenem", "r_Ticarcillin.Clavulanate", "r_Vancomycin", "r_Ampicillin",
                       "r_Cefotaxime","r_TOTAL","Human_Milk","Maternal_Milk", "Donor_Milk", "Formula","Fortification","Vitamin_A",
                       "Caffeine","Iron","Furosemide_Lasix","m_ampicillin","m_ceftriaxone","m_azithromycin",
                       "m_amoxicillin", "m_cefazolin","m_erythromycin","m_gentamicin","m_penicillin","m_vancomycin",
                       "m_clindamycin","m_cefotaxime", "dur_membrane_rupture","Total.Antibiotic.Days", "Cohort", "CRIB.II.Score"]

  threshold = int(self.parameters["threshold"])#20
  #method="top_largest"
  method= self.parameters["method"]#"top_present"

  df = pd.read_csv(amr)
  amr_dict={"amr":[],"value":[], "nvalues":[]}

  for col in df.columns:
    if "gb|" in col:
        amr_dict["value"].append(df[col].mean())
        amr_dict["amr"].append(col)
        amr_dict["nvalues"].append(len(df[df[col]>0]))


  # Threshold. Top 200
  amr_sorted = sorted(amr_dict["value"])
  amr_df = pd.DataFrame(amr_dict)
  print(len(amr_df))

  if method=="top_largest":
    # Top Top 10 %:
    top_200 = amr_sorted[-threshold]
    #top_ten = amr_df["value"].quantile(.80)
    amr_df = amr_df[amr_df["value"]>top_200]
    amr_list = list(amr_df["amr"])

  if method=="top_present":
    amr_df = amr_df[amr_df["nvalues"]>threshold]
    amr_list = list(amr_df["amr"])

  print(len(amr_df))

  new_amr_df = df[["sample"]+amr_list]
  new_amr_df.index = new_amr_df["sample"]
  new_amr_df = new_amr_df.drop("sample", axis=1)
  new_amr_df.rename(columns=lambda x: x.replace("\n","").split("|")[-2] + x.replace("\n","").split("|")[-1], inplace=True)

  #new_amr_df["sample"] = new_amr_df.index
  # Merge with major df:
  major_df = pd.read_csv(major_file)

  major_df = major_df.merge(new_amr_df, on="sample", how="left")
  major_df = major_df.fillna(1)

  cols_to_drop = ["Clindamycin","Cefotaxime"]
  for col in cols_to_drop:
    major_df = major_df.drop(col, axis=1)

  major_df = major_df.drop("sample", axis=1)
  major_df = major_df.drop("Individual", axis=1)
  major_df = major_df.drop("AveragePTR", axis=1)
  #major_df = major_df.drop("Cohort", axis=1)
  major_df = major_df.drop("Antibiotic_Treatment", axis=1)
  major_df = major_df.drop("Trimethoprim-Sulfamathoxazole", axis=1)
  major_df = major_df.drop("Antibiotic_Treatment_unfiltered", axis=1)



  #major_df.to_csv(out_file, index=False)
  major_df.to_csv(out_file, index=False)


  # In[2]:


  # Plot ABR profile

  amr_df = pd.read_csv(amr)#"A-out/ptr_amr.csv")
  amr_cols = []
  for col in amr_df.columns:
    if "ARO" in col:
        amr_cols.append(col)
        
  amr_df.index = amr_df['Cohort']
  amr_df = amr_df[amr_cols]
  amr_df


  # In[3]:


  # check significance;
  from scipy import stats
  amr_df_tmp = amr_df.copy()
  amr_df_tmp['cohort'] = amr_df_tmp.index

  control_df = amr_df_tmp[amr_df_tmp['cohort']=='Control']
  ab_df = amr_df_tmp[amr_df_tmp['cohort']=='Antibiotic']

  cc_group = []
  ac_group = []

  for amr in amr_cols:
    U, p = stats.mannwhitneyu(control_df[amr], ab_df[amr])
    if p<0.05:
        if control_df[amr].mean()>ab_df[amr].mean():
            cc_group.append(amr)
        else:
            ac_group.append(amr)
  print("{} genes are associated with CC".format(len(cc_group)))
  print("{} genes are associated with AC".format(len(ac_group)))


  # In[4]:


  amr_df


  # In[5]:




  amr_df = amr_df.T

  amr_aro_df = amr_df
  amr_aro_df['ARO Accession'] = amr_df.index
  amr_aro_df['ARO Accession'] = amr_df['ARO Accession'].apply(lambda x: x.split('|')[-2])
  amr_aro_df = amr_aro_df[['ARO Accession']]

  amr_df
  #   Merge ARO accession with gene family

  aro_index=out_dir+"/aro_index.tsv"#"analysis-out/4-Select_major_ABR/aro_index.tsv"
  aro_df = pd.read_csv(aro_index, sep='\t')
  aro_df

  merged_df = amr_aro_df.merge(aro_df, how='left', on = "ARO Accession")
  amr_df.index = merged_df['AMR Gene Family']
  amr_df = amr_df.drop("ARO Accession", axis=1)
  amr_df['AMR_gene_family'] = amr_df.index

  amr_df = amr_df.groupby('AMR_gene_family').sum()
  amr_df = amr_df.T
  amr_df.to_csv(out_dir+"/gene_family.csv")#("analysis-out/4-Select_major_ABR/gene_family.csv")
  amr_df


  # In[6]:


  # check significance;
  from scipy import stats

  major_gene_families=['resistance-nodulation-cell division (RND) antibiotic efflux pump',
  'major facilitator superfamily (MFS) antibiotic efflux pump',
  'TEM beta-lactamase',
  'ATP-binding cassette (ABC) antibiotic efflux pump',
  'pmr phosphoethanolamine transferase',
  'ABC-F ATP-binding cassette ribosomal protection protein',
  'ampC-type beta-lactamase',
  'tetracycline-resistant ribosomal protection protein',
  'OKP beta-lactamase',
  'ACT beta-lactamase']

  amr_df_tmp = amr_df[major_gene_families]
  amr_df_tmp['cohort'] = amr_df_tmp.index

  control_df = amr_df_tmp[amr_df_tmp['cohort']=='Control']
  ab_df = amr_df_tmp[amr_df_tmp['cohort']=='Antibiotic']
  for amr in amr_df_tmp.columns:
    U, p = stats.mannwhitneyu(control_df[amr], ab_df[amr])
    if p<0.05:
        print(amr, p)


  # In[7]:


  merged_df.to_csv(out_dir+"/gene_family_out.csv")#'analysis-out/4-Select_major_ABR/gene_family_out.csv')


  # In[8]:


  major_gene_families=['resistance-nodulation-cell division (RND) antibiotic efflux pump',
  'major facilitator superfamily (MFS) antibiotic efflux pump',
  'ATP-binding cassette (ABC) antibiotic efflux pump',
  'TEM beta-lactamase',
  'ampC-type beta-lactamase',
  'OKP beta-lactamase',
  'ACT beta-lactamase',                    
  'pmr phosphoethanolamine transferase',
  'ABC-F ATP-binding cassette ribosomal protection protein',
  'tetracycline-resistant ribosomal protection protein']


  # In[9]:


  amr_df.index


  # In[10]:


  #print(amr_df_tmp.columns)
  #amr_df_tmp = amr_df.drop(["cohort"], axis=1)
  new_df={"CPM":[], "Gene Family":[], 'Cohort':[]}
  for amr in amr_df_tmp.columns:
    if amr!='cohort':
        new_df["CPM"]+=list(amr_df_tmp[amr])
        new_df['Gene Family']+=[amr for x in range(len(amr_df_tmp))]
        new_df["Cohort"]+=list(amr_df_tmp.index)
    
  new_df = pd.DataFrame(new_df)
  new_df


  # In[11]:


  import seaborn as sns

  sns.set_theme()
  sns.set(rc={'figure.figsize':(7,11.7)})

  g = sns.barplot(x = 'CPM', y = 'Gene Family', hue = 'Cohort', data = new_df,
            palette = 'hls', order=major_gene_families)
  #g.set(xlim=(1, 3))
  g

