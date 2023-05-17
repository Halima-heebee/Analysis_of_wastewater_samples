import os
import pandas as pd
from pathlib import Path

#_removedbiasedSNPs.tsv

directory = '/home/mbxha18/allcombinations/kit14/10.4_SUP_PROM_VARSCAN/ill_match'

All_control = pd.read_csv(f'/home/mbxha18/OvsI_prom/ONT_match_prom/All_control.tsv',sep="\t")
a_conc = pd.read_csv("/home/mbxha18/allcombinations/TWIST_synthetic_mixes_sample_sheet_150222_edited.txt", sep="\t",index_col=0)

unfiltered=[]
all_subdirs = [d for d in os.listdir(f'{directory}') if os.path.isdir(d)]
for dir in all_subdirs: 
    print (f"{dir}")
    os.chdir(f'{directory}/{dir}')
    sample_ont = f"{dir}_ont"
    sample_ill = f"{dir}_ill"
    
    ont = pd.read_csv(f"{directory}/{dir}/{dir}_removedbiasedSNPs.tsv",sep="\t")
    av_cov_ont = pd.read_csv(f"{directory}/{dir}/{dir}.ont.mosdepth.summary.txt",sep="\t")
    ont = ont.rename(columns={"VarAllele": "ALT", "Position": "POS", "Ref": "REF", "Reads1": "REF_DP", "Reads2": "ALT_DP", "VarFreq" : "ALT_F", "Chrom" : "REGION"})
    ont['ALT'] = ont['ALT'].map(lambda x: x.lstrip('+-'))
    ont["ALT_F"] = ont["ALT_F"].str.rstrip('%').astype('float')
    ont = ont.sort_values('ALT_F').drop_duplicates(subset='POS', keep='last', ignore_index= True).sort_values('POS').reset_index(drop=True)
    ont = ont[["POS","REF","ALT","REF_DP","ALT_DP", "ALT_F"]]
    
    if ont.empty == True:
        new_row = pd.DataFrame(data= {'POS':0, 'REF':0, 'ALT':0, 'REF_DP':0, "ALT_DP":0, "ALT_F":0 }, index=[0])
        ont = pd.concat([ont, new_row]).reset_index(drop=True)

    ont_depth = pd.read_csv(f"{directory}/{dir}/{dir}.ont.samtools.cov",sep="\t",
                                    header = None, names = ["REGION","POS","DEPTH"], usecols = ["POS","DEPTH"])

    try: 
        for index, row in a_conc.iterrows():
            if (row["plate_well"]  == dir.split('-')[0]):
                print(row["lineage"])
                ont["lineage"] = row["lineage"]
                sample_lineage = ont.lineage.values.tolist()[0]
    except IndexError:
        pass

    control_for_sample = All_control[All_control["lineage"].eq(sample_lineage)].drop_duplicates(ignore_index= True)

    ont_true = pd.merge(control_for_sample, ont, how= "left", on=["POS","REF","ALT","lineage"])
    ont_true_depth = pd.merge(ont_true, ont_depth, how= "left", on=["POS"])
    ont_table = ont_true_depth.fillna(0)

    def categorise(row):  
        if row['lineage'] == "neg/blank":
            return 'neg'
        if row['ALT_DP'] == 0 and row["DEPTH"] <= 7: 
            return 'no_coverage'
        elif row['ALT_DP']== 0 and row["DEPTH"] >= 8 :
            return 'not_called'
        else: 
            return 'called'
    ont_table["REASON"] = ont_table.apply(lambda row: categorise(row), axis=1)

    A = ont_table[["POS","REF","ALT","REF_DP","ALT_DP","ALT_F","DEPTH","REASON","lineage"]]
    
    ill = pd.read_csv(f"{directory}/{dir}/{dir}.ill.mincov8.minvar0.03.varscan.tsv",sep="\t")
    av_cov_ill = pd.read_csv(f"{directory}/{dir}/{dir}.ill.mosdepth.summary.txt",sep="\t")
    ill = ill.rename(columns={"VarAllele": "ALT", "Position": "POS", "Ref": "REF", "Reads1": "REF_DP", "Reads2": "ALT_DP", "VarFreq" : "ALT_F"})
    ill['ALT'] = ill['ALT'].map(lambda x: x.lstrip('+-'))
    ill["ALT_F"] = ill["ALT_F"].str.rstrip('%').astype('float')
    ill = ill.sort_values('ALT_F').drop_duplicates(subset='POS', keep='last', ignore_index= True).sort_values('POS').reset_index(drop=True)
    ill = ill[["POS","REF","ALT","REF_DP","ALT_DP", "ALT_F"]]
    
    if ill.empty == True:
        new_row = pd.DataFrame(data= {'POS':0, 'REF':0, 'ALT':0, 'REF_DP':0, "ALT_DP":0, "ALT_F":0 }, index=[0])
        ill = pd.concat([ill, new_row]).reset_index(drop=True)

    ill_depth = pd.read_csv(f"{directory}/{dir}/{dir}.ill.samtools.cov",sep="\t", 
                            header = None, names = ["REGION","POS","DEPTH"], usecols = ["POS","DEPTH"])

    try: 
        for index, row in a_conc.iterrows():
            if (row["plate_well"]  == dir.split('-')[0]):
                ill["lineage"] = row["lineage"]
                sample_lineage = ill.lineage.values.tolist()[0]
    except IndexError:
        pass
    
    control_for_sample = All_control[All_control["lineage"].eq(sample_lineage)].drop_duplicates(ignore_index= True)
    

    ill_true = pd.merge(control_for_sample, ill, how= "left", on=["POS","REF","ALT","lineage"], suffixes=("_ill", "_true"))
    ill_true_depth = pd.merge(ill_true, ill_depth, how= "left", on=["POS"])
    ill_table = ill_true_depth.fillna(0)

    def categorise(row):  
        if row['lineage'] == "neg/blank":
            return 'neg'
        if row['ALT_DP'] == 0 and row["DEPTH"] <= 7: 
            return 'no_coverage'
        elif row['ALT_DP']== 0 and row["DEPTH"] >= 8:
            return 'not_called'
        else: 
            return 'called'
    ill_table["REASON"] = ill_table.apply(lambda row: categorise(row), axis=1)

    B = ill_table[["POS","REF","ALT","REF_DP","ALT_DP","ALT_F","DEPTH","REASON","lineage"]]
    
    All = pd.merge(A, B, how="outer", on=["POS","REF","ALT","lineage"], suffixes=("_ont", "_ill"))
    All.to_csv(f'{directory}/{dir}/{dir}.variant.br.mincov8.minvar0.03.varscan.tsv',sep="\t",index = False)
    
    
    ont_FP = ont[~ont['POS'].isin(control_for_sample['POS'])]
    ill_FP = ill[~ill['POS'].isin(control_for_sample['POS'])]
    
    
    ont_only = ont_FP[~ont_FP['POS'].isin(ill_FP['POS'])]
    ill_only = ill_FP[~ill_FP['POS'].isin(ont_FP['POS'])]
    

    ont_only.to_csv(f'{directory}/{dir}/{dir}.FP_ont_only',sep="\t",index = False)
    ill_only.to_csv(f'{directory}/{dir}/{dir}.FP_ill_only',sep="\t",index = False)
    
    
    #Using the truth set to calculate probability
    try:
        covg8I = All[All["DEPTH_ill"].gt(7)]
        covg8O = All[All["DEPTH_ont"].gt(7)]
        IP = len(covg8I)/len(All)
        OP = len(covg8O)/len(All)
    
    #Using all positions in the genome
        ont_depthg8 = ont_depth[ont_depth["DEPTH"].gt(7)]
        ill_depthg8 = ill_depth[ill_depth["DEPTH"].gt(7)]

        OPAll= len(ont_depthg8)/len(ont_depth)
        IPAll= len(ill_depthg8)/len(ill_depth)
    except ZeroDivisionError:
        OPAll=0
    #To calculate Sensitivity_WC
    All_Cov = All[~All['REASON_ont'].isin(["no_coverage"]) ]# df[~df["column"].isin(["value"])]
    All_Cov2 = All_Cov[~All_Cov['REASON_ill'].isin(["no_coverage"]) ]


        # TP = candidates identified by truth and ONT as true
        # FP = candidates identified by ONT as true but not true
        # FN = candidates not identified by ONT as true
        # TN = candidates not identified by both truth and ONT
        # FNwc= candidates not identified by ONT and has coverage
    for index, row in ont_table.iterrows():
        if (row["lineage"]  == "neg"):
            True_variant = 0
            TP_ont = 0
            FP_ont = 0
            FN_ont = 0
            TN_ont = 0
            FNwc_ont = 0
        else:
            True_variant = len(control_for_sample["POS"]) 
            TP_ont = len(control_for_sample["POS"]) - len(A[A["ALT_DP"].eq(0)])
            FP_ont = len(ont["POS"].unique()) - len(A[A["ALT_DP"].ne(0)])
            FN_ont = True_variant - len(A[A["ALT_DP"].ne(0)])
            TN_ont = 29903 - len(ont["POS"].unique()) - FN_ont
            FNwc_ont = len(All_Cov2[All_Cov2["REASON_ont"].eq("not_called")]) 



    try: 
        Sensitivity_ONT = TP_ont/(TP_ont + FN_ont)
        Precision_ONT = TP_ont/(TP_ont + FP_ont)
        Jaccard_similarity_ONT = TP_ont/(TP_ont + FP_ont + FN_ont)

        Sensitivity_ONT_WC = len(All_Cov2[All_Cov2["REASON_ont"].eq("called")])/len(All_Cov2)
        Jaccard_similarity_ONT_WC = len(All_Cov2[All_Cov2["REASON_ont"].eq("called")])/(len(All_Cov2) + FNwc_ont)

    except ZeroDivisionError:
        Precision_ONT = 0
        Sensitivity_ONT = 0
        Jaccard_similarity_ONT = 0
        Sensitivity_ONT_WC = 0
        Jaccard_similarity_ONT_WC = 0

    for index, row in ont_table.iterrows():
        if (row["lineage"]  == "neg"):
            True_variant = 0
            TP_ill = 0
            FP_ill = 0
            FN_ill = 0
            TN_ill = 0
            FNwc_ill = 0
        else:
            True_variant = len(control_for_sample["POS"])
            TP_ill = len(control_for_sample["POS"]) - len(B[B["ALT_DP"].eq(0)])
            FP_ill = len(ill["POS"].unique()) - len(B[B["ALT_DP"].ne(0)])
            FN_ill = True_variant - len(B[B["ALT_DP"].ne(0)])
            TN_ill = 29903 - len(ill["POS"].unique()) - FN_ill
            FNwc_ill = len(All_Cov2[All_Cov2["REASON_ill"].eq("not_called")]) 
            

    try:
        Sensitivity_ILL = TP_ill/(TP_ill + FN_ill)
        Precision_ILL = TP_ill/(TP_ill + FP_ill)
        Jaccard_similarity_ILL = TP_ill/(TP_ill + FP_ill + FN_ill)

        Sensitivity_ILL_WC = len(All_Cov2[All_Cov2["REASON_ill"].eq("called")])/len(All_Cov2)
        Jaccard_similarity_ILL_WC = len(All_Cov2[All_Cov2["REASON_ill"].eq("called")])/(len(All_Cov2) + FNwc_ill)

    except ZeroDivisionError:
        Precision_ILL = 0
        Sensitivity_ILL = 0
        Jaccard_similarity_ILL = 0
        Sensitivity_ILL_WC = 0
        Jaccard_similarity_ILL_WC = 0


    for index, row in a_conc.iterrows():
        if (row["plate_well"]  == dir.split('-')[0]):
            alpha_conc = row["Alpha-C15"]
            alpha_exp_frequency = row ["alpha_exp_frequency"]
            beta_conc = row["Beta-C16"]
            beta_exp_frequency = row["beta_exp_frequency"]
            delta_conc = row["Delta-C23"]
            delta_exp_frequency = row["deltaC23_exp_frequency"]
            deltaAY2_conc = row["DeltaAY.2-C29"]
            deltaAY2_exp_frequency = row["deltaAY2_exp_frequency"]
            omicron_conc = row["Omicron-C48"]
            omicron_frequency = row["omicron_exp_frequency"]
            lineage = row["lineage"]
            
    d = {'sample': [sample_ont, sample_ill],'True_variant': [True_variant, True_variant], 'TP':[TP_ont,TP_ill], 'FP':[FP_ont, FP_ill], 'FN':[FN_ont,FN_ill], 'TN':[TN_ont,TN_ill], 'FN_wc':[FNwc_ont,FNwc_ill], 'sensitivity': [Sensitivity_ONT, Sensitivity_ILL], 'precision': [Precision_ONT, Precision_ILL], 'Jaccard_similarity': [Jaccard_similarity_ONT, Jaccard_similarity_ILL],'sensitivity_wc': [Sensitivity_ONT_WC, Sensitivity_ILL_WC], 'Jaccard_similarity_wc': [Jaccard_similarity_ONT_WC, Jaccard_similarity_ILL_WC], 'Truth_proportion': [OP, IP], 'All_proportion': [OPAll, IPAll],'lineage' : [lineage, lineage], 'alpha_conc': [alpha_conc, alpha_conc], 'beta_conc': [beta_conc,beta_conc], 'delta_conc': [delta_conc,delta_conc],'deltaAY.2_conc':[deltaAY2_conc,deltaAY2_conc], 'omicron_conc': [omicron_conc,omicron_conc], 'alpha_exp_frequency': [alpha_exp_frequency,alpha_exp_frequency],'beta_exp_frequency': [beta_exp_frequency,beta_exp_frequency], 'delta_exp_frequency':[delta_exp_frequency,delta_exp_frequency], 'deltaAY2_exp_frequency':[deltaAY2_exp_frequency,deltaAY2_exp_frequency], 'omicron_frequency' : [omicron_frequency,omicron_frequency], 'mean_cov': [av_cov_ont.loc[0]['mean'],av_cov_ill.loc[0]['mean']]}
    df = pd.DataFrame(data=d)
    unfiltered.append(df)

unfiltered_table = pd.concat(unfiltered).reset_index(drop=True)
print (unfiltered_table)
unfiltered_table.to_csv(f'{directory}/bigtable.mincov8.minvar0.03.br_kit14.10.4_SUP_PROM_VARSCAN.tsv',sep="\t") 














