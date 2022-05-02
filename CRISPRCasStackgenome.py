from Bio import SeqIO
import pandas as pd
import os
from CRISPRCasStackFastaPreprocess import findSubStrIndex
import random
import numpy as np
import joblib
from jax_unirep import get_reps


def start_CRISPRidentify(input_address,output_address):
    os.system(r'python CRISPRidentify.py '+'--file '+input_address+' --result_folder '+output_address)
    Complete_summary_path=os.path.join(output_address,'Complete_summary.csv')
    fp=pd.read_csv(Complete_summary_path)
    if 'Global ID' in fp.columns:
        new_csv=fp.drop('Global ID', 1)
    else:
        columns = ['Name','ID','Region index','Start','End','Length','Consensus repeat','Repeat Length',
                   'Average Spacer Length','Number of spacers','Strand','Category']
        new_csv=pd.DataFrame(columns=columns)
    result_path=os.path.join(output_address,'CRISPR_result.csv')
    new_csv.to_csv(result_path,index=False)
    return 0

def delete_file(output_address):
    delete_Complete_Cas_summarypath = os.path.join(output_address, 'Complete_Cas_summary.csv')
    os.system('rm -f '+delete_Complete_Cas_summarypath)
    delete_Complete_summarypath = os.path.join(output_address, 'Complete_summary.csv')
    os.system('rm -f '+delete_Complete_summarypath)
    return 0


def start_CRISPR(input_address,output_address):
    start_CRISPRidentify(input_address,output_address)
    delete_file(output_address)
    return 0

def gene_prediction(input_gene_address,output_protein_path,output_gene_path,output_gbk_path):
    os.system(r'prodigal -i '+input_gene_address+' -o '+output_gbk_path+' -a '+output_protein_path+' -d '+output_gene_path)
    return 0

def extract_info(singlefastapath):
    fp=open(singlefastapath)
    for record in SeqIO.parse(fp, "fasta"):
        gene_description = record.description
        protein_lenth=len(str(record.seq))
        p2,p3,p4,p5=findSubStrIndex('@',gene_description,2),findSubStrIndex('@',gene_description,3),findSubStrIndex('@',gene_description,4),findSubStrIndex('@',gene_description,5)
        gene_start=gene_description[p2+1:p3]
        gene_end=gene_description[p4+1:p5]
        name = record.name[:findSubStrIndex('@',gene_description,1)]
    fp.close()
    return name,gene_start,gene_end,protein_lenth

def judge_cas(fastafolderpath,psefolderAddress,aatpfolderAddress,trifolderAddress,csv_save_address,threshold):
    psefile_list = os.listdir(psefolderAddress)
    aatpfile_list = os.listdir(aatpfolderAddress)
    trifile_list = os.listdir(trifolderAddress)
    psefile_list.sort(key=lambda x: int(x.split('.')[0]))
    aatpfile_list.sort(key=lambda x: int(x.split('.')[0]))
    trifile_list.sort(key=lambda x: int(x.split('.')[0]))
    num = 0
    columns_name=['num','name','location','gene_start','gene_end','protein_lenth','probability','ID','E-value','description']
    protein_information_csv = pd.DataFrame(columns=columns_name)
    for file in psefile_list:
        if (file in aatpfile_list) and (file in trifile_list):
            file_ext = os.path.splitext(file)
            file_subname, file_type = file_ext
            pse_single_csv_address = os.path.join(psefolderAddress, file)
            aatp_single_csv_address = os.path.join(aatpfolderAddress, file)
            tri_single_csv_address = os.path.join(trifolderAddress, file)
            single_fasta_address = os.path.join(fastafolderpath, file_subname + '.fasta')
            psesingle_csv = pd.read_csv(pse_single_csv_address)
            psefeature = np.array(psesingle_csv.iloc[0:, :psesingle_csv.shape[1]])
            aatpsingle_csv = pd.read_csv(aatp_single_csv_address)
            aatpfeature = np.array(aatpsingle_csv.iloc[0:, :aatpsingle_csv.shape[1]])
            trisingle_csv = pd.read_csv(tri_single_csv_address)
            trifeature = np.array(trisingle_csv.iloc[0:, :trisingle_csv.shape[1]])
            fp = open(single_fasta_address)
            for record in SeqIO.parse(fp, "fasta"):
                seq = record.seq
            fp.close()
            h_avg = get_reps(seq)[0]
            new=pd.DataFrame({'UniRef_809':h_avg[0][808], 'aatp384':aatpfeature[0][384], 'UniRef_1164':h_avg[0][1163],
                              'tri_gram_pssm3629':trifeature[0][3629], 'aatp390':aatpfeature[0][390], 'aatp38':aatpfeature[0][38],
                              'pse_pssm4':psefeature[0][4], 'UniRef_1473':h_avg[0][1472], 'UniRef_15':h_avg[0][14],
                              'aatp397':aatpfeature[0][397], 'UniRef_308':h_avg[0][307], 'aatp158':aatpfeature[0][158],
                              'pse_pssm38':psefeature[0][38], 'aatp62':aatpfeature[0][62], 'UniRef_222':h_avg[0][221],
                              'aatp402':aatpfeature[0][402], 'tri_gram_pssm3998':trifeature[0][3998], 'tri_gram_pssm7218':trifeature[0][7218],
                              'aatp4':aatpfeature[0][4], 'UniRef_1553':h_avg[0][1552], 'UniRef_608':h_avg[0][607],
                              'tri_gram_pssm7490':trifeature[0][7490], 'aatp55':aatpfeature[0][55]}, index=[1])
            x_test = np.array(new.iloc[0:, :new.shape[1]])
            LGBMpath=os.path.join(os.path.split(os.path.realpath(__file__))[0],'model','LGBM.pkl')
            LGBN = joblib.load(LGBMpath)
            LGBM_feature=LGBN.predict_proba(x_test)
            RFpath=os.path.join(os.path.split(os.path.realpath(__file__))[0],'model','RF.pkl')
            RF = joblib.load(RFpath)
            RF_feature = RF.predict_proba(x_test)
            ERTpath = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'model', 'ERT.pkl')
            ERT = joblib.load(ERTpath)
            ERT_feature = ERT.predict_proba(x_test)
            GBDTpath = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'model', 'GBDT.pkl')
            GBDT = joblib.load(GBDTpath)
            GBDT_feature = GBDT.predict_proba(x_test)
            probability_feature = pd.DataFrame({'LGBM_0':LGBM_feature[0][0], 'LGBM_1':LGBM_feature[0][1], 'GBDT_0':GBDT_feature[0][0],
                              'GBDT_1':GBDT_feature[0][1], 'RF_0':RF_feature[0][0], 'RF_1':RF_feature[0][1],
                              'ET_0':ERT_feature[0][0], 'ET_1':ERT_feature[0][1]}, index=[1])
            CRISPRCasStackpath = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'model', 'CRISPRCasStack.pkl')
            CRISPRCasStack = joblib.load(CRISPRCasStackpath)
            x_test_probability = np.array(probability_feature.iloc[0:, :probability_feature.shape[1]])
            CRISPRCasStack_feature = CRISPRCasStack.predict_proba(x_test_probability)
            cas_probability = round(CRISPRCasStack_feature[0][1], 3)
            if cas_probability > threshold:
                num += 1
                name, gene_start, gene_end, protein_lenth = extract_info(single_fasta_address)
                location=file_subname
                new_simple = pd.DataFrame({'num': num,
                                            'name': name,
                                           'location': location,
                                            'gene_start': gene_start,
                                           'gene_end': gene_end,
                                           'protein_lenth': protein_lenth,
                                           'probability':cas_probability,
                                           'ID':'Null',
                                           'E-value':'Null',
                                           'description':'Null'}, index=[1])
                protein_information_csv = protein_information_csv.append(new_simple)
    protein_information_csv.to_csv(csv_save_address,index=False)
    return 0

def screening(mlcas_csv,hmmcas_csv,allfasta_path,finalresult_path):
    all_fasta=open(allfasta_path)
    fasta_name=[]
    num = 0
    num_list = []
    ml = pd.read_csv(mlcas_csv)
    ml_empty=ml.empty
    hmm = pd.read_csv(hmmcas_csv)
    hmm_empty = hmm.empty
    ml_namelist=list(ml['name'])
    hmm_namelist = list(hmm['name'])
    columns_name = ['num', 'name','location', 'gene_start', 'gene_end', 'protein_lenth', 'probability', 'ID', 'E-value',
                    'description']
    combine_mlhmm=pd.DataFrame(columns=columns_name)
    if ml_empty and hmm_empty:
        combine_mlhmm.to_csv(finalresult_path, index=False)
    else:
        for record in SeqIO.parse(all_fasta, "fasta"):
            name = record.name  #
            if '|' in name:
                name=name.replace('|','_')
            if '.' in name:
                name=name.replace('.','_')
            if ' ' in name:
                name=name.replace(' ','@')
            fasta_name.append(name)
        for fasta in fasta_name:
            if fasta in hmm_namelist:
                num+=1
                num_list.append(num)
                new_simple=hmm[hmm['name'].isin([fasta])]
                combine_mlhmm=combine_mlhmm.append(new_simple)
                continue
            if fasta in ml_namelist:
                num+=1
                num_list.append(num)
                new_simple = ml[ml['name'].isin([fasta])]
                combine_mlhmm=combine_mlhmm.append(new_simple)
                continue

        combine_mlhmm['num']=num_list
        combine_mlhmm.to_csv(finalresult_path,index=False)
    all_fasta.close()
    return 0