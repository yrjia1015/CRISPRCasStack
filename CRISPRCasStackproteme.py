import os
from Bio import SeqIO
import pandas as pd
import numpy as np
import joblib
import random
from jax_unirep import get_reps

def extract_single_fasta(fasta_address, singlefasta_folder):
    all_fasta = open(fasta_address)
    fasta_name = []
    file_name = []
    num = 0
    for record in SeqIO.parse(all_fasta, "fasta"):
        num += 1
        name = record.name
        description=record.description
        seq = str(record.seq)[:len(str(record.seq)) - 1]
        if '|' in description:
            description=description.replace('|','_')
        if '.' in description:
            description=description.replace('.','_')
        if ' ' in description:
            description=description.replace(' ','@')
        fasta_name.append(name)
        single_file_name = str(num) + ".fasta"
        file_name.append(single_file_name)
        path=os.path.join(singlefasta_folder,single_file_name)
        single_fasta = open(path, 'w')
        single_fasta.write('>' + description + '\n')
        single_fasta.write(seq)
        single_fasta.close()
    all_fasta.close()
    return 0

def screen(fastafolder_path,dbfolder_path,pssmfolder_path):
    fasta_list=os.listdir(fastafolder_path)
    fasta_list.sort(key=lambda x: int(x.split('.')[0]))
    for i in fasta_list:
        file_ext = os.path.splitext(i)
        file_subname, file_type = file_ext
        db_path=os.path.join(dbfolder_path,'uniref50.fasta')
        fasta_path=os.path.join(fastafolder_path,i)
        pssm_name=file_subname+'.pssm'
        pssm_path=os.path.join(pssmfolder_path,pssm_name)
        os.system("psiblast -in_msa "+fasta_path+ " -db "+db_path+" -num_iterations 3 -out_ascii_pssm "+pssm_path)
    return 0

def single_pssmfeature(fastapath,pssmfolder,outpsecsvpath,outaatpcsvpath,outtricsvpath):
    os.system('perl possum_standalone.pl -i '+fastapath+' -o '+outpsecsvpath+' -t pse_pssm -p '+pssmfolder+' -h T')
    os.system('perl possum_standalone.pl -i ' + fastapath + ' -o ' + outaatpcsvpath + ' -t aatp -p ' + pssmfolder + ' -h T')
    os.system('perl possum_standalone.pl -i ' + fastapath + ' -o ' + outtricsvpath + ' -t tri_gram_pssm -p ' + pssmfolder + ' -h T')
    return 0

def all_feature(fastafolder,pssmfolder,psecsvfolder,aatpcsvfolder,tricsvfolder):
    fasta_list=os.listdir(fastafolder)
    fasta_list.sort(key=lambda x: int(x.split('.')[0]))
    for file in fasta_list:
        file_ext = os.path.splitext(file)
        file_subname, file_type = file_ext
        fasta_path=os.path.join(fastafolder,file)
        csv_name=file_subname+'.csv'
        psecsv_path=os.path.join(psecsvfolder,csv_name)
        aatpcsv_path = os.path.join(aatpcsvfolder, csv_name)
        tricsv_path = os.path.join(tricsvfolder, csv_name)
        single_pssmfeature(fasta_path,pssmfolder,psecsv_path,aatpcsv_path,tricsv_path)
    return 0

def Proteome_extract_info(singlefastapath):
    fp = open(singlefastapath)
    for record in SeqIO.parse(fp, "fasta"):
        name = record.description
        protein_lenth = len(str(record.seq))
    fp.close()
    return name, protein_lenth

def Proteome_judge_cas(fastafolderpath,psefolderAddress,aatpfolderAddress,trifolderAddress,csv_save_address,threshold):
    psefile_list = os.listdir(psefolderAddress)
    aatpfile_list = os.listdir(aatpfolderAddress)
    trifile_list = os.listdir(trifolderAddress)
    psefile_list.sort(key=lambda x: int(x.split('.')[0]))
    aatpfile_list.sort(key=lambda x: int(x.split('.')[0]))
    trifile_list.sort(key=lambda x: int(x.split('.')[0]))
    num = 0
    columns_name=['num','name','location','protein_lenth','probability','ID','E-value','description']
    protein_information_csv = pd.DataFrame(columns=columns_name)

    for file in psefile_list:
        if (file in aatpfile_list) and (file in trifile_list):
            file_ext = os.path.splitext(file)  # 将文件名和后缀分开
            file_subname, file_type = file_ext  # 前缀，后缀
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
            x_test_probability = np.array(probability_feature.iloc[0:, :probability_feature.shape[1]])  # 每个样本的向量形式
            CRISPRCasStack_feature = CRISPRCasStack.predict_proba(x_test_probability)
            cas_probability = round(CRISPRCasStack_feature[0][1], 3)
            if cas_probability > threshold:
                location=file_subname
                num += 1
                name, protein_lenth = Proteome_extract_info(single_fasta_address)
                new_simple = pd.DataFrame({'num': num,
                                           'name': name,
                                           'location':location,
                                           'protein_lenth': protein_lenth,
                                           'probability':cas_probability,
                                           'ID':'Null',
                                           'E-value':'Null',
                                           'description':'Null'}, index=[1])
                protein_information_csv = protein_information_csv.append(new_simple)
    protein_information_csv.to_csv(csv_save_address,index=False)
    return 0

def Proteome_screening(mlcas_csv,hmmcas_csv,allfasta_path,finalresult_path):
    all_fasta=open(allfasta_path)
    fasta_name=[]
    num = 0
    num_list = []
    ml=pd.read_csv(mlcas_csv)
    ml_empty=ml.empty
    ml_namelist=list(ml['name'])
    hmm=pd.read_csv(hmmcas_csv)
    hmm_empty=hmm.empty
    hmm_namelist = list(hmm['name'])
    columns_name = ['num', 'name', 'location','protein_lenth', 'probability', 'ID', 'E-value',
                    'description']
    combine_mlhmm=pd.DataFrame(columns=columns_name)
    if ml_empty and hmm_empty:
        combine_mlhmm.to_csv(finalresult_path, index=False)
    else:
        for record in SeqIO.parse(all_fasta, "fasta"):
            name = record.name
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
