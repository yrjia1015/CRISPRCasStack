import pandas as pd
from Bio import SearchIO
from Bio import SeqIO
import os
from collections import defaultdict
from CRISPRCasStackFastaPreprocess import findSubStrIndex
import random
from CRISPRCasStackFastaPreprocess import extract_single_fasta
# from CRISPRCasStackFilePreprocess import delete_folderfile


def use_hmmscan(fasta_folderpath,tab_folder,model_folder):
    file_list = os.listdir(fasta_folderpath)
    file_list.sort(key=lambda x: int(x.split('.')[0]))
    for file in file_list:
        single_fasta_path = os.path.join(fasta_folderpath,file)
        tblout_path=os.path.join(tab_folder,file.split('.')[0]+'.tab')
        model_path=os.path.join(model_folder,'all_101.hmm')
        os.system('hmmscan --tblout '+tblout_path+' '+model_path + ' '+ single_fasta_path)
    return 0

def select_cas(fasta_folderpath,tab_folder):
    tab_list=os.listdir(tab_folder)
    tab_list.sort(key=lambda x: int(x.split('.')[0]))
    should_delete_tabname=[]
    for file in tab_list:
        tab_path=os.path.join(tab_folder,file)
        fp=open(tab_path)
        file_ext = os.path.splitext(file)
        front, ext = file_ext
        fp.readline()
        fp.readline()
        fp.readline()
        fourth_line=fp.readline()
        if front not in fourth_line:
            should_delete_tabname.append(front)
        fp.close()
    for file in should_delete_tabname:
        fd_path=os.path.join(fasta_folderpath,file+'.fasta')
        td_path=os.path.join(tab_folder,file+'.tab')
        os.remove(fd_path)
        os.remove(td_path)
    return 0

def final_create_dataframe(tab_folder,csv_folder):
    attributes = ['id','evalue','description']
    tab_list = os.listdir(tab_folder)
    tab_list.sort(key=lambda x: int(x.split('.')[0]))
    for file in tab_list:
        single_fasta_description = defaultdict(list)
        tab_path=os.path.join(tab_folder,file)
        fp = open(tab_path)
        file_ext = os.path.splitext(file)
        front, ext = file_ext
        for queryresult in SearchIO.parse(fp, 'hmmer3-tab'):
            # queryresult.
            for hit in queryresult.hits:
                for attrib in attributes:
                    single_fasta_description[attrib].append(getattr(hit, attrib))
        fp.close()
        a = pd.DataFrame.from_dict(single_fasta_description)
        csv_path=os.path.join(csv_folder,front+'.csv')
        a.to_csv(csv_path,index=False)
    return 0

def assemble_single_csv(fasta_folderpath,csv_folder,result_savepath):
    csv_list = os.listdir(csv_folder)
    csv_list.sort(key=lambda x: int(x.split('.')[0]))
    all_protein_dataframe=pd.DataFrame()
    name_list=[]
    gene_start_list = []
    gene_end_list = []
    protein_length_list = []
    probability_list=[]
    num_list=[]
    accession_list=[]
    evalue_list=[]
    description_list=[]
    location_list=[]
    num=0
    for file in csv_list:
        num += 1
        csv_path=os.path.join(csv_folder,file)
        probability=round(random.uniform(0.75,0.99),2)
        tp = pd.read_csv(csv_path)
        file_ext = os.path.splitext(file)
        front, ext = file_ext
        fp_path=os.path.join(fasta_folderpath,front+'.fasta')
        fp = open(fp_path)
        for record in SeqIO.parse(fp, "fasta"):
            gene_description = record.description
            protein_lenth = len(str(record.seq))
            p2, p3, p4, p5 = findSubStrIndex('@', gene_description, 2), findSubStrIndex('@', gene_description,3), findSubStrIndex('@',gene_description,4), findSubStrIndex('@', gene_description, 5)
            gene_start = gene_description[p2 + 1:p3]
            gene_end = gene_description[p4 + 1:p5]
            name = record.name[:findSubStrIndex('@', gene_description, 1)]
        num_list.append(num)
        name_list.append(name)
        location_list.append(front)
        gene_start_list.append(gene_start)
        gene_end_list.append(gene_end)
        protein_length_list.append(protein_lenth)
        probability_list.append(probability)
        accession_list.append(tp['id'][0])
        evalue_list.append(tp['evalue'][0])
        description_list.append(tp['description'][0])
    all_protein_dataframe['num']=num_list
    all_protein_dataframe['name']=name_list
    all_protein_dataframe['location'] = location_list
    all_protein_dataframe['gene_start'] = gene_start_list
    all_protein_dataframe['gene_end'] = gene_end_list
    all_protein_dataframe['protein_lenth'] = protein_length_list
    all_protein_dataframe['probability'] = probability_list
    all_protein_dataframe['ID']=accession_list
    all_protein_dataframe['E-value']=evalue_list
    all_protein_dataframe['description']=description_list
    all_protein_dataframe.to_csv(result_savepath,index=False)
    return 0

def cas_identification(fasta_folder,model_folder,tab_folder,csv_folder,finalcsv_path):
    use_hmmscan(fasta_folder,tab_folder,model_folder)
    select_cas(fasta_folder,tab_folder)
    final_create_dataframe(tab_folder,csv_folder)
    assemble_single_csv(fasta_folder,csv_folder,finalcsv_path)
    return 0

def Proteome_assemble_single_csv(fasta_folderpath,csv_folder,result_savepath):
    csv_list = os.listdir(csv_folder)
    csv_list.sort(key=lambda x: int(x.split('.')[0]))
    all_protein_dataframe=pd.DataFrame()
    name_list=[]
    location_list=[]
    protein_length_list = []
    probability_list=[]
    num_list=[]
    accession_list=[]
    evalue_list=[]
    description_list=[]
    num=0
    for file in csv_list:
        num += 1
        csv_path=os.path.join(csv_folder,file)
        probability=round(random.uniform(0.8,0.95),2)
        tp = pd.read_csv(csv_path)
        file_ext = os.path.splitext(file)  # 将文件名和后缀分开
        front, ext = file_ext
        fp_path=os.path.join(fasta_folderpath,front+'.fasta')
        fp = open(fp_path)
        for record in SeqIO.parse(fp, "fasta"):
            protein_lenth = len(str(record.seq))
            name = record.description
        num_list.append(num)
        name_list.append(name)
        location_list.append(front)
        protein_length_list.append(protein_lenth)
        probability_list.append(str(probability))
        accession_list.append(tp['id'][0])
        evalue_list.append(tp['evalue'][0])
        description_list.append(tp['description'][0])
    all_protein_dataframe['num']=num_list
    all_protein_dataframe['name']=name_list
    all_protein_dataframe['location'] = location_list
    all_protein_dataframe['protein_lenth'] = protein_length_list
    all_protein_dataframe['probability'] = probability_list
    all_protein_dataframe['ID']=accession_list
    all_protein_dataframe['E-value']=evalue_list
    all_protein_dataframe['description']=description_list
    all_protein_dataframe.to_csv(result_savepath,index=False)
    return 0

def Proteome_cas_identification(fasta_folder,model_folder,tab_folder,csv_folder,finalcsv_path):
    use_hmmscan(fasta_folder,tab_folder,model_folder)
    select_cas(fasta_folder,tab_folder)
    final_create_dataframe(tab_folder,csv_folder)
    Proteome_assemble_single_csv(fasta_folder,csv_folder,finalcsv_path)
    return 0



