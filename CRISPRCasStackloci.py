# coding=utf-8 #
import pandas as pd

def classfication_casloci(single_loci_list):
    lower_single_loci_list = []
    for sample in single_loci_list:
        lower_single_loci_list.append(sample.lower())
    loci_type = 'Cas_loci'
    for gene in lower_single_loci_list:
        if 'cas3' in gene:
            loci_type = 'Cas_loci_typeI'
            break
        if 'cas10' in gene:
            loci_type = 'Cas_loci_typeIII'
            break
        if 'csf1' in gene:
            loci_type = 'Cas_loci_typeIV'
            break
        if 'cas9' in gene:
            loci_type = 'Cas_loci_typeII'
            break
        if 'cas12' in gene:
            loci_type = 'Cas_loci_typeV'
            break
        if 'cas13' in gene:
            loci_type = 'Cas_loci_typeVI'
            break
    if loci_type == 'Cas_loci_typeI':
        num_Ia = 0
        for gene in lower_single_loci_list:
            if 'cas8a' in gene:
                num_Ia += 1
                break
        for gene in lower_single_loci_list:
            if 'csa5' in gene:
                num_Ia += 1
                break
        if num_Ia == 2:
            loci_type = 'Cas_loci_typeI-A'
            return loci_type
        if len(lower_single_loci_list) > 5:
            for gene in lower_single_loci_list:
                if 'cas8b' in gene:
                    loci_type = 'Cas_loci_typeI-B'
                    return loci_type
        num_Ic = 0
        for gene in lower_single_loci_list:
            if 'cas5' in gene:
                num_Ic += 1
                break
        for gene in lower_single_loci_list:
            if 'cas8c' in gene:
                num_Ic += 1
                break
        if num_Ic == 2:
            loci_type = 'Cas_loci_typeI-C'
            return loci_type
        for gene in lower_single_loci_list:
            if 'cas10d' in gene:
                loci_type = 'Cas_loci_typeI-D'
                return loci_type
        num_Ie = 0
        for gene in lower_single_loci_list:
            if 'cse1' in gene:
                num_Ie += 1
                break
        for gene in lower_single_loci_list:
            if 'cse2' in gene:
                num_Ie += 1
                break
        for gene in lower_single_loci_list:
            if 'cas4' in gene:
                num_Ie -= 1
                break
        if num_Ie == 2:
            loci_type = 'Cas_loci_typeI-E'
            return loci_type
        num_If = 0
        for gene in lower_single_loci_list:
            if 'csy1' in gene:
                num_If += 1
                break
        for gene in lower_single_loci_list:
            if 'csy2' in gene:
                num_If += 1
                break
        for gene in lower_single_loci_list:
            if 'csy3' in gene:
                num_If += 1
                break
        for gene in lower_single_loci_list:
            if 'cas6' in gene:
                num_If += 1
                break
        if num_If == 4:
            loci_type = 'Cas_loci_typeI-F1'
            return loci_type
        return loci_type
    if loci_type == 'Cas_loci_typeIII':
        for gene in lower_single_loci_list:
            if 'csm2' in gene:
                loci_type = 'Cas_loci_typeIII-A'
                return loci_type
        for gene in lower_single_loci_list:
            if 'cmr5' in gene:
                loci_type = 'Cas_loci_typeIII-B'
                return loci_type
        for gene in lower_single_loci_list:
            if 'csx10' in gene:
                loci_type = 'Cas_loci_typeIII-D'
                return loci_type
        return loci_type
    if loci_type == 'Cas_loci_typeIV':
        for gene in lower_single_loci_list:
            if 'csf3' in gene:
                loci_type = 'Cas_loci_typeIV-C'
                return loci_type
        return loci_type
    if loci_type == 'Cas_loci_typeII':
        for gene in lower_single_loci_list:
            if 'csn2' in gene:
                loci_type = 'Cas_loci_typeII-A'
                return loci_type
        for gene in lower_single_loci_list:
            if 'cas4' in gene:
                loci_type = 'Cas_loci_typeII-B'
                return loci_type
        num_IIc = 0
        for gene in lower_single_loci_list:
            if 'cas1' in gene:
                num_IIc += 1
                break
        for gene in lower_single_loci_list:
            if 'cas2' in gene:
                num_IIc += 1
                break
        if num_IIc == 2:
            loci_type = 'Cas_loci_typeII-B'
            return loci_type
        return loci_type
    if loci_type == 'Cas_loci_typeV':
        for gene in lower_single_loci_list:
            if 'cas12a' in gene:
                loci_type = 'Cas_loci_typeV-A'
                return loci_type
        for gene in lower_single_loci_list:
            if 'cas12b' in gene:
                loci_type = 'Cas_loci_typeV-B'
                return loci_type
        for gene in lower_single_loci_list:
            if 'cas12c' in gene:
                loci_type = 'Cas_loci_typeV-C'
                return loci_type
        for gene in lower_single_loci_list:
            if 'cas12d' in gene:
                loci_type = 'Cas_loci_typeV-D'
                return loci_type
        for gene in lower_single_loci_list:
            if 'cas12e' in gene:
                loci_type = 'Cas_loci_typeV-E'
                return loci_type
        for gene in lower_single_loci_list:
            if 'cas12f' in gene:
                loci_type = 'Cas_loci_typeV-F'
                return loci_type
        for gene in lower_single_loci_list:
            if 'cas12g' in gene:
                loci_type = 'Cas_loci_typeV-G'
                return loci_type
        for gene in lower_single_loci_list:
            if 'cas12h' in gene:
                loci_type = 'Cas_loci_typeV-H'
                return loci_type
        for gene in lower_single_loci_list:
            if 'cas12i' in gene:
                loci_type = 'Cas_loci_typeV-I'
                return loci_type
        for gene in lower_single_loci_list:
            if 'cas12k' in gene:
                loci_type = 'Cas_loci_typeV-K'
                return loci_type
        return loci_type
    if loci_type == 'Cas_loci_typeVI':
        for gene in lower_single_loci_list:
            if 'cas13a' in gene:
                loci_type = 'Cas_loci_typeVI-A'
                return loci_type
        for gene in lower_single_loci_list:
            if 'cas13b' in gene:
                if 'B1' in gene:
                    loci_type = 'Cas_loci_typeVI-B1'
                    return loci_type
        for gene in lower_single_loci_list:
            if 'cas13b' in gene:
                if 'B2' in gene:
                    loci_type = 'Cas_loci_typeVI-B2'
                    return loci_type
        for gene in lower_single_loci_list:
            if 'cas13c' in gene:
                loci_type = 'Cas_loci_typeVI-C'
                return loci_type
        for gene in lower_single_loci_list:
            if 'cas13d' in gene:
                loci_type = 'Cas_loci_typeVI-D'
                return loci_type
        return loci_type
    return loci_type

def find_cas_loci(casprotein_csvpath,casloci_savepath):
    cas_protein_csv=pd.read_csv(casprotein_csvpath)
    cas_protein_csv_empty=cas_protein_csv.empty
    location=list(cas_protein_csv['location'])
    eachloci = []
    eachloci_3 = []
    loci_description = ['loci_ID', 'loci_type', 'loci_start', 'loci_end', 'each_gene']
    cas_loci_csv = pd.DataFrame(columns=loci_description)
    if cas_protein_csv_empty:
        cas_loci_csv.to_csv(casloci_savepath, index=False)
    else:
        for x in sorted(set(location)):
            eachloci.append(x)
            if x + 1 not in location:
                if len(eachloci) > 2:
                    # print(eachloci)
                    eachloci_3.append(eachloci)
                eachloci = []
        if len(eachloci_3) != 0:
            num = 0
            for loci in eachloci_3:
                num += 1
                loci_first_gene=cas_protein_csv[cas_protein_csv['location'].isin([loci[0]])]
                loci_start=loci_first_gene['gene_start'].values[0]
                loci_last_gene = cas_protein_csv[cas_protein_csv['location'].isin([loci[len(loci)-1]])]
                loci_end = loci_last_gene['gene_end'].values[0]
                loci_casgene_name=[]
                for cas_gene_location in loci:
                    cas_gene_infomation = cas_protein_csv[cas_protein_csv['location'].isin([cas_gene_location])]
                    description=cas_gene_infomation['ID'].values[0]
                    loci_casgene_name.append(description)
                loci_type=classfication_casloci(loci_casgene_name)
                each_gene = ' , '.join(loci_casgene_name)
                new_loci=pd.DataFrame({'loci_ID':num,
                                       'loci_type':loci_type,
                                       'loci_start':loci_start,
                                       'loci_end':loci_end,
                                       'each_gene':each_gene}, index=[1])
                cas_loci_csv = cas_loci_csv.append(new_loci)
            cas_loci_csv.to_csv(casloci_savepath,index=False)
        else:
            cas_loci_csv.to_csv(casloci_savepath, index=False)

    return 0

def detect_CRISPRCAS(casloci_csvpath,crispr_csvpath,CRISPRCAS_savepath):
    casloci=pd.read_csv(casloci_csvpath)
    casloci_empty = casloci.empty
    CRISPR=pd.read_csv(crispr_csvpath)
    CRISPR_empty = CRISPR.empty
    num = 0
    num_list = []
    CRISPRCAS_description = ['num', 'loci_type', 'loci_start', 'loci_end', 'each_gene',
                             'CRISPR_start','CRISPR_end','CRISPR_lengh','CRISPR_Consensus_repeat',
                             'Number_of_spacers','CRISPR_strand']
    CRISPRCAS_csv = pd.DataFrame(columns=CRISPRCAS_description)

    if CRISPR_empty or casloci_empty:
        CRISPRCAS_csv.to_csv(CRISPRCAS_savepath, index=False)
    else:
        for index,loci in casloci.iterrows():
            loci_start = loci['loci_start']
            loci_end = loci['loci_end']
            for index,crispr in CRISPR.iterrows():
                crispr_start = crispr['Start']
                crispr_end = crispr['End']
                if abs(crispr_start - loci_end) <= 5000 or abs(loci_start - crispr_end) <= 5000 :
                    num += 1
                    newCRISPRCAS=pd.DataFrame({'num':num,
                                               'loci_type':loci['loci_type'],
                                               'loci_start':loci_start,
                                               'loci_end':loci_end,
                                               'each_gene':loci['each_gene'],
                                               'CRISPR_start':crispr_start,
                                               'CRISPR_end':crispr_end,
                                               'CRISPR_lengh':crispr['Length'],
                                               'CRISPR_Consensus_repeat':crispr['Consensus repeat'],
                                               'Number_of_spacers':crispr['Number of spacers'],
                                               'CRISPR_strand':crispr['Strand']}, index=[1])
                    CRISPRCAS_csv = CRISPRCAS_csv.append(newCRISPRCAS)
                    # break
        CRISPRCAS_csv.to_csv(CRISPRCAS_savepath,index=False)
    return 0

