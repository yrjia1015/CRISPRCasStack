# CRISPRCasStack
CRISPRCasStack is a toolkit capable of accurately identifying Cas proteins and comprehensively predicting CRISPR-Cas-related components, with the goal of accurately identifying potential Cas proteins that cannot currently be identified based on homology through a machine learning-based approach. At the same time, CRISPRCasStack integrates a set of state-of-the-art tools that can accurately identify CRISPR-Cas locus on prokaryotic genomes. Two user input modes are provided in the CRISPRCasStack toolkit: Genome sequence input mode and Proteins sequence input mode. If you want to detect CRISPR-Cas locus related information, please enter the genome sequence. If you want to detect only Cas proteins, please enter the proteome sequence.

## Prerequisite

### System requirements

First you need to use CRISPRCasStack on the Linux system and in a conda environment.

### Install requirements.txt

```
conda install --yes --file requirements.txt
```

### Decompress the hmm file 

Due to file upload restrictions, you will need to extract the zip file from the hmm folder to the hmm folder.The database is then formatted in the hmm folder with the following command.

```
hmmpress AllProfiles.hmm
```

### Preparing the uniref50 database 

Due to file upload restrictions, you will need to download the uniref50(Download address:https://www.uniprot.org/downloads#unireflink) database to the database folder and makeblastdb it with the following command

```
makeblastdb -in uniref50.fasta -parse_seqids -hash_index -dbtype prot
```

## Using CRISPRCasStack

### help parameter

You can see the help by using the `-h` option

```

python CRISPRCasStack.py -h

```

### Mandatory parameters
#### Given an input file in fasta format ( single genome sequence or proteome sequence)
* `-i <Path of the input file>`
#### Give the folder where the predictions will be stored 
* `-o <Path of the output folder>`
#### Give the selected calculation mode (select `g` if the input is a genomic sequence, or `p` if the input is a protein sequence) 
* `-m <Path of the output folder>`

### Additional parameters
#### Threshold (default is 0.5)
* `-t <Threshold values for determining Cas proteins>`

### Example
#### Take test.txt and test_g.fasta in the testinput folder as an example
```

python CRISPRCasStack.py -i testinput/test.txt -o testoutput -m p

```

```

python CRISPRCasStack.py -i testinput/test_g.fasta -o testoutput -m g

```
