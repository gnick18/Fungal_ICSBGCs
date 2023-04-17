

# Reproducible Code

[More information on CHTC at UW-Madison](https://chtc.cs.wisc.edu/)

## Flowplan

| **Step** | What this step does                                          | Where it was ran                                       |
| -------- | ------------------------------------------------------------ | ------------------------------------------------------ |
| 1        | Obtain and annotate all of the genomes from NCBI (Ran on: Thu Dec 9 15:49:00 CST 2021) | CHTC supercomputer at UW-Madison                       |
| 2        | Running CASSIS on every ICS domain to obtain gene cluster predictions | CHTC supercomputer at UW-Madison                       |
| 3        | Creating gbk, fna, and gff files for each ICS BGC prediction | CHTC supercomputer at UW-Madison                       |
| 4        | Running BigScape, with a custom anchor domain script         | CHTC supercomputer at UW-Madison                       |
| 5        | Dereplicating the BigScape output                            | My local computer                                      |
| 6        | Creating metadata table needed for downstream processes      | My local computer                                      |
| 7        | Running clinker on dereplicated folders (GCFs, GCCs)         | CHTC supercomputer at UW-Madison                       |
| 8        | Investigating gene conservation and domain conservation      | CHTC supercomputer at UW-Madison and my local computer |
| Misc     | Expansion and Contraction analysis and scripts               | My local computer                                      |
| Misc     | Kendall's Tau Correlation Analysis                           | My local computer                                      |

*<u>Important note!</u>* Many of the scripts found in this reproducible-script, aka `.sub` and parts of `.sh`, are only pertentent to supercomputer users at UW-Madison (CHTC). They are copied here for reproducibility and transparency but were not designed to be run locally on your own personal computer. **To avoid the need to re-run the entire pipeline, all of the raw predictions, groupings, and meta-data tables are available in the supplemental data and publication GitHub repo! Additionally, we have designed a publically available [web-server](https://isocyanides.fungi.wisc.edu/) to aid in the exploration of all data generated with this publication** To make skimming this document easier, key lines of code with pertanenet parameters have been seperated and placed above the raw scripts for each step detailed in this document. If trying to replicate part or all of the code at your own University or Institutation, the raw scripts would need to be adapted or edited depending on your computer's specifications/dependencies. 

### Overall constant dependencies, versions, etc

- Python 3.7.9 was used throughout the entire project
- The Anaconda package manager was used to install and mainintain many of the related python and R-packages
- All programs run locally were done so on Mac-OS operating systems and thus are not guarenteed to work on the Windows operating systems!
- The CHTC supercomputer runs on a Linux distribution ([CentOS Stream 8](https://chtc.cs.wisc.edu/uw-research-computing/os-transition-htc)) 

## Step 1: Obtaining and annotating the fungal genomes

### List of dependent packages in conda:

```
# packages in environment at /home/gnickles/miniconda3/envs/default:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main  
_openmp_mutex             4.5                       1_gnu  
biopython                 1.78             py39h7f8727e_0  
blas                      1.0                         mkl  
bottleneck                1.3.2            py39hdd57654_1  
ca-certificates           2020.10.14                    0    anaconda
certifi                   2021.10.8        py39h06a4308_0  
intel-openmp              2021.4.0          h06a4308_3561  
ld_impl_linux-64          2.35.1               h7274673_9  
libffi                    3.3                  he6710b0_2  
libgcc-ng                 9.3.0               h5101ec6_17  
libgomp                   9.3.0               h5101ec6_17  
libstdcxx-ng              9.3.0               hd4cf53a_17  
mkl                       2021.4.0           h06a4308_640  
mkl-service               2.4.0            py39h7f8727e_0  
mkl_fft                   1.3.1            py39hd3c417c_0  
mkl_random                1.2.2            py39h51133e4_0  
ncurses                   6.3                  h7f8727e_2  
numexpr                   2.7.3            py39h22e1b3c_1  
numpy                     1.21.2           py39h20f2e39_0  
numpy-base                1.21.2           py39h79a1101_0  
openssl                   1.1.1l               h7f8727e_0  
pandas                    1.3.4            py39h8c16a72_0  
pip                       21.2.4           py39h06a4308_0  
python                    3.9.7                h12debd9_1  
python-dateutil           2.8.2              pyhd3eb1b0_0  
pytz                      2021.3             pyhd3eb1b0_0  
readline                  8.1                  h27cfd23_0  
setuptools                58.0.4           py39h06a4308_0  
six                       1.16.0             pyhd3eb1b0_0  
sqlite                    3.36.0               hc218d9a_0  
tk                        8.6.11               h1ccaba5_0  
tzdata                    2021e                hda174b7_0  
wheel                     0.37.0             pyhd3eb1b0_1  
xz                        5.2.5                h7b6447c_0  
zlib                      1.2.11               h7b6447c_3  
```

*All of the most recent versions as of Dec, 2021 were used*

### Condensed summary of the key programs and parameters 

1. Downloading the genomes using the ncbi datasets tool

```sh
#getting summary of available fungal genomes
datasets summary genome --annotated taxon fungi
#downloading the genomes using the datasets tool
datasets download genome accession
OR 
datasets download genome --annotated taxon fungi
```

2. Getting protein domain predictions for every genome

```sh
hmmpress Pfam-A.hmm
hmmsearch -E 1e-4 --cpu 1 --tblout "$Name"_hmmer.out Pfam-A.hmm "$Name"_prot.faa
```

3. Running Antismash

```sh
antismash --taxon fungi --genefinding-gff "$Name".gff "$Name".fna --output-dir "$Name"_anti --cpus 2
```

### Raw scripts

*Note!* All of these scripts are designed to work specifically on UW-Madison's supercomputer cluster and would fail if trying to run locally. They are copied here for reproducibility and transparency but were not designed to be run locally on your own personal computer. If trying to replicate part or all of the code, some editing will be needed before running it locally.

**ConvertSummaryToJson.sub**

```sub
# Specify the HTCondor Universe (vanilla is the default and is used
universe = vanilla
#
log = accessionSummary.log
error = accessionSummary.err
#
executable = /home/gnickles/FungalGenomes/Scripts/ConvertSummaryToJson.sh
#
#need to go in this order: [m/d/y aka 09/01/2020] (grab genomes released since this date, put nan if not using this filtering), mF to start from (1 is the defult if nothing provided)
arguments = nan 1
output = accessionSummary.out
#
#make sure it has gluster
requirements = (Target.HasCHTCStaging == true) && (OpSysMajorVer == 7)
#
# Specify that HTCondor should transfer files to and from the
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
#Transfering the files
transfer_input_files =/staging/gnickles/default.tar.gz, /home/gnickles/FungalGenomes/Scripts/ConvertSummaryToJson.py, /home/gnickles/FungalGenomes/Programs/datasets
transfer_output_files = AccessionOut.tar.gz
#
# Tell HTCondor what amount of compute resources
request_cpus = 1
request_memory = 10GB
request_disk = 10GB
#
#
# Getting the instances of the job
queue
```

**ConvertSummaryToJson.sh**

```sh
#!/bin/sh
#
#This line is super important otherwise the datasets will never run!!
unset http_proxy
# have job exit if any command returns with non-zero exit status (aka failure)
set -e
#
# replace env-name on the right hand side of this line with the name of your conda environment
ENVNAME=default
# if you need the environment directory to be named something other than the environment name, change this line
ENVDIR=$ENVNAME
#
# these lines handle setting up the environment; you shouldn't have to modify them
export PATH
mkdir $ENVDIR
tar -xzf $ENVNAME.tar.gz -C $ENVDIR
. $ENVDIR/bin/activate
#
#
outputDir=$(pwd)

#if there is an older accession file that is getting passed, put it as the 4th argument on this command
python3 ConvertSummaryToJson.py "$outputDir" "$1" "$2"

mkdir = AccessionOut
mv genomes.txt AccessionOut/
mv *_GeneomesMetaData.tsv AccessionOut/
mv all_access.txt AccessionOut/
mv *.json AccessionOut/
tar -czvf AccessionOut.tar.gz AccessionOut

```

**ConvertSummaryToJson.py**

```python
#! /usr/bin/env python
##RUN THIS PROGRAM FROM THE FOLDER YOU WANT EVERYTHING TO END UP IN, otherwise it will need some editing
#Questions or issues: email gnickles@wisc.edu

import os
from sys import argv
import json
import pandas as pd
from datetime import date
from Bio import Entrez
import ssl

######
# GLOBAL MF VARIABLE
######
today=date.today()
todaysDate = today.strftime("%m_%d_%y")
#bypassing the certificate verification step
ssl._create_default_https_context = ssl._create_unverified_context

#gets the json summary
def GetJsonSummary(outputDirectory, releasedSince):
    # Textual month, day and year
    print("Running the program and getting the json file:\n")
    print("Today's Date: ", todaysDate)
    #setting the output file for the os command
    outFile = ""
    cmd = ""
    if releasedSince.lower() == "nan":
        outFile = os.path.join(outputDirectory, "AllFungi.json")
        cmd = "./datasets summary genome --annotated taxon fungi > " + outFile
    else:
        alternateDate = "_".join(releasedSince.split('/'))
        outFile = os.path.join(outputDirectory, "FungiSince-" + alternateDate + ".json")
        cmd = "./datasets summary genome --annotated taxon fungi --released-since " + releasedSince + " > " + outFile
    os.system(cmd)
    return outFile

def ParseJson(oldAccession, outFile, mF_count):
    metaDF = pd.DataFrame(columns = ['Accession','tax_id','Contig_n50','seq_length','submission_date','submitter','assembly_level'])
    accessionDF = pd.DataFrame(columns = ["mF",'accession','assembly_level','gene_count','n50','strain',"taxonomy"])
    with open(outFile) as j:
        data = json.load(j)
        assemblies = data['assemblies']
        totalCount = data['total_count']
        print("Total number of genomes to parse: " + str(totalCount))
        numLeft = totalCount
        for assembly in assemblies:
            if numLeft % 100 == 0:
                print('Genomes left: ' + str(numLeft))
            #getting all the variables and stuff
            assemblyData = assembly['assembly']

            accession = assemblyData['assembly_accession'].strip()
            mfString = ""
            seqLength = ""
            subDate = ""
            assembly_level = ""
            taxID = taxID = ""
            n50 = ""
            submitter = ""
            geneCount = ""
            strain = ""
            #moving on if this accession is already in the all_access file
            try:
                if accession in oldAccession['accession']:
                    continue
                    numLeft -= 1
            except:
                pass
            try:
                seqLength = assemblyData['seq_length'].strip()
            except:
                seqLength = ""
            try:
                subDate = assemblyData['submission_date'].strip()
            except:
                subDate = ""
            try:
                assembly_level = assemblyData['assembly_level'].strip()
            except:
                assembly_level = ""
            try:
                taxID = assemblyData['org']['tax_id'].strip()
            except:
                taxID = taxID = ""
            try:
                n50 = assemblyData['contig_n50']
            except:
                n50 = ""
            try:
                submitter = assemblyData['submitter'].strip()
            except:
                submitter = ""
            try:
                geneCount = assemblyData['annotation_metadata']['stats']['gene_counts']['total']
            except:
                geneCount = ""
            try:
                strain = assemblyData['org']['strain'].strip()
            except:
                strain = ""
            try:
                mfString = GetMF(mF_count)
            except:
                mfString = ""
            mF_count += 1
            ####
            # Getting the taxonomy
            ####
            taxonomyString = GetTax(taxID)
            ####
            # Adding to dataframes
            ####
            newDataRow_Meta = {"Accession": accession,
            "tax_id": taxID,
            "Contig_n50": n50,
            "seq_length": seqLength,
            "submission_date": subDate,
            "submitter": submitter,
            "assembly_level": assembly_level
            }
            metaDF = metaDF.append(newDataRow_Meta,ignore_index=True)
            #adding the new data to the accession Table
            newDataRow_Acc = {"mF":mfString,
            'accession':accession,
            'assembly_level': assembly_level,
            'gene_count': geneCount,
            'n50':n50,
            'strain':strain,
            'taxonomy':taxonomyString
            }
            accessionDF = accessionDF.append(newDataRow_Acc, ignore_index=True)
            numLeft -= 1
    return [metaDF,accessionDF]

def GetMF(mF_count):
    mfString = str(mF_count) + "mF"
    return mfString

#Gets the taxonomy data using the entrez tool from ncbi
def GetTax(taxID):
    #Running the Entrez search
    search = Entrez.efetch(id = taxID, db = "taxonomy", retmode = "xml")
    results = Entrez.read(search)[0]
    linEx = results['LineageEx']

    #Filling the stuff into the correct column
    kingdom = ''
    phylum=''
    classs = ''
    order=''
    family=''
    genus=''
    species=''
    genome_info=results['ScientificName']
    for l in linEx:
        rank = l['Rank']
        if rank.lower() == 'kingdom':
            kingdom = l['ScientificName']
        elif rank.lower() == 'phylum':
            phylum = l['ScientificName']
        elif rank.lower() == 'class':
            classs = l['ScientificName']
        elif rank.lower() == 'order':
            order = l['ScientificName']
        elif rank.lower() == 'family':
            family = l['ScientificName']
        elif rank.lower() == 'genus':
            genus = l['ScientificName']
        elif rank.lower() == 'species':
            species = l['ScientificName']
    #converting everything to a nice string for the all_access file
    newTax = [kingdom,phylum,classs,order,family,genus,species,genome_info]
    newTax2 = []
    for n in newTax:
        n2 = n.strip()
        newTax2.append(n2)
    returnThis = ",".join(newTax2)
    return returnThis

if __name__ == '__main__':
    #change if you're not Grant! But honestly it's just a formality and is only for the NCBI's data collection
    Entrez.email = "gnickles@wisc.edu"
    #put the full path!!!
    outputDirectory = argv[1]
    #use the following format: [m/d/y] aka 09/01/2020 else put 'nan'
    releasedSince = argv[2]
    mF_count = int(argv[3])
    oldAccession = ""
    if len(argv) > 4:
        oldAccessionFile = argv[4]
        oldAccession = pd.read_csv(oldAccessionFile, sep='\t',names=["mF",'accession','assembly_level','gene_count','n50','strain',"taxonomy"])
    #gets the json using datasets
    outFile = GetJsonSummary(outputDirectory,releasedSince)
    #makes the table results
    results = ParseJson(oldAccession, outFile, mF_count)
    metaDF = results[0]
    accessionDF = results[1]
    #writing out the files
    metaOutPath = os.path.join(outputDirectory, todaysDate + "_GeneomesMetaData.tsv")
    accOutPath = os.path.join(outputDirectory, "all_access.txt")
    metaDF.to_csv(metaOutPath, sep='\t',header=True,index=False)
    accessionDF.to_csv(accOutPath,sep='\t',header=True,index=False)

    genomesList = []
    for index, row in accessionDF.iterrows():
        try:
            genomesList.append(row['accession'].strip())
        except:
            print(row + " didn't work!")
    #writing out the queue file
    with open(os.path.join(outputDirectory,"genomes.txt"),'w+') as genomesTxt:
        for g in genomesList:
            # if g == genomesList[0]:
            #     print(g)
            #     genomesTxt.write(g)
            # else:
            print(g)
            genomesTxt.write(g + "\n")
```

**ncbi-download.sub**

```
#Specify the HTCondor Universe (vanilla is the default and is used
universe = vanilla
#
arguments = $(genome)
#
#Specify your executable (single binary or a script that runs several
executable = /home/gnickles/FungalGenomes/Scripts/ncbi-download.sh
#
output = $(genome)_ncbi.out
log = $(genome)_ncbi.log
error = $(genome)_ncbi.err
#
#make sure it has gluster
requirements = (OpSysMajorVer =?= 7)
#
#Specify that HTCondor should transfer files to and from the
#should_transfer_files = YES
#when_to_transfer_output = ON_EXIT
transfer_input_files = /home/gnickles/FungalGenomes/Programs/datasets, /home/gnickles/FungalGenomes/all_access.txt
transfer_output_files = $(genome).tar.gz
#
#Tell HTCondor what amount of compute resources
request_cpus = 1
request_memory = 1GB
request_disk = 5GB
max_idle = 2000
#
#
#In order for the code to work without edits the queue file should always be named genomes.txt and be placed in the QueueFile folder
queue genome from /home/gnickles/FungalGenomes/QueueFile/genomes.txt
```

**ncbi-download.sh**

```sh
#!/bin/sh
#
baseDir=$(pwd)
#
#This line is super important otherwise the program will never run!!
unset http_proxy
#
accession="$1"
#making the output directory
mkdir "$accession"
#setting the mF number
#getting the MF
num=$(cat *all_access.txt | grep "$accession" | awk '{ print $1 }')

NewName="$num"-Z
#downloading the genome; this needs a -g for non-linux machines
./datasets download genome accession "$accession"
unzip ncbi_dataset.zip
#moving the gff
mv ./ncbi_dataset/data/"$accession"/genomic.gff ./"$accession".gff
#these won't alwasy exist but we want to remove them just in case
rm ./ncbi_dataset/data/"$accession"/rna.fna
rm ./ncbi_dataset/data/"$accession"/chrMT.fna
rm ./ncbi_dataset/data/"$accession"/cds_from_genomic.fna
#
#Finding the names of fna file for debugging purposes
for file in ./*.fna; do
  echo $file > fnaFileNames.txt
done
#
#fixing the chromosomes
cat ./ncbi_dataset/data/"$accession"/*.fna   | grep ">" | awk -v var="$NewName" '{print $1"\t"var""NR}' | sed 's/>//' > ./"$accession"/names
cat ./ncbi_dataset/data/"$accession"/*.fna  | awk  '/^>/{print ">zzz" ++i; next}{print}' | sed "s/zzz/$NewName/" > ./"$accession"/"$accession".fna
#renaming the gff file
awk -v OFS="\t" -F '\t' 'NR==FNR{a[$1]=$2; next}{$1=a[$1]; print}' ./"$accession"/names "$accession".gff | grep . > zemp
mv zemp ./"$accession"/"$accession".gff
#removing the older version
rm ./"$accession".gff
#removing uneccesary info on the protein fasta file so it is left with just the protein name
awk '{if ($1 ~ /^>/) {print $1} else {print}}' ./ncbi_dataset/data/"$accession"/protein.faa > ./"$accession"/"$accession"_prot.faa
# moving over the sequence report
mv ./ncbi_dataset/data/"$accession"/*.jsonl ./"$accession"/"$accession"_sequence_report.jsonl
#compressing the folder and some files
gzip ./"$accession"/"$accession".fna
gzip ./"$accession"/"$accession".gff
gzip ./"$accession"/"$accession"_prot.faa
#
tar -C "$baseDir" -czvf "$accession".tar.gz ./"$accession"
```

**hmmer.sub**

```
# Specify the HTCondor Universe (vanilla is the default)
universe = vanilla
#
log = $(name)_hmmer.log
error = $(name)_hmmer.err
#
# Specify your executable (single binary or a script that runs several
executable = hmmer.sh
arguments = $(name)
#the hmmer outfiles are hugeeeeeeeeeee
#output = $(NewProcess)_hmmer.out
#
#make sure it has gluster
requirements = (Target.HasCHTCStaging == true) && (OpSysMajorVer == 7)
#
# Specify that HTCondor should transfer files to and from the
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files =/home/gnickles/FungalGenomes/Programs/hmmer-3.3.2_compiled.tar.gz,/home/gnickles/FungalGenomes/$(name).tar.gz
#
# Tell HTCondor what amount of compute resources
request_cpus = 1
request_memory = 3GB
request_disk = 5GB
#
queue name from /home/gnickles/FungalGenomes/QueueFile/genomes.txt
```

**hmmer.sh**

```sh
#!/bin/sh
#
cp /staging/gnickles/Pfam-A.hmm.gz .
#
tar -xzvf hmmer-3.3.2_compiled.tar.gz
tar -xzvf "$1".tar.gz
gunzip Pfam-A.hmm.gz
gunzip ./"$1"/"$1"_prot.faa.gz
#
export PATH=$(pwd)/hmmer-3.3.2/src:$PATH
#
hmmpress Pfam-A.hmm
hmmsearch -E 1e-4 --cpu 1 --tblout "$1"_hmmer.out Pfam-A.hmm ./"$1"/"$1"_prot.faa
#
#Grant: added the sed statement to remove any lines that have --- in them
cat "$1"_hmmer.out | sed '/^#/d' | awk 'NR>3 {print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | sort -u -k1,2> temp
echo -e "Target_Name\tDomain\tDomain_Accession\tE-value\tscore\tbias" > topLine
cat topLine temp > "$1"_hmmer.out
rm temp topLine
#
rm *.fa *.gz *.tsv *.txt *h3f *h3i *h3m *h3p *hmm
gzip "$1"_hmmer.out
```

**gio_ncbi.sub**

```
# Specify the HTCondor Universe (vanilla is the default and is used
universe = vanilla
#
log = $(name)_gio.log
error = $(name)_gio.err
#
# exectuables and whatnot
executable = gio_ncbi.sh
arguments = $(name)
output = $(name)_gio.out
#
#make sure it has gluster
requirements = (Target.HasCHTCStaging == true) && (OpSysMajorVer == 7)
#
#Specify that HTCondor should transfer files to and from the
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
#
transfer_input_files =/home/gnickles/FungalGenomes/$(name).tar.gz
transfer_output_files = $(name)_GIO_ncbi_parent.txt.gz
#
# Tell HTCondor what amount of compute resources
request_cpus = 1
request_memory = 1GB
request_disk = 5GB
#
#
queue name from /home/gnickles/FungalGenomes/QueueFile/genomes.txt
```

**gio_ncbi.sh**

```sh
#!/bin/sh
#
#extracting the folder and moving everything to the base dir
tar -xzf "$1".tar.gz
mv ./"$1"/* ./
#decompressing everything
gunzip "$1"_prot.faa.gz
gunzip "$1".gff.gz
gunzip "$1".fna.gz
#
#grab out all prot names
grep ">" "$1"_prot.faa | sed 's/>//' | awk '{print $1}'> tempy.txt
#iterate through each prot name
while read -r REPLY;do
#the note and orig_prot id were causing issues for a very small number of isoaltes so i removed them. they have not been super well
#tested yet though
grep -w "$REPLY" "$1".gff | sed 's/Note=.*;//'| sed 's/orig_protein_id.*;//'| grep -w $REPLY > zempy
#
Parent_id=`grep -w "$REPLY" zempy | grep Parent |grep -v Parent="$REPLY"\; | sed 's/.*Parent=//' | sed 's/;.*//' | sort -u | tr '\n' "," | sed 's/,$//'`
#the above catestrophically breaks the loop if parent id ends the loop as a blank variable.. no idea why.
if [[ $Parent_id == "" ]];then
Parent_id=`echo $REPLY`
fi
#
grep -w $Parent_id "$1".gff >> zempy
cat zempy | sort -u > temp
mv temp zempy
#
echo $REPLY > search_terms
echo $Parent_id >>search_terms
#
#use both the parent id and the protein name to find lines with the corresponding gff entries
#and then process them to figure out where the earliest start is and latest end is for that gene's entries
grep -w -f search_terms zempy| awk -v var="$REPLY" -v parent=$Parent_id -F '\t' '{print $1"\t"var"\t"parent"\t"$4}' | sort -n -k4,4 > start_temp
grep -w -f search_terms zempy | awk -F '\t' '{print $5"\t"$7"\t"$9}'| sort -n -r -k1,1 >end_temp
start=`head -1 start_temp`
end=`head -1 end_temp`
echo -e "$start""\t""$end" >>  "$1"_GIO_ncbi_parent.txt
#
done < tempy.txt
# changed the nr to an n to prevent reverse sorting 
cat  "$1"_GIO_ncbi_parent.txt | sort -k1,1 -k4,4n >temp
mv temp  "$1"_GIO_ncbi_parent.txt
gzip "$1"_GIO_ncbi_parent.txt
#
# rm *fna *.gff *.fa *json tempy.txt zempy *.tar.gz search_terms end_temp start_temp
```

**Antismash.sub**

```sub
# Specify the HTCondor Universe (vanilla is the default and is used
universe = docker
docker_image=aflatoxing/antismash5
#
arguments = $(name)
#
# Specify your executable (single binary or a script that runs several
executable = Antismash.sh
#
output = $(name)_anti.out
log = $(name)_anti.log
error = $(name)_anti.err
#
#make sure it has gluster
requirements = (Target.HasCHTCStaging == true) && (OpSysMajorVer == 7)
#
# Specify that HTCondor should transfer files to and from the
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = /home/groups/keller_group/FungalGenomes/$(name).tar.gz
transfer_output_files = $(name)_anti.tar.gz
#
#Tell HTCondor what amount of compute resources
request_cpus = 2
request_memory = 15GB
request_disk = 10GB
#
queue name from /home/groups/keller_group/FungalGenomes/QueueFile/genomes.txt
```

**Antismash.sh**

```sh
#!/bin/bash
#
export PATH=/opt/conda/bin:$PATH
#
Name="$1"
#extracting the folder and moving everything to the base dir
tar -xzf "$Name".tar.gz
mv ./"$Name"/* ./
#decompressing everything
gunzip "$Name"_prot.faa.gz
gunzip "$Name".gff.gz
gunzip "$Name".fna.gz
#Run antismash
#
. activate antismash
#
antismash --taxon fungi --genefinding-gff "$Name".gff "$Name".fna --output-dir "$Name"_anti --cpus 2
#
#tidy up the antismash output
rm ./"$Name"_anti/"$Name".gbk
rm ./"$Name"_anti/*.zip
# rm ./"$Name"_anti/*.json
gzip ./"$Name"_anti/regions.js
gzip ./"$Name"_anti/"$Name".json
#
tar -czvf "$Name"_anti.tar.gz ./"$Name"_anti
#
```





## Step 2: Running CASSIS on the ICS containing genes

There are three inputs for CASSIS. The following is from the authors of the algorithm:

https://sbi.hki-jena.de/cassis/Help.php

*Date that I grabbed this table from the website: 12/15/2021*

**CASSIS requires two input files and an anchor gene name (not the protein name)**

| Input                  |                                                              |
| ---------------------- | ------------------------------------------------------------ |
| Genomic sequence file  | A multiFASTA file (`.fasta`) containing the DNA sequences of all contigs of the species. Contig numbers have to coincide with the annotation file. |
| Genome annotation file | A simple text file (no Excel!) in tabular format with at least five columns. These are: `gene <string> | contig <string> | start position <int> | stop position <int> | strand <+ or->` The column separator must be tab (`\t`). The annotation file must be sorted ascending by contig, start and stop. Start positions must be smaller than stop positions. Contig names have to coincide with the genome sequence file. |
| Anchor Gene            | Feature ID of the clusters anchor gene. The ID has to coincide with (i.e. be part of) the annotation file. This gene will be the starting point of the cluster prediction. |

There are a lot of dependencies needed to run CASSIS, and it is not friendly to the MacOS as the exectuable is a binary file. To get this working properly I created a Docker container that has all of the neccesary dependencies and programs to get it running. As long as the user inputs the correct files this will work normally. It is also publicly available.

[Link to the docker container](https://hub.docker.com/r/gnick18/cassis2)

Additionally, because my code adapts NCBI files to the CASSIS friendly format it can be leverage for other use cases with NCBI annotated genomes. You just need to generate the GIO file with the script that was included in **Step 1.**

### Condensed summary of parameters used with CASSIS

```sh
cassis --dir CASSIS_"$genome"_"$counter" --annotation "$genome".gff --genome "$genome".fna --meme 1.0e+005 --anchor "$anchor" --cluster "$anchor" --fimo 0.00006 --frequency 14 --gap-length 2 --prediction --num-cpus 1 --verbose
```

### Raw scripts

**GetAnchor_RunCassis.sub**

```
# Specify the HTCondor Universe (vanilla is the default)
universe = docker
docker_image=gnick18/cassis2
# docker_image=antismash/standalone-lite
#
log = $(genome)_anchor.log
error = $(genome)_anchor.err
#
# Specify your executable (single binary or a script that runs several
executable = GetAnchor_RunCassis.sh
arguments = $(genome)
output = $(genome)_anchor.out
#
#make sure it has gluster
requirements = (Target.HasCHTCStaging == true) && (OpSysMajorVer == 7)
#
# Specify that HTCondor should transfer files to and from the
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files =/home/gnickles/FungalGenomes/$(genome)/$(genome).gff.gz, /home/gnickles/FungalGenomes/$(genome)/$(genome).fna.gz, /home/gnickles/FungalGenomes/$(genome)/$(genome)_hmmer.out.gz, /home/gnickles/FungalGenomes/$(genome)/$(genome)_GIO_ncbi_parent.txt.gz
# /home/gnickles/Cassis/cassis
transfer_output_files = CASSIS_$(genome).tar.gz
#
# Tell HTCondor what amount of compute resources
request_cpus = 1
request_memory = 2GB
request_disk = 4GB
#
queue genome from /home/gnickles/FungalGenomes/QueueFile/genomes.txt
```

**GetAnchor_RunCassis.sh**

```sh
#!/bin/bash
#extracting the files
gunzip *.gz
#
genome="$1"
#setting up the gff file in the cassis compatible format
sed 's/;/\t/g' "$genome".gff | sed 's/ID=//g' | awk '{print $9,$1,$4,$5,$7}' | grep gene | awk -v OFS='\t' '{ $1=$1; print }' > temp
cat temp > "$genome".gff
rm temp
#changing the hmmer out files to list the gene name next to the protein domain

#getting all of the lines with ICS domains
if grep -q PF05141 "$genome"_hmmer.out; then
  #this means there was an ICS gene predicted by hmmer
  ###
  # Getting the gene list that has ICS domains in them
  ###
  #grabbing the proteins
  ICSProteins=$(cat "$genome"_hmmer.out | grep PF05141 | awk '{ print $1 }')
  #converting the proteins to a list of genes with the GIO file
  for line in $ICSProteins; do
    cat "$genome"_GIO_ncbi_parent.txt | grep $line | awk '{ print $7 }' | sed 's/;/\n/g' | grep Parent | sed 's/Parent=//g' >> ICSGenes.txt
  done
  ICSGenes=$(cat ICSGenes.txt)
  # looping over all the anchor genes running CASSIS each time
  counter=1
  mkdir CASSIS_"$genome"
  for anchor in $ICSGenes; do
    #make the output folder for the specific anchor gene
    mkdir CASSIS_"$genome"_"$counter"
    cassis --dir CASSIS_"$genome"_"$counter" --annotation "$genome".gff --genome "$genome".fna --meme 1.0e+005 --anchor "$anchor" --cluster "$anchor" --fimo 0.00006 --frequency 14 --gap-length 2 --prediction --num-cpus 1 --verbose
    #copying the one file I need to the folder that will be outputted
    cp CASSIS_"$genome"_"$counter"/"$anchor"/CLUSTER/*_all_predictions* CASSIS_"$genome"/
  done
  #compressing the outputFile
  tar -czvf CASSIS_"$genome".tar.gz CASSIS_"$genome"
else
  #this means there were no ICS genes predicted by hmmer, the program will exit
  echo "There were no ICS genes" > NoICS.txt
  tar -czvf CASSIS_"$genome".tar.gz NoICS.txt
  exit 1
fi
```

Requesting 4 GB of disk was sufficient for a vast majority of the genomes, however 103 runs on the supercomputer ran out of disk space. Once the rest of the genomes completed, l re-ran the one's that failed with 8 GB of disk. 

**Running the last 103**

```sh
# CHANGES IN THE SUBMIT SCRIPT
# Tell HTCondor what amount of compute resources
request_cpus = 3
request_memory = 4GB
request_disk = 8GB

##########

# TERMINAL RUNNING THE LAST 103
queue genome from /home/gnickles/FungalGenomes/QueueFile/runTheseCassis.txt
(base) [gnickles@submit-1 Cassis]$ condor_submit GetAnchor_RunCassis.sub 
Submitting job(s).......................................................................................................
103 job(s) submitted to cluster 14518870.
```

11 of the genomes ran out of disk usage so I ran one last batch with a lot of resources (overkill for a vast majority of the genomes).

```sh
# CHANGES IN THE SUBMIT SCRIPT
# Tell HTCondor what amount of compute resources
request_cpus = 4
request_memory = 5GB
request_disk = 15GB

##########

# TERMINAL RUNNING THE LAST 11
-- Schedd: submit-1.chtc.wisc.edu : <128.105.244.191:9618?... @ 12/16/21 11:40:07
OWNER    BATCH_NAME      SUBMITTED   DONE   RUN    IDLE  TOTAL JOB_IDS
gnickles ID: 14521422  12/16 11:35      _     11      _     11 14521422.0-10

11 jobs; 0 completed, 0 removed, 0 idle, 11 running, 0 held, 0 suspended
```

## Step 3: Creating gbk, fna, and gff files for each ICS BGC prediction

The following is the largest step in terms of custom code scripted for this publication. All of the raw gbk, fna, and gff files that were outputed can be found and download from: 

https://github.com/gnick18/Fungal_ICSBGCs

### Obataining the dependencies with miniconda or anaconda

```sh
#compiling the miniconda environment
#
conda create -n ICSCreator
conda activate ICSCreator
conda install pandas
conda install -c bioconda emboss
conda install -c anaconda biopython
#List of packages this installed 
# packages in environment at /home/gnickles/miniconda3/envs/ICSCreator:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main  
_openmp_mutex             4.5                       1_gnu  
biopython                 1.78             py39h7f8727e_0  
blas                      1.0                         mkl  
bottleneck                1.3.2            py39hdd57654_1  
ca-certificates           2020.10.14                    0    anaconda
certifi                   2021.10.8        py39h06a4308_0  
emboss                    6.6.0                h6debe1e_0    bioconda
expat                     2.4.1                h2531618_2  
fontconfig                2.13.1               h6c09931_0  
freetype                  2.11.0               h70c0345_0  
icu                       58.2                 he6710b0_3  
intel-openmp              2021.4.0          h06a4308_3561  
jpeg                      9d                   h7f8727e_0  
ld_impl_linux-64          2.35.1               h7274673_9  
libffi                    3.3                  he6710b0_2  
libgcc-ng                 9.3.0               h5101ec6_17  
libgd                     2.3.3                h695aa2c_0  
libgomp                   9.3.0               h5101ec6_17  
libpng                    1.6.37               hbc83047_0  
libstdcxx-ng              9.3.0               hd4cf53a_17  
libtiff                   4.2.0                h85742a9_0  
libuuid                   1.0.3                h7f8727e_2  
libwebp-base              1.2.0                h27cfd23_0  
libxml2                   2.9.12               h03d6c58_0  
lz4-c                     1.9.3                h295c915_1  
mkl                       2021.4.0           h06a4308_640  
mkl-service               2.4.0            py39h7f8727e_0  
mkl_fft                   1.3.1            py39hd3c417c_0  
mkl_random                1.2.2            py39h51133e4_0  
ncurses                   6.3                  h7f8727e_2  
numexpr                   2.7.3            py39h22e1b3c_1  
numpy                     1.21.2           py39h20f2e39_0  
numpy-base                1.21.2           py39h79a1101_0  
openssl                   1.1.1l               h7f8727e_0  
pandas                    1.3.4            py39h8c16a72_0  
pip                       21.2.4           py39h06a4308_0  
python                    3.9.7                h12debd9_1  
python-dateutil           2.8.2              pyhd3eb1b0_0  
pytz                      2021.3             pyhd3eb1b0_0  
readline                  8.1                  h27cfd23_0  
setuptools                58.0.4           py39h06a4308_0  
six                       1.16.0             pyhd3eb1b0_0  
sqlite                    3.36.0               hc218d9a_0  
tk                        8.6.11               h1ccaba5_0  
tzdata                    2021e                hda174b7_0  
wheel                     0.37.0             pyhd3eb1b0_1  
xz                        5.2.5                h7b6447c_0  
zlib                      1.2.11               h7b6447c_3  
zstd                      1.4.9                haebb681_0  
#making the environment portable
conda install -c conda-forge conda-pack
conda pack -n ICSCreator
chmod 644 ICSCreator.tar.gz
(base) [gnickles@submit-1 ICSCreator]$ ls -sh ICSCreator.tar.gz 
520M ICSCreator.tar.gz
```

### Test showing adding the 'biosynthetic flag' impacts the GCF classifications

When antismash normally runs, it adds flags in the `.gbk` files that look like this (example of normal antismash output CDS region)

```
     CDS             complement(10001..13180)
                     /Dbxref="NCBI_GP:EAY30987.1"
                     /ID="cds-EAY30987.1"
                     /Name="EAY30987.1"
                     /gbkey="CDS"
                     /gene="M23134_07394"
                     /gene_functions="biosynthetic (rule-based-clusters)
                     lanthipeptide: Lant_dehydr_N"
                     /gene_functions="biosynthetic (rule-based-clusters)
                     lanthipeptide: Lant_dehydr_C"
                     /gene_functions="biosynthetic-additional (smcogs)
                     SMCOG1155:Lantibiotic dehydratase domain protein (Score:
                     970.6; E-value: 1.3e-293)"
                     /gene_kind="biosynthetic"
                     /locus_tag="M23134_07394"
                     /phase="0"
                     /product="lantibiotic dehydratase, C terminus family"
                     /protein_id="EAY30987.1"
                     /sec_met_domain="Lant_dehydr_N (E-value: 4.2e-136,
                     bitscore: 448.7, seeds: 73, tool: rule-based-clusters)"
                     /sec_met_domain="Lant_dehydr_C (E-value: 8.3e-36, bitscore:
                     118.2, seeds: 35, tool: rule-based-clusters)"
                     /source="Genbank"
                     /transl_table=11
                     /translation="VDTVRQHLQAAYRLPALQQALYIASPDLFAQLDRWHALAALPAPN
                     KKQRKDLQNLEQKLWMFLARAAYRCTPFGTFAAVTVGEFRLQRHSSLALDKLGNMRTSS
                     RLDMVYLLRLRASLEAEPKIREQLRYFPNNTLYRNGDDMLCSFFDEQTLQHVASSKEYH
                     EVAWELFIACGQAGGMRIGEMIAHLMTGPELTEADCREFVEYLINHHLLISELYPSLNG
                     EEYLEQLHASIAQMPAAQGYLAIINHIRRQLALIDATPDPAQSLGLYQQLAGYLSEQPF
                     YQAQETAITAHPNEKTASKKELIQVDAFRDENGQALSKTVPRHLVRQIATLIDTVGEVP
                     TTKADLSTFAQAFFKRYENEAVPLLLAMDEVTGVGYPYGKRVATLPTYLQELGVVVGEG
                     GAMPQVADPKAGQTSALDALLEQKYVEALQTGAPAIHLSKGEVQSLKASQLPTKEKPRS
                     LPPAMVAMGKLLYQAPDDLEHNPKQDNPRNFRFQLYHATGPVVGNLLGRFCYGSPLLQQ
                     QLSATIAQVEPDDNSVVYAEIDHLPGQRVGNVIIRPRLRDYSIALGGANTNDVHQIPVT
                     DLWLSAQDGGKRLVLHSRKLNKQVVPRMSNAHNYAHNAHPLYKFLCDLQYQQEYRQLPP
                     FFRAQRWQQRSYLPRLVYDGTLILSPRQWRVNYKALEAVQLQKKEDMHAVSTWTHKIER
                     LRQQLRLPQAVLLVEADNKLYIDFAHPLSIGVFLQKVFDRKELVLEEVIEQHLPSQPYS
                     NEMIFPIINTAAPSSKQSPAFDAQTWVPDVQKNFSFGSEWVYFKIYLGTLGSDRLIAED
                     LHPLSKELRDSGHIDQWFFIRYTDPEAYGTHLRWRLHINHYPKEFGEIVQAFSQKLLGL
                     KQQGKIIHYVPDTYQRELNRYLPYNMALSEEWFALDSDFVCHALALLAKDEATNDTLKS
                     ALALQRIKALLQAFDLPEETRQQITSQGKTNFSREFGLDVHPQAQRLKKALAKQALPPW
                     LPPAYRLLFDQLTEQEAPLAQQIKANLQQHGASETQLRSLLSSYIHMFYNRFFASGQRL
                     EEMWGYWRLEG"
```

Adittionally I was getting this error everytime I tried running BigScape on my files oringial files:

```
  Warning: Input set has files with no Biosynthetic Genes (affects alignment mode)
   See no_biosynthetic_genes_list.txt
```

After looking into the BigScape notes I deduced what line was determing that call was:

- Version of BigScape: `BiG-SCAPE 1.0.1 (2020-01-27)`
- `/gene_kind="biosynthetic`

To test for if this label significantly alters the networks I took 100 normal AntiSmash runs and ran it on BigScape with a cuttoff of 0.5. I then made a small python script that deletes all of the biosynthetic labels.

**DeleteBiosynthetic.py** 

```python
#! /usr/bin/env python
#Author: Grant Nickles

#this script will loop over the antismash file and delete the lines that flag it has having a gene cluster
import pdb
import os

rootDirectory = r'/Volumes/T7/antismash/100gbks'

for gbk in os.listdir(rootDirectory):
    if os.path.isfile(os.path.join(rootDirectory, gbk)) and gbk.endswith(".gbk"):
        gbkPath = os.path.join(rootDirectory, gbk)
        #looping over the file
        editedFile = []
        with open(gbkPath) as gbk:
            lines = gbk.readlines()
            for index, line in enumerate(lines):
                if '/gene_kind="biosynthetic' not in line:
                    editedFile.append(line)
        #writing out the new file
        os.remove(gbkPath)
        with open(gbkPath, "w") as gbkFinal:
            for l in editedFile:
                gbkFinal.write(l)

```

When examining the outputs, **the networks that had the biosynthetic tags grouped more things together.** By putting an emphasis on the biosynthetic genes in the prediction the program puts more emphasis on these alignments. Since antismash was making these calls, I knew somewhere in their code repo must be the information to determine if something should be called biosynthetic or not. After doing a lot of digging I found the file which is used by the code to color the known bioynsthetic genes on the html file. This is perfect as it can be levereged to label the gbk files in this publication. 

[biosynthetic pfam file from AntiSmash](./BiosyntheticDomains/pfam_biosynthetic.txt)

As a result of this test, I created a final script that loops over each gbk file, and using the HMMer generated protein domain predictions, adds the `/gene_kind="biosynthetic` flag to the relavent CDS regions.

### The scripts

**ICSCreator_Supercomputer.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles
#usage: The main script, this is only used to take user inputs for key variables that the other programs need to run properly. It also calls the other python scripts in the proper order to run the pipeline.
#Argument Notes: User can enter the inputs from the terminal or comment out this section and manually type them into the script themselves. If one option is being used the other MUST be commented out or it could result in the wrong inputs being used.
from CassisToClusters import CassisToCluster
from SeqretLooper import SeqretLooper
from GBKEditor import GBKEditor
import os
import pdb
import sys
from multiprocessing import Pool, freeze_support

if __name__ == "__main__":
    try:
        import pyfiglet
        ascii_art = pyfiglet.print_figlet('ICS', colors='WHITE', font="isometric1")
        print(ascii_art)
    except:
        pass
    #the root directory that has all the genomes files and the cluster prediction files
    rootDirectory = sys.argv[1]
    windows = ['CASSIS'] #DON'T CHANGE THIS PARAMETER IF USING ON CASSIS, if you add other parameters for the window the code will likley break
    # these bottom three parameters can still be altered
    annotation = "PF05141"
    deleteIntermetiates = False
    clusterType = "ICS"
################
# CREATING THE GFF AND CHROMOSOME FASTA FILES TO BE LATER USED IN GBK CREATION
# Output Explained
#File Structure: Root/GenomeFolder_1/Annotation_Windows/filtered_chromosomes + bp_Windows/indiviual_gffs
# example
# ~ ls GCA_000002855.2_new_fixed.tar.gz_folder/Annotation_Windows/
# ~     1_297mF-Z3.fna  Window10bp/     Window15bp/     Window5bp/      Window7bp/
# ~ ls GCA_000002855.2_new_fixed.tar.gz_folder/Annotation_Windows/Window5bp/
# ~     1_GCA_000002855.2_5bp_.gff
    print("Running CassisToCluster on all the genomes")
    #listing out the subfolder
    genomeFolder = []
    for genomeFolders in os.listdir(rootDirectory): #stepping through the genome folders
        if os.path.isdir(os.path.join(rootDirectory, genomeFolders)) and str(genomeFolders).startswith("GC"):
            genomeFolder = os.path.join(rootDirectory, genomeFolders)
    # setting up the multiprocessing, it will use the max cores as defined by os.cpu_count() by default
    # with Pool() as pool:
    #     pool.starmap(CassisToCluster, [(genomeListPath, rootDirectory) for genomeListPath in genomeListPaths])
    CassisToCluster(genomeFolder, rootDirectory)
################
# USING SEQRET TO CREATE THE GBK SCAFFOLD
# GBK files will be stored instide the Annotation_Windows/Window#bp folders based on the gff files that were made by the ProteinWindow script
    print("\n#########################\n")
    print("Finished creating the annotation(gff) and chromosome(fasta) files based on the annotation provided. \nUsing seqret to create the gbk scaffolds.\n")
    with open(os.path.join(rootDirectory, "results.txt"), "a+") as result:
        result.write("\n#########################\n\n")
    SeqretLooper(rootDirectory, deleteIntermetiates)
    if deleteIntermetiates == True:
        YoN = ""
    else:
        YoN = "not"
    print("\nRunning seqret has finished. As requested the intermediate files were " + YoN + " deleted.")
################
# USING THE CUSTOM PROGRAM TO EDIT THE GBK SCAFFOLD AS BIGSCAPE READABLE
# The fixed gbk files will be the only ones in the folders once this program has finished
    print("\n#########################\n")
    print("Now starting the editing of the gbk files.")
    with open(os.path.join(rootDirectory, "results.txt"), "a+") as result:
        result.write("\n#########################\n\n")
    GBKEditor(rootDirectory, annotation)
    print("\nEditing finished of all the gbk files! All the the outputs that were on the terminal can be revisited on the results.txt file stored withing the rootDirectory provided.")
################
# Copying the results into a folder by the window size that was grabbed
    # print("\n#########################\n")
    # print("The program is copying all of the results into their own folder. Each folder contains all of the gbks for that given window size. This is found in the rootDirectory provided!")
    # GBKCopier(rootDirectory, windows)
################
# # Creating the parsed files from the copied window files.
#     print("\n#########################\n")
#     print("The program is starting the creation of the simulated Antismash parsed files")
#     antiFolder = os.path.join(rootDirectory, "anti_parsed_new")
#     os.mkdir(antiFolder)
#     pdb.set_trace()
#     AntiParsedMakerFaster(genomeFolder, clusterType, antiFolder)
################
# Program finished message
    print("\n#########################\n")
    try:
        import pyfiglet
        ascii_art = pyfiglet.print_figlet('The program is complete :)', colors='PINK')
        print(ascii_art)
    except:
        pass
```

**CassisToClusters.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles
#usage: Will create an Annotation_Windows folder while editing the fasta and gff files
#assumptions: The genomes in the root need to have the all predictions file with the other fasta gff GFF etc.
#argument notes:

import os
import sys
import pdb
import shutil # using this to remove full directories
from pandas import read_csv
from GffFilterer import GffFilterer
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord


def CassisToCluster(genomeFolder, rootDirectory):
    genome = str(genomeFolder).split("/")[-1][:15] #this removes all but the unique identifier of the genome
    # print("Starting analysis on " + genome)
    allRanges = {}
    allRanges["Cluster"] = [] #dictionary of lists, containing all the window as the key and a list of lists with each annotation [[chromosome, start, end, backboneGene][chromosome, start, end, backboneGene]]
    geneClusterBorders = {} #dictionary with anchorGene as key and list of start and stop gene as the value. This is populated from the CASSIS results
    ############
    # Getting the gene cluster boarders
    ############
    annotationFolderAlreadyMade = False
    allPredictionsFileFound = False
    for files in os.listdir(genomeFolder):
        if os.path.isfile(os.path.join(genomeFolder, files)) and "all_predictions" in str(files):
            #Reading in the file
            predictionPath = os.path.join(genomeFolder, files)
            anchorGene = files.split("_all_predictions")[0]
            predictionTSV = read_csv(predictionPath, sep='\t', header=0)
            if predictionTSV.empty == False:
                allPredictionsFileFound = True
                #making the Annotation folder, if there isn't an all_predictions file there wasn't any ICS hits
                if annotationFolderAlreadyMade == False:
                    windowsFolder = os.path.join(genomeFolder, "Annotation_Windows")
                    os.mkdir(windowsFolder)
                    annotationFolderAlreadyMade = True
                # This section grabs the genes and gets ride of the cobined genes that CASSIS sometimes has. gene-AFUA_4G01370+gene-AFUA_4G01380 becomes gene-AFUA_4G01370
                # Grabbing the first row as that is the best hit
                firstGene = str(predictionTSV.iloc[0,4]).split('+')[0]
                lastGene = str(predictionTSV.iloc[0,5]).split('+')[-1]
                geneClusterBorders[anchorGene] = [firstGene, lastGene]
    if allPredictionsFileFound == False:
        print(genome + " did not contain any hits for the annotation: ICS")
        with open(os.path.join(rootDirectory, "results.txt"), "a+") as result:
            result.write(genome + " did not contain any hits for the annotation: ICS\n")
        return #exiting out of the function early
    ############
    # Populating the allRanges Dictionary
    ############
    allClusters = []
    newFolderPaths = []
    for files in os.listdir(genomeFolder):
        if os.path.isfile(os.path.join(genomeFolder, files)) and ".gff" in str(files) and "unedited" not in str(files) and not files.startswith("."):
            gffPath = os.path.join(genomeFolder, files)
            gffFull = read_csv(gffPath, sep='\t', names=["contig","n", "annotation type", "start", "stop", "n2" ,"direction", "n3", "info"])
            gffGenesOnly = gffFull[gffFull['annotation type'].str.contains('gene')]
            GFF = gffGenesOnly.loc[:,["contig","start", "stop", "info"]]
            for anchor, boundries in geneClusterBorders.items():
                #grabbing the rows that contain the boundry genes
                firstGeneRow = GFF.loc[GFF['info'].str.contains(boundries[0] + ";",case=False)]
                lastGeneRow = GFF.loc[GFF['info'].str.contains(boundries[-1] + ";",case=False)]
                # There are some situations where two rows are pulled back due to two mRNA transcripts on the same gene being annotated
                # in this case I only take the one with the larger boundry which should match the parent gene on the gff file
                if len(firstGeneRow) > 1:
                    firstGeneRow = FindLargstRow(firstGeneRow)
                if len(lastGeneRow) > 1:
                    lastGeneRow = FindLargstRow(lastGeneRow)
                ###########
                # MAKING SURE THE CONTIGS ARE THE SAME, if they are not the program fixes it by making a fake new contig
                ###########
                if firstGeneRow['contig'].to_string(index = False) != lastGeneRow['contig'].to_string(index = False):
                    #getting the path of the fasta file
                    pathOfFasta = ""
                    for files in os.listdir(genomeFolder):
                        if os.path.isfile(os.path.join(genomeFolder, files)) and ".fna" in str(files) and "ncbi_prot" not in str(files):
                            pathOfFasta = os.path.join(genomeFolder, files)
                    newFolderPath = os.path.join(genomeFolder, "COMBINED_" + firstGeneRow['contig'].to_string(index=False).strip() + "_" + lastGeneRow['contig'].to_string(index=False).strip())
                    # sometimes multiple gene clusters fall withing the same contigs, if this happens the program won't remake the edited fasta and gff files for that folder
                    folderAlreadyMade = False
                    try:
                        os.mkdir(newFolderPath)
                        newFolderPaths.append(newFolderPath)
                    except:
                        folderAlreadyMade = True
                    editedRows = ContigsDontMatch(firstGeneRow, lastGeneRow, pathOfFasta, gffPath, newFolderPath, genome, folderAlreadyMade)
                    firstGeneRow = editedRows[0]
                    lastGeneRow = editedRows[-1]
                try:
                    start = int(firstGeneRow["start"])
                    end = int(lastGeneRow["stop"])
                    contig = firstGeneRow['contig'].to_string(index = False).strip()
                except:
                    print(genome + " needs to mannually edited!")
                    with open(os.path.join(rootDirectory, "lookIntoThese.txt"), "a+") as result:
                        result.write(genome + " needs to mannually edited!\n")
                    continue
                #now I have all the pieces to add to allClusters in the order [chromosome, start, end, backboneGene]
                allClusters.append([contig,start,end,anchor])
            # Now that the allClusters list of lists has been populated, I'll add this to the allRanges dictionary with the 'Window' key being set as 'Clusters'
    if len(allClusters) == 0:
        return # exiting out the function early
    alreadyIn = []
    for cluster in allClusters:
        if cluster[3] not in alreadyIn:
            start = cluster[1]
            stop = cluster[2]
            contig = cluster[0]
            matchingAnchors = [cluster[3]]
            alreadyIn.append(cluster[3])
            for cluster2 in allClusters:
                start2 = cluster2[1]
                stop2 = cluster2[2]
                contig2 = cluster2[0]
                if cluster2 != cluster and start == start2 and stop == stop2 and contig == contig2:
                    #they are the same cluster but multiple anchor genes
                    matchingAnchors.append(cluster2[3])
                    alreadyIn.append(cluster2[3])
            allRanges['Cluster'].append([contig, start, stop, matchingAnchors])
    #Checking for multiple backbone genes in same clusters
    ############
    # Filtering the GFF and Fasta Files
    ############
    GffFilterer(allRanges, genomeFolder, genome)
    totalHits = None
    for element in allRanges:
        totalHits = len(allRanges[element])
    print("Filtering complete for " + genome + ", total annotation hits: " + str(totalHits))
    with open(os.path.join(rootDirectory, "results.txt"), "a+") as result:
        result.write("Filtering complete for " + genome + ", total annotation hits: " + str(totalHits) + "\n")
    # deleting the intermediate combined contig folder
    for folder in newFolderPaths:
        shutil.rmtree(folder)

def ContigsDontMatch(firstGeneRow, lastGeneRow, pathOfFasta, pathOfGFF, newFolderPath, genome, folderAlreadyMade):
    editedContigSeq = ""
    metaData = ""
    firstGeneContig = firstGeneRow['contig'].to_string(index=False).strip()
    orderedContigLengths = []
    firstContigLength = 0
    contigsCombined = 1
    lastGeneContig = lastGeneRow['contig'].to_string(index=False).strip()
    #########
    # Making the edited fasta file
    #########
    foundFirstContig = False
    for record in SeqIO.parse(pathOfFasta, "fasta"):
        if firstGeneContig in record.id:
            editedContigSeq = editedContigSeq + record.seq
            foundFirstContig = True
            orderedContigLengths.append([record.id,len(record.seq)])
        elif foundFirstContig == True and lastGeneContig not in record.id:
            orderedContigLengths.append([record.id,len(record.seq)])
            contigsCombined += 1
            editedContigSeq = editedContigSeq + record.seq
        elif lastGeneContig in record.id and foundFirstContig == True:
            orderedContigLengths.append([record.id,len(record.seq)])
            editedContigSeq = editedContigSeq + record.seq
            contigsCombined += 1
            break
    #########
    # Figuring out how much to add to each start and stop location
    #########
    counter = 1
    sequentialLength = 0
    amountToAddToStartStop = {}
    for contig in orderedContigLengths:
        #it's the first contig combined
        if counter == 1:
            amountToAddToStartStop[contig[0]] = 0
            sequentialLength = contig[1]
            counter += 1
        #it's the last contig combined
        elif counter == contigsCombined:
            amountToAddToStartStop[contig[0]] = sequentialLength
        #it's one of the middle contigs that was added
        else:
            amountToAddToStartStop[contig[0]] = sequentialLength
            sequentialLength = sequentialLength + contig[1]
            counter += 1
    #########
    # Looping through the fasta file to get the descriptive data (not the easiest way to do this but it works so I didn't bother editing)
    #########
    if folderAlreadyMade == False:
        with open(pathOfFasta) as fna:
            lines = fna.readlines()
            for index, line in enumerate(lines):
                strippedLine = line.rstrip()
                # grabbing all of the meta information after the contig
                if index == 0:
                    metaData = " ".join(line.split()[1:])
                    break
        editedContigSeq = str(editedContigSeq) #converting the seq object to a string
        newRecord = SeqRecord(Seq(editedContigSeq), id=firstGeneContig + "_" + lastGeneContig, description = metaData)
        #Writing the edited fasta file
        newFastaPath = os.path.join(newFolderPath, "EDITED_" + genome + ".fna")
        SeqIO.write(newRecord, newFastaPath, "fasta")
    #########
    # Making the edited gff file
    #########
    gffFull = read_csv(pathOfGFF, sep='\t', names=["contig","n", "annotation type", "start", "stop", "n2" ,"direction", "n3", "info"])
    for key, amountToAdd in amountToAddToStartStop.items():
        #changing all of the start and stop columns
        gffFull.loc[gffFull['contig'] == key, 'start'] = gffFull['start'] + amountToAdd
        gffFull.loc[gffFull['contig'] == key, 'stop'] = gffFull['stop'] + amountToAdd
        #changing the contig row
        gffFull.loc[gffFull['contig'] == key, 'contig'] = firstGeneContig + "_" + lastGeneContig
    if folderAlreadyMade == False:
        # writing out the new gff file
        newGFFPath = os.path.join(newFolderPath, "EDITED_" + genome + ".gff")
        gffFull.to_csv(newGFFPath, index = False, sep='\t')
    #########
    # Editing and returning the first and last gene row to match the new boundries
    #########
    # pandas is giving me a warning message here but the code works fine, I'm supressing printing out the warning with this block
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        firstGeneRow.contig = firstGeneContig + "_" + lastGeneContig
        lastGeneRow.contig = firstGeneContig + "_" + lastGeneContig
        lastGeneRow.start = lastGeneRow['start'].astype(int).copy() + sequentialLength
        lastGeneRow.stop = lastGeneRow['stop'].astype(int).copy() + sequentialLength
        return [firstGeneRow, lastGeneRow]


def FindLargstRow(rows):
    rows['length'] = rows['stop'] - rows['start']
    largestRowIndex = rows['length'].idxmax() #finding the row with the largest length
    largestRow = rows.loc[[largestRowIndex]]
    return largestRow
```

**GffFilterer.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles
#usage: GffFilterer [dictionary from ProteinWindows] [genome directory] [name of the genome]
#This code is used inside GffFilterer and will store the filtered chromosomes inside the Annotation_Windows folder
#Argument Notes: allRanges = dictionary of lists, containing all the window as the key and a list of lists with each annotation [[chromosome, start, end][chromosome, start, end]]
#This code is used inside ProteinWindows.py and will store the filtered gff files inside the corrosponding windows folder within the Annotation_Windows folder
import os
import sys
import pdb
from FastaFilterer import FastaFilterer
from FastaFilterer import CombinedContigsFastaFilterer

# the subDirectory is the path to the genome folder, and genome is the name of said genome
def GffFilterer(allRanges, subDirectory, genome):
    codeDirectory = os.getcwd()
    #Saving the path of the annotation window
    windowsFolder = None
    for annoDir in os.listdir(subDirectory):
        if os.path.isdir(os.path.join(subDirectory, annoDir)) and "Annotation" in str(annoDir):
            windowsFolder = os.path.join(subDirectory, annoDir)
            break
    #######
    # CHECKING FOR COMBINED CONTIG FOLDERS
    #######
    combinedContigs = False
    for folders in os.listdir(subDirectory):
        if os.path.isdir(os.path.join(subDirectory, folders)) and "COMBINED" in str(folders):
            combinedContigs = True
    # if there is combined gffs the program will store them as a list of lists
    # List of list structure: [[basename of gff, path][basename of gff, path]]
    # all edited gffs will have the work EDITED in the front of them
    if combinedContigs == True:
        allGFFPaths = CheckForEditedGffs(subDirectory)
    #finding and storing the path of the gffFile
    else:
        gffPath = None
        for file in os.listdir(subDirectory): #this part of the code assumes only one gff is present in each subFolder
            if os.path.isfile(os.path.join(subDirectory, file)):
                if str(file).endswith('.gff') and "unedited" not in str(file): #only opens the gff file
                    gffPath = os.path.join(subDirectory, file)
                    break
    for window in allRanges:
        windowSubFolder = os.path.join(windowsFolder, "CASSIS")
        os.mkdir(windowSubFolder)
        #creating the framework of the dictionary which will store all the indexes to save for each BGC annotation hit
        indexSave = {} #dictionary: Key represents every BGC hit, the value is an index of every line to save for that range of bp and if it is a backbone gene [index, True/False]
        bgcStartStop = {} #dictionary: Key is same as indexSave, the value is a list with the start bp of the first gene and the last bp of the last gene
        for annotation in allRanges[window]:
            chromosome = annotation[0]
            start = annotation[1]
            end = annotation[2]
            string = str(chromosome) + "_" + str(start) + "_" + str(end)
            indexSave[string] = []
            bgcStartStop[string] = [-1,-1]
        #########
        # This version of the gff opening will only happen if there are multiple gffs as it is less effienct (BigO time is exponential)
        #########
        if combinedContigs == True:
            results = CombinedContigsIndexSave(allGFFPaths, allRanges, indexSave, bgcStartStop, window)
            indexSave = results[0]
            bgcStartStop = results[1]
            RunFastaGFFEditing(indexSave, bgcStartStop, genome, window, windowSubFolder, allGFFPaths, subDirectory, codeDirectory)
        ########
        # This normal version will use the same gff for all of the annotations and is computationally more efficient (BigO time is linear)
        ########
        #saving all of the indexes of genes falling in each annotation range into indexSave
        else:
            with open(gffPath) as gff:
                lines = gff.readlines()
                for index, line in enumerate(lines):
                    line = line.strip('\n')
                    #Searching for the lines the contain only the locus tags within the largest files
                    splitLine = line.split('\t')
                    lineChromosome = splitLine[0]
                    lineStart = int(splitLine[3])
                    lineEnd = int(splitLine[4])
                    annotationType = splitLine[2]
                    for annotation in allRanges[window]:
                        chromosome = annotation[0]
                        start = annotation[1]
                        end = annotation[2]
                        backboneNames = annotation[3]
                        firstFound = False
                        lastFound = False
                        if lineChromosome == chromosome and annotationType.lower() != 'region': #if the chromosomes match and it's not a region annotation...
                            s = str(chromosome) + "_" + str(start) + "_" + str(end)
                            if start <= lineStart <= end: #if the start of the feature is in the range..
                                indexSave[s].append([index,IsBackboneLine(backboneNames,line)]) #save the line
                                if lineStart < bgcStartStop[s][0] or bgcStartStop[s][0] < 0: #if this annotation is the first or lower than other found annotations
                                    bgcStartStop[s][0] = lineStart
                                elif lineEnd > bgcStartStop[s][1] or bgcStartStop[s][1] < 0: #if this annotation is the first or higher than other found annotations
                                    bgcStartStop[s][1] = lineEnd
                            elif start <= lineEnd <= end: #or if the end of the feature is in the range..
                                indexSave[s].append([index,IsBackboneLine(backboneNames,line)]) #save the line
                                if lineStart < bgcStartStop[s][0] or bgcStartStop[s][0] < 0: #if this annotation is the first or lower than other found annotations
                                    bgcStartStop[s][0] = lineStart
                                elif lineEnd > bgcStartStop[s][1] or bgcStartStop[s][1] < 0: #if this annotation is the first or higher than other found annotations
                                    bgcStartStop[s][1] = lineEnd
            #indexSave will now have all of the indexes for each gff file. This next section of the code will write out all of the unique gff files
            #Counter is a simple counting varaiable that will count up after each subsequent annoation in indexSave
            #Naming for each gff file will be: genome + "_" + window + "bp" + "_" + counter + ".gff"
            #Example: GCA_0000028552_5bp_1.gff, GCA_0000028552_5bp_2.gff ...
            counter = 1
            for BGC in indexSave:
                start = bgcStartStop[BGC][0]
                end = bgcStartStop[BGC][1]
                #creates the filtered fasta files to be used by seqret later on
                try:
                    FastaFilterer(BGC, subDirectory, counter, start, end, str(window))
                except:
                    print("There was an issue with using FastaFilterer on " + genome)
                #opening and creating the gff file
                fileName = "region"+str(counter)+"_"+genome+"_"+str(window)+".gff"
                Temp1 = os.path.join(windowSubFolder, "Temp1_" + fileName)
                with open(Temp1, "a+") as gffFile:
                    filteredLines = []
                    for i in indexSave[BGC]:
                        if i[1] == False:
                            filteredLines.append(lines[i[0]])
                        else: #if the line is a backboneGene the code will add this info to the key
                            editedLine = lines[i[0]][:-2] + ";backbone=True\n"
                            filteredLines.append(editedLine)
                    for l in filteredLines:
                        gffFile.write(l)
                counter += 1
                #Editing the start and stop columns
                os.chdir(windowSubFolder)
                subtractAmount = str(start - 1)
                Temp2 = os.path.join(windowSubFolder, "Temp2_" + fileName)
                cmdStartEdit = """awk 'BEGIN { FS=OFS="\t" } {$4 = $4 - """ + subtractAmount + "; print}' " + "Temp1_" + fileName + " > Temp2_" + fileName
                cmdEndEdit = """awk 'BEGIN { FS=OFS="\t" } {$5 = $5 - """ + subtractAmount + "; print}' " + "Temp2_" + fileName + " > " + fileName
                os.system(cmdStartEdit)
                os.system(cmdEndEdit)
                os.remove(Temp1)
                os.remove(Temp2)
                os.chdir(codeDirectory)

def RunFastaGFFEditing(indexSave, bgcStartStop, genome, window, windowSubFolder, allGFFPaths, subDirectory, codeDirectory):
    #indexSave will now have all of the indexes for each gff file. This next section of the code will write out all of the unique gff files
    #Counter is a simple counting varaiable that will count up after each subsequent annoation in indexSave
    #Naming for each gff file will be: genome + "_" + window + "bp" + "_" + counter + ".gff"
    #Example: GCA_0000028552_5bp_1.gff, GCA_0000028552_5bp_2.gff ...
    counter = 1
    for BGC in indexSave:
        start = bgcStartStop[BGC][0]
        end = bgcStartStop[BGC][1]
        contigName = "_".join(BGC.split("_")[0:2])
        ## FINDING THE CORRECT FASTA FILE
        allFAs = GetTheEditedFastas(subDirectory)
        fastaPath = ""
        uneditedPath = ""
        normalFasta = False
        for fastaFiles in allFAs:
            # if the chromosome is in the name of the file it is an edited file, this will never be the case for gene clusters that should need the full unedited fasta
            if contigName in fastaFiles:
                fastaPath = fastaFiles
            elif "EDITED_" not in fastaFiles:
                uneditedPath = fastaFiles
        if len(fastaPath) == 0:
            fastaPath = uneditedPath
            normalFasta = True
        #creates the filtered fasta files to be used by seqret later on
        try:
            CombinedContigsFastaFilterer(BGC, subDirectory, fastaPath, counter, start, end, str(window), normalFasta)
        except:
            print("There was an issue with using CombinedContigsFastaFilterer on " + genome)
        #opening and creating the gff file
        fileName = "region"+str(counter)+"_"+genome+"_"+str(window)+".gff"
        Temp1 = os.path.join(windowSubFolder, "Temp1_" + fileName)
        #########
        # POPULATING THE LINES VARIABLE WITH CORRECT GFF FILE
        #########
        gffPath = ""
        lines = ""
        uneditedPath = ""
        for gffFiles in allGFFPaths:
            # if the chromosome is in the name of the file it is an edited file, this will never be the case for gene clusters that should need the full unedited gff
            if contigName in gffFiles:
                gffPath = gffFiles
            elif "EDITED_" not in gffFiles:
                uneditedPath = gffFiles
        if len(gffPath) == 0:
            gffPath = uneditedPath
        with open(gffPath) as gff:
            lines = gff.readlines()
        #########
        # WRITING OUT THE GFF FILE
        #########
        with open(Temp1, "a+") as gffFile:
            filteredLines = []
            for i in indexSave[BGC]:
                if i[1] == False:
                    filteredLines.append(lines[i[0]])
                else: #if the line is a backboneGene the code will add this info to the key
                    editedLine = lines[i[0]][:-2] + ";backbone=True\n"
                    filteredLines.append(editedLine)
            for l in filteredLines:
                gffFile.write(l)
        counter += 1
        #Editing the start and stop columns
        os.chdir(windowSubFolder)
        subtractAmount = str(start - 1)
        Temp2 = os.path.join(windowSubFolder, "Temp2_" + fileName)
        cmdStartEdit = """awk 'BEGIN { FS=OFS="\t" } {$4 = $4 - """ + subtractAmount + "; print}' " + "Temp1_" + fileName + " > Temp2_" + fileName
        cmdEndEdit = """awk 'BEGIN { FS=OFS="\t" } {$5 = $5 - """ + subtractAmount + "; print}' " + "Temp2_" + fileName + " > " + fileName
        os.system(cmdStartEdit)
        os.system(cmdEndEdit)
        os.remove(Temp1)
        os.remove(Temp2)
        os.chdir(codeDirectory)


def CombinedContigsIndexSave(allGFFPaths, allRanges, indexSave, bgcStartStop, window):
    for annotation in allRanges[window]:
        chromosome = annotation[0]
        start = annotation[1]
        end = annotation[2]
        backboneNames = annotation[3]
        firstFound = False
        lastFound = False
        gffPath = ""
        #######
        # Finding the correct gff file to open
        #######
        uneditedPath = ""
        for gffFiles in allGFFPaths:
            # if the chromosome is in the name of the file it is an edited file, this will never be the case for gene clusters that should need the full unedited gff
            if chromosome in gffFiles:
                gffPath = gffFiles
            elif "EDITED_" not in gffFiles:
                uneditedPath = gffFiles
        if len(gffPath) == 0:
            gffPath = uneditedPath
        #saving all of the indexes of genes falling in each annotation range into indexSave
        with open(gffPath) as gff:
            lines = gff.readlines()
            for index, line in enumerate(lines):
                if index != 0:
                    line = line.strip('\n')
                    #Searching for the lines the contain only the locus tags within the largest files
                    splitLine = line.split('\t')
                    lineChromosome = splitLine[0]
                    lineStart = int(splitLine[3])
                    lineEnd = int(splitLine[4])
                    annotationType = splitLine[2]
                    if lineChromosome == chromosome and annotationType.lower() != 'region': #if the chromosomes match and its not a region annotation...
                        s = str(chromosome) + "_" + str(start) + "_" + str(end)
                        if start <= lineStart <= end: #if the start of the feature is in the range..
                            indexSave[s].append([index,IsBackboneLine(backboneNames,line)]) #save the line
                            if lineStart < bgcStartStop[s][0] or bgcStartStop[s][0] < 0: #if this annotation is the first or lower than other found annotations
                                bgcStartStop[s][0] = lineStart
                            elif lineEnd > bgcStartStop[s][1] or bgcStartStop[s][1] < 0: #if this annotation is the first or higher than other found annotations
                                bgcStartStop[s][1] = lineEnd
                        elif start <= lineEnd <= end: #or if the end of the feature is in the range..
                            indexSave[s].append([index,IsBackboneLine(backboneNames,line)]) #save the line
                            if lineStart < bgcStartStop[s][0] or bgcStartStop[s][0] < 0: #if this annotation is the first or lower than other found annotations
                                bgcStartStop[s][0] = lineStart
                            elif lineEnd > bgcStartStop[s][1] or bgcStartStop[s][1] < 0: #if this annotation is the first or higher than other found annotations
                                bgcStartStop[s][1] = lineEnd
    return [indexSave, bgcStartStop]


def GetTheEditedFastas(subDirectory):
    allFAs = []
    # first storing the normal fasta file
    for files in os.listdir(subDirectory):
        if os.path.isfile(os.path.join(subDirectory, files)) and str(files).endswith('.fna') and "ncbi_prot" not in str(files):
            allFAs.append(os.path.join(subDirectory, files))
            break #there should only be one occurance so break after finding for efficiency
    # now going through and storing all of the other fasta files
    for folders in os.listdir(subDirectory):
        if os.path.isdir(os.path.join(subDirectory, folders)) and "COMBINED" in str(folders):
            combinedFolder = os.path.join(subDirectory, folders)
            for editedFiles in os.listdir(combinedFolder):
                if os.path.isfile(os.path.join(combinedFolder, editedFiles)) and str(editedFiles).endswith('.fna'):
                    allFAs.append(os.path.join(combinedFolder, editedFiles))
    return allFAs

def CheckForEditedGffs(subDirectory):
    allGffs = []
    # first storing the normal gff file
    for files in os.listdir(subDirectory):
        if os.path.isfile(os.path.join(subDirectory, files)) and str(files).endswith('.gff') and "unedited" not in str(files):
            allGffs.append(os.path.join(subDirectory, files))
            break #there should only be one occurance so break after finding for efficiency
    # now going through and storing all of the other gff files
    for folders in os.listdir(subDirectory):
        if os.path.isdir(os.path.join(subDirectory, folders)) and "COMBINED" in str(folders):
            combinedFolder = os.path.join(subDirectory, folders)
            for editedFiles in os.listdir(combinedFolder):
                if os.path.isfile(os.path.join(combinedFolder, editedFiles)) and str(editedFiles).endswith('.gff'):
                    allGffs.append(os.path.join(combinedFolder, editedFiles))
    return allGffs

def IsBackboneLine(backboneNames, line):
    key = line.split("\t")[8]
    if "ID=gene" in key:
        name = "-".join(key.split(";")[0].split("=")[-1].split("-")[1:])
        for b in backboneNames:
            if name in b:
                return True
    return False
#key.split(";")[1].startswith("Name") or
```

**GBKEditor.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles
#usage: GBKEditor [rootDirectory which is taken in my GBKCreator] [Annotation that was used to make the gbk files]
#This is used to edit the GBK files that were made using Emboss's Seqret tool. The editing is so the programs can be read into BigScape and other network programs.
import os
import sys
import pdb

def removeString(stringToRemove, line):
    editedLine = "".join(line.split(stringToRemove))
    return editedLine

def GBKEditor(rootDirectory, annotation):
    for genomeFolders in os.listdir(rootDirectory): #stepping through the genome folders
        if os.path.isdir(os.path.join(rootDirectory, genomeFolders)) and genomeFolders.startswith("GC"):
            genomeFolder = os.path.join(rootDirectory, genomeFolders)
            genome = str(genomeFolders)[:15] #this removes all but the unique identifier of the genome
            print("Editing the gbk files of " + genome)
            with open(os.path.join(rootDirectory, "results.txt"), "a+") as result:
                result.write("Editing the gbk files of " + genome + "\n")
            for annotationFolders in os.listdir(genomeFolder): #stepping into the Annotation_Windows folder inside the genome folder
                if os.path.isdir(os.path.join(genomeFolder, annotationFolders)):
                    annoFolder = os.path.join(genomeFolder, annotationFolders)
                    #this process will occur for every gbk file in every window folder. First the program will loop through all of the window folders.
                    for windowFolders in os.listdir(annoFolder):
                        if os.path.isdir(os.path.join(annoFolder, windowFolders)):
                            windowFolder = os.path.join(annoFolder, windowFolders)
                            windowSize = windowFolders[6:]
                            #this process will loop through the window folder and repeat for every gbk file inside the folder
                            for file in os.listdir(windowFolder):
                                #will only open the gbk files in the folder with this line
                                if os.path.isfile(os.path.join(windowFolder, file)) and str(file).endswith('.gbk'):
                                    gbkPath = os.path.join(windowFolder, file)
                                    gbkEdited = []
                                    totalLength = 0
                                    with open(gbkPath) as gbk:
                                        #all the the gbk editing
                                        lines = gbk.readlines()
                                        insideOrigin = False
                                        insideFeatures = False
                                        #Some sequences have weird formatting where the seq locations include the name of the sequence:
                                        #   gene            complement(NW_024065878.1:1..882)
                                        # and not
                                        #   gene            complement(1..2256)
                                        # if the program finds an occurance of this it sets locatedRemovalString to be true and will remove all occurences of said string
                                        # this is ok to remove becuase it is the sequence name and is never anywhere in the annotation descriptions
                                        stringToRemove = ""
                                        locatedRemovalString = False
                                        for index, line in enumerate(lines):
                                            #editing the header in the script
                                            begingingOfLine = line[0:8]
                                            lettersInBegining = any(c.isalpha() for c in begingingOfLine)
                                            if index == 0:
                                                totalLength = line.split()[2]
                                                gbkEdited.append(line)
                                                locus = line.split()[1] #grabbing the locus in the script
                                                #adding the def line
                                                definition = "DEFINITION  " + genome + " results based on the annotation " + annotation + " and " + windowSize + " inputs.\n"
                                                gbkEdited.append(definition)
                                                #adding the accession
                                                accession = "ACCESSION   " + genome + "\n"
                                                gbkEdited.append(accession)
                                                #adding the VERSION of the software
                                                version = "VERSION     1\n"
                                                gbkEdited.append(version)
                                                #adding the kewwords line
                                                kewwords = "KEYWORDS    .\n"
                                                gbkEdited.append(kewwords)
                                                #adding the source line
                                                source = "SOURCE      ICSCreator.py written by Grant Nickles (gnickles@wisc.edu)\n"
                                                gbkEdited.append(source)
                                                #adding the organism line
                                                orgamism = "ORGANISM    " + genome + "\n"
                                                gbkEdited.append(orgamism)
                                            #checking if it's inside the features:
                                            elif line.startswith("FEATURES"):
                                                insideFeatures = True
                                                gbkEdited.append(line)
                                            # if the line shows the location of the gene, mRNA etc. and it's not in the origin and the program hasn't already located the removal string if there is any
                                            elif lettersInBegining and not insideOrigin and not locatedRemovalString and insideFeatures:
                                                #checking if it needs to be edited
                                                if ":" in line:
                                                    locatedRemovalString = True
                                                    #checking if the line is a complement line or not
                                                    stringToRemove = ""
                                                    if "complement" in line:
                                                        stringToRemove = str(line.split("(")[-1].split(":")[0]) + ":"
                                                    else:
                                                        stringToRemove = str(line.split("     gene            ")[-1].split(":")[0]) + ":"
                                                    gbkEdited.append(removeString(stringToRemove, line))
                                                else:
                                                    gbkEdited.append(line)
                                            #making sure the program isn't insid the origin
                                            elif line.rstrip() == 'ORIGIN' and insideFeatures:
                                                insideOrigin = True
                                                gbkEdited.append(line)
                                            elif locatedRemovalString == True:
                                                gbkEdited.append(removeString(stringToRemove, line))
                                            #for all other lines
                                            else:
                                                #saving all other lines as is
                                                gbkEdited.append(line)
                                    #overwritting the gbk file with the edited lines
                                    with open(gbkPath, "w") as gbkFinal:
                                        for l in gbkEdited:
                                            gbkFinal.write(l)

# #testing
# rootDirectory = r'/Volumes/HardDrive/test'
# annotation = "PF05141"
# GBKEditor(rootDirectory, annotation)
```

**FastaFilterer.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles
#usage: FastaFilterer [biosynthetic gene cluster as a list: [chromosome, start location, stop location]] [genome directory] [annotation number]
#This code is used inside GffFilterer and will store the filtered chromosomes inside the Annotation_Windows folder
import os
import sys
import pdb
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SeqIO

def FastaFilterer(BGC, subDirectory, count, start, end, window):
    chromosome = BGC.split("_")[0]
    fnaPath = None
    for file in os.listdir(subDirectory): #this part of the code assumes only one gff is present in each subFolder
        if os.path.isfile(os.path.join(subDirectory, file)):
            if str(file).endswith('.fna'): #only opens the fna scaffold file, this part of the code must be changed if the file ending is not fna and is fa or fasta
                fnaPath = os.path.join(subDirectory, file)
                break
    recordToSave = None
    for record in SeqIO.parse(fnaPath, "fasta"):
        if chromosome in record.description:
            recordToSave = record
            break
    annotatedPath = os.path.join(subDirectory, "Annotation_Windows")
    fileName = os.path.join(annotatedPath, str(count) + "_" + window + "_" + chromosome + ".fna")
    SeqIO.write(recordToSave[(start-1):end], fileName, "fasta")

#will only run this if there are combined contigs
def CombinedContigsFastaFilterer(BGC, subDirectory, fnaPath, count, start, end, window, normalFasta):
    #TODO: should I change this join statement to a - instead of _? might screw up the seqret part
    if normalFasta == True:
        chromosome = BGC.split("_")[0]
    else:
        chromosome = "_".join(BGC.split("_")[0:2])
    recordToSave = None
    # there should be only one record to save for this section but I'll recycle the code from above as it will still work
    recordToSave = None
    for record in SeqIO.parse(fnaPath, "fasta"):
        if chromosome in record.description:
            recordToSave = record
            break
    annotatedPath = os.path.join(subDirectory, "Annotation_Windows")
    fileName = os.path.join(annotatedPath, str(count) + "_" + window + "_" + chromosome + ".fna")
    SeqIO.write(recordToSave[(start-1):end], fileName, "fasta")


# #test script
# subDirectory = "/Volumes/HardDrive/Test/GCA_000002855.2.tar.gz_folder/Annotation_Windows/"
# count = "1"
# BGC = "297mF-Z1_00010101_00002492"
# FastaFilterer(BGC, subDirectory, count)
```

**RunICSCreator.sub**

```sh
# Specify the HTCondor Universe (vanilla is the default)
universe = vanilla
#
log = $(genome)_ICS.log
error = $(genome)_ICS.err
#
# Specify your executable (single binary or a script that runs several
executable = RunICSCreator.sh
arguments = $(genome)
output = $(genome)_ICS.out
#
#make sure it has gluster
requirements = (Target.HasCHTCStaging == true) && (OpSysMajorVer == 7)
#
# Specify that HTCondor should transfer files to and from the
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files =/home/gnickles/FungalGenomes/$(genome)/$(genome).gff.gz, /home/gnickles/FungalGenomes/$(genome)/$(genome).fna.gz, /home/gnickles/FungalGenomes/$(genome)/$(genome)_hmmer.out.gz, /home/gnickles/FungalGenomes/$(genome)/$(genome)_GIO_ncbi_parent.txt.gz, /home/gnickles/Cassis/CassisResults/CASSIS_$(genome).tar.gz, /home/gnickles/ICSCreator/CassisToClusters.py, /home/gnickles/ICSCreator/FastaFilterer.py, /home/gnickles/ICSCreator/GBKEditor.py, /home/gnickles/ICSCreator/GffFilterer.py, /home/gnickles/ICSCreator/ICSCreator_Supercomputer.py, /home/gnickles/ICSCreator/SeqretLooper.py, /home/gnickles/ICSCreator/ICSCreator.tar.gz, /home/gnickles/ICSCreator/RunICSCreator.sh
transfer_output_files = $(genome)_ICSPredictions.tar.gz
#
# Tell HTCondor what amount of compute resources
request_cpus = 1
request_memory = 5GB
request_disk = 5GB
#
queue genome from /home/gnickles/FungalGenomes/QueueFile/LastOne.txt
```

**RunICSCreator.sh**

```sh
#!/bin/bash
#####
# Unpacking the miniconda environment
#####
# # have job exit if any command returns with non-zero exit status (aka failure)
# set -e

# replace env-name on the right hand side of this line with the name of your conda environment
ENVNAME=ICSCreator
# if you need the environment directory to be named something other than the environment name, change this line
ENVDIR=$ENVNAME

# these lines handle setting up the environment; you shouldn't have to modify them
export PATH
mkdir $ENVDIR
tar -xzf $ENVNAME.tar.gz -C $ENVDIR
. $ENVDIR/bin/activate

#####
# Setting up the Folder
#####
genome="$1"
mkdir "$genome"
mv *"$genome".* *"$genome"_* ./"$genome"
cd "$genome"

#####
# Unzipping everything
#####
tar -xzvf CASSIS_"$genome".tar.gz
rm CASSIS_"$genome".tar.gz
gunzip *.gz
#moving the Cassis predictions into the base folder
mv ./CASSIS_"$genome"/* ./
rm -r ./CASSIS_"$genome"/
cd ..
#####
# Running the program
#####
rootDir=$(pwd)
python3 ICSCreator_Supercomputer.py $rootDir

#####
# Moving around the results if there was a BGC prediction
#####
mkdir "$genome"_ICSPredictions
if [ -d ./"$genome"/Annotation_Windows ]; then
  cp ./"$genome"/Annotation_Windows/*.fna ./"$genome"_ICSPredictions
  cp ./"$genome"/Annotation_Windows/CASSIS/*.gbk ./"$genome"_ICSPredictions
  cp ./"$genome"/Annotation_Windows/CASSIS/*.gff ./"$genome"_ICSPredictions
fi

#####
# Adding in the info on the number of annotation hits for this genome
#####
cat results.txt | grep Filtering > ./GCF_010015735.1_ICSPredictions/AnnotationHits.txt

#####
# Compressing the output
#####
tar -czvf "$genome"_ICSPredictions.tar.gz "$genome"_ICSPredictions
```

**SeqretLooper.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles
#usage: SeqretLooper [rootDirectory which is taken in my GBKCreator] [Boolean: Delete the intermediate files?]
#This code is used after the ProteinWindows program has finished. It will loop through every genome folder to create gbk files for each gff file. Uses the filtered fna files to complete this.
import os
import sys
import pdb

def SeqretLooper(rootDirectory, deleteIntermetiates):
    codeDirectory = os.getcwd()
    for genomeFolders in os.listdir(rootDirectory): #stepping through the genome folders
        if os.path.isdir(os.path.join(rootDirectory, genomeFolders)) and genomeFolders.startswith("GC"):
            genome = str(genomeFolders)[:15] #this removes all but the unique identifier of the genome
            print(genome + " running on seqret")
            with open(os.path.join(rootDirectory, "results.txt"), "a+") as result:
                result.write(genome + " running on seqret" + "\n")
            genomeFolder = os.path.join(rootDirectory, genomeFolders)
            for annotationFolders in os.listdir(genomeFolder): #stepping into the Annotation_Windows folder inside the genome folder
                if os.path.isdir(os.path.join(genomeFolder, annotationFolders)):
                    annoFolder = os.path.join(genomeFolder, annotationFolders)
                    chromosomes = {} #dictionary storing the paths of the fna files
                    for file in os.listdir(annoFolder): #finding and storing the paths of the chromosome files, key is the number of the annotation and the value is the full file path
                        if os.path.isfile(os.path.join(annoFolder, file)) and file.endswith(".fna"):
                            fastaPath = os.path.join(annoFolder, file)
                            number = file.split("_")[0]
                            windowSize = file.split("_")[1]
                            numberAndWindow = number + "_" + windowSize
                            chromosomes[numberAndWindow] =  fastaPath
                    for windowFolders in os.listdir(annoFolder): #stepping into each Windows subfolder
                        if os.path.isdir(os.path.join(annoFolder, windowFolders)):
                            windowFolder = os.path.join(annoFolder, windowFolders)
                            windowSize = windowFolders[6:]
                            for gff in os.listdir(windowFolder): #listing all of the gff files
                                name = gff[:-4]
                                number = (name.split("_")[0].split("region")[-1])
                                gffPath = os.path.join(windowFolder, gff)
                                #finding what chromosome scaffold to use
                                try:
                                    chromosomePath = chromosomes[number + "_" + "Cluster"]
                                except:
                                    print(genome + " needs to mannually edited!")
                                    with open(os.path.join(rootDirectory, "lookIntoThese.txt"), "a+") as result:
                                        result.write(genome + " needs to mannually edited!/n")
                                    continue
                                os.chdir(windowFolder)
                                #running seqret and pipping the std.error and std.out into the results.txt file for future reference
                                cmd = "seqret -sequence " + chromosomePath + " -feature --fformat gff -fopenfile " + gffPath + " -osformat genbank -auto -osname2 " + name + " -osextension2 gbk >> ../../../results.txt 2>&1"
                                os.system(cmd)
                                os.chdir(codeDirectory)
                                if deleteIntermetiates == True:
                                    #delete all of the intermediate gff files
                                    os.remove(gffPath)
                    if deleteIntermetiates == True:
                        #delete the intermediate chromosome fasta file
                        for fnaPath in chromosomes.values():
                            os.remove(fnaPath)

# #testing
# rootDirectory = r'/Volumes/HardDrive/TestingCassisOutput/'
# deleteIntermetiates = False
# SeqretLooper(rootDirectory, deleteIntermetiates)
```

**AddBiosyntheticFlags.sub**

```
# Specify the HTCondor Universe (vanilla is the default)
universe = vanilla
#
log = $(genome)_biosyn.log
error = $(genome)_biosyn.err
#
# Specify your executable (single binary or a script that runs several
executable = AddBiosyntheticFlags.sh
arguments = $(genome)
output = $(genome)_biosyn.out
#
#make sure it has gluster
requirements = (Target.HasCHTCStaging == true) && (OpSysMajorVer == 7)
#
# Specify that HTCondor should transfer files to and from the
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = /home/gnickles/FungalGenomes/$(genome)/$(genome)_hmmer.out.gz, /home/gnickles/BiosyntheticGenomeTag/gbks/$(genome).tar.gz, /home/gnickles/BiosyntheticGenomeTag/pfam_biosynthetic.txt, /home/gnickles/ICSCreator/ICSCreator.tar.gz, /home/gnickles/BiosyntheticGenomeTag/AddBiosyntheticFlags.py
transfer_output_files = $(genome)_biosynGBKs.tar.gz
#
# Tell HTCondor what amount of compute resources
request_cpus = 1
request_memory = 1GB
request_disk = 3GB
#
queue genome from /home/gnickles/BiosyntheticGenomeTag/runThese.txt
```

**AddBiosyntheticFlags.sh**

```sh
#!/bin/bash
#####
# Unpacking the miniconda environment
#####
# replace env-name on the right hand side of this line with the name of your conda environment
ENVNAME=ICSCreator
# if you need the environment directory to be named something other than the environment name, change this line
ENVDIR=$ENVNAME

# these lines handle setting up the environment; you shouldn't have to modify them
export PATH
mkdir $ENVDIR
tar -xzf $ENVNAME.tar.gz -C $ENVDIR
. $ENVDIR/bin/activate
#extracting the files
genome="$1"
gunzip "$genome"_hmmer.out.gz
tar -xzvf "$genome".tar.gz
mv ./"$genome"/* ./
rm "$genome".tar.gz
#running the python program
rootDir=$(pwd)
python3 AddBiosyntheticFlags.py "$rootDir"
#compressing the results
mkdir "$genome"_biosynGBKs
mv *.gbk "$genome"_biosynGBKs
tar -czvf "$genome"_biosynGBKs.tar.gz "$genome"_biosynGBKs
```

**AddBiosyntheticFlags.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles
# in order to allow for better alignments when running BigScape, the code needs to add biosynthetic flags to
# assumes the protein domain, gbks to be changed, and the text file listing the pfam domains are all in the base directory
import os
import sys
import pandas as pd
import numpy as np

def GetFiles(rootDirectory):
    #getting the file paths
    gbks = []
    biosynDF = np.nan
    pfamsDF = np.nan
    for file in os.listdir(rootDirectory):
        #gbks
        if os.path.isfile(os.path.join(rootDirectory, file)) and not file.startswith("."):
            if file.endswith(".gbk") and not file.startswith('.'):
                gbks.append(os.path.join(rootDirectory, file))
            ### CHANGE THIS LINE IF THE FILE IS NAMED SOMETHING DIFFERENT
            elif file == "pfam_biosynthetic.txt":
                biosyntheticFilePath = os.path.join(rootDirectory, file)
                biosynDF = pd.read_csv(biosyntheticFilePath, sep = '\t')
            elif file.endswith("hmmer.out"):
                hmmerFilePath = os.path.join(rootDirectory, file)
                pfamsDF = pd.read_csv(hmmerFilePath, sep = '\t')
    results = [gbks, biosynDF, pfamsDF]
    return results

def AddBiosyntheticLabel(gbks, biosynDF, pfamsDF):
    ######
    # opening and parsing over the gbk file
    ######
    for gbkPath in gbks:
        editedGBK = []
        with open(gbkPath) as gbk:
            lines = gbk.readlines()
            #parsing over the gbk file
            insideCDS = False
            for index, line in enumerate(lines):
                #inside the CDS region if this is true
                if line.startswith("     CDS"):
                    insideCDS = True
                    editedGBK.append(line)
                    continue
                elif insideCDS == True:
                    editedGBK.append(line)
                    if line.startswith("                     /protein_id="):
                        protein = line.strip().split('=')[1][1:-1]
                        isBiosyn = IsBiosyntheticProtein(protein, pfamsDF, biosynDF)
                        if isBiosyn == True:
                            newLine = r'                     /gene_kind="biosynthetic"' + '\n'
                            editedGBK.append(newLine)
                    elif line.startswith("          "):
                        continue
                    #if this is true continue while checking for another protein ID
                    elif line.startswith("     CDS"):
                        continue
                    #the parsing has hit a new feature
                    else:
                        insideCDS = False
                else:
                    editedGBK.append(line)
        #overwritting the gbk file with the edited lines
        with open(gbkPath, "w") as gbkFinal:
            for l in editedGBK:
                gbkFinal.write(l)


def IsBiosyntheticProtein(protein, pfamsDF, biosynDF):
    pfamsOfProtein = pfamsDF.loc[pfamsDF.Target_Name == protein]
    #the protein had multiple hits
    if len(pfamsOfProtein) > 1:
        for index, row in pfamsOfProtein.iterrows():
            pfamTrimmed = row["Domain_Accession"].strip().split(".")[0]
            #checking if it is in the biosynDF
            if pfamTrimmed in biosynDF.Accession.values:
                return True
    #the protein only had one hit
    else:
        pfamTrimmed = pfamsOfProtein["Domain_Accession"].to_string(index=False).strip().split(".")[0]
        isBiosyn = pfamTrimmed in biosynDF.Accession.values
        return isBiosyn


if __name__ == '__main__':
    rootDirectory = sys.argv[1]
    results = GetFiles(rootDirectory)
    gbks = results[0]
    biosynDF = results[1]
    pfamsDF = results[2]
    AddBiosyntheticLabel(gbks, biosynDF, pfamsDF)
```

### Scripts that adds in biosynthetic tag to gbk files

**AddBiosyntheticFlags.sub**

```sh
# Specify the HTCondor Universe (vanilla is the default)
universe = vanilla
#
log = $(genome)_biosyn.log
error = $(genome)_biosyn.err
#
# Specify your executable (single binary or a script that runs several
executable = AddBiosyntheticFlags.sh
arguments = $(genome)
output = $(genome)_biosyn.out
#
#make sure it has gluster
requirements = (Target.HasCHTCStaging == true) && (OpSysMajorVer == 7)
#
# Specify that HTCondor should transfer files to and from the
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = /home/gnickles/FungalGenomes/$(genome)/$(genome)_hmmer.out.gz, /home/gnickles/BiosyntheticGenomeTag/gbks/$(genome).tar.gz, /home/gnickles/BiosyntheticGenomeTag/pfam_biosynthetic.txt, /home/gnickles/ICSCreator/ICSCreator.tar.gz, /home/gnickles/BiosyntheticGenomeTag/AddBiosyntheticFlags.py
transfer_output_files = $(genome)_biosynGBKs.tar.gz
#
# Tell HTCondor what amount of compute resources
request_cpus = 1
request_memory = 1GB
request_disk = 3GB
#
queue genome from /home/gnickles/BiosyntheticGenomeTag/runThese.txt

```

**AddBiosyntheticFlags.sh**

```sh
#!/bin/bash
#####
# Unpacking the miniconda environment
#####
# replace env-name on the right hand side of this line with the name of your conda environment
ENVNAME=ICSCreator
# if you need the environment directory to be named something other than the environment name, change this line
ENVDIR=$ENVNAME

# these lines handle setting up the environment; you shouldn't have to modify them
export PATH
mkdir $ENVDIR
tar -xzf $ENVNAME.tar.gz -C $ENVDIR
. $ENVDIR/bin/activate
#extracting the files
genome="$1"
gunzip "$genome"_hmmer.out.gz
tar -xzvf "$genome".tar.gz
mv ./"$genome"/* ./
rm "$genome".tar.gz
#running the python program
rootDir=$(pwd)
python3 AddBiosyntheticFlags.py "$rootDir"
#compressing the results
mkdir "$genome"_biosynGBKs
mv *.gbk "$genome"_biosynGBKs
tar -czvf "$genome"_biosynGBKs.tar.gz "$genome"_biosynGBKs
```

**AddBiosyntheticFlags.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles
# in order to allow for better alignments when running BigScape, the code needs to add biosynthetic flags to
# assumes the protein domain, gbks to be changed, and the text file listing the pfam domains are all in the base directory
import os
import sys
import pandas as pd
import numpy as np

def GetFiles(rootDirectory):
    #getting the file paths
    gbks = []
    biosynDF = np.nan
    pfamsDF = np.nan
    for file in os.listdir(rootDirectory):
        #gbks
        if os.path.isfile(os.path.join(rootDirectory, file)) and not file.startswith("."):
            if file.endswith(".gbk") and not file.startswith('.'):
                gbks.append(os.path.join(rootDirectory, file))
            ### CHANGE THIS LINE IF THE FILE IS NAMED SOMETHING DIFFERENT
            elif file == "pfam_biosynthetic.txt":
                biosyntheticFilePath = os.path.join(rootDirectory, file)
                biosynDF = pd.read_csv(biosyntheticFilePath, sep = '\t')
            elif file.endswith("hmmer.out"):
                hmmerFilePath = os.path.join(rootDirectory, file)
                pfamsDF = pd.read_csv(hmmerFilePath, sep = '\t')
    results = [gbks, biosynDF, pfamsDF]
    return results

def AddBiosyntheticLabel(gbks, biosynDF, pfamsDF):
    ######
    # opening and parsing over the gbk file
    ######
    for gbkPath in gbks:
        editedGBK = []
        with open(gbkPath) as gbk:
            lines = gbk.readlines()
            #parsing over the gbk file
            insideCDS = False
            for index, line in enumerate(lines):
                #inside the CDS region if this is true
                if line.startswith("     CDS"):
                    insideCDS = True
                    editedGBK.append(line)
                    continue
                elif insideCDS == True:
                    editedGBK.append(line)
                    if line.startswith("                     /protein_id="):
                        protein = line.strip().split('=')[1][1:-1]
                        isBiosyn = IsBiosyntheticProtein(protein, pfamsDF, biosynDF)
                        if isBiosyn == True:
                            newLine = r'                     /gene_kind="biosynthetic"' + '\n'
                            editedGBK.append(newLine)
                    elif line.startswith("          "):
                        continue
                    #if this is true continue while checking for another protein ID
                    elif line.startswith("     CDS"):
                        continue
                    #the parsing has hit a new feature
                    else:
                        insideCDS = False
                else:
                    editedGBK.append(line)
        #overwritting the gbk file with the edited lines
        with open(gbkPath, "w") as gbkFinal:
            for l in editedGBK:
                gbkFinal.write(l)

def IsBiosyntheticProtein(protein, pfamsDF, biosynDF):
    pfamsOfProtein = pfamsDF.loc[pfamsDF.Target_Name == protein]
    #the protein had multiple hits
    if len(pfamsOfProtein) > 1:
        for index, row in pfamsOfProtein.iterrows():
            pfamTrimmed = row["Domain_Accession"].strip().split(".")[0]
            #checking if it is in the biosynDF
            if pfamTrimmed in biosynDF.Accession.values:
                return True
    #the protein only had one hit
    else:
        pfamTrimmed = pfamsOfProtein["Domain_Accession"].to_string(index=False).strip().split(".")[0]
        isBiosyn = pfamTrimmed in biosynDF.Accession.values
        return isBiosyn

if __name__ == '__main__':
    rootDirectory = sys.argv[1]
    results = GetFiles(rootDirectory)
    gbks = results[0]
    biosynDF = results[1]
    pfamsDF = results[2]
    AddBiosyntheticLabel(gbks, biosynDF, pfamsDF)
```

This ran without any issues, the gbk results were then used to run my program on BigScape. 

## Step 4: Running BiG-SCAPE to generate the GCF classifications

### Condensed summary of the parameters used with BiG-SCAPE

```python
bigscape.py -i ./AllResults/ -o ./BigScape_ICS --cores 4 --anchorfile anchor_domainsNEW.txt --cutoffs 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 --clan_cutoff 0.35 0.5 --include_gbk_str .gbk
```

### The scripts

#### Removing the 9 files that contained errors

Before running this analysis I needed to check if there were any totally blank gbk/gff files. This occurs from a small number of misformatted files on NCBI breaking the code in Step 3. Because this is such an extreme minority of the total cluster predictions I will simply remove them from further analysis before moving forward.

**Files that had 0B gff files:**

```sh
(base) gnickles@Nymeria-Keller-Lab:/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ICSGenomes/clusterGffs$ ls -l | grep 0B
-rwxrwxrwx  1 gnickles  staff     0B Jan 26 12:12 GCA_000769735.1_289_Penicillium_expansum.gff*
-rwxrwxrwx  1 gnickles  staff     0B Jan 26 12:12 GCA_000769745.1_297_Penicillium_expansum.gff*
-rwxrwxrwx  1 gnickles  staff     0B Jan 26 12:12 GCA_012656225.1_1285_Aspergillus_fumigatiaffinis.gff*
-rwxrwxrwx  1 gnickles  staff     0B Jan 26 12:13 GCA_020500565.1_2043_Aspergillus_fumigatus.gff*
-rwxrwxrwx  1 gnickles  staff     0B Jan 26 12:13 GCA_020501885.1_2277_Aspergillus_fumigatus.gff*
-rwxrwxrwx  1 gnickles  staff     0B Jan 26 12:13 GCA_020503165.1_2464_Aspergillus_fumigatus.gff*
-rwxrwxrwx  1 gnickles  staff     0B Jan 26 12:13 GCA_020578095.1_3007_Parastagonospora_nodorum.gff*
-rwxrwxrwx  1 gnickles  staff     0B Jan 26 12:13 GCA_020580775.1_3123_Parastagonospora_nodorum.gff*
-rwxrwxrwx  1 gnickles  staff     0B Jan 26 12:13 GCF_000769745.1_3526_Penicillium_expansum.gff*
```

These will be removed before moving forward. The rest were uploaded to CHTC for use in the analysis

#### Resolving the 130+ genomes that failed on the supercomputer

In the initial run of BiG-SCAPE there were 130+ genomes with truncated gbk files due to prematurely ended runs on the UW-Madison supercomputer.

Example of the errors:

```sh
Importing GenBank files
   Error with file /home/input/./GCA_000001985.1/Cluster_1/1_Talaromyces_marneffei.gbk: 
    'Premature end of line during features table'
    (This file will be excluded from the analysis)
   Error with file /home/input/./GCA_000001985.1/Cluster_3/3_Talaromyces_marneffei.gbk: 
    'Premature end of line during features table'
    (This file will be excluded from the analysis)
   Error with file /home/input/./GCA_000001985.1/Cluster_4/4_Talaromyces_marneffei.gbk: 
    'Premature end of line during features table'
    (This file will be excluded from the analysis)
   Error with file /home/input/./GCA_000002715.1/Cluster_10/10_Aspergillus_clavatus.gbk: 
    'Premature end of line during features table'
    (This file will be excluded from the analysis)
```

The error ended up being from the supercomputer. Perplexingly if I individually downloaded a truncated genome and re-ran Step 3 locally the output was non-truncated. 

In the end, to resolve this niche issue I had to generate a bunch of small scripts that delt with organizing, re-naming, and downloading for these 130 problem genomes. I then ran this small subset locally. To prevent confusion from duplacte genomes of the same species having similar BGC predictions, **this script also changed the name of each BGC to a unique identifying number (1-3800), which can be found in the Supplmental Tables and GitHub repository.** It is unlikley that a reader of these Supplemental Methods will find these raw scripts useful as they are designed for my personal laptop and how the data was organized on my personal hard-drive (which is why they aren't extensively documented with comment lines like other code I made), but parts might be good to reference for reproducibility.

**3300GenomeOrganization.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles
#This script is going to be used to make some small functions that can organize everything in my folders

import os
import pdb
import shutil
import pandas as pd

keyFilePath = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ICSGenomes/KeysForClusters.txt'

def FindAndReplace(rootDirectory, gbkFolder, copyFrom):
    #super inefficient but I'm just finding and replacing any of the things that have issues and I re-ran on my local computer
    listOfGenomesDownloaded = []
    for genome in os.listdir(rootDirectory):
        if os.path.isdir(os.path.join(rootDirectory, genome)) and genome.startswith("GC"):
            genomePath = os.path.join(rootDirectory, genome)
            moveToNextGenome = False
            #stepping into the cluster folder
            for cluster in os.listdir(genomePath):
                if moveToNextGenome == True:
                    break
                if os.path.isdir(os.path.join(genomePath, cluster)) and cluster.startswith("Cluster"):
                    clusterPath = os.path.join(genomePath, cluster)
                    for file in os.listdir(clusterPath):
                        if os.path.isfile(os.path.join(clusterPath, file)) and file.endswith(".gbk") and not file.startswith("."):
                            gbkPath = os.path.join(clusterPath, file)
                            #grabbing the final line in the file
                            last_line = ""
                            with open(gbkPath) as gbk:
                                for line in gbk:
                                    pass
                                last_line = line
                            if "//" not in last_line:
                                print(genome)
                                listOfGenomesDownloaded.append(genome)
                                moveToNextGenome = True
                                break
    #######
    # Looping through the rootDirectory to find and replace the correct files
    #######
    keyFile = pd.read_csv(keyFilePath, sep = '\t')
    for genome in os.listdir(copyFrom):
        if os.path.isdir(os.path.join(copyFrom, genome)) and genome.startswith("GC") and genome in listOfGenomesDownloaded:
            genomePath = os.path.join(copyFrom, genome)
            #populating the dictionary that stores all of the paths for the files to be copied
            annotationFolder = os.path.join(genomePath, "Annotation_Windows")
            copyingDic = {}
            for file in os.listdir(annotationFolder):
                if os.path.isfile(os.path.join(annotationFolder, file)) and not file.startswith(".") and file.endswith(".fna"):
                    fnaNumber = str(file.split("_")[0])
                    fnaPath = os.path.join(annotationFolder, file)
                    copyingDic[fnaNumber] = [fnaPath]
            #adding in the gff and gbk files to the dictionary
            cassisFolder = os.path.join(annotationFolder, "CASSIS")
            for file in os.listdir(cassisFolder):
                if not file.startswith('.'):
                    filePath = os.path.join(cassisFolder, file)
                    clusterNumber = str(file.split("ion")[-1].split("_")[0])
                    copyingDic[clusterNumber].append(filePath)
            ######
            # Doing the deleting and copying
            ######
            for clusterNumber, filePaths in copyingDic.items():
                #getting the copyingFrom Paths
                copyingFromFNA = ""
                copyingFromGFF = ""
                copyingFromGBK = ""
                for path in filePaths:
                    if path.endswith(".fna"):
                        copyingFromFNA = path
                    elif path.endswith(".gbk"):
                        copyingFromGBK = path
                    else:
                        copyingFromGFF = path

                baseName = ".".join(copyingFromGFF.split("/")[-1].split(".")[0:-1])
                clusterName = keyFile.loc[keyFile.Full_Cluster_Name == baseName, "Number_Species"].to_string(index=False).strip()
                clusterNumber = clusterName.split("_")[0]
                #setting the new names
                fnaNewName = clusterNumber + "_" + "_".join(copyingFromFNA.split("/")[-1].split("_")[1:])
                gffNewName = clusterName + ".gff"
                gbkNewName = clusterName + ".gbk"
                #setting the copy to paths
                actualGenomePath = os.path.join(rootDirectory, genome)
                actualClusterPath = os.path.join(actualGenomePath, "Cluster_" + clusterNumber)
                newGffPath = os.path.join(actualClusterPath, gffNewName)
                newGbkPath = os.path.join(actualClusterPath, gbkNewName)
                consolidatedFolderGBK = os.path.join(gbkFolder, gbkNewName)
                newFnaPath = os.path.join(actualClusterPath, fnaNewName)
                #deleting all of the files in the folder stuff is getting copied to:
                cmd = "rm " + actualClusterPath + "/*.fna " + actualClusterPath + "/*.gff " + actualClusterPath + "/*.gbk"
                os.system(cmd)
                #copying over the fasta file
                cmd = "cp " + copyingFromFNA + " " + newFnaPath
                os.system(cmd)
                #copying over the gff file
                cmd = "cp " + copyingFromGFF + " " + newGffPath
                os.system(cmd)
                #copying the gbk file to the two places
                cmd = "cp " + copyingFromGBK + " " + newGbkPath
                os.system(cmd)
                cmd = "cp " + copyingFromGBK + " " + consolidatedFolderGBK
                os.system(cmd)


def DownloadBrokenBitches(rootDirectory, saveTo):
    listGenomesToDownload = []
    for genome in os.listdir(rootDirectory):
        if os.path.isdir(os.path.join(rootDirectory, genome)) and genome.startswith("GC"):
            genomePath = os.path.join(rootDirectory, genome)
            moveToNextGenome = False
            #stepping into the cluster folder
            for cluster in os.listdir(genomePath):
                if moveToNextGenome == True:
                    break
                if os.path.isdir(os.path.join(genomePath, cluster)) and cluster.startswith("Cluster"):
                    clusterPath = os.path.join(genomePath, cluster)
                    for file in os.listdir(clusterPath):
                        if os.path.isfile(os.path.join(clusterPath, file)) and file.endswith(".gbk") and not file.startswith("."):
                            gbkPath = os.path.join(clusterPath, file)
                            #grabbing the final line in the file
                            last_line = ""
                            with open(gbkPath) as gbk:
                                for line in gbk:
                                    pass
                                last_line = line
                            if "//" not in last_line:
                                print(genome)
                                listGenomesToDownload.append(genome)
                                moveToNextGenome = True
                                break
    #Downloading the correct files for the folders that need to be re-run from the supercomputer
    for genome in listGenomesToDownload:
        newFolder = os.path.join(saveTo, genome)
        os.mkdir(newFolder)
        #downloading the files
        gioFile = genome + "_GIO_ncbi_parent.txt.gz"
        hmmerFile = genome + "_hmmer.out.gz"
        gffFile = genome + ".gff.gz"
        fnaFile = genome + ".fna.gz"
        filesToCopy = [gioFile, hmmerFile, gffFile, fnaFile]
        for file in filesToCopy:
            cmd = "scp gnickles@submit1.chtc.wisc.edu:/home/gnickles/FungalGenomes/" + genome + "/" + file + " " + newFolder
            os.system(cmd)
        #unzipping the files
        for file in os.listdir(newFolder):
            if file.endswith(".gz") and not file.startswith("."):
                filePath = os.path.join(newFolder, file)
                cmd = "gunzip " + filePath
                os.system(cmd)
        #copying the cassis results from the HardDrive
        cassisFolder = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/CassisResults'
        genomeCassisResults = os.path.join(cassisFolder, "CASSIS_" + genome)
        cmd = "cp " + genomeCassisResults + "/gene* " + newFolder
        os.system(cmd)

def CopyFromSuperComputer(rootDirectory):
    for genome in os.listdir(rootDirectory):
        genomePath = os.path.join(rootDirectory, genome)
        newFolder = os.path.join(genomePath, "RawData")
        os.mkdir(newFolder)
        gioFile = genome + "_GIO_ncbi_parent.txt.gz"
        hmmerFile = genome + "_hmmer.out.gz"
        protFile = genome + "_prot.faa.gz"
        gffFile = genome + ".gff.gz"
        filesToCopy = [gioFile, hmmerFile, protFile, gffFile]
        for file in filesToCopy:
            cmd = "scp gnickles@submit1.chtc.wisc.edu:/home/gnickles/FungalGenomes/" + genome + "/" + file + " " + newFolder
            os.system(cmd)

#the function that does all of the magic. Creates the key file, organizes everything, and re-names every BGC to have a unique number along with the species name, i.e. 18_Candida_auris
def CreateKeyNames(rootDirectory, allAccessPath, gbkFolder):
    allAccess = pd.read_csv(allAccessPath, sep = '\t')
    counter = 1
    keyFile = pd.DataFrame(columns = ["Number_Only","Number_Species","Full_Cluster_Name","NCBI_Accession","Kingdom","Phylum","Class","Order","Family","Genus","Species","Genome_Information"])
    for genome in os.listdir(rootDirectory):
        print("Starting " + genome)
        genomePath = os.path.join(rootDirectory, genome)
        tSplit = GetKeyInfo(allAccess, genome)
        #looping over all of the cluster hits
        for cluster in os.listdir(genomePath):
            if os.path.isfile(os.path.join(genomePath, cluster)) and cluster.endswith(".fna") and not cluster.startswith("."):
                clusterNumber = cluster.split("_")[0]
                #making the folder to put all of the cluster stuff into it
                newFolder = os.path.join(genomePath, "Cluster_" + str(counter))
                os.mkdir(newFolder)
                newName = ""
                if tSplit[6] != "":
                    newName = str(counter) + "_" + tSplit[6]
                elif tSplit[7] != "":
                    newName = str(counter) + "_" + tSplit[7]
                else:
                    pdb.set_trace()
                if len(newName.split("_")) > 3 or "(" in newName or ")" in newName or "[" in newName or "]" in newName:
                    newName = input("The genome is " + genome + " and the ugly name is " + newName)
                newFastaName = str(counter) + "_" + "_".join(cluster.split("_")[1:])

                #######
                # Moving over the fna, gff, and gbk files to the cluster folder
                #######
                #GBK
                gbkName = "region" + clusterNumber + "_" + genome + "_Cluster.gbk"
                copyFrom = os.path.join(genomePath, gbkName)
                copyTo = os.path.join(gbkFolder, newName + ".gbk")
                cmd = "cp " + copyFrom + " " + copyTo
                os.system(cmd)
                moveTo = os.path.join(newFolder, newName + ".gbk")
                shutil.move(copyFrom, moveTo)
                #GFF
                gffName = "region" + clusterNumber + "_" + genome + "_Cluster.gff"
                gffPath = os.path.join(genomePath, gffName)
                moveTo = os.path.join(newFolder, newName + ".gff")
                shutil.move(gffPath, moveTo)
                #FNA
                fnaPath = os.path.join(genomePath, cluster)
                moveTo = os.path.join(newFolder, newFastaName)
                shutil.move(fnaPath, moveTo)
                #######
                # Adding in this entry into the key file
                #######
                newDataRow = {"Number_Only": counter,
                "Number_Species": newName,
                "Full_Cluster_Name": "region" + clusterNumber + "_" + genome + "_Cluster",
                "NCBI_Accession": genome,
                "Kingdom": tSplit[0],
                "Phylum": tSplit[1],
                "Class": tSplit[2],
                "Order": tSplit[3],
                "Family": tSplit[4],
                "Genus": tSplit[5],
                "Species": tSplit[6],
                "Genome_Information": tSplit[7]
                }
                keyFile = keyFile.append(newDataRow, ignore_index = True)
                counter += 1
    keyFilePath = os.path.join(rootDirectory, "KeysForClusters.txt")
    pdb.set_trace()
    keyFile.to_csv(keyFilePath, sep='\t', header = True, index = False)
def GetKeyInfo(allAccess, genome):
    genomeRow = allAccess.loc[allAccess.accession == genome]
    taxonomy = genomeRow['taxonomy'].to_string(index = False).strip()
    tSplit = taxonomy.split(",")
    return tSplit

#this function removes any of the cluster files, and re copyies over the raw ICS predictions. Allows the user to reset everything if something broke when running the CreateKeyNames function
def FixTheFolders(rootDirectory, predictionsDir):
    for genome in os.listdir(rootDirectory):
        if genome.startswith("GC"):
            clusterFound = False
            #Removing the cluster folders
            genomePath = os.path.join(rootDirectory, genome)
            for clusters in os.listdir(genomePath):
                if os.path.isdir(os.path.join(genomePath, clusters)) and clusters.startswith("Cluster_"):
                    clusterFound = True
                    #removing that folder
                    clusterPath = os.path.join(genomePath, clusters)
                    shutil.rmtree(clusterPath)
            if clusterFound == True:
                print("starting the cleaning of " + genome)
                for predictions in os.listdir(predictionsDir):
                    if os.path.isdir(os.path.join(predictionsDir, predictions)) and predictions.startswith(genome):
                        predictionPath = os.path.join(predictionsDir, predictions)
                        for file in os.listdir(predictionPath):
                            if not file.startswith("."):
                                copyFrom = os.path.join(predictionPath, file)
                                copyTo = os.path.join(genomePath, file)
                                cmd = "cp " + copyFrom + " " + copyTo
                                os.system(cmd)


if __name__ == '__main__':
    pd.options.display.max_colwidth = 5000
    rootDirectory = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ICSGenomes'
    allAccessPath = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/all_access.txt'
    gbkFolder = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ICSGenomes/gbks'
    predictionsDir = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ICSPredictions'
    saveTo = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ReRun'
    # CopyFromSuperComputer(rootDirectory)
    # CreateKeyNames(rootDirectory, allAccessPath, gbkFolder)
    # FixTheFolders(rootDirectory, predictionsDir)
    # DownloadBrokenBitches(rootDirectory, saveTo)
    FindAndReplace(rootDirectory, gbkFolder, saveTo)
# def MakeListOfHits():
```

#### Running Bigscape on the supercomputer

Once this issue with the genomes was resolved BigScape was able to read and understand the rest of the hits. However, it kept failing on my local computer as I didn't have nearly enough RAM to generate the output. 

Because of this I generated some supercomputer scripts that were able to do the same thing, but allowed me to have 100 GB of RAM.

**BigScape_ICS.sub**

```sh
# Specify the HTCondor Universe (vanilla is the default and is used
universe = docker
docker_image=aflatoxing/bigscape
#
+LongJob = true
#
log = bigscape.log
error = bigscape.err
#
# Specify your executable (single binary or a script that runs several
executable = BigScape_ICS.sh
output = bigscape.out
#
#make sure it has gluster
requirements = (Target.HasCHTCStaging == true) && (OpSysMajorVer == 7)
#
# Specify that HTCondor should transfer files to and from the
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = /home/gnickles/BigScape/anchor_domainsNEW.txt
transfer_output_files = BigScape_results.tar.gz
#
# Tell HTCondor what amount of compute resources
request_cpus = 4
request_memory = 100GB
request_disk = 10GB
#
# Tell HTCondor to run 3 instances of our job:
queue 1
```

- originally I tried running this while asking for 16 GB of memory (my computer has 16) but it still failed. It probably was overkill to ask for 100 GB but this is honestly nothing for CHTC to allocate. CPU's are much harder to come by hence why I still only asked for the same amount my laptop has. 

**BigScape_ICS.sh**

```sh
#!/bin/bash
#
cp /staging/gnickles/gbksBio.tar.gz ./
tar -xvzf gbksBio.tar.gz

bigscape.py -i ./AllResults/ -o ./BigScape_ICS --cores 4 --anchorfile anchor_domainsNEW.txt --cutoffs 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 --clan_cutoff 0.35 0.5 --include_gbk_str .gbk

tar -czvf BigScape_results.tar.gz ./BigScape_ICS
```

**In the end I found the best cuttoff to be 0.3, and will move forward using it in my analysis**

- Higher cuttoffs merged known ICS BGCS together when they should not be
- Lower cuttoffs only contained highly identical species and thus was not informative 
- Goal was to identify BGCs that are highly likley to produce the same or similar metabolites, while allow for enough leway to see what species might have it
- 0.3 strongly seemed like the sweet spot, although likley various cuttoffs would be more optimal for different individual GCFs, since there is little to no information to use as sanity checks outside of crmA, xanB, and ICSA/B this will be sufficient for my study

#### Adding metadata of taxonomy onto the bigscape network

I generated a script called **OverlayTax_BigscapeNetwork.py** whose only inputs are the accession file, which is automatically generated by the downloading and annotation scripts in step 1, and the network file from bigscape. This allows for better visualization of said network.

**OverlayTax_BigscapeNetwork.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles
#the network and accession represent the location on the computer for these files
import os
import sys
import pdb
import pandas as pd

def addMetaData(network, accession):
    #loading in the two files from their location on the computer
    networkDF = pd.read_csv(network, sep='\t', header=0)
    networkDF["TempName_S"] = None
    networkDF["TempName_T"] = None
    accessionDF = pd.read_csv(accession, sep="\t", names = ['number','TempName','sequence type', 'n1', 'n2', 'misc', 'SpeciesName'])
    TempAndSpecies = accessionDF[['TempName','SpeciesName']] #grabbing only these two columns
    #the species column that Mickey made uses spaces to separtate all of the information as follows
    #Kingdom, Division, Class, Order, Family, genus, species, genome information
    #ex: Fungi,Ascomycota,Eurotiomycetes,Eurotiales,Aspergillaceae,Aspergillus,Aspergillus_fumigatus,Aspergillus_fumigatus_A1163
    #source node
    metaDataDF_S = pd.DataFrame(columns = ["TempName_S","Kingdom_S","Phylum_S","Class_S","Order_S","Family_S","Genus_S","Species_S","Genome_Information_S"])
    for index, row in TempAndSpecies.iterrows():
        try:
            infoSplit = row["SpeciesName"].split(',')
        except:
            if str(row['SpeciesName']) == 'nan':
                print(row['TempName'] + " had no taxonomy information on it.")
                continue
        newDataRow = {"TempName_S": str(row["TempName"]),
        "Kingdom_S": str(infoSplit[0]).strip(),
        "Phylum_S": str(infoSplit[1]).strip(),
        "Class_S": str(infoSplit[2]).strip(),
        "Order_S": str(infoSplit[3]).strip(),
        "Family_S": str(infoSplit[4]).strip(),
        "Genus_S": str(infoSplit[5]).strip(),
        "Species_S": str(infoSplit[6]).strip(),
        "Genome_Information_S": str(infoSplit[7]).strip()
        }
        metaDataDF_S = metaDataDF_S.append(newDataRow, ignore_index=True)
    #target node
    metaDataDF_T = metaDataDF_S.copy()
    metaDataDF_T = metaDataDF_T.rename(columns=lambda x: '_'.join([x.split("_")[0],"T"]))

    #going through the network df and finding the row in df that matches up based on the temp name
    for index2, rowNetwork in networkDF.iterrows():
        #grabbing the temp name out of the network file, remove the region an kb parts in the name
        networkName_S = "_".join(rowNetwork["Clustername 1"].split("_")[1:3])
        networkDF.at[index2, "TempName_S"] = networkName_S
        #doing the same with the target nodes
        networkName_T = "_".join(rowNetwork["Clustername 2"].split("_")[1:3])
        networkDF.at[index2, "TempName_T"] = networkName_T
    #this is the line that does the magic, combines based on the TempData Column with the merge command
    # First adding metadata based on Source Node
    mergedTargetDF = networkDF.merge(right=metaDataDF_T, how='left', on="TempName_T")
    finalDF = mergedTargetDF.merge(right=metaDataDF_S, how='left', on="TempName_S")

    networkName = "MetaData_"+network.split(r"/")[-1]
    folderPathSplit = "/".join(network.split("/")[:-1]) #Removes file on the path
    networkPath = os.path.join(folderPathSplit, networkName)
    finalDF.to_csv(path_or_buf = networkPath, sep = "\t", header = True, index=False) #writing out the file

if __name__ == '__main__':
    networks = [r'/Volumes/HardDrive/ALLFungalGenomes_ICS/BigScapeResults/network_files/2021-07-15_12-36-44_hybrids_glocal/Others/Others_c0.30.network',r'/Volumes/HardDrive/ALLFungalGenomes_ICS/BigScapeResults/network_files/2021-07-15_12-36-44_hybrids_glocal/Others/Others_c0.35.network',r'/Volumes/HardDrive/ALLFungalGenomes_ICS/BigScapeResults/network_files/2021-07-15_12-36-44_hybrids_glocal/Others/Others_c0.40.network',r'/Volumes/HardDrive/ALLFungalGenomes_ICS/BigScapeResults/network_files/2021-07-15_12-36-44_hybrids_glocal/Others/Others_c0.45.network',r'/Volumes/HardDrive/ALLFungalGenomes_ICS/BigScapeResults/network_files/2021-07-15_12-36-44_hybrids_glocal/Others/Others_c0.50.network',r'/Volumes/HardDrive/ALLFungalGenomes_ICS/BigScapeResults/network_files/2021-07-15_12-36-44_hybrids_glocal/Others/Others_c0.55.network',r'/Volumes/HardDrive/ALLFungalGenomes_ICS/BigScapeResults/network_files/2021-07-15_12-36-44_hybrids_glocal/Others/Others_c0.60.network']

    accession = r"/Volumes/HardDrive/ALLFungalGenomes_ICS/all_access.txt"
    networkNumber = 0.30
    for n in networks:
        print("starting the network with cutoff value: " + str(networkNumber))
        addMetaData(n, accession)
        networkNumber += 0.05
```




## Step 5 Dereplicating the BiG-SCAPE clan output

The following scripts could be useful to anyone who wants to similary dereplicate a BiG-SCAPE generated network to account for bridge clusters. 

At a high level there are two major steps

1. Create what we decided to call an "Antiparsed" files, which is essentially a table summary of each BGC prediction. Makes it much easier to keep track of each BGC in an efficient manner that is bash friendly. 
2. Dereplicate the networks by checking for the connectivity between nodes, all while accounting for the bias that can occur from duplicate genomes of the same species.
   1. This is essentially done by looking at which nodes a given node is connected to, and if the neighbors of that node contain similar connection patterns
   2. If an extreme bridge cluster connection is identified, it is systematically analyzed either kept as a single GCF, or broken into multiple GCFs

**After these steps the user ends up with a folder called dereplicated that contains the final GCF classifications, with bridge clusters removed.**

### Step 5.1 Scripts for generating the Antiparsed files

Input files to make the python script work: 

**AntiParsedMaker.sub**

```sh
# Specify the HTCondor Universe (vanilla is the default)
universe = vanilla
#
log = $(genome)_anti.log
error = $(genome)_anti.err
#
# Specify your executable (single binary or a script that runs several
executable = AntiParsedMaker.sh
arguments = $(genome)
output = $(genome)_anti.out
#
#make sure it has gluster
requirements = (Target.HasCHTCStaging == true) && (OpSysMajorVer == 7)
#
# Specify that HTCondor should transfer files to and from the
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = /home/gnickles/AntiParsed/compressedGffs/$(genome)_gffs.tar.gz,/home/gnickles/FungalGenomes/$(genome)/$(genome)_GIO_ncbi_parent.txt.gz, /home/gnickles/ICSCreator/ICSCreator.tar.gz, /home/gnickles/AntiParsed/AntiParsedMaker_CHTC.py
transfer_output_files = $(genome)_anti.tar.gz
#
# Tell HTCondor what amount of compute resources
request_cpus = 1
request_memory = 1GB
request_disk = 3GB
#
queue genome from /home/gnickles/AntiParsed/antiQueue.txt
```

**AntiParsedMaker.sh**

```sh
#!/bin/bash
#####
# Unpacking the miniconda environment
#####
# replace env-name on the right hand side of this line with the name of your conda environment
ENVNAME=ICSCreator
# if you need the environment directory to be named something other than the environment name, change this line
ENVDIR=$ENVNAME

# these lines handle setting up the environment; you shouldn't have to modify them
export PATH
mkdir $ENVDIR
tar -xzf $ENVNAME.tar.gz -C $ENVDIR
. $ENVDIR/bin/activate
#####
# Extracting the files
#####
genome="$1"
tar -xzvf "$genome"_gffs.tar.gz
mv ./"$genome"_gffs/*.gff ./
gunzip "$genome"_GIO_ncbi_parent.txt.gz
#####
# Running the program
#####
rootDir=$(pwd)
python3 AntiParsedMaker_CHTC.py "$rootDir" "$genome"
#
tar -czvf "$genome"_anti.tar.gz ./anti_parsed_new
```

**AntiParsedMaker_CHTC.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles
#usage: AntiParsedMaker [directory containing genomes] [Type of cluster ie. NRPS, Other, ICS etc.]
#Argument Notes: The rootDirectory and clusterType are taken from the GBKCreator file.
#If running on its own all of the complete GBK files need to be in their own folder preferably called Window[bp amount]All
#and there needs to be GIO files with the same genome name inside of them as the cluster files
import os
import sys
import pandas as pd

#Files I need GIO for genome, gff of the cluster predictions
def GetTheProteinName(rootDirectory, geneName):
    for files in os.listdir(rootDirectory):
        if os.path.isfile(os.path.join(rootDirectory, files)) and "GIO_ncbi_parent.txt" in files:
            try:
                gioPath = os.path.join(rootDirectory, files)
                gioDF = pd.read_csv(gioPath, sep = "\t", names =['contig', 'protein', 'rna','start','stop','direction','key'], engine = 'c')
                # geneLoci = geneName.split("gene-")[-1]
                proteinNameRows = gioDF[gioDF['key'].str.contains(geneName)]
                # this happens when there are more than one match, we need to select the row that has the EXACT match
                if len(proteinNameRows) > 1:
                    for index, row in proteinNameRows.iterrows():
                        keyInfo = row['key']
                        geneOnly = keyInfo.split('gene-')[-1].split(";")[0].strip()
                        checkThisGene = "gene-" + geneOnly
                        if checkThisGene == geneName:
                            proteinNameRow = row
                            break
                else:
                    proteinNameRow = proteinNameRows
                try:
                    protein = proteinNameRow['protein'].to_string(index = False).strip()
                    return protein
                except:
                    protein = proteinNameRow['protein']
                    return protein
            except:
                pdb.set_trace()
                return "ERROR IN FILE"

def AntiParsedMaker(rootDirectory, clusterType, genome):
    #now the program will step through all of the genome folders
    print(genome)
    errorInFile = False
    #Opening the parsed dataframe
    parsed = pd.DataFrame(columns = ["count","cluster","contig","protein","gene","start","stop","+/-","key","cluster type"])
    count = 1
    for files in os.listdir(rootDirectory):
        if errorInFile == True:
            break
        if os.path.isfile(os.path.join(rootDirectory, files)) and files.endswith('.gff'):
            gffPath = os.path.join(rootDirectory, files)
            # there was a single genome that I had to add that encoding line to get it to work, not really sure why but it was different than the other gff files in its encoding
            gffDF = pd.read_csv(gffPath, names = ['contig', 'source', 'annotation type', 'start', 'stop', 'n', 'direction', 'n2', 'key'], encoding='cp1252', sep = '\t')
            gffGenesOnly = gffDF[gffDF['annotation type'] == 'gene']
            for index, row in gffGenesOnly.iterrows():
                keyInfo = row['key']
                geneName = keyInfo.split('ID=')[-1].split(';')[0]
                protein = GetTheProteinName(rootDirectory, geneName)
                if protein == "ERROR IN FILE":
                    print("there is an issure with: " + genome)
                    errorInFile = True
                    break
                #checking if this is the backbone gene
                clusterTypeAddToRow = clusterType
                if 'backbone=True' in keyInfo:
                    clusterTypeAddToRow = clusterType + " " + clusterType
                #getting the cluster number
                clusterName = files.split("_")[2:3][0]
                #adding in the new data row
                newDataRow = {"count": count,
                "cluster": clusterName,
                "contig": row['contig'],
                'protein': protein,
                'gene' : geneName,
                'start' : row['start'],
                'stop' : row['stop'],
                '+/-' : row['direction'],
                'key' : row['key'],
                'cluster type': clusterTypeAddToRow
                }
                parsed = parsed.append(newDataRow, ignore_index=True)
        count += 1
    #Generating the output
    parsedFileName = genome + "_parsed_anti"
    antiFolder = os.path.join(rootDirectory, "anti_parsed_new")
    if not os.path.isdir(os.path.join(rootDirectory, "anti_parsed_new")):
        os.mkdir(antiFolder)

    newFilePath = os.path.join(antiFolder, parsedFileName)
    # pandas to csv overwrites existing files, this is not the most efficient way to do this but it will work
    parsed.to_csv(newFilePath, sep = "\t", header = False, index=False)


if __name__ == '__main__':
    #If running on its own
    rootDirectory = sys.argv[1]
    genome = sys.argv[2] #TODO on sh script
    clusterType = "ICS"

    AntiParsedMaker(rootDirectory, clusterType, genome)
```

### Step 5.2: Scripts for dereplicating the BigScape clan output

The derepication pipeline, which was written by Dr. Milton Drott, was leveraged for this section. It takes the bigscape network and clan predictions, to see what is connected to what. Based on these connections it tries to remove outlier connections (i.e. one node connection linking two real GCFs). It was edited slightly to work on my file structure and laptop.

To run this the user has to have a folder called `anti_parsed-new` (this can be changed if the file path is changed in the scripts) in the same directory as these scripts. You the user tell it which network and clan file to use.

**Step1.sh**

```sh
#Step 1
#assumes that you have the parsed_new files in a folder called anti_parsed_new within the current directory
#get gene counts
#GN-EDIT: ./anti_parsed_new/*_parsed_new_hybrids_anti_fix --> *_parsed_anti
cat ./anti_parsed_new/*_parsed_anti | awk -F'\t' '{print $2"\t"$10}' | awk '{if ($3=="") print $1"\t"$2"\t"1; else print $1"\t"$2"\t"1"\t"1}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | sort -k1,1> temp2
#
#get bb gene counts
#GN-EDIT: ./anti_parsed_new/*_parsed_new_hybrids_anti_fix --> *_parsed_anti
cat ./anti_parsed_new/*_parsed_anti | awk -F'\t' '{print $2"\t"$10}' | awk '{if ($3=="") print $1"\t"$2"\t"1; else print $1"\t"$2"\t"1"\t"1}' | awk -F '\t' '{a[$1"\t"$2"\t"]+= $4} END{for (i in a) print i, a[i]}' | sort -k1,1 > temp
#
join temp temp2 | awk '{print $1"\t"$2"\t"$4"\t"$3}'> all_bgc_types_with_counts
#
#GN-EDIT: Changed this to the directory of the network file
Network='/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/BigScapeResults/Run_BiosynGBKs/EditedNetwork_Clans/EDITED_Others_clans_0.35_0.50.tsv'
#Network=$1
#get all bgcs
#
#note the “grep mF” may have ot be changed for other projects
#GN-EDIT:changed grep mf to grep GC
awk 'NR>1 {print $1"\n"$2}' "$Network" | grep GC | sort -u > allbgcs
#
Counter="1"
REPLY=`head -1 allbgcs`
#
Rlength=`echo -e "$REPLY" | wc -c| sed 's/ //g'`
#
while [ "$Rlength" -gt 2 ]; do
#
#get nearest neighbors and find their nearest neighbors
grep -w "$REPLY" "$Network" | awk '{print $1"\n"$2}' | sort -u > round1
#
grep -w -f round1 "$Network" | awk '{print $1"\n"$2}' | sort -u > round2
#
#determine if nearest neighbors have more clusters than second nearest neighbors
Count1=`cat round1 | wc -l | sed 's/ //g'`
Count2=`cat round2 | wc -l | sed 's/ //g'`
#
echo -e "$Counter""\t""$Count1""\t""$Count2"
#
#if nearest and second nearest are diff continue on until nth neighbor is same as n+1
while [ "$Count1" -lt "$Count2" ]; do
#
grep -w -f round2 "$Network"| awk '{print $1"\n"$2}' | sort -u > round1
grep -w -f round1 "$Network"| awk '{print $1"\n"$2}' | sort -u > round2
#
Count1=`cat round1 | wc -l | sed 's/ //g'`
Count2=`cat round2 | wc -l | sed 's/ //g'`
#
echo -e "$Counter""\t""$Count1""\t""$Count2"
#
done
#
cat round2 > "$Counter"family
#
#
#remove the BGCs that just got nearest neighbored above
#note the “grep mF” may have to be changed for other projects
#GN-EDIT: grep mF --> grep GC
cat "$Counter"family "$Counter"family allbgcs | sort | uniq -u | grep GC> temp
mv temp allbgcs
#
#grab out the backbone gene types for counterfamily
#
grep -w -f "$Counter"family all_bgc_types_with_counts > "$Counter"family_with_type
#
#lets identify rare bgc combos:
awk '{print $2}' "$Counter"family_with_type | awk '{count[$0]++} END {for (type in count) print type, count[type]}' > temp
#
#and count their occurances
awk 'FNR==NR{s+=$2;next;} {printf "%s\t%s\t%s\n",$1,$2,100*$2/s}' temp temp > zemp
#
#if the total count is 5 or bigger and a bgc type is present in fewer than 20% add a warning
awk 'FNR==NR{SUM+=$2;next} {if (SUM>=5 && $3<20){print $1}}' zemp zemp >warn_bb
while read;do cat "$Counter"family_with_type | awk -v var="$REPLY" '{if ($2==var) print $1"\t"$2"\t"$3"\t"$4"\t""DANGER!"; else print}'> zemp2;mv zemp2 "$Counter"family_with_type;done <warn_bb
rm "$Counter"family
#
echo -e "$Counter""\t""$Count1""\t""$Count2"
#
REPLY=`head -1 allbgcs`
Rlength=`echo -e "$REPLY" | wc -c| sed 's/ //g'`
Counter=$[Counter + 1]
done
mkdir family_with_type
rm temp* warn_bb zemp round*
mv *family_with_type ./family_with_type/
```

**Step2.sh**

```sh
#Step 2
#this script will print same_gcfs and finalsecond rm errors every time, this was done out of an abundance of caution as
#on the first run if these are present they could be problematic.
#
#Network="/Volumes/HardDrive/ICSProject/Eurotiales_all/BigScape/7500AllCopy/network/Others_c0.30.network"
#Clans="/Volumes/HardDrive/ICSProject/Eurotiales_all/BigScape/7500AllCopy/clans/Others_clans_0.30_0.70.tsv"
Network='/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/BigScapeResults/Run_BiosynGBKs/EditedNetwork_Clans/EDITED_Others_c0.30.network'
Clans='/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/BigScapeResults/Run_BiosynGBKs/EditedNetwork_Clans/EDITED_Others_clans_0.35_0.50.tsv'
#
mkdir neighbor_info
#where I is the number of dereplicated clusters if you want to do all fo them)
#GN-EDIT: works with the LStep2.sh script for automation
# for i in $(seq 1 "$3");do
#### THIS IS THE LINE TO EDIT
for i in {1..210};do
echo "$i"
#
rm *nn
#
awk '{print $1}' ./family_with_type/"$i"family_with_type > iso_in_network
#
grep -w -f iso_in_network "$Clans" | awk '{print $3,$1}'| sed 's/_alignment.fasta.*>/ /' > iso_with_gcf
#
awk '{print $1}' iso_with_gcf | sort -u > all_gcfs
#
grep -w -f iso_in_network "$Network" > reduced_net
#
#grab out all the isos in a gcf and place them and their nearest neighbors into a file
while read;do grep -w "$REPLY" iso_with_gcf | awk '{print $2}' >tempz; grep -w -f tempz reduced_net | awk '{print $1"\n"$2}' | sort -u > "$REPLY"nn;done < <(cat all_gcfs)
#
#
while read;do
for f in *nn;do
target=`echo -e "$f"|sed 's/nn//'`
#
grep -w "$REPLY" iso_with_gcf| awk '{print $2}' > source_iso
grep -w "$target" iso_with_gcf| awk '{print $2}' > target_iso
#
sourcecount=`cat source_iso | wc -l | sed 's/ //g'`;
targetcount=`cat target_iso | wc -l | sed 's/ //g'`;
#
common=`cat "$REPLY"nn "$f" | sort | uniq -d | wc -l | sed 's/ //g'`;
common_source=`cat "$REPLY"nn "$f" | sort | uniq -d | grep -w -f source_iso |wc -l | sed 's/ //g'`;
common_target=`cat "$REPLY"nn "$f" | sort | uniq -d | grep -w -f target_iso |wc -l | sed 's/ //g'`;
source_nn_count=`cat "$REPLY"nn | wc -l | sed 's/ //g'`;
target_nn_count=`cat "$f" |  wc -l | sed 's/ //g'`;
echo -e "$REPLY""nn""\t""$f""\t""$sourcecount""\t""$source_nn_count""\t""$common_source""\t""$common""\t""$common_target""\t""$targetcount""\t""$target_nn_count";
done;
done < <(cat all_gcfs) > neighbors
#
#
#Filter neighbors to require proportions of all columns to be greater than 15% -- now changed to 50% for various reasons
#the higher this filter is the the less likley the program will concider something as like
#
awk '$7>0 && $1!=$2 {print}' neighbors | awk '{portion_source=$5/$3}{portion_target=$7/$8}{portion_nn_source=$6/$4}{portion_nn_target=$6/$9} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"|",portion_source,portion_nn_source,portion_target,portion_nn_target}' | awk '$11>0.15 && $12>0.15 && $13>0.15 && $14>0.15 {print}' | column -t > filtered_neighbors
#
awk '$7>0 && $1!=$2 {print}' neighbors | awk '{portion_source=$5/$3}{portion_target=$7/$8}{portion_nn_source=$6/$4}{portion_nn_target=$6/$9} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"|",portion_source,portion_nn_source,portion_target,portion_nn_target}' | column -t > z
mv z neighbors
#
rm same_gcfs
#
cat all_gcfs >temp_gcfs
gcf_count=`cat all_gcfs|  wc -l | sed 's/ //g'`
#
while [ $gcf_count -gt 0 ];do
#
gcf=`head -1 temp_gcfs`
#
grep -w "$gcf""nn" filtered_neighbors | awk '{print $1"\n"$2}'| sort -u > same
echo -e "$gcf""nn" >> same
#
#allow for there to be none in the same group but still add back in without increasing count
echo -e "$gcf""nn" >> same
cat same | sort -u > tempz
mv tempz same
#
grep -w -f same filtered_neighbors | awk '{print $1"\n"$2}'| sort -u > secondsame
#
Count1=`cat same | wc -l | sed 's/ //g'`
Count2=`cat secondsame | wc -l | sed 's/ //g'`
#
while [ "$Count1" -lt "$Count2" ]; do
#
grep -w -f secondsame filtered_neighbors | awk '{print $1"\n"$2}'| sort -u > same
grep -w -f same filtered_neighbors | awk '{print $1"\n"$2}'| sort -u > secondsame
#
Count1=`cat same | wc -l | sed 's/ //g'`
Count2=`cat secondsame | wc -l | sed 's/ //g'`
#
done
#
cat same secondsame| sed 's/nn//' | sort -u > tempz
mv tempz same
#
firstsame=`cat same | sed 's/$/,/' | tr -d '\n'`
#
echo -e "$firstsame" >> same_gcfs
#
cat same temp_gcfs | sort | uniq -u > tempz
mv tempz temp_gcfs
gcf_count=`cat temp_gcfs|  wc -l | sed 's/ //g'`
done
#tidy1
rm finalsecond secondsame same temp_gcfs
#
#
mv neighbors ./neighbor_info/"$i"neighbors
mv filtered_neighbors ./neighbor_info/"$i"filtered_neighbors
mv same_gcfs ./neighbor_info/"$i"same_gcfs
#
#tidy2
rm reduced_net all_gcfs iso_with_gcf iso_in_network
done
```

**Step3.sh**

```sh
#STEP3 WITH DEREPLICATION OF SPECIES DURING/BEFORE MODAL BGC DETERMINATION
#ASSUMES YOU HAVE HYRBIDS ANTI PARSED
#will print various errors about awk not finding zzemp2 and failing to remove various directories – haven’t gone over/tidied this yet.
mkdir dereplicated
#
grep hybrid ./anti_parsed_new/*parsed_anti | awk -F'\t' '{print $2,$10}' |sort -u> all_bgc_hybrids
#Network="/Volumes/HardDrive/ICSProject/Eurotiales_all/BigScape/Window7500All/bigscape/network_files/2021-03-09_09-37-08_hybrids_glocal/Others/Others_c0.30.network"
#Clans="/Volumes/HardDrive/ICSProject/Eurotiales_all/BigScape/Window7500All/bigscape/network_files/2021-03-09_09-37-08_hybrids_glocal/Others/Others_clans_0.30_0.70.tsv"
Network='/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/BigScapeResults/Run_BiosynGBKs/EditedNetwork_Clans/EDITED_Others_c0.30.network'
Clans='/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/BigScapeResults/Run_BiosynGBKs/EditedNetwork_Clans/EDITED_Others_clans_0.35_0.50.tsv'
master_access="/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/all_access.txt"
#
#while read parent_gcf;do
#GN-EDIT: changed for the LStep3.sh script
# for parent_gcf in $(seq 1 "$3");do
for parent_gcf in {1..210};do
#for parent_gcf in {1..24};do
echo -e starting "$parent_gcf"
#
awk '{print $1}' ./family_with_type/"$parent_gcf"family_with_type > iso_in_network
grep -w -f iso_in_network "$Clans"| awk '{print $3,$1}'| sed 's/_alignment.fasta.*>/ /' > iso_with_gcf
#
#output the combined gcfs into a single file
Counter="1"
while read;do
echo "$REPLY" | tr ',' '\n'|grep . > gcfs
grep -w -f gcfs iso_with_gcf | awk '{print $2}' > ttemp
grep -w -f ttemp ./family_with_type/"$parent_gcf"family_with_type | awk '{print $1"\t"$2"\t"$3"\t"$4}' > "$Counter"_combo"$parent_gcf"_withtype
    Counter=$[Counter + 1]
done < <(cat ./neighbor_info/"$parent_gcf"same_gcfs)
#
#
#then iterate through..
ls -lh *_combo"$parent_gcf"_withtype | awk '{print $9}' | sed 's/_combo.*//' > combo_gcfs
while read -r combo_gcf;do
#
#add in hybrid designation
#
while read -r bgc type gene_count bb_count; do
myhybrid=`grep -w "$bgc" ./all_bgc_hybrids  | grep hybrid| wc -l | sed 's/ //g'`
hybrid=`grep "$bgc" ./all_bgc_hybrids | wc -l | sed 's/ //g'`
#
echo -e "$myhybrid""\t""$hybrid" | awk -F'\t' -v bgc="$bgc" -v type="$type" -v gene_count="$gene_count" -v bb_count="$bb_count" '{if ($1==1) print bgc"\t"type"\t"gene_count"\t"bb_count"\t""hybrid"; if ($1==0 && $2==1) print bgc"\t"type"\t"gene_count"\t"bb_count"\t""maybe_hybrid"; if ($1==0 && $2==0) print bgc"\t"type"\t"gene_count"\t"bb_count"\t""not_hybrid"}'
done < <(cat "$combo_gcf"_combo"$parent_gcf"_withtype)> bemp
mv bemp "$combo_gcf"_combo"$parent_gcf"_withtype
#
#
#lets identify rare bgc combos:
#count bb type occurances
awk '{print $2}' "$combo_gcf"_combo"$parent_gcf"_withtype | awk '{count[$0]++} END {for (type in count) print type, count[type]}' > ttemp
#
#and count their occurrences in percent
awk 'FNR==NR{s+=$2;next;} {printf "%s\t%s\t%s\n",$1,$2,100*$2/s}' ttemp ttemp > zzemp
#
#if the total count is 5 or bigger and a bgc type is present in fewer than 20% add a warning
awk 'FNR==NR{SUM+=$2;next} {if (SUM>=5 && $3<20){print $1}}' zzemp zzemp >warn_bb
#
# below will add warnings to the ones that have weird bbgene content etc. based on 20% above.
while read;do cat "$combo_gcf"_combo"$parent_gcf"_withtype| awk -v var="$REPLY" '{if ($2==var) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""DANGER_type!!"}'>> zzemp2;done <warn_bb
awk '{print $1}' zzemp2 > zzemp3
awk '{print $1}' "$combo_gcf"_combo"$parent_gcf"_withtype > names
cat names zzemp3 | sort | uniq -u > add_back
#at this point add_back contains all of the isolates that did not get a type warning
#
#grab out gene counts and add WARNINGS
grep -w -f add_back "$combo_gcf"_combo"$parent_gcf"_withtype | sort -n -k3,3 | awk '{print $3}' > ttemp
#add warnings to those that have 50% more or less gene than the mean gene count.
Mean=`cat ttemp | awk '{sum+=$1}END {print sum/NR}' `
#
UpperLim=`echo -e "$Mean" | awk '{half=$1/2} {print $1+half}'`
LowerLim=`echo -e "$Mean" | awk '{half=$1/2} {print $1-half}'`
#
grep -w -f add_back "$combo_gcf"_combo"$parent_gcf"_withtype | awk -v var="$UpperLim" '{if ($3>var) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""DANGER_count!!"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5}'| awk -v var="$LowerLim" '{if ($3<var) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""DANGER_count!!"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' >>zzemp2
rm ttemp
#
# GN-EDIT: this change was made so the grep -v lines were unable to remove any of the warning files. This is becuase i want to not remove any of the warnings files.
# mv zzemp2 "$parent_gcf"combo_"$combo_gcf"_withtype_warnings
cat zzemp2 | sed 's/DANGER/MANGER/g' > "$parent_gcf"combo_"$combo_gcf"_withtype_warnings
#
# and for backbone gene counts… too many may be multiple clusters together, too few may be fragmented
grep DANGER "$parent_gcf"combo_"$combo_gcf"_withtype_warnings > add_back
#
Modebb_count=`grep -v DANGER "$parent_gcf"combo_"$combo_gcf"_withtype_warnings | awk '{a[$4]++} END {for (i in a) if (a[i] > freq) {most=i; freq=a[i]} printf("%s\n", most)}'`
#
grep -v DANGER "$parent_gcf"combo_"$combo_gcf"_withtype_warnings | awk -v var="$Modebb_count" '{if ($4!=var) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""DANGER_bb_count!!"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > zzemp2
cat add_back >>zzemp2
mv zzemp2 "$parent_gcf"combo_"$combo_gcf"_withtype_warnings
#
#Okay now what about dereplicating ones in the same species..
grep -v DANGER "$parent_gcf"combo_"$combo_gcf"_withtype_warnings  | awk '{print $1}' | sed 's/_/ /g' | awk '{print $2"_"$3}' > isos
's/GC.*/GC/'
#
rm "$parent_gcf"iso_spp"$combo_gcf"
while read -r iso;do
access=`grep -w "$iso" "$master_access" |awk '{print $2}'`;
#
Access_check=`echo -e "$access" | wc -c | sed 's/ //g'`
#
#
if [ $Access_check -gt 2 ];then
#
species=`grep -w "$access" "$master_access"| awk -F'\t' '{print $7}' | awk -F',' '{print $7}'`;
#GN-EDITS: removed the -w in the grep command because using accesssions and not mF nomenclature
# removing the danger grep section
grep "$iso" "$parent_gcf"combo_"$combo_gcf"_withtype_warnings | grep -v DANGER| awk -v var="$species" '{print var"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' >> "$parent_gcf"iso_spp"$combo_gcf"
# grep "$iso" "$parent_gcf"combo_"$combo_gcf"_withtype_warnings | awk -v var="$species" '{print var"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' >> "$parent_gcf"iso_spp"$combo_gcf"
#
else
#this would indicate that the same accession is occurring more than once.
echo -e  something odd with "$iso" in "$combo_gcf" from parent gcf "$parent_gcf"
#
fi
done< <(cat isos)
#
awk '{print $1}' "$parent_gcf"iso_spp"$combo_gcf" | sort -u> species
#
#loop through each set 11 times to generate a modal bbgene value for each species
mode_counter="1"
while [ $mode_counter -lt 11 ]; do
#
#this randomly selects one representative of eachspecies and then calculates modal bb type across species (repeating many times to make sure we get the #right answer even if multiple reps of same species.
#
while read -r species;do grep -w "$species" "$parent_gcf"iso_spp"$combo_gcf" | shuf -n 1 ;done< <(cat species) | awk '{a[$3]++} END {for (i in a) if (a[i] > freq) {most=i; freq=a[i]} printf("%s\n", most)}'
    mode_counter=$[mode_counter + 1]
#
done > modal_bb
#
#problem here is that if frequencies are even this will not take a random one..
final_mode_bb=`awk '{a[$1]++} END {for (i in a) if (a[i] > freq) {most=i; freq=a[i]} printf("%s\n", most)}' modal_bb`
#
awk -v var="$final_mode_bb" '$3==var {print}' "$parent_gcf"iso_spp"$combo_gcf" > "$parent_gcf"iso_spp"$combo_gcf"_modal
#
#
#The following will take the iso spp and select a single representative from each species containing modal bb
#identify the unique isolates
awk '{print $1}' "$parent_gcf"iso_spp"$combo_gcf"_modal | sort -u > "$parent_gcf"iso_spp_uniq
while read;do awk -v var="$REPLY" '$1==var {print}' "$parent_gcf"iso_spp"$combo_gcf"_modal | shuf -n 1;done < <(cat "$parent_gcf"iso_spp_uniq) > "$parent_gcf"cluster_spp_dereplicated"$combo_gcf"
#
echo -e "$combo_gcf""\t""$final_mode_bb" >> gcf_summary_ttemp
#
#tidy
rm "$combo_gcf"_combo"$parent_gcf"_withtype "$parent_gcf"iso_spp"$combo_gcf"_modal
#
#
mv "$parent_gcf"iso_spp"$combo_gcf" ./dereplicated/"$parent_gcf"iso_spp_NEW"$combo_gcf"
mv "$parent_gcf"cluster_spp_dereplicated*"$combo_gcf" ./dereplicated/"$parent_gcf"cluster_spp_NEW_dereplicated"$combo_gcf"
mv "$parent_gcf"combo_"$combo_gcf"_withtype_warnings ./dereplicated/"$parent_gcf"combo_"$combo_gcf"_withtype_warnings_NEW
#
done < <(cat combo_gcfs)
#
#add warnings if multiple files share a modal backbone as these maybeshould be combined…
while read -r gcfnum bb;do
echo -e "$bb" | tr ',' '\n' |grep . > bbs
awk -v var="$gcfnum" '$2!=var {print}' gcf_summary_ttemp > others
species_count=`cat ./dereplicated/"$parent_gcf"cluster_spp_NEW_dereplicated"$gcfnum" | wc -l | sed 's/ //g' `
#warning=`grep -w -f bbs others | wc -l | sed 's/ //g'|awk '{if  ($1>=2) print "WARN"; else print "fine"}'`
echo -e "$gcfnum""\t""$bb""\t""$species_count"
done < <(cat gcf_summary_ttemp) > "$parent_gcf"gcf_summary
#
#
mv "$parent_gcf"gcf_summary ./dereplicated/"$parent_gcf"gcf_summary_NEW
#
#tidy
rm spp_dereplicated isos add_back names zzemp3 warn_bb zzemp ttemp "$parent_gcf"iso_spp_uniq bbs others gcf_summary_ttemp "$parent_gcf"same_gcfs isos
done
```

## Step 6: Getting protein domain statistics on the files

In order to do this in the most efficient way, I wanted to generate a large tsv table that contained a row for each unique protein domain in every gene in every cluster. This allowed for very easy subsetting, sorting, and analyzing in some of the more downstream analysis.

**GenerateLargeTable.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles (gnickles@wisc.edu)
#Usage: this program is designed to work with the output from Dr. Milton Drott's depreplication pipeline. This MUST be run before using this script for it to work seamlessly.
#Assumptions: there is a folder called dereplication that contains all of the input data, and the user will provide where the output should be saved to (the folder does not need to be made before running the code)

import os
import sys
import pdb
import pandas as pd
import numpy as np

##########
# GLOBAL VARIABLES
##########
# These are run here and not the main() since they need to be called in so many subfunctions
keyFile = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ICSGenomes/KeysForClusters.txt'
keyDF = pd.read_csv(keyFile, sep='\t')
# generating the large dataframe that will store all of the information
col_names = ['GCF', 'GCC', 'Cluster_Name', 'NCBI_Accession', 'Protein', 'Domain', 'Domain_Accession', 'E-value_hmmer', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Genome_Information']
megaDF = pd.DataFrame(columns = col_names)


#Subfunction for GenerateLargeTable that gets the taxonomic information based on a cluster key name
def GetTaxInfo(clusterName):
    clusterRow = keyDF.loc[keyDF.Number_Species == clusterName]
    neededColumns = clusterRow[['NCBI_Accession', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Genome_Information']]
    neededSeries = neededColumns.squeeze()
    resultDic = neededSeries.to_dict()
    return resultDic

#Subfunction for GenerateLargeTable to locates the raw data path for this ncbi accession; assumes folder stucture is Accession/RawData/ i.e. GCF_000002655/RawData
def FindGenomePath(ncbiAccession, genomesFolder):
    for genomes in os.listdir(genomesFolder):
        if os.path.isdir(os.path.join(genomesFolder, genomes)) and genomes == ncbiAccession:
            genomePath = os.path.join(genomesFolder, genomes)
            dataPath = os.path.join(genomePath, "RawData")
            return [genomePath, dataPath]

#loops over all of the cluster tables to populate the mega informational table, this uses the combo_with_warnings file to do so
def GenerateLargeTable(dereplicatedFolder, genomesFolder, outputDir):
    print("Starting the population of the megaDF")
    #looping over the dereplication folder
    counter = 1
    for combo in os.listdir(dereplicatedFolder):
        #only stepping in to with warnings files, as these have not removed any of the clusters, this will be done by my own code later on
        if os.path.isfile(os.path.join(dereplicatedFolder, combo)) and not combo.startswith(".") and "withtype_warnings_NEW" in combo and "combo" in combo:
            comboPath = os.path.join(dereplicatedFolder, combo)
            clanNum = combo.split("_")[1]
            famNum = combo.split("combo")[0]
            print("GCF: " + famNum + " GCC: " + clanNum)
            #reading in the combo file as a tsv, the variable is called GCC as it represents a gene cluster clan
            gccDF = pd.read_csv(comboPath, sep='\t', names=['Cluster_Name',"BGC_Type","Num_Genes", "Num_Backbones","Hyrid","Warnings"])
            for index, row in gccDF.iterrows():
                #making the row dictionary, which will later be added to the large DF once fully populated
                clusterName = row['Cluster_Name']
                ## CHECKING IF THERE IS AN ISSUE IN THE KEY FILE FOR THIS ENTRY
                numsInCluster = sum(c.isdigit() for c in clusterName)
                if clusterName[numsInCluster:(numsInCluster+1)] != "_":
                    print("ISSUE: " + combo + " " + clusterName)
                    continue
                taxInfo = GetTaxInfo(clusterName) # returns a dictionary to be merged later on
                ncbiAccession = taxInfo['NCBI_Accession']

                newDataRow = {'GCF': famNum,
                'GCC': clanNum,
                'Cluster_Name': clusterName
                }
                newDataRow.update(taxInfo)

                #Getting the proteins, protein domains, and adding the rows
                results = FindGenomePath(ncbiAccession, genomesFolder)
                genomePath = results[0]
                dataPath = results[1]
                AddProteinInformation(newDataRow, genomePath, dataPath)
        if counter == 5:
            megaOutPath = os.path.join(outputDir,  "IntermediateICSPredictions.tsv")
            megaDF.to_csv(megaOutPath, index=False, sep='\t')
            counter = 1
        else:
            counter += 1

    #####
    # Writing out the results
    #####
    megaOutPath = os.path.join(outputDir,  "AllICSPredictions.tsv")
    megaDF.to_csv(megaOutPath, index=False, sep='\t')

#program that does the heavy lifting for adding in the protein and pfam information to the megaDF
def AddProteinInformation(newDataRow, genomePath, dataPath):
    #####
    # Reading in the hmmer file as a df
    #####
    hmmerDf = np.nan
    for file in os.listdir(dataPath):
        if os.path.isfile(os.path.join(dataPath, file)) and file.endswith("hmmer.out"):
            hmmerPath = os.path.join(dataPath, file)
            hmmerDf = pd.read_csv(hmmerPath, sep='\t')
    #####
    # Finding the gff for this cluster and also reading it as a df
    #####
    clusterName = newDataRow['Cluster_Name']
    clusterNumber = clusterName.split("_")[0]
    gffDf = np.nan
    for folder in os.listdir(genomePath):
        if os.path.isdir(os.path.join(genomePath, folder)) and ("Cluster_" + str(clusterNumber)) == folder:
            clusterPath = os.path.join(genomePath, folder)
            #now that we are in the correct cluster folder, we need to get the path to the gff file
            #there should only be one in this folder and it should correspond to the BGC prediction
            for file in os.listdir(clusterPath):
                if os.path.isfile(os.path.join(clusterPath, file)) and file.endswith(".gff") and not file.startswith("."):
                    gffPath = os.path.join(clusterPath, file)
                    col_names = ['contig','source','annotation_type','start','stop','n/a','direction','na','info']
                    gffDf = pd.read_csv(gffPath, names=col_names, sep='\t')
                    break
    #####
    # Looping over the gff for each protein
    #####
    cdsRegions = gffDf.loc[gffDf.annotation_type == "CDS",'info'].tolist()
    proteinsAlreadyAdded = []
    for region in cdsRegions:
        splitInfo = region.split(";")
        for i in splitInfo:
            if "protein_id" in i:
                protein = i.split("=")[-1].strip()
                if protein not in proteinsAlreadyAdded:
                    #####
                    # Getting the domains that associate with this protein
                    #####
                    GetDomains(protein, hmmerDf, newDataRow)
                    proteinsAlreadyAdded.append(protein)

def GetDomains(protein, hmmerDf, newDataRow):
    #this needs to be added to update the global variable in the function
    global megaDF
    domainsDf = hmmerDf.loc[hmmerDf.Target_Name == protein]
    #if this conditional is true, there were no prediction protein domains for this protein
    if len(domainsDf) == 0:
        appendThis = {"Domain": np.nan,
        "Domain_Accession": np.nan,
        'E-value_hmmer': np.nan,
        "Protein": protein}
        newDataRow.update(appendThis)
        # Adding in the complete data row to the megaTable
        megaDF = megaDF.append(newDataRow,ignore_index=True)
    else:
        for index, row in domainsDf.iterrows():
            domain = row['Domain']
            domainAccession = row['Domain_Accession']
            evalue = row['E-value']
            appendThis = {"Domain": domain,
            "Domain_Accession": domainAccession,
            'E-value_hmmer': evalue,
            "Protein": protein}
            #adding in this domain information
            newDataRow.update(appendThis)
            # Adding in the complete data row to the megaTable
            megaDF = megaDF.append(newDataRow,ignore_index=True)

    return megaDF

if __name__ == '__main__':
    dereplicatedFolder = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/DereplicationOfHits/dereplicated'
    genomesFolder = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ICSGenomes'
    outputDir = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ConservationAnalysis'
    if not os.path.isdir(outputDir):
        os.mkdir(outputDir)

    GenerateLargeTable(dereplicatedFolder, genomesFolder, outputDir)

```

The output of this is the "megaDf", aka the table that contains ALL of the information. Then I needed to do some data cleaning. So I created a script called `CleanTheTable.py` that can filter out the table depending on how many unique species or genomes are in a given GCF or GCC. 

**CleanTheTable.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles (gnickles@wisc.edu)
#Usage: This program is designed to be run after generating the generate large table script
# It will remove small groups or groups without enough unique species, and will remove outlier clusters

import os
import sys
import pdb
import pandas as pd
import numpy as np
import re

##########
# GLOBAL VARIABLES
##########
col_names = ['FilteringBy', 'Name', 'NumGenomes', 'NumSpecies', 'pORf', 'ListGenomes', "ListSpecies"]
cleaningTable = pd.DataFrame(columns = col_names)

megaDFPath = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ConservationAnalysis/AllICSPredictions.tsv'
megaDF = pd.read_csv(megaDFPath, sep='\t')
#this means that any clan or family with less than 5 genomes or 3 unique species will be removed from the analysis; these values can be changed based on the user's discretion
genomeMin = 0
uniqSpeciesMin = 3

def CountGenomesAndSpecies(columnParameter, filteringList, outputDir):
    global cleaningTable
    passedList = []
    for f in filteringList:
        filteredDF = megaDF.loc[megaDF[columnParameter] == f]
        uniqGenomes = filteredDF.NCBI_Accession.unique()
        numUniqGenomes = len(uniqGenomes)
        #doing the same for the number of unique species
        filteredDF['RealSpecies'] = filteredDF['Cluster_Name'].str.replace(r'^[0-9]+_', '', case = False, regex=True)
        uniqSpecies = filteredDF.RealSpecies.unique()
        numUniqSpecies = len(uniqSpecies)
        ###
        # Adding in the results to the cleaningTable
        ###
        pORf = ""
        if numUniqGenomes <= genomeMin or numUniqSpecies <= uniqSpeciesMin:
            pORf = "Failed"
        else:
            pORf = "Passed"
            passedList.append(f)
            print(str(f) + " passed!")
        newDataRow = {"FilteringBy": columnParameter,
        "Name": f,
        "NumGenomes": numUniqGenomes,
        "NumSpecies": numUniqSpecies,
        "pORf": pORf,
        "ListGenomes": ",".join(uniqGenomes),
        "ListSpecies": ",".join(uniqSpecies)
        }
        cleaningTable = cleaningTable.append(newDataRow, ignore_index = True)
    return passedList

def RemoveSmallGroups(outputDir):
    global cleaningTable
    global megaDF

    #constructing the filtered table based on the families, i.e removing only families that have less than the cuttoff criteria
    familiesList = megaDF.GCF.unique()
    print("FAMILIES FILTERING")
    passedFamilies = CountGenomesAndSpecies('GCF', familiesList, outputDir)

    #fist I have to convert the two columns to string type
    megaDF['GCF'] = megaDF['GCF'].astype('str')
    megaDF['GCC'] = megaDF['GCC'].astype('str')
    #constructing the filtered table based clans, i.e removing only clans that have less than the cutoff criteria
    megaDF['GCF_GCC'] = megaDF.GCF.str.cat(megaDF.GCC, sep="_")
    clansList = megaDF.GCF_GCC.unique()
    print("\nCLANS FILTERING")
    passedClans = CountGenomesAndSpecies('GCF_GCC', clansList, outputDir)
    ###
    # Writing out the cleaningTable
    ###
    cleaningPath = os.path.join(outputDir, "FilteringTable.tsv")
    cleaningTable.to_csv(cleaningPath, sep='\t', index=False)

    # Returning the list of passed groups
    return [passedFamilies, passedClans]

def RemoveOutlierClusters(passedFamilies, passedClans, allFamilyDF, allClanDF):
    #FAMILY FIRST
    print("\nREMOVING THE OUTLIER CLUSTERS IF THEY ONLY HAVE ONE GENE IN THE PREDICTION\n")
    print("FAMILIES")
    #first a simple loop that will remove any of the rows with only a single protein in them. I don't want them to skew the domain percents later on
    allClusters = allFamilyDF.Cluster_Name.unique()
    for cluster in allClusters:
        clusterDF = allFamilyDF.loc[allFamilyDF.Cluster_Name == cluster]
        proteinsInCluster = clusterDF.Protein.unique()
        #this condition should never be satisfied but I have a set trace here to make sure nothing has gone wrong
        if len(proteinsInCluster) == 0:
            pdb.set_trace()
        elif len(proteinsInCluster) == 1:
            #there is only one protein in this cluster, so remove it from the df
            print(cluster + " only had 1 protein")
            allFamilyDF = allFamilyDF.loc[allFamilyDF.Cluster_Name != cluster]
        else:
            continue
    #adding in the column that contains the number of proteins in each cluster
    allFamilyDF = AddProteinNumColumn(allClusters, allFamilyDF)

    #NOW CLANS
    print("\nCLANS")
    #first a simple loop that will remove any of the rows with only a single protein in them. I don't want them to skew the domain percents later on
    allClusters = allClanDF.Cluster_Name.unique()
    for cluster in allClusters:
        clusterDF = allClanDF.loc[allClanDF.Cluster_Name == cluster]
        proteinsInCluster = clusterDF.Protein.unique()
        #this condition should never be satisfied but I have a set trace here to make sure nothing has gone wrong
        if len(proteinsInCluster) == 0:
            pdb.set_trace()
        elif len(proteinsInCluster) == 1:
            #there is only one protein in this cluster, so remove it from the df
            print(cluster + " only had 1 protein")
            allClanDF = allClanDF.loc[allClanDF.Cluster_Name != cluster]
        else:
            continue
    #adding in the column that contains the number of proteins in each cluster
    allClanDF = AddProteinNumColumn(allClusters, allClanDF)

    #now a more complicated loop, this one will remove any clusters that are more than 2 standard deviations outside the mean for that given GCC or GCF
    #doing the family groupings first
    print("\nREMOVING CLUSTERS WITH >= 2 STDS OF GENES OUTSIDE OF THE MEAN GENE COUNT\n")
    print("GCFs:")
    for family in passedFamilies:
        allFamilyDF = GetStdDeviations(family, 'GCF', allFamilyDF)
    #now for the clan groupings
    print("\nGCCs:")
    for clan in passedClans:
        allClanDF = GetStdDeviations(clan, 'GCF_GCC', allClanDF)

    # RETURNING THE RESULTS
    return [allFamilyDF, allClanDF]

#subfunction for RemoveOutlierClusters that generates the stds for each group, and removes the outliers from the table
def GetStdDeviations(groupName, filteringColumn, fullDF):
    groupDF = fullDF[fullDF[filteringColumn] == str(groupName)]
    clustersInGroup = groupDF.Cluster_Name.unique()

    numGenesList = []
    for cluster in clustersInGroup:
        clusterDF = groupDF.loc[groupDF.Cluster_Name == cluster]
        numGenes = clusterDF.Num_Genes.unique()
        #this condition should never be fulfilled but is here as a safety check
        if len(numGenes) == 0 or len(numGenes) > 1:
            print("Something went wrong!")
            pdb.set_trace()
        else:
            numGenesList.append(numGenes[0])

    #Now getting the standard deviations for this list
    std = np.std(numGenesList)
    meanCount = np.mean(numGenesList)
    twoStdsAbove = meanCount + (std * 2)
    twoStdsBelow = meanCount - (std * 2)

    #finding and removing columns that deviate from this range
    for cluster in clustersInGroup:
        clusterDF = groupDF.loc[groupDF.Cluster_Name == cluster]
        numGenes = clusterDF.Num_Genes.unique()
        if numGenes >= twoStdsAbove or numGenes <= twoStdsBelow:
            print(str(groupName) + "  Mean: " + str(round(meanCount, 2)) + "  " + cluster + ": " + str(int(numGenes[0])))
            fullDF = fullDF.loc[fullDF.Cluster_Name != cluster]
    return fullDF


#subfunction for RemoveOutlierClusters that adds a new column to the megaDF with the number of genes in a given cluster
def AddProteinNumColumn(allClusters, groupDF):
    for cluster in allClusters:
        clusterDF = groupDF.loc[groupDF.Cluster_Name == cluster]
        proteinsInCluster = clusterDF.Protein.unique()
        numGenes = len(proteinsInCluster)
        groupDF.loc[groupDF.Cluster_Name == cluster, "Num_Genes"] = numGenes
    return groupDF

def RemoveFailedGroups(passedFamilies, passedClans):
    passedFamilies = [str(i) for i in passedFamilies]
    allFamilyDF = megaDF[megaDF['GCF'].isin(passedFamilies)]

    allClanDF = megaDF[megaDF['GCF_GCC'].isin(passedClans)]
    return [allFamilyDF, allClanDF]

if __name__ == '__main__':
    #computer specific path variables! Must be set by the user before running
    outputDir = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ConservationAnalysis'
    gbks = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ICSGenomes/gbks'

    print("Filtering based on genome amount and the number of unique species:\n")
    #The passed families is first spot, clans is second spot in the list
    passed = RemoveSmallGroups(outputDir)
    passedFamilies = passed[0]
    passedClans = passed[1]
    #removing the failed groups from the dataframes and splitting it into two versions, one for the family classification and the other the clan
    results = RemoveFailedGroups(passedFamilies, passedClans)
    allFamilyDF = results[0]
    allClanDF = results[1]
    famLenBefore = len(allFamilyDF.Cluster_Name.unique())
    clanLenBefore = len(allClanDF.Cluster_Name.unique())
    #removing all of the outliers from each dataframe seperatly and returns both edited dataframes
    results = RemoveOutlierClusters(passedFamilies, passedClans, allFamilyDF, allClanDF)
    allFamilyDF = results[0]
    allClanDF = results[1]
    famLenAfter = len(allFamilyDF.Cluster_Name.unique())
    clanLenAfter = len(allClanDF.Cluster_Name.unique())
    print("\nNUM OF CLUSTERS ELIMINATED FROM REMOVING THE OUTLIERS\nFamily dataframe\n Before: " + str(famLenBefore) + "  After: " + str(famLenAfter) + "\nClan dataframe\n Before: " + str(clanLenBefore) + "  After: " + str(clanLenAfter))
    #saving the edited dataframes as their own tables
    famPath = os.path.join(outputDir, "FilteredFamilies.tsv")
    clanPath = os.path.join(outputDir, "FilteredClans.tsv")

    allFamilyDF.to_csv(famPath, sep='\t', index=False)
    allClanDF.to_csv(clanPath, sep='\t', index=False)
```

And finally this program just organized the gbk files into folders based on the dereplicated calls which sets it up for step 10.

**Organize GBKsFromTables.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles (gnickles@wisc.edu)
#Usage: This program should be run after running CleanTheTable.py. The user provides it the path to the gbks

import os
import pdb
import pandas as pd
import numpy as np
import shutil

def MoveTheCluster(gbkFolder, outputDir, cluster):
    for file in os.listdir(gbkFolder):
        if os.path.isfile(os.path.join(gbkFolder, file)) and cluster in file and not file.startswith("."):
            moveFrom = os.path.join(gbkFolder, file)
            moveTo = os.path.join(outputDir, file)
            cmd = "cp " + moveFrom + " " + moveTo
            os.system(cmd)
            break

def OrganizeTheGBKs(gbkFolder, clanDF, familyDF, outputDir):
    #Family First
    print("moving the GCFs")
    famOutDir = os.path.join(outputDir, "GCFs")
    if os.path.isdir(famOutDir):
        shutil.rmtree(famOutDir)
    os.mkdir(famOutDir)

    famList = familyDF.GCF.unique()
    for fam in famList:
        gcfDir = os.path.join(famOutDir, str(fam))
        os.mkdir(gcfDir)
        print(fam)
        gcfDF = familyDF.loc[familyDF.GCF == fam]
        clustersList = gcfDF.Cluster_Name.unique()
        for cluster in clustersList:
            MoveTheCluster(gbkFolder, gcfDir, cluster)

    #Now the Clans
    print("moving the GCCs")
    clanOutDir = os.path.join(outputDir, "GCCs")
    if os.path.isdir(clanOutDir):
        shutil.rmtree(clanOutDir)
    os.mkdir(clanOutDir)

    clanList = clanDF.GCF_GCC.unique()
    for clan in clanList:
        gccDir = os.path.join(clanOutDir, clan)
        os.mkdir(gccDir)
        print(clan)
        gccDF = clanDF.loc[clanDF.GCF_GCC == clan]
        clustersList = gccDF.Cluster_Name.unique()
        for cluster in clustersList:
            MoveTheCluster(gbkFolder, gccDir, cluster)


if __name__ == '__main__':
    gbkFolder = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ICSGenomes/gbks'
    clanPath = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ConservationAnalysis/FilteredClans.tsv'
    clanDF = pd.read_csv(clanPath, sep='\t')
    familyPath = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ConservationAnalysis/FilteredFamilies.tsv'
    familyDF = pd.read_csv(familyPath, sep='\t')
    #if this is not made the program will make it for you
    outputDir = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ConservationAnalysis/GCF_GCCs'
    if not os.path.isdir(outputDir):
        os.mkdir(outputDir)

    OrganizeTheGBKs(gbkFolder, clanDF, familyDF, outputDir)
```

## Step 7: Running clinker on the outputs

This step will be run on the CHTC supercomputer, as there are too many to resonably run on my local computer. And the large one's require a heavy amount of RAM (100GB +). The following is the submit script (the memory was adjusted up and down depending on the size of the GCF or GCC).

**RunClinker.sub**

```sh
# Specify the HTCondor Universe (vanilla is the default and is used
universe = vanilla
#
log = $(group)_clinker.log
error = $(group)_clinker.err
#log = small2_clinker.log
#error = small2_clinker.err
#
# Specify your executable (single binary or a script that runs several
executable = RunClinker.sh
arguments = $(group)
output = $(group)_clinker.out
# output = small2_clinker.out
#
#make sure it has gluster
requirements = (Target.HasCHTCStaging == true) && (OpSysMajorVer == 7)
#
# Specify that HTCondor should transfer files to and from the
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
#may also need to do one with the fixed… although maybe the answer to that is to write a script that figures out if the fixed has a better busco than the #other.. FILE IS all_fungal_busco.txt
transfer_input_files = /home/gnickles/Clinker/RawData/largeFiles/$(group).tar.gz
transfer_output_files = $(group)_clinker.tar.gz
#
# Tell HTCondor what amount of compute resources
request_cpus = 4
request_memory = 70GB
request_disk = 5GB
#
# Tell HTCondor to get names from txt file:
queue group from /home/gnickles/Clinker/RawData/largeFiles/largeFiles.txt
```

And this is the bash script that did all of the magic. 

**RunClinker.sh**

```sh
#!/bin/bash
#####
# Unpacking the miniconda environment
#####
cp /staging/gnickles/clinker.tar.gz ./
# replace env-name on the right hand side of this line with the name of your conda environment
ENVNAME=clinker
# if you need the environment directory to be named something other than the environment name, change this line
ENVDIR=$ENVNAME

# these lines handle setting up the environment; you shouldn't have to modify them
export PATH
mkdir $ENVDIR
tar -xzf $ENVNAME.tar.gz -C $ENVDIR
. $ENVDIR/bin/activate

mkdir gitFiles
cd gitFiles/
git clone https://github.com/gamcil/clinker.git
cd clinker/
pip install .
cd ..
cd ..

#####
# Extracting my files
#####
tar -xzvf "$1".tar.gz

#####
# Extracting my files
#####
clinker ./"$1"/*.gbk -s "$1".json -p "$1".html -j 4

mkdir "$1"_clinker
mv "$1".json "$1".html "$1"_clinker
tar -czvf "$1"_clinker.tar.gz "$1"_clinker
```

After this, the outputs were downloaded and placed onto my computer.

## Step 8: Investigating gene conservation and domain conservation

In order to garner confidence in my predictions I need to rely heavily on gene conservation of the BGC across evolutionary time. However, this output needs to be interpreted and anlyzed by human eyes, and different species are more related to each other than others, certain domains are key indicators of gene clusters, etc. So I created a script that allows the user to generate summary tables of any GCC or GCF they are curious about. I also made a script that finds what GCC or GCF a given species or BGC is in. This will better inform others in the lab if they want to dive deeper into a specific species, GCF/GCC, etc.

**FindTheGCF.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles (gnickles@wisc.edu)

#this program will loop over the consolidated folder to grab out the ICS backbone genes

import os
import sys
import pdb
import codecs
import pandas as pd


# "region1_GCF_000002655.1"
# "region1_GCA_000002855.2"

##########
# THIS MUST BE SET BEFORE RUNNING THE CODE ON SPECIES SEARCHES
##########
keyFile = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ICSGenomes/KeysForClusters.txt'

#subfunction made for species searching
def CheckIfClusterInCombo(comboFilePath, searchingFor, speciesOrAccession):
    #getting a list of the species that have a gene cluster prediction
    keyFileDF = pd.read_csv(keyFile, sep='\t')
    #filtering by the searchFor string
    clusterPredictions = ""
    if speciesOrAccession.lower() == "species":
        clusterPredictions = keyFileDF.loc[keyFileDF.Number_Species.str.contains(searchingFor), "Number_Species"].tolist()
    elif speciesOrAccession.lower() == "accession":
        clusterPredictions = keyFileDF.loc[keyFileDF.NCBI_Accession.str.contains(searchingFor), "Number_Species"].tolist()
    found = False
    returnedResults = []
    comboDF = pd.read_csv(comboFilePath, names=['BGC', "type","1", "2", "3", "4"], sep='\t')
    for cluster in clusterPredictions:
        checkDF = comboDF.loc[comboDF.BGC == cluster.strip()]
        if len(checkDF) > 0:
            returnedResults.append([cluster, True])
        else:
            returnedResults.append([cluster, False])
    #if it gets to this point the clusters were not in the file
    return returnedResults

#this subfunction just returns a boolean depending on if the cluster is in the file
def comboCheckIfClusterInComboAccession(comboFilePath, searchingFor):
    #getting a list of the species that have a gene cluster prediction
    keyFileDF = pd.read_csv(keyFile, sep='\t')
    #filtering by the searchFor string
    clusterPredictions = keyFileDF.loc[keyFileDF.NCBI_Accession.str.contains(searchingFor), "Number_Species"].tolist()
    returnedResults = []
    comboDF = pd.read_csv(comboFilePath, names=['BGC', "type","1", "2", "3", "4"], sep='\t')
    accesionGFCs = comboDF.loc
    #I added this codecs line since I was getting encoding errors, I found how to do it here:
    # https://stackoverflow.com/questions/12468179/unicodedecodeerror-utf8-codec-cant-decode-byte-0x9c
    with codecs.open(comboFilePath, 'r', encoding='utf-8', errors='ignore') as combo:
        if searchingFor in combo.read():
            return True
        #otherwise
        return False

#the main function that finds if a search term in is what gene cluster family/clan predictions
def FindTheGCF(dereplicatedFolder, searchingFor, speciesOrAccession):
    printedTopLine = False
    for files in os.listdir(dereplicatedFolder):
        if os.path.isfile(os.path.join(dereplicatedFolder, files)) and 'withtype_warnings_NEW' in files:
            comboFilePath = os.path.join(dereplicatedFolder, files)
            isInFile = CheckIfClusterInCombo(comboFilePath, searchingFor, speciesOrAccession)
            for searchTerm in isInFile:
                cluster = searchTerm[0]
                found = searchTerm[1]
                if printedTopLine == False:
                    print("Search hits for: " + searchingFor)
                    printedTopLine = True
                elif found == True and printedTopLine == True:
                    print(cluster + " was found in this combo file:\n" + comboFilePath)


#this searches through all of the files, only should be used if checking for bugs
def FindTheGCF_ALL(dereplicatedFolder, searchingFor):
    if speciesOrAccession.lower() == "accession":
        for files in os.listdir(dereplicatedFolder):
            if os.path.isfile(os.path.join(dereplicatedFolder, files)):
                comboFilePath = os.path.join(dereplicatedFolder, files)
                isInFile = CheckIfClusterIncombo(comboFilePath, searchingFor)
                if isInFile == True:
                    print("The file(s) containing " + searchingFor + " is in: \n" + comboFilePath)
    elif speciesOrAccession.lower() == "species":
        print("This function is not complete!")
        #TODO

dereplicatedFolder = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/DereplicationOfHits/dereplicated'
# "Aspergillus_fumigatus"
# "Penicillium_expansum"
searchingFor = "GCA_000011425.1"
speciesOrAccession = 'accession'
#can add in ability to search at any taxonomic rank if someone feels like that would be useful; just chat to me
FindTheGCF(dereplicatedFolder, searchingFor, speciesOrAccession)
```

And this is the script that conducts the similarity analysis based on the clinker results:

**SimilarityAnalysis.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles
#usage: Parses the output of clinker, assumes this was already run and the output file has the same name as the containing folder with a tsv extension, and filters out things that don't fit the given criteria
#THIS IS MADE TO BE RUN ON THE SUPERCOMPUTER
import os
from sys import argv
import pandas as pd
import shutil
from collections import Counter #needed for finding the duplicate species
import random # needed for grabbing random subsamples from lists
import numpy as np #used in the GetCombinations function
from operator import mul
from functools import reduce
import json

############
# GLOBAL VARIABLES
############
cuttoff = 0.5
numGenesThreshold = 3

############
# FUNCTIONS
############
# Main function: assumes the tsv files from clinker have the same name as the containing folder with a tsv extension
def AnalyzeClan(clan, clanPath, outputDirectory):
    #creating the dataframes
    results = CreateTheDataFrames(clanPath)
    final_combos = results[0]
    fullDF = results[1]
    fileFound = False
    #Reading in the json file
    for file in os.listdir(clanPath):
        #this condition should only be true for one file in the folder!!
        if os.path.isfile(os.path.join(clanPath, file)) and file.endswith(".json") and not file.startswith("."):
            jsonPath = os.path.join(clanPath, file)
            fileFound = True
            with open(jsonPath) as j:
                data = json.load(j)
                #generate the uid map
                map = MapTheJson(data)
                #match the uids to the correct place in the dataframe
                fullDF = PopulateGroups(map, data['groups'], fullDF)
                #writing out the fullDF for referencing later
                fullPath = os.path.join(outputDirectory, clan + "_fullDF.tsv")
                fullDF.to_csv(fullPath, sep='\t', index=False)
    #clinker needs to be run on this folder, if everything goes smoothly this condition will never be met
    if fileFound == False:
        print("Clinker still needs to be run on the following folder: " + clan)
    #Running the analysis on the dataframes to determine if it passes the cuttoff threshold
    pORf_results = PassOrFail(fullDF, final_combos, outputDirectory, clan)
    pORf = pORf_results[0]
    dfToWrite = pORf_results[-1]
    #writing out the files with a print statement on if they passed
    if pORf == True:
        print(clan + " passed!")
    else:
        print(clan + " failed :'(")
    WriteOutFiles(outputDirectory, clanPath, clan, dfToWrite, pORf)

#this function writes out all of the things to the output directory, and makes the sub clan folder
def WriteOutFiles(outputDirectory, clanPath, clan, clanDF, pORf):
    #writing out the tsv file with the gene comparisons
    if pORf:
        tsvOutputPath = os.path.join(outputDirectory, clan + "_" + str(cuttoff) + "Similarity_PASSED.tsv")
        clanDF.to_csv(tsvOutputPath, sep='\t',header=True, index=True, na_rep = "NaN")
    else:
        tsvOutputPath = os.path.join(outputDirectory, clan + "_" + str(cuttoff) + "Similarity_FALIED.tsv")
        clanDF.to_csv(tsvOutputPath, sep='\t',header=True, index=True, na_rep = "NaN")

#Run the analysis on clansDFs to see if it will pass
def PassOrFail(fullDF, final_combos, outputDirectory, clan):
    #as long as every pairwise combination returns the same conclusion this condition will be set to true
    allCombosTheSame = False
    resultsFromEachCombo = []
    returnThisClanDF = ""
    returnClan_alreadyEdited = False
    #setting the BGC_index as the index column
    fullDF = fullDF.set_index('BGC_names')
    fullDF.index = fullDF.index.astype('str')
    ####
    # setting the df for returning if the whole thing passes
    ####
    returnThisClanDF = fullDF.filter(items = final_combos[0], axis = 0)
    #converting everything to a True False matrix
    #The [* 1] makes all of the False values a 0 and the True values a 1
    fullDF = fullDF.notna() * 1
    ####
    # subsetting the dataframe into the groups
    ####
    clanDFs = []
    for combo in final_combos:
        clanDF = fullDF.filter(items = combo, axis = 0)
        clanDFs.append(clanDF)
    #Getting the mean values of each column, and converting them to True for pass or False for fail based on the cuttoff
    for i in range(len(clanDFs)):
        clanDF = clanDFs[i]
        clanDF.loc['mean'] = clanDF.mean()

    # fullDF.loc['PorF'] = fullDF.loc['mean'] >= cuttoff
        pORf = clanDF.loc['mean'] >= cuttoff
        numFalse = 0
        numTrue = 0
        try:
            numFalse = pORf.value_counts().loc[False]
        except KeyError:
            numFalse = 0
        try:
            numTrue = pORf.value_counts().loc[True]
        except KeyError:
            numTrue = 0
        #this clanDF passed! 🤩
        if numTrue >= numGenesThreshold:
            resultsFromEachCombo.append(True)

            if returnClan_alreadyEdited == False:
                #adding in the mean data
                #you need the .T to make it add it as a row and not a column!
                returnThisClanDF = pd.concat([returnThisClanDF, pORf.to_frame().T], ignore_index=False)
                #this next line works becuase the new series is alwasy added to the bottom
                #editing the mean index string
                as_list = returnThisClanDF.index.tolist()
                idx = as_list.index('mean')
                as_list[idx] = 'Passed'
                returnThisClanDF.index = as_list
                returnClan_alreadyEdited = True

        #This clanDF failed 😩
        else:
            resultsFromEachCombo.append(False)
            if returnClan_alreadyEdited == False:
                returnThisClanDF = pd.concat([returnThisClanDF, pORf.to_frame().T], ignore_index=False)
                as_list = returnThisClanDF.index.tolist()
                idx = as_list.index('mean')
                as_list[idx] = 'Passed'
                returnThisClanDF.index = as_list
                returnClan_alreadyEdited = True
    #making sure that all of the result combos came to the same conclusion
    theSame = all(element == resultsFromEachCombo[0] for element in resultsFromEachCombo)
    if theSame:
        pORf = resultsFromEachCombo[0]
        #copy over the files and write out one of the clanDFs
        return [pORf, returnThisClanDF]
    else:
        #the conclusions are not agreeing, look into this group by hand to determine why
        print("##################\n" + clan + " had some combos that passed the filtering and some that failed. Look into this group by hand.\n##################")
        return [False, returnThisClanDF]

#simple funciton that checks if something is already at a given spot in the dataframe, empty spots much be populated with np.nan for this to work
def DoesElementExistHere(fullDF, columnHeader, BGC):
    #checks if spot is null, resets the index, and grabs the true of false value only
    checkIfNull = pd.isnull(fullDF.loc[fullDF.BGC_names == BGC, columnHeader]).reset_index(drop=True)[0]
    if checkIfNull:
        return False
    else:
        return True
#populated the dataframes based on the groups from the json file and the mappings!
def PopulateGroups(map, groups, fullDF):
    bgcsInDF = fullDF['BGC_names'].tolist()
    for group in groups:
        groupName = "_".join(group['label'].split())
        #adding this group as a column on the dataframe
        fullDF[groupName] = np.nan
        #going through the genes in the group and adding them to the DF
        genesInGroup = group['genes']
        for uid in genesInGroup:
            mappedUID = map[uid]
            protein_id = mappedUID[0]
            bgcName = "_".join(mappedUID[-1].split())
            #making sure it wasn't one of the outliers that got removed
            if bgcName in bgcsInDF:
                check = DoesElementExistHere(fullDF, groupName, bgcName)
                #if there isn't anything at this location add it this way:
                if not check:
                    fullDF.loc[fullDF.BGC_names == bgcName, groupName] = protein_id
                else: #append it to the existing spot
                    fullDF.loc[fullDF.BGC_names == bgcName, groupName] = fullDF.loc[fullDF.BGC_names == bgcName, groupName].to_string(index=False).strip() + "," + protein_id
    return fullDF


def MapTheJson(data):
    map = {}
    clusters = data['clusters']
    for cluster in clusters:
        bgcName = clusters[cluster]['name']
        #there should only be one element in this list
        loci = clusters[cluster]['loci'][0]
        geneList = data['loci'][loci]['genes']
        for gene in geneList:
            gene_uid = gene
            protein_id = data['genes'][gene_uid]['label']
            map[gene_uid] = [protein_id, bgcName]
    return map

# def MapTheJsonFromHTML(clusters):
#     map = {}
#     for cluster in clusters:
#         bgcName = clusters[cluster]['name']
#         loci = cluster['loci']
#         for group in loci:
#             genes = group['genes']
#             for gene in genes:
#                 gene_uid = gene['uid']
#                 protein_id = gene['label']
#                 map[gene_uid] = [protein_id, bgcName]
#     return map


#adds the no duplicate list to every combo
def AddUnfiltered(all_combos, noDuplicates_BGCS):
    final_combos = []
    for combo in all_combos:
        final_combos.append(combo + noDuplicates_BGCS)
    return final_combos

#permutations of the duplicate species
def GetCombinations(allDuplicates_BGCS):
    num_genes = len(allDuplicates_BGCS)
    # np.prod([len(gene) for gene in allDuplicates_BGCS])
    num_coords = reduce(mul, [len(gene) for gene in allDuplicates_BGCS], 1)
    all_combos = []
    #there are less than 10^4 pairwise combinations
    if num_coords < 10000:
        for coord_num in range(num_coords):
            this_coord = []
            for i in range(num_genes):
                this_len = len(allDuplicates_BGCS[i])
                this_coord.append(coord_num % this_len)
                coord_num //= this_len
                # all_coords.append(this_coord)
            all_combos.append([allDuplicates_BGCS[i][coord] for i, coord in enumerate(this_coord)])
    #if this situation arrises, grab 100 random samplings and return that instead
    else:
        #getting random elements from a list
        counter = 1
        while counter <= 100:
            this_sampling = []
            for i in range(num_genes):
                sppList = allDuplicates_BGCS[i]
                spp = random.sample(sppList, 1)[0]
                this_sampling.append(spp)
            all_combos.append(this_sampling)
            counter += 1
    return all_combos

#subfunction that returns bgcs based on a species
def FindTheBGC_FromSpecies(bgcList, species, duplicates):
    if duplicates == False:
        for bgc in bgcList:
            speciesName = "_".join(bgc.split("_")[1:])
            if speciesName == species:
                return bgc # returns the one bgc
    elif duplicates == True:
        sppList = []
        for bgc in bgcList:
            speciesName = "_".join(bgc.split("_")[1:])
            if speciesName == species:
                sppList.append(bgc)
        #Assuming the list is >3, the program will grab three random species and return that sublist
        if len(sppList) <= 3:
            return sppList # returns a list of the species duplicate bgcs
        else:
            randomSubset = random.sample(sppList, 3)
            return randomSubset #random sampling of three from the duplicate species
    #something went wrong if this condition is fulfilled
    else:
        print('There was an error with this clan. Double check the output.')

#subfunction that creates the dataframe(s) with pairwise combinations for duplicate species
def CreateTheDataFrames(clanPath):
    unfilteredListOfBGCS = []
    unfilteredListOfSpecies = []
    for gbks in os.listdir(clanPath):
        if os.path.isfile(os.path.join(clanPath,gbks)) and gbks.endswith(".gbk") and not gbks.startswith("."):
            bgcName = gbks.split(".gb")[0]
            speciesName = "_".join(bgcName.split("_")[1:])
            unfilteredListOfBGCS.append(bgcName)
            unfilteredListOfSpecies.append(speciesName)
    #making the full dataframe
    colNames = ['BGC_names']
    xtra = {'BGC_names': unfilteredListOfBGCS}
    fullDF = pd.DataFrame(columns = colNames)
    fullDF = fullDF.append(pd.DataFrame(xtra))

    #now that the unfiltered lists are created I need to make the sublists
    noDuplicates_BGCS = []
    #any elements in this list have duplicates
    duplicateSpecies = [key for key in Counter(unfilteredListOfSpecies).keys() if Counter(unfilteredListOfSpecies)[key]>1]
    #Step 1: first I want to populate the no duplicateSpecies list
    for spp in unfilteredListOfSpecies:
        if spp not in duplicateSpecies:
            spp_BGC = FindTheBGC_FromSpecies(unfilteredListOfBGCS, spp, False)
            noDuplicates_BGCS.append(spp_BGC)
    #Step 2: making the duplicate species lists if there are duplicates
    #list of lists, of the duplicate species
    allDuplicates_BGCS = []
    if len(duplicateSpecies) > 0:
        for spp in duplicateSpecies:
            duplicateList = FindTheBGC_FromSpecies(unfilteredListOfBGCS, spp, True)
            #adding it to the larger list of lists
            allDuplicates_BGCS.append(duplicateList)
    #Step 3: Making permutations of duplicate scecies samplings from step 2 IF there were duplicate species
    final_combos = []
    if len(allDuplicates_BGCS) > 0:
        all_combos = GetCombinations(allDuplicates_BGCS)
        final_combos = AddUnfiltered(all_combos, noDuplicates_BGCS)
    else:
        final_combos = [noDuplicates_BGCS]
    return [final_combos, fullDF]

############
# MAIN
############
if __name__ == "__main__":
    outputDirectory = argv[1]
    #making the outputDirectory if it doesn't already exist
    if not os.path.isdir(outputDirectory):
        os.mkdir(outputDirectory)
    clanPath = argv[2]
    clan = argv[3]
    ####
    # Setting the threshold for filtering: what percent of the genomes in the dataframes must share 3 or more genes
    # NumGenes is the number of genes that must hit that threshold for a pass to be counted
    ###
    AnalyzeClan(clan, clanPath, outputDirectory)
```

**SimilarityAnalysis.sh**

```sh
#!/bin/sh
#
# have job exit if any command returns with non-zero exit status (aka failure)
cp /staging/gnickles/default.tar.gz ./
#
# replace env-name on the right hand side of this line with the name of your conda environment
ENVNAME=default
# if you need the environment directory to be named something other than the environment name, change this line
ENVDIR=$ENVNAME
#
# these lines handle setting up the environment; you shouldn't have to modify them
export PATH
mkdir $ENVDIR
tar -xzf $ENVNAME.tar.gz -C $ENVDIR
. $ENVDIR/bin/activate
#
#
#opening the files
tar -xzvf "$1".tar.gz
outputDir=$(pwd)
#back to the base dir. Now I'll grab the path of the clan files
cd "$1"
clanPath=$(pwd)
cd ..
python3 SimilarityAnalysis.py "$outputDir" "$clanPath" "$1"
```

**SimilarityAnalysis.sub**

```sh
# Specify the HTCondor Universe (vanilla is the default and is used
universe = vanilla
#
log = $(name)_0.5similarity.log
error = $(name)_0.5similarity.err
#
# Specify your executable (single binary or a script that runs several
executable = SimilarityAnalysis.sh
arguments = $(name)
output = $(name)_0.5similarity.out
#
#make sure it has gluster
requirements = (Target.HasCHTCStaging == true) && (OpSysMajorVer == 7)
#
# Specify that HTCondor should transfer files to and from the
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
#Transfering the files
transfer_input_files = /home/gnickles/SyntenyAnalysis/$(name).tar.gz, /home/gnickles/SyntenyAnalysis/SimilarityAnalysis.py
#transfer_output_files = $(name)_similarityResults.tsv
#
# Tell HTCondor what amount of compute resources
request_cpus = 1
request_memory = 10GB
request_disk = 10GB
#
#
# Getting the instances of the job
queue name from /home/gnickles/SyntenyAnalysis/groups.txt
```

After this, the results were downloaded back on my computer. Then the final step is to conduct a protein domain analysis, and combine it with the gene cluster conservation analysis that was done by clinker. The core file that is needed from the clinker runs is the json file saved when the -s extension is used. 

**CreateSummaryReport.py**

```python
#! /usr/bin/env python
#Author: Grant Nickles

import os
from sys import argv
import pandas as pd
import numpy as np #used in the GetCombinations function
import pdb
from array import *
import re

#the main function that calls all of the subfunctions to generate the summary report document
def CreateSummaryReport(fullTable, groupName, groupTable, outputFilePath):
    ## making the top part of the summary document
    with open(outputFilePath, "w") as result:
        result.write("################\nAnalysis of " + groupName + "\n################")

    ## the first step is to generate the conserved protein domains section
    fullTable['Accession_pfam'] = fullTable.Domain_Accession.str.cat(fullTable.Domain, sep="_")
    fullTable['GCF'] = fullTable['GCF'].astype('str')
    fullTable['GCC'] = fullTable['GCC'].astype('str')
    GetConservedDomains(fullTable, groupName, outputFilePath)

    # the second step is to write out the conserved gene table
    GetConservedGenes(groupTable, fullTable)


#the function that will take the information from the clinker run to get the conserved groups
def GetConservedGenes(groupTable, fullTable):
    with open(outputFilePath, "a+") as result:
        result.write("################\nConserved groups of genes and their domains (to see specific protein names in a given group find that group in the fullDF file):\n")
    # groupTable.replace(np.nan, 0, inplace=True)
    groupNames = groupTable.columns
    numClusters = len(groupTable)
    for group in groupNames:
        if group.startswith("Group"):
            groupColumn = groupTable[group].dropna()
            numInGroup = len(groupColumn)
            percent = numInGroup / numClusters
            if percent >= 0.2:
                fixedPercent = int(round(percent, 2) * 100)
                #getting the protein domains associated with this group
                keepTrying = True
                counter = 0
                exProtein = ""
                while keepTrying == True:
                    exProtein = groupColumn.iloc[counter]
                    if len(exProtein.split(",")) > 1:
                        counter += 1
                        continue
                    else:
                        keepTrying = False
                if keepTrying == True:
                    pdb.set_trace()
                proteinDF = fullTable.loc[fullTable.Protein == exProtein]
                proteinPfams = proteinDF[['Accession_pfam','E-value_hmmer']]
                filteredPfams = proteinPfams[proteinPfams['E-value_hmmer'] <= (1 * pow(10,-10))]
                filteredPfams['Accession_pfam'] = filteredPfams['Accession_pfam'].astype('str')
                pfamList = ", ".join(filteredPfams['Accession_pfam'].tolist())
                with open(outputFilePath, "a+") as result:
                    result.write(group + "      " + str(fixedPercent) + "% conservation      " + pfamList + "\n")




#the function that writes the first part of the outputFile; assumes the groupName will either be a single number (GCF) or two numbers seperated by an underscore (GCC)
def GetConservedDomains(fullTable, groupName, outputFilePath):
    #GCF
    filteringColumn = ""
    if len(groupName.split("_")) == 1:
        filteringColumn = "GCF"
    elif len(groupName.split("_")) == 2:
        filteringColumn = "GCF_GCC"
    else:
        print("There was an error with the groupName!!")
        pdb.set_trace()
    #subsetting the dataframe
    groupDF = fullTable[fullTable[filteringColumn] == groupName]
    pfamArray = groupDF.Accession_pfam.unique()
    clustersList = groupDF.Cluster_Name.unique()
    ## Printing the clusters and the unique species
    numGenomes = len(clustersList)
    groupDF['RealSpecies'] = groupDF['Cluster_Name'].str.replace(r'^[0-9]+_', '', case = False, regex=True)
    uniqSpecies = groupDF.RealSpecies.unique()
    numSpecies = len(uniqSpecies)
    with open(outputFilePath, "a+") as result:
        result.write("\nNumber of genomes in group: " + str(numGenomes) + "\nNumber of species in group: " + str(numSpecies) + "\n\n")
        result.write("Species Names: " + ", ".join(uniqSpecies.tolist()) + "\n")
        result.write("Cluster Names: " + ", ".join(clustersList.tolist()) + "\n\n")

    #making the new dataframe
    col_names = np.insert(pfamArray, 0, "Cluster_Name")
    col_names = np.insert(col_names, 0, "Species")
    domainsDF = pd.DataFrame(columns = col_names)
    for cluster in clustersList:
        clusterDF = groupDF[groupDF['Cluster_Name'] ==  cluster]
        species = re.sub(r'^[0-9]+_', "", cluster)
        addThis = {"Cluster_Name": cluster, "Species": species}
        domainsDF = domainsDF.append(addThis, ignore_index=True)
        # getting the list pfams in this cluster
        clusterPfams = clusterDF.Accession_pfam.unique()
        #updating the dataframe for each cluster
        for pfam in clusterPfams:
            domainsDF.loc[domainsDF.Cluster_Name == cluster, pfam] = 1
    ## Replacing the nan values with zeros so we can calculate the mean values
    domainsDF.replace(np.nan, 0, inplace=True)
    ##Get the means for each species
    domainsDF = MeanValuesSpecies(domainsDF, uniqSpecies.tolist())
    ## Getting the mean values
    # pfamsOnly = domainsDF.drop('Cluster_Name', 1)
    pfamMeans = domainsDF.mean().sort_values(ascending=False)
    ## Keeping only the rows that have more than 20% conservation
    pfams20 = pfamMeans[pfamMeans >= 0.2]
    ##Printing out these pfams
    with open(outputFilePath, "a+") as result:
        result.write("Protein domains that were present in 20% of the clusters (duplicate species are averaged):\n")
        #looping over the series
        for item in pfams20.iteritems():
            domain = item[0]
            percent = str(item[1] * 100)
            if str(domain).lower() == "nan":
                continue
            else:
                outString = domain + ": " + percent + r'%'
                result.write(outString + "\n")


def MeanValuesSpecies(domainsDF, uniqSpecies):
    editedDF = pd.DataFrame(columns = domainsDF.columns)
    for sp in uniqSpecies:
        speciesSubset = domainsDF[domainsDF['Species'] == sp]
        #there were no duplicates, add it to the edited dataframe as is
        if len(speciesSubset) == 1:
            editedDF = editedDF.append(speciesSubset, ignore_index = True)
        else:
            meansOfSpecies = speciesSubset.mean().to_frame().transpose()
            meansOfSpecies['Species'] = sp
            meansOfSpecies['Cluster_Name'] = "Average_For_Species"
            editedDF = pd.concat([editedDF, meansOfSpecies], axis=0)
    return editedDF

if __name__ == '__main__':
    pd.set_option('display.max_rows', None)

    # groupFolder = argv[1]
    groupFullDF = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ConservationAnalysis/1_Analysis/1_2_0.5Similarity_PASSED.tsv'
    groupTable = pd.read_csv(groupFullDF, sep='\t')
    groupName = "1_2"
    #put the direction to either the filtered family file or the filtered clan file
    tablePath = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ConservationAnalysis/FilteredFamilies.tsv'
    fullTable = pd.read_csv(tablePath, sep='\t')
    #outputFile path: the program will overwrite this file if it already exists!
    outputFilePath = r'/Volumes/T7/ICSProject/ConsolidatedResults_3300Genomes/ConservationAnalysis/1_Analysis/1_2_GCFAnalysis.txt'
    CreateSummaryReport(fullTable, groupName, groupTable, outputFilePath)
```

This output file looks like the following:

```
################
Analysis of 1
################
Number of genomes in group: 47
Number of species in group: 17

Species Names: Talaromyces_marneffei, Talaromyces_stipitatus, Talaromyces_verruculosus, Penicillium_occitanis, Talaromyces_pinophilus, Penicillium_roqueforti, Penicillium_expansum, Penicillium_italicum, Penicillium_griseofulvum,
...
...
...
Protein domains that were present in 20% of the clusters (duplicate species are averaged):
PF05141.13_DIT1_PvcA: 100.0%
PF01965.25_DJ-1_PfpI: 80.88235294117648%
PF07690.17_MFS_1: 75.16339869281046%
PF10558.10_MTP18: 35.78431372549019%
PF01042.22_Ribonuc_L-PSP: 34.80392156862745%
PF11022.9_ATP19: 31.372549019607842%
PF01648.21_ACPS: 29.411764705882355%
PF08282.13_Hydrolase_3: 25.816993464052292%
...
...
...
################
Conserved groups of genes and their domains (to see specific protein names in a given group find that group in the fullDF file):
Group_0      22% conservation      PF00004.30_AAA, PF13191.7_AAA_16, PF07724.15_AAA_2, PF13671.7_AAA_33, PF07728.15_AAA_5, PF17862.2_AAA_lid_3, PF02359.19_CDC48_N, PF05496.13_RuvB_N, PF00004.30_AAA, PF13191.7_AAA_16, PF07724.15_AAA_2, PF13671.7_AAA_33, PF07728.15_AAA_5, PF17862.2_AAA_lid_3, PF02359.19_CDC48_N, PF05496.13_RuvB_N
Group_1      33% conservation      PF00122.21_E1-E2_ATPase, PF00403.27_HMA, PF00702.27_Hydrolase
```



## Expansion and Contraction analysis and scripts

This script could be modified for use on other projects. The only two inputs are 

1. the rooted tree to loop over 
2. A two column table, one column with the leaf names and the second with the copy number of what you're studying

You also have to tell the tree where to start iterating over the tree from. This can be done by passing two distant species, and the program will find the MRCA between those leaf tips. It will create the visualization automattically for you marking where significant expansions and contractions have occured.

**IdentifyExpansions_Contractions.py**

```python
#!/usr/bin/env python
#Author: Grant Nickles (gnickles@wisc.edu)
#Usage: This script will systematically locate expansions and contractions in the tree

#Loading the packages
import pandas as pd
#http://etetoolkit.org/download/
from ete3 import Tree, TreeStyle, TextFace, NodeStyle, faces, Face, AttrFace
import pdb
import os
import researchpy as rp # conda install -c researchpy researchpy
import scipy.stats as stats
import numpy as np
import warnings
import colorsys

## Adding in the heatmap to the leafs of the tree
def AddHeatMapFace(BGC_num, node, color):
    heatmapFace = Face()
    heatmapFace.margin_left=200
    heatmapFace.margin_right=200
    heatmapFace.margin_top=100
    heatmapFace.margin_bottom=100
    heatmapFace.rotable=True
    heatmapFace.background.color = color
    node.add_face(heatmapFace, column=0, position="aligned")

    nameFace = AttrFace("name", fsize=50)
    node.add_face(nameFace, column=0, position='branch-top')
    return

#The main function that looks over an ete3 inputed tree to find locations with statistically significant expansions
def LoopOverTree(tree, startingNode, icsBGC_DF, significanceCuttoff, sizeOfDif, heatmapColors):
    statisticallyHigherNodes = []
    #Step 1: looping over the tree
    for node in startingNode.traverse("preorder"):
        if node.is_leaf():
            BGC_num = icsBGC_DF.loc[icsBGC_DF['ID'] == node.name, "NumICSBGCs"].iloc[0]
            node.name = node.name + "_" + str(BGC_num)
            color = heatmapColors[BGC_num]
            AddHeatMapFace(BGC_num, node, color)
        else:
            children = node.get_children()
            #making sure the children are not leafs, if they are, move on!
            if children[0].is_leaf() or children[1].is_leaf():
                continue
            leftLineage = []
            leftNode = ""
            rightLineage = []
            rightNode = ""
            counter = 0
            for c in children:
                if counter == 0:
                    leftLineage = c.get_leaf_names()
                    leftNode = c
                    counter += 1
                else:
                    rightLineage = c.get_leaf_names()
                    rightNode = c

            #Step 2: Adding categorical label to ics dataframe with this information
            allLeafs = leftLineage + rightLineage
            lineageICSs = icsBGC_DF.loc[icsBGC_DF.ID.isin(allLeafs)]
            if len(lineageICSs) != len(allLeafs):
                dfIDs = lineageICSs['ID'].unique().tolist()
                notShared = set(allLeafs) ^ set(dfIDs)
                pdb.set_trace()

            lineageICSs.loc[lineageICSs['ID'].isin(leftLineage), "Lineage"] = "Left"
            lineageICSs.loc[lineageICSs['ID'].isin(rightLineage), "Lineage"] = "Right"
            #making sure there are more than one value in each lineage

            #Step 3: Running the t-test to see if there is a statistically sig dif. between the mean ICS BGC values of the two sister lineages
            try:
                summary, results = rp.ttest(group1 = lineageICSs['NumICSBGCs'][lineageICSs['Lineage'] == "Right"], group1_name = "Right_Lineage", group2 = lineageICSs['NumICSBGCs'][lineageICSs['Lineage'] == "Left"], group2_name= "Left_Lineage")
            except:
                continue
            #Step 4: Some filtering step, I want the difference to be larger than 2 at a MIN. And the pvalue needs to be lower than 0.05
            dif = results.iloc[0,1]
            pValue = results.iloc[3,1]

            if dif > 1 and pValue <= 0.05:
            #Step 5: Adding the larger lineage to the list
                print("The lineage on the left: " + ",".join(leftLineage))
                print("The lineage on the right: " + ",".join(rightLineage))
                print(summary)
                if dif >= 0:
                    ModifyTreeNodes_Significance(rightNode, leftNode)
                elif dif < 0:
                    ModifyTreeNodes_Significance(leftNode, rightNode)

#The helper function that edits the node IDs on the Tree
def ModifyTreeNodes_Significance(expansionNode, contractionNode):
    expansionNode.name = expansionNode.name + "**"
    contractionNode.name = contractionNode.name + "*"

    #adding in the node styles
    nStyle_exp = NodeStyle()
    nStyle_exp["fgcolor"] = "DarkGreen"
    nStyle_exp["hz_line_color"] = "Green"
    nStyle_exp["vt_line_color"] = "Green"
    nStyle_exp["size"] = 200
    nStyle_exp["hz_line_width"] = 125
    nStyle_exp["vt_line_width"] = 125
    nStyle_exp["bgcolor"] = 'Honeydew'

    nStyle_con = NodeStyle()
    nStyle_con["fgcolor"] = "FireBrick"
    nStyle_con["hz_line_color"] = "Red"
    nStyle_con["vt_line_color"] = "Red"
    nStyle_con["size"] = 200
    nStyle_con["hz_line_width"] = 125
    nStyle_con["vt_line_width"] = 125
    nStyle_con["bgcolor"] = "MistyRose"

    expansionNode.set_style(nStyle_exp)
    contractionNode.set_style(nStyle_con)

    return

#creates the visualization of the tree expansions and contractions
def MakeVisualization(tree):
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = True
    t.show()

def GenerateHeatMapDictionary():
    print("TODO")


#Where the code is run
if __name__ == '__main__':
    ######
    # USER DEFINED VARIABLES
    ######
    #The tree file for the program to use, MUST BE OUTGROUPED PROPERLY FOR THIS TO WORK
    treeFile = r'/Users/gnickles/Desktop/ICS_Fungi_Project/DataAnalysis/ExpansionContractions/TTest_Input/AscomyceteTree_Subset.tre'
    #num icsBGCs, two columns, first is called ID with the tree leafs, second column is number of ICS BGCs (headers included)
    numICSBGCs_file = r'./TTest_Input/ICSBGCS_AstralSpeciesTree.tsv'
    #specific to the ete3 formats for reading and writing newwick files (see documentation)
    fmt = 0
    #defining where the program should start parsing through the tree
    RootOrNode = "Node"
    findBetween = ['Aspergillus_fumigatus', 'Saccharomyces_cerevisiae'] #leave blank if using root
    # pvalue must be <= this value
    significanceCuttoff = 0.05
    #min size difference between the group means to display, lower will include more significant hits, higher will return only the best expansion candidates
    sizeOfDif = 1
    saveTreeAs = r'/Users/gnickles/Desktop/ICS_Fungi_Project/DataAnalysis/ExpansionContractions/TTest_Input/AscomyceteTree_EditedTree.tre'
    #Color dictionary
    # heatmapColors = {0:"white",1:"#ff0404",2:"#ed6a00",3:"#d79800",4:"#b8bd00",5:"#8adf07",6:"#00ff15"}
    heatmapColors = {0:"white",1:"Crimson",2:"Tomato",3:"Yellow",4:"Orange",5:"YellowGreen",6:"Green"}
    ######
    # READING IN THE FILES
    ######
    tree = Tree(treeFile, format=fmt)
    if RootOrNode == "Node":
        node = tree.get_common_ancestor(findBetween[0], findBetween[1])
    elif RootOrNode == "Root":
        node = tree.get_tree_root()
    icsBGC_DF = pd.read_csv(numICSBGCs_file, sep='\t')

    for n in tree.traverse():
        if n.support < 1:
            pdb.set_trace()
        nstyle = NodeStyle()
        nstyle["fgcolor"] = "black"
        nstyle["size"] = 1
        nstyle["vt_line_width"] = 25
        nstyle["hz_line_width"] = 25
        n.set_style(nstyle)

    ######
    # Step 1. Locating the statistically significant nodes on the tree using T-Tests
    ######
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        statisticallyHigherNodes = LoopOverTree(tree, node, icsBGC_DF, significanceCuttoff, sizeOfDif, heatmapColors)

    ######
    # Step 2. Visualizing the tree
    ######
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = True
    ts.title.add_face(TextFace("Expansions and Contractions", fsize=20), column=0)
    ts.mode = "c"
    # ts.branch_vertical_margin = 3
    ts.draw_guiding_lines = True
    # ts.root_opening_factor = 0.1
    ts.force_topology
    ts.allow_face_overlap = True
    #saving and displaying the results
    # tree.render(r'/Users/gnickles/Desktop/ICS_Fungi_Project/DataAnalysis/ExpansionContractions/TTest_Output/AscomyceteExpansions_Contractions.png', w=300, units="mm", tree_style=ts, dpi=700)
    cmd = r"open /Users/gnickles/Desktop/ICS_Fungi_Project/DataAnalysis/ExpansionContractions/TTest_Output/AscomyceteExpansions_Contractions.png"
    os.system(cmd)
    tree.show(tree_style=ts)
    #tree.write(format = 0, outfile=saveTreeAs)
    #tree.render(r'/Users/gnickles/Desktop/ICS_Fungi_Project/DataAnalysis/ExpansionContractions/TTest_Output/TTest_ExpContractions.svg', w=300, units="mm", tree_style=ts)


```

