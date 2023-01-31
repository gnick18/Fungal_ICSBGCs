# Fungal_ICSBGCs
Reproducible code and data repository for the fungal ICS BGC prediction publication. All of the core code used in the publication in addition to the raw prediction files can be found here. 

The metadata for the genomes in addition to tables that include more information on select BGCs can be found in the **Supplemental Information**. 

![ICS_FungalTreeFigure](./Images/ICS_FungalTreeFigure.png)

*Explanation of folder organization:*

### **ICSBGCs**

This folder stores all of the raw ICS BGC predictions in genbank, gff, and fasta format. The BGC predictions for each genome is organized in a nested folder structure by the the NCBI accession of the genome it was found in.

- The genbank files are compatible with popular natural product programs such as BiG-SCAPE and clinker
- The gff file contains the annotations for the BGC region
- The fasta file contains the raw nucleotide sequence of the loci encapsulating all of the genes predicted to be in the ICS BGC
- NOTE: These are the unedited and unfiltered predictions generated from the ICSBGC prediction pipeline. Information on specific gene cluster familes, or refined BGC predictions can be found in other folders

