import csv
with open('/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME_reefgenomics/01_Coral/Stylophora_pistillata_genome_reefgenomics/Spist_Annotation_with_gene_name_GO.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader: 
        #print(row[1],row[2])
        a=row[2]
        b=a.split(",")
        for GO in b:
            print(row[1],GO) 


