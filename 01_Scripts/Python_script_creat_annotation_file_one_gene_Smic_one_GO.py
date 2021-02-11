import csv
with open('/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME/02_Symbiodinium/01_CLADE_A_Smicroadriaticum/Python_first_script/Smicro_genename_gene_newname_GO.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader: 
        #print(row[1],row[3])
        a=row[3]
        b=a.split(",")
        for GO in b:
            print(row[1],GO) 


