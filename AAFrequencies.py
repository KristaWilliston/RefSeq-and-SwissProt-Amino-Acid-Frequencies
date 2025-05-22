#0.  Read through and try the examples from Chapters 2-5 of BioPython's Tutorial.

#1a. Download human proteins from RefSeq and compute amino-acid frequencies for the (RefSeq) human proteome.
#       Which amino-acid occurs the most? The least?
#       Hint: access RefSeq human proteins in human.protein.fasta.gz from the course data-folder.

#1b. Download human proteins from SwissProt and compute amino-acid frequencies for the SwissProt human proteome.
#       Which amino-acid occurs the most? The least?
#       Hint: access UniProt XML format SwissProt human proteins from 	http://www.uniprot.org/downloads -> “Taxonomic divisions”

#1c. How similar are the human amino-acid frequencies of in RefSeq and SwissProt? 
#       Which amino-acids show the biggest difference in frequency?

import Bio.SeqIO
import sys
import gzip

# Check the input
if len(sys.argv)<3:
    print >>sys.stderr, "Please provide a fasta file and an XML file"
    sys.exit(1)

# Get the sequence filename
#human.protein.fasta.gz
fasta_filename = sys.argv[1]
#uniprot_sprot_human.xml.gz
xml_filename = sys.argv[2]

# Open the FASTA file and iterate through its sequences
fasta_counts = {}
fasta_total = 0
fasta_file = gzip.open(fasta_filename, mode="rt")
for seq_record in Bio.SeqIO.parse(fasta_file, "fasta"):
    for aa in seq_record.seq:
        #1a. Compute Frequencies
        if aa not in fasta_counts:
            fasta_counts[aa] = 0
        fasta_counts[aa] += 1     #add 1 to the FASTA counts
        fasta_total += 1          #add1 to the FASTA total

#FASTA Frequencies
print("FASTA File Amino Acid Frequencies:")
for k, v in sorted(fasta_counts.items(), key=lambda p: p[1], reverse = True):
    fasta_percentage = (v / fasta_total) * 100  #change values from numbers to percentages
    print(round(fasta_percentage, 3), "%", k)   #round percentages to 3 decimal places
print()

# Open the XML file and iterate through its sequences
xml_counts = {}
xml_total = 0
xml_file = gzip.open(xml_filename, mode="rt")
for seq_record in Bio.SeqIO.parse(xml_file, "uniprot-xml"):
    for aa in seq_record.seq:
        #1b. Compute Frequencies
        if aa not in xml_counts:
            xml_counts[aa] = 0
        xml_counts[aa] += 1         #add 1 to the XML counts
        xml_total += 1              #add 1 to the XML total

#XML Frequencies
print("XML File Amino Acid Frequencies:")
for k, v in sorted(xml_counts.items(), key = lambda p: p[1], reverse = True):
    xml_percentage = (v / xml_total) * 100      #change values from numbers to percentages
    print(round(xml_percentage, 3), "%", k)     #round percentages to 3 decimal places
print()


#1c. Difference In Frequency
differences = {}                    #to store abs difference in percentages

#check all amino acids in FASTA file
for aa, fasta_count in fasta_counts.items():
    fasta_percentage = (fasta_count / fasta_total) * 100                                #calcs FASTA percentage of amino acids
    xml_percentage = xml_counts.get(aa,0) / xml_total * 100 if xml_total > 0 else 0     #gets count of amino acids from XML dictionary & calcs percentage avoiding 0's
    difference = abs(fasta_percentage - xml_percentage)                                 #calcs abs difference between FASTA and XML percentages
    differences[aa] = round(difference, 3)                                              #stores value in differences dictionary
                                    
#check all amino acids in XML file
for aa in xml_counts:
    if aa not in differences:                                                           #if amino acid isnt already in differences dictionary
        xml_percentage = (xml_counts[aa] / xml_total) * 100 if xml_total > 0 else 0     #calcs frequency of amino acid that isn't already present in the XML file
        fasta_percentage = (fasta_counts.get(aa, 0) / fasta_total) * 100                #calcs abs difference between FASTA and XML percentages
        difference = abs(fasta_percentage - xml_percentage)                             #stores value in differences dictionary
        differences[aa] = round(difference, 3)

#sort results from biggest difference to least difference in percentage
sorted_differences = sorted(differences.items(), key = lambda item: item[1], reverse = True)
print("Amino Acid Differences in Percentages:")
for aa, diff in sorted_differences:
    print(aa, diff, "%")
