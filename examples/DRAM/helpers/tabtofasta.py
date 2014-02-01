import re, sys, csv, argparse
import Bio.Seq
import Bio.SeqRecord
import Bio.Alphabet
from Bio import SeqIO

nuc_ambigs = {"A" : "A",
              "C" : "C",
              "G" : "G",
              "T" : "T", 
              "R" : "AG",
              "Y" : "CT",
              "S" : "CG",
              "W" : "AT",
              "K" : "GT",
              "M" : "AC",
              "B" : "CGT",
              "D" : "AGT",
              "H" : "ACT",
              "V" : "ACG"}
              
codon_table = {"AAA":"K",
            "AAC":"N",
            "AAG":"K",
            "AAT":"N",
            "ACA":"T",
            "ACC":"T",
            "ACG":"T",
            "ACT":"T",
            "AGA":"R",
            "AGC":"S",
            "AGG":"R",
            "AGT":"S",
            "ATA":"I",
            "ATC":"I",
            "ATG":"M",
            "ATT":"I",
            "CAA":"Q",
            "CAC":"H",
            "CAG":"Q",
            "CAT":"H",
            "CCA":"P",
            "CCC":"P",
            "CCG":"P",
            "CCT":"P",
            "CGA":"R",
            "CGC":"R",
            "CGG":"R",
            "CGT":"R",
            "CTA":"L",
            "CTC":"L",
            "CTG":"L",
            "CTT":"L",
            "GAA":"E",
            "GAC":"D",
            "GAG":"E",
            "GAT":"D",
            "GCA":"A",
            "GCC":"A",
            "GCG":"A",
            "GCT":"A",
            "GGA":"G",
            "GGC":"G",
            "GGG":"G",
            "GGT":"G",
            "GTA":"V",
            "GTC":"V",
            "GTG":"V",
            "GTT":"V",
            "TAA":"X",
            "TAC":"Y",
            "TAG":"X",
            "TAT":"Y",
            "TCA":"S",
            "TCC":"S",
            "TCG":"S",
            "TCT":"S",
            "TGA":"X",
            "TGC":"C",
            "TGG":"W",
            "TGT":"C",
            "TTA":"L",
            "TTC":"F",
            "TTG":"L",
            "TTT":"F"}

translations = {}

for n1 in nuc_ambigs:
    for n2 in nuc_ambigs:
        for n3 in nuc_ambigs:
            prots = set ()
            for l1 in nuc_ambigs[n1]:
                for l2 in nuc_ambigs[n2]:
                    for l3 in nuc_ambigs[n3]:
                        prots.add (codon_table[''.join ([l1,l2,l3])])
                        
            key = ''.join(prots)
            if key not in translations:
                translations[key] = []
            
            translations[key].append (''.join ([n1,n2,n3]))
                          
            



position_number_re = re.compile ('^P([0-9]+)$')
drug_score         = re.compile ('^([^\s]+)\sFold$')

rt = "PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDEDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGLTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNDIQKLVGKLNWASQIYPGIKVRQLCKLLRGTKALTEVIPLTEEAELELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEPFKNLKTGKYARMRGAHTNDVKQLTEAVQKITTESIVIWGKTPKFKLPIQKETWETWWTEYWQATWIPEWEFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRETKLGKAGYVTNRGRQKVVTLTDTTNQKTELQAIYLALQDSGLEVNIVTDSQYALGIIQAQPDQSESELVNQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDKLVSAGIRKVLF"
pr = "PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Convert a tab-separated drug resistance phenotype file to a sequence file and the CSV labels file'
    )
 
    parser.add_argument(
        '-r', '--reference',
        type=str,
	choices = ('pr','rt'),
        help='aminoacid sequence to fill consensus positions from',
        required = True,
    )
    parser.add_argument (
        '-p', '--pheno',
        type=argparse.FileType('r'),
        help='the tab file with empirical drug data are stored',
        required = True,
    )
 
    parser.add_argument (
        '-f', '--fasta',
        type=argparse.FileType('w'),
        help='write the FASTA file to',
        required = True,
    )
    
    parser.add_argument (
        '-l', '--labels',
        type=argparse.FileType('w'),
        help='write the labels file to',
        required = True,
    ) 
 
    args = None
    retcode = -1
    args = parser.parse_args()
    
    resistance_reader = csv.reader (args.pheno, delimiter = '\t')
    col_names = next (resistance_reader)
    
    if args.reference == 'rt':
        args.reference = rt
    elif args.reference == 'pr':
        args.reference = pr
        
    position_mapper = {}
    drug_mapper     = {}
    drug_list       = ['ID']
    
    for idx, col_name in enumerate(col_names):
        m = position_number_re.match (col_name)
        if m:
            position_mapper [idx] = int (m.group(1)) - 1
        else:
            m = drug_score.match (col_name)
            if m:
                drug_list.append (m.group(1))
                drug_mapper [idx] = len (drug_list)-1
                
                
    neuts = {}
    
    label_writer = csv.writer(args.labels)
    label_writer.writerow (drug_list)

    for line in resistance_reader: 
        seq = []
        drugz = []
        drugz.append (line[0])
        
        for k in position_mapper:
            if line[k] == '-':
                seq.append (translations[args.reference[position_mapper[k]]][0])
            elif line[k] != '.' and line [k] != 'X':
                key = ''.join(set (line[k]))
                if key in translations:
                    seq.append (translations[key][0])
            
        for d in drug_mapper:
            drugz.append (line[d] if line[d] != 'NA' else None)
                
        prot_seq = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq (''.join (seq), Bio.Alphabet.IUPAC.protein), id = line[0], name = line[0], description = '')
        SeqIO.write(prot_seq, args.fasta, "fasta")
        label_writer.writerow (drugz)
      
    
    sys.exit(retcode)



