import json, sys, csv, argparse, operator, re

s_res = re.compile ('^(High-level resistance|Intermediate resistance)$')
tag_re = re.compile ('^(.+)\\|[0-9]+$')

parser = argparse.ArgumentParser(
    description='Compare IDEPI and Stanford DRAM classifications'
)

parser.add_argument (
    '-i', '--idepi',
    type=argparse.FileType('r'),
    help='idepi predict JSON output file',
    required = True,
)

parser.add_argument (
    '-s', '--stanford',
    type=argparse.FileType('r'),
    help='read the Stanford algorithm classification file from',
    required = True,
)

parser.add_argument (
    '-d', '--drug',
    help='which antiretroviral drug',
    type=str,
    required = True,
) 


def Cohen_kappa (cfm):
    N = sum ([sum (k) for k in cfm])
    agree  = (cfm[0][0] + cfm[1][1]) / N
    idepi_yes = cfm[1][0] + cfm[1][1]
    stanford_yes = cfm[0][1] + cfm[1][1]
    random = (idepi_yes*stanford_yes + (N-idepi_yes)*(N-stanford_yes))/N**2
    return (agree-random)/(1-random)

args = None
retcode = -1
args = parser.parse_args()

labels = json.load (args.idepi)['predictions']
idepi_ids = {}

confusion_matrix = [[0,0],[0,0]]

for s in labels:
    tag = tag_re.match (s["id"])
    if tag:
        tag = tag.group(1)
    else:
        tag = s["id"]
    if tag in idepi_ids:
        idepi_ids [tag] = max(s["value"],idepi_ids [tag])
    else:
        idepi_ids [tag] = s["value"]

stanford = csv.reader (args.stanford, delimiter = '\t')
column_names = next (stanford)

idx = column_names.index (args.drug)

for s in stanford:
    if len (s):
        row = 1 if idepi_ids[s[0]] > 0 else 0
        column = 1 if s_res.match (s[idx]) else 0
        
        #if row!=column:
        #    print (s[0], idepi_ids[s[0]], s[idx])
        
        confusion_matrix[row][column] += 1
        
N = sum ([sum (r) for r in confusion_matrix])        
print ("N = %d" % N)
print ("IDEPI         Stanford     Count   Percentage",
       "\nSusceptible   Susceptible %6d    %.3g" % (confusion_matrix[0][0], confusion_matrix[0][0]/N*100.),
       "\nResistant     Resistant   %6d    %.3g" % (confusion_matrix[1][1], confusion_matrix[1][1]/N*100.),
       "\nResistant     Susceptible %6d    %.3g" % (confusion_matrix[1][0], confusion_matrix[1][0]/N*100.),
       "\nSusceptible   Resistant   %6d    %.3g" % (confusion_matrix[0][1], confusion_matrix[0][1]/N*100.))
       
        
print ("\nCohen's kappa = %g" % Cohen_kappa(confusion_matrix))
        



