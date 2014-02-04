import json, sys, csv, argparse, operator, re
from math import sqrt

has_HAD = re.compile ('^.+\_HAD$')
patient_ID = re.compile ('^([^\_]+).+$')

parser = argparse.ArgumentParser(
    description='tabulate validation results for the HAD example'
)

parser.add_argument (
    '-j', '--json',
    type=argparse.FileType('r'),
    help='IDEPI prediction file',
    required = True,
)


def report (cfm):
    N = sum ([sum (k) for k in cfm])
    accu  = (cfm[0][0] + cfm[1][1]) / N
    specificity = cfm[1][1] / (cfm [0][1] + cfm[1][1])
    sensitivity = cfm[1][1] / (cfm[1][0] + cfm[1][1])

    return [specificity, sensitivity, accu]
    

args = None
retcode = -1
args = parser.parse_args()

labels = json.load (args.json)['predictions']
idepi_ids = {}

confusion_matrix = [[0,0],[0,0]]

by_individual = {}

for s in labels:
    row = 1 if has_HAD.match (s["id"]) is not None else 0
    column = 1 if s["value"] > 0 else 0
    
    pid = patient_ID.match (s["id"]).group (1)
    
    if pid not in by_individual:
        by_individual [pid] = [row, [0,0]]
        
    by_individual [pid][1][column] += 1


for id, prediction in by_individual.items():
    prop_had = prediction [1][1] / sum (prediction [1])
    print ('--Individual %s (%s). Prediction: %g%% HAD' % (id, 'HAD' if prediction [0] else 'Non-HAD', prop_had * 100.))
    
    row    = prediction[0]
    column = 1 if prop_had >= 0.5 else 0
    
    confusion_matrix[row][column] += 1
    
       
 
N = sum ([sum (r) for r in confusion_matrix])        

print ("\nN = %d\n" % N)

print (  "IDEPI     Measured     Count   Percentage",
       "\nNon-HAD   Non-HAD      %6d    %.3g" % (confusion_matrix[0][0], confusion_matrix[0][0]/N*100.),
       "\nHAD       HAD          %6d    %.3g" % (confusion_matrix[1][1], confusion_matrix[1][1]/N*100.),
       "\nNon-HAD   HAD          %6d    %.3g" % (confusion_matrix[1][0], confusion_matrix[1][0]/N*100.),
       "\nHAD       Non-HAD      %6d    %.3g" % (confusion_matrix[0][1], confusion_matrix[0][1]/N*100.))
       
metrics = report (confusion_matrix)
 
print ("\nSpecificity %.2g" %   metrics[0],
       "\nSensitivity %.2g" %   metrics[1],
       "\nAccuracy    %.2g" %   metrics[2])
       
