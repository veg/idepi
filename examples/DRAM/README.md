# Drug resistance association mutations (DRAM) in HIV-1

In this example, we will use IDEPI to infer the features in 
HIV-1 reverse transcriptase (RT) associated with nevirapine (NVP)
resistance. All the code examples assume that the current 
working directory is **data/DRAM** (inside the IDEPI
directory).

The data for this analysis are Vicro Antivirogram phenotypic
measurements available at the [Stanford HIV drug resistance
database](http://hivdb.stanford.edu/modules/lookUpFiles/phenoGenoDataSet/RT_Antivirogram_DataSet.txt)
The exact used we ran for the IDEPI manuscript is available in 
**raw/RT.txt**.

## Preliminaries

It is first necessary to convert the raw phenotypic measurements,
reported as a list of amino-acid mutations from the HIV-1 reference
sequence (HXB2), into actual RT nucleotide sequences, and also 
generate a CSV file with the expected format of:

    SeqID1, Drug1_IC50, ..., DrugM_IC50
    ...
    SeqIDN, Drug1_IC50, ..., DrugM_IC50
    
A Python script for performing this conversion is provided
in **helpers/tabtofasta.py**. 

    usage: tabtofasta.py [-h] -r {pr,rt} -p PHENO -f FASTA -l LABELS

    Convert a tab-separated drug resistance phenotype file to a sequence file and
    the CSV labels file

    optional arguments:
      -h, --help            show this help message and exit
      -r {pr,rt}, --reference {pr,rt}
                            aminoacid sequence to fill consensus positions from
      -p PHENO, --pheno PHENO
                            the tab file with empirical drug data are stored
      -f FASTA, --fasta FASTA
                            write the FASTA file to
      -l LABELS, --labels LABELS
                            write the labels file to

For instance, make the following call to generate an (unaligned) FASTA file, and the 
corresponding phenotype labels (fold reduction in susceptibility). 

    $python3 helpers/tabtofasta.py -r rt -p raw/RT.txt -f training/rt.fas -l training/phenotypes.csv
    
Because we are using HIV-1 pol, there is no prebuilt HMMER model in IDEPI to 
generate the necessary multiple sequence alignment (MSA), so we will need to align
the sequences in the **rt.fasta** file. Such an MSA can be found in **training/rt.msa**.

## Top DRAM identification

The first example uses the **idepi discrete** module to perform a cross-validation 
based model performance evaluation, and extract top **N** features. 

### Fit the model with a single feature (**N=1**).

    $idepi discrete --csv training/rt.msa training/phenotypes.csv --refmsa training/rt.msa --refseq training/hxb2_rt.fas --label "max(IC50)>=5" --mcc --numfeats 1 NVP
    
We'll go through the arguments in a second, but first let's look at the output. IDEPI
reports (to stderr) a number of warning messages of the form:
    
    /opt/python-3.3.1/lib/python3.3/site-packages/idepi/datasource/__init__.py:248: UserWarning: skipping sequence '9576', VALUE not found
    warn("skipping sequence '%s', VALUE not found" % r.id)

These messages simply mean that a subset of sequences in the phenotype label file
do not have any data for nevirapine (NVP), and will not be used for model training.
IDEPI will also create a Stockholm (.sto) file with the HMMER alignment for the
sequence data in the current working directory.

The output of the program is written to stdout as a JSON file (which is human readable,
but also machine readable into Python, JavaScript, and almost all other widely used
languages).

    {
      "metadata": {
        "antibodies":  [ "NVP" ],
        "balance":     0.622861,
        "features":    1,
        "folds":       5,
        "label":      "max(IC50)>=5",
        "sequences":   1461
      },
      "statistics": {
        "accuracy":    { "mean": 0.652941, "std": 0.069201 },
        "f1score":     { "mean": 0.621673, "std": 0.098219 },
        "mcc":         { "mean": 0.457632, "std": 0.100586 },
        "npv":         { "mean": 0.522352, "std": 0.053964 },
        "ppv":         { "mean": 0.963492, "std": 0.055948 },
        "sensitivity": { "mean": 0.460440, "std": 0.106938 },
        "specificity": { "mean": 0.970958, "std": 0.043419 }
      },
      "weights": [
        { "position": "K103K", "N": 5, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } }
      ]
    }

1. **_metadata_** reports the information about the settings and data used by IDEPI
    + **antibodies** (NVP) should really be called "phenotypes", and leftover from the 
    the original IDEPI design, focusing only on antibody epitopes
    + **balance** (0.622) describes the proportion of sequences defined as having the
    "resistant" phenotype (defined by the 'label' field)
    + **features** (1) reports the number of features used by IDEPI to build the model
    + **folds** (5) the (default) number of cross-validation folds used by IDEPI to 
    compute model performance 
    + **label** ("max(IC50)>=5") the Python formula used by IDEPI to convert phenotype
    measures (IC50) into a binary phenotype. In this case, any sequence with 5x or greater
    reduction in susceptibility will be labeled resistant.
    + **sequences** (1461) the total number of (sequence, phenotype) pairs used by IDEPI
    to build the model

2. **_statistics_** records all the standard model performance metrics computed with nested
    cross validation. The mean and standard deviation for each metric (in this case over
    5 CV folds) are reported. E.g. the mean Matthew's Correlation Coefficient (mcc) is
    0.46 (the value may differ from run to run because of the stochastic nature of CV).

3. **_weights_** is the combined report of all the features found by the classification 
procedure in any of the CV folds. 
    + **position** the feature selected by the model (in this case the presence of a lysine
    at reference position 103 where the reference sequence also has a lysine).
    + **N** in how many of the cross validation folds was this feature selected (the greater
    the value, the more robust the feature)
    + **rank** the rank of this feature (mean and st.dev) in the models derived by the CV 
    folds (see below for more interesting examples). Here, in all 5 folds, the same feature
    (_K103K_) was selected.
    + **value** does this feature correspond to susceptibility (value -1) or resistance
    (value 1) in each of the CV fold. In this example, the _K103K_ feature was universally 
    associated with susceptibility, i.e. a sequence with anything other than
    a lysine at reference position 103 would be classified resistant.
    
Now we describe all the command line arguments.
    
    $idepi discrete --csv training/rt.msa training/phenotypes.csv --refmsa training/rt.msa --refseq training/hxb2_rt.fas --label "max(IC50)>=5" --mcc --numfeats 1 NVP

-  **--csv training/rt.msa training/phenotypes.csv** read the sequences and matched phenotypes from 
this pair of files
+  **--refmsa training/rt.msa** train the HMMER alignment model from this alignment. For this example,
the sequences are already aligned, but the HMMER model still needs to be built so that we can align
new sequences (and the reference sequence) to it. 
+  **--refseq training/hxb2_rt.fas** use this sequence to define reference coordinates
+  **--label "max(IC50)>=5"** define sequences as resistant if their IC50 change is 5-fold or greater. 
The reason we use "max" is because IC50 is encoded as a vector (of possibly multiple phenotypes, see 
later examples).
+  **--mcc** tune the model to maximize the Matthews correlation coefficient
+  **--numfeats 1** build models with a single feature
+  **NVP** read phenotypic information from the NVP column in the phenotypes file.
	
### Fit the model with five features (**N=5**).

Only the **--numfeats** argument is changed.
    
    $idepi discrete --csv training/rt.msa training/phenotypes.csv --refmsa training/rt.msa --refseq training/hxb2_rt.fas --label "max(IC50)>=5" --mcc --numfeats 5 NVP

IDEPI reports the following model performance and features

    {
      "metadata": {
        "antibodies":  [ "NVP" ],
        "balance":     0.622861,
        "features":    5,
        "folds":       5,
        "label":      "max(IC50)>=5",
        "sequences":   1461
      },
      "statistics": {
        "accuracy":    { "mean": 0.858972, "std": 0.067390 },
        "f1score":     { "mean": 0.874432, "std": 0.067357 },
        "mcc":         { "mean": 0.738751, "std": 0.107121 },
        "npv":         { "mean": 0.743124, "std": 0.106687 },
        "ppv":         { "mean": 0.974937, "std": 0.036442 },
        "sensitivity": { "mean": 0.794505, "std": 0.114886 },
        "specificity": { "mean": 0.965536, "std": 0.050365 }
      },
      "weights": [
        { "position": "L74L",   "N": 4, "rank": { "mean": 4.75, "std": 0.87 }, "value": { "mean": -0.50, "std": 1.73 } },
        { "position": "K103K",  "N": 5, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "Y181Y",  "N": 5, "rank": { "mean": 2.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "Y188Y",  "N": 1, "rank": { "mean": 5.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "G190G",  "N": 5, "rank": { "mean": 3.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "T215T",  "N": 1, "rank": { "mean": 4.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "G231[]", "N": 1, "rank": { "mean": 4.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "E233[]", "N": 1, "rank": { "mean": 5.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "H235[]", "N": 1, "rank": { "mean": 4.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "P236P",  "N": 1, "rank": { "mean": 4.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } }
      ]
    }
    
The key difference from the **N=1** example is how the **_weights_** record is interpreted here.
*K103K* (susceptible) is the consistent top feature (rank 1) among all the CV folds, followed by *Y181Y* (susceptible, rank 2),
and *G190G* (susceptible, rank 3). For the 4th and 5th ranked features there is now noticeable variation among the folds, with 
several features, like *T215T* (resistant), *P236P* (resistant), *G231[]* (resistant, i.e. a deletion/missing data at position 231), etc
appearing as the 4th ranked feature in one CV fold only. This type of behavior suggests that either the model is starting to 
overfit the data, or that we haven't added enough features to encompass the complexity of the DRAM profile.  

### Direct IDEPI to perform a grid search over the number of features

It is possible to supply a range of values to the **--numfeats** option, and ask IDEPI to select the model with the 
best performance according to the selected criterion.

    idepi discrete --csv training/rt.msa training/phenotypes.csv --refmsa training/rt.msa --refseq training/hxb2_rt.fas --label "max(IC50)>=5" --mcc --numfeats 1:5,10:50:5,60:100:10 NVP

The grid is specified using a comma separated list of **from:to:step**, where **step** is an optional argument (defaulting to one).
E.g. the _1:5,10:50:5,60:100:10_ argument specifies the 1,2,3,4,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100 grid.

**Note** this analysis make take a long time to run.

The output will look much like in the previous two examples, except that the feature list will
be quite long (abbreviated here to features ranked between 1 and 10 on average). 
IDEPI selected **N=80** as the number of features yielding the best mean MCC (0.83).  

    {
      "metadata": {
        "antibodies":  [ "NVP" ],
        "balance":     0.622861,
        "features":    80,
        "folds":       5,
        "label":      "max(IC50)>=5",
        "sequences":   1461
      },
      "statistics": {
        "accuracy":    { "mean": 0.915805, "std": 0.030167 },
        "f1score":     { "mean": 0.928642, "std": 0.026165 },
        "mcc":         { "mean": 0.834433, "std": 0.059813 },
        "npv":         { "mean": 0.831624, "std": 0.045891 },
        "ppv":         { "mean": 0.982918, "std": 0.026773 },
        "sensitivity": { "mean": 0.880220, "std": 0.037588 },
        "specificity": { "mean": 0.974595, "std": 0.039407 }
      },
      "weights": [
      ...
        { "position": "L74L",           "N": 5, "rank": { "mean":  5.60, "std":  3.90 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "K103K",          "N": 5, "rank": { "mean":  1.00, "std":  0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "K103N",          "N": 5, "rank": { "mean":  6.20, "std":  0.89 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "Y181Y",          "N": 5, "rank": { "mean":  2.00, "std":  0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "Y188Y",          "N": 5, "rank": { "mean":  6.60, "std":  1.79 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "G190G",          "N": 5, "rank": { "mean":  3.00, "std":  0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
       ...
       ]
    }




## Learning the model of NVP resistance

Unlike **idepi discrete** which uses nested cross-validation to study the robustness of the 
features, **idepi learn** uses the complete data set (although it still uses cross-validation
to tune internal model, e.g. LSVM parameters) to build a predictive model and save it to file
so that sequences of unknown phenotypes can be classified later.

The arguments are exactly the same as those passed to **idepi discrete**, except the final
positional argument is required to specify the name of the file (in this case **models/NVP.model**)
to save the IDEPI model to. We build a model on **N=80** features here, based on the grid 
search from the previous step.

    $idepi learn --csv training/rt.msa training/phenotypes.csv --refmsa training/rt.msa --refseq training/hxb2_rt.fas --label "max(IC50)>=5" --mcc --numfeats 80 NVP models/NVP.model

The output to stdout is exactly the same as for **idepi discrete**, except that there is only 
one "fold". For brevity, only the top 5 ranked features are retained in the listing below.

    {
      "metadata": {
        "antibodies":  [ "NVP" ],
        "balance":     0.622861,
        "features":    80,
        "folds":       1,
        "label":      "max(IC50)>=5",
        "sequences":   1461
      },
      "statistics": {
        "accuracy":    { "mean": 0.924025, "std": 0.000000 },
        "f1score":     { "mean": 0.935875, "std": 0.000000 },
        "mcc":         { "mean": 0.850016, "std": 0.000000 },
        "npv":         { "mean": 0.843750, "std": 0.000000 },
        "ppv":         { "mean": 0.986602, "std": 0.000000 },
        "sensitivity": { "mean": 0.890110, "std": 0.000000 },
        "specificity": { "mean": 0.980036, "std": 0.000000 }
      },
      "weights": [
        ...
        { "position": "L74L",           "N": 1, "rank": { "mean":  5.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "K103K",          "N": 1, "rank": { "mean":  1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "Y181Y",          "N": 1, "rank": { "mean":  2.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "G190G",          "N": 1, "rank": { "mean":  3.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "E233[]",         "N": 1, "rank": { "mean":  4.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } }
        ...
       ]
    }

## Predicting the phenotype of sequences

Given the model constructed with **idepi learn** and a set of sequences with an unknown
phenotype, we will use **idepi predict** to classify the NVP resistance status of each
sequence. The 1655 sequences used for this purpose come from the paper by
[Avila-Rios et al](http://www.ncbi.nlm.nih.gov/pubmed/22110765), included in the 
IDEPI distribution in the **validation/PMID22110765.msa** file.

**Note** If ambiguous nucleotides were found in key DRAM positions in a sequence,
multiple possible resolutions were created (e.g. sequences __M000MA6FTK_2008-06|0__ and 
__M000MA6FTK_2008-06|1__ represent two possible resolutions). A sequence was classified
resistant if **any** of its resolutions were. This is a standard practice in the field.

Execute:

    $idepi predict models/NVP.model validation/PMID22110765.msa > validation/idepi.json
    
This **idepi predict** script expects to be given the model (saved by **idepi learn**) and
a collection of sequences (homologous to those used to train the model) as a FASTA file (
it does not have to be aligned).

The output, written to stdout, and redirected to **validation/idepi.json** in this example,
is a JSON file, where, for each input sequence, the JSON object contains a record of the 
form (in the _predictions_ field)

    { "id": "JC00MC6GRU_2007-11", 
    "value": -1, 
    "features": [ "L74L", "V75V", "A98A", "L100L", "K101K", "V1
    06V", "V108V", "F116F", "I135V", "R172R", "K173K", "V179V", "Y181Y", "Y188Y", "V189V", 
    "G190G", "E203E", "L210L", "T215T", "K219K", "H221H", "K223K", "L228L", "G231G", "Y232Y", "E233E", 
    "H235H", "P236P" ] }
    
+ **id** is the name of the sequence   
+ **value** reports the phenotype classification, defined according the the label **max(IC50)>=5**
used to train the model. A value of "1" means the label is True, a value of "-1" means the value is 
False (susceptible in our case).
+ **features** lists all the model features present in the current sequence; for LSVM classifiers,
a weighted sum of these features defines the final score.

## Comparing to the gold standard

We use the [Stanford HIVdb program](http://sierra2.stanford.edu/sierra/servlet/JSierra?action=sequenceInput)
to infer the resistance phenotype of the sequences in **validation/PMID22110765.msa**
The raw output of the program is stored in **validation/stanford_db.txt**

Use the utility script **helpers/idepi2stanford.py** to generate a 2x2 table of 
classification agreement/disagreement between the two algorithms
    
    $python3 helpers/idepi2stanford.py -i validation/idepi.json -s validation/stanford_db.txt -d NVP
    
Output:

    N = 1639
    IDEPI         Stanford     Count   Percentage 
    Susceptible   Susceptible   1604    97.9 
    Resistant     Resistant       26    1.59 
    Resistant     Susceptible      3    0.183 
    Susceptible   Resistant        6    0.366

    Cohen's kappa = 0.849668