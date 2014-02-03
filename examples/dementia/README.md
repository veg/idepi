# Finding envelope features predictive of HIV-associated dementia [HAD] 

This example is concerned with identify the set of genetic features
(also commonly called "signatures") of HIV-associated dementia (HAD),
based on sequences isolated from HAD-positive and HAD-negative
HIV infected individuals. 

The data used for this exercise have been obtained from the
[HIV Brain Sequence Database](http://www.hivbrainseqdb.org/about/), 
and previously analyzed in a 2012 paper by [Holman and Gabuzda](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0049538).

All the code examples assume that the current 
working directory is **data/dementia** (inside the IDEPI
directory). 


## Preliminaries

The sequences extracted from the HIV Brain Sequence Database
and used to train the HAD classifier and extract predictive
features are located in **training/HIVSeqDB.fasta**. The database
also contained phenotype information (HAD status) which was 
used to generate sequence labels for training.

The phenotype information in this classification problem is
binary, namely whether or not the patient had HAD.
To conform to IDEPI phenotype specifications, designed to handle
IC50 data, one needs to create a CSV file connecting each 
sequence ID with its phenotype as follows

    ID,HAD
    DQ976387,50
    ...
    DQ288646,0    
    ...
    
The default sequence labeling algorithm in IDEPI will label 
a sequence "resistant" (HAD positive in this case) if its phenotype is >= 20, hence the value
of 50 for the strain isolated from a HAD+ individual in the above example.


## Identifying the genetic feature most predictive of HAD.

We use the **idepi discrete** module identify the top features in 
the training data associated with the development of HAD.

### Fit the model with a single feature (**N=1**).

    $idepi discrete --csv training/HIVSeqDB.fasta training/HAD.csv --numfeats 1 --mcc HAD 
    
The arguments above tell IDEPI to 

+ read paired sequence/phenotype data from training/HIVSeqDB.fasta and training/HAD.csv, respectively
+ train models with a single feature
+ optimize Matthews Correlation Coefficient (mcc) when tuning LSVM parameters
+ use the _HAD_ column to define sequence phenotypes

> The warning messages written to stderr are there because not 
> every sequence in the .fasta file has an associated phenotype in 
> the .csv file.

Because it is hypothesized here that HAD is largely determined by HIV-1 envelope,
IDEPI's existing HMMER model is used to align training
sequences to the standard HXB2 coordinates. 

The output of the program is written to stdout as a JSON file.

    {
      "metadata": {
        "antibodies":  [ "HAD" ],
        "balance":     0.702671,
        "features":    1,
        "folds":       5,
        "label":      "max(IC50) > 20",
        "sequences":   861
      },
      "statistics": {
        "accuracy":    { "mean": 0.828102, "std": 0.038312 },
        "f1score":     { "mean": 0.884866, "std": 0.036295 },
        "mcc":         { "mean": 0.574884, "std": 0.073256 },
        "npv":         { "mean": 0.836237, "std": 0.213031 },
        "ppv":         { "mean": 0.833701, "std": 0.035967 },
        "sensitivity": { "mean": 0.945455, "std": 0.114873 },
        "specificity": { "mean": 0.550905, "std": 0.154718 }
      },
      "weights": [
        { "position": "T297K", "N": 4, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "R298D", "N": 1, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } }
      ]
    }

    
1. **_metadata_** reports the information about the settings and data used by IDEPI
    + **antibodies** (HAD) should really be called "phenotypes", and leftover from the 
    the original IDEPI design, focusing only on antibody epitopes
    + **balance** (0.703) describes the proportion of sequences defined as having the
    "resistant" phenotype (defined by the 'label' field); here these are the HAD+ 
    sequences
    + **features** (1) reports the number of features used by IDEPI to build the model
    + **folds** (5) the (default) number of cross-validation folds used by IDEPI to 
    compute model performance 
    + **label** ("max(IC50)>=20") the Python formula used by IDEPI to convert phenotype
    measures (IC50) into a binary phenotype. 
    + **sequences** (861) the total number of (sequence, phenotype) pairs used by IDEPI
    to build the model

2. **_statistics_** records all the standard model performance metrics computed with nested
    cross validation. The mean and standard deviation for each metric (in this case over
    5 CV folds) are reported. E.g. the mean Matthew's Correlation Coefficient (mcc) is
    0.57 (the value may differ from run to run because of the stochastic nature of CV).

3. **_weights_** is the combined report of all the features found by the classification 
procedure in any of the CV folds. 
    + **position** the feature selected by the model. In this case there is some ambiguity
    about the identity of the top model feature. 4 of the folds selecte the presence of a 
    lysine in reference position 297 (HXB2 has a threonine in this position), i.e. _T297K_, 
    while the remaining fold identifies a neighboring position _R298D_. 
    + **N** in how many of the cross validation folds was this feature selected (the greater
    the value, the more robust the feature)
    + **rank** the rank of this feature (mean and st.dev) in the models derived by the CV 
    folds (see below for more interesting examples). H
    + **value** does this feature correspond to susceptibility (value -1) or resistance
    (value 1) in each of the CV fold. In this example, the both features were universally 
    associated with being HAD-, i.e. a sequence that has the requisite substitution is
    classified as HAD-.
    

### Perform a grid search over the number of features

Specify a grid of possible values to search over for the optimal number of features. 

    $idepi discrete --csv training/HIVSeqDB.fasta training/HAD.csv --numfeats 1:5,10:50:5,60:100:10 --mcc HAD 

The grid is specified using a comma separated list of **from:to:step**, where **step** is an optional argument (defaulting to one).
E.g. the _1:5,10:50:5,60:100:10_ argument specifies the 1,2,3,4,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100 grid.

**Note** this analysis may take a long time to run.

The output will look much like in the previous two examples, except that the feature list will
be quite long (abbreviated here to features ranked between 1 and 5 on average). 
IDEPI selected **N=90** as the number of features yielding the best mean MCC (0.89).  
> PNGS is an abbreviation for Potential N-linked Glycosylation Site.
> See the _tropism_ or _bnab_ examples for a more detailed explanation.

    {
      "metadata": {
        "antibodies":  [ "HAD" ],
        "balance":     0.702671,
        "features":    90,
        "folds":       5,
        "label":      "max(IC50) > 20",
        "sequences":   861
      },
      "statistics": {
        "accuracy":    { "mean": 0.952373, "std": 0.012868 },
        "f1score":     { "mean": 0.965986, "std": 0.009780 },
        "mcc":         { "mean": 0.887940, "std": 0.028415 },
        "npv":         { "mean": 0.917860, "std": 0.088268 },
        "ppv":         { "mean": 0.968902, "std": 0.031513 },
        "sensitivity": { "mean": 0.963636, "std": 0.043102 },
        "specificity": { "mean": 0.925943, "std": 0.079170 }
      },
      "weights": [
        { "position": "T297K",           "N": 5, "rank": { "mean":  1.60, "std":  2.68 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "R298D",           "N": 5, "rank": { "mean":  2.80, "std":  2.19 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "PNGS(N448)",      "N": 5, "rank": { "mean":  2.80, "std":  1.67 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "PNGS(T188)",      "N": 5, "rank": { "mean":  4.00, "std":  2.45 }, "value": { "mean":  1.00, "std": 0.00 } }
        ...
      ]
    }

## Learning the model for HAD.

Unlike **idepi discrete** which uses nested cross-validation to evaluate the robustness of the 
features, **idepi learn** uses the complete data set (although it still uses cross-validation
to tune internal model, e.g. LSVM parameters) to build a predictive model and save it to file
so that sequences of unknown phenotypes can be classified later.

The arguments are exactly the same as those passed to **idepi discrete**, except the final
positional argument is required to specify the name of the file (in this case **models/HAD.model**)
to save the IDEPI model to. We build a model on **N=90** features here, based on the grid 
search from the previous step.

    $idepi learn --csv training/HIVSeqDB.fasta training/HAD.csv --numfeats 90 --mcc HAD models/HAD.model 

The output to stdout has the same content as for **idepi discrete**, except that there is only 
one "fold" (all of the training data). For brevity, only the top 5 ranked features are retained in the listing below.

    {
      "metadata": {
        "antibodies":  [ "HAD" ],
        "balance":     0.702671,
        "features":    90,
        "folds":       1,
        "label":      "max(IC50) > 20",
        "sequences":   861
      },
      "statistics": {
        "accuracy":    { "mean": 0.953542, "std": 0.000000 },
        "f1score":     { "mean": 0.967213, "std": 0.000000 },
        "mcc":         { "mean": 0.887901, "std": 0.000000 },
        "npv":         { "mean": 0.939024, "std": 0.000000 },
        "ppv":         { "mean": 0.959350, "std": 0.000000 },
        "sensitivity": { "mean": 0.975207, "std": 0.000000 },
        "specificity": { "mean": 0.902344, "std": 0.000000 }
      },
      "weights": [
        { "position": "T297K",           "N": 1, "rank": { "mean":  1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "I320[]",          "N": 1, "rank": { "mean":  2.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "R298D",           "N": 1, "rank": { "mean":  3.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "PNGS(N448)",      "N": 1, "rank": { "mean":  4.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "PNGS(T188)",      "N": 1, "rank": { "mean":  5.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } }
     ]
    }


## Predicting the HAD status of _validation_ sequences

In this section, we will apply the model derived by  **idepi learn** 
to the set of sequences from **10** individuals not used for training the model 
from the [Holman and Gabuzda](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0049538) paper, 
included in the IDEPI distribution in the **validation/HAD.fas** file.

> All of the validation samples have HAD (i.e are positive), and can contain multiple sequences from 
> the same individual. Following Holman and Gabuzda, we classify an individual as HAD+
> if the majority of sequences from that individual are classified as HAD+.

**Note** As noted previously, because the current version of IDEPI 
expects codon sequences as input, the original amino-acid sequences 
have been reverse translated into codon sequences. 

Execute:

    $idepi predict models/HAD.model validation/HAD.fas > validation/idepi.json
    
This **idepi predict** script expects to be given the model (saved by **idepi learn**) and
a collection of sequences (homologous to those used to train the model) as a FASTA file, not
necessarily.

The output, written to stdout, and redirected to **validation/idepi.json** in this example,
is a JSON file, where, for each input sequence, the JSON object contains a record (in the _predictions_ field) of the 
following form 
 
    { "id": "E21_2_HAD", 
      "value": 1, 
      "features": [ "L288L", "T297T", "R298R", "G314G", "R315G", 
     "F317F", "PNGS(N276)", "PNGS(N295)", "PNGS(N448)", "PNGS(K130+N276)", 
     "PNGS(K130+N295)", "PNGS(K130+N301)", "PNGS(K130+N448)", "PNGS(N160+N339)", 
     "PNGS(N197+K362)", "PNGS(N276+N289)", "PNGS(N289+N386)", "PNGS(N289+N392)", 
     "PNGS(N301+N356)", "PNGS(N301+N386)", "PNGS(N339+N637)", "PNGS(N356+N386)", 
     "PNGS(N356+N448)", "PNGS(K362+N448)" ] }
 
    
+ **id** is the name of the sequence   
+ **value** reports the phenotype classification, defined according the the label **max(IC50)>=20**
used to train the model. A value of "1" means the label is True, a value of "-1" means the value is 
False (susceptible in our case).
+ **features** lists all the model features present in the current sequence; for LSVM classifiers,
a weighted sum of these features defines the final score.

To evaluate IDEPI's predictions, we can use the auxiliary script found in **helpers/tabulate_validation.py**.
The script compares the phenotypes stored in **validation/idepi.json** to those encoded in the 
sequence names themselves (_HAD_), and calculates the confusion table and various performance
metrics. The script also collapses all the different isolates from the same individual 
and makes the majority call.

    $python3 helpers/tabulate_validation.py -j validation/idepi.json
    
The output (compare to Table 4 in the manuscript) is as follows:
    
    --Individual 6568 (HAD). Prediction: 85.7143% HAD
    --Individual 55 (HAD). Prediction: 100% HAD
    --Individual 7766 (HAD). Prediction: 85.7143% HAD
    --Individual 47 (HAD). Prediction: 100% HAD
    --Individual E21 (HAD). Prediction: 100% HAD
    --Individual 61 (HAD). Prediction: 100% HAD
    --Individual 60 (HAD). Prediction: 100% HAD
    --Individual 62 (HAD). Prediction: 100% HAD
    --Individual CA110 (HAD). Prediction: 100% HAD
    --Individual 59 (HAD). Prediction: 100% HAD

    N = 10

    IDEPI     Measured     Count   Percentage 
    Non-HAD   Non-HAD           0    0 
    HAD       HAD              10    100 
    Non-HAD   HAD               0    0 
    HAD       Non-HAD           0    0

    Specificity 1 
    Sensitivity 1 
    Accuracy    1


