# Identifying epitopes and predicting neutralization by anti-HIV broadly neutralizing antibodies (bNab) 

In this tutorial, we will explain how IDEPI can be used for its **primary purpose**: 
identifying the epitopes of (monoclonal) broadly neutralizing antibodies (bNab).

The data used for the training and validation come from three separate sources.

1. The 181 (we were able to find sequence data from 178 of them) HIV-1 strains and IC50 titers from 
the _reference_ screening panel, used in [a 2013 J Virol paper by Chuang et al](http://jvi.asm.org/content/87/18/10047.short)
2. The 247 sequences used to train (and 55 -- to validate) a b12 bNab classifier from [a 2010 PLoS Comp Biol paper by Gnanakaran et al](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000955)
3. Publicly available sequences for the PGT-121 bNab (mostly from [a 2011 Nature paper by Walker et al] (http://www.ncbi.nlm.nih.gov/pubmed/21849977)

These data are available either in **training** and **validation** directories, or
as a part of the built-in neutralization database distributed with IDEPI.

All the code examples assume that the current 
working directory is **data/bnab** (inside the IDEPI
directory). 

## Identifying the single most predictive feature of neutralization/susceptibility 

We use the **idepi discrete** module identify the single top features in 
the training data associated with the neutralization (or resistance) to 
neutralization by a bNab. 

### PGT-121.

Issue the following command:

    $idepi discrete --numfeats 1 --mcc --clonal PGT-121 
    
The arguments above tell IDEPI to 

+ train models with a single feature
+ optimize Matthews Correlation Coefficient (mcc) when tuning LSVM parameters
+ use the IC50 data for the _PGT-121_ antibody to define sequence phenotypes

Note that unlike other examples, this call does not explicitly tell IDEPI where 
to find the sequences or their associated IC50 titers. This is because IDEPI ships with 
a prebuilt [SQLite](http://www.sqlite.org) database (which can be found 
at **idepi/data/allneuts.sqlite3** in the IDEPI source directory). You
can type some random string instead of PGT-121 to have IDEPI return the list of 
known antibodies (not all of which are actual bNab), currently 

>'10E8', '17b', '2F5', '2G12', '3BNC117', '447-52D', '447D', '4E10', '5-MAB COCKTAIL', '58.2', '8ANC131', '8ANC195', 'B6', 'C1', 'C10', 'C4', 'CD4-IGG2', 'CH103', 'E17', 'E3', 'E6', 'F105', 'F425', 'HIVIG', 'HIVIG-C', 'HMB-H-P1', 'HMB-H-P2', 'HMB-H-P3', 'LTNP2', 'N16', 'N16 CLADE B HIV+ PLA', 'N16-TX*', 'NAC17', 'NAC18', 'NIH45-46', 'PG16', 'PG9', 'PGT-121', 'PGT-122', 'PGT-123', 'PGT-125', 'PGT-126', 'PGT-127', 'PGT-128', 'PGT-130', 'PGT-131', 'PGT-135', 'PGT-136', 'PGT-137', 'PGT-141', 'PGT-142', 'PGT-143', 'PGT-144', 'PGT-145', 'PGT-151', 'PGT-152', 'PGV04', 'PLASMA 1', 'PLASMA 2', 'PLASMA 3', 'PLASMA 4', 'PLASMA 5', 'SCD4', 'TRIMAB', 'VRC-CH31', 'VRC-PG04', 'VRC-PG20', 'VRC01', 'VRC03', 'VRC06', 'VRC23', 'X5', 'Z1648', 'Z1652', 'Z1679', 'Z1682', 'Z1684', 'Z2', 'Z23', 'Z85', 'Z87', 'b12', 'b13'


Because bNab neutralization is determined by HIV-1 envelope,
IDEPI's existing HMMER model is used to align training
sequences to the standard HXB2 coordinates. 

The output of the program is written to stdout as a JSON file.

    {
      "metadata": {
        "antibodies":  [ "PGT-121" ],
        "balance":     0.359375,
        "features":    1,
        "folds":       5,
        "label":      "max(IC50) > 20",
        "sequences":   128
      },
      "statistics": {
        "accuracy":    { "mean": 0.796615, "std": 0.172615 },
        "f1score":     { "mean": 0.727302, "std": 0.260070 },
        "mcc":         { "mean": 0.572386, "std": 0.400268 },
        "npv":         { "mean": 0.874706, "std": 0.193316 },
        "ppv":         { "mean": 0.685354, "std": 0.203517 },
        "sensitivity": { "mean": 0.780000, "std": 0.355000 },
        "specificity": { "mean": 0.805147, "std": 0.099740 }
      },
      "weights": [
        { "position": "PNGS(N301+N332)", "N": 5, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } }
      ]
    }


    
1. **_metadata_** reports the information about the settings and data used by IDEPI
    + **antibodies** (PGT-121) is the antibody being studied
    + **balance** (0.36) describes the proportion of sequences defined as having the
    "resistant" phenotype (defined by the 'label' field); here a sequenced is resistant
    if its neutralization IC50 is >20 micrograms/milliliter. 
    + **features** (1) reports the number of features used by IDEPI to build the model
    + **folds** (5) the (default) number of cross-validation folds used by IDEPI to 
    compute model performance 
    + **label** ("max(IC50)>=20") the Python formula used by IDEPI to convert phenotype
    measures (IC50) into a binary phenotype. 
    + **sequences** (128) the total number of (sequence, phenotype) pairs used by IDEPI
    to build the model

2. **_statistics_** records all the standard model performance metrics computed with nested
    cross validation. The mean and standard deviation for each metric (in this case over
    5 CV folds) are reported. E.g. the mean Matthew's Correlation Coefficient (mcc) is
    0.57 (the value may differ from run to run because of the stochastic nature of CV).

3. **_weights_** is the combined report of all the features found by the classification 
procedure in any of the CV folds. 
    + **position** the feature selected by the model. In this case the single feature 
    selected is a pair of potential N-linked glycosylation sites, anchored at reference 
    positions 301 and 332. **Any** sequence change which disrupts **either** of the two 
    PNGS motifs will be interpreted as conferring resistance.
    + **N** in how many of the cross validation folds was this feature selected (the greater
    the value, the more robust the feature)
    + **rank** the rank of this feature (mean and st.dev) in the models derived by the CV 
    folds (see below for more interesting examples). H
    + **value** does this feature correspond to susceptibility (value -1) or resistance
    (value 1) in each of the CV fold. In this example, a sequence has both PNGS motifs
    is classified as susceptible.

### PG-9

Issue the following command:

    $idepi discrete --mcc --numfeats 1 --clonal PG9
    
Differences from the PGT-121 case:

+ We use the **--clonal** argument to tell IDEPI to only use sequences labeled as clonal (as opposed to "bulk")
+ The name of the antibody has changed

Output:
    
    {
      "metadata": {
        "antibodies":  [ "NAC17", "PG9" ],
        "balance":     0.262458,
        "features":    1,
        "folds":       5,
        "label":      "max(IC50) > 20",
        "sequences":   301
      },
      "statistics": {
        "accuracy":    { "mean": 0.780601, "std": 0.073823 },
        "f1score":     { "mean": 0.421271, "std": 0.113114 },
        "mcc":         { "mean": 0.364691, "std": 0.214812 },
        "npv":         { "mean": 0.792969, "std": 0.028678 },
        "ppv":         { "mean": 0.742540, "std": 0.420131 },
        "sensitivity": { "mean": 0.303333, "std": 0.097432 },
        "specificity": { "mean": 0.950202, "std": 0.098824 }
      },
      "weights": [
        { "position": "PNGS(N160)",      "N": 3, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "PNGS(N160+N197)", "N": 1, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "PNGS(N160+N301)", "N": 1, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } }
      ]
    }

> The **antibodies** field lists _two_ entries; this is because the two refer to the same antibody.
  
       
### b12

Issue the following command:

    $idepi discrete --csv training/lanl_b12.fas training/lanl_b12.csv  --mcc --numfeats 1 --label "max(IC50)>=50" B12
    
Differences from the PGT-121 case:

+ We specify a pair of external files, through the **--csv** argument, which contain the sequences and matching IC50 titers
+ The definition of a resistant sequence now requires IC50 to be >= 50 (default is >20) (the **--label** argument); this is
done to match the definition in [Gnanakaran et al](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000955).
+ The antibody (B12) matches the column name available in **training/lanl_b12.csv**

Output:
    
    {
      "metadata": {
        "antibodies":  [ "B12" ],
        "balance":     0.643725,
        "features":    1,
        "folds":       5,
        "label":      "max(IC50)>=50",
        "sequences":   247
      },
      "statistics": {
        "accuracy":    { "mean": 0.632735, "std": 0.290497 },
        "f1score":     { "mean": 0.663271, "std": 0.438427 },
        "mcc":         { "mean": 0.255810, "std": 0.412483 },
        "npv":         { "mean": 0.527773, "std": 0.262263 },
        "ppv":         { "mean": 0.728007, "std": 0.164654 },
        "sensitivity": { "mean": 0.642540, "std": 0.531378 },
        "specificity": { "mean": 0.613072, "std": 0.202857 }
      },
      "weights": [
        { "position": "I161M", "N": 1, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "D185D", "N": 4, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } }
      ]
    }

### 8ANC131

Issue the following command:

    $idepi discrete --csv training/181.fas training/181.csv --numfeats 1 --mcc 8ANC131
    
Differences from the PGT-121 case:

+ We specify a pair of external files, through the **--csv** argument, which contain the sequences and matching IC50 titers
+ The antibody (8ANC131) matches the column name available in **training/181.csv**

Output:
    
    {
      "metadata": {
        "antibodies":  [ "8ANC131" ],
        "balance":     0.308989,
        "features":    1,
        "folds":       5,
        "label":      "max(IC50) > 20",
        "sequences":   178
      },
      "statistics": {
        "accuracy":    { "mean": 0.617778, "std": 0.261952 },
        "f1score":     { "mean": 0.411729, "std": 0.228513 },
        "mcc":         { "mean": 0.150087, "std": 0.369735 },
        "npv":         { "mean": 0.729125, "std": 0.108926 },
        "ppv":         { "mean": 0.438333, "std": 0.294717 },
        "sensitivity": { "mean": 0.436364, "std": 0.325246 },
        "specificity": { "mean": 0.698667, "std": 0.456748 }
      },
      "weights": [
        { "position": "N92N",            "N": 1, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "I646L",           "N": 1, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "PNGS(N339+Q442)", "N": 3, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } }
      ]
    }
      
### 8ANC195

Issue the following command:

    $idepi discrete --csv training/181.fas training/181.csv --numfeats 1 --mcc 8ANC195
    
Differences from the PGT-121 case:

+ We specify a pair of external files, through the **--csv** argument, which contain the sequences and matching IC50 titers
+ The antibody (8ANC195) matches the column name available in **training/181.csv**

Output:
    
    {
      "metadata": {
        "antibodies":  [ "8ANC195" ],
        "balance":     0.426966,
        "features":    1,
        "folds":       5,
        "label":      "max(IC50) > 20",
        "sequences":   178
      },
      "statistics": {
        "accuracy":    { "mean": 0.797460, "std": 0.187565 },
        "f1score":     { "mean": 0.735678, "std": 0.265119 },
        "mcc":         { "mean": 0.585957, "std": 0.382458 },
        "npv":         { "mean": 0.788595, "std": 0.181354 },
        "ppv":         { "mean": 0.819534, "std": 0.207757 },
        "sensitivity": { "mean": 0.674167, "std": 0.310935 },
        "specificity": { "mean": 0.890952, "std": 0.132583 }
      },
      "weights": [
        { "position": "PNGS(N234+N276)", "N": 5, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } }
      ]
    }

### 10E8

Issue the following command:

    $idepi discrete --csv training/181.fas training/181.csv --numfeats 1 --mcc --label "max(IC50)>=5" 10E8
    
Differences from the PGT-121 case:

+ We specify a pair of external files, through the **--csv** argument, which contain the sequences and matching IC50 titers
+ The definition of a resistant sequence now requires IC50 to be >= 5 (default is >20) (the **--label** argument); this is
done because 10E8 is unusually potent. 
+ The antibody (10E8) matches the column name available in **training/181.csv**

Output (not the complete lack of predictive power here: MCC is negative, ppv is 0, and no two replicates agree on the top feature, this is because
only 4% of the training set are 'positive', i.e. resistant):
    
     {
      "metadata": {
        "antibodies":  [ "10E8" ],
        "balance":     0.0393258,
        "features":    1,
        "folds":       5,
        "label":      "max(IC50)>=5",
        "sequences":   178
      },
      "statistics": {
        "accuracy":    { "mean":  0.841746, "std": 0.325113 },
        "f1score":     { "mean":  0.000000, "std": 0.000000 },
        "mcc":         { "mean": -0.059083, "std": 0.115424 },
        "npv":         { "mean":  0.954865, "std": 0.031811 },
        "ppv":         { "mean":  0.000000, "std": 0.000000 },
        "sensitivity": { "mean":  0.000000, "std": 0.000000 },
        "specificity": { "mean":  0.876639, "std": 0.341713 }
      },
      "weights": [
        { "position": "E153Q",           "N": 1, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "K171E",           "N": 1, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "K305K",           "N": 1, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "H564H",           "N": 1, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "PNGS(T413+E824)", "N": 1, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } }
      ]
    }

    
## Perform a grid search over the number of features

For a more detailed description, please see the **DRAM**, **tropism**, or **dementia** examples.

### PGT-121

Command

    $idepi discrete --mcc --numfeats 1:5,10:50:5,60:100:10 PGT-121

Output (limited to top 5 ranking features)

    {
      "metadata": {
        "antibodies":  [ "PGT-121" ],
        "balance":     0.359375,
        "features":    1,
        "folds":       5,
        "label":      "max(IC50) > 20",
        "sequences":   128
      },
      "statistics": {
        "accuracy":    { "mean": 0.796615, "std": 0.172615 },
        "f1score":     { "mean": 0.727302, "std": 0.260070 },
        "mcc":         { "mean": 0.572386, "std": 0.400268 },
        "npv":         { "mean": 0.874706, "std": 0.193316 },
        "ppv":         { "mean": 0.685354, "std": 0.203517 },
        "sensitivity": { "mean": 0.780000, "std": 0.355000 },
        "specificity": { "mean": 0.805147, "std": 0.099740 }
      },
      "weights": [
        { "position": "PNGS(N301+N332)", "N": 5, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } }
      ]
    }
    
### PG-9

Command

    $idepi discrete --mcc --numfeats 1:5,10:50:5,60:100:10 --clonal PG9

Output (limited to top 10 ranking features)

    {
      "metadata": {
        "antibodies":  [ "NAC17", "PG9" ],
        "balance":     0.262458,
        "features":    60,
        "folds":       5,
        "label":      "max(IC50) > 20",
        "sequences":   301
      },
      "statistics": {
        "accuracy":    { "mean": 0.780765, "std": 0.094978 },
        "f1score":     { "mean": 0.572046, "std": 0.168912 },
        "mcc":         { "mean": 0.427343, "std": 0.230270 },
        "npv":         { "mean": 0.845248, "std": 0.050251 },
        "ppv":         { "mean": 0.592918, "std": 0.206827 },
        "sensitivity": { "mean": 0.556667, "std": 0.155300 },
        "specificity": { "mean": 0.860707, "std": 0.099602 }
      },
      "weights": [
      ...
        { "position": "E47E",            "N": 1, "rank": { "mean":  2.00, "std":  0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "V84V",            "N": 1, "rank": { "mean":  2.00, "std":  0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "PNGS(N160)",      "N": 5, "rank": { "mean":  7.00, "std": 16.49 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "V169E",           "N": 5, "rank": { "mean":  5.00, "std":  5.66 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "A336V",           "N": 1, "rank": { "mean":  8.00, "std":  0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "W395T",           "N": 1, "rank": { "mean":  6.00, "std":  0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "N677Q",           "N": 1, "rank": { "mean":  7.00, "std":  0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "N809I",           "N": 1, "rank": { "mean":  9.00, "std":  0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "PNGS(T138+N816)", "N": 1, "rank": { "mean":  6.00, "std":  0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "PNGS(N160+N301)", "N": 5, "rank": { "mean":  6.20, "std":  6.84 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "PNGS(N234+N332)", "N": 2, "rank": { "mean":  9.00, "std":  7.07 }, "value": { "mean":  0.00, "std": 1.41 } },
        { "position": "PNGS(N332+G404)", "N": 1, "rank": { "mean":  6.00, "std":  0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "PNGS(H625+N637)", "N": 1, "rank": { "mean":  4.00, "std":  0.00 }, "value": { "mean":  1.00, "std": 0.00 } }
       ...
      ]
    }
    
### 8ANC131

Command
    
    $idepi discrete --csv training/181.fas training/181.csv --numfeats 1:5,10:50:5,60:100:10  --mcc 8ANC131

Output (limited to top 5 ranking features)

    {
      "metadata": {
        "antibodies":  [ "8ANC131" ],
        "balance":     0.308989,
        "features":    15,
        "folds":       5,
        "label":      "max(IC50) > 20",
        "sequences":   178
      },
      "statistics": {
        "accuracy":    { "mean": 0.634603, "std": 0.147954 },
        "f1score":     { "mean": 0.445763, "std": 0.282741 },
        "mcc":         { "mean": 0.193922, "std": 0.334608 },
        "npv":         { "mean": 0.765776, "std": 0.146224 },
        "ppv":         { "mean": 0.423284, "std": 0.201359 },
        "sensitivity": { "mean": 0.509091, "std": 0.456360 },
        "specificity": { "mean": 0.691000, "std": 0.275031 }
      },
      "weights": [
        ...
        { "position": "PNGS(142a)",      "N": 1, "rank": { "mean":  2.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "K151G",           "N": 5, "rank": { "mean":  5.00, "std": 5.66 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "I201V",           "N": 1, "rank": { "mean":  5.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "V255V",           "N": 1, "rank": { "mean":  4.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "N463K",           "N": 3, "rank": { "mean":  3.67, "std": 2.94 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "D474N",           "N": 3, "rank": { "mean":  3.33, "std": 1.63 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "R476K",           "N": 1, "rank": { "mean":  5.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "I646L",           "N": 1, "rank": { "mean":  1.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "H720L",           "N": 1, "rank": { "mean":  5.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "Q805R",           "N": 1, "rank": { "mean":  5.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "PNGS(T138+Q442)", "N": 2, "rank": { "mean":  4.50, "std": 2.12 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "PNGS(142a+N234)", "N": 3, "rank": { "mean":  3.33, "std": 1.63 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "PNGS(N332+N816)", "N": 1, "rank": { "mean":  2.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "PNGS(N339+Q442)", "N": 4, "rank": { "mean":  3.75, "std": 9.53 }, "value": { "mean":  1.00, "std": 0.00 } }
        ...
      ]
    }

### 8ANC195

Command
    
    $idepi discrete --csv training/181.fas training/181.csv --numfeats 1:5,10:50:5,60:100:10  --mcc 8ANC195

Output

    {
      "metadata": {
        "antibodies":  [ "8ANC195" ],
        "balance":     0.426966,
        "features":    2,
        "folds":       5,
        "label":      "max(IC50) > 20",
        "sequences":   178
      },
      "statistics": {
        "accuracy":    { "mean": 0.826667, "std": 0.156547 },
        "f1score":     { "mean": 0.822995, "std": 0.152105 },
        "mcc":         { "mean": 0.677084, "std": 0.304943 },
        "npv":         { "mean": 0.937402, "std": 0.153194 },
        "ppv":         { "mean": 0.735263, "std": 0.156634 },
        "sensitivity": { "mean": 0.935833, "std": 0.153161 },
        "specificity": { "mean": 0.745714, "std": 0.171739 }
      },
      "weights": [
        { "position": "PNGS(N160+N230)", "N": 5, "rank": { "mean": 2.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "PNGS(N234+N276)", "N": 5, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } }
      ]
    }
    
### 10E8

Command
    
    $idepi discrete --csv training/181.fas training/181.csv --numfeats 1:5,10:50:5,60:100:10  --mcc --label "max(IC50)>=5" 10E8

Output

    {
      "metadata": {
        "antibodies":  [ "10E8" ],
        "balance":     0.0393258,
        "features":    5,
        "folds":       5,
        "label":      "max(IC50)>=5",
        "sequences":   178
      },
      "statistics": {
        "accuracy":    { "mean": 0.932857, "std": 0.128226 },
        "f1score":     { "mean": 0.250000, "std": 0.866025 },
        "mcc":         { "mean": 0.234942, "std": 0.879120 },
        "npv":         { "mean": 0.970606, "std": 0.039590 },
        "ppv":         { "mean": 0.233333, "std": 0.869227 },
        "sensitivity": { "mean": 0.300000, "std": 0.894427 },
        "specificity": { "mean": 0.958824, "std": 0.128876 }
      },
      "weights": [
        { "position": "E153Q",           "N": 4, "rank": { "mean": 3.25, "std": 3.57 }, "value": { "mean":  0.50, "std": 1.73 } },
        { "position": "K171E",           "N": 4, "rank": { "mean": 2.50, "std": 2.24 }, "value": { "mean":  0.50, "std": 1.73 } },
        { "position": "K305K",           "N": 1, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "Q363S",           "N": 1, "rank": { "mean": 5.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "P369V",           "N": 2, "rank": { "mean": 3.50, "std": 0.71 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "P417P",           "N": 1, "rank": { "mean": 3.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "H564H",           "N": 1, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "T676S",           "N": 1, "rank": { "mean": 3.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "T676T",           "N": 3, "rank": { "mean": 3.33, "std": 2.16 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "A792A",           "N": 1, "rank": { "mean": 4.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "PNGS(D137+K362)", "N": 1, "rank": { "mean": 2.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "PNGS(N234+N356)", "N": 1, "rank": { "mean": 4.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "PNGS(K362+E824)", "N": 1, "rank": { "mean": 5.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "PNGS(T413+E824)", "N": 3, "rank": { "mean": 2.33, "std": 2.16 }, "value": { "mean":  1.00, "std": 0.00 } }
      ]
    }

    
### b12

Command

    $idepi discrete --csv training/lanl_b12.fas training/lanl_b12.csv  --mcc --numfeats 1:5,10:50:5,60:100:10 --label "max(IC50)>=50" B12
    
Output (limited to top 5 ranking features, note the considerable variability among the selected features)

    {
      "metadata": {
        "antibodies":  [ "B12" ],
        "balance":     0.643725,
        "features":    5,
        "folds":       5,
        "label":      "max(IC50)>=50",
        "sequences":   247
      },
      "statistics": {
        "accuracy":    { "mean": 0.696327, "std": 0.040963 },
        "f1score":     { "mean": 0.756736, "std": 0.044952 },
        "mcc":         { "mean": 0.355502, "std": 0.092011 },
        "npv":         { "mean": 0.568480, "std": 0.049800 },
        "ppv":         { "mean": 0.781864, "std": 0.064054 },
        "sensitivity": { "mean": 0.736089, "std": 0.101784 },
        "specificity": { "mean": 0.624837, "std": 0.169267 }
      },
      "weights": [
        { "position": "N80S",            "N": 1, "rank": { "mean": 2.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "M149N",           "N": 1, "rank": { "mean": 5.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "I161M",           "N": 1, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "Y173Y",           "N": 1, "rank": { "mean": 5.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "D185D",           "N": 5, "rank": { "mean": 1.80, "std": 3.58 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "T202K",           "N": 1, "rank": { "mean": 4.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "T399Y",           "N": 1, "rank": { "mean": 4.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "I424I",           "N": 1, "rank": { "mean": 2.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "I424V",           "N": 1, "rank": { "mean": 4.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "N425R",           "N": 1, "rank": { "mean": 3.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "I467T",           "N": 1, "rank": { "mean": 2.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "K502R",           "N": 1, "rank": { "mean": 4.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "V749A",           "N": 1, "rank": { "mean": 2.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "N750S",           "N": 1, "rank": { "mean": 3.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "L774F",           "N": 3, "rank": { "mean": 3.00, "std": 1.41 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "L774L",           "N": 1, "rank": { "mean": 5.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "A792G",           "N": 1, "rank": { "mean": 5.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "S804G",           "N": 1, "rank": { "mean": 3.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "PNGS(M149+N392)", "N": 1, "rank": { "mean": 3.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } }
      ]
    }

## Learning the model for b12.

For a more detailed description, please see the **DRAM**, **tropism**, or **dementia** examples.

    $idepi learn --csv training/lanl_b12.fas training/lanl_b12.csv --mcc --numfeats 5 --label "max(IC50)>=50" B12 models/B12.json

Output

    {
      "metadata": {
        "antibodies":  [ "B12" ],
        "balance":     0.643725,
        "features":    5,
        "folds":       1,
        "label":      "max(IC50)>=50",
        "sequences":   247
      },
      "statistics": {
        "accuracy":    { "mean": 0.748988, "std": 0.000000 },
        "f1score":     { "mean": 0.797386, "std": 0.000000 },
        "mcc":         { "mean": 0.471425, "std": 0.000000 },
        "npv":         { "mean": 0.630000, "std": 0.000000 },
        "ppv":         { "mean": 0.829932, "std": 0.000000 },
        "sensitivity": { "mean": 0.767296, "std": 0.000000 },
        "specificity": { "mean": 0.715909, "std": 0.000000 }
      },
      "weights": [
        { "position": "M149N", "N": 1, "rank": { "mean": 5.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "I161M", "N": 1, "rank": { "mean": 4.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "D185D", "N": 1, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "I424V", "N": 1, "rank": { "mean": 2.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "V812V", "N": 1, "rank": { "mean": 3.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } }
      ]
    }
   


## Predicting the neutralization of _validation_ sequences by b12

In this section, we will apply the model derived by  **idepi learn** 
to the set of **xx** sequences not used for training the model 
from the [Holman and Gabuzda](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0049538) paper, 
included in the IDEPI distribution in the **validation/lanl_b12_test.fas** file.


Execute:

    $idepi predict models/B12.json validation/lanl_b12_test.fas > validation/idepi.json
    

The output, written to stdout, and redirected to **validation/idepi.json** in this example,
is a JSON file, where, for each input sequence, the JSON object contains a record (in the _predictions_ field) of the 
following form 
 
    { "id": "HM215367_50",    "value":  1, "features": [ "M149N", "I161M",          "I424V"          ] },
  
    
+ **id** is the name of the sequence   
+ **value** reports the phenotype classification, defined according the the label **max(IC50)>=20**
used to train the model. A value of "1" means the label is True, a value of "-1" means the value is 
False (susceptible in our case).
+ **features** lists all the model features present in the current sequence; for LSVM classifiers,
a weighted sum of these features defines the final score.

To evaluate IDEPI's predictions, we can use the auxiliary script found in **helpers/tabulate_validation.py**.
The script compares the phenotypes stored in **validation/idepi.json** to those encoded in the 
sequence names themselves (IC50 values following an underscore), and calculates the confusion table and various performance
metrics. The script also collapses all the different isolates from the same individual 
and makes the majority call.

    $python3 helpers/tabulate_validation.py -j validation/idepi.json
    
The output (compare to Table 4 in the manuscript) is as follows:
    
    N = 55

    IDEPI         Measured     Count   Percentage 
    Susceptible   Susceptible     12    21.8 
    Resistant     Resistant       28    50.9 
    Resistant     Susceptible      7    12.7 
    Susceptible   Resistant        8    14.5

    Specificity 0.8 
    Sensitivity 0.78 
    Accuracy    0.73 
    MCC         0.4



