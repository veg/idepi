# Co-receptor usage/cell tropism prediction for HIV-1

In this tutorial, we apply IDEPI to the training and validation
data from a 2010 paper by [Dybowski et al](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000743)
on computational prediction of co-receptor usage from V3 sequence data.

> We would like to express our gratitude to Dr Hoffman, the senior author
> of this study, who was very kind to send us the sequences
> comprising the validation dataset.

All the code examples assume that the current 
working directory is **data/tropism** (inside the IDEPI
directory). 


## Preliminaries


The current version of IDEPI expects nucleotide sequences,
and because V3 data in published manuscripts is most frequently
provided as amino-acid data, it was necessary to "untranslate"
the protein data into IDEPI-compatible codon data (in a future
release IDEPI will support a flag to handle protein data directly).

The phenotype information in this classification problem is
binary, namely which coreceptor does a particular strain use.
To conform to IDEPI phenotype specifications, designed to handle
IC50 data, one needs to create a CSV file connecting each 
sequence ID with its phenotype as follows

    seqID, CXCR4, DUAL
    CCR5_C.ZA.1999.99ZATM7.12320,0,0
    CXCR4_B.NL.-.ACH168_C23.12936,50,0
    
The default sequence labeling algorithm in IDEPI will label 
a sequence "resistant" if its phenotype is >= 20, hence the value
of 50 for the CXCR4-using strain in the above example.


## Top two CXCR4 (X4) usage feature identification

We use the **idepi discrete** module identify the top two features in 
the training data associated with X4 tropism. 

### Fit the model with a single feature (**N=1**).

    $idepi discrete --csv training/V3.fas training/tropism.csv --numfeats 1 --mcc CXCR4 
    
The arguments above tell IDEPI to 

+ read paired sequence/phenotype data from training/V3.fas and training/tropism.csv, respectively
+ train models with a single feature
+ optimize Matthews Correlation Coefficient (mcc) when tuning LSVM parameters
+ use the _CXCR4_ column to define sequence phenotypes

Because co-receptor usage is determined by HIV-1 envelope,
IDEPI's existing HMMER model is used to align training
sequences to the standard HXB2 coordinates. 

The output of the program is written to stdout as a JSON file.

    {
      "metadata": {
        "antibodies":  [ "CXCR4" ],
        "balance":     0.15118,
        "features":    1,
        "folds":       5,
        "label":      "max(IC50) > 20",
        "sequences":   1356
      },
      "statistics": {
        "accuracy":    { "mean": 0.897498, "std": 0.014909 },
        "f1score":     { "mean": 0.593582, "std": 0.085067 },
        "mcc":         { "mean": 0.552214, "std": 0.080537 },
        "npv":         { "mean": 0.915507, "std": 0.015995 },
        "ppv":         { "mean": 0.739700, "std": 0.063764 },
        "sensitivity": { "mean": 0.497561, "std": 0.106873 },
        "specificity": { "mean": 0.968730, "std": 0.011221 }
      },
      "weights": [
        { "position": "PNGS(N301)", "N": 5, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } }
      ]
    }
    
1. **_metadata_** reports the information about the settings and data used by IDEPI
    + **antibodies** (CXCR4) should really be called "phenotypes", and leftover from the 
    the original IDEPI design, focusing only on antibody epitopes
    + **balance** (0.151) describes the proportion of sequences defined as having the
    "resistant" phenotype (defined by the 'label' field); here these are the CXCR4 (or Dual
    tropic) sequences
    + **features** (1) reports the number of features used by IDEPI to build the model
    + **folds** (5) the (default) number of cross-validation folds used by IDEPI to 
    compute model performance 
    + **label** ("max(IC50)>=20") the Python formula used by IDEPI to convert phenotype
    measures (IC50) into a binary phenotype. 
    + **sequences** (1356) the total number of (sequence, phenotype) pairs used by IDEPI
    to build the model

2. **_statistics_** records all the standard model performance metrics computed with nested
    cross validation. The mean and standard deviation for each metric (in this case over
    5 CV folds) are reported. E.g. the mean Matthew's Correlation Coefficient (mcc) is
    0.55 (the value may differ from run to run because of the stochastic nature of CV).

3. **_weights_** is the combined report of all the features found by the classification 
procedure in any of the CV folds. 
    + **position** the feature selected by the model (in this case the presence of potential
    N-linked glycosylation site at HXB2 coordinate 301).
    + **N** in how many of the cross validation folds was this feature selected (the greater
    the value, the more robust the feature)
    + **rank** the rank of this feature (mean and st.dev) in the models derived by the CV 
    folds (see below for more interesting examples). Here, in all 5 folds, the same feature
    (_PNGS301_) was selected.
    + **value** does this feature correspond to susceptibility (value -1) or resistance
    (value 1) in each of the CV fold. In this example, the _PNGS301_ feature was universally 
    associated with **not** using X4, i.e. a sequence that does not have a PNGS anchored 
    at 301 will be classified as X4 using.
    
### Fit the model with two features (**N=2**).

Only the **--numfeats** argument is changed.
    
    $ idepi discrete --csv training/V3.fas training/tropism.csv --numfeats 2 --mcc CXCR4
    
IDEPI reports the following model performance and features

    {
      "metadata": {
        "antibodies":  [ "CXCR4" ],
        "balance":     0.15118,
        "features":    2,
        "folds":       5,
        "label":      "max(IC50) > 20",
        "sequences":   1356
      },
      "statistics": {
        "accuracy":    { "mean": 0.917406, "std": 0.011134 },
        "f1score":     { "mean": 0.712167, "std": 0.055287 },
        "mcc":         { "mean": 0.665772, "std": 0.057383 },
        "npv":         { "mean": 0.943704, "std": 0.015122 },
        "ppv":         { "mean": 0.751448, "std": 0.025119 },
        "sensitivity": { "mean": 0.678049, "std": 0.093831 },
        "specificity": { "mean": 0.960034, "std": 0.007297 }
      },
      "weights": [
        { "position": "PNGS(N301)", "N": 5, "rank": { "mean": 1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "R306R",      "N": 5, "rank": { "mean": 2.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } }
      ]
    }
    
By adding a feature associated with X4 usage (_R306R_), IDEPI improves the MCC of the model from 0.55 to 0.67.
This feature is also a part of the [11/25 rule](http://www.prn.org/index.php/management/article/hiv_tropism_1002) for tropism detection.

### Direct IDEPI to perform a grid search over the number of features

Specify a grid of possible values to search over for the optimal number of features. 

    $idepi discrete --csv training/V3.fas training/tropism.csv --numfeats 1:5,10:50:5,60:100:10 --mcc CXCR4 

The grid is specified using a comma separated list of **from:to:step**, where **step** is an optional argument (defaulting to one).
E.g. the _1:5,10:50:5,60:100:10_ argument specifies the 1,2,3,4,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100 grid.

**Note** this analysis may take a long time to run.

The output will look much like in the previous two examples, except that the feature list will
be quite long (abbreviated here to features ranked between 1 and 10 on average). 
IDEPI selected **N=90** as the number of features yielding the best mean MCC (0.78).  

    {
      "metadata": {
        "antibodies":  [ "CXCR4" ],
        "balance":     0.15118,
        "features":    90,
        "folds":       5,
        "label":      "max(IC50) > 20",
        "sequences":   1356
      },
      "statistics": {
        "accuracy":    { "mean": 0.936567, "std": 0.025884 },
        "f1score":     { "mean": 0.810358, "std": 0.068809 },
        "mcc":         { "mean": 0.777690, "std": 0.080545 },
        "npv":         { "mean": 0.980145, "std": 0.008304 },
        "ppv":         { "mean": 0.742599, "std": 0.093563 },
        "sensitivity": { "mean": 0.892683, "std": 0.043631 },
        "specificity": { "mean": 0.944382, "std": 0.025013 }
      },
      "weights": [
        { "position": "PNGS(N301)", "N": 5, "rank": { "mean":  1.00, "std":  0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "R306R",      "N": 5, "rank": { "mean":  2.00, "std":  0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "N300Y",      "N": 5, "rank": { "mean":  3.20, "std":  0.89 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "N302N",      "N": 5, "rank": { "mean":  3.80, "std":  0.89 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "I323I",      "N": 5, "rank": { "mean":  5.60, "std":  1.79 }, "value": { "mean": -0.60, "std": 1.79 } },
        { "position": "T303T",      "N": 5, "rank": { "mean":  6.20, "std":  1.67 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "322e[]",     "N": 5, "rank": { "mean":  7.00, "std":  3.16 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "A316V",      "N": 5, "rank": { "mean":  8.80, "std":  3.85 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "T303I",      "N": 5, "rank": { "mean":  9.40, "std":  1.79 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "G321K",      "N": 5, "rank": { "mean":  9.40, "std":  5.02 }, "value": { "mean":  1.00, "std": 0.00 } }
      ]
    }


## Learning the model for X4 usage

Unlike **idepi discrete** which uses nested cross-validation to evaluate the robustness of the 
features, **idepi learn** uses the complete data set (although it still uses cross-validation
to tune internal model, e.g. LSVM parameters) to build a predictive model and save it to file
so that sequences of unknown phenotypes can be classified later.

The arguments are exactly the same as those passed to **idepi discrete**, except the final
positional argument is required to specify the name of the file (in this case **models/V3.model**)
to save the IDEPI model to. We build a model on **N=90** features here, based on the grid 
search from the previous step.

    $idepi learn --csv training/V3.fas training/tropism.csv --numfeats 90 --mcc CXCR4 models/V3.model 

The output to stdout has the same content as for **idepi discrete**, except that there is only 
one "fold" (all of the training data). For brevity, only the top 5 ranked features are retained in the listing below.

    {
      "metadata": {
        "antibodies":  [ "CXCR4" ],
        "balance":     0.15118,
        "features":    90,
        "folds":       1,
        "label":      "max(IC50) > 20",
        "sequences":   1356
      },
      "statistics": {
        "accuracy":    { "mean": 0.954277, "std": 0.000000 },
        "f1score":     { "mean": 0.857798, "std": 0.000000 },
        "mcc":         { "mean": 0.832775, "std": 0.000000 },
        "npv":         { "mean": 0.984000, "std": 0.000000 },
        "ppv":         { "mean": 0.809524, "std": 0.000000 },
        "sensitivity": { "mean": 0.912195, "std": 0.000000 },
        "specificity": { "mean": 0.961772, "std": 0.000000 }
      },
      "weights": [
        ...
        { "position": "N300Y",      "N": 1, "rank": { "mean":  3.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "PNGS(N301)", "N": 1, "rank": { "mean":  1.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } },
        { "position": "N302N",      "N": 1, "rank": { "mean":  4.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "R306R",      "N": 1, "rank": { "mean":  2.00, "std": 0.00 }, "value": { "mean":  1.00, "std": 0.00 } },
        { "position": "I323I",      "N": 1, "rank": { "mean":  5.00, "std": 0.00 }, "value": { "mean": -1.00, "std": 0.00 } }
        ...
      ]
    }


## Predicting the phenotype of _validation_ sequences

In this section, we will apply the model derived by  **idepi learn** 
to the set of **74** sequences not used for training the model (and 
whose phenotypes have been measured experimentally), from the 
[Dybowski et al](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000743)
paper, included in the  IDEPI distribution in the **validation/Dybowski.fas** file.

**Note** As noted previously, because the current version of IDEPI 
expects codon sequences as input, the original amino-acid sequences 
have been reverse translated into codon sequences. 

Execute:

    $idepi predict models/V3.model validation/Dybowski.fas > validation/idepi.json
    
This **idepi predict** script expects to be given the model (saved by **idepi learn**) and
a collection of sequences (homologous to those used to train the model) as a FASTA file, not
necessarily.

The output, written to stdout, and redirected to **validation/idepi.json** in this example,
is a JSON file, where, for each input sequence, the JSON object contains a record (in the _predictions_ field) of the 
following form 

    { "id": "V369_X4_", 
      "value": 1, 
      "features": [ "N300Y", "N301T", "N302N", "T303V", "I307T", "R308H", 
                    "P313P", "R315R", "F317Y", "V318Y", "I320T", "G321K", "322eI", "I323T", "M326I", "R327R", "Q328Q", "A329A",
                    "H330H" ] }
    
+ **id** is the name of the sequence   
+ **value** reports the phenotype classification, defined according the the label **max(IC50)>=20**
used to train the model. A value of "1" means the label is True, a value of "-1" means the value is 
False (susceptible in our case).
+ **features** lists all the model features present in the current sequence; for LSVM classifiers,
a weighted sum of these features defines the final score.

To evaluate IDEPI's predictions, we can use the auxiliary script found in **helpers/tabulate_validation.py**.
The script compares the phenotypes stored in **validation/idepi.json** to those encoded in the 
sequence names themselves (R5 vs X4), and calculates the confusion table and various performance
metrics.

    $python3 helpers/tabulate_validation.py -j validation/idepi.json
    
The output (compare to Table 4 in the manuscript) is as follows:
    
    N = 74

    IDEPI     Measured     Count   Percentage 
    CCR5      CCR5             53    71.6 
    CXCR4     CXCR4            14    18.9 
    CCR5      CXCR4             5    6.76 
    CXCR4     CCR5              2    2.7

    Specificity 0.74 
    Sensitivity 0.88 
    Accuracy    0.91 
    MCC         0.74
 
