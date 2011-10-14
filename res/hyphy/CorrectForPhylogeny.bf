MESSAGE_LOGGING = 0;

LoadFunctionLibrary ("ReadDelimitedFiles.bf");
LoadFunctionLibrary ("NJ.bf");
LoadFunctionLibrary ("p_Distance_aa.bf");

SetDialogPrompt ("Please specify an amino-acid file:");

DataSet       inputAlignment = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData   = CreateFilter (inputAlignment,1);

// check to see if this is a protein alignment 
GetDataInfo (_AncestralFilterChars, filteredData, "CHARACTERS");
assert (Columns (_AncestralFilterChars) == 20, "Expected a protein alignment");

// clean up sequence names

newNamesToOldNames = {};
usedNames          = {};

fprintf (stdout, "\nConverting sequence names to valid HyPhy identifiers...\n");

for (s = 0; s < filteredData.species; s += 1)
{
    GetString (thisSequenceName, filteredData, s);
    newName                     = normalizeSequenceID (thisSequenceName, "usedNames");
    newNamesToOldNames[newName] = thisSequenceName;
    if (thisSequenceName != newName)
    {
        SetParameter (inputAlignment, s, newName);
        fprintf (stdout, thisSequenceName, "=>", newName, "\n");
    }
}

fprintf (stdout, "\nBuilding a NJ tree...\n");

treeString = InferTreeTopology (1);
fprintf (stdout, "\nFitting an evolutionary model for HIV-1 sequences (this could take a bit of time)...\n");

LoadFunctionLibrary ("HIVbetween+F.bf", {"00":"Fixed Rates"});
Tree aTree = treeString;

LikelihoodFunction lf = (filteredData, aTree);
Optimize (res, lf);

fprintf (stdout, lf);

fprintf (stdout, "\nReconstructing ancestors (this could also take a bit of time)...\n");

DataSet       _marginalAncestors       = ReconstructAncestors (lf,MARGINAL,DOLEAVES);
DataSetFilter _marginalAncestorsFilter = CreateFilter (_marginalAncestors, 1);

GetDataInfo (_marginalFilterSiteToPatternMap, filteredData);
GetString (_AncestralNodeNames, _marginalAncestorsFilter, -1);

_characterDimension = Columns (_AncestralFilterChars);

_strLen = 0;
for (_idx_1 = 0; _idx_1 < _marginalAncestorsFilter.species; _idx_1 += 1)
{
    // +22 for the ' <UNKNOWN DESCRIPTION>' that's appended stupidly by HyPhy
    _accum = 1 + 22;
    for (_idx_2 = 1; (1. + _idx_1) / _idx_2 > 1.; _idx_2 = _idx_2 * 10)
    {
        _accum += 1;
    }
    _strLen += _accum;
}
// -1 to remove space for the leading comma
// is handled by the 0 above

_msm = {_marginalAncestorsFilter.species,_marginalAncestorsFilter.sites*_characterDimension};
_bim = {_marginalAncestorsFilter.species,_marginalAncestorsFilter.sites*_characterDimension};
_ids = "";
_ids * _strLen;

for (_idx_1 = 0; _idx_1 < _marginalAncestorsFilter.species; _idx_1 += 1)
{
       if (_idx_1 > 0)
       {
            _ids * ",";
       }
       _ids * ("" + newNamesToOldNames[_AncestralNodeNames[_idx_1]]);

       for (_idx_2 = 0; _idx_2 < _marginalAncestorsFilter.sites; _idx_2 += 1)
       {
               _patternIndex = _marginalFilterSiteToPatternMap[_idx_2];

               _charProbs = _marginalAncestors.marginal_support_matrix[{{_idx_1,_patternIndex*_characterDimension}}][{{_idx_1,(1+_patternIndex)*_characterDimension-1}}];
               // _charProbs = _charProbs["(1-_MATRIX_ELEMENT_VALUE_)"];

               GetDataInfo (_charsPresent, filteredData, _idx_1, _patternIndex);

               for (_idx_3 = 0;  _idx_3 < _characterDimension; _idx_3 += 1)
               {
                       _msm[_idx_1][_idx_2 * _characterDimension + _idx_3] = _charProbs[_idx_3];
                       _bim[_idx_1][_idx_2 * _characterDimension + _idx_3] = _charsPresent[_idx_3];
               }
       }

}

function _THyPhyAskFor(key)
{
    if (key == "ids")
    {
        return _ids;
    }
    if (key == "order")
    {
        return "" + Join(",",_AncestralFilterChars);
    }
    if (key == "marginalSupportMatrix")
    {
        return _msm;
    }
    if (key == "binaryIdentitiesMatrix")
    {
        return _bim;
    }
    if (key == "numSites")
    {
        return _marginalAncestorsFilter.sites;
    }
    if (key == "numChars")
    {
        return _characterDimension;
    }
    if (key == "numSpecies")
    {
        return  _marginalAncestorsFilter.species;
    }
    return "_THyPhy_NOT_HANDLED_";
}
