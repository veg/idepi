LoadFunctionLibrary ("ReadDelimitedFiles");
LoadFunctionLibrary ("NJ");
LoadFunctionLibrary ("p_Distance_aa");

SetDialogPrompt 	("Please specify an amino-acid file:");

DataSet				inputAlignment     = ReadDataFile (PROMPT_FOR_FILE);

dataFilePath		=  LAST_FILE_PATH;

DataSetFilter		filteredData	   = CreateFilter (inputAlignment,1);

// check to see if this is a protein alignment 
GetDataInfo (_AncestalFilterChars, filteredData, "CHARACTERS");
assert (Columns (_AncestalFilterChars) == 20, "Expected a protein alignment");

// clean up sequence names

newNamesToOldNames = {};
usedNames		   = {};

fprintf (stdout, "\nConverting sequence names to valid HyPhy identifiers...\n");

for (s = 0; s < filteredData.species; s += 1)
{
	GetString 					  (thisSequenceName, filteredData, s);
	newName 					= normalizeSequenceID (thisSequenceName, "usedNames");
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

LoadFunctionLibrary ("HIVbetween+F", {"00":"Fixed Rates"});
Tree aTree = treeString;

LikelihoodFunction lf = (filteredData, aTree);
Optimize (res, lf);

fprintf (stdout, lf);

fprintf (stdout, "\nReconstructing ancestors (this could also take a bit of time)...\n");

DataSet	 		_marginalAncestors 			= ReconstructAncestors (lf,MARGINAL,DOLEAVES);
DataSetFilter	_marginalAncestorsFilter	= CreateFilter 		   (_marginalAncestors, 1);
GetDataInfo 	(_marginalFilterSiteToPatternMap, filteredData);
GetString		(_AncestralNodeNames, _marginalAncestorsFilter, -1);

_idx_3 = 0;
_characterDimension 	= Columns (_AncestalFilterChars);


_marginalInformation	= {};
/* [(i,j)] -> {chars,1} - marginal support for each character */

for (_idx_1 = 0; _idx_1 < _marginalAncestorsFilter.species; _idx_1 = _idx_1 + 1)
{
	for (_idx_2 = 0; _idx_2 < _marginalAncestorsFilter.sites; _idx_2 = _idx_2 + 1)
	{
		_patternIndex 				 = _marginalFilterSiteToPatternMap[_idx_2];
		_marginalInformation[_idx_3] = _marginalAncestors.marginal_support_matrix[{{_idx_1,_patternIndex*_characterDimension}}][{{_idx_1,(1+_patternIndex)*_characterDimension-1}}];
		_idx_3 						 = _idx_3+1;
	}

}

_outputCSV = ""; _outputCSV * 8192; 
/* _outputCSV * ("Sequence,[Site:["+Join(",",_AncestalFilterChars)+"]]\n"); */

_outputCSV * ("[\n");

_idx_3 = 0;
for (_idx_1 = 0; _idx_1 < _marginalAncestorsFilter.species; _idx_1 = _idx_1 + 1)
{
    if (_idx_1 > 0)
    {
        _outputCSV * (",\n");
    }

    _outputCSV * ("{\"id\":\"" + newNamesToOldNames[_AncestralNodeNames[_idx_1]] + "\",\"values\":[");

	for (_idx_2 = 0; _idx_2 < _marginalAncestorsFilter.sites; _idx_2 = _idx_2 + 1)
	{
        /* trailing commas */
        if (_idx_2 > 0)
        {
            _outputCSV * (",");
        }

        _outputCSV * ("[");
        /* _outputCSV * ("\n" + newNamesToOldNames[_AncestralNodeNames[_idx_1]] + "," + (1+_idx_2)); */

		_maxValue = 0;
		_maxIndex = 0;

		for (_idx_4 = 0; _idx_4 < _characterDimension; _idx_4 = _idx_4 + 1)
		{
			_thisCharacter = (_marginalInformation[_idx_3])[_idx_4];
            
            if (_idx_4 > 0)
            {
                _outputCSV * (",");
            }
			
            _outputCSV * (""+_thisCharacter);
			
            if (_thisCharacter > _maxValue)
			{
				_maxValue = _thisCharacter;
				_maxIndex = _idx_4;
			}
		}
		
		_idx_3 = _idx_3 + 1;
        
        _outputCSV * ("]");
	}
    
    _outputCSV * ("]}");
}

_outputCSV * ("\n]\n");

_outputCSV * 0;

dataFilePath += "_correction.csv";

fprintf (dataFilePath, CLEAR_FILE, _outputCSV);
