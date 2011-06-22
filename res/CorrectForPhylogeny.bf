MESSAGE_LOGGING = 1;

LoadFunctionLibrary ("ReadDelimitedFiles");
LoadFunctionLibrary ("NJ");
LoadFunctionLibrary ("p_Distance_aa");

SetDialogPrompt 	("Please specify an amino-acid file:");

DataSet				inputAlignment     = ReadDataFile (PROMPT_FOR_FILE);

DataSetFilter		filteredData	   = CreateFilter (inputAlignment,1);

// check to see if this is a protein alignment 
GetDataInfo (_AncestralFilterChars, filteredData, "CHARACTERS");
assert (Columns (_AncestralFilterChars) == 20, "Expected a protein alignment");

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
_characterDimension 	= Columns (_AncestralFilterChars);


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

_ids = {_marginalAncestorsFilter.species};
_output = {_marginalAncestorsFilter.species,(_marginalAncestorsFilter.sites*_characterDimension)};

// _outputCSV = ""; _outputCSV * 8192; 
/* _outputCSV * ("Sequence,[Site:["+Join(",",_AncestralFilterChars)+"]]\n"); */

// _outputCSV * ("{\"order\":[\""+Join("\",\"",_AncestralFilterChars)+"\"],\n\"data\":[\n");

_idx_3 = 0;
for (_idx_1 = 0; _idx_1 < _marginalAncestorsFilter.species; _idx_1 = _idx_1 + 1)
{
    if (_idx_1 > 0)
    {
        _outputCSV * (",\n");
    }

    _ids[_idx_1] = newNamesToOldNames[_AncestralNodeNames[_idx_1]];

    // _outputCSV * ("{\"id\":\"" + newNamesToOldNames[_AncestralNodeNames[_idx_1]] + "\",\"values\":[");

	for (_idx_2 = 0; _idx_2 < _marginalAncestorsFilter.sites; _idx_2 = _idx_2 + 1)
	{
        /* trailing commas */
        if (_idx_2 > 0)
        {
            _outputCSV * (",");
        }

        // _outputCSV * ("[");
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
		
            _output[_idx_1][((_characterDimension * _idx_2) + _idx_4)] = _thisCharacter;

            // _outputCSV * (""+_thisCharacter);
			
            if (_thisCharacter > _maxValue)
			{
				_maxValue = _thisCharacter;
				_maxIndex = _idx_4;
			}
		}
		
		_idx_3 = _idx_3 + 1;
        
        // _outputCSV * ("]");
	}
    
    // _outputCSV * ("]}");
}

// _outputCSV * ("\n]\n}\n");
// _outputCSV * 0;
// fprintf ("/dev/null", CLEAR_FILE, _outputCSV);

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
    if (key == "data")
    {
        return _output;
    }
    return "_THyPhy_NOT_HANDLED_";
}
