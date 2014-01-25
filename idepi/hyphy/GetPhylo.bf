MESSAGE_LOGGING = 0;

LoadFunctionLibrary ("ReadDelimitedFiles.bf");
LoadFunctionLibrary ("NJ.bf");
LoadFunctionLibrary ("p_Distance_aa");

SetDialogPrompt ("Please specify an amino-acid file:");

DataSet       inputAlignment = ReadDataFile (_inputFile);
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

LoadFunctionLibrary ("HIVbetween+F.mdl", {"00":"Fixed Rates"});
Tree aTree = treeString;

newTreeString = Format(aTree, 1, 1);

LikelihoodFunction lf = (filteredData, aTree);
Optimize (res, lf);

fprintf (stdout, lf);

fprintf (stdout, "\nReconstructing ancestors (this could also take a bit of time)...\n");

DataSet       _marginalAncestors       = ReconstructAncestors (lf);
DataSetFilter _marginalAncestorsFilter = CreateFilter (_marginalAncestors, 1);

_output = "";
_output * ((_marginalAncestorsFilter.sites + 50) * (_marginalAncestorsFilter.species + filteredData.species));

for (_idx_1 = 0; _idx_1 < filteredData.species; _idx_1 += 1) {
    GetDataInfo(_seq, filteredData, _idx_1);
    GetString(_name, filteredData, _idx_1);
    _output * (">" + _name + "\n" + _seq + "\n");
}

for (_idx_1 = 0; _idx_1 < _marginalAncestorsFilter.species; _idx_1 += 1) {
    GetDataInfo(_seq, _marginalAncestorsFilter, _idx_1);
    GetString(_name, _marginalAncestorsFilter, _idx_1);
    _output * (">" + _name + "\n" + _seq + "\n");
}

_output * 0;

function _THyPhyAskFor(key)
{
    if (key == "ancestors")
    {
        return _output;
    }
    if (key == "tree")
    {
        return newTreeString;
    }
    return "_THyPhy_NOT_HANDLED_";
}
