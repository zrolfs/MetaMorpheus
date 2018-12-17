using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.DeNovoSequencing
{
    public class DeNovoSequencingEngine : MetaMorpheusEngine
    {
        protected const int FragmentBinsPerDalton = 1000;
        private const int MaximumMassDifference = 500 * FragmentBinsPerDalton;
        private const int MaxNodesAtOnce = 1000;
        private int MinimumMassDifference = 0;
        private readonly List<Modification> FixedModifications;
        private readonly List<Modification> VariableModifications;
        private readonly PeptideSpectralMatch[] PeptideSpectralMatches;
        private readonly Ms2ScanWithSpecificMass[] ArrayOfSortedMS2Scans;
        private List<string>[] AminoAcidCombinationDictionary;
        private readonly bool NTerminus;

        public DeNovoSequencingEngine(PeptideSpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans, List<Modification> variableModifications, List<Modification> fixedModifications, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            PeptideSpectralMatches = globalPsms;
            ArrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
            VariableModifications = variableModifications;
            FixedModifications = fixedModifications;
            NTerminus = commonParameters.DigestionParams.FragmentationTerminus == FragmentationTerminus.N;
        }


        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Preparing de novo index...");
            //get a list of all mass combinations as a dictionary (if takes a long time, may do indexing)
            //Get list of amino acid residues (TODO: mzlib public get)

            AminoAcidCombinationDictionary = new List<string>[MaximumMassDifference + 1];

            //DEBUG
            char[] alphabet = "ABCDEFGHJKLMNOPQRSTUVWXYZ".ToCharArray();
            //char[] alphabet = "DKLPTVY".ToCharArray();
            List<Residue> allResidues = new List<Residue>();
            foreach (char AA in alphabet)
            {
                if (Residue.TryGetResidue(AA, out Residue residue))
                {
                    allResidues.Add(residue);
                }
            }


            //TODO: Modifications
            //foreach (Modification fixedMod in FixedModifications)
            //{
            //    int modificationMass = (int)Math.Round(fixedMod.MonoisotopicMass.Value * FragmentBinsPerDalton);
            //    char AA = fixedMod.mot
            //}

            //convert all fixed mod residues to modified residues

            //add additional residues for variable mods

            //get all combinations below the specified maximum mass and populate the dictionary.

            //Order by mass, low to high
            allResidues = allResidues.OrderBy(x => x.MonoisotopicMass).ToList();
            char[] aminoAcids = allResidues.Select(x => x.Letter).ToArray();
            //get masses
            //int[] masses = allResidues.Select(x => (int)Math.Round(x.MonoisotopicMass * FragmentBinsPerDalton)).ToArray();
            double[] masses = allResidues.Select(x => x.MonoisotopicMass * FragmentBinsPerDalton).ToArray();

            //record the minimum mass difference to search for
            MinimumMassDifference = (int)Math.Floor(commonParameters.ProductMassTolerance.GetMinimumValue(masses[0]));


            List<int> indexList = new List<int>(); //indexes for each "amino acid" in the current combination
            string sequenceTag = ""; //the string of combined residues for the current combination
            int lastAminoAcidIndex = aminoAcids.Length - 1; //get the last index for combinatorics
            int smallestMass = (int)Math.Floor(masses[0]); //smallest mass, which is different from the smallest mass difference by some small ppm
            double mass = 0; //total mass
            int massInt = 0;

            do
            {
                //Determine the longest length possible for the current indexes
                double changingMass = masses[0]; //need to reset here as we try to add on the smallest mass
                mass += changingMass;
                massInt = (int)Math.Round(mass);
                while (massInt < MaximumMassDifference)
                {
                    sequenceTag += aminoAcids[0];
                    indexList.Add(0);

                    if (AminoAcidCombinationDictionary[massInt] == null)
                    {
                        AminoAcidCombinationDictionary[massInt] = new List<string> { sequenceTag };
                    }
                    else
                    {
                        AminoAcidCombinationDictionary[massInt].Add(sequenceTag);
                    }
                    mass += changingMass;
                    massInt = (int)Math.Round(mass);
                }

                //Wrap up the leftovers from the while loop
                mass -= (changingMass + masses[indexList[indexList.Count - 1]]); //undo last add that added an extra "amino acid" and the last amino acid so we can add a new one onto the tag
                sequenceTag = sequenceTag.Substring(0, sequenceTag.Length - 1); //clip off the c-terminal amino acid so we can add a different one

                //done determining the longest length. Most likely looks like GGGGGGGG


                //increase the mass of the terminal amino acid, so we'll get GGGGGGGA, GGGGGGGS, etc
                int currentLastIndex = indexList[indexList.Count - 1];
                while (currentLastIndex != lastAminoAcidIndex) //go until we cycle through all amino acids, or we get past the highest mass
                {
                    currentLastIndex++;
                    changingMass = masses[currentLastIndex];
                    mass += changingMass;
                    massInt = (int)Math.Round(mass);

                    if (massInt < MaximumMassDifference) //if the mass is lower than the maximum cutoff
                    {
                        sequenceTag += aminoAcids[currentLastIndex];

                        if (AminoAcidCombinationDictionary[massInt] == null)
                        {
                            AminoAcidCombinationDictionary[massInt] = new List<string> { sequenceTag };
                        }
                        else
                        {
                            AminoAcidCombinationDictionary[massInt].Add(sequenceTag);
                        }
                        sequenceTag = sequenceTag.Substring(0, sequenceTag.Length - 1); //clip off the c-terminal amino acid so we can add a different one
                    }
                    else //oops! we're too big now, we should leave
                    {
                        currentLastIndex = lastAminoAcidIndex; //get out of this loop, we're done here.
                    }
                    mass -= changingMass; //remove the mass we added so that we can add a new mass
                }

                //We cycled through the C-terminus, so now we need to change something else. Move in one residue and see if we can update that. If not, keep moving in residues
                do
                {
                    indexList.RemoveAt(indexList.Count - 1);
                    if (indexList.Count != 0)
                    {
                        mass -= masses[indexList[indexList.Count - 1]];
                        sequenceTag = sequenceTag.Substring(0, sequenceTag.Length - 1); //clip off the c-terminal amino acid so we can add a different one
                    }
                } while (indexList.Count != 0 && indexList[indexList.Count - 1] == lastAminoAcidIndex);

                //remove the terminal amino acid
                if (indexList.Count != 0)
                {
                    //++the new index
                    indexList[indexList.Count - 1]++;
                    int index = indexList[indexList.Count - 1];
                    mass += masses[index];
                    massInt = (int)Math.Round(mass);

                    //check if we're at a valid mass
                    if (massInt < MaximumMassDifference)
                    {
                        sequenceTag += aminoAcids[index];

                        if (AminoAcidCombinationDictionary[massInt] == null)
                        {
                            AminoAcidCombinationDictionary[massInt] = new List<string> { sequenceTag };
                        }
                        else
                        {
                            AminoAcidCombinationDictionary[massInt].Add(sequenceTag);
                        }
                    }
                    else
                    {
                        mass -= masses[index];
                        indexList[indexList.Count - 1] = lastAminoAcidIndex;
                    }
                }
            } while (indexList.Count != 0);

            //We have nice theoretical masses, but we need to get the actual fragment masses if we want to compare some of the lower mass experimental peaks.
            //We should only search for ions of a single fragmentation type and supplement the other with complementary ions.
            var ProductTypesToSearch = DissociationTypeCollection.ProductsFromDissociationType[commonParameters.DissociationType].Intersect(TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[commonParameters.DigestionParams.FragmentationTerminus]).ToList();
            List<int> productShifts = new List<int>();
            var productShift = (int)Math.Round(DissociationTypeCollection.GetMassShiftFromProductType(ProductTypesToSearch[0]) * FragmentBinsPerDalton);





            Status("Performing de novo sequencing...");

            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing de novo sequencing... ", nestedIds));

            Parallel.ForEach(Partitioner.Create(0, ArrayOfSortedMS2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, (range, loopState) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops)
                    {
                        loopState.Stop();
                        return;
                    }

                    //1. PREPARE THE EXPERIMENTAL SPECTRA

                    //Get spectrum
                    Ms2ScanWithSpecificMass spectrum = ArrayOfSortedMS2Scans[i];

                    //get masses
                    //int precursorMass = (int)Math.Round(spectrum.PrecursorMass * FragmentBinsPerDalton);
                    //int[] experimentalMasses = spectrum.ExperimentalFragments.Select(x => (int)Math.Round(x.monoisotopicMass * FragmentBinsPerDalton)).Where(x => x < precursorMass - productShift).ToArray();
                    //int[] experimentalIntensities = spectrum.ExperimentalFragments.Select(x=>x.totalIntensity)
                    //List<int> complimentaryMasses = experimentalMasses.ToList();
                    ////addCompIons
                    //if (complementaryIonConversionDictionary.TryGetValue(commonParameters.DissociationType, out double protonMassShift)) //TODO: this is broken for EThcD because that method needs two conversions
                    //{
                    //    int shift = precursorMass + (int)Math.Round(Chemistry.ClassExtensions.ToMass(protonMassShift, 1) * FragmentBinsPerDalton);
                    //    for (int index = 0; index < experimentalMasses.Length; index++)
                    //    {
                    //        complimentaryMasses.Add(shift - experimentalMasses[index]);
                    //    }
                    //}
                    //else
                    //{
                    //    throw new NotImplementedException();
                    //}
                    //experimentalMasses = complimentaryMasses.Select(x => x - productShift).OrderBy(x => x).ToArray(); //subtract the product mass to search for and sort.

                    int precursorMass = (int)Math.Round(spectrum.PrecursorMass * FragmentBinsPerDalton);
                    int[] experimentalMasses = spectrum.ExperimentalFragments.Select(x => (int)Math.Round(x.monoisotopicMass * FragmentBinsPerDalton)).Where(x => x < precursorMass - productShift).ToArray();
                    double[] experimentalIntensities = spectrum.ExperimentalFragments.Where(x => x.monoisotopicMass * FragmentBinsPerDalton < precursorMass - productShift).Select(x => x.totalIntensity).ToArray();

                    List<MzPeak> complementaryPeaks = new List<MzPeak>();

                    if (complementaryIonConversionDictionary.TryGetValue(commonParameters.DissociationType, out double protonMassShift)) //TODO: this is broken for EThcD because that method needs two conversions
                    {
                        double shift = precursorMass + (int)Math.Round(Chemistry.ClassExtensions.ToMass(protonMassShift, 1) * FragmentBinsPerDalton);
                        for (int index = 0; index < experimentalMasses.Length; index++)
                        {
                            int experimentalMass = experimentalMasses[index];
                            double experimentalIntensity = experimentalIntensities[index];
                            complementaryPeaks.Add(new MzPeak(experimentalMass, experimentalIntensity));
                            if (experimentalMass > 200 * FragmentBinsPerDalton) //limit the mass range so that immonium ions don't get doubled. Issue with K if K isn't on the terminus
                            {
                                complementaryPeaks.Add(new MzPeak(shift - experimentalMass, experimentalIntensity));
                            }
                        }
                    }
                    else
                    {
                        throw new NotImplementedException();
                    }
                    complementaryPeaks = complementaryPeaks.OrderBy(x => x.Mz).ToList();
                    experimentalMasses = complementaryPeaks.Select(x => (int)x.Mz - productShift).ToArray();
                    experimentalIntensities = complementaryPeaks.Select(x => x.Intensity).ToArray();


                    //2. IDENTIFY THE POSSIBLE BEGINNINGS OF EACH SEQUENCE

                    //we need to identify "seed" masses, or we're screwed from the beginning.
                    List<DeNovoTree> temporaryNodeTrees = new List<DeNovoTree>(); //save the seeds
                    int maxSeedIndex = BinarySearchBin(experimentalMasses, MaximumMassDifference, false); //get the index of the largest experimental mass to look for

                    int lowestMassSeen = 0; //record the lowest mass used to prevent double counting of the same mass shift (happens when multiple experimental peaks are close in mass)
                    for (int index = 0; index <= maxSeedIndex; index++) //for all masses in the range of interest
                    {
                        int experimentalMass = experimentalMasses[index];

                        //get the range of theoretical bins for this mass
                        int minMass = Math.Max((int)Math.Floor(commonParameters.ProductMassTolerance.GetMinimumValue(experimentalMass)), lowestMassSeen);
                        int maxMass = Math.Min((int)Math.Ceiling(commonParameters.ProductMassTolerance.GetMaximumValue(experimentalMass)), MaximumMassDifference);

                        //get combinations with the exact same mass
                        for (; minMass <= maxMass; minMass++)
                        {
                            List<string> combinations = AminoAcidCombinationDictionary[minMass];
                            if (combinations != null)
                            {
                                temporaryNodeTrees.Add(new DeNovoTree(combinations, minMass, experimentalIntensities[index])); //save the theoretical mass, since that is what we claim is there
                            }
                        }
                        lowestMassSeen = maxMass;
                    }


                    //3. EXPAND THE SEQUENCE NETWORKS TO OBTAIN THE PRECURSOR MASS

                    //We now, hopefully, have lots of seed masses to build off of. But which are correct? We scream, for we do not know.
                    //Grow off of each seed in a tree structure. If we can't complete the tree to obtain the precursor mass, remove it.
                    List<DeNovoTree> finalForest = new List<DeNovoTree>();

                    //get the range of theoretical bins for this mass
                    int minPrecursorMass = (int)Math.Floor(commonParameters.PrecursorMassTolerance.GetMinimumValue(precursorMass));
                    int maxPrecursorMass = (int)Math.Ceiling(commonParameters.PrecursorMassTolerance.GetMaximumValue(precursorMass));

                    do
                    {
                        List<DeNovoTree> nextForest = new List<DeNovoTree>();
                        foreach (DeNovoTree tree in temporaryNodeTrees)
                        {
                            //check to see if we can complete the tree by reaching the precursor mass
                            int currentMass = tree.CurrentMass;

                            //get the range of theoretical bins for this mass
                            int minMass = minPrecursorMass - productShift - currentMass;
                            int maxMass = Math.Min(maxPrecursorMass - productShift - currentMass, MaximumMassDifference);

                            //subtract precursor from current
                            if (minMass < MaximumMassDifference)
                            {
                                bool deadEnd = true; //check if the mass difference is even worth looking for new nodes

                                for (; minMass <= maxMass; minMass++)
                                {
                                    List<string> combinations = AminoAcidCombinationDictionary[minMass];
                                    if (combinations != null)
                                    {
                                        finalForest.Add(new DeNovoTree(tree.Nodes, combinations, 0, tree.NumberOfPossibleSequences * combinations.Count, tree.IntensitySum)); //current mass "0" doesn't matter anymore
                                        deadEnd = false;
                                    }
                                }
                                if (deadEnd) //we can't hit the precursor, so no reason to keep looking at this tree
                                {
                                    continue;
                                }
                            }

                            //get indexes for experimental peaks
                            int minIndex = BinarySearchBin(experimentalMasses, MinimumMassDifference + currentMass, true);
                            int maxIndex = BinarySearchBin(experimentalMasses, MaximumMassDifference + currentMass, false);

                            lowestMassSeen = 0; //reset for the new tree

                            for (; minIndex <= maxIndex; minIndex++)
                            {
                                int experimentalMass = experimentalMasses[minIndex];

                                minMass = Math.Max((int)Math.Floor(commonParameters.ProductMassTolerance.GetMinimumValue(experimentalMass)) - currentMass, lowestMassSeen);
                                maxMass = Math.Min((int)Math.Ceiling(commonParameters.ProductMassTolerance.GetMaximumValue(experimentalMass)) - currentMass, MaximumMassDifference);

                                for (; minMass <= maxMass; minMass++)
                                {
                                    List<string> combinations = AminoAcidCombinationDictionary[minMass];
                                    if (combinations != null)
                                    {
                                        double debugInt = experimentalIntensities[minIndex];
                                        nextForest.Add(new DeNovoTree(tree.Nodes, combinations, currentMass + minMass, tree.NumberOfPossibleSequences * combinations.Count, tree.IntensitySum + experimentalIntensities[minIndex]));
                                    }
                                }
                                lowestMassSeen = maxMass;
                            }
                        }
                        temporaryNodeTrees.Clear();

                        //save the list of trees for the next iteration
                        //often, there are too many combinations generated to effectively interrgate each. We need to crop the list to reduce computational demand
                        if (nextForest.Count > MaxNodesAtOnce) //if there are too many combinations
                        {
                            //nextForest = nextForest.OrderBy(x => x.NumberOfPossibleSequences).ToList();
                            nextForest = nextForest.OrderBy(x => x.CurrentMass).ToList();
                            for (int treeIndex = 0; treeIndex < MaxNodesAtOnce; treeIndex++)
                            {
                                temporaryNodeTrees.Add(nextForest[treeIndex]);
                            }
                        }
                        else
                        {
                            temporaryNodeTrees = nextForest;
                        }
                    } while (temporaryNodeTrees.Count != 0);


                    //4. SCORE THE NETWORKS AND SAVE THE HIGHEST SCORING SEQUENCES

                    //need to consider if it matters that the nodes themselves are unambiguous if there's another tree with the same number of nodes. 
                    //Isn't number of nodes more important?
                    //But the number of combinations is important for the actual probability
                    //...Isn't it?
                    //Need to relate score with probability, and the relaion is unlikely to be linear.

                    //Nodes should always triumph
                    //After that, it's a toss up between all of the things
                    //Could do the number of (nodes^3) * the number of combinations for each divided by the sum of all of these?
                    //everything with the same number of nodes within a single scan will have the same score

                    //e-value could be useful?


                    //if we found anything at all
                    if (finalForest.Count != 0)
                    {
                        HashSet<string> deNovoSequences = new HashSet<string>();
                        List<double> deNovoScores = new List<double>();

                        //score each tree
                        AssignDeNovoScores(finalForest, spectrum.TotalIonCurrent);

                        //sort by score
                        finalForest = finalForest.OrderByDescending(x => x.Score).ToList();

                        int currentTreeIndex = 0;
                        DeNovoTree currentTree = finalForest[currentTreeIndex];

                        //grab the nodes from the highest scoring tree to determine the sequences
                        int[] maxCombinationIndexes = new int[currentTree.Nodes.Count];
                        for (int index = 0; index < maxCombinationIndexes.Length; index++)
                        {
                            maxCombinationIndexes[index] = currentTree.Nodes[index].Count;
                        }
                        //set the index for each node to track the combinations
                        int[] currentCombinationIndexes = new int[maxCombinationIndexes.Length];

                        int currentNodeIteration = 0;
                        int psmIndex = i * commonParameters.NumberOfSequencesPerPrecursor; //use to figure out where to place the psm in the psm array
                                                                                           //get as many sequences as desired, prioritizing the highest scoring sequences
                        for (int psmNumber = 0; psmNumber < commonParameters.NumberOfSequencesPerPrecursor; psmNumber++)
                        {
                            //check if we already recorded all of the combinations from the previous tree
                            if (currentNodeIteration == currentTree.NumberOfPossibleSequences)
                            {
                                //switch to the next most confident tree
                                currentTreeIndex++;
                                if (currentTreeIndex != finalForest.Count) //if there are still trees left to explore
                                {
                                    currentNodeIteration = 0; //reset
                                    currentTree = finalForest[currentTreeIndex];

                                    maxCombinationIndexes = new int[currentTree.Nodes.Count];
                                    for (int index = 0; index < maxCombinationIndexes.Length; index++)
                                    {
                                        maxCombinationIndexes[index] = currentTree.Nodes[index].Count;
                                    }
                                    currentCombinationIndexes = new int[maxCombinationIndexes.Length];
                                }
                                else
                                {
                                    break;
                                }
                            }


                            //get nodes and combinations, changing the last available combination each time
                            string sequence = "";
                            bool updateNeeded = true;

                            for (int index = 0; index < currentCombinationIndexes.Length; index++)
                            {
                                //get the current index
                                int currentCombinationIndex = currentCombinationIndexes[index];
                                //determine where to add the combination
                                if (NTerminus)
                                {
                                    sequence += currentTree.Nodes[index][currentCombinationIndex];
                                }
                                else
                                {
                                    sequence = currentTree.Nodes[index][currentCombinationIndex] + sequence;
                                }

                                //update the index if necessary
                                if (updateNeeded)
                                {
                                    if (currentCombinationIndex == maxCombinationIndexes[index] - 1)
                                    {
                                        currentCombinationIndexes[index] = 0;
                                    }
                                    else
                                    {
                                        currentCombinationIndexes[index]++;
                                        updateNeeded = false;
                                    }
                                }
                            }

                            //if it's already been added, ignore it
                            if (deNovoSequences.Contains(sequence))
                            {
                                psmNumber--;
                            }
                            else
                            {
                                if (currentTree.Score > 0)
                                {
                                    PeptideWithSetModifications peptide = new PeptideWithSetModifications(sequence, new Dictionary<string, Modification>(), digestionParams: commonParameters.DigestionParams);

                                    List<Product> peptideTheorProducts = peptide.Fragment(commonParameters.DissociationType, commonParameters.DigestionParams.FragmentationTerminus).ToList();
                                    List<MatchedFragmentIon> matchedIons = MatchFragmentIons(spectrum, peptideTheorProducts, commonParameters);
                                    PeptideSpectralMatch psm = new PeptideSpectralMatch(peptide, 0, currentTree.Score, i, spectrum, commonParameters.DigestionParams, matchedIons);
                                    PeptideSpectralMatches[psmIndex + psmNumber] = psm;
                                    deNovoSequences.Add(sequence);
                                }
                            }
                            currentNodeIteration++;
                        }
                    }
                    //reduce RAM efficiently
                    finalForest.Clear();
                    //Score based on the unconfident residues is more important. Rather than average local confidence, should be multiplying each local confidence (100, 100, 70 should be 70, not 80). If you also have a 100, 100, 30, then you'd have a final of 100 combined. 
                    //"N" can never be over 50%, because of "GG". "Q" can never be over 33%, because of "AG" and "GA". Interesting question on if "Gg" and "gG" are different, and should be 66% or 50%...

                    //would be interesting to post search adjust probabilities based on frequency of assignments                 

                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / ArrayOfSortedMS2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing de novo sequencing... ", nestedIds));
                    }
                }
            });

            foreach (PeptideSpectralMatch psm in PeptideSpectralMatches.Where(p => p != null))
            {
                psm.ResolveAllAmbiguities();
            }

            return new MetaMorpheusEngineResults(this);
        }

        private static int BinarySearchBin(int[] experimentalMasses, int massOfInterest, bool lowRange)
        {
            int m = 0;
            int l = 0;
            int h = experimentalMasses.Length - 1;
            if (h == -1)
            {
                return h;
            }

            // binary search in the fragment bin for precursor mass
            while (l <= h)
            {
                m = (h + l) / 2;

                if (h - l < 2)
                    break;
                if (experimentalMasses[m] < massOfInterest)
                    l = m + 1;
                else
                    h = m - 1;
            }

            if (lowRange)
            {
                while (m != 0 && experimentalMasses[m - 1] > massOfInterest)
                {
                    m--;
                }
                while (m != experimentalMasses.Length && experimentalMasses[m] < massOfInterest) //out of index is an option to be output, since it is inclusive
                {
                    m++;
                }
            }
            else
            {
                while (m != experimentalMasses.Length - 1 && experimentalMasses[m + 1] < massOfInterest)
                {
                    m++;
                }
                while (m != -1 && experimentalMasses[m] > massOfInterest) //-1 is an option to be output, since it is inclusive
                {
                    m--;
                }
            }
            return m;
        }

        private void AssignDeNovoScores(List<DeNovoTree> forest, double TIC)
        {
            int maxNodes = forest.Max(x => x.Nodes.Count);
            int minNodes = forest.Min(x => x.Nodes.Count);

            long[] combinationTotalsPerNodeCount = new long[maxNodes + 1];
            for (int node = minNodes; node <= maxNodes; node++)
            {
                List<DeNovoTree> treesAtThisNodeCount = forest.Where(x => x.Nodes.Count == node).ToList();
                long sum = 0;
                foreach (DeNovoTree tree in treesAtThisNodeCount)
                {
                    sum += tree.NumberOfPossibleSequences;
                }
                combinationTotalsPerNodeCount[node] = sum;
            }

            long totalProbabilitySum = 0;
            for (int node = minNodes; node <= maxNodes; node++)
            {
                totalProbabilitySum += (int)Math.Floor(Math.Pow(1d * node / maxNodes, 2) * combinationTotalsPerNodeCount[node]);
            }

            foreach (DeNovoTree tree in forest)
            {
                tree.SetDeNovoScore(totalProbabilitySum, TIC, maxNodes);
            }
        }
    }
}
