using Chemistry;
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
        private const int MaxNodesAtOnce = 100;
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
            char[] alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".ToCharArray();
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
            int[] masses = allResidues.Select(x => (int)Math.Round(x.MonoisotopicMass * FragmentBinsPerDalton)).ToArray();

            //record the minimum mass difference to search for
            MinimumMassDifference = (int)Math.Floor(commonParameters.ProductMassTolerance.GetMinimumValue(masses[0]));


            List<int> indexList = new List<int>(); //indexes for each "amino acid" in the current combination
            string sequenceTag = ""; //the string of combined residues for the current combination
            int lastAminoAcidIndex = aminoAcids.Length - 1; //get the last index for combinatorics
            int smallestMass = masses[0]; //smallest mass, which is different from the smallest mass difference by some small ppm
            int mass = 0; //total mass

            HashSet<string> sequencesAddedDEBUG = new HashSet<string>();

            do
            {
                //Determine the longest length possible for the current indexes
                int changingMass = smallestMass; //need to reset here as we try to add on the smallest mass
                mass += changingMass;
                while (mass < MaximumMassDifference)
                {
                    sequenceTag += aminoAcids[0];
                    indexList.Add(0);
                    int minMass = (int)Math.Floor(commonParameters.ProductMassTolerance.GetMinimumValue(mass));
                    int maxMass = Math.Min((int)Math.Ceiling(commonParameters.ProductMassTolerance.GetMaximumValue(mass)), MaximumMassDifference);

                    sequencesAddedDEBUG.Add(sequenceTag);
                    for (int currentMass = minMass; currentMass <= maxMass; currentMass++)
                    {
                        if (AminoAcidCombinationDictionary[currentMass] == null)
                        {
                            AminoAcidCombinationDictionary[currentMass] = new List<string> { sequenceTag };
                        }
                        else
                        {
                            AminoAcidCombinationDictionary[currentMass].Add(sequenceTag);
                        }
                    }
                    mass += changingMass;
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
                    if (mass < MaximumMassDifference) //if the mass is lower than the maximum cutoff
                    {
                        sequenceTag += aminoAcids[currentLastIndex];

                        int minMass = (int)Math.Floor(commonParameters.ProductMassTolerance.GetMinimumValue(mass));
                        int maxMass = Math.Min((int)Math.Ceiling(commonParameters.ProductMassTolerance.GetMaximumValue(mass)), MaximumMassDifference);

                        for (int currentMass = minMass; currentMass <= maxMass; currentMass++)
                        {
                            if (AminoAcidCombinationDictionary[currentMass] == null)
                            {
                                AminoAcidCombinationDictionary[currentMass] = new List<string> { sequenceTag };
                            }
                            else
                            {
                                AminoAcidCombinationDictionary[currentMass].Add(sequenceTag);
                            }
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
                    //check if we're at a valid mass
                    if (mass < MaximumMassDifference)
                    {
                        sequenceTag += aminoAcids[index];
                        int minMass = (int)Math.Floor(commonParameters.ProductMassTolerance.GetMinimumValue(mass));
                        int maxMass = Math.Min((int)Math.Ceiling(commonParameters.ProductMassTolerance.GetMaximumValue(mass)), MaximumMassDifference);
                        sequencesAddedDEBUG.Add(sequenceTag);
                        for (int currentMass = minMass; currentMass <= maxMass; currentMass++)
                        {
                            if (AminoAcidCombinationDictionary[currentMass] == null)
                            {
                                AminoAcidCombinationDictionary[currentMass] = new List<string> { sequenceTag };
                            }
                            else
                            {
                                AminoAcidCombinationDictionary[currentMass].Add(sequenceTag);
                            }
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

                    //Get spectrum
                    Ms2ScanWithSpecificMass spectrum = ArrayOfSortedMS2Scans[i];
                    //get masses
                    int precursorMass = (int)Math.Round(spectrum.PrecursorMass * FragmentBinsPerDalton);
                    int[] experimentalMasses = spectrum.ExperimentalFragments.Select(x => (int)Math.Round(x.monoisotopicMass * FragmentBinsPerDalton)).Where(x=>x<precursorMass-productShift).ToArray();

                    List<int> complimentaryMasses = experimentalMasses.ToList();
                    //addCompIons
                    if (complementaryIonConversionDictionary.TryGetValue(commonParameters.DissociationType, out double protonMassShift)) //TODO: this is broken for EThcD because that method needs two conversions
                    {
                        int shift = precursorMass + (int)Math.Round(Chemistry.ClassExtensions.ToMass(protonMassShift, 1) * FragmentBinsPerDalton);
                        for (int index = 0; index < experimentalMasses.Length; index++)
                        {
                            complimentaryMasses.Add(shift - experimentalMasses[index]);
                        }
                    }
                    else
                    {
                        throw new NotImplementedException();
                    }
                    experimentalMasses = complimentaryMasses.Select(x => x - productShift).OrderBy(x => x).ToArray(); //subtract the product mass to search for and sort.

                    //we need to identify "seed" masses, or we're screwed from the beginning.
                    List<DeNovoTree> temporaryNodeTrees = new List<DeNovoTree>();
                    int maxSeedIndex = BinarySearchBin(experimentalMasses, MaximumMassDifference, false);

                    for (int index = 0; index <= maxSeedIndex; index++)
                    {
                        int experimentalMass = experimentalMasses[index];
                        List<string> combinations = AminoAcidCombinationDictionary[experimentalMass];
                        if (combinations != null)
                        {
                            //nodeTrees.Add(new DeNovoTree(combinations, experimentalMass + productShift, productShift));
                            temporaryNodeTrees.Add(new DeNovoTree(combinations, experimentalMass));
                        }
                    }

                    //We now, hopefully, have lots of seed masses to build off of. But which are correct? We scream, for we do not know.
                    //Grow off of each seed in a tree structure. If we can't complete the tree to obtain the precursor mass, chop the tree down.
                    bool doneGrowing = true;
                    List<DeNovoTree> finalForest = new List<DeNovoTree>();
                    do
                    {
                        doneGrowing = true;
                        List<DeNovoTree> nextForest = new List<DeNovoTree>();
                        foreach (DeNovoTree tree in temporaryNodeTrees)
                        {
                            //check we actually haven't...
                            int currentMass = tree.CurrentMass;
                            //subtract precursor from current
                            int finalMassDifference = precursorMass - productShift - currentMass;
                            if (finalMassDifference < MaximumMassDifference)
                            {
                                bool deadEnd = true; //check if the mass difference is even worth looking for new nodes
                                List<string> combinations = AminoAcidCombinationDictionary[finalMassDifference];
                                if (combinations != null)
                                {
                                    finalForest.Add(new DeNovoTree(tree.Nodes, combinations, finalMassDifference, precursorMass));
                                    deadEnd = false;
                                }
                                if (deadEnd) //we can't hit the precursor, so no reason to keep looking at this tree
                                {
                                    continue;
                                }
                            }

                            int minIndex = BinarySearchBin(experimentalMasses, MinimumMassDifference + currentMass, true);
                            int maxIndex = BinarySearchBin(experimentalMasses, MaximumMassDifference + currentMass, false);
                            for (; minIndex <= maxIndex; minIndex++)
                            {
                                int experimentalMass = experimentalMasses[minIndex];
                                int experimentalMassDifference = experimentalMass - currentMass;
                                List<string> combinations = AminoAcidCombinationDictionary[experimentalMassDifference];
                                if (combinations != null)
                                {
                                    nextForest.Add(new DeNovoTree(tree.Nodes, combinations, experimentalMassDifference, experimentalMass));
                                    doneGrowing = false;
                                }
                            }

                            //Just need to keep it efficient. All the good ones are near the beginning anyway...
                            if (nextForest.Count > MaxNodesAtOnce)
                            {
                                break;
                            }
                        }
                        temporaryNodeTrees.Clear();
                        temporaryNodeTrees = nextForest;
                    } while (!doneGrowing);

                    //if we found anything at all
                    if (finalForest.Count != 0)
                    {
                        HashSet<string> deNovoSequences = new HashSet<string>();
                        List<double> deNovoScores = new List<double>();

                        //we need to figure out the probabilities for these trees.
                        foreach (DeNovoTree tree in finalForest)
                        {
                            tree.CalculateNumberOfPossibleSequences();
                        }

                        finalForest = finalForest.OrderBy(x => x.NumberOfPossibleSequences).ToList(); //we now have the number of possible sequences for each tree

                        //Start with 100% probability
                        //Divide by the number of trees (three trees each have 33% probability)
                        //Divide the remaining probability (33%) within each tree based on the amount of ambiguity in said tree

                        int currentTreeIndex = 0;
                        DeNovoTree currentTree = finalForest[currentTreeIndex];
                        double score = (100d) / (finalForest.Count * currentTree.NumberOfPossibleSequences); //THIS IS THE SCORE


                        //grab the nodes from the highest scoring tree to determine the sequences
                        int[] maxCombinationIndexes = new int[currentTree.Nodes.Count];
                        for (int index = 0; index < maxCombinationIndexes.Length; index++)
                        {
                            maxCombinationIndexes[index] = currentTree.Nodes[index].Count;
                        }
                        //set the index for each node to track the combinations
                        int[] currentCombinationIndexes = new int[maxCombinationIndexes.Length];

                        int currentNodeIteration = 0;
                        int psmIndex = i * 5;
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
                                    currentNodeIteration=0; //reset
                                    currentTree = finalForest[currentTreeIndex];
                                    score = (100d) / (finalForest.Count * currentTree.NumberOfPossibleSequences); //THIS IS THE SCORE

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
                                if (score > 0)
                                {
                                    PeptideWithSetModifications peptide = new PeptideWithSetModifications(sequence, new Dictionary<string, Modification>(), digestionParams: commonParameters.DigestionParams);

                                    List<Product> peptideTheorProducts = peptide.Fragment(commonParameters.DissociationType, commonParameters.DigestionParams.FragmentationTerminus).ToList();
                                    List<MatchedFragmentIon> matchedIons = MatchFragmentIons(spectrum, peptideTheorProducts, commonParameters);
                                    PeptideSpectralMatch psm = new PeptideSpectralMatch(peptide, 0, score, i, spectrum, commonParameters.DigestionParams, matchedIons);
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
                try
                {
                    psm.ResolveAllAmbiguities();
                }
                catch
                {

                }
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
                while (experimentalMasses[m] < massOfInterest && m != experimentalMasses.Length - 1)
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
    }
}
