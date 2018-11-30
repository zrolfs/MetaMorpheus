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
        private readonly MassDiffAcceptor SearchMode;
        private readonly List<Modification> FixedModifications;
        private readonly List<Modification> VariableModifications;
        private readonly PeptideSpectralMatch[] PeptideSpectralMatches;
        private readonly Ms2ScanWithSpecificMass[] ArrayOfSortedMS2Scans;
        private static readonly Dictionary<char, int> ResidueDictionary;


        static DeNovoSequencingEngine() //populate residues and masses
        {
            char[] alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".ToCharArray();
            ResidueDictionary = new Dictionary<char, int>();
            foreach (char AA in alphabet)
            {
                if(Residue.TryGetResidue(AA, out Residue residue))
                {
                    ResidueDictionary[AA] = (int)Math.Round(residue.MonoisotopicMass * FragmentBinsPerDalton); 
                }
            }
        }

        public DeNovoSequencingEngine(PeptideSpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans, List<Modification> variableModifications, List<Modification> fixedModifications, List<Protein> proteinList, MassDiffAcceptor searchMode, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            PeptideSpectralMatches = globalPsms;
            ArrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
            VariableModifications = variableModifications;
            FixedModifications = fixedModifications;
            SearchMode = searchMode;
        }


        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Preparing de novo index...");
            foreach(Modification fixedMod in FixedModifications)
            {
                int modificationMass = (int)Math.Round(fixedMod.MonoisotopicMass.Value * FragmentBinsPerDalton);
                char AA = fixedMod.mot
            }


            Status("Performing de novo sequencing...");

            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing de novo sequencing... ", nestedIds));

            int maxAminoAcidMass = 
            Residue[] residueIndex = new Residue[];

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
                    double[] masses = spectrum.ExperimentalFragments.Select(x => x.monoisotopicMass).ToArray(); //assume sorted low to high.

                    //Score based on the unconfident residues is more important. Rather than average local confidence, should be multiplying each local confidence (100, 100, 70 should be 70, not 80). If you also have a 100, 100, 30, then you'd have a final of 100 combined. 
                    //"N" can never be over 50%, because of "GG". "Q" can never be over 50%, because of "AG".


                    // digest each protein into peptides and search for each peptide in all spectra within precursor mass tolerance
                    foreach (PeptideWithSetModifications peptide in Proteins[i].Digest(commonParameters.DigestionParams, FixedModifications, VariableModifications))
                    {
                        List<Product> peptideTheorProducts = peptide.Fragment(commonParameters.DissociationType, commonParameters.DigestionParams.FragmentationTerminus).ToList();

                        foreach (ScanWithIndexAndNotchInfo scan in GetAcceptableScans(peptide.MonoisotopicMass, SearchMode))
                        {
                            List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan.TheScan, peptideTheorProducts, commonParameters);

                            double thisScore = CalculatePeptideScore(scan.TheScan.TheScan, matchedIons, 0);

                            bool meetsScoreCutoff = thisScore >= commonParameters.ScoreCutoff;

                            // this is thread-safe because even if the score improves from another thread writing to this PSM,
                            // the lock combined with AddOrReplace method will ensure thread safety
                            if (meetsScoreCutoff || commonParameters.CalculateEValue)
                            {
                                // valid hit (met the cutoff score); lock the scan to prevent other threads from accessing it
                                lock (myLocks[scan.ScanIndex])
                                {
                                    bool scoreImprovement = PeptideSpectralMatches[scan.ScanIndex] == null || (thisScore - PeptideSpectralMatches[scan.ScanIndex].RunnerUpScore) > -PeptideSpectralMatch.ToleranceForScoreDifferentiation;

                                    if (scoreImprovement)
                                    {
                                        if (PeptideSpectralMatches[scan.ScanIndex] == null)
                                        {
                                            PeptideSpectralMatches[scan.ScanIndex] = new PeptideSpectralMatch(peptide, scan.Notch, thisScore, scan.ScanIndex, scan.TheScan, commonParameters.DigestionParams, matchedIons);
                                        }
                                        else
                                        {
                                            PeptideSpectralMatches[scan.ScanIndex].AddOrReplace(peptide, thisScore, scan.Notch, commonParameters.ReportAllAmbiguity, matchedIons);
                                        }
                                    }

                                    if (commonParameters.CalculateEValue)
                                    {
                                        PeptideSpectralMatches[scan.ScanIndex].AllScores.Add(thisScore);
                                    }
                                }
                            }
                        }
                    }

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


    }
}
