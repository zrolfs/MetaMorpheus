using MassSpectrometry;
using SharpLearning.Common.Interfaces;
using SharpLearning.Containers.Matrices;
using SharpLearning.CrossValidation.TrainingTestSplitters;
using SharpLearning.GradientBoost.Learners;
using SharpLearning.GradientBoost.Models;
using SharpLearning.Metrics.Regression;
using SharpLearning.Optimization;
using SharpLearning.RandomForest.Learners;
using SharpLearning.RandomForest.Models;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MzLibUtil;

namespace EngineLayer.Aggregation
{
    public class AggregationEngine : MetaMorpheusEngine
    {
        private readonly double MaxRetentionTimeDifferenceAllowedInMinutes;
        private readonly double MinCosineScoreAllowed;
        private readonly int NumberOfMS1SpectraToAverage;
        private readonly Ms2ScanWithSpecificMass[] MS2Scans;

        private readonly MsDataFile OriginalFile;
        public MsDataFile AggregatedDataFile { get; private set; }
        public Tolerance SuggestedPrecursorTolerance { get; private set; }
        public Tolerance SuggestedProductTolerance { get; private set; }


        public AggregationEngine(MsDataFile originalFile, Ms2ScanWithSpecificMass[] ms2scans, CommonParameters commonParameters, List<string> nestedIds, double maxRetentionTimeDifferenceAllowedInMinutes, double minCosineScoreAllowed, int numberOfMS1SpectraToAverage) : base(commonParameters, nestedIds)
        {
            OriginalFile = originalFile;
            MS2Scans = ms2scans;
            MaxRetentionTimeDifferenceAllowedInMinutes = maxRetentionTimeDifferenceAllowedInMinutes;
            MinCosineScoreAllowed = minCosineScoreAllowed;
            NumberOfMS1SpectraToAverage = numberOfMS1SpectraToAverage;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Identifying MS2 groups");

            //need to group together which scans to compare
            //index ms2scans by precursor masses
            //Could change this to the max MS1 range observed
            int maxMass = 30000;
            int binsPerDalton = 1000;
            List<int>[] massIndex = new List<int>[maxMass * binsPerDalton]; //scan number possesing a precursor mass of the index
            for (int i = 0; i < MS2Scans.Length; i++)
            {
                int mass = (int)Math.Round(MS2Scans[i].PrecursorMonoisotopicPeakMz * binsPerDalton);
                if (massIndex[mass] == null)
                {
                    massIndex[mass] = new List<int> { i };
                }
                else
                {
                    massIndex[mass].Add(i);
                }
            }
            //somewhat tricky. We want to group a scan with all the other scans, but there's a chance that the tolerance falls outside for some but not for others.
            //we can just expand the range for tolerance whenever another is added.

            //want to group all scans to compare and later subgroup those where the comparison is above a threshold
            bool[] seen = new bool[MS2Scans.Length]; //don't double group stuff
            List<List<Ms2ScanWithSpecificMass>> groups = new List<List<Ms2ScanWithSpecificMass>>();
            for (int i = 0; i < MS2Scans.Length; i++)
            {
                //get indexes to compare
                if (!seen[i])
                {
                    seen[i] = true; //we've seen it, so don't use it again
                    var scan = MS2Scans[i]; //get the scan
                    int obsFragmentFloorMass = (int)Math.Floor((commonParameters.PrecursorMassTolerance.GetMinimumValue(scan.PrecursorMonoisotopicPeakMz)) * binsPerDalton);
                    int obsFragmentCeilingMass = (int)Math.Ceiling((commonParameters.PrecursorMassTolerance.GetMaximumValue(scan.PrecursorMonoisotopicPeakMz)) * binsPerDalton);
                    int minBinObs = -1; //save outer bounds so we can expand tolerances if needed
                    int maxBinObs = -1;
                    List<Ms2ScanWithSpecificMass> groupToAdd = new List<Ms2ScanWithSpecificMass>(); //current group
                    //foreach mz in range, expand range if necessary
                    for (int bin = obsFragmentFloorMass; bin <= obsFragmentCeilingMass; bin++) //go through bins and add to the group
                    {
                        List<int> scans = massIndex[bin];
                        if (scans != null) //FIXME: It's still possible to group things twice since large mz ppm can hit a smaller mz where the small ppm wouldn't
                        {
                            if (minBinObs == -1)
                            {
                                minBinObs = bin;
                            }
                            maxBinObs = bin;
                            foreach (int scanIndex in scans)
                            {
                                seen[scanIndex] = true;
                                groupToAdd.Add(MS2Scans[scanIndex]);
                            }
                        }
                    }

                    //get lower bound if it's expanded
                    int decreasingValues = obsFragmentFloorMass - 1;
                    int minimumValue = (int)Math.Floor(commonParameters.PrecursorMassTolerance.GetMinimumValue(minBinObs));
                    while (decreasingValues > minimumValue)
                    {
                        List<int> scans = massIndex[decreasingValues];
                        if (scans != null)
                        {
                            minimumValue = (int)Math.Floor(commonParameters.PrecursorMassTolerance.GetMinimumValue(decreasingValues));

                            foreach (int scanIndex in scans)
                            {
                                seen[scanIndex] = true;
                                groupToAdd.Add(MS2Scans[scanIndex]);
                            }
                        }
                        decreasingValues--;
                    }

                    //get upper bound if it's expanded
                    int increasingValues = obsFragmentCeilingMass + 1;
                    int maximumValue = (int)Math.Ceiling(commonParameters.PrecursorMassTolerance.GetMaximumValue(maxBinObs));
                    while (increasingValues < maximumValue)
                    {
                        List<int> scans = massIndex[increasingValues];
                        if (scans != null)
                        {
                            maximumValue = (int)Math.Ceiling(commonParameters.PrecursorMassTolerance.GetMaximumValue(increasingValues));

                            foreach (int scanIndex in scans)
                            {
                                seen[scanIndex] = true;
                                groupToAdd.Add(MS2Scans[scanIndex]);
                            }
                        }
                        increasingValues++;
                    }
                    groups.Add(groupToAdd);
                }
            }

            //Now that we've separated groups by mass, let's try to separate based on retention time
            List<List<Ms2ScanWithSpecificMass>> retGroups = new List<List<Ms2ScanWithSpecificMass>>();
            foreach (List<Ms2ScanWithSpecificMass> group in groups) //go over all the previously made groups
            {
                List<List<Ms2ScanWithSpecificMass>> subGroups = new List<List<Ms2ScanWithSpecificMass>>(); //local subgroups go here, needed so you don't regroup previous classifications
                foreach (Ms2ScanWithSpecificMass scan in group) //iterate through each scan in the previous group
                {
                    bool foundSpot = false;
                    foreach (List<Ms2ScanWithSpecificMass> subGroup in subGroups) //see if the scan fits somewhere, if not make a new subgroup
                    {
                        Ms2ScanWithSpecificMass scanInSubGroup = subGroup.Last(); //only need to check the lawst
                        {
                            if (scan.RetentionTime > scanInSubGroup.RetentionTime - MaxRetentionTimeDifferenceAllowedInMinutes
                                && scan.RetentionTime < scanInSubGroup.RetentionTime + MaxRetentionTimeDifferenceAllowedInMinutes) //if a match, add it
                            {
                                subGroup.Add(scan);
                                foundSpot = true;
                                break;
                            }
                        }
                        if (foundSpot)
                        {
                            break;
                        }
                    }
                    if (!foundSpot) //if not a match, create a new category
                    {
                        subGroups.Add(new List<Ms2ScanWithSpecificMass> { scan });
                    }
                }
                retGroups.AddRange(subGroups);
            }
            groups = retGroups; //save reassignments

            //Now let's separate based on cosine score
            List<List<Ms2ScanWithSpecificMass>> scoredGroups = new List<List<Ms2ScanWithSpecificMass>>();
            foreach (List<Ms2ScanWithSpecificMass> group in groups) //go over all the previously made groups
            {
                List<List<Ms2ScanWithSpecificMass>> subGroups = new List<List<Ms2ScanWithSpecificMass>>(); //local subgroups go here, needed so you don't regroup previous classifications
                foreach (Ms2ScanWithSpecificMass scan in group) //iterate through each scan in the previous group
                {
                    bool foundSpot = false;
                    foreach (List<Ms2ScanWithSpecificMass> subGroup in subGroups) //see if the scan fits somewhere, if not make a new subgroup
                    {
                        foreach (Ms2ScanWithSpecificMass scanInSubGroup in subGroup) //iterate through each member of each previously found group
                        {
                            if (MinCosineScoreAllowed <= CosineScore(scan.TheScan.MassSpectrum, scanInSubGroup.TheScan.MassSpectrum, commonParameters.ProductMassTolerance)) //if a match, add it
                            {
                                subGroup.Add(scan);
                                foundSpot = true;
                                break;
                            }
                        }
                        if (foundSpot)
                        {
                            break;
                        }
                    }
                    if (!foundSpot) //if not a match, create a new category
                    {
                        subGroups.Add(new List<Ms2ScanWithSpecificMass> { scan });
                    }
                }
                scoredGroups.AddRange(subGroups);
            }
            groups = scoredGroups; //save

            AggregatedDataFile = new MsDataFile(OriginalFile.GetAllScansList().ToArray(), OriginalFile.SourceFile);
            return new MetaMorpheusEngineResults(this);
        }


        public double CosineScore(MzSpectrum scan1, MzSpectrum scan2, Tolerance tolerance)
        {
            double[] mz1 = scan1.XArray;
            double[] intensity1 = scan1.YArray;
            double[] mz2 = scan2.XArray;
            double[] intensity2 = scan2.YArray;
            return CosineScore(mz1, mz2, intensity1, intensity2, tolerance);
        }

        public double CosineScore(MsDataScan scan1, MsDataScan scan2, Tolerance tolerance)
        {
            double[] mz1 = scan1.MassSpectrum.XArray;
            double[] intensity1 = scan1.MassSpectrum.YArray;
            double[] mz2 = scan2.MassSpectrum.XArray;
            double[] intensity2 = scan2.MassSpectrum.YArray;
            return CosineScore(mz1, mz2, intensity1, intensity2, tolerance);
        }

        public double CosineScore(double[] mz1, double[] mz2, double[] intensity1, double[] intensity2, Tolerance tolerance)
        {
            //convert spectra to vectors
            List<double> vector1 = new List<double>();
            List<double> vector2 = new List<double>();
            int i = 0;
            int j = 0;

            while (i != mz1.Length && j != mz2.Length)
            {
                double one = mz1[i];
                double two = mz2[j];
                if (tolerance.Within(one, two))
                {
                    vector1.Add(intensity1[i]);
                    vector2.Add(intensity2[j]);
                    i++;
                    j++;
                }
                else if (one > two)
                {
                    vector1.Add(0);
                    vector2.Add(intensity2[j]);
                    j++;
                }
                else //two>one
                {
                    vector1.Add(intensity1[i]);
                    vector2.Add(0);
                    i++;
                }
            }
            //wrap up leftover peaks
            for (; i < mz1.Length; i++)
            {
                vector1.Add(intensity1[i]);
                vector2.Add(0);
            }
            for (; j < mz2.Length; j++)
            {
                vector1.Add(0);
                vector2.Add(intensity2[j]);
            }

            //cosine score of vectors
            //numerator
            double numerator = 0;
            for (i = 0; i < vector1.Count; i++)
            {
                numerator += vector1[i] * vector2[i];
            }

            //denominator
            double denominator = Math.Sqrt(vector1.Sum(x => x * x)) * Math.Sqrt(vector2.Sum(x => x * x));

            //calculate cosine score
            return Math.Round(numerator / denominator * 1000) / 1000;
        }

    }
}