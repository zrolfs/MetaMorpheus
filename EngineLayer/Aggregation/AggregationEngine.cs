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
        private readonly Ms2ScanWithSpecificMass[] MS2Scans;

        private readonly MsDataFile OriginalFile;
        public MsDataFile AggregatedDataFile { get; private set; }
        public Tolerance SuggestedPrecursorTolerance { get; private set; }
        public Tolerance SuggestedProductTolerance { get; private set; }


        public AggregationEngine(MsDataFile originalFile, Ms2ScanWithSpecificMass[] ms2scans, CommonParameters commonParameters, List<string> nestedIds, double maxRetentionTimeDifferenceAllowedInMinutes, double minCosineScoreAllowed) : base(commonParameters, nestedIds)
        {
            OriginalFile = originalFile;
            MS2Scans = ms2scans;
            MaxRetentionTimeDifferenceAllowedInMinutes = maxRetentionTimeDifferenceAllowedInMinutes;
            MinCosineScoreAllowed = minCosineScoreAllowed;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Averaging MS1 spectra");
            //we are ONLY averaging m/z and NOT intensity. Intensity averaging would convolute downstream quantification
            MsDataScan[] ms1scans = OriginalFile.GetAllScansList().Where(x => x.MsnOrder == 1).ToArray();
            //we have a set of peaks in the ms1 scan, and we'll cycle through until:
            //-we have two consecutive ms1 scans that do not contain a peak, 
            //-we reach the end of the file
            //-we're told to stop by numberOfMs1SpectraToAverage
            List<double>[] peaksThatWereFound = new List<double>[ms1scans.Length]; //index is the scan index that tracks the peaks (doubles) previously found to prevent duplicate comparisons
            double[][] ms1mzs = ms1scans.Select(x => x.MassSpectrum.XArray).ToArray();

            for (int index = 0; index < ms1mzs.Length; index++)
            {
                double[] currentMzs = ms1mzs[index]; //grab current ms1scan
                List<double>[] allPeaksFound = currentMzs.Select(x => new List<double> { x }).ToArray(); //convert each mz into a seed for the list.
                bool[] peaksRecentlyFound = Enumerable.Repeat(true, currentMzs.Length).ToArray(); //set default to true, since they're all in this spectrum.

                //start cycling
                for(int secondaryIndex = index+1; secondaryIndex<ms1mzs.Length; secondaryIndex++)
                {

                }
            }


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


            Status("Averaging MS2 spectra");
            //Each MS2 spectra is going to be assigned a new scan number that's placed at the earliest occurance of the group.
            List<Ms2ScanWithSpecificMass> syntheticSpectra = new List<Ms2ScanWithSpecificMass>();
            List<int> minimumScanNumbersForSyntheticSpectra = new List<int>();
            foreach (List<Ms2ScanWithSpecificMass> group in groups)
            {
                if (group.Count == 1) //if nothing to aggregate, just save it
                {
                    syntheticSpectra.Add(group[0]);
                    minimumScanNumbersForSyntheticSpectra.Add(group[0].OneBasedScanNumber);
                }
                else
                {
                    minimumScanNumbersForSyntheticSpectra.Add(group.Select(x => x.OneBasedScanNumber).Min());

                    List<double> synethicMZs = new List<double>();
                    List<double> syntheticIntensities = new List<double>();

                    //get masses
                    double[][] groupedMzs = group.Select(x => x.TheScan.MassSpectrum.XArray).ToArray();
                    double[][] groupedIntensities = group.Select(x => x.TheScan.MassSpectrum.YArray).ToArray();


                    //require a peak to be present in at least half of all grouped spectra for it to be saved

                }
            }

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