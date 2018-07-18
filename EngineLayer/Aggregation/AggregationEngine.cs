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
using MathNet.Numerics.Statistics;
using System.Linq;

namespace EngineLayer.Aggregation
{
    public class AggregationEngine : MetaMorpheusEngine
    {
        private readonly double MaxRetentionTimeDifferenceAllowedInMinutes;
        private readonly double MinCosineScoreAllowed;
        private readonly Tolerance PrecursorTolerance;
        private readonly Tolerance ProductTolerance;
        private MsDataFile originalFile;
        private readonly string OriginalFilePath;
        private const int numberOfStrikesBeforeOut = 2;

        public MsDataFile AggregatedDataFile { get; private set; }
        public List<double> elutionProfileWidthsInScans { get; private set; }
        public Tolerance SuggestedPrecursorTolerance { get; private set; }
        public Tolerance SuggestedProductTolerance { get; private set; }


        public AggregationEngine(MsDataFile originalFile, string originalFilePath, CommonParameters commonParameters, List<string> nestedIds, double maxRetentionTimeDifferenceAllowedInMinutes, double minCosineScoreAllowed) : base(commonParameters, nestedIds)
        {
            this.originalFile = originalFile;
            MaxRetentionTimeDifferenceAllowedInMinutes = maxRetentionTimeDifferenceAllowedInMinutes;
            MinCosineScoreAllowed = minCosineScoreAllowed;
            PrecursorTolerance = commonParameters.PrecursorMassTolerance;
            ProductTolerance = commonParameters.ProductMassTolerance;
            elutionProfileWidthsInScans = new List<double>();
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Averaging MS1 spectra");
            //we are ONLY averaging m/z and NOT intensity. Intensity averaging would convolute downstream quantification
            List<MsDataScan> originalScans = originalFile.GetAllScansList();
            MsDataScan[] ms1scans = originalScans.Where(x => x.MsnOrder == 1).ToArray();
            //we have a set of peaks in the ms1 scan, and we'll cycle through until:
            //-we have two consecutive ms1 scans that do not contain a peak, 
            //-we reach the end of the file
            double[][] ms1mzs = ms1scans.Select(x => x.MassSpectrum.XArray).ToArray();
            List<double>[][] allMs1PeaksFound = new List<double>[ms1mzs.Length][]; //this is a tricky index. Each ms1 scan has an index of List<double>[], where each peak of the ms1 scan has a List<double> that contains all the grouped mzs for that peak.

            for (int seedScanIndex = 0; seedScanIndex < ms1mzs.Length; seedScanIndex++)
            {
                double[] seedMzs = ms1mzs[seedScanIndex]; //grab current ms1scan
                int[] numberOfStrikesForEachPeak = new int[seedMzs.Length]; //Pseudo boolean where 0 is found, 1 is not found, 2 (numberOfStrikesBeforeOut) is out.

                //see what peaks have been found already
                List<double>[] seedPeaksFound = allMs1PeaksFound[seedScanIndex];
                for (int seedPeakIndex = 0; seedPeakIndex < seedMzs.Length; seedPeakIndex++)
                {
                    if (seedPeaksFound[seedPeakIndex] == null) //this is a new peak for us that was not previously seen
                    {
                        seedPeaksFound[seedPeakIndex] = new List<double> { seedMzs[seedPeakIndex] };
                    }
                    else //we saw this already, just use that existing list
                    {
                        numberOfStrikesForEachPeak[seedPeakIndex] = numberOfStrikesBeforeOut;
                    }
                }

                //start cycling through other ms1 scans
                for (int branchScanIndex = seedScanIndex + 1; branchScanIndex < ms1mzs.Length; branchScanIndex++)
                {
                    bool done = true; //this is used to determine if all of the seed peaks have disappeared from the branches
                    double[] branchMzs = ms1mzs[branchScanIndex];

                    int seedPeakIndex = 0;
                    int branchPeakIndex = 0;
                    //foreach peak in the seed spectrum, find matches in the branch spectrum
                    for (; seedPeakIndex < seedMzs.Length; seedPeakIndex++)
                    {
                        if (numberOfStrikesForEachPeak[seedPeakIndex] != numberOfStrikesBeforeOut) //if this seed is not dead to us
                        {
                            double seedMz = seedMzs[seedPeakIndex];
                            //see if it has a buddy! Move the branch until we pass it
                            for (; branchPeakIndex < branchMzs.Length; branchPeakIndex++)
                            {
                                if (branchMzs[branchPeakIndex] > seedMz) //we went too far! or did we?
                                {
                                    break;
                                }
                            }

                            //backup the branchIndex if needed
                            if (branchPeakIndex == branchMzs.Length
                                || (branchPeakIndex != 0 && branchMzs[branchPeakIndex] - seedMz > seedMz - branchMzs[branchPeakIndex - 1]))
                            {
                                branchPeakIndex--;
                            }


                            if (ProductTolerance.Within(seedMz, branchMzs[branchPeakIndex])) //if a match
                            {
                                done = false;
                                numberOfStrikesForEachPeak[seedPeakIndex] = 0; //reset
                                seedPeaksFound[seedPeakIndex].Add(branchMzs[branchPeakIndex]);
                                allMs1PeaksFound[branchScanIndex][branchPeakIndex] = seedPeaksFound[seedPeakIndex]; //update index so future seeds don't have to look
                            }
                            else
                            {
                                if (++numberOfStrikesForEachPeak[seedPeakIndex] != numberOfStrikesBeforeOut)
                                {
                                    done = false;
                                }
                                else //hey, this is dead to us now. Let's do some post processing. (ie average)
                                {
                                    //average the grouped mzs
                                    AverageMzs(seedPeaksFound[seedPeakIndex]);
                                    //get elution profile
                                    int scanIndexesBetweenSeedAndLastBranch = branchScanIndex - numberOfStrikesBeforeOut - seedScanIndex;
                                    if (scanIndexesBetweenSeedAndLastBranch != 0)
                                    {
                                        elutionProfileWidthsInScans.Add(scanIndexesBetweenSeedAndLastBranch);
                                    }
                                }
                            }
                        }
                    }

                    if (done) //if none of the peaks are still viable. 
                    {
                        break; //let's move on to the next seed
                    }
                }

                //Finish residual post-processing for peaks that made it to the end of the file
                for (int peakIndex = 0; peakIndex < numberOfStrikesForEachPeak.Length; peakIndex++)
                {
                    if (numberOfStrikesForEachPeak[peakIndex] != numberOfStrikesBeforeOut) //if we haven't finished this already
                    {
                        AverageMzs(seedPeaksFound[peakIndex]);
                    }
                }

                //All grouping is done for this seed, now update the scan with the aggregated mzs
                List<double> mzsToAdd = new List<double>();
                List<double> intensitiesToAdd = new List<double>();
                double[] unchangedIntensities = ms1scans[seedScanIndex].MassSpectrum.YArray; //get the old intensities
                for (int peakIndex = 0; peakIndex < seedPeaksFound.Length; peakIndex++)
                {
                    List<double> mz = seedPeaksFound[peakIndex];
                    if (mz.Count != 0) // if we didn't remove the peak
                    {
                        mzsToAdd.Add(mz[0]);
                        intensitiesToAdd.Add(unchangedIntensities[peakIndex]);
                    }
                }

                //overwrite the existing spectrum
                ms1scans[seedScanIndex] = CloneDataScanWithNewSpectrum(ms1scans[seedScanIndex], new MzSpectrum(mzsToAdd.ToArray(), intensitiesToAdd.ToArray(), false));
            }

            //update the datafile
            originalFile = new MsDataFile(originalScans.ToArray(), originalFile.SourceFile); //don't need to update ms2 precursor info, because GetMS2Scans does it for us

            //estimate optimal retention time tolerance
            //we previously recorded the difference between scan indexes for elutions, now use those
            double totalRunTime = ms1scans.Last().RetentionTime;
            double estimatedTimeBetweenMs1Scans = totalRunTime / ms1scans.Length;
            double averageElution = elutionProfileWidthsInScans.Average() * estimatedTimeBetweenMs1Scans;
            double innerQuartile = Statistics.InterquartileRange(elutionProfileWidthsInScans) * estimatedTimeBetweenMs1Scans;
            //currently not doing anything with this info

            Status("Getting ms2 scans...");
            Ms2ScanWithSpecificMass[] MS2Scans = GetMs2Scans(originalFile, OriginalFilePath, commonParameters.DoPrecursorDeconvolution, commonParameters.UseProvidedPrecursorInfo, commonParameters.DeconvolutionIntensityRatio, commonParameters.DeconvolutionMaxAssumedChargeState, commonParameters.DeconvolutionMassTolerance).ToArray();

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
                    int maxBinObs = -1; //save outer bounds so we can expand tolerances if needed
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
            List<int> scanNumbersForSyntheticSpectra = new List<int>();
            foreach (List<Ms2ScanWithSpecificMass> group in groups)
            {
                if (group.Count == 1) //if nothing to aggregate, just save it
                {
                    syntheticSpectra.Add(group[0]);
                    scanNumbersForSyntheticSpectra.Add(group[0].OneBasedScanNumber);
                }
                else
                {
                    scanNumbersForSyntheticSpectra.Add(group[group.Count/2].OneBasedScanNumber);

                    List<double> synethicMZs = new List<double>();
                    List<double> syntheticIntensities = new List<double>();

                    //get values
                    double[][] groupedMzs = group.Select(x => x.TheScan.MassSpectrum.XArray).ToArray();
                    double[][] groupedIntensities = group.Select(x => x.TheScan.MassSpectrum.YArray).ToArray();

                    List<double> allMzs = new List<double>();
                    List<double> groupedMzMaxRanges = new List<double>();

                    //Sort array of all mzs sorted low to high
                    allMzs.Sort();

                    //group mzs
                    int seedIndex = 0;
                    for (int i = 0; i < allMzs.Count; i++)
                    {
                        seedIndex = i; //the actual seed
                        double maxValue = ProductTolerance.GetMaximumValue(allMzs[i]);
                        i++;
                        for (; i < allMzs.Count; i++)
                        {
                            double currentMz = allMzs[i];
                            if (currentMz > maxValue) //stop if out of range
                            {
                                break;
                            }
                            else //update range
                            {
                                maxValue = ProductTolerance.GetMaximumValue(currentMz);
                            }
                        }
                        //record
                        groupedMzMaxRanges.Add(maxValue);
                    }

                    //we're done grouping peaks, see if there's enough info (half of all grouped Ms2s must contain this peak group)
                    int numScansNeeded = group.Count / 2;
                    int[] peakPositionArray = new int[group.Count];
                    foreach (double maxValue in groupedMzMaxRanges) //foreach peak group
                    {
                        for(int i=0; i<group.Count; i++) //foreach scan we're looking at here
                        {
                            double[] currentScanMzs = groupedMzs[i]; //get mzs
                            double[] currentScanIntensities = groupedIntensities[i]; //get intensities
                            int j = peakPositionArray[i]; //previous position we left off at
                            for (; j < currentScanMzs.Length; j++)
                            {
                                if(currentScanMzs[j]<maxValue)
                                {

                                }
                            }
                            peakPositionArray[i] = j;
                        }
                    }
                    ///////////
                }
            }

            AggregatedDataFile = new MsDataFile(originalFile.GetAllScansList().ToArray(), originalFile.SourceFile);
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

        public void AverageMzs(List<double> referenceListOfMzs)
        {
            //Currently NOT using intensity for weighting. Reason being that it's more computationally intensive to save those values.
            if (referenceListOfMzs.Count != 1) //if it's worth averaging
            {
                double averageMZ = referenceListOfMzs.Average();
                referenceListOfMzs.Clear(); //need to clear to keep the reference
                referenceListOfMzs.Add(averageMZ); //add the average as the only peak
            }
            else //this peak was only found once... isn't that a little odd? like maybe it's noise?
            {
                referenceListOfMzs.Clear(); //Let's wipe it, but keep the entry to know that it's blank now.
            }
        }

        public MsDataScan CloneDataScanWithNewSpectrum(MsDataScan oldScan, MzSpectrum updatedSpectrum)
        {
            return new MsDataScan(
                updatedSpectrum,
                oldScan.OneBasedScanNumber,
                oldScan.MsnOrder,
                oldScan.IsCentroid,
                oldScan.Polarity,
                oldScan.RetentionTime,
                oldScan.ScanWindowRange,
                oldScan.ScanFilter,
                oldScan.MzAnalyzer,
                oldScan.TotalIonCurrent,
                oldScan.InjectionTime,
                oldScan.NoiseData,
                oldScan.NativeId,
                oldScan.SelectedIonMZ,
                oldScan.SelectedIonChargeStateGuess,
                oldScan.SelectedIonIntensity,
                oldScan.IsolationMz,
                oldScan.IsolationWidth,
                oldScan.DissociationType,
                oldScan.OneBasedPrecursorScanNumber,
                oldScan.SelectedIonMonoisotopicGuessMz
                );
        }
    }
}