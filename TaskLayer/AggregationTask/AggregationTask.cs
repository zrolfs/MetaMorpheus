using EngineLayer;
using EngineLayer.Aggregation;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class AggregationTask : MetaMorpheusTask
    {
        private const string aggregateSuffix = "-agg";

        public AggregationTask() : base(MyTask.Aggregate)
        {
            CommonParameters = new CommonParameters(
                productMassTolerance: new PpmTolerance(25),
                precursorMassTolerance: new PpmTolerance(15),
                trimMsMsPeaks: false,
                doPrecursorDeconvolution: false,
                scoreCutoff: 10);

            AggregationParameters = new AggregationParameters();
        }

        public AggregationParameters AggregationParameters { get; set; }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            // write prose settings
            ProseCreatedWhileRunning.Append("The following Aggregation settings were used: ");
            ProseCreatedWhileRunning.Append("Precursor mass tolerance = " + CommonParameters.PrecursorMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("Product mass tolerance = " + CommonParameters.ProductMassTolerance + ". ");
            ProseCreatedWhileRunning.Append("Max retention time difference allowed (min) = " + CommonParameters.ProductMassTolerance + ". ");
            ProseCreatedWhileRunning.Append("Min cosine score allowed = " + CommonParameters.ProductMassTolerance + ". ");
            ProseCreatedWhileRunning.Append("Number of MS1 spectra to average = " + CommonParameters.ProductMassTolerance + ". ");

            // start the Aggregation task
            Status("Aggregating...", new List<string> { taskId });
            MyTaskResults = new MyTaskResults(this)
            {
                NewSpectra = new List<string>(),
                NewFileSpecificTomls = new List<string>()
            };

            object lock1 = new object();

            var myFileManager = new MyFileManager(true);

            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                if (GlobalVariables.StopLoops) { break; }

                // get filename stuff
                var originalUnAggregatedFilePath = currentRawFileList[spectraFileIndex];
                var originalUnAggregatedFilenameWithoutExtension = Path.GetFileNameWithoutExtension(originalUnAggregatedFilePath);
                string aggregatedFilePath = Path.Combine(OutputFolder, originalUnAggregatedFilenameWithoutExtension + aggregateSuffix + ".mzML");

                // mark the file as in-progress
                StartingDataFile(originalUnAggregatedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUnAggregatedFilePath });

                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                MsDataFile myMsDataFile;

                // load the file
                Status("Loading spectra file...", new List<string> { taskId, "Individual Spectra Files" });
                lock (lock1)
                {
                    myMsDataFile = myFileManager.LoadFile(originalUnAggregatedFilePath, CommonParameters.TopNpeaks, CommonParameters.MinRatio, CommonParameters.TrimMs1Peaks, CommonParameters.TrimMsMsPeaks);
                }

                //normalize the ms2 spectra to facilitate 
                //Status("Normalizing MS2 scans...", new List<string> { taskId, "Individual Spectra Files" });
                //const int numHighestIntensityPeaksToSumAndNormalizeTo = 10;
                //List<MsDataScan> scans = new List<MsDataScan>();
                AggregationEngine engine = new AggregationEngine(myMsDataFile, CommonParameters, new List<string> { taskId, "Individual Spectra Files", originalUnAggregatedFilenameWithoutExtension }, AggregationParameters.MaxRetentionTimeDifferenceAllowedInMinutes, AggregationParameters.MinCosineScoreAllowed, AggregationParameters.NumberOfMS1SpectraToAverage);

                // get datapoints to fit Aggregation function to
                Status("Aggregating data points...", new List<string> { taskId, "Individual Spectra Files" });

                //Toml.WriteFile(fileSpecificParams, newTomlFileName, tomlConfig);

               // SucessfullyFinishedWritingFile(newTomlFileName, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, aggregatedFilePath, false);
            }
            return MyTaskResults;
        }
    }
}