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
                productMassTolerance: new PpmTolerance(15),
                precursorMassTolerance: new PpmTolerance(4),
                trimMsMsPeaks: false,
                doPrecursorDeconvolution: false
                );

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
                var origDataFile = currentRawFileList[spectraFileIndex];
                var origDataFileWithoutExtension = Path.GetFileNameWithoutExtension(origDataFile);
                string aggregatedFilePath = Path.Combine(OutputFolder, origDataFileWithoutExtension + aggregateSuffix + ".mzML");

                // mark the file as in-progress
                StartingDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });

                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                Status("Loading spectra file...", thisId);
                MsDataFile myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks);

                // get datapoints to fit Aggregation function to
                Status("Aggregating data points...", thisId);
                AggregationEngine engine = new AggregationEngine(myMsDataFile, origDataFile, combinedParams, new List<string> { taskId, "Individual Spectra Files", origDataFileWithoutExtension }, AggregationParameters.MaxRetentionTimeDifferenceAllowedInMinutes, AggregationParameters.MinCosineScoreAllowed);
                engine.Run();



                //Toml.WriteFile(fileSpecificParams, newTomlFileName, tomlConfig);

                // SucessfullyFinishedWritingFile(newTomlFileName, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, aggregatedFilePath, false);
            }
            return MyTaskResults;
        }
    }
}