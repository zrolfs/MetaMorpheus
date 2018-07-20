using EngineLayer;
using EngineLayer.Aggregation;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using System.Collections.Generic;
using System.IO;

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
                trimMsMsPeaks: false
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

                Status("Writing aggregated spectra...", thisId);
                //write ms file
                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(engine.AggregatedDataFile, aggregatedFilePath, false);
                MyTaskResults.NewSpectra.Add(aggregatedFilePath);

                //write a copy of the old file specific parameters (if applicable)
                if (fileSettingsList[spectraFileIndex] != null)
                {
                    string newTomlFileName = Path.Combine(OutputFolder, origDataFileWithoutExtension + aggregateSuffix + ".toml");
                    Toml.WriteFile(fileSettingsList[spectraFileIndex], newTomlFileName, tomlConfig);
                    MyTaskResults.NewFileSpecificTomls.Add(newTomlFileName);
                }

                //write diagnostics = 
                string diagnosticFileName = Path.Combine(OutputFolder, origDataFileWithoutExtension + "-Diagnostics" + ".txt");
                using (StreamWriter file = new StreamWriter(diagnosticFileName))
                {
                    engine.diagnosticLines.ForEach(x => file.WriteLine(x));
                }
                    Status("Done!", thisId);
            }
            return MyTaskResults;
        }
    }
}