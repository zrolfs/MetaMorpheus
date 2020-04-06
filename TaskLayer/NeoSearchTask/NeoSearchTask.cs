using EngineLayer;
using EngineLayer.Neo;
using MassSpectrometry;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class NeoSearchTask : MetaMorpheusTask
    {
        #region Private Fields

        private List<DbForTask> StoredDatabases = new List<DbForTask>();

        #endregion Private Fields

        #region Public Constructors

        public NeoSearchTask() : base(MyTask.Neo)
        {
            NeoParameters = new NeoParameters();

            IDigestionParams tempDigParams = new DigestionParams
            {
                MinPeptideLength = 8,
                MaxPeptideLength = 16,
                Protease = GlobalVariables.ProteaseDictionary["non-specific"],
                MaxMissedCleavages = 15
            };
            CommonParameters = new CommonParameters
            {
                DigestionParams = tempDigParams,
                DoPrecursorDeconvolution = false,
                PrecursorMassTolerance = null,
                ProductMassTolerance = null
            };
        }

        #endregion Public Constructors

        #region Public Enums

        public enum NeoTaskType { AggregateTargetDecoyFiles, GenerateSplicedPeptides, AggregateNormalSplicedFiles, SearchTransDb, SearchTranslatedDb };

        #endregion Public Enums

        #region Public Properties

        public NeoTaskType NeoType { get; set; }

        public NeoParameters NeoParameters { get; set; }

        #endregion Public Properties

        #region Public Methods

        public NeoSearchTask Clone()
        {
            return (NeoSearchTask)this.MemberwiseClone();
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificSettings[] fileSettingsList)
        {
            myTaskResults = new MyTaskResults(this);
            NeoFindAmbiguity.numInterveningResidues = NeoParameters.MaxDistanceAllowed;
            NeoFindAmbiguity.normalCis = NeoParameters.NormalCis;
            NeoFindAmbiguity.reverseCis = NeoParameters.ReverseCis;
            NeoFindAmbiguity.maxMissingConsecutivePeaks = NeoParameters.MaxMissedConsecutiveFragments;
            NeoFindAmbiguity.maxNumPossibleSequences = NeoParameters.MaxCandidatesPerSpectrum;

            //order list of files so we don't look for the products of one raw file in a different raw file (will lead to inconsistent and frustrating crashes)
            currentRawFileList.Sort();
            NeoParameters.TargetFilePath.Sort();
            NeoParameters.DecoyFilePath.Sort();
            NeoParameters.NFilePath.Sort();
            NeoParameters.CFilePath.Sort();

            List<ProductType> ionTypes = new List<ProductType>();
            if (CommonParameters.BIons)
                ionTypes.Add(ProductType.B);
            if (CommonParameters.YIons)
                ionTypes.Add(ProductType.Y);
            if (CommonParameters.ZdotIons)
                ionTypes.Add(ProductType.Zdot);
            if (CommonParameters.CIons)
                ionTypes.Add(ProductType.C);
            NeoFindAmbiguity.ionsUsed = ionTypes;

            if (NeoType.Equals(NeoTaskType.AggregateTargetDecoyFiles))
            {
                //getfolders
                for (int i = 0; i < currentRawFileList.Count; i++)
                {
                    if (NeoParameters.DecoySearch)
                    {
                        string decoyPath = new DirectoryInfo(OutputFolder).Name;
                        string taskString = decoyPath.Split('-')[0];
                        int taskNum = Convert.ToInt32(taskString.Substring(4, taskString.Length - 4));
                        taskNum--;
                        NeoParameters.DecoyFilePath.Add(OutputFolder.Substring(0, OutputFolder.Length - decoyPath.Length) + "Task" + taskNum + "-SearchTask\\" + Path.GetFileNameWithoutExtension(currentRawFileList[i]) + "_PSMs.psmtsv");
                        if (NeoParameters.TargetSearch)
                        {
                            string targetPath = new DirectoryInfo(OutputFolder).Name;
                            taskNum--;
                            NeoParameters.TargetFilePath.Add(OutputFolder.Substring(0, OutputFolder.Length - targetPath.Length) + "Task" + taskNum + "-SearchTask\\" + Path.GetFileNameWithoutExtension(currentRawFileList[i]) + "_PSMs.psmtsv");
                        }
                    }
                    else if (NeoParameters.TargetSearch)
                    {
                        string targetPath = new DirectoryInfo(OutputFolder).Name;
                        string taskString = targetPath.Split('-')[0];
                        int taskNum = Convert.ToInt32(taskString.Substring(4, taskString.Length - 4));
                        taskNum--;
                        NeoParameters.TargetFilePath.Add(OutputFolder.Substring(0, OutputFolder.Length - targetPath.Length) + "Task" + taskNum + "-SearchTask\\" + Path.GetFileNameWithoutExtension(currentRawFileList[i]) + "_PSMs.psmtsv");
                    }
                    else
                    {
                        //do nothing, both exist
                    }
                    AggregateSearchFiles.Combine(NeoParameters.TargetFilePath[i], NeoParameters.DecoyFilePath[i], OutputFolder + "\\" + Path.GetFileNameWithoutExtension(currentRawFileList[i]));
                }
            }
            else if (NeoType.Equals(NeoTaskType.AggregateNormalSplicedFiles))
            {
                //get folder info
                string currentPath = new DirectoryInfo(OutputFolder).Name;
                string taskString = currentPath.Split('-')[0];
                string outputFolderSubstring = OutputFolder.Substring(0, OutputFolder.Length - currentPath.Length);
                for (int i = 0; i < currentRawFileList.Count; i++)
                {
                    int taskNum = Convert.ToInt32(taskString.Substring(4, taskString.Length - 4)) - 1;
                    string translatedPath = outputFolderSubstring + "Task" + taskNum + "-SearchTask\\" + Path.GetFileNameWithoutExtension(currentRawFileList[i]) + "_PSMs.psmtsv";

                    //combine TL and traditional search
                    //getfolders
                    AggregateSearchFiles.Combine(NeoParameters.TargetFilePath[i], NeoParameters.DecoyFilePath[i], OutputFolder + "\\" + Path.GetFileNameWithoutExtension(currentRawFileList[i]) + "_NonsplicedAggregate", translatedPath);

                    //reset database
                    myTaskResults.newDatabases = StoredDatabases;

                    //assign all file paths using folder info
                    taskNum -= 2;
                    string transPath = outputFolderSubstring + "Task" + taskNum + "-SearchTask\\" + Path.GetFileNameWithoutExtension(currentRawFileList[i]) + "_PSMs.psmtsv";
                    taskNum -= 2;
                    string cisPath = outputFolderSubstring + "Task" + taskNum + "-SearchTask\\" + Path.GetFileNameWithoutExtension(currentRawFileList[i]) + "_PSMs.psmtsv";
                    string normalPath = OutputFolder + "\\" + Path.GetFileNameWithoutExtension(currentRawFileList[i]) + "_NonsplicedAggregate_Targets.psmtsv";
                    string transAggPath = OutputFolder + "\\" + Path.GetFileNameWithoutExtension(currentRawFileList[i]) + "_TransAggregate_Targets.psmtsv";

                    //combine cis and trans for trans analysis
                    AggregateSearchFiles.Combine(cisPath, NeoParameters.DecoyFilePath[i], OutputFolder + "\\" + Path.GetFileNameWithoutExtension(currentRawFileList[i]) + "_TransAggregate", transPath); //this is trans and cis thrown together for trans analysis

                    //do cis and trans (separate) analysis
                    AggregateSearchFiles.RecursiveNeoAggregation(normalPath, cisPath, OutputFolder, Path.GetFileNameWithoutExtension(currentRawFileList[i]) + "_CisResults.psmtsv");
                    AggregateSearchFiles.RecursiveNeoAggregation(normalPath, transAggPath, OutputFolder, Path.GetFileNameWithoutExtension(currentRawFileList[i]) + "_TransResults.psmtsv");
                }
            }
            else if (NeoType.Equals(NeoTaskType.GenerateSplicedPeptides))
            {
                NeoMassCalculator.ImportMasses();

                ParallelOptions parallelOptions = new ParallelOptions();
                if (CommonParameters.MaxParallelFilesToAnalyze.HasValue)
                    parallelOptions.MaxDegreeOfParallelism = CommonParameters.MaxParallelFilesToAnalyze.Value;
                MyFileManager myFileManager = new MyFileManager(true);
                myTaskResults.newDatabases = new List<DbForTask>();

                //Import Database
                Status("Loading modifications...", taskId);

                #region Load modifications

                List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
                List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();
                List<string> localizeableModificationTypes = CommonParameters.ListOfModTypesLocalize == null ? new List<string>() : CommonParameters.ListOfModTypesLocalize.ToList();
                if (CommonParameters.LocalizeAll)
                    localizeableModificationTypes = GlobalVariables.AllModTypesKnown.ToList();
                else
                    localizeableModificationTypes = GlobalVariables.AllModTypesKnown.Where(b => localizeableModificationTypes.Contains(b)).ToList();

                #endregion Load modifications

                NeoFindAmbiguity.theoreticalProteins = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, true, DecoyType.None, localizeableModificationTypes, b.IsContaminant, out Dictionary<string, Modification> unknownModifications)).ToList();

                List<string> dbPathList = dbFilenameList.Select(x => x.FilePath).OrderBy(x => x).ToList();
                string concatonatedDatabaseString = dbPathList[0];
                dbPathList.RemoveAt(0);
                foreach(string path in dbPathList)
                {
                    concatonatedDatabaseString += "_" + Path.GetFileNameWithoutExtension(path);
                }
                    NeoFindAmbiguity.PopulateSequenceLookUpDictionaries(concatonatedDatabaseString);

                //Import Spectra
                for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
                {
                    var origDataFile = currentRawFileList[spectraFileIndex];
                    ICommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);
                    //get parameters updated
                    NeoFindAmbiguity.precursorMassTolerancePpm = combinedParams.PrecursorMassTolerance.Value;
                    NeoFindAmbiguity.productMassTolerancePpm = combinedParams.ProductMassTolerance.Value;

                    var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                    NewCollection(Path.GetFileName(origDataFile), thisId);
                    Status("Loading spectra file...", thisId);
                    IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks);
                    Status("Getting ms2 scans...", thisId);
                    Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams.DoPrecursorDeconvolution, combinedParams.UseProvidedPrecursorInfo, combinedParams.DeconvolutionIntensityRatio, combinedParams.DeconvolutionMaxAssumedChargeState, combinedParams.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
                    Ms2ScanWithSpecificMass[] indexedSpectra = new Ms2ScanWithSpecificMass[arrayOfMs2ScansSortedByMass.Max(x => x.OneBasedScanNumber) + 1];
                    foreach (Ms2ScanWithSpecificMass scan in arrayOfMs2ScansSortedByMass)
                        indexedSpectra[scan.OneBasedScanNumber] = scan;

                    //if termini input

                    if (NeoParameters.SearchNTerminus || NeoParameters.SearchCTerminus)
                    {
                        string nPath = "";
                        string cPath = "";
                        //if no termini input
                        string taskHeader = "Task";
                        string[] pathArray = OutputFolder.Split('\\');
                        string basePath = "";
                        for (int i = 0; i < pathArray.Length - 1; i++)
                            basePath += pathArray[i] + '\\';
                        string currentTaskNumber = pathArray[pathArray.Length - 1].Split('-')[0];
                        currentTaskNumber = currentTaskNumber.Substring(taskHeader.Length, currentTaskNumber.Length - taskHeader.Length);
                        string NHeader = "";
                        string CHeader = "";
                        if (NeoParameters.SearchCTerminus)
                        {
                            CHeader = taskHeader + (Convert.ToInt16(currentTaskNumber) - 1);
                            if (NeoParameters.SearchNTerminus)
                                NHeader = taskHeader + (Convert.ToInt16(currentTaskNumber) - 2);
                        }
                        else
                            NHeader = taskHeader + (Convert.ToInt16(currentTaskNumber) - 1);
                        foreach (string s in Directory.GetDirectories(basePath))
                        {
                            if (s.Contains(NHeader))
                                nPath = s;
                            else if (s.Contains(CHeader))
                                cPath = s;
                        }
                        string fileName = Path.GetFileNameWithoutExtension(currentRawFileList[spectraFileIndex]) + "_PSMs.psmtsv";
                        nPath += "\\" + fileName;
                        cPath += "\\" + fileName;

                        //Read N and C files
                        lock (NeoParameters.NFilePath)
                        {
                            while (NeoParameters.NFilePath.Count <= spectraFileIndex)
                            {
                                NeoParameters.NFilePath.Add("");
                                NeoParameters.CFilePath.Add("");
                            }
                            NeoParameters.NFilePath[spectraFileIndex] = nPath;
                            NeoParameters.CFilePath[spectraFileIndex] = cPath;
                        }
                    }

                    Status("Importing Search Results...", taskId);
                    List<NeoPsm> psms = ImportPsmtsv.ImportNeoPsms(NeoParameters.NFilePath[spectraFileIndex], NeoParameters.CFilePath[spectraFileIndex]);

                    //Splice
                    Status("Splicing Fragments...", taskId);
                    List<NeoPsm> candidates = NeoSplicePeptides.SplicePeptides(psms);

                    //Find Ambiguity
                    Status("Identifying Ambiguity...", taskId);
                    NeoFindAmbiguity.FindAmbiguity(candidates, indexedSpectra);

                    //Export Results
                    Status("Exporting Results...", taskId);
                    NeoExport.ExportAll(candidates, indexedSpectra, OutputFolder, Path.GetFileNameWithoutExtension(currentRawFileList[spectraFileIndex]));

                    //Switch databases

                    string outputFolder = NeoExport.path + NeoExport.folder + @"\" + NeoExport.folder + "_" + Path.GetFileNameWithoutExtension(currentRawFileList[spectraFileIndex]) + "_FusionDatabaseAppendixNC.fasta";
                    myTaskResults.newDatabases.Add(new DbForTask(outputFolder, false));
                }
                StoredDatabases = dbFilenameList;
            }
            else if (NeoType.Equals(NeoTaskType.SearchTransDb))
            {
                myTaskResults.newDatabases = new List<DbForTask>();
                for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
                {
                    string outputFolder = NeoExport.path + NeoExport.folder + @"\" + NeoExport.folder + "_" + Path.GetFileNameWithoutExtension(currentRawFileList[spectraFileIndex]) + "_FusionDatabaseAppendixTS.fasta";
                    myTaskResults.newDatabases.Add(new DbForTask(outputFolder, false));
                }
            }
            else //if SearchTranslatedDb
            {
                myTaskResults.newDatabases = new List<DbForTask>();
                for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
                {
                    string outputFolder = NeoExport.path + NeoExport.folder + @"\" + NeoExport.folder + "_" + Path.GetFileNameWithoutExtension(currentRawFileList[spectraFileIndex]) + "_FusionDatabaseAppendixTL.fasta";
                    myTaskResults.newDatabases.Add(new DbForTask(outputFolder, false));
                }
            }

            return myTaskResults;
        }

        #endregion Protected Methods
    }
}