﻿using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.HistogramAnalysis;
using EngineLayer.Localization;
using EngineLayer.ModificationAnalysis;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class PostSearchAnalysisTask : MetaMorpheusTask
    {
        public PostSearchAnalysisParameters Parameters { get; set; }
        private List<EngineLayer.ProteinGroup> ProteinGroups { get; set; }
        private IEnumerable<IGrouping<string, PeptideSpectralMatch>> PsmsGroupedByFile { get; set; }

        public PostSearchAnalysisTask()
            : base(MyTask.Search)
        {
        }

        public MyTaskResults Run()
        {
            // Stop loop if canceled
            if (GlobalVariables.StopLoops) { return Parameters.SearchTaskResults; }

            if (Parameters.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.ModOpen
                || Parameters.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.Open
                || Parameters.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.Custom
                )
            {
                // This only makes sense if there is a mass difference that you want to localize. No use for exact and missed monoisotopic mass searches.
                Parameters.SearchParameters.DoLocalizationAnalysis = true;
            }
            else
            {
                Parameters.SearchParameters.DoLocalizationAnalysis = false;
            }

            //update all psms with peptide info
            if (Parameters.SearchParameters.SearchType != SearchType.NonSpecific) //if it hasn't been done already
            {
                Parameters.AllPsms = Parameters.AllPsms.Where(psm => psm != null).ToList();
                Parameters.AllPsms.ForEach(psm => psm.ResolveAllAmbiguities());

                Parameters.AllPsms = Parameters.AllPsms.OrderByDescending(b => b.Score)
                   .ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue)
                   .GroupBy(b => (b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

                CalculatePsmFdr();
            }

            DoMassDifferenceLocalizationAnalysis();
            ProteinAnalysis();
            QuantificationAnalysis();

            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files" }));

            HistogramAnalysis();
            WritePsmResults();
            WriteProteinResults();
            WriteQuantificationResults();
            WritePrunedDatabase();
            WritePeptideResults(); // modifies the FDR results for PSMs, so do this last

            return Parameters.SearchTaskResults;
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            return null;
        }

        /// <summary>
        /// Calculate estimated false-discovery rate (FDR) for peptide spectral matches (PSMs)
        /// </summary>
        private void CalculatePsmFdr()
        {
            // TODO: because FDR is done before parsimony, if a PSM matches to a target and a decoy protein, there may be conflicts between how it's handled in parsimony and the FDR engine here
            // for example, here it may be treated as a decoy PSM, where as in parsimony it will be determined by the parsimony algorithm which is agnostic of target/decoy assignments
            // this could cause weird PSM FDR issues

            Status("Estimating PSM FDR...", Parameters.SearchTaskId);
            new FdrAnalysisEngine(Parameters.AllPsms, Parameters.NumNotches, CommonParameters, new List<string> { Parameters.SearchTaskId }).Run();

            // sort by q-value because of group FDR stuff
            // e.g. multiprotease FDR, non/semi-specific protease, etc
            Parameters.AllPsms = Parameters.AllPsms
                .OrderBy(p => p.FdrInfo.QValue)
                .ThenByDescending(p => p.Score)
                .ThenBy(p => p.FdrInfo.CumulativeTarget)
                .ToList();

            Status("Done estimating PSM FDR!", Parameters.SearchTaskId);
        }

        private void ProteinAnalysis()
        {
            if (!Parameters.SearchParameters.DoParsimony)
            {
                return;
            }

            Status("Constructing protein groups...", Parameters.SearchTaskId);

            //if SILAC
            List<PeptideSpectralMatch> psmsForProteinParsimony = Parameters.AllPsms;
            if (Parameters.SilacLabels != null)
            {
                psmsForProteinParsimony = SilacConversions.UpdatePsmsForParsimony(Parameters.SilacLabels, psmsForProteinParsimony);
            }

            // run parsimony
            ProteinParsimonyResults proteinAnalysisResults = (ProteinParsimonyResults)(new ProteinParsimonyEngine(psmsForProteinParsimony, Parameters.SearchParameters.ModPeptidesAreDifferent, CommonParameters, new List<string> { Parameters.SearchTaskId }).Run());

            // score protein groups and calculate FDR
            ProteinScoringAndFdrResults proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psmsForProteinParsimony,
                Parameters.SearchParameters.NoOneHitWonders, Parameters.SearchParameters.ModPeptidesAreDifferent, true, CommonParameters, new List<string> { Parameters.SearchTaskId }).Run();

            ProteinGroups = proteinScoringAndFdrResults.SortedAndScoredProteinGroups;

            foreach (PeptideSpectralMatch psm in Parameters.AllPsms)
            {
                psm.ResolveAllAmbiguities();
            }

            Status("Done constructing protein groups!", Parameters.SearchTaskId);
        }

        private void DoMassDifferenceLocalizationAnalysis()
        {
            if (Parameters.SearchParameters.DoLocalizationAnalysis)
            {
                Status("Running mass-difference localization analysis...", Parameters.SearchTaskId);
                for (int spectraFileIndex = 0; spectraFileIndex < Parameters.CurrentRawFileList.Count; spectraFileIndex++)
                {
                    CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, Parameters.FileSettingsList[spectraFileIndex]);

                    var origDataFile = Parameters.CurrentRawFileList[spectraFileIndex];
                    Status("Running mass-difference localization analysis...", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", origDataFile });
                    MsDataFile myMsDataFile = Parameters.MyFileManager.LoadFile(origDataFile, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks, combinedParams);
                    new LocalizationEngine(Parameters.AllPsms.Where(b => b.FullFilePath.Equals(origDataFile)).ToList(),
                        myMsDataFile, combinedParams, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", origDataFile }).Run();
                    Parameters.MyFileManager.DoneWithFile(origDataFile);
                    ReportProgress(new ProgressEventArgs(100, "Done with localization analysis!", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", origDataFile }));
                }
            }

            // count different modifications observed
            new ModificationAnalysisEngine(Parameters.AllPsms, CommonParameters, new List<string> { Parameters.SearchTaskId }).Run();
        }

        private void QuantificationAnalysis()
        {
            if (!Parameters.SearchParameters.DoQuantification)
            {
                return;
            }

            // pass quantification parameters to FlashLFQ
            Status("Quantifying...", Parameters.SearchTaskId);

            // construct file info for FlashLFQ
            var spectraFileInfo = new List<SpectraFileInfo>();

            // get experimental design info for normalization
            if (Parameters.SearchParameters.Normalize)
            {
                string assumedExperimentalDesignPath = Directory.GetParent(Parameters.CurrentRawFileList.First()).FullName;
                assumedExperimentalDesignPath = Path.Combine(assumedExperimentalDesignPath, GlobalVariables.ExperimentalDesignFileName);

                if (File.Exists(assumedExperimentalDesignPath))
                {
                    var experimentalDesign = File.ReadAllLines(assumedExperimentalDesignPath)
                        .ToDictionary(p => p.Split('\t')[0], p => p);

                    foreach (var file in Parameters.CurrentRawFileList)
                    {
                        string filename = Path.GetFileNameWithoutExtension(file);

                        var expDesignForThisFile = experimentalDesign[filename];
                        var split = expDesignForThisFile.Split('\t');

                        string condition = split[1];
                        int biorep = int.Parse(split[2]);
                        int fraction = int.Parse(split[3]);
                        int techrep = int.Parse(split[4]);

                        // experimental design info passed in here for each spectra file
                        spectraFileInfo.Add(new SpectraFileInfo(fullFilePathWithExtension: file,
                            condition: condition,
                            biorep: biorep - 1,
                            fraction: fraction - 1,
                            techrep: techrep - 1));

                        Parameters.MyFileManager.DoneWithFile(file);
                    }
                }
                else
                {
                    throw new MetaMorpheusException("Could not find experimental design file at location:\n" + assumedExperimentalDesignPath);
                }
            }
            else
            {
                foreach (var file in Parameters.CurrentRawFileList)
                {
                    // experimental design info passed in here for each spectra file
                    spectraFileInfo.Add(new SpectraFileInfo(fullFilePathWithExtension: file, condition: "", biorep: 0, fraction: 0, techrep: 0));
                    Parameters.MyFileManager.DoneWithFile(file);
                }
            }

            // get PSMs to pass to FlashLFQ
            var unambiguousPsmsBelowOnePercentFdr = Parameters.AllPsms.Where(p =>
                p.FdrInfo.QValue <= 0.01
                && p.FdrInfo.QValueNotch <= 0.01
                && !p.IsDecoy
                && p.FullSequence != null).ToList(); //if ambiguous, there's no full sequence

            //If SILAC (Pre-Quantification) preparation
            Dictionary<PeptideSpectralMatch, List<PeptideSpectralMatch>> silacPsmMatcher = new Dictionary<PeptideSpectralMatch, List<PeptideSpectralMatch>>(); //use for getting back info once we've separated light/heavy. It is the light psm for the key, and the list of heavy psm replicas for the value.
            Dictionary<string, List<string>> silacProteinGroupMatcher = new Dictionary<string, List<string>>(); //use for getting back info once we've separated light/heavy. It is the light protein group for the key, and a list of heavy protein group replicas for the value.
            if (Parameters.SilacLabels != null)
            {
                //PEPTIDE LEVEL CHANGES
                //go through all the psms and duplicate them until a psm copy exists for the unlabeled and labeled proteins
                //The number of psms should roughly increase by a factor of N, where N is the number of labels.
                //It may not increase exactly by a factor of N if the amino acid that gets labeled doesn't exist in the peptide
                List<SilacLabel> allSilacLabels = Parameters.SilacLabels;
                List<PeptideSpectralMatch> silacPsms = new List<PeptideSpectralMatch>();
                foreach (PeptideSpectralMatch psm in unambiguousPsmsBelowOnePercentFdr)
                {
                    //see which label, if any, this peptide has
                    string peptideBaseSequence = psm.BaseSequence;
                    SilacLabel observedLabel = SilacConversions.GetRelevantLabelFromBaseSequence(peptideBaseSequence, allSilacLabels); //returns null if no label

                    //if it's not the light form, make a light form and add it
                    PeptideSpectralMatch lightPsm = observedLabel == null ? psm : SilacConversions.GetSilacPsm(psm, observedLabel, true);
                    silacPsms.Add(lightPsm);

                    //see which labels, if any, this peptide can have
                    peptideBaseSequence = SilacConversions.GetSilacLightBaseSequence(peptideBaseSequence, observedLabel); //update to the light sequence, if applicable
                    List<SilacLabel> labelsToModify = allSilacLabels.Where(x => peptideBaseSequence.Contains(x.OriginalAminoAcid)).ToList(); //see which residues can be labeled
                    //foreach of the possible labels this peptide can have, add a new psm for it
                    List<PeptideSpectralMatch> addedSilacPsms = labelsToModify.Select(x => SilacConversions.GetSilacPsm(lightPsm, x, false)).ToList(); //if we had a heavy psm to start with, this will create an identical heavy
                    addedSilacPsms.ForEach(x => silacPsms.Add(x)); //add to psm list
                    silacPsmMatcher[lightPsm] = addedSilacPsms; //add to match list
                }

                //update the list for FlashLFQ
                silacPsms.ForEach(x => x.ResolveAllAmbiguities()); //update the base/full sequences
                unambiguousPsmsBelowOnePercentFdr = silacPsms;


                //PROTEIN LEVEL CHANGES
                //Need to create the heavy proteins in our parsimony list and assign the heavy peptides to these heavy proteins
                //During the search and protein parsimony, heavy and light peptides were considered to come from the same protein (light)
                //Now, we want to quantify the light and heavy proteins separately and we need to differentiate them.
                if (ProteinGroups != null) //if we did parsimony
                {
                    List<EngineLayer.ProteinGroup> silacProteinGroups = new List<EngineLayer.ProteinGroup>();
                    foreach (EngineLayer.ProteinGroup proteinGroup in ProteinGroups)
                    {
                        //add the unlabeled protein group
                        //this method removes the heavy psms from the protein groups and adds their light replica
                        EngineLayer.ProteinGroup unlabeledProteinGroup = SilacConversions.GetSilacProteinGroups(unambiguousPsmsBelowOnePercentFdr, proteinGroup);
                        silacProteinGroups.Add(unlabeledProteinGroup);
                        //add the labeled protein group(s)
                        List<EngineLayer.ProteinGroup> addedProteinGroups = allSilacLabels.Select(x => SilacConversions.GetSilacProteinGroups(unambiguousPsmsBelowOnePercentFdr, proteinGroup, x)).ToList(); //foreach label, create a new heavy protein group
                        addedProteinGroups.ForEach(x => silacProteinGroups.Add(x)); //add to psm list
                        silacProteinGroupMatcher[proteinGroup.ProteinGroupName] = addedProteinGroups.Select(x => x.ProteinGroupName).ToList(); //add to match list
                    }
                    ProteinGroups = silacProteinGroups;
                }
            }


            // pass protein group info for each PSM
            var psmToProteinGroups = new Dictionary<PeptideSpectralMatch, List<FlashLFQ.ProteinGroup>>();
            if (ProteinGroups != null)
            {
                foreach (var proteinGroup in ProteinGroups)
                {
                    var proteinsOrderedByAccession = proteinGroup.Proteins.OrderBy(p => p.Accession);

                    var flashLfqProteinGroup = new FlashLFQ.ProteinGroup(proteinGroup.ProteinGroupName,
                        string.Join("|", proteinsOrderedByAccession.Select(p => p.GeneNames.Select(x => x.Item2).FirstOrDefault())),
                        string.Join("|", proteinsOrderedByAccession.Select(p => p.Organism).Distinct()));

                    foreach (var psm in proteinGroup.AllPsmsBelowOnePercentFDR.Where(v => v.FullSequence != null))
                    {
                        if (psmToProteinGroups.TryGetValue(psm, out var flashLfqProteinGroups))
                        {
                            flashLfqProteinGroups.Add(flashLfqProteinGroup);
                        }
                        else
                        {
                            psmToProteinGroups.Add(psm, new List<FlashLFQ.ProteinGroup> { flashLfqProteinGroup });
                        }
                    }
                }
            }
            else
            {
                // if protein groups were not constructed, just use accession numbers
                var accessionToPg = new Dictionary<string, FlashLFQ.ProteinGroup>();
                foreach (var psm in unambiguousPsmsBelowOnePercentFdr)
                {
                    var proteins = psm.BestMatchingPeptides.Select(b => b.Peptide.Protein).Distinct();

                    foreach (var protein in proteins)
                    {
                        if (!accessionToPg.ContainsKey(protein.Accession))
                        {
                            accessionToPg.Add(protein.Accession, new FlashLFQ.ProteinGroup(protein.Accession, string.Join("|", protein.GeneNames.Select(p => p.Item2).Distinct()), protein.Organism));
                        }

                        if (psmToProteinGroups.TryGetValue(psm, out var proteinGroups))
                        {
                            proteinGroups.Add(accessionToPg[protein.Accession]);
                        }
                        else
                        {
                            psmToProteinGroups.Add(psm, new List<FlashLFQ.ProteinGroup> { accessionToPg[protein.Accession] });
                        }
                    }
                }
            }

            //group psms by file
            var psmsGroupedByFile = unambiguousPsmsBelowOnePercentFdr.GroupBy(p => p.FullFilePath);

            // some PSMs may not have protein groups (if 2 peptides are required to construct a protein group, some PSMs will be left over)
            // the peptides should still be quantified but not considered for protein quantification
            var undefinedPg = new FlashLFQ.ProteinGroup("UNDEFINED", "", "");
            //sort the unambiguous psms by protease to make MBR compatible with multiple proteases
            Dictionary<Protease, List<PeptideSpectralMatch>> proteaseSortedPsms = new Dictionary<Protease, List<PeptideSpectralMatch>>();
            Dictionary<Protease, FlashLfqResults> proteaseSortedFlashLFQResults = new Dictionary<Protease, FlashLfqResults>();

            foreach (DigestionParams dp in Parameters.ListOfDigestionParams)
            {
                if (!proteaseSortedPsms.ContainsKey(dp.Protease))
                {
                    proteaseSortedPsms.Add(dp.Protease, new List<PeptideSpectralMatch>());
                }
            }
            foreach (var psm in unambiguousPsmsBelowOnePercentFdr)
            {
                if (!psmToProteinGroups.ContainsKey(psm))
                {
                    psmToProteinGroups.Add(psm, new List<FlashLFQ.ProteinGroup> { undefinedPg });
                }

                proteaseSortedPsms[psm.DigestionParams.Protease].Add(psm);
            }

            // pass PSM info to FlashLFQ
            var flashLFQIdentifications = new List<Identification>();
            foreach (var spectraFile in psmsGroupedByFile)
            {
                var rawfileinfo = spectraFileInfo.Where(p => p.FullFilePathWithExtension.Equals(spectraFile.Key)).First();

                foreach (var psm in spectraFile)
                {
                    flashLFQIdentifications.Add(new Identification(rawfileinfo, psm.BaseSequence, psm.FullSequence,
                        psm.PeptideMonisotopicMass.Value, psm.ScanRetentionTime, psm.ScanPrecursorCharge, psmToProteinGroups[psm]));
                }
            }

            // run FlashLFQ
            var FlashLfqEngine = new FlashLfqEngine(
                allIdentifications: flashLFQIdentifications,
                normalize: Parameters.SearchParameters.Normalize,
                ppmTolerance: Parameters.SearchParameters.QuantifyPpmTol,
                matchBetweenRuns: Parameters.SearchParameters.MatchBetweenRuns,
                silent: true,
                optionalPeriodicTablePath: GlobalVariables.ElementsLocation,
                maxThreads: CommonParameters.MaxThreadsToUsePerFile);

            if (flashLFQIdentifications.Any())
            {
                Parameters.FlashLfqResults = FlashLfqEngine.Run();
            }

            //MultiProtease MBR capability code
            //Parameters.FlashLfqResults = null;

            //foreach (var proteasePsms in proteaseSortedPsms)
            //{
            //    var flashLFQIdentifications = new List<Identification>();
            //    var proteasePsmsGroupedByFile = proteasePsms.Value.GroupBy(p => p.FullFilePath);
            //    foreach (var spectraFile in proteasePsmsGroupedByFile)
            //    {
            //        var rawfileinfo = spectraFileInfo.Where(p => p.FullFilePathWithExtension.Equals(spectraFile.Key)).First();

            //        foreach (var psm in spectraFile)
            //        {
            //            flashLFQIdentifications.Add(new Identification(rawfileinfo, psm.BaseSequence, psm.FullSequence,
            //                psm.PeptideMonisotopicMass.Value, psm.ScanRetentionTime, psm.ScanPrecursorCharge, psmToProteinGroups[psm]));
            //        }
            //    }

            //    // run FlashLFQ
            //    var FlashLfqEngine = new FlashLFQEngine(
            //        allIdentifications: flashLFQIdentifications,
            //        normalize: Parameters.SearchParameters.Normalize,
            //        ppmTolerance: Parameters.SearchParameters.QuantifyPpmTol,
            //        matchBetweenRuns: Parameters.SearchParameters.MatchBetweenRuns,
            //        silent: true,
            //        optionalPeriodicTablePath: GlobalVariables.ElementsLocation);

            //    if (flashLFQIdentifications.Any())
            //    {
            //        //make specific to protease
            //        var results = FlashLfqEngine.Run();

            //        if (Parameters.FlashLfqResults == null)
            //        {
            //            Parameters.FlashLfqResults = results;
            //        }
            //        else
            //        {
            //            Parameters.FlashLfqResults.MergeResultsWith(results);
            //        }
            //    }
            //}

            // get protein intensity back from FlashLFQ
            if (ProteinGroups != null && Parameters.FlashLfqResults != null)
            {
                foreach (var proteinGroup in ProteinGroups)
                {
                    proteinGroup.FilesForQuantification = spectraFileInfo;
                    proteinGroup.IntensitiesByFile = new Dictionary<SpectraFileInfo, double>();

                    foreach (var spectraFile in proteinGroup.FilesForQuantification)
                    {
                        if (Parameters.FlashLfqResults.ProteinGroups.TryGetValue(proteinGroup.ProteinGroupName, out var flashLfqProteinGroup))
                        {
                            proteinGroup.IntensitiesByFile.Add(spectraFile, flashLfqProteinGroup.GetIntensity(spectraFile));
                        }
                        else
                        {
                            proteinGroup.IntensitiesByFile.Add(spectraFile, 0);
                        }
                    }
                }
            }

            //If SILAC (Post-Quantification), compress the light/heavy protein group pairs into the same light protein group but different files
            //Create new files for each silac label and file so that "file 1" now becomes "file 1 (light)" and "file 1 (heavy)"
            //Change heavy residue into the light residue plus a string label ("PEPTIDEa" -> "PEPTIDEK(+8.014)")
            //This light to heavy conversion needs to happen for the flashLFQ peptides here, but can't for the psm peptides, which are constrained to the protein
            //i.e. pwsms currently don't have sequences; they have start/end residues and a protein sequence. We have to change the output sequences when they're created.
            if (Parameters.SilacLabels != null)
            {
                List<SilacLabel> silacLabels = Parameters.SilacLabels;

                //MAKE NEW RAW FILES
                //update number of spectra files to include a new file for each label*condition
                Dictionary<SpectraFileInfo, string> fileToLabelDictionary = new Dictionary<SpectraFileInfo, string>(); //figure out which file is which label, since some files will be only light and others only heavy. Key is file, value is the label string (label.MassDifference)
                Dictionary<SpectraFileInfo, SpectraFileInfo> labeledToUnlabeledFile = new Dictionary<SpectraFileInfo, SpectraFileInfo>(); //keep track of the heavy-to-light pairs. If multiple, looks like 3-1 and 2-1, but no 3-2 (only heavy to light, no heavy to heavy)
                List<SpectraFileInfo> silacSpectraFileInfo = new List<SpectraFileInfo>(); //new files

                //foreach existing file
                foreach (SpectraFileInfo originalFile in spectraFileInfo)
                {
                    //add the existing file as the light
                    silacSpectraFileInfo.Add(originalFile);
                    //foreach label, add a new file with the label
                    foreach (SilacLabel label in silacLabels)
                    {
                        SpectraFileInfo silacFile = new SpectraFileInfo(originalFile.FilenameWithoutExtension + "(" + label.OriginalAminoAcid + label.MassDifference + ")." + originalFile.FullFilePathWithExtension.Split('.').Last(), originalFile.Condition, originalFile.BiologicalReplicate, originalFile.TechnicalReplicate, originalFile.Fraction);
                        silacSpectraFileInfo.Add(silacFile);
                        fileToLabelDictionary[silacFile] = label.MassDifference;
                        labeledToUnlabeledFile[silacFile] = originalFile;
                    }
                }


                //UPDATE PROTEIN GROUPS
                //remove the heavy protein groups so that there are only light ones
                //add the intensities of the heavy groups into the newly created heavy SpectraFileInfos
                if (ProteinGroups != null) //if we did parsimony
                {
                    List<EngineLayer.ProteinGroup> silacProteinGroups = new List<EngineLayer.ProteinGroup>();

                    //foreach protein group (which has its own quant for each file)
                    foreach (EngineLayer.ProteinGroup proteinGroup in ProteinGroups)
                    {
                        proteinGroup.FilesForQuantification = silacSpectraFileInfo; //update fileinfo for the group
                                                                                    //grab the light groups. Using these light groups, find their heavy group pair(s), add them to the light group quant info, and then remove the heavy groups
                        if (silacProteinGroupMatcher.TryGetValue(proteinGroup.ProteinGroupName, out List<string> silacSubGroupNames)) //try to find the light protein groups. If it's not light, ignore it
                        {
                            //the out variable contains all the other heavy protein groups that were generated for this light protein group
                            //go through the files and see if any of them contain the same label. If not, put zeroes for those missing "files"
                            Dictionary<SpectraFileInfo, double> updatedIntensitiesByFile = proteinGroup.IntensitiesByFile; //light intensities

                            //go through all files (including "silac" files)
                            List<EngineLayer.ProteinGroup> subGroup = ProteinGroups.Where(x => silacSubGroupNames.Contains(x.ProteinGroupName)).ToList(); //find the protein groups where the accession contains "light" accession of the current protein group
                            foreach (SpectraFileInfo fileInfo in silacSpectraFileInfo) //for every file (light and heavy)
                            {
                                //if it doesn't have a value, then it's a silac file (light missing values still have a value "0")
                                if (!updatedIntensitiesByFile.ContainsKey(fileInfo))
                                {
                                    string labelSignature = fileToLabelDictionary[fileInfo]; //a string associated with a silac label
                                    EngineLayer.ProteinGroup foundGroup = subGroup.Where(x => x.Proteins.Any(y => y.Accession.Contains(labelSignature))).FirstOrDefault(); //get the protein groups containing this label
                                    updatedIntensitiesByFile[fileInfo] = foundGroup == null ? 0 : foundGroup.IntensitiesByFile[labeledToUnlabeledFile[fileInfo]]; //update the intensity for that label in the light group
                                }
                                //else do nothing. The light version is already in the dictionary
                            }

                            proteinGroup.IntensitiesByFile = updatedIntensitiesByFile;
                            silacProteinGroups.Add(proteinGroup);
                        }
                    }
                    ProteinGroups = silacProteinGroups; //update

                    //UPDATE FLASHLFQ PROTEINS
                    if (Parameters.FlashLfqResults != null) //can be null if nothing was quantified (all peptides are ambiguous)
                    {
                        Dictionary<string, FlashLFQ.ProteinGroup> flashLfqProteins = Parameters.FlashLfqResults.ProteinGroups; //dictionary of protein group names to protein groups
                                                                                                                               //if the protein group is a heavy protein group, get rid of it. We already accounted for it above.
                        var keys = flashLfqProteins.Keys.ToList();
                        foreach (string key in keys)
                        {
                            if (silacLabels.Any(x => key.Contains(x.MassDifference)))
                            {
                                flashLfqProteins.Remove(key);
                            }
                        }
                    }
                }

                //UPDATE FLASHLFQ SPECTRA FILES
                if (Parameters.FlashLfqResults != null) //can be null if nothing was quantified (all peptides are ambiguous)
                {
                    List<SpectraFileInfo> originalFiles = Parameters.FlashLfqResults.SpectraFiles;
                    foreach (SpectraFileInfo info in silacSpectraFileInfo)
                    {
                        if (!originalFiles.Contains(info))
                        {
                            originalFiles.Add(info);
                        }
                    }
                }

                //UPDATE PEPTIDE INFO
                //convert all psm/peptide/proteingroup sequences from the heavy label to the light label for output
                //We can do this for all of the FlashLFQ peptides/peaks, because they use string sequences.
                //We are unable to do this for Parameters.AllPsms, because they store proteins and start/end residues instead
                //for Psms, we need to convert during the writing.
                List<PeptideSpectralMatch> allPsms = Parameters.AllPsms;
                for (int i = 0; i < allPsms.Count; i++)
                {
                    allPsms[i].ResolveHeavySilacLabel(silacLabels, Parameters.SearchParameters.ModsToWriteSelection);
                }

                //Convert all lfqpeaks from heavy (a) to light (K+8.014) for output
                if (Parameters.FlashLfqResults != null) //can be null if nothing was quantified (all peptides are ambiguous)
                {
                    var lfqPeaks = Parameters.FlashLfqResults.Peaks;
                    List<SpectraFileInfo> peakKeys = lfqPeaks.Keys.ToList();
                    foreach (SpectraFileInfo key in peakKeys)
                    {
                        List<FlashLFQ.ChromatographicPeak> peaks = lfqPeaks[key];
                        for (int i = 0; i < peaks.Count; i++)
                        {
                            var peak = peaks[i];
                            List<Identification> identifications = new List<Identification>();
                            foreach (var id in peak.Identifications)
                            {
                                SilacLabel label = SilacConversions.GetRelevantLabelFromBaseSequence(id.BaseSequence, silacLabels);
                                HashSet<FlashLFQ.ProteinGroup> originalGroups = id.proteinGroups;
                                List<FlashLFQ.ProteinGroup> updatedGroups = new List<FlashLFQ.ProteinGroup>();
                                foreach (FlashLFQ.ProteinGroup group in originalGroups)
                                {
                                    string groupName = group.ProteinGroupName;
                                    if (label == null) //if light
                                    {
                                        updatedGroups.Add(group);
                                    }
                                    else
                                    {
                                        groupName = groupName.Replace(label.MassDifference, "");
                                        updatedGroups.Add(new FlashLFQ.ProteinGroup(groupName, group.GeneName, group.Organism));
                                    }
                                }

                                Identification updatedId = new Identification(
                                    id.fileInfo,
                                    SilacConversions.GetSilacLightBaseSequence(id.BaseSequence, label),
                                    SilacConversions.GetSilacLightFullSequence(id.ModifiedSequence, label),
                                    id.monoisotopicMass,
                                    id.ms2RetentionTimeInMinutes,
                                    id.precursorChargeState,
                                    updatedGroups,
                                    id.OptionalChemicalFormula,
                                    id.UseForProteinQuant
                                    );
                                identifications.Add(updatedId);
                            }
                            FlashLFQ.ChromatographicPeak updatedPeak = new FlashLFQ.ChromatographicPeak(identifications.First(), peak.IsMbrPeak, peak.SpectraFileInfo);
                            for (int j = 1; j < identifications.Count; j++) //add all the original identification
                            {
                                updatedPeak.MergeFeatureWith(new FlashLFQ.ChromatographicPeak(identifications[j], peak.IsMbrPeak, peak.SpectraFileInfo), FlashLfqEngine.Integrate);
                            }
                            updatedPeak.IsotopicEnvelopes = peak.IsotopicEnvelopes; //need to set isotopicEnevelopes, since the new identifications didn't have them.
                            updatedPeak.CalculateIntensityForThisFeature(FlashLfqEngine.Integrate); //needed to update info
                            peaks[i] = updatedPeak;
                        }
                    }

                    //convert all lfq peptides from heavy to light for output
                    Dictionary<string, Peptide> lfqPwsms = Parameters.FlashLfqResults.PeptideModifiedSequences;
                    List<string> pwsmKeys = lfqPwsms.Keys.ToList();
                    foreach (string key in pwsmKeys)
                    {
                        Peptide currentPeptide = lfqPwsms[key];
                        SilacLabel label = SilacConversions.GetRelevantLabelFromFullSequence(currentPeptide.Sequence, silacLabels);
                        if (label != null) //if it's a heavy peptide
                        {
                            lfqPwsms.Remove(key); //get rid of it
                                                  //update the light version
                            string lightSequence = SilacConversions.GetSilacLightFullSequence(currentPeptide.Sequence, label, false); //get the light sequence
                            List<SpectraFileInfo> heavyFiles = silacSpectraFileInfo.Where(x => x.FilenameWithoutExtension.Contains(label.MassDifference)).ToList(); //these are the heavy raw file names

                            //Find the light peptide (which has a value for the light datafile) and set the intensity for the heavy datafile from the current peptide
                            if (lfqPwsms.TryGetValue(lightSequence, out Peptide lightPeptide)) //this should always have a value, since we made replicas earlier, and yet it sometimes doesn't...
                            {
                                foreach (SpectraFileInfo heavyFile in heavyFiles)
                                {
                                    SpectraFileInfo lightFile = labeledToUnlabeledFile[heavyFile];
                                    lightPeptide.SetIntensity(heavyFile, currentPeptide.GetIntensity(lightFile));
                                    lightPeptide.SetDetectionType(heavyFile, currentPeptide.GetDetectionType(lightFile));
                                }
                            }
                            else //if there's no light, create a new entry for the heavy
                            {
                                //new peptide
                                Peptide updatedPeptide = new Peptide(lightSequence, currentPeptide.UseForProteinQuant);
                                //update the heavy info, set the light values to zero
                                foreach (SpectraFileInfo info in heavyFiles)
                                {
                                    updatedPeptide.SetIntensity(info, currentPeptide.GetIntensity(info));
                                    updatedPeptide.SetDetectionType(info, currentPeptide.GetDetectionType(info));
                                }

                                //set the other values to zero
                                List<SpectraFileInfo> otherInfo = silacSpectraFileInfo.Where(x => !heavyFiles.Contains(x)).ToList();
                                foreach (SpectraFileInfo info in otherInfo)
                                {
                                    updatedPeptide.SetIntensity(info, 0);
                                    updatedPeptide.SetDetectionType(info, DetectionType.NotDetected);
                                }
                                HashSet<FlashLFQ.ProteinGroup> originalGroups = currentPeptide.proteinGroups;
                                HashSet<FlashLFQ.ProteinGroup> updatedGroups = new HashSet<FlashLFQ.ProteinGroup>();
                                foreach (FlashLFQ.ProteinGroup group in originalGroups)
                                {
                                    string groupName = group.ProteinGroupName;
                                    groupName = groupName.Replace(label.MassDifference, "");
                                    updatedGroups.Add(new FlashLFQ.ProteinGroup(groupName, group.GeneName, group.Organism));
                                }
                                updatedPeptide.proteinGroups = updatedGroups;
                                lfqPwsms[updatedPeptide.Sequence] = updatedPeptide;
                            }
                        }
                    }
                }
            }
        }

        private void HistogramAnalysis()
        {
            if (Parameters.SearchParameters.DoHistogramAnalysis)
            {
                var limitedpsms_with_fdr = Parameters.AllPsms.Where(b => (b.FdrInfo.QValue <= 0.01)).ToList();
                if (limitedpsms_with_fdr.Any(b => !b.IsDecoy))
                {
                    Status("Running histogram analysis...", new List<string> { Parameters.SearchTaskId });
                    var myTreeStructure = new BinTreeStructure();
                    myTreeStructure.GenerateBins(limitedpsms_with_fdr, Parameters.SearchParameters.HistogramBinTolInDaltons);
                    var writtenFile = Path.Combine(Parameters.OutputFolder, "MassDifferenceHistogram.tsv");
                    WriteTree(myTreeStructure, writtenFile);
                    FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });
                }
            }
        }

        private void WritePsmResults()
        {
            Status("Writing PSM results...", Parameters.SearchTaskId);
            var FilteredPsmListForOutput = Parameters.AllPsms
                .Where(p => p.FdrInfo.QValue <= CommonParameters.QValueOutputFilter
                && p.FdrInfo.QValueNotch <= CommonParameters.QValueOutputFilter).ToList();

            if (!Parameters.SearchParameters.WriteDecoys)
            {
                FilteredPsmListForOutput.RemoveAll(b => b.IsDecoy);
            }
            if (!Parameters.SearchParameters.WriteContaminants)
            {
                FilteredPsmListForOutput.RemoveAll(b => b.IsContaminant);
            }

            // write PSMs
            string writtenFile = Path.Combine(Parameters.OutputFolder, "AllPSMs.psmtsv");
            WritePsmsToTsv(FilteredPsmListForOutput, writtenFile, Parameters.SearchParameters.ModsToWriteSelection);
            FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });

            // write PSMs for percolator
            writtenFile = Path.Combine(Parameters.OutputFolder, "AllPSMs_FormattedForPercolator.tsv");
            WritePsmsForPercolator(FilteredPsmListForOutput, writtenFile, CommonParameters.QValueOutputFilter);
            FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });

            // write summary text
            Parameters.SearchTaskResults.AddNiceText("All target PSMS within 1% FDR: " + Parameters.AllPsms.Count(a => a.FdrInfo.QValue <= 0.01 && !a.IsDecoy));
            if (Parameters.SearchParameters.DoParsimony)
            {
                Parameters.SearchTaskResults.AddNiceText("All target protein groups within 1% FDR: " + ProteinGroups.Count(b => b.QValue <= 0.01 && !b.IsDecoy)
                    + Environment.NewLine);
            }

            PsmsGroupedByFile = FilteredPsmListForOutput.GroupBy(p => p.FullFilePath);

            foreach (var file in PsmsGroupedByFile)
            {
                // write summary text
                var psmsForThisFile = file.ToList();
                string strippedFileName = Path.GetFileNameWithoutExtension(file.First().FullFilePath);

                Parameters.SearchTaskResults.AddNiceText("MS2 spectra in " + strippedFileName + ": " + Parameters.NumMs2SpectraPerFile[strippedFileName][0]);
                Parameters.SearchTaskResults.AddNiceText("Precursors fragmented in " + strippedFileName + ": " + Parameters.NumMs2SpectraPerFile[strippedFileName][1]);
                Parameters.SearchTaskResults.AddNiceText("Target PSMs within 1% FDR in " + strippedFileName + ": " + psmsForThisFile.Count(a => a.FdrInfo.QValue <= 0.01 && !a.IsDecoy));

                // writes all individual spectra file search results to subdirectory
                if (Parameters.CurrentRawFileList.Count > 1)
                {
                    // create individual files subdirectory
                    Directory.CreateDirectory(Parameters.IndividualResultsOutputFolder);

                    // write PSMs
                    writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + "_PSMs.psmtsv");
                    WritePsmsToTsv(psmsForThisFile, writtenFile, Parameters.SearchParameters.ModsToWriteSelection);
                    FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", file.First().FullFilePath });

                    // write PSMs for percolator
                    writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + "_PSMsFormattedForPercolator.tsv");
                    WritePsmsForPercolator(psmsForThisFile, writtenFile, CommonParameters.QValueOutputFilter);
                    FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", file.First().FullFilePath });
                }
            }
        }

        private void WriteProteinResults()
        {
            if (Parameters.SearchParameters.DoParsimony)
            {
                // write protein groups to tsv
                string writtenFile = Path.Combine(Parameters.OutputFolder, "AllProteinGroups.tsv");
                WriteProteinGroupsToTsv(ProteinGroups, writtenFile, new List<string> { Parameters.SearchTaskId }, CommonParameters.QValueOutputFilter);

                // write all individual file results to subdirectory
                // local protein fdr, global parsimony, global psm fdr
                if (Parameters.CurrentRawFileList.Count > 1 || Parameters.SearchParameters.WriteMzId || Parameters.SearchParameters.WritePepXml)
                {
                    Directory.CreateDirectory(Parameters.IndividualResultsOutputFolder);

                    foreach (var fullFilePath in PsmsGroupedByFile.Select(v => v.Key))
                    {
                        string strippedFileName = Path.GetFileNameWithoutExtension(fullFilePath);

                        List<PeptideSpectralMatch> psmsForThisFile = PsmsGroupedByFile.Where(p => p.Key == fullFilePath).SelectMany(g => g).ToList();
                        var subsetProteinGroupsForThisFile = ProteinGroups.Select(p => p.ConstructSubsetProteinGroup(fullFilePath)).ToList();

                        ProteinScoringAndFdrResults subsetProteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(subsetProteinGroupsForThisFile, psmsForThisFile,
                            Parameters.SearchParameters.NoOneHitWonders, Parameters.SearchParameters.ModPeptidesAreDifferent,
                            false, CommonParameters, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }).Run();

                        subsetProteinGroupsForThisFile = subsetProteinScoringAndFdrResults.SortedAndScoredProteinGroups;

                        Parameters.SearchTaskResults.AddNiceText("Target protein groups within 1 % FDR in " + strippedFileName + ": " + subsetProteinGroupsForThisFile.Count(b => b.QValue <= 0.01 && !b.IsDecoy));

                        // write individual spectra file protein groups results to tsv
                        if (Parameters.CurrentRawFileList.Count > 1)
                        {
                            writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + "_ProteinGroups.tsv");
                            WriteProteinGroupsToTsv(subsetProteinGroupsForThisFile, writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }, CommonParameters.QValueOutputFilter);
                        }

                        // write mzID
                        if (Parameters.SearchParameters.WriteMzId)
                        {
                            Status("Writing mzID...", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });

                            var mzidFilePath = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + ".mzID");
                            MzIdentMLWriter.WriteMzIdentMl(psmsForThisFile, subsetProteinGroupsForThisFile, Parameters.VariableModifications, Parameters.FixedModifications, Parameters.SilacLabels,
                                new List<Protease> { CommonParameters.DigestionParams.Protease }, CommonParameters.QValueOutputFilter, CommonParameters.ProductMassTolerance,
                                CommonParameters.PrecursorMassTolerance, CommonParameters.DigestionParams.MaxMissedCleavages, mzidFilePath);

                            FinishedWritingFile(mzidFilePath, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });
                        }

                        // write pepXML
                        if (Parameters.SearchParameters.WritePepXml)
                        {
                            Status("Writing pepXML...", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });

                            var pepXMLFilePath = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + ".pep.XML");
                            PepXMLWriter.WritePepXml(psmsForThisFile, Parameters.DatabaseFilenameList, Parameters.VariableModifications, Parameters.FixedModifications,
                                CommonParameters, pepXMLFilePath, CommonParameters.QValueOutputFilter);

                            FinishedWritingFile(pepXMLFilePath, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });
                        }

                        ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }));
                    }
                }
            }
        }

        private void WriteQuantificationResults()
        {
            if (Parameters.SearchParameters.DoQuantification && Parameters.FlashLfqResults != null)
            {
                // write peaks
                WritePeakQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, "AllQuantifiedPeaks", new List<string> { Parameters.SearchTaskId });

                // write peptide quant results
                WritePeptideQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, "AllQuantifiedPeptides", new List<string> { Parameters.SearchTaskId });

                // write individual results
                if (Parameters.CurrentRawFileList.Count > 1)
                {
                    foreach (var file in Parameters.FlashLfqResults.Peaks)
                    {
                        WritePeakQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.IndividualResultsOutputFolder,
                            file.Key.FilenameWithoutExtension + "_QuantifiedPeaks", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", file.Key.FullFilePathWithExtension });
                    }
                }
            }
        }

        private void WritePrunedDatabase()
        {
            if (Parameters.SearchParameters.WritePrunedDatabase)
            {
                Status("Writing Pruned Database...", new List<string> { Parameters.SearchTaskId });
                HashSet<Modification> modificationsToWriteIfBoth = new HashSet<Modification>();
                HashSet<Modification> modificationsToWriteIfInDatabase = new HashSet<Modification>();
                HashSet<Modification> modificationsToWriteIfObserved = new HashSet<Modification>();

                var confidentPsms = Parameters.AllPsms.Where(b => b.FdrInfo.QValueNotch <= 0.01 && b.FdrInfo.QValue <= 0.01 && !b.IsDecoy && b.BaseSequence != null).ToList();
                var proteinToConfidentBaseSequences = new Dictionary<Protein, List<PeptideWithSetModifications>>();

                // associate all confident PSMs with all possible proteins they could be digest products of (before or after parsimony)
                foreach (PeptideSpectralMatch psm in confidentPsms)
                {
                    var myPepsWithSetMods = psm.BestMatchingPeptides.Select(p => p.Peptide);

                    foreach (var peptide in myPepsWithSetMods)
                    {
                        if (proteinToConfidentBaseSequences.TryGetValue(peptide.Protein.NonVariantProtein, out var myPepList))
                        {
                            myPepList.Add(peptide);
                        }
                        else
                        {
                            proteinToConfidentBaseSequences.Add(peptide.Protein.NonVariantProtein, new List<PeptideWithSetModifications> { peptide });
                        }
                    }
                }

                // Add user mod selection behavours to Pruned DB
                foreach (var modType in Parameters.SearchParameters.ModsToWriteSelection)
                {
                    foreach (Modification mod in GlobalVariables.AllModsKnown.Where(b => b.ModificationType.Equals(modType.Key)))
                    {
                        if (modType.Value == 1) // Write if observed and in database
                        {
                            modificationsToWriteIfBoth.Add(mod);
                        }
                        if (modType.Value == 2) // Write if in database
                        {
                            modificationsToWriteIfInDatabase.Add(mod);
                        }
                        if (modType.Value == 3) // Write if observed
                        {
                            modificationsToWriteIfObserved.Add(mod);
                        }
                    }
                }

                //generates dictionary of proteins with only localized modifications
                var ModPsms = Parameters.AllPsms.Where(b => b.FdrInfo.QValueNotch <= 0.01 && b.FdrInfo.QValue <= 0.01 && !b.IsDecoy && b.FullSequence != null).ToList();
                var proteinToConfidentModifiedSequences = new Dictionary<Protein, List<PeptideWithSetModifications>>();

                foreach (PeptideSpectralMatch psm in ModPsms)
                {
                    var myPepsWithSetMods = psm.BestMatchingPeptides.Select(p => p.Peptide);

                    foreach (var peptide in myPepsWithSetMods)
                    {
                        if (proteinToConfidentModifiedSequences.TryGetValue(peptide.Protein.NonVariantProtein, out var myPepList))
                        {
                            myPepList.Add(peptide);
                        }
                        else
                        {
                            proteinToConfidentModifiedSequences.Add(peptide.Protein.NonVariantProtein, new List<PeptideWithSetModifications> { peptide });
                        }
                    }
                }

                // mods included in pruned database will only be confidently localized mods (peptide's FullSequence != null)
                foreach (var nonVariantProtein in Parameters.ProteinList.Select(p => p.NonVariantProtein).Distinct())
                {
                    if (!nonVariantProtein.IsDecoy)
                    {
                        proteinToConfidentModifiedSequences.TryGetValue(nonVariantProtein, out var psms);
                        HashSet<(int, Modification, SequenceVariation)> modsObservedOnThisProtein = new HashSet<(int, Modification, SequenceVariation)>(); // sequence variant is null if mod is not on a variant
                        foreach (PeptideWithSetModifications psm in psms ?? new List<PeptideWithSetModifications>())
                        {
                            foreach (var idxModKV in psm.AllModsOneIsNterminus)
                            {
                                int proteinIdx = GetOneBasedIndexInProtein(idxModKV.Key, psm);
                                SequenceVariation relevantVariant = psm.Protein.AppliedSequenceVariations.FirstOrDefault(sv => VariantApplication.IsSequenceVariantModification(sv, proteinIdx));
                                SequenceVariation unappliedVariant =
                                    relevantVariant == null ? null : // it's not a sequence variant mod
                                    psm.Protein.SequenceVariations.FirstOrDefault(sv => sv.Description != null && sv.Description.Equals(relevantVariant.Description));
                                modsObservedOnThisProtein.Add((VariantApplication.RestoreModificationIndex(psm.Protein, proteinIdx), idxModKV.Value, unappliedVariant));
                            }
                        }

                        IDictionary<(SequenceVariation, int), List<Modification>> modsToWrite = new Dictionary<(SequenceVariation, int), List<Modification>>();

                        //Add if observed (regardless if in database)
                        foreach (var observedMod in modsObservedOnThisProtein)
                        {
                            var tempMod = observedMod.Item2;

                            if (modificationsToWriteIfObserved.Contains(tempMod))
                            {
                                var svIdxKey = (observedMod.Item3, observedMod.Item1);
                                if (!modsToWrite.ContainsKey(svIdxKey))
                                {
                                    modsToWrite.Add(svIdxKey, new List<Modification> { observedMod.Item2 });
                                }
                                else
                                {
                                    modsToWrite[svIdxKey].Add(observedMod.Item2);
                                }
                            }
                        }

                        // Add modification if in database (two cases: always or if observed)
                        foreach (var modkv in nonVariantProtein.OneBasedPossibleLocalizedModifications)
                        {
                            foreach (var mod in modkv.Value)
                            {
                                //Add if always In Database or if was observed and in database and not set to not include
                                if (modificationsToWriteIfInDatabase.Contains(mod) ||
                                    (modificationsToWriteIfBoth.Contains(mod) && modsObservedOnThisProtein.Contains((modkv.Key, mod, null))))
                                {
                                    if (!modsToWrite.ContainsKey((null, modkv.Key)))
                                    {
                                        modsToWrite.Add((null, modkv.Key), new List<Modification> { mod });
                                    }
                                    else
                                    {
                                        modsToWrite[(null, modkv.Key)].Add(mod);
                                    }
                                }
                            }
                        }

                        // Add variant modification if in database (two cases: always or if observed)
                        foreach (SequenceVariation sv in nonVariantProtein.SequenceVariations)
                        {
                            foreach (var modkv in sv.OneBasedModifications)
                            {
                                foreach (var mod in modkv.Value)
                                {
                                    //Add if always In Database or if was observed and in database and not set to not include
                                    if (modificationsToWriteIfInDatabase.Contains(mod) ||
                                        (modificationsToWriteIfBoth.Contains(mod) && modsObservedOnThisProtein.Contains((modkv.Key, mod, sv))))
                                    {
                                        if (!modsToWrite.ContainsKey((sv, modkv.Key)))
                                        {
                                            modsToWrite.Add((sv, modkv.Key), new List<Modification> { mod });
                                        }
                                        else
                                        {
                                            modsToWrite[(sv, modkv.Key)].Add(mod);
                                        }
                                    }
                                }
                            }
                        }

                        if (proteinToConfidentBaseSequences.ContainsKey(nonVariantProtein.NonVariantProtein))
                        {
                            // adds confidently localized and identified mods
                            nonVariantProtein.OneBasedPossibleLocalizedModifications.Clear();
                            foreach (var kvp in modsToWrite.Where(kv => kv.Key.Item1 == null))
                            {
                                nonVariantProtein.OneBasedPossibleLocalizedModifications.Add(kvp.Key.Item2, kvp.Value);
                            }
                            foreach (var sv in nonVariantProtein.SequenceVariations)
                            {
                                sv.OneBasedModifications.Clear();
                                foreach (var kvp in modsToWrite.Where(kv => kv.Key.Item1 != null && kv.Key.Item1.Equals(sv)))
                                {
                                    sv.OneBasedModifications.Add(kvp.Key.Item2, kvp.Value);
                                }
                            }
                        }
                    }
                }

                //writes all proteins
                if (Parameters.DatabaseFilenameList.Any(b => !b.IsContaminant))
                {
                    string outputXMLdbFullName = Path.Combine(Parameters.OutputFolder, string.Join("-", Parameters.DatabaseFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "pruned.xml");
                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), Parameters.ProteinList.Select(p => p.NonVariantProtein).Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);
                    FinishedWritingFile(outputXMLdbFullName, new List<string> { Parameters.SearchTaskId });
                }
                if (Parameters.DatabaseFilenameList.Any(b => b.IsContaminant))
                {
                    string outputXMLdbFullNameContaminants = Path.Combine(Parameters.OutputFolder, string.Join("-", Parameters.DatabaseFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "pruned.xml");
                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), Parameters.ProteinList.Select(p => p.NonVariantProtein).Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);
                    FinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { Parameters.SearchTaskId });
                }

                //writes only detected proteins
                if (Parameters.DatabaseFilenameList.Any(b => !b.IsContaminant))
                {
                    string outputXMLdbFullName = Path.Combine(Parameters.OutputFolder, string.Join("-", Parameters.DatabaseFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "proteinPruned.xml");
                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinToConfidentBaseSequences.Keys.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);
                    FinishedWritingFile(outputXMLdbFullName, new List<string> { Parameters.SearchTaskId });
                }
                if (Parameters.DatabaseFilenameList.Any(b => b.IsContaminant))
                {
                    string outputXMLdbFullNameContaminants = Path.Combine(Parameters.OutputFolder, string.Join("-", Parameters.DatabaseFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "proteinPruned.xml");
                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinToConfidentBaseSequences.Keys.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);
                    FinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { Parameters.SearchTaskId });
                }
            }
        }

        private void WritePeptideResults()
        {
            Status("Writing peptide results...", Parameters.SearchTaskId);

            // write best (highest-scoring) PSM per peptide
            string writtenFile = Path.Combine(Parameters.OutputFolder, "AllPeptides.psmtsv");
            List<PeptideSpectralMatch> peptides = Parameters.AllPsms.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();
            
            new FdrAnalysisEngine(peptides, Parameters.NumNotches, CommonParameters, new List<string> { Parameters.SearchTaskId }, "Peptide").Run();

            if (!Parameters.SearchParameters.WriteDecoys)
            {
                peptides.RemoveAll(b => b.IsDecoy);
            }
            if (!Parameters.SearchParameters.WriteContaminants)
            {
                peptides.RemoveAll(b => b.IsContaminant);
            }
            peptides.RemoveAll(p => p.FdrInfo.QValue > CommonParameters.QValueOutputFilter);

            WritePsmsToTsv(peptides, writtenFile, Parameters.SearchParameters.ModsToWriteSelection);
            FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });

            Parameters.SearchTaskResults.AddNiceText("All target peptides within 1% FDR: " + peptides.Count(a => a.FdrInfo.QValue <= 0.01 && !a.IsDecoy));

            foreach (var file in PsmsGroupedByFile)
            {
                // write summary text
                var psmsForThisFile = file.ToList();
                string strippedFileName = Path.GetFileNameWithoutExtension(file.First().FullFilePath);
                var peptidesForFile = psmsForThisFile.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();
                new FdrAnalysisEngine(peptidesForFile, Parameters.NumNotches, CommonParameters, new List<string> { Parameters.SearchTaskId }, "Peptide").Run();
                Parameters.SearchTaskResults.AddNiceText("Target peptides within 1% FDR in " + strippedFileName + ": " + peptidesForFile.Count(a => a.FdrInfo.QValue <= 0.01 && !a.IsDecoy) + Environment.NewLine);

                // writes all individual spectra file search results to subdirectory
                if (Parameters.CurrentRawFileList.Count > 1)
                {
                    // create individual files subdirectory
                    Directory.CreateDirectory(Parameters.IndividualResultsOutputFolder);

                    // write best (highest-scoring) PSM per peptide
                    writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + "_Peptides.psmtsv");
                    WritePsmsToTsv(peptidesForFile, writtenFile, Parameters.SearchParameters.ModsToWriteSelection);
                    FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", file.First().FullFilePath });
                }
            }
        }

        private static int GetOneBasedIndexInProtein(int oneIsNterminus, PeptideWithSetModifications peptideWithSetModifications)
        {
            if (oneIsNterminus == 1)
            {
                return peptideWithSetModifications.OneBasedStartResidueInProtein;
            }
            if (oneIsNterminus == peptideWithSetModifications.Length + 2)
            {
                return peptideWithSetModifications.OneBasedEndResidueInProtein;
            }
            return peptideWithSetModifications.OneBasedStartResidueInProtein + oneIsNterminus - 2;
        }

        private static void WriteTree(BinTreeStructure myTreeStructure, string writtenFile)
        {
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tArea 0.01t\tArea 0.255\tFracLocalizeableTarget\tMine\tUnimodID\tUnimodFormulas\tUnimodDiffs\tAA\tCombos\tModsInCommon\tAAsInCommon\tResidues\tprotNtermLocFrac\tpepNtermLocFrac\tpepCtermLocFrac\tprotCtermLocFrac\tFracWithSingle\tOverlappingFrac\tMedianLength\tUniprot");
                foreach (Bin bin in myTreeStructure.FinalBins.OrderByDescending(b => b.Count))
                {
                    output.WriteLine(bin.MassShift.ToString("F4", CultureInfo.InvariantCulture)
                        + "\t" + bin.Count.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountDecoy.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.LocalizeableTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget - bin.LocalizeableTarget).ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.Count == 0 ? double.NaN : (double)bin.CountDecoy / bin.Count).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.01))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.255))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget == 0 ? double.NaN : (double)bin.LocalizeableTarget / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.Mine
                        + "\t" + bin.UnimodId
                        + "\t" + bin.UnimodFormulas
                        + "\t" + bin.UnimodDiffs
                        + "\t" + bin.AA
                        + "\t" + bin.Combos
                        + "\t" + string.Join(",", bin.ModsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.AAsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.ResidueCount.OrderByDescending(b => b.Value).Select(b => b.Key + ":" + b.Value))
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.ProtNlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.PepNlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.PepClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.ProtClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.FracWithSingle).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + ((double)bin.Overlapping / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.MedianLength).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.UniprotID);
                }
            }
        }

        private void WritePsmsForPercolator(List<PeptideSpectralMatch> psmList, string writtenFileForPercolator, double qValueCutoff)
        {
            using (StreamWriter output = new StreamWriter(writtenFileForPercolator))
            {
                output.WriteLine("SpecId\tLabel\tScanNr\tF1\tF2\tPeptide\tProteins");
                output.WriteLine("DefaultDirection\t-\t-\t1\t1\t\t");
                for (int i = 0; i < psmList.Count; i++)
                {
                    var psm = psmList[i];

                    if (psm.FdrInfo.QValue > qValueCutoff || psm.FdrInfo.QValueNotch > qValueCutoff)
                    {
                        continue;
                    }

                    output.Write(i.ToString());
                    output.Write('\t' + (psm.IsDecoy ? -1 : 1).ToString());
                    output.Write('\t' + psm.ScanNumber.ToString());

                    // Features
                    output.Write('\t' + string.Join("\t", psm.Features));

                    // HACKY: Ignores all ambiguity
                    var pwsm = psm.BestMatchingPeptides.First().Peptide;

                    output.Write('\t' + (pwsm.PreviousAminoAcid + "." + pwsm.FullSequence + "." + pwsm.NextAminoAcid).ToString());
                    output.Write('\t' + (pwsm.Protein.Accession).ToString());
                    output.WriteLine();
                }
            }
        }

        private void WriteProteinGroupsToTsv(List<EngineLayer.ProteinGroup> proteinGroups, string filePath, List<string> nestedIds, double qValueCutoff)
        {
            if (proteinGroups != null && proteinGroups.Any())
            {
                using (StreamWriter output = new StreamWriter(filePath))
                {
                    output.WriteLine(proteinGroups.First().GetTabSeparatedHeader());
                    for (int i = 0; i < proteinGroups.Count; i++)
                    {
                        if ((!Parameters.SearchParameters.WriteDecoys && proteinGroups[i].IsDecoy) ||
                            (!Parameters.SearchParameters.WriteContaminants && proteinGroups[i].IsContaminant))
                        {
                            continue;
                        }
                        else if (proteinGroups[i].QValue <= qValueCutoff)
                        {
                            output.WriteLine(proteinGroups[i]);
                        }
                    }
                }

                FinishedWritingFile(filePath, nestedIds);
            }
        }

        private void WritePeptideQuantificationResultsToTsv(FlashLfqResults flashLFQResults, string outputFolder, string fileName, List<string> nestedIds)
        {
            var fullSeqPath = Path.Combine(outputFolder, fileName + ".tsv");

            flashLFQResults.WriteResults(null, fullSeqPath, null);

            FinishedWritingFile(fullSeqPath, nestedIds);
        }

        private void WritePeakQuantificationResultsToTsv(FlashLfqResults flashLFQResults, string outputFolder, string fileName, List<string> nestedIds)
        {
            var peaksPath = Path.Combine(outputFolder, fileName + ".tsv");

            flashLFQResults.WriteResults(peaksPath, null, null);

            FinishedWritingFile(peaksPath, nestedIds);
        }
    }
}