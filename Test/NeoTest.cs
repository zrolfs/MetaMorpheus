using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    internal static class NeoTest
    {
        #region Public Methods

        [Test]
        public static void TestNeoAll()
        {
            NeoSearchTask ye5 = new NeoSearchTask()
            {
                CommonParameters = new CommonParameters
                {
                    ConserveMemory = false,
                    ScoreCutoff = 1,
                    TotalPartitions = 2,
                    TaskDescriptor = "NeoSearchTask",
                    DigestionParams = new DigestionParams
                    {
                        MaxMissedCleavages = 1,
                        MinPeptideLength = 1,
                        InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain
                    }
                },
                NeoParameters = new NeoParameters
                {
                    Calibrate = false
                }
            };
            List<MetaMorpheusTask> taskList = NeoLoadTomls.LoadTomls(ye5);

            #region Setup tasks

            foreach (var modFile in Directory.GetFiles(@"Mods"))
                GlobalVariables.AddMods(PtmListLoader.ReadModsFromFile(modFile));

            #endregion Setup tasks

            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => taskList[0].CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => taskList[0].CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();

            // Generate data for files
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1");
            Protein DecoyProtein = new Protein("MEHTNAK", "accessiond");
            Protein CisSpliceProtein = new Protein("PEPTIANTHE", "spliced");
            Protein TransSpliceProtein = new Protein("AACNNPEPTIDE", "spliced");
            var digestedList = ParentProtein.Digest(taskList[0].CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
            PeptideWithSetModifications decoyPep = DecoyProtein.Digest(taskList[0].CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0];
            PeptideWithSetModifications cissplicePep = CisSpliceProtein.Digest(taskList[0].CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0];
            PeptideWithSetModifications transsplicePep = TransSpliceProtein.Digest(taskList[0].CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0];
            Assert.AreEqual(5, digestedList.Count);

            var dictHere = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            dictHere.Add(3, new List<Modification> { new ModificationWithMass("21", null, motif, TerminusLocalization.Any, 21.981943) });
            Protein ParentProteinToNotInclude = new Protein("MPEPTIDEK", "accession2", "organism", new List<Tuple<string, string>>(), dictHere);
            digestedList = ParentProteinToNotInclude.Digest(taskList[0].CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { digestedList[0], digestedList[1], digestedList[2], digestedList[3], decoyPep, cissplicePep, transsplicePep });

            Protein proteinWithChain = new Protein("MAACNNNCAA", "accession3", "organism", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new List<ProteolysisProduct> { new ProteolysisProduct(4, 8, "chain") }, "name2", "fullname2");

            #region Write the files

            string mzmlName = @"ok.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
            string xmlName = "okk.xml";
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { ParentProtein, proteinWithChain }, xmlName);
            DbForTask db = new DbForTask(xmlName, false);
            #endregion Write the files

            List<(string, MetaMorpheusTask)> taskList2 = new List<(string, MetaMorpheusTask)>();
            foreach(MetaMorpheusTask task in taskList)
            {
                taskList2.Add(("Task" + (taskList2.Count + 1) + "-" + task.CommonParameters.TaskDescriptor, task));
            }
            var engine = new EverythingRunnerEngine(taskList2, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false) }, Environment.CurrentDirectory);
            engine.Run();
            Assert.IsTrue(Directory.Exists(@"Task13-NeoSearchTask"));
        }

        [Test]
        public static void TestNeoMultipleRawFiles()
        {
            NeoSearchTask ye5 = new NeoSearchTask()
            {
                CommonParameters = new CommonParameters
                {
                    ConserveMemory = false,
                    ScoreCutoff = 1,
                    TotalPartitions = 2,
                    TaskDescriptor = "NeoSearchTask",
                    DigestionParams = new DigestionParams
                    {
                        MaxMissedCleavages = 1,
                        MinPeptideLength = 1,
                        InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain
                    }
                },
                NeoParameters = new NeoParameters
                {
                    Calibrate = false
                }
            };
            List<MetaMorpheusTask> taskList = NeoLoadTomls.LoadTomls(ye5);

            #region Setup tasks

            foreach (var modFile in Directory.GetFiles(@"Mods"))
                GlobalVariables.AddMods(PtmListLoader.ReadModsFromFile(modFile));

            #endregion Setup tasks

            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => taskList[0].CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => taskList[0].CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();

            // Generate data for files
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1");
            Protein DecoyProtein = new Protein("MEHTNAK", "accessiond");
            Protein CisSpliceProtein = new Protein("PEPTIANTHE", "spliced");
            Protein TransSpliceProtein = new Protein("AACNNPEPTIDE", "spliced");
            var digestedList = ParentProtein.Digest(taskList[0].CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
            PeptideWithSetModifications decoyPep = DecoyProtein.Digest(taskList[0].CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0];
            PeptideWithSetModifications cissplicePep = CisSpliceProtein.Digest(taskList[0].CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0];
            PeptideWithSetModifications transsplicePep = TransSpliceProtein.Digest(taskList[0].CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0];
            Assert.AreEqual(5, digestedList.Count);

            var dictHere = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            dictHere.Add(3, new List<Modification> { new ModificationWithMass("21", null, motif, TerminusLocalization.Any, 21.981943) });
            Protein ParentProteinToNotInclude = new Protein("MPEPTIDEK", "accession2", "organism", new List<Tuple<string, string>>(), dictHere);
            digestedList = ParentProteinToNotInclude.Digest(taskList[0].CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { digestedList[0], digestedList[1], digestedList[2], digestedList[3], decoyPep, cissplicePep, transsplicePep });

            Protein proteinWithChain = new Protein("MAACNNNCAA", "accession3", "organism", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new List<ProteolysisProduct> { new ProteolysisProduct(4, 8, "chain") }, "name2", "fullname2");

            #region Write the files

            string mzmlName = @"ok.mzML"; 
            string mzmlName2 = @"ok2.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName2, false);
            string xmlName = "okk.xml";
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { ParentProtein, proteinWithChain }, xmlName);
            DbForTask db = new DbForTask(xmlName, false);
            #endregion Write the files

            List<(string, MetaMorpheusTask)> taskList2 = new List<(string, MetaMorpheusTask)>();
            foreach (MetaMorpheusTask task in taskList)
            {
                taskList2.Add(("Task" + (taskList2.Count + 1) + "-" + task.CommonParameters.TaskDescriptor, task));
            }
            var engine = new EverythingRunnerEngine(taskList2, new List<string> { mzmlName, mzmlName2 }, new List<DbForTask> { new DbForTask(xmlName, false) }, Environment.CurrentDirectory);
            engine.Run();
            Assert.IsTrue(Directory.Exists(@"Task13-NeoSearchTask"));
        }

        [Test]
        public static void TestNeoLoadTomls()
        {
            const int numMaxLengthToBeFollowed = 572;
            const int numMinLengthToBeFollowed = 34;
            const int numMaxMissedCleavages = 62;
            const string proteaseName = "chymotrypsin (don't cleave before proline)";
            const int scoreCutoff = -3;
            IDigestionParams tempDigParams = new DigestionParams
            {
                MinPeptideLength = numMinLengthToBeFollowed,
                MaxPeptideLength = numMaxLengthToBeFollowed,
                Protease = GlobalVariables.ProteaseDictionary[proteaseName],
                MaxMissedCleavages = numMaxMissedCleavages
            };
            NeoSearchTask ye5 = new NeoSearchTask()
            {

                CommonParameters = new CommonParameters
                {
                    DigestionParams = tempDigParams,
                    DoPrecursorDeconvolution = true,
                    PrecursorMassTolerance = null,
                    ProductMassTolerance = null,
                    BIons = true,
                    CIons = true,
                    YIons = true,
                    ZdotIons = true,
                    ScoreCutoff = scoreCutoff
                }
            };
            List<MetaMorpheusTask> taskList = NeoLoadTomls.LoadTomls(ye5);
            Assert.IsTrue(taskList.Count == 14);
            for(int i=0; i<taskList.Count; i++)
            {
                MetaMorpheusTask task = taskList[i];
                Assert.IsTrue(task.CommonParameters.DigestionParams.MinPeptideLength == numMinLengthToBeFollowed 
                    || task.CommonParameters.DigestionParams.Protease.Name.Contains("single")
                    || task.CommonParameters.DigestionParams.Protease.Name == "top-down");
                Assert.IsTrue(task.CommonParameters.DigestionParams.MaxPeptideLength == numMaxLengthToBeFollowed);
                Assert.IsTrue(task.CommonParameters.DigestionParams.MaxMissedCleavages == numMaxMissedCleavages);
                Assert.IsTrue(task.CommonParameters.DigestionParams.Protease.Name == proteaseName
                    || task.CommonParameters.DigestionParams.Protease.Name.Contains("single")
                    || task.CommonParameters.DigestionParams.Protease.Name == "top-down");
                Assert.IsTrue(task.CommonParameters.BIons || task.CommonParameters.YIons);
                Assert.IsTrue(task.CommonParameters.CIons || task.CommonParameters.ZdotIons);
                Assert.IsTrue(task.CommonParameters.ScoreCutoff == scoreCutoff);
            }
        }

        [Test]
        public static void TestNeoLoadNC()
        {

        }



        #endregion Public Methods

    }
}
