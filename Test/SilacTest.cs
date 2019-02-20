﻿using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class SilacTest
    {
        [Test]
        public static void TestSilacQuantification()
        {
            Residue heavyLysine = new Residue("a", 'a', "a", Chemistry.ChemicalFormula.ParseFormula("C{13}6H12N{15}2O"), ModificationSites.All); //+8 lysine
            Residue lightLysine = Residue.GetResidue('K');
            List<SilacLabel> silacLabels = new List<Proteomics.SilacLabel>
            {
                new SilacLabel(lightLysine.Letter, lightLysine.Letter, heavyLysine.ThisChemicalFormula.Formula, heavyLysine.MonoisotopicMass - lightLysine.MonoisotopicMass)
            };
            SearchTask task = new SearchTask
            {
                CommonParameters = new CommonParameters(),
                SearchParameters = new SearchParameters
                {
                    SilacLabels = silacLabels
                }
            };

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("1", task) };

            PeptideWithSetModifications lightPeptide = new PeptideWithSetModifications("PEPTIDEK", new Dictionary<string, Proteomics.Modification>());

            List<double> massDifferences = new List<double> { heavyLysine.MonoisotopicMass - lightLysine.MonoisotopicMass };
            MsDataFile myMsDataFile1 = new TestDataFile(lightPeptide, massDifferences);
            string mzmlName = @"silac.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlName, false);

            string xmlName = "SilacDb.xml";
            {
                Protein theProtein = new Protein("PEPTIDEK", "accession1");
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { theProtein }, xmlName);
            }

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSilac");
            Directory.CreateDirectory(outputFolder);
            var theStringResult = task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1").ToString();

            Assert.IsTrue(theStringResult.Contains("All target PSMS within 1% FDR: 1")); //it's not a psm, it's a MBR feature

            //test proteins
            string[] output = File.ReadAllLines(TestContext.CurrentContext.TestDirectory + @"/TestSilac/AllProteinGroups.tsv");
            Assert.AreEqual(output.Length, 2);
            Assert.IsTrue(output[0].Contains("Intensity_silac\tIntensity_silac(K+8.014)")); //test that two files were made
            Assert.IsTrue(output[1].Contains("875000\t437500")); //test the heavy intensity is half that of the light (per the raw file)

            //test peptides
            output = File.ReadAllLines(TestContext.CurrentContext.TestDirectory + @"/TestSilac/AllQuantifiedPeptides.tsv");
            Assert.AreEqual(output.Length, 3);
            Assert.IsTrue(output[1].Contains("PEPTIDEK\taccession1\t"));//test the accession was not modified
            Assert.IsTrue(output[1].Contains("875000")); //test intensity
            Assert.IsTrue(output[2].Contains("PEPTIDEK(+8.014)\taccession1\t")); //test the accession was not modified
            Assert.IsTrue(output[2].Contains("437500")); //test intensity

            //test peaks
            output = File.ReadAllLines(TestContext.CurrentContext.TestDirectory + @"/TestSilac/AllQuantifiedPeaks.tsv");
            Assert.AreEqual(output.Length, 3);
            Assert.IsTrue(output[1].Contains("silac\t")); //test the filename was NOT modified (it was for proteins, but we don't want it for peptides)
            Assert.IsTrue(output[2].Contains("silac\t"));//test the filename was NOT modified (it was for proteins, but we don't want it for peptides)
            Assert.IsTrue(output[1].Contains("PEPTIDEK\t")); //test light sequence was not modified
            Assert.IsTrue(output[2].Contains("PEPTIDEK(+8.014)\t")); //test heavy sequence was output correctly (do NOT want "PEPTIDEa")
            Assert.IsTrue(output[1].Contains("927.45")); //test light mass
            Assert.IsTrue(output[2].Contains("935.46")); //test heavy mass

            //delete files
            Directory.Delete(outputFolder, true);
            File.Delete(xmlName);
            File.Delete(mzmlName);
        }
    }
}