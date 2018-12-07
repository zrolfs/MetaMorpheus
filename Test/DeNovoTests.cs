using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.DeNovoSequencing;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    internal static class DeNovoTests
    {
        [Test]
        public static void TestDeNovo()
        {
            //theoretical spectrum that doesn't resemble real spectra at all.
            var msFile = new TestDataFile(new PeptideWithSetModifications("PEPTIDE", null), 2, 8000, 52);
            CommonParameters commonParameters = new CommonParameters(doPrecursorDeconvolution:false, digestionParams: new DigestionParams(fragmentationTerminus:FragmentationTerminus.C));
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(msFile, null, commonParameters).ToArray();
            PeptideSpectralMatch[] globalPsms = new PeptideSpectralMatch[listOfSortedms2Scans.Length*commonParameters.NumberOfSequencesPerPrecursor];


            DeNovoSequencingEngine deNovoEngine = new DeNovoSequencingEngine(globalPsms, listOfSortedms2Scans, null, null, commonParameters, new List<string>());
            var results = deNovoEngine.Run();

            int psms = globalPsms.Count(x => x != null);
            Assert.AreEqual(psms, commonParameters.NumberOfSequencesPerPrecursor);
            Assert.AreEqual(globalPsms[0].BaseSequence, "PEPTIDE");
            Assert.IsTrue(globalPsms[0].Score > globalPsms[4].Score);
        }
    }
}
