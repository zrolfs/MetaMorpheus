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
        private readonly int NumberOfMS1SpectraToAverage;

        private readonly MsDataFile OriginalFile;
        public MsDataFile AggregatedDataFile { get; private set; }
        public Tolerance SuggestedPrecursorTolerance { get; private set; }
        public Tolerance SuggestedProductTolerance { get; private set; }


        public AggregationEngine(MsDataFile originalFile, CommonParameters commonParameters, List<string> nestedIds, double maxRetentionTimeDifferenceAllowedInMinutes, double minCosineScoreAllowed, int numberOfMS1SpectraToAverage) : base(commonParameters, nestedIds)
        {
            OriginalFile = originalFile;
            MaxRetentionTimeDifferenceAllowedInMinutes = maxRetentionTimeDifferenceAllowedInMinutes;
            MinCosineScoreAllowed = minCosineScoreAllowed;
            NumberOfMS1SpectraToAverage = numberOfMS1SpectraToAverage;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Identifying MS2 groups");
            AggregatedDataFile = new MsDataFile(OriginalFile.GetAllScansList().ToArray(), OriginalFile.SourceFile); 
            return new MetaMorpheusEngineResults(this);
        }

    }
}