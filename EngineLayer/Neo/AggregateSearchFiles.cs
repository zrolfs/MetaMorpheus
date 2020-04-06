using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer.Neo
{
    public static class AggregateSearchFiles
    {
        #region Public Methods

        public static void Combine(string primaryFilePath, string secondaryFilePath, string outputFilePath)
        {
            string[] primaryLines = File.Exists(@primaryFilePath) ? (File.ReadAllLines(@primaryFilePath)): new string[] { "" };
            string[] secondaryLines = File.Exists(@secondaryFilePath) ? (File.ReadAllLines(@secondaryFilePath)): new string[] { "" };
            List<PsmTsvLine> aggregatedLines = AggregateDifferentDatabaseSearches(primaryLines, secondaryLines);
            aggregatedLines = CalculateFDR(aggregatedLines);

            using (StreamWriter file = new StreamWriter(@outputFilePath + "_TargetsAndDecoys.psmtsv"))
            {
                file.WriteLine(primaryLines[0]);
                foreach (PsmTsvLine psm in aggregatedLines)
                    file.WriteLine(psm.ToString());
            }

            List<PsmTsvLine> targets = AssignFDRToTarget(primaryLines, secondaryLines);
            using (StreamWriter file = new StreamWriter(@outputFilePath + "_Targets.psmtsv"))
            {
                file.WriteLine(primaryLines[0]);
                foreach (PsmTsvLine psm in targets)
                    file.WriteLine(psm.ToString());
            }
        }

        public static void Combine(string primaryFilePath, string secondaryFilePath, string outputFilePath, string translatedFilePath)
        {
            string[] primaryLines = File.Exists(@primaryFilePath) ? (File.ReadAllLines(@primaryFilePath)) : new string[] { "" };
            string[] secondaryLines = File.Exists(@secondaryFilePath) ? (File.ReadAllLines(@secondaryFilePath)) : new string[] { "" };
            string[] translatedLines = File.Exists(@translatedFilePath) ? File.ReadAllLines(@translatedFilePath): new string[] { ""};
            List<PsmTsvLine> targetLines = AggregateDifferentDatabaseSearches(primaryLines, translatedLines);
            List<string> targetList = targetLines.Select(x => x.ToString()).ToList();
            targetList.Insert(0, primaryLines[0]); //header
            primaryLines = targetList.ToArray();
            List<PsmTsvLine> aggregatedLines = AggregateDifferentDatabaseSearches(primaryLines, secondaryLines);
            aggregatedLines = CalculateFDR(aggregatedLines);

            using (StreamWriter file = new StreamWriter(@outputFilePath + "_TargetsAndDecoys.psmtsv"))
            {
                file.WriteLine(primaryLines[0]);
                foreach (PsmTsvLine psm in aggregatedLines)
                    file.WriteLine(psm.ToString());
            }

            List<PsmTsvLine> targets = AssignFDRToTarget(primaryLines, secondaryLines);
            using (StreamWriter file = new StreamWriter(@outputFilePath + "_Targets.psmtsv"))
            {
                file.WriteLine(primaryLines[0]);
                foreach (PsmTsvLine psm in targets)
                    file.WriteLine(psm.ToString());
            }
        }

        public static void RecursiveNeoAggregation(string standardFilePath, string neoResultFilePath, string outputFolder, string identifier)
        {
            //This method determines the optimum cutoff for gold standard identification and the minimum score difference required for a splice to outscore a normal
            double desiredConfidence = 0.01;
            double bestQ = -1; //q Cutoff at highest score
            double bestThreshold = -1; //highest number of splice assignments at a 1% local FDR
            int numSplicedHighScore = 0; //highest number of spliced peptides found
            const double thresholdStep = 0.5;
            const int minimumScoreDifference = 1;
            const int maxmimumScoreDifference = 5;

            string[] primaryLines = (File.ReadAllLines(@standardFilePath)); //this one is written by a Neo-Task and should never be missing
            string[] secondaryLines = File.Exists(@neoResultFilePath) ? (File.ReadAllLines(@neoResultFilePath)) : new string[] { "" };
            List<PsmTsvLine> primaryPsms = ImportPsmtsv.ImportLinesToAggregate(primaryLines);
            primaryPsms.ForEach(x => x.neoType = PsmTsvLine.NeoType.Normal);
            List<PsmTsvLine> secondaryPsms = ImportPsmtsv.ImportLinesToAggregate(secondaryLines);

            //fringe case where the qThreshold can't get up(down) to the desired confidence because not enough decoys were found.
            //if qThrehold isn't changing, axe it.
            double previousQThrehold = 0;
            for (double qThreshold = 0; qThreshold <= desiredConfidence;)
            {
                qThreshold = UpdateQThreshold(primaryPsms, qThreshold, true);
                for (double scoreDifferenceThreshold = minimumScoreDifference; scoreDifferenceThreshold < maxmimumScoreDifference; scoreDifferenceThreshold += thresholdStep)
                {
                    List<PsmTsvLine> aggregatedLines = Percolate(primaryPsms, secondaryPsms, qThreshold, scoreDifferenceThreshold);
                    int numSplicedScore = CalculateNumberOfConfidentSpliced(aggregatedLines, desiredConfidence);
                    if (numSplicedHighScore < numSplicedScore)
                    {
                        numSplicedHighScore = numSplicedScore;
                        bestQ = qThreshold;
                        bestThreshold = scoreDifferenceThreshold;
                    }
                }
                if (previousQThrehold == qThreshold) //if we haven't changed (because we can't), get out of the loop
                {
                    break;
                }
                else //update the previous
                {
                    previousQThrehold = qThreshold;
                }
            }
            List<PsmTsvLine> finalAggregatedLines = Percolate(primaryPsms, secondaryPsms, bestQ, bestThreshold);
            finalAggregatedLines = finalAggregatedLines.OrderByDescending(x => x.score).ToList();
            using (StreamWriter file = new StreamWriter(Path.Combine(outputFolder, identifier)))
            {
                file.WriteLine(primaryLines[0] + " \t FusionType" + " \t IsFusion" + " \t IsDecoy"); //header
                foreach (PsmTsvLine line in finalAggregatedLines)
                    file.WriteLine(line.ToString());
            }
            using (StreamWriter file = new StreamWriter(Path.Combine(outputFolder, "PercolatorInfo_" + identifier)))
            {
                file.WriteLine("Maxmimum q-Value of Gold Standards: " + bestQ);
                file.WriteLine("Minimum Score Difference for Splice Selection Over Normal: " + bestThreshold);
                file.WriteLine("FDR for spliced peptides: " + desiredConfidence);
                file.WriteLine("Number of spliced peptides identified at " + desiredConfidence + ": " + numSplicedHighScore);
            }
        }

        public static double UpdateQThreshold(List<PsmTsvLine> primaryLines, double qThreshold, bool increaseQ)
        {
            List<double> qValues = primaryLines.Select(x => Convert.ToDouble(x.q)).ToList(); //grab all q values
            if (increaseQ) //get next highest qValue
            {
                qValues = qValues.OrderBy(x => x).ToList();
                foreach (double qValue in qValues)
                {
                    if (qValue > qThreshold)
                    {
                        return qValue;
                    }
                }
            }
            else //get next lowest qValue
            {
                qValues = qValues.OrderByDescending(x => x).ToList();
                foreach (double qValue in qValues)
                {
                    if (qValue < qThreshold)
                    {
                        return qValue;
                    }
                }
            }
            return qThreshold; //if nothing, return the same qValue and catch it later
        }

        public static int CalculateNumberOfConfidentSpliced(List<PsmTsvLine> aggregatedLines, double desiredConfidence)
        {
            int numConfidentSpliced = 0;
            int numDecoySpliced = 0;
            aggregatedLines = aggregatedLines.OrderByDescending(x => x.score).ToList();
            foreach (PsmTsvLine line in aggregatedLines)
            {
                if (numConfidentSpliced != 0 && (1.0 * numDecoySpliced) / numConfidentSpliced >= desiredConfidence)
                    break;
                if (line.neoType.Equals(PsmTsvLine.NeoType.Spliced))
                    numConfidentSpliced++;
                else if (line.neoType.Equals(PsmTsvLine.NeoType.DecoySpliced))
                    numDecoySpliced++;
            }
            return numConfidentSpliced;
        }

        public static List<PsmTsvLine> Percolate(List<PsmTsvLine> primaryPsms, List<PsmTsvLine> secondaryPsms, double qThreshold, double scoreDifferenceThreshold)
        {
            //get minimum score of qThreshold
            List<PsmTsvLine> primaryAtQ = primaryPsms.Where(x => Convert.ToDouble(x.q) + 0.000001 > qThreshold && Convert.ToDouble(x.q) - 0.000001 < qThreshold).ToList();
            double minScoreAllowed = primaryAtQ.Count == 0 ? 0:primaryAtQ.Min(x => x.score);

            List<PsmTsvLine> aggregatedLines = new List<PsmTsvLine>();
            int p = 0;
            int s = 0;
            //While loop generates a combined list of primary (normal) and secondary (spliced) psms based on scoreDifferenceThreshold
            while (p < primaryPsms.Count && s < secondaryPsms.Count)
            {
                PsmTsvLine psmP = primaryPsms[p];
                PsmTsvLine psmS = secondaryPsms[s];
                if (psmP.scanNumber < psmS.scanNumber)
                {
                    aggregatedLines.Add(psmP);
                    p++;
                }
                else if (psmP.scanNumber > psmS.scanNumber)
                {
                    if (psmS.score > minScoreAllowed)
                    {
                        psmS.neoType = PsmTsvLine.NeoType.Spliced;
                        aggregatedLines.Add(psmS);
                    }
                    s++;
                }
                else
                {
                    if (psmP.score > psmS.score - scoreDifferenceThreshold || psmS.score < minScoreAllowed)
                    {
                        aggregatedLines.Add(psmP);
                    }
                    else
                    {
                        psmS.neoType = (Convert.ToDouble(psmP.q) <= qThreshold) ? PsmTsvLine.NeoType.DecoySpliced : PsmTsvLine.NeoType.Spliced; //if the beaten score belonged to a gold standard, this is a false discovery
                        aggregatedLines.Add(psmS);
                    }
                    p++;
                    s++;
                }
            }
            //wrap up any leftover psms without scanCounts
            for (; p < primaryPsms.Count; p++)
            {
                aggregatedLines.Add(primaryPsms[p]);
            }

            for (; s < secondaryPsms.Count; s++)
            {
                PsmTsvLine psmS = secondaryPsms[s];
                if (psmS.score > minScoreAllowed)
                {
                    psmS.neoType = PsmTsvLine.NeoType.Spliced;
                    aggregatedLines.Add(psmS);
                }
            }

            return aggregatedLines;
        }

        public static int CalculateNumberOfConfidentSpliced(List<PsmTsvLine> aggregatedLines)
        {
            int numConfidentSpliced = 0;
            int numDecoySpliced = 0;
            aggregatedLines = aggregatedLines.OrderByDescending(x => x.score).ToList();
            foreach (PsmTsvLine line in aggregatedLines)
            {
                if (numConfidentSpliced != 0 && (1.0 * numDecoySpliced) / numConfidentSpliced >= 0.05)
                    break;
                if (line.neoType.Equals(PsmTsvLine.NeoType.Spliced))
                    numConfidentSpliced++;
                else if (line.neoType.Equals(PsmTsvLine.NeoType.DecoySpliced))
                    numDecoySpliced++;
            }
            return numConfidentSpliced;
        }

        public static List<PsmTsvLine> AggregateDifferentDatabaseSearches(string[] primaryLines, string[] secondaryLines)
        {
            List<PsmTsvLine> primaryPsms = ImportPsmtsv.ImportLinesToAggregate(primaryLines);
            List<PsmTsvLine> secondaryPsms = ImportPsmtsv.ImportLinesToAggregate(secondaryLines);
            List<PsmTsvLine> aggregatedLines = new List<PsmTsvLine>();
            int p = 0;
            int s = 0;
            while (p < primaryPsms.Count && s < secondaryPsms.Count)
            {
                PsmTsvLine psmP = primaryPsms[p];
                PsmTsvLine psmS = secondaryPsms[s];
                if (psmP.scanNumber < psmS.scanNumber)
                {
                    aggregatedLines.Add(psmP);
                    p++;
                }
                else if (psmP.scanNumber > psmS.scanNumber)
                {
                    aggregatedLines.Add(psmS);
                    s++;
                }
                else
                {
                    if (psmP.score > psmS.score - 0.001 && psmP.score < psmS.score + 0.001)
                        aggregatedLines.Add(psmP.AggregateLine(psmS));
                    else if (psmP.score > psmS.score)
                        aggregatedLines.Add(psmP);
                    else
                        aggregatedLines.Add(psmS);
                    p++;
                    s++;
                }
            }
            for (; p < primaryPsms.Count; p++)
                aggregatedLines.Add(primaryPsms[p]);

            for (; s < secondaryPsms.Count; s++)
                aggregatedLines.Add(secondaryPsms[s]);

            return aggregatedLines;
        }

        #endregion Public Methods

        #region Private Methods

        private static List<PsmTsvLine> CalculateFDR(List<PsmTsvLine> aggregatedLines)
        {
            aggregatedLines = aggregatedLines.OrderByDescending(x => x.score).ToList();
            int targets = 0;
            int decoys = 0;
            foreach (PsmTsvLine line in aggregatedLines)
            {
                if (line.DCT.Contains("T") || line.DCT.Contains("C"))
                    targets++;
                else
                    decoys++;

                line.target = targets.ToString();
                line.decoy = decoys.ToString();
                line.q = ((1.0d * decoys) / (targets)).ToString();
            }
            return aggregatedLines;
        }

        private static List<PsmTsvLine> AssignFDRToTarget(string[] primaryLines, string[] secondaryLines)
        {
            List<PsmTsvLine> primaryPsms = ImportPsmtsv.ImportLinesToAggregate(primaryLines);
            List<PsmTsvLine> secondaryPsms = ImportPsmtsv.ImportLinesToAggregate(secondaryLines);
            primaryPsms = primaryPsms.OrderByDescending(x => x.score).ToList();
            secondaryPsms = secondaryPsms.OrderByDescending(x => x.score).ToList();
            int p = 0;
            int s = 0;
            int target = 0;
            int decoy = 0;
            double qMax = 0;
            double decoyScore = secondaryPsms.Count != 0 ? secondaryPsms[s].score : 0;

            while (p < primaryPsms.Count)
            {
                PsmTsvLine targetLine = primaryPsms[p];
                if (targetLine.score > decoyScore || s == secondaryPsms.Count - 1)
                {
                    target++;
                    targetLine.target = target.ToString();
                    targetLine.decoy = decoy.ToString();
                    double qValue = (1.0d * decoy / target);
                    qMax = (qMax > qValue) ? qMax : qValue;
                    targetLine.q = qMax.ToString();
                    p++;
                }
                else
                {
                    decoy++;
                    s++;
                    decoyScore = secondaryPsms[s].score;
                }
            }
            return primaryPsms;
        }

        private static void SubtractScoresFromFusions(List<PsmTsvLine> finalList, double scoreCutOff)
        {
            foreach (PsmTsvLine line in finalList)
                if (line.neoType != PsmTsvLine.NeoType.Normal)
                    line.score -= scoreCutOff;
        }

        #endregion Private Methods
    }
}