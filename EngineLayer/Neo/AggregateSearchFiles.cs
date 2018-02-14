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
            string[] primaryLines = (System.IO.File.ReadAllLines(@primaryFilePath));
            string[] secondaryLines = (System.IO.File.ReadAllLines(@secondaryFilePath));
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
            string[] primaryLines = (System.IO.File.ReadAllLines(@primaryFilePath));
            string[] secondaryLines = (System.IO.File.ReadAllLines(@secondaryFilePath));
            string[] translatedLines = (System.IO.File.ReadAllLines(@translatedFilePath));
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
            double qThreshold = 0;
            double oldQThreshold = 0;
            double scoreDifferenceThreshold = 1;
            double oldScoreDifferenceThreshold = 1;
            int numSplicedHighScoreQ = -1; //highest number of splice assignments at a 1% local FDR
            int numSplicedHighScoreThreshold = -1; //highest number of splice assignments at a 1% local FDR
            int numSplicedScore = 0; //current number of splic assignments at a 1% local FDR
            bool increaseQ = true;
            bool increaseScoreDifference = true;

            string[] primaryLines = (System.IO.File.ReadAllLines(@standardFilePath));
            string[] secondaryLines = (System.IO.File.ReadAllLines(@neoResultFilePath));
            List<PsmTsvLine> primaryPsms = ImportPsmtsv.ImportLinesToAggregate(primaryLines);
            primaryPsms.ForEach(x => x.neoType = PsmTsvLine.NeoType.Normal);
            List<PsmTsvLine> secondaryPsms = ImportPsmtsv.ImportLinesToAggregate(secondaryLines);

            do //determine if score difference should be changed
            {
                if (numSplicedScore < numSplicedHighScoreThreshold || (numSplicedScore == numSplicedHighScoreThreshold && increaseScoreDifference && UpdateScoreDifferenceThreshold(scoreDifferenceThreshold, increaseScoreDifference) > 20)) //check second time around if move score Threhold the other way
                    increaseScoreDifference = false;
                else
                {
                    numSplicedHighScoreThreshold = numSplicedScore; //update highscore
                    oldScoreDifferenceThreshold = scoreDifferenceThreshold; //update score difference
                }
                scoreDifferenceThreshold = UpdateScoreDifferenceThreshold(scoreDifferenceThreshold, increaseScoreDifference); //update score Threshold
                double qStart = qThreshold;
                List<PsmTsvLine> aggregatedLines = Percolate(primaryPsms, secondaryPsms, qThreshold, scoreDifferenceThreshold);
                numSplicedScore = CalculateNumberOfConfidentSpliced(aggregatedLines);
                do //determine gold standards to use
                {
                    oldQThreshold = qThreshold; //updateQ

                    if (numSplicedScore < numSplicedHighScoreQ || (numSplicedScore == numSplicedHighScoreQ && increaseQ && UpdateQThreshold(primaryPsms, qThreshold, increaseQ) > 0.05)) //check second time around if move qValue the other way
                    {
                        if (numSplicedHighScoreThreshold == numSplicedScore && qStart==0)
                            break;
                        increaseQ = false;
                    }
                    else
                    {
                        numSplicedHighScoreQ = numSplicedScore; //update highscore
                    }

                    qThreshold = UpdateQThreshold(primaryPsms, qThreshold, increaseQ); //get qValue
                    aggregatedLines = Percolate(primaryPsms, secondaryPsms, qThreshold, scoreDifferenceThreshold);
                    numSplicedScore = CalculateNumberOfConfidentSpliced(aggregatedLines);
                } while ((numSplicedScore >= numSplicedHighScoreQ || increaseQ) 
                    && (numSplicedScore > numSplicedHighScoreQ || qThreshold < 0.05 || !increaseQ) && qThreshold > 0); //do again the otherway if done increasing
                if (qThreshold == 0 && numSplicedScore>=numSplicedHighScoreQ)
                    oldQThreshold = qThreshold;
                List<PsmTsvLine> oldAggregatedLines = Percolate(primaryPsms, secondaryPsms, oldQThreshold, scoreDifferenceThreshold);
                numSplicedScore = CalculateNumberOfConfidentSpliced(oldAggregatedLines);
                increaseQ = true;
                qThreshold = oldQThreshold;
            } while ((numSplicedScore >= numSplicedHighScoreThreshold || increaseScoreDifference) 
                && (numSplicedScore > numSplicedHighScoreThreshold || ((scoreDifferenceThreshold < 20 || !increaseScoreDifference) && scoreDifferenceThreshold > -20)));
            List<PsmTsvLine> finalAggregatedLines = Percolate(primaryPsms, secondaryPsms, oldQThreshold, oldScoreDifferenceThreshold);
            numSplicedScore = CalculateNumberOfConfidentSpliced(finalAggregatedLines);
            SubtractScoresFromFusions(finalAggregatedLines, oldScoreDifferenceThreshold);
            finalAggregatedLines = finalAggregatedLines.OrderByDescending(x => x.score).ToList();
            using (StreamWriter file = new StreamWriter(Path.Combine(outputFolder, identifier)))
            {
                file.WriteLine(primaryLines[0]+" \t FusionType"); //header
                foreach (PsmTsvLine line in finalAggregatedLines)
                    file.WriteLine(line.ToString());
            }
            using (StreamWriter file = new StreamWriter(Path.Combine(outputFolder, "PercolatorInfo_" + identifier)))
            {
                file.WriteLine("Maxmimum q-Value of Gold Standards: " + oldQThreshold);
                file.WriteLine("Minimum Score Difference for Splice Selection Over Normal: " + oldScoreDifferenceThreshold);
            }
        }

        public static double UpdateQThreshold(List<PsmTsvLine> primaryLines, double qThreshold, bool increaseQ)
        {
            List<double> qValues = primaryLines.Select(x => Convert.ToDouble(x.q)).ToList(); //grab all q values
            if (increaseQ) //get next highest qValue
            {
                qValues = qValues.OrderBy(x => x).ToList();
                foreach (double qValue in qValues)
                    if (qValue > qThreshold)
                        return qValue;
            }
            else //get next lowest qValue
            {
                qValues = qValues.OrderByDescending(x => x).ToList();
                foreach (double qValue in qValues)
                    if (qValue < qThreshold)
                        return qValue;
            }
            return qThreshold; //if nothing, return the same qValue
        }

        public static double UpdateScoreDifferenceThreshold(double scoreDifferenceThreshold, bool increaseScoreDifference)
        {
            if (increaseScoreDifference)
                return scoreDifferenceThreshold += 0.05;
            else
                return scoreDifferenceThreshold -= 0.05;
        }

        public static List<PsmTsvLine> Percolate(List<PsmTsvLine> primaryPsms, List<PsmTsvLine> secondaryPsms, double qThreshold, double scoreDifferenceThreshold)
        {
            //get minimum score of qThreshold
            List<PsmTsvLine> primaryAtQ = primaryPsms.Where(x => Convert.ToDouble(x.q) + 0.000001 > qThreshold && Convert.ToDouble(x.q) - 0.000001 < qThreshold).ToList();
            double minScoreAllowed = primaryAtQ.Min(x => x.score);

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
                else if (psmP.scanNumber > psmS.scanNumber && psmS.score-scoreDifferenceThreshold > minScoreAllowed)
                {
                    psmS.neoType = PsmTsvLine.NeoType.Spliced;
                    aggregatedLines.Add(psmS);
                    s++;
                }
                else
                {
                    if (psmP.score > psmS.score - scoreDifferenceThreshold || psmS.score - scoreDifferenceThreshold < minScoreAllowed)
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
                aggregatedLines.Add(primaryPsms[p]);

            for (; s < secondaryPsms.Count; s++)
                aggregatedLines.Add(secondaryPsms[s]);

            return aggregatedLines;
        }

        public static int CalculateNumberOfConfidentSpliced(List<PsmTsvLine> aggregatedLines)
        {
            int numConfidentSpliced = 0;
            int numDecoySpliced = 0;
            aggregatedLines = aggregatedLines.OrderByDescending(x => x.score).ToList();
            foreach (PsmTsvLine line in aggregatedLines)
            {
                if (numConfidentSpliced!=0 && (1.0 * numDecoySpliced) / numConfidentSpliced >= 0.05)
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
                aggregatedLines.Add(primaryPsms[s]);

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
            PsmTsvLine decoyLine = secondaryPsms[s];

            while (p < primaryPsms.Count)
            {
                PsmTsvLine targetLine = primaryPsms[p];
                if (targetLine.score > decoyLine.score || s == secondaryPsms.Count)
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
                    decoyLine = secondaryPsms[s];
                }
            }
            return primaryPsms;
        }

        private static void SubtractScoresFromFusions(List<PsmTsvLine> finalList, double scoreCutOff)
        {
            foreach(PsmTsvLine line in finalList)
                if(line.neoType!=PsmTsvLine.NeoType.Normal)
                    line.score -= scoreCutOff;
        }

        #endregion Private Methods
    }
}