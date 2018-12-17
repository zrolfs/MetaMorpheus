using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.DeNovoSequencing
{
    public class DeNovoTree
    {
        public List<List<string>> Nodes { get; }
        public int CurrentMass { get; }
        public int NumberOfPossibleSequences { get; }
        public double Score { get; private set; }
        public double IntensitySum { get; private set; }

        public DeNovoTree(List<string> combinations, int currentMass, double intensity)
        {
            Nodes = new List<List<string>> { combinations };
            CurrentMass = currentMass;
            NumberOfPossibleSequences = combinations.Count;
            IntensitySum = intensity;
            //NumberOfPossibleSequences = 1;
        }

        public DeNovoTree(List<List<string>> nodes, List<string> addedCombinations, int currentMass, int numberOfPossibleSequences, double intensity)
        {
            Nodes = new List<List<string>>();
            nodes.ForEach(x => Nodes.Add(x));
            Nodes.Add(addedCombinations);
            CurrentMass = currentMass;
            NumberOfPossibleSequences = numberOfPossibleSequences;
            IntensitySum = intensity;
            //NumberOfPossibleSequences = 1;
        }

        public void SetDeNovoScore(long totalProbabilities, double TIC, int maxNodes)
        {
            //Start with 100% probability
            //Divide by the number of trees (three trees each have 33% probability)
            //Divide the remaining probability (33%) within each tree based on the amount of ambiguity in said tree

            //we have the number of possible sequences for each tree to use as a metric for the quality of the assignment
            //however, there are sneaky amino acid combinations that indicate false confidence
            //Example: "PPP" is the only mass for 309.158, but there's not a lot of confidence in "PPP/PPP/PPP" (two peaks)
            //Only one sequence is possible for this tree (PPPPPPPPP), but we only have two peaks!
            //Especially if there's another awesome tree "P/E/P/T/I/D/E", but there are two sequences possible because of I and L (2>1)
            //It doesn't seem fair that "PPPPPPPPP" should outscore "PEPTIDE" in this scenario

            //We're going to divide the number of possible sequences by the number of nodes SQUARED to account for the added confidence that comes from having multiple nodes.

            //Score = (100d * Nodes.Count * Nodes.Count) / (numCompetingTrees * NumberOfPossibleSequences);
            Score = ((100000000d * Math.Pow(1d*Nodes.Count/maxNodes,2)) / totalProbabilities)*IntensitySum/TIC;
        }
    }
}
