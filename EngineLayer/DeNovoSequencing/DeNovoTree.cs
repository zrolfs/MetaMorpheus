using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.DeNovoSequencing
{
    public class DeNovoTree
    {
        public List<List<string>> Nodes { get; }
        public int CurrentMass { get; }
        //public int ProductMassShift { get; }
        public int NumberOfPossibleSequences { get; private set; }

        public DeNovoTree(List<string> combinations, int currentMass)//, int productMassShift)
        {
            //Nodes = new List<DeNovoNode> { new DeNovoNode(combinations, currentMass) };
            Nodes = new List<List<string>> { combinations };
            CurrentMass = currentMass;
           // ProductMassShift = productMassShift;
        }

        public DeNovoTree(List<List<string>> nodes, List<string> addedCombinations, int massDifference, int currentMass)
        {
            Nodes = new List<List<string>>();
            nodes.ForEach(x => Nodes.Add(x));
            //Nodes.Add(new DeNovoNode(addedCombinations, massDifference));
            Nodes.Add(addedCombinations);
            CurrentMass = currentMass;
        }

        public void CalculateNumberOfPossibleSequences()
        {
            NumberOfPossibleSequences = Nodes.Count != 0 ? 1 : 0;
            foreach(List<string> node in Nodes)
            {
                NumberOfPossibleSequences *= node.Count;
            }
        }
    }
}
