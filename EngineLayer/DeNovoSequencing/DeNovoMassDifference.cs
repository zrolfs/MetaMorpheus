using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.DeNovoSequencing
{
    public class DeNovoMassDifference
    {
        public Residue Residue { get; }
        public double MassDifference { get; }
        public double LowestMass { get; }

        public DeNovoMassDifference(Residue residue, double massDifference, double lowestMass)
        {
            Residue = residue;
            MassDifference = massDifference;
            LowestMass = lowestMass;
        }
    }
}
