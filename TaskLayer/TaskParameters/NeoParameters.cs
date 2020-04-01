using System.Collections.Generic;

namespace TaskLayer
{
    public class NeoParameters
    {
        #region Public Constructors

        public NeoParameters()
        {
            Calibrate = true;
            GPTMD = true;
            TargetSearch = true;
            DecoySearch = true;
            SearchNTerminus = true;
            SearchCTerminus = true;
            MaxMissedConsecutiveFragments = 2;
            //MaxMissedTotalFragments = 5;
            MaxCandidatesPerSpectrum = 2000;
            MaxDistanceAllowed = 25;
            NormalCis = true;
            ReverseCis = true;
            TargetFilePath = new List<string>();
            DecoyFilePath = new List<string>();
            NFilePath = new List<string>();
            CFilePath = new List<string>();
        }

        #endregion Public Constructors

        #region Public Properties

        public bool Calibrate { get; set; }
        public double? PrecursorTolerancePPM { get; set; }
        public double? ProductTolerancePPM { get; set; }
        public bool GPTMD { get; set; }
        public bool TargetSearch { get; set; }
        public List<string> TargetFilePath { get; set; }
        public bool DecoySearch { get; set; }
        public List<string> DecoyFilePath { get; set; }

        public bool SearchNTerminus { get; set; }
        public List<string> NFilePath { get; set; }
        public bool SearchCTerminus { get; set; }
        public List<string> CFilePath { get; set; }

        public int MaxMissedConsecutiveFragments { get; set; }
        //public int MaxMissedTotalFragments { get; set; }
        public int MaxCandidatesPerSpectrum { get; set; }

        public int MinDistanceAllowed { get; set; }
        public int MaxDistanceAllowed { get; set; }
        public bool NormalCis { get; set; }
        public bool ReverseCis { get; set; }

        #endregion Public Properties
    }
}