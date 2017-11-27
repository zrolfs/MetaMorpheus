using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public class SequencesToActualProteinPeptidesEngineResults : MetaMorpheusEngineResults
    {
        #region Public Constructors

        public SequencesToActualProteinPeptidesEngineResults(MetaMorpheusEngine s, Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, Dictionary<string, bool> globalIsDecoy) : base(s)
        {
            this.CompactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
        }
        public SequencesToActualProteinPeptidesEngineResults(MetaMorpheusEngine s, Dictionary<CompactPeptideBase, HashSet<string>> compactPeptideToProteinPeptideMatching, Dictionary<string, bool> globalIsDecoy) : base(s)
        {
            this.CompactPeptideToProteinPeptideMatchingString = compactPeptideToProteinPeptideMatching;
            this.globalIsDecoy = globalIsDecoy;
        }
        #endregion Public Constructors

        #region Public Properties

        public Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> CompactPeptideToProteinPeptideMatching { get; }
        public Dictionary<CompactPeptideBase, HashSet<string>> CompactPeptideToProteinPeptideMatchingString { get; }
        public Dictionary<string, bool> globalIsDecoy { get; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.Append("CompactPeptideToProteinPeptideMatching.Count: " + CompactPeptideToProteinPeptideMatchingString.Count);
            return sb.ToString();
        }

        #endregion Public Methods
    }
}