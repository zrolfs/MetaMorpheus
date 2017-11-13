using Proteomics;

namespace EngineLayer
{
    public abstract class Peptide
    {
        #region Private Fields

        private string baseSequence;

        #endregion Private Fields

        #region Protected Constructors

        protected Peptide(Protein protein, int oneBasedStartResidueInProtein, int oneBasedEndResidueInProtein, string peptideDescription = null)
        {
            Protein = protein;
            OneBasedStartResidueInProtein = oneBasedStartResidueInProtein;
            OneBasedEndResidueInProtein = oneBasedEndResidueInProtein;
            Length = OneBasedEndResidueInProtein - OneBasedStartResidueInProtein + 1;
            PeptideDescription = peptideDescription;
            cis = false;
        }
        protected Peptide(int startTwo, int endTwo,Protein protein, int oneBasedStartResidueInProtein, int oneBasedEndResidueInProtein, string peptideDescription = null)
        {
            Protein = protein;
            OneBasedStartResidueInProtein = oneBasedStartResidueInProtein;
            OneBasedEndResidueInProtein = oneBasedEndResidueInProtein;
            Length = OneBasedEndResidueInProtein - OneBasedStartResidueInProtein + 1;
            PeptideDescription = peptideDescription;
            this.startTwo = startTwo;
            this.endTwo = endTwo;
            cis = true;
        }


        #endregion Protected Constructors

        #region Public Properties

        public Protein Protein { get; }
        public int OneBasedStartResidueInProtein { get; }
        public int OneBasedEndResidueInProtein { get; }
        public bool cis { get; }
        public string PeptideDescription { get; }
        public int startTwo { get; }
        public int endTwo { get; }

        public int Length { get; }

        public virtual char PreviousAminoAcid
        {
            get
            {
                return OneBasedStartResidueInProtein > 1 ? Protein[OneBasedStartResidueInProtein - 2] : '-';
            }
        }

        public virtual char NextAminoAcid
        {
            get
            {
                return endTwo < Protein.Length ? Protein[endTwo] : '-';
            }
        }

        public string BaseSequence
        {
            get
            {
                if (baseSequence == null)
                    baseSequence = Protein.BaseSequence.Substring(OneBasedStartResidueInProtein - 1, Length);
                return baseSequence;
            }
        }

        #endregion Public Properties

        #region Public Indexers

        public char this[int zeroBasedIndex]
        {
            get
            {
                if (cis && zeroBasedIndex + OneBasedStartResidueInProtein - 1 > OneBasedEndResidueInProtein)
                    return Protein.BaseSequence[zeroBasedIndex + OneBasedStartResidueInProtein - OneBasedEndResidueInProtein + startTwo - 1];
                else
                    return Protein.BaseSequence[zeroBasedIndex + OneBasedStartResidueInProtein - 1];
            }
        }

        #endregion Public Indexers
    }
}