namespace TaskLayer
{
    public class AggregationParameters
    {
        public AggregationParameters()
        {
            MaxRetentionTimeDifferenceAllowedInMinutes = 10;
            MinCosineScoreAllowed = 0.075;
            NumberOfMS1SpectraToAverage = int.MaxValue;
        }

        public double MaxRetentionTimeDifferenceAllowedInMinutes { get; set; }
        public double MinCosineScoreAllowed { get; set; }
        public int NumberOfMS1SpectraToAverage { get; set; }
    }
}