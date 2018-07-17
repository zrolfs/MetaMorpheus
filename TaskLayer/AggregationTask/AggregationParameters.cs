namespace TaskLayer
{
    public class AggregationParameters
    {
        public AggregationParameters()
        {
            MaxRetentionTimeDifferenceAllowedInMinutes = 10;
            MinCosineScoreAllowed = 0.075;
        }

        public double MaxRetentionTimeDifferenceAllowedInMinutes { get; set; }
        public double MinCosineScoreAllowed { get; set; }
    }
}