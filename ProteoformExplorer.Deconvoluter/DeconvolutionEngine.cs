using BayesianEstimation;
using Chemistry;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using ThermoFisher.CommonCore.Data.Business;
using UsefulProteomicsDatabases;
using Polarity = MassSpectrometry.Polarity;
using Tolerance = MzLibUtil.Tolerance;

namespace ProteoformExplorer.Deconvoluter
{
    public class DeconvolutionEngine
    {
        private static readonly double[][] allMasses;
        private static readonly double[][] allIntensities;
        private static readonly double[] mostIntenseMasses;
        private static readonly double[] diffToMonoisotopic;

        static DeconvolutionEngine()
        {
            Loaders.LoadElements();

            // AVERAGINE
            const double averageC = 5.0359;
            const double averageH = 7.9273;
            const double averageO = 1.52996;
            const double averageN = 1.3608;
            const double averageS = 0.0342;

            const double fineRes = 0.125;
            const double minRes = 1e-8;

            double maxMass = 100000;
            double averagineDaltonResolution = 50.0;

            double averagineMass = PeriodicTable.GetElement("C").PrincipalIsotope.AtomicMass * averageC
                + PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * averageH
                + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass * averageO
                + PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass * averageN
                + PeriodicTable.GetElement("S").PrincipalIsotope.AtomicMass * averageS;

            int numAveragines = (int)Math.Ceiling(maxMass / averagineDaltonResolution) + 1;

            allMasses = new double[numAveragines][];
            allIntensities = new double[numAveragines][];
            mostIntenseMasses = new double[numAveragines];
            diffToMonoisotopic = new double[numAveragines];

            for (int i = 1; i < numAveragines; i++)
            {
                double mass = i * averagineDaltonResolution;

                double averagines = mass / averagineMass;

                if (mass < 50)
                {
                    continue;
                }

                ChemicalFormula chemicalFormula = new ChemicalFormula();
                chemicalFormula.Add("C", Convert.ToInt32(averageC * averagines));
                chemicalFormula.Add("H", Convert.ToInt32(averageH * averagines));
                chemicalFormula.Add("O", Convert.ToInt32(averageO * averagines));
                chemicalFormula.Add("N", Convert.ToInt32(averageN * averagines));
                chemicalFormula.Add("S", Convert.ToInt32(averageS * averagines));

                var chemicalFormulaReg = chemicalFormula;
                IsotopicDistribution isotopicDistribution = IsotopicDistribution.GetDistribution(chemicalFormulaReg, fineRes, minRes);
                var masses = isotopicDistribution.Masses.ToArray();
                var intensities = isotopicDistribution.Intensities.ToArray();
                var mostIntense = intensities.Max();
                int indOfMostIntensePeak = Array.IndexOf(intensities, mostIntense);

                for (int j = 0; j < intensities.Length; j++)
                {
                    intensities[j] /= mostIntense;
                }

                mostIntenseMasses[i] = masses[indOfMostIntensePeak];
                diffToMonoisotopic[i] = masses[indOfMostIntensePeak] - chemicalFormulaReg.MonoisotopicMass;
                allMasses[i] = masses;
                allIntensities[i] = intensities;
            }
        }

        public double MinMass { get; private set; }
        public double MinFractionIntensityRequired { get; private set; }
        public double IntensityRatioLimit { get; private set; }
        public double PearsonCorrelationRequired { get; private set; }
        public double SignalToNoiseRequired { get; private set; }
        public Tolerance PpmTolerance { get; private set; }
        public int MinCharge { get; private set; }
        public int MaxCharge { get; private set; }
        public int MinPeaks { get; private set; }

        public DeconvolutionEngine(double minMass, double minFractionIntensityRequired, double intensityRatio, double pearsonCorrelationRequired, double signalToNoiseRequired,
            double ppmTolerance, int minCharge, int maxCharge, int minPeaks)
        {
            MinMass = minMass;
            this.MinFractionIntensityRequired = minFractionIntensityRequired;
            this.IntensityRatioLimit = intensityRatio;
            this.PearsonCorrelationRequired = pearsonCorrelationRequired;
            this.SignalToNoiseRequired = signalToNoiseRequired;
            this.PpmTolerance = new PpmTolerance(ppmTolerance);
            this.MinCharge = minCharge;
            this.MaxCharge = maxCharge;
            this.MinPeaks = minPeaks;
        }

        public IEnumerable<DeconvolutedEnvelope> Deconvolute(MsDataFile file, string fileName = null, int threads = -1)
        {
            var items = new List<DeconvolutedEnvelope>();

            var scans = file.GetMS1Scans().ToList();
            int scansComplete = 0;

            Parallel.ForEach(Partitioner.Create(0, scans.Count),
               new ParallelOptions { MaxDegreeOfParallelism = threads },
               (range, loopState) =>
               {
                   for (int i = range.Item1; i < range.Item2; i++)
                   {
                       var scan = scans[i];

                       var deconvolutedEnvs = Deconvolute(scan).ToList();

                       foreach (var env in deconvolutedEnvs)
                       {
                           env.SpectraFileName = fileName;
                       }

                       lock (items)
                       {
                           items.AddRange(deconvolutedEnvs);
                           scansComplete++;
                       }
                   }
               });

            foreach (var item in items.OrderByDescending(p => p.Score))
            {
                // TODO: feature finding
                yield return item;
            }
        }

        public IEnumerable<DeconvolutedEnvelope> Deconvolute(MsDataScan scan)
        {
            var deconvolutedEnvs = Deconvolute(scan.MassSpectrum, scan.MassSpectrum.Range, scan.Polarity).ToList();

            foreach (var env in deconvolutedEnvs)
            {
                env.RetentionTime = scan.RetentionTime;
                env.OneBasedScan = scan.OneBasedScanNumber;
                env.TotalScanDeconvolutedIntensity = deconvolutedEnvs.Sum(p => p.Peaks.Sum(v => v.ExperimentalIntensity));

                yield return env;
            }
        }

        public IEnumerable<DeconvolutedEnvelope> Deconvolute(MzSpectrum spectrum, MzRange mzRange, Polarity polarity)
        {
            // if no peaks in the scan, stop
            if (spectrum.Size == 0)
            {
                yield break;
            }

            // get list of envelope candidates for this scan
            var indicies = GetPeaksThatPassSignalToNoiseFilter(spectrum).ToList();
            var candidateEnvelopes = GetEnvelopeCandidates(spectrum, mzRange, indicies, polarity);
            var parsimoniousEnvelopes = RunEnvelopeParsimony(candidateEnvelopes, spectrum);

            // return deconvoluted envelopes
            foreach (DeconvolutedEnvelope envelope in parsimoniousEnvelopes.Where(p =>
                p.SignalToNoise >= SignalToNoiseRequired
                && p.MonoisotopicMass >= MinMass
                && (polarity == Polarity.Positive && p.Charge >= MinCharge || polarity == Polarity.Negative && p.Charge <= MinCharge)
                && p.Peaks.Count >= MinPeaks))
            {
                yield return envelope;
            }
        }

        public List<DeconvolutedEnvelope> GetEnvelopeCandidates(MzSpectrum spectrum, MzRange mzRange, List<int> optionalIndicies = null, Polarity polarity = Polarity.Positive)
        {
            List<DeconvolutedEnvelope> envelopeCandidates = new List<DeconvolutedEnvelope>();
            List<DeconvolutedPeak> peaksBuffer = new List<DeconvolutedPeak>();
            HashSet<double> mzsClaimed = new HashSet<double>(); // this is empty, no m/z peaks have been claimed yet
            HashSet<int> potentialChargeStates = new HashSet<int>();
            List<(double, double)> intensitiesBuffer = new List<(double, double)>();

            // get list of envelope candidates for this scan
            for (int p = 0; p < spectrum.XArray.Length; p++)
            {
                double mz = spectrum.XArray[p];

                if (optionalIndicies != null && !optionalIndicies.Contains(p))
                {
                    continue;
                }

                // check to see if this peak is in the m/z deconvolution range
                if (mz < mzRange.Minimum)
                {
                    continue;
                }
                else if (mz > mzRange.Maximum)
                {
                    break;
                }

                // get rough list of charge states to check for based on m/z peaks around this peak
                potentialChargeStates = GetPotentialChargeStates(potentialChargeStates, spectrum, p, polarity);

                // examine different charge state possibilities and get corresponding envelope candidates
                foreach (int z in potentialChargeStates)
                {
                    DeconvolutedEnvelope candidateEnvelope = GetIsotopicEnvelope(spectrum, p, z, peaksBuffer, mzsClaimed, intensitiesBuffer);

                    if (candidateEnvelope != null)
                    {
                        envelopeCandidates.Add(candidateEnvelope);
                    }
                }
            }

            // determing which envelopes have neighboring charges
            //Dictionary<int, List<DeconvolutedEnvelope>> envelopesGroupedByCharge = envelopeCandidates.GroupBy(p => p.Charge).ToDictionary(p => p.Key, v => v.ToList());

            //foreach (var candidate in envelopeCandidates)
            //{
            //    for (int z = 1; z <= MaxCharge; z++)
            //    {
            //        double mass = candidate.Peaks.First().ExperimentalMz.ToMass(candidate.Charge);

            //        if (!envelopesGroupedByCharge.TryGetValue(z, out var envelopesWithThisCharge))
            //        {
            //            continue;
            //        }

            //        bool chargeAndMassObserved = envelopesWithThisCharge.Any(p => p.Peaks.Any(v => PpmTolerance.Within(v.ExperimentalMz.ToMass(p.Charge), mass)));

            //        if (chargeAndMassObserved)
            //        {
            //            candidate.NeighboringCharges++;
            //        }
            //    }
            //}

            return envelopeCandidates;
        }

        public List<DeconvolutedEnvelope> RunEnvelopeParsimony(List<DeconvolutedEnvelope> envelopeCandidates, MzSpectrum spectrum)
        {
            var parsimoniousEnvelopes = new List<DeconvolutedEnvelope>();
            var peaksBuffer = new List<DeconvolutedPeak>();
            HashSet<double> mzsClaimed = new HashSet<double>();
            List<DeconvolutedEnvelope> overlappingEnvelopes = new List<DeconvolutedEnvelope>();
            List<(double, double)> intensitiesBuffer = new List<(double, double)>();

            // greedy algorithm. pick a set of mutually-exclusive isotopic envelopes, ordered by score,
            // until no more envelopes can be found in the spectrum
            DeconvolutedEnvelope nextEnvelope = GetNextBestEnvelope(envelopeCandidates, null);
            //bool harmonicsRemoved = false;

            while (nextEnvelope != null)
            {
                DeconvolutedEnvelope chosenEnvelope = nextEnvelope;
                parsimoniousEnvelopes.Add(chosenEnvelope);

                foreach (var peak in chosenEnvelope.Peaks)
                {
                    mzsClaimed.Add(peak.ExperimentalMz);
                }

                overlappingEnvelopes.Clear();
                overlappingEnvelopes.AddRange(envelopeCandidates.Where(p => p.Peaks.Any(v => mzsClaimed.Contains(v.ExperimentalMz))));

                //TODO: remove trailing edges on left/right side of envelope
                // iterate through isotopes, if SN < 1 or 2, add to list of mzsClaimed

                foreach (DeconvolutedEnvelope overlappingEnvelope in overlappingEnvelopes.OrderByDescending(p => p.Score))
                {
                    int mzIndex = spectrum.GetClosestPeakIndex(overlappingEnvelope.Peaks.First().ExperimentalMz);
                    int candidateEnvelopeIndex = envelopeCandidates.IndexOf(overlappingEnvelope);

                    DeconvolutedEnvelope redeconEnv = GetIsotopicEnvelope(spectrum, mzIndex, overlappingEnvelope.Charge, peaksBuffer, mzsClaimed,
                        intensitiesBuffer);

                    if (redeconEnv != null)
                    {
                        redeconEnv.NoiseFwhm = overlappingEnvelope.NoiseFwhm;
                        //redeconEnv.SignalToNoise = overlappingEnvelope.SignalToNoise;
                        redeconEnv.Baseline = overlappingEnvelope.Baseline;
                    }

                    envelopeCandidates[candidateEnvelopeIndex] = redeconEnv;
                }

                // remove all invalid isotopic envelopes
                envelopeCandidates.RemoveAll(p => p == null);

                //remove harmonics
                //if (chosenEnvelope.MonoisotopicMass > 10000 && !harmonicsRemoved)
                //{
                //    RemoveHarmonics(envelopeCandidates, chosenEnvelope, spectrum);
                //    harmonicsRemoved = true;
                //}

                // get the next-best envelope
                nextEnvelope = GetNextBestEnvelope(envelopeCandidates, chosenEnvelope);
            }

            return parsimoniousEnvelopes;
        }

        public void RemoveHarmonics(List<DeconvolutedEnvelope> envelopeCandidates, DeconvolutedEnvelope previousBestEnvelope, MzSpectrum spectrum)
        {


            //HashSet<DeconvolutedEnvelope> harmonicEnvelopes = new HashSet<DeconvolutedEnvelope>();
            //HashSet<double> mzs = new HashSet<double>();
            //List<DeconvolutedEnvelope> temp = new List<DeconvolutedEnvelope>();

            //Dictionary<double, List<DeconvolutedEnvelope>> envelopesGroupedByMz = new Dictionary<double, List<DeconvolutedEnvelope>>();

            //foreach (var envelope in envelopeCandidates)
            //{
            //    foreach (var peak in envelope.Peaks)
            //    {
            //        if (!envelopesGroupedByMz.ContainsKey(peak.ExperimentalMz))
            //        {
            //            envelopesGroupedByMz.Add(peak.ExperimentalMz, new List<DeconvolutedEnvelope>());
            //        }

            //        envelopesGroupedByMz[peak.ExperimentalMz].Add(envelope);
            //    }
            //}

            //foreach (DeconvolutedEnvelope envelope in envelopeCandidates)
            //{
            //    mzs.Clear();
            //    temp.Clear();



            //    //int neighboringChargeCount = neighboringCharge.Count();

            //    foreach (var peak in envelope.Peaks)
            //    {
            //        mzs.Add(peak.ExperimentalMz);

            //        temp.AddRange(envelopesGroupedByMz[peak.ExperimentalMz].Where(p =>
            //            envelope.Charge % p.Charge == 0 &&
            //            envelope.Charge != p.Charge));
            //    }

            //    if (!temp.Any())
            //    {
            //        continue;
            //    }

            //    var neighboringCharge = envelopeCandidates.Where(p => (p.Charge == envelope.Charge + 1 || p.Charge == envelope.Charge - 1)
            //        && PpmTolerance.Within(p.MonoisotopicMass, envelope.MonoisotopicMass)).ToList();

            //    if (neighboringCharge.Count == 0)
            //    {
            //        continue;
            //    }

            //    var potentialHarmonicEnvelopes = temp.Where(p =>
            //        p.Peaks.Count(v => mzs.Contains(v.ExperimentalMz)) >= 2);
            //    //.ToList();

            //    foreach (var env in potentialHarmonicEnvelopes)
            //    {
            //        harmonicEnvelopes.Add(env);
            //    }

            //    //for (int i = 2; i < MaxCharge / envelope.Charge; i++)
            //    //{
            //    //    int harmonicZ = envelope.Charge * i;
            //    //    double mainMz = envelope.Peaks.First().ExperimentalMz;

            //    //    var theoreticalHarmonicMz = (mainMz.ToMass(harmonicZ) + Constants.C13MinusC12).ToMz(harmonicZ);
            //    //    int ind = spectrum.GetClosestPeakIndex(theoreticalHarmonicMz);
            //    //    double experimentalHarmonicMz = spectrum.XArray[ind];

            //    //    if (PpmTolerance.Within(experimentalHarmonicMz.ToMass(harmonicZ), theoreticalHarmonicMz.ToMass(harmonicZ))
            //    //        && experimentalHarmonicMz != mainMz)
            //    //    {
            //    //        var harmonicCandidates = envelopeCandidates.Where(p => p.Charge == harmonicZ
            //    //            && p.Peaks.Any(v => v.ExperimentalMz == experimentalHarmonicMz || v.ExperimentalMz == mainMz))
            //    //            .ToList(); // DEBUG

            //    //        if (harmonicCandidates.Any())
            //    //        {
            //    //            harmonicEnvelopes.Add(envelope);
            //    //        }
            //    //    }

            //    //    // TODO: look for -1 dalton peak
            //    //}
            //}

            //foreach (DeconvolutedEnvelope harmonicEnvelope in harmonicEnvelopes)
            //{
            //    envelopeCandidates.Remove(harmonicEnvelope);
            //}
        }

        public DeconvolutedEnvelope GetIsotopicEnvelope(MzSpectrum spectrum, int p, int z, List<DeconvolutedPeak> deconvolutedPeaks,
            HashSet<double> alreadyClaimedMzs, List<(double, double)> intensitiesBuffer)
        {
            double mz = spectrum.XArray[p];
            double intensity = spectrum.YArray[p];

            if (alreadyClaimedMzs.Contains(mz))
            {
                return null;
            }

            double mass = mz.ToMass(z);

            // get the index of an averagine envelope close in mass
            int massIndex = GetMassIndex(mass) + 1;

            if (massIndex >= mostIntenseMasses.Length)
            {
                return null;
            }

            double[] averagineEnvelopeMasses = allMasses[massIndex];
            double[] averagineEnvelopeIntensities = allIntensities[massIndex];
            double monoMass = mass - diffToMonoisotopic[massIndex];

            GetPotentialIsotopePeaks(deconvolutedPeaks, massIndex, alreadyClaimedMzs, spectrum, z, mass, intensity);

            if (deconvolutedPeaks.Count < 2)
            {
                return null;
            }

            // calculate % intensity missing
            double sumIntensity = deconvolutedPeaks.Sum(p => Math.Min(p.ExperimentalIntensity, p.TheoreticalIntensity));
            double expectedTotalIntensity = 0;

            double maxTheorIntensity = deconvolutedPeaks.First().TheoreticalIntensity;
            for (int i = 0; i < averagineEnvelopeIntensities.Length; i++)
            {
                double expectedIsotopeIntensity = averagineEnvelopeIntensities[i] * maxTheorIntensity;
                expectedTotalIntensity += expectedIsotopeIntensity;
            }

            double fracIntensityObserved = sumIntensity / expectedTotalIntensity;

            if (fracIntensityObserved < MinFractionIntensityRequired)
            {
                return null;
            }

            // calculate correlation to averagine
            double corr = Correlation.Pearson(deconvolutedPeaks.Select(p => p.ExperimentalIntensity), deconvolutedPeaks.Select(p => p.TheoreticalIntensity));

            DeconvolutedEnvelope env = null;

            // this is just to save memory, but quality filtering can happen outside of this method after the envelope has been returned, if desired
            if (corr >= PearsonCorrelationRequired)
            {
                // create + return the isotopic envelope object
                env = new DeconvolutedEnvelope(deconvolutedPeaks.ToList(), monoMass, z, corr, fracIntensityObserved);
            }
            else
            {
                // see if a subset of the peaks is a valid envelope
                List<DeconvolutedEnvelope> subsetEnvelopeCandidates = new List<DeconvolutedEnvelope>();

                //TODO: calculate subsets of peaks that have been gathered. return the best one that meets the filtering criteria,
                // or all of the ones that meet the filtering criteria?
                var sortedPeaks = deconvolutedPeaks.OrderBy(p => p.TheoreticalMz).ToList();

                List<DeconvolutedPeak> subsetPeaks = new List<DeconvolutedPeak>();
                for (int start = 0; start < sortedPeaks.Count; start++)
                {
                    for (int end = sortedPeaks.Count - 1; end >= start + 1; end--)
                    {
                        subsetPeaks.Clear();

                        for (int k = start; k <= end; k++)
                        {
                            subsetPeaks.Add(sortedPeaks[k]);
                        }

                        //TODO: make this more efficient by changing start + end
                        if (!subsetPeaks.Any(p => p.ExperimentalMz == mz))
                        {
                            continue;
                        }

                        corr = Correlation.Pearson(subsetPeaks.Select(p => p.ExperimentalIntensity), subsetPeaks.Select(p => p.TheoreticalIntensity));
                        sumIntensity = subsetPeaks.Sum(p => Math.Min(p.ExperimentalIntensity, p.TheoreticalIntensity));
                        fracIntensityObserved = sumIntensity / expectedTotalIntensity;

                        if (corr >= PearsonCorrelationRequired && fracIntensityObserved >= MinFractionIntensityRequired && subsetPeaks.Count >= 2)
                        {
                            var subsetEnvelope = new DeconvolutedEnvelope(subsetPeaks.OrderBy(p => Math.Abs(p.ExperimentalMz - mz)).ToList(),
                                monoMass, z, corr, fracIntensityObserved);

                            subsetEnvelopeCandidates.Add(subsetEnvelope);
                        }
                    }
                }

                // use the best subset isotopic envelope
                env = subsetEnvelopeCandidates.OrderByDescending(p => p.Score).FirstOrDefault();
            }

            if (env != null)
            {
                HashSet<double> originalEnvelopeMzs = new HashSet<double>(env.Peaks.Select(p => p.ExperimentalMz));

                var sn = GetBaselineAndNoise(p, intensitiesBuffer, spectrum);
                env.Baseline = sn.baseline;
                env.NoiseFwhm = sn.noiseFwhm;

                var spectralAngle = env.GetNormalizedSpectralAngle(spectrum, averagineEnvelopeMasses, averagineEnvelopeIntensities, PpmTolerance, null);
                env.NormalizedSpectralAngle = spectralAngle;

                for (int i = 2; i < 6; i++)
                {
                    int harmonicCharge = i * env.Charge;
                    double harmonicMass = mz.ToMass(harmonicCharge);
                    int harmonicMassIndex = GetMassIndex(harmonicMass);

                    if (harmonicMassIndex >= allMasses.Length)
                    {
                        break;
                    }

                    var harmonicIsotopes = GetPotentialIsotopePeaks(deconvolutedPeaks, harmonicMassIndex, null, spectrum, harmonicCharge, harmonicMass, intensity);

                    double[] harmonicAveragineEnvelopeMasses = allMasses[harmonicMassIndex];
                    double[] harmonicAveragineEnvelopeIntensities = allIntensities[harmonicMassIndex];
                    double harmonicMonoMass = mass - diffToMonoisotopic[harmonicMassIndex];

                    var harmonicEnvelope = new DeconvolutedEnvelope(deconvolutedPeaks, harmonicMass, harmonicCharge, 0, 0);
                    harmonicEnvelope.NoiseFwhm = env.NoiseFwhm;
                    harmonicEnvelope.Baseline = env.Baseline;
                    double harmonicSpectralAngle = harmonicEnvelope.GetNormalizedSpectralAngle(spectrum, harmonicAveragineEnvelopeMasses, harmonicAveragineEnvelopeIntensities, PpmTolerance,
                        originalEnvelopeMzs, true);

                    if (!double.IsNaN(harmonicSpectralAngle))
                    {
                        env.InterstitialSpectralAngle = Math.Max(harmonicSpectralAngle, env.InterstitialSpectralAngle);
                    }
                }

                //if (env.SignalToNoise > 5 && env.FractionIntensityMissing > 0.3)
                //{
                //    return null;
                //}
            }

            return env;
        }

        public HashSet<int> GetPotentialChargeStates(HashSet<int> potentialChargeStates, MzSpectrum scan, int xArrayIndex, Polarity polarity)
        {
            // get rough list of charge states to check for based on m/z peaks around this peak
            potentialChargeStates.Clear();
            double mz = scan.XArray[xArrayIndex];

            for (int i = xArrayIndex + 1; i < scan.XArray.Length; i++)
            {
                double potentialIsotopeMz = scan.XArray[i];

                if (potentialIsotopeMz > mz + 1.2)
                {
                    break;
                }

                if (polarity == Polarity.Negative)
                {
                    for (int z = -MinCharge; z >= -MaxCharge; z--)
                    {
                        if (PpmTolerance.Within(potentialIsotopeMz.ToMass(z), mz.ToMass(z) + Constants.C13MinusC12))
                        {
                            potentialChargeStates.Add(z);
                        }
                    }
                }
                else
                {
                    for (int z = MinCharge; z <= MaxCharge; z++)
                    {
                        if (PpmTolerance.Within(potentialIsotopeMz.ToMass(z), mz.ToMass(z) + Constants.C13MinusC12))
                        {
                            potentialChargeStates.Add(z);
                        }
                    }
                }
            }

            return potentialChargeStates;
        }

        public IEnumerable<int> GetPeaksThatPassSignalToNoiseFilter(MzSpectrum spectrum)
        {
            HashSet<int> indexes = new HashSet<int>();
            List<double> intensities = new List<double>();

            double start = spectrum.Range.Minimum;
            double end = start + 100;

            for (; start < spectrum.Range.Maximum; start += 50)
            {
                intensities.Clear();
                end = start + 100;

                int ind = spectrum.GetClosestPeakIndex(start);
                int ind2 = spectrum.GetClosestPeakIndex(end);

                for (int i = ind; i < ind2; i++)
                {
                    intensities.Add(spectrum.YArray[i]);
                }

                var hdi = Util.GetHighestDensityInterval(intensities.ToArray(), 0.2);
                double noiseSigma = (hdi.hdi_end - hdi.hdi_start) * 1.9685;
                double noiseFwhm = noiseSigma * 2.3558;
                double baseline = (hdi.hdi_end + hdi.hdi_start) / 2;

                for (int i = ind; i < ind2; i++)
                {
                    double signalToNoise = (spectrum.YArray[i] - baseline) / noiseFwhm;

                    if (signalToNoise > SignalToNoiseRequired)
                    {
                        indexes.Add(i);
                    }
                }
            }

            return indexes;
        }

        public List<DeconvolutedPeak> GetPotentialIsotopePeaks(List<DeconvolutedPeak> deconvolutedPeaks, int massIndex, HashSet<double> alreadyClaimedMzs,
            MzSpectrum spectrum, int z, double modeMass, double modeIntensity)
        {
            deconvolutedPeaks.Clear();

            double[] averagineEnvelopeMasses = allMasses[massIndex];
            double[] averagineEnvelopeIntensities = allIntensities[massIndex];
            double monoMass = modeMass - diffToMonoisotopic[massIndex];

            int indOfMostIntense = Array.IndexOf(averagineEnvelopeIntensities, 1);

            // 1 is to the right, -1 is to the left in the envelope
            int isotopeDirection = 1;
            for (int i = indOfMostIntense; i < averagineEnvelopeMasses.Length && i >= 0; i += isotopeDirection)
            {
                double isotopeMassShift = averagineEnvelopeMasses[i] - averagineEnvelopeMasses[indOfMostIntense];
                double isotopeTheoreticalMass = modeMass + isotopeMassShift;
                double theoreticalIsotopeMz = isotopeTheoreticalMass.ToMz(z);
                double theoreticalIsotopeIntensity = averagineEnvelopeIntensities[i] * modeIntensity;

                var peakIndex = spectrum.GetClosestPeakIndex(theoreticalIsotopeMz);

                //TODO: look for other peaks in the scan that could be this isotope that meet the m/z tolerance
                var isotopeExperMz = spectrum.XArray[peakIndex];
                var isotopeExperIntensity = spectrum.YArray[peakIndex];

                double intensityRatio = isotopeExperIntensity / theoreticalIsotopeIntensity;

                double isotopeExperimentalMass = isotopeExperMz.ToMass(z);
                bool withinMassTol = PpmTolerance.Within(isotopeExperMz.ToMass(z), isotopeTheoreticalMass);
                bool withinIntensityTol = intensityRatio < IntensityRatioLimit && intensityRatio > 1 / IntensityRatioLimit;
                bool unclaimedMz = (alreadyClaimedMzs == null || !alreadyClaimedMzs.Contains(isotopeExperMz))
                    && !deconvolutedPeaks.Select(p => p.ExperimentalMz).Contains(isotopeExperMz);

                if (withinMassTol // check mass tolerance
                    && withinIntensityTol // check intensity tolerance
                    && unclaimedMz) // check to see if this peak has already been claimed by another envelope or this envelope
                {
                    deconvolutedPeaks.Add(new DeconvolutedPeak(isotopeExperMz, theoreticalIsotopeMz, z, isotopeExperIntensity,
                        theoreticalIsotopeIntensity, i, averagineEnvelopeIntensities[i]));
                }
                else
                {
                    if (isotopeDirection == 1)
                    {
                        isotopeDirection = -1;
                        i = indOfMostIntense;
                    }
                    else
                    {
                        break;
                    }
                }
            }

            return deconvolutedPeaks;
        }

        public double GetModeMassFromMonoisotopicMass(double monoMass)
        {
            int index = GetMassIndex(monoMass) + 1;

            if (index >= allMasses.Length)
            {
                index = allMasses.Length - 1;
            }

            var masses = allMasses[index];
            int indexOfMostAbundant = Array.IndexOf(masses, 1);

            double modeMass = monoMass + diffToMonoisotopic[index];

            return modeMass;
        }

        private int GetMassIndex(double mass)
        {
            int massIndex = Array.BinarySearch(mostIntenseMasses, mass);

            if (massIndex < 0)
            {
                massIndex = ~massIndex;
            }

            return massIndex;
        }

        private DeconvolutedEnvelope GetNextBestEnvelope(List<DeconvolutedEnvelope> envelopeCandidates, DeconvolutedEnvelope previousBestEnvelope)
        {
            if (previousBestEnvelope == null)
            {
                return envelopeCandidates.OrderByDescending(p => p.Score).FirstOrDefault();
            }

            List<DeconvolutedEnvelope> sameMassEnvelopes = new List<DeconvolutedEnvelope>();

            for (int i = -3; i <= 3; i++)
            {
                double monoMass = i * Constants.C13MinusC12 + previousBestEnvelope.MonoisotopicMass;
                sameMassEnvelopes.AddRange(envelopeCandidates.Where(p => PpmTolerance.Within(p.MonoisotopicMass, monoMass)));
            }

            if (sameMassEnvelopes.Any())
            {
                return sameMassEnvelopes.OrderByDescending(p => p.Score).FirstOrDefault();
            }

            return envelopeCandidates.OrderByDescending(p => p.Score).FirstOrDefault();
        }

        public (double baseline, double noiseFwhm) GetBaselineAndNoise(int index, List<(double, double)> intensitiesBuffer, MzSpectrum spectrum)
        {
            double noiseHalfWindowWidth = 50;

            double mz = spectrum.XArray[index];
            intensitiesBuffer.Clear();

            for (int j = index; j < spectrum.XArray.Length; j++)
            {
                if (spectrum.XArray[j] > mz + noiseHalfWindowWidth)
                {
                    break;
                }

                intensitiesBuffer.Add((spectrum.XArray[j], spectrum.YArray[j]));
            }
            for (int j = index - 1; j >= 0; j--)
            {
                if (spectrum.XArray[j] < mz - noiseHalfWindowWidth)
                {
                    break;
                }

                intensitiesBuffer.Add((spectrum.XArray[j], spectrum.YArray[j]));
            }

            var hdi = Util.GetHighestDensityInterval(intensitiesBuffer.Select(p => p.Item2).ToArray(), 0.2);
            double noiseSigma = (hdi.hdi_end - hdi.hdi_start) * 1.9685;
            double noiseFwhm = noiseSigma * 2.3558;
            double baseline = (hdi.hdi_end + hdi.hdi_start) / 2;

            return (baseline, noiseFwhm);
        }
    }
}
