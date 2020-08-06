using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;


namespace DSPLib
{
    /// <summary>
    /// Digital signal processing class
    /// </summary>
    public class Dsp
    {
        private List<Complex> oldComplexSamples;
        private List<Complex> complexSignal;
        private Complex[] overlapSamples;
        public List<Complex> baseBandSignal;
        public List<double> timeStamp;

        /// <summary>
        /// Factorial of number
        /// </summary>
        /// <param name="k">= Integer number to compute factorial</param>
        /// <returns></returns>
        private double Factorial(double k)
        {
            if (k == 0)
            {
                return 1;
            }
            else
            {
                return k * Factorial(k - 1);
            }
        }

        /// <summary>
        /// Fast fourier transform of a complex sequence
        /// </summary>
        /// <param name="x">= Complex sequence</param>
        public void FastFourierTransform(ref Complex[] x)
        {
            // Find the power of two for the total FFT size up to 2^32
            bool foundIt = false;
            for (int powerTwo = 0; powerTwo <= 32; powerTwo++)
            {
                double n = Math.Pow(2.0, powerTwo);
                if ((x.Length) == n)
                {
                    foundIt = true;
                    break;
                }
            }
            if (foundIt == false)
                throw new ArgumentOutOfRangeException("Length of sequence was not an even power of 2! FFT cannot continue.");
            int N = x.Length;
            if (N <= 1)
            {
                return;
            }
            // divide
            Complex[] even = new Complex[N / 2];
            Complex[] odd = new Complex[N / 2];
            for (int i = 0; i < N / 2; i++)
            {
                even[i] = x[2 * i];
                odd[i] = x[2 * i + 1];
            }
            // Conquer
            FastFourierTransform(ref even);
            FastFourierTransform(ref odd);
            // Combine
            for (int k = 0; k < N / 2; ++k)
            {
                Complex t = Complex.FromPolarCoordinates(1.0, -2 * Math.PI * k / N) * odd[k];
                x[k] = even[k] + t;
                x[k + N / 2] = even[k] - t;
            }
        }

        /// <summary>
        /// Complex conjugate of array of complex numbers
        /// </summary>
        /// <param name="x">= Array of complex numbers</param>
        /// <returns></returns>
        private Complex[] ArrayConjugate(Complex[] x)
        {
            Complex[] conjugateX = new Complex[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                conjugateX[i] = Complex.Conjugate(x[i]);
            }
            return conjugateX;
        }

        /// <summary>
        /// Multiplication of two array of complex numbers
        /// </summary>
        /// <param name="array1">= Array of complex numbers</param>
        /// <param name="array2">= Array of complex numbers</param>
        /// <returns></returns>
        private Complex[] ArrayMultiply(Complex[] array1, Complex[] array2)
        {
            if (array1.Length != array2.Length)
            {
                throw new ArgumentOutOfRangeException("Vectors of complex numbers sizes must match");
            }
            Complex[] complexMult = new Complex[array1.Length];
            for (int i = 0; i < array1.Length; i++)
            {
                complexMult[i] = array1[i] * array2[i];
            }
            return complexMult;
        }

        /// <summary>
        /// Computes the IFFT of a frequency domain sequence
        /// </summary>
        /// <param name="x">= Sequence in the frequency domain</param>
        private void InverseFFT(ref Complex[] x)
        {
            // Conjugate the complex numbers
            x = ArrayConjugate(x);
            // forward fft
            FastFourierTransform(ref x);
            // Conjugate the complex numbers again
            x = ArrayConjugate(x);
            // Scale the numbers
            for (int i = 0; i < x.Length; i++)
            {
                x[i] = x[i] / x.Length;
            }
        }

        /// <summary>
        /// Circular convolution
        /// </summary>
        /// <param name="firFFT">= FIR filter discrete time fourier transform</param>
        /// <param name="sequence">= Received sequence block</param>
        /// <param name="signalOffset">= Offset in received signal</param>
        /// <param name="L"></param>
        /// <returns></returns>
        private List<Complex> CircularConvolution(
            ref List<Complex> firFFT,
            ref List<Complex> sequence,
            ref Complex[] overlapSequence,
            List<Complex> appendSequence = null
            )
        {
            int N = firFFT.Count;
            Complex[] sequenceFFT = ComputePaddedFFT(ref sequence, N, appendSequence, overlapSequence.ToList());
            Complex[] circConv = ArrayMultiply(sequenceFFT, firFFT.ToArray());
            InverseFFT(ref circConv);
            return circConv.ToList();
        }

        /// <summary>
        /// Constructor for digital signal processing
        /// </summary>
        public Dsp()
        {
            oldComplexSamples = new List<Complex>();
            complexSignal = new List<Complex>();
            baseBandSignal = new List<Complex>();
            overlapSamples = null;
            timeStamp = new List<double>();
        }

        /// <summary>
        /// Computes convolution for remaining sequence samples
        /// </summary>
        /// <param name="firFilter">= Filter finite impulse response</param>
        /// <param name="L">= Length of window sequence</param>
        public void LastSamplesConvolution(ref Complex[] firFilter, int L)
        {
            if (complexSignal.Count > 0)
            {
                int M = firFilter.Length - L + 1;
                overlapSamples = oldComplexSamples.GetRange(0, M - 1).ToArray();
                Complex[] appendSamples = new Complex[L - complexSignal.Count];
                Array.Clear(appendSamples, 0, L - complexSignal.Count);
                List<Complex> firFilterList = firFilter.ToList();
                List<Complex> baseBandSig = CircularConvolution(ref firFilterList, ref complexSignal, ref overlapSamples, appendSamples.ToList());
                baseBandSignal.AddRange(baseBandSig.GetRange(M - 1, L));
                // Remove all collected samples once the data points are collected and processed
                complexSignal.Clear();
            }
        }

        /// <summary>
        /// Find the minimum length of block of data to perform circular convolution
        /// </summary>
        /// <param name="firFilterLength">= Length of FIR filter impulse response</param>
        /// <returns></returns>
        public int MinBlockLength(int firFilterLength)
        {
            int N = 0;
            int i = 0;
            while (true)
            {
                if (N > firFilterLength)
                {
                    int M = N - firFilterLength + 1;
                    if (M >= firFilterLength)
                    {
                        break;
                    }
                }
                i = i + 1;
                N = (int)Math.Pow(2, i);
            }
            return N;
        }

        public Complex[] ComputePaddedFFT(
            ref List<Complex> sequence,
            int totalLength,
            List<Complex> appendSequence = null,
            List<Complex> prependSequence = null
            )
        {
            int origSeqLength = sequence.Count;
            int beforeSeqLength = 0;
            Complex[] sequenceFFT = new Complex[totalLength];
            Array.Clear(sequenceFFT, 0, totalLength);
            if (prependSequence != null)
            {
                beforeSeqLength = prependSequence.Count;
                prependSequence.CopyTo(sequenceFFT, 0);
            }
            if (appendSequence != null)
            {
                appendSequence.CopyTo(sequenceFFT, beforeSeqLength + origSeqLength);
            }
            sequence.CopyTo(sequenceFFT, beforeSeqLength);
            FastFourierTransform(ref sequenceFFT);
            return sequenceFFT;
        }

        /// <summary>
        /// Coherent signal acquisition from received modulated signal
        /// </summary>
        /// <param name="dSample">= Incoming signal sample</param>
        /// <param name="dPhase">= Phase of incoming singal sample</param>
        /// <param name="dTimeStamp">= Time stamp of incoming signal sample</param>
        /// <param name="dCarrier">= Carrier of incoming signal sample</param>
        /// <param name="L">= Block size of signal</param>
        /// <param name="M">= Length of FIR filter impulse response</param>
        /// <param name="firFilter">FFT transform of FIR impulse response</param>
        public bool AcquireSignal(
            double dSample,
            double dPhase,
            double dTimeStamp,
            double dCarrier,
            int L,
            ref Complex[] firFilter
            )
        {
            bool bDemodulated = false;
            // Get In-phase and quadrature phase components of mixed signal for demodulation
            double inPhase = dSample * Math.Cos(2 * Math.PI * dCarrier * dTimeStamp + dPhase);
            double quadPhase = dSample * Math.Sin(2 * Math.PI * dCarrier * dTimeStamp + dPhase);
            int M = firFilter.Length - L + 1;
            complexSignal.Add(new Complex(inPhase, quadPhase));
            timeStamp.Add(dTimeStamp);
            // Perform a real time convolution FFT
            if (complexSignal.Count % L == 0)
            {
                // First time append zeros
                if (overlapSamples == null)
                {
                    overlapSamples = new Complex[M - 1];
                    Array.Clear(overlapSamples, 0, M - 1);
                }
                else
                {
                    int lowIndex = L - M + 1;
                    overlapSamples = complexSignal.GetRange(lowIndex, M - 1).ToArray();
                }
                List<Complex> firFilterList = firFilter.ToList();
                List<Complex> baseBandSig = CircularConvolution(ref firFilterList, ref complexSignal, ref overlapSamples);
                baseBandSignal.AddRange(baseBandSig.GetRange(M - 1, L));
                // Remove all collected samples once the data points are collected and processed
                oldComplexSamples.Clear();
                oldComplexSamples.AddRange(complexSignal.GetRange(L - M + 1, M - 1));
                complexSignal.Clear();
                bDemodulated = true;
            }
            return bDemodulated;
        }

        /// <summary>
        /// Modified bessel function of the first kind used in Kaiser FIR filters
        /// The sum is truncated when the increment is very small compared with
        /// the magnitude of the Bessel value
        /// </summary>
        /// <param name="x">=Argument to Bessel function</param>
        /// <param name="nu">=Order of the Bessel function</param>
        /// <returns></returns>
        public double ModifiedBesselFirstKind(double x, double nu, double tol = 1e-5)
        {
            double sum = 0;
            double den;
            double num;
            int k = 0;
            double bessel = Math.Pow(0.5 * x, nu);
            while (true)
            {
                num = Math.Pow((0.25 * Math.Pow(x, 2)), k);
                den = Factorial(k) * Factorial((int)nu + k);
                // Stop when the increments are small enough
                if (Math.Abs(num / den) < sum * tol)
                {
                    break;
                }
                sum += (num / den);
                k += 1;
            }
            bessel *= sum;
            return bessel;
        }

        /// <summary>
        /// Kaiser window design for FIR band pass filter
        /// </summary>
        /// <param name="delta">= Gain deviation from unity</param>
        /// <param name="omegap">= Low pass band frequency</param>
        /// <param name="omegas">= High stop band frequency</param>
        /// <param name="fs">= Sampling frequency</param>
        /// <returns></returns>
        public List<Complex> KaiserWindowFir(double delta, double omegap, double omegas, double fs)
        {
            List<Complex> h = new List<Complex>();
            double bandWidth = 2 * Math.PI * (omegas - omegap) / fs;
            double omegac = Math.PI * (omegap + omegas) / fs;
            double A = -20 * Math.Log10(delta);
            double beta = 0.0;
            // Length of FIR impulse response
            double fracM = (((A - 8) / (2.285 * bandWidth)) - (int)((A - 8) / (2.285 * bandWidth)));
            int plusFracM = (int)(Math.Ceiling(fracM));
            int M = (int)((A - 8) / (2.285 * bandWidth)) + plusFracM;
            if (A > 50)
            {
                beta = 0.1102 * (A - 8.7);
            }
            else if (21 <= A && 50 >= A)
            {
                beta = 0.5842 * Math.Pow(A - 21, 0.4) + 0.07886 * (A - 21);
            }
            double besselDen = ModifiedBesselFirstKind(beta, 0);
            // Impulse response of the FIR filter
            for (int n = 0; n <= M; n++)
            {
                double alpha = (double)M / 2.0;
                double besselArg = beta * Math.Pow(1 - Math.Pow((n - alpha) / alpha, 2), 0.5);
                double besselNum = ModifiedBesselFirstKind(besselArg, 0);
                double hdNum = Math.Sin(omegac * (n - alpha));
                double hdDen = Math.PI * (n - alpha);
                // Limit case of sinc function
                if (n == alpha)
                {
                    hdNum = omegac;
                    hdDen = Math.PI;
                }
                h.Add(new Complex(hdNum / hdDen * besselNum / besselDen, 0.0));
            }
            return h;
        }
    }

    public class PeakDetectionOutput
    {
        public List<double> basebandSamples;
        public List<int> peakIndexes;
        public List<double> peakTimeStamps;
        public List<double> peakTimeIntervals;
        public List<double> avgFilter;
        public List<double> stdDevFilter;

        /// <summary>
        /// Gets the magnitude of the I + jQ base band signal samples
        /// </summary>
        /// <param name="baseBandSignal">= Complex base band demodulated signal</param>
        /// <returns>= Magnitude of base band demodulated signal</returns>
        public static List<double> MagnitudeOfComplex(List<Complex> baseBandSignal)
        {
            List<double> baseBandSamples = new List<double>();
            for (int i=0; i<baseBandSignal.Count; i++)
            {
                baseBandSamples.Add(baseBandSignal[i].Magnitude);
            }
            return baseBandSamples;
        }
    }

    public static class PeakDetection
    {
        /// <summary>
        /// Compute the mean of a window list
        /// </summary>
        /// <param name="windowList">= data window</param>
        /// <returns></returns>
        private static double Mean(List<double> windowList)
        {
            return windowList.Average();
        }

        /// <summary>
        /// Compute the standard deviation of a window list
        /// </summary>
        /// <param name="windowValues">= data window</param>
        /// <returns></returns>
        private static double StdDev(List<double> windowValues)
        {
            double ret = 0;
            if (windowValues.Count() > 1)
            {
                double avg = windowValues.Average();
                double sum = windowValues.Sum(d => Math.Pow(d - avg, 2));
                ret = Math.Sqrt((sum) / (windowValues.Count() - 1));
            }
            return ret;
        }

        /// <summary>
        /// Algorithm to detect peak values in a sequence
        /// From robust peak detection algorithm forum in stack overflow (www.stackoverflow.com)
        /// based on the principle of dispersion where the data standard deviation is modified
        /// for detection of outliers
        /// </summary>
        /// <param name="inputSignal">= Complex value base band signal</param>
        /// <param name="timeStamps">= Time stamp of the sequence</param>
        /// <param name="lag">= Window lag</param>
        /// <param name="threshold">= Used to detect a peak in sequence in standard deviations</param>
        /// <param name="influence">= between 0 and 1, where 1 is normal influence, 0.5 is half</param>
        /// <returns>Class PeakDetectionOutput with time stamps where the peak happens</returns>
        public static PeakDetectionOutput PeakDetectionAlgorithm(
            List<Complex> inputSignal,
            List<double> timeStamps,
            int lag,
            double threshold,
            double influence,
            double lastPeakTimeStamp = 0,
            int sequenceOffset=0
            )
        {
            var result = new PeakDetectionOutput();
            List<double> peakTimeStamps = new List<double>();
            List<double> peakTimeIntervals = new List<double>();
            List<int> peakIndexes = new List<int>();
            List<int> peakRecords = new List<int>();
            double[] inputSamples = PeakDetectionOutput.MagnitudeOfComplex(inputSignal).ToArray();
            double[] filteredY = new List<double>(inputSamples.ToList()).ToArray();
            double[] avgFilter = new double[inputSignal.Count];
            double[] stdFilter = new double[inputSignal.Count];
            double peakRecord = 0;
            int peakIndex = -1;
            List<double> slidingWindow = inputSamples.ToList().GetRange(0, lag);
            avgFilter[lag - 1] = Mean(slidingWindow);
            stdFilter[lag - 1] = StdDev(slidingWindow);
            // Real time peak detection algorithm
            for (int i=lag; i<inputSignal.Count; i++)
            {
                // Evaluated sample represents a potential peak
                if (Math.Abs(inputSamples[i] - avgFilter[i - 1]) > threshold * stdFilter[i - 1])
                {
                    peakIndexes.Add(i);
                    // Find the maximum value of the peak in contiguous recorded peaks
                    if (peakIndexes.Count > 1)
                    {
                        // Condition for contiguous potential detected peaks
                        if ( ((peakIndexes.ElementAt(peakIndexes.Count - 1) - peakIndexes.ElementAt(peakIndexes.Count - 2)) <= 1) )
                        {
                            // Record the largest peak value and its index
                            if (inputSamples[peakIndexes.ElementAt(peakIndexes.Count - 1)] >= peakRecord)
                            {
                                peakRecord = inputSamples[peakIndexes.ElementAt(peakIndexes.Count - 1)];
                                peakIndex = i;
                            }
                        }
                        else
                        {
                            peakTimeStamps.Add(timeStamps.ElementAt(peakIndex));
                            peakRecords.Add(peakIndex + sequenceOffset);
                            // Compute time interval difference including from the block of data from before
                            if (peakTimeStamps.Count > 1)
                            {
                                peakTimeIntervals.Add(
                                    peakTimeStamps.ElementAt(peakTimeStamps.Count - 1) - 
                                    peakTimeStamps.ElementAt(peakTimeStamps.Count - 2)
                                    );
                            }
                            else
                            {
                                peakTimeIntervals.Add(
                                    peakTimeStamps.ElementAt(peakTimeStamps.Count - 1) - lastPeakTimeStamp
                                    );
                            }
                            peakRecord = inputSamples[i];
                            peakIndex = i;
                        }
                    }
                    else
                    {
                        peakIndex = i;
                    }
                    filteredY[i] = influence * inputSamples[i] + (1 - influence) * filteredY[i - 1];
                }
                else
                {
                    filteredY[i] = inputSamples[i];
                }
                slidingWindow = filteredY.ToList().GetRange(i - lag, lag + 1);
                avgFilter[i] = Mean(slidingWindow);
                stdFilter[i] = StdDev(slidingWindow);
            }
            // Last peak record
            if (peakIndex > 0)
            {
                peakTimeStamps.Add(timeStamps.ElementAt(peakIndex));
                peakRecords.Add(peakIndex + sequenceOffset);
                // Compute last peak interval
                if (peakTimeStamps.Count > 1)
                {
                    peakTimeIntervals.Add(
                        peakTimeStamps.ElementAt(peakTimeStamps.Count - 1) -
                        peakTimeStamps.ElementAt(peakTimeStamps.Count - 2)
                        );
                }
                else
                {
                    peakTimeIntervals.Add(
                        peakTimeStamps.ElementAt(peakTimeStamps.Count - 1) - lastPeakTimeStamp
                        );
                }
            }
            // Output class data
            result.basebandSamples = inputSamples.ToList();
            result.avgFilter = avgFilter.ToList();
            result.stdDevFilter = stdFilter.ToList();
            result.peakTimeStamps = peakTimeStamps;
            result.peakTimeIntervals = peakTimeIntervals;
            result.peakIndexes = peakRecords;
            return result;
        }

        /// <summary>
        /// Detect whether a time interval is an outlier
        /// </summary>
        /// <param name="timeIntervals">= Data with time intervals in between detected peaks</param>
        /// <param name="cumulativeNumSamples">= Total number of samples accumulated</param>
        /// <param name="timeMoveAverage">= Old time intervals moving average</param>
        /// <param name="timeMoveVariance">= Old time intervals moving variance</param>
        /// <param name="threshold">= Number of standard deviations for considering a data point an outlier</param>
        /// <returns></returns>
        public static List<int> TimeIntervalOutlierDetection(
            List<double> timeIntervals,
            int cumulativeNumSamples,
            double timeMoveAverage,
            ref double timeMoveVariance,
            double threshold
            )
        {
            List<int> outliersIndexList = new List<int>();
            double sumVector = timeIntervals.Sum();
            double sumSquaresVector = timeIntervals.Sum(x => x * x);
            int n = cumulativeNumSamples;
            int M = timeIntervals.Count;
            // Compute the moving variance of the sequence of time interval data
            if (cumulativeNumSamples >= timeIntervals.Count + 1)
            {
                timeMoveVariance = (n - M - 1) / (n - 1) * timeMoveVariance +
                    (M * (n - M) * Math.Pow(timeMoveAverage, 2) -
                     2 * (n - M) * timeMoveAverage * sumVector +
                     n * sumSquaresVector - Math.Pow(sumVector, 2)) / (n * (n - 1));
                // Search for outlier indexes in the current block of data
                for (int i = 0; i < timeIntervals.Count; i++)
                {
                    if (timeIntervals.ElementAt(i) > (timeMoveAverage + threshold * Math.Sqrt(timeMoveVariance)) ||
                         timeIntervals.ElementAt(i) < (timeMoveAverage - threshold * Math.Sqrt(timeMoveVariance)))
                    {
                        outliersIndexList.Add(i);
                    }
                }
            }
            else
            {
                timeMoveVariance = Math.Pow(StdDev(timeIntervals), 2);
            }
            return outliersIndexList;
        }
    }
}
