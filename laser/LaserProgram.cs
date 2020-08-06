using System;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Numerics;
using DSPLib;

namespace laser
{
    /// <summary>
    /// Class for processing raw data to detect pulses
    /// Program implemented in Main
    /// 
    /// Arguments:
    /// -inputfile:     Input to the raw data
    /// -outputfile:    Output CSV to the base band signal
    /// -carrier:       Carrier modulation frequency of the collected data
    /// -sampling:      Sampling rate of the modulated signal
    /// -passband:      Kaiser window FIR low pass filter band
    /// -cutoffband:    Kaiser window FIR low pass filter cutoff band to -60dB
    /// -gaindev:       Kaiser window FIR low pass filter gain deviation from unity in the pass band
    /// -peaklag:       Window length for peak detection statistics
    /// -peakthreshold: Threshold used to detect a peak in the baseband signal
    /// -peakinfluence: Value used to detect a peak in the baseband signal as new measurements arrive
    /// -tifile:        Output CSV file where the time intervals in between peaks are stored
    /// -tithreshold:   Threshold used to detect an outlier in the time intervals series
    /// 
    ///  The format of the baseband signal output file is:
    ///  (sample, baseband value)
    ///  
    ///  The format of the time intervals output file is:
    ///  (time stamp of peak, time interval between successive peaks, upper limit, lower limit, is time interval an outlier Y/N)
    /// </summary>
    class LaserProgram
    {
        private const string sDisable = "disable";
        private const string sEnable = "enable";

        /// <summary>
        /// Class encapsulating writing to the output file
        /// </summary>
        class WriteToFile
        {
            public StreamWriter outFileHandler { get; }
            public WriteToFile(string outputFile)
            {
                outFileHandler = new StreamWriter(outputFile);
            }

            public void Data(double timeStamp, double timePulseDifference)
            {
                outFileHandler.WriteLine("{0},{1}", timeStamp, timePulseDifference);
            }

            public void SummaryOutliers(PeakDetectionOutput detector, List<int> outliers, double avg, double stdDev, double thrs)
            {
                for (int i=0; i<detector.peakTimeStamps.Count; i++)
                {
                    if (outliers.Contains(i))
                    {
                        outFileHandler.WriteLine(
                            "{0},{1},{2},{3},Y",
                            detector.peakTimeStamps.ElementAt(i).ToString("n5"),
                            detector.peakTimeIntervals.ElementAt(i).ToString("n5"),
                            (avg + thrs * stdDev).ToString("n5"),
                            (avg - thrs * stdDev).ToString("n5")
                            );
                    }
                    else
                    {
                        outFileHandler.WriteLine(
                            "{0},{1},{2},{3},N",
                            detector.peakTimeStamps.ElementAt(i).ToString("n5"),
                            detector.peakTimeIntervals.ElementAt(i).ToString("n5"),
                            (avg + thrs * stdDev).ToString("n5"),
                            (avg - thrs * stdDev).ToString("n5")
                            );
                    }
                }
            }
        }

        /// <summary>
        /// Parse program arguments
        /// </summary>
        /// <param name="args">= Input program arguments</param>
        /// <returns></returns>
        static Dictionary<string, string> Parse(string[] args)
        {
            // List the program flags
            Dictionary<string, string> progArgs = new Dictionary<string, string>
            {
                { "inputfile", sDisable },
                { "outputfile", sDisable },
                { "tifile", sDisable },
                { "carrier", sDisable },
                { "sampling", sDisable },
                { "passband", sDisable },
                { "cutoffband", sDisable },
                { "gaindev", sDisable },
                { "peaklag", sDisable },
                { "peakthreshold", sDisable },
                { "peakinfluence", sDisable },
                { "tithreshold", sDisable }
            };

            Console.WriteLine($"Parameter count = {args.Length}");
            string sLastUsedOption = "";
            for (int i = 0; i < args.Length; i++)
            {
                string[] argParse = args[i].Replace("-/", "-").Split('-');
                string sOption = argParse[argParse.Length - 1];
                // Argument directive
                if (progArgs.ContainsKey(sOption))
                {
                    // Invalidate last option
                    if (sLastUsedOption != "")
                    {
                        if (progArgs[sLastUsedOption] == sEnable)
                        {
                            progArgs[sLastUsedOption] = sDisable;
                        }
                    }
                    progArgs[sOption] = sEnable;
                    sLastUsedOption = sOption;
                    continue;
                }
                // Argument value
                else
                {
                    if (progArgs.ContainsKey(sLastUsedOption))
                    {
                        progArgs[sLastUsedOption] = sOption;
                    }
                }
            }
            return progArgs;
        }

        /// <summary>
        /// Print base band signal
        /// </summary>
        /// <param name="fHandler"></param>
        /// <param name="baseBandSignal"></param>
        static void SignalOutput(WriteToFile fHandler, List<Complex> ordinate)
        {
            for (int i = 0; i < ordinate.Count; i++)
            {
                fHandler.Data(i, ordinate.ElementAt(i).Magnitude);
            }
        }

        /// <summary>
        /// Main entry point
        /// </summary>
        /// <param name="args">= Input program arguments</param>
        static void Main(string[] args)
        {
            try
            {
                Dictionary<string, string> progArgs = Parse(args);
                // Program arguments
                double fs = Convert.ToDouble(progArgs["sampling"]);
                double Ts = 1.0 / fs;
                double fc = Convert.ToDouble(progArgs["carrier"]);
                double passBand = Convert.ToDouble(progArgs["passband"]);
                double stopBand = Convert.ToDouble(progArgs["cutoffband"]);
                double gainDeviation = Convert.ToDouble(progArgs["gaindev"]);
                int peakLag = Convert.ToInt32(progArgs["peaklag"]);
                double peakThreshold = Convert.ToDouble(progArgs["peakthreshold"]);
                double peakInfluence = Convert.ToDouble(progArgs["peakinfluence"]);
                double tiThreshold = Convert.ToDouble(progArgs["tithreshold"]);
                Dsp Signal = new Dsp();
                
                // Low pass FIR filter 
                List<Complex> firImpulseResponse = Signal.KaiserWindowFir(gainDeviation, passBand, stopBand, fs);
                // FIR filter length
                int M = firImpulseResponse.Count;
                // Block length
                int N = Signal.MinBlockLength(M);
                // FFT input sequence length
                int L = N - M + 1;
                Complex[] appendSequence = new Complex[L - 1];
                Array.Clear(appendSequence, 0, L - 1);
                // FIR filter frequency response
                Complex[] firFFT = Signal.ComputePaddedFFT(ref firImpulseResponse, N, appendSequence.ToList());
                // Main class output values
                List<double> timeStamps = new List<double>();
                List<double> peakTimeStamps = new List<double>();
                List<double> peakTimeIntervals = new List<double>();
                List<double> timeIntervalsMovingAverage = new List<double>();
                List<int> peakIndexes = new List<int>();
                List<int> outliersTimeIntervals = null;

                // Read data from input file
                using StreamReader inputFile = new StreamReader(progArgs["inputfile"]);
                WriteToFile outFile = new WriteToFile(progArgs["outputfile"]);
                WriteToFile timeIntervalFile = new WriteToFile(progArgs["tifile"]);
                inputFile.BaseStream.Seek(0, SeekOrigin.Begin);
                string str = inputFile.ReadLine();
                double dLSBValue = -1;
                double dVoltageValue = -1;
                double dValue = 0.0;
                double dAverage = 0.0;
                double dTimeStamp = 0.0;
                double dPhase = 0.0;
                double timeIntervalsMoveAverage = 0.0;
                double timeIntervalsMoveVariance = 0.0;
                double lastPeakTimeStamp = 0.0;
                int iNumSamples = 1;
                int numDetections = 0;
                int numMovingAverageSamples = 0;
                while (str != null)
                {
                    string[] sDataValues = str.Split(',');
                    try
                    {
                        dTimeStamp = ((double)iNumSamples - 1.0) * Ts;
                        dLSBValue = Convert.ToDouble(sDataValues[0]);
                        dVoltageValue = Convert.ToDouble(sDataValues[1]);
                        dAverage = ((double)iNumSamples - 1.0) / ((double)iNumSamples) * dAverage + dVoltageValue / ((double)iNumSamples);
                        // Remove DC offset
                        dValue = dVoltageValue - dAverage;
                        // Extract I and Q from signal and get base band signal
                        bool bPeakDetection = Signal.AcquireSignal(dValue, dPhase, dTimeStamp, fc, L, ref firFFT);
                        // Perform peak detection every L samples of the input sequence found from FFT transform
                        if (bPeakDetection)
                        {
                            PeakDetectionOutput detector = PeakDetection.PeakDetectionAlgorithm(
                                Signal.baseBandSignal.GetRange(numDetections * L, L),
                                Signal.timeStamp.GetRange(numDetections * L, L),
                                peakLag,
                                peakThreshold,
                                peakInfluence,
                                lastPeakTimeStamp,
                                numDetections * L
                                );
                            // Shift by L samples
                            numDetections += 1;
                            // Compute moving average of time intervals
                            numMovingAverageSamples += detector.peakTimeIntervals.Count;
                            outliersTimeIntervals = PeakDetection.TimeIntervalOutlierDetection(
                                detector.peakTimeIntervals,
                                numMovingAverageSamples,
                                timeIntervalsMoveAverage,
                                ref timeIntervalsMoveVariance,
                                tiThreshold
                                );
                            double currentWindowSum = detector.peakTimeIntervals.Sum();
                            timeIntervalsMoveAverage =
                                (numMovingAverageSamples - detector.peakTimeIntervals.Count) * 
                                timeIntervalsMoveAverage / numMovingAverageSamples + 
                                currentWindowSum / numMovingAverageSamples;
                            // Time stamp for the next block of data to find the time intervals
                            if (detector.peakTimeStamps.Count > 0) { lastPeakTimeStamp = detector.peakTimeStamps.Last(); }
                            // Store outputs either list or write to a file
                            timeIntervalFile.SummaryOutliers(
                                detector, 
                                outliersTimeIntervals,
                                timeIntervalsMoveAverage,
                                Math.Sqrt(timeIntervalsMoveVariance),
                                tiThreshold
                                );
                            timeIntervalsMovingAverage.Add(timeIntervalsMoveAverage);
                            peakTimeStamps.AddRange(detector.peakTimeStamps);
                            peakIndexes.AddRange(detector.peakIndexes);
                            peakTimeIntervals.AddRange(detector.peakTimeIntervals);
                        }
                        timeStamps.Add(dTimeStamp);
                        iNumSamples += 1;
                    }
                    catch (FormatException)
                    {

                    }
                    catch (IndexOutOfRangeException)
                    {

                    }
                    str = inputFile.ReadLine();
                }
                Signal.LastSamplesConvolution(ref firFFT, L);
                SignalOutput(outFile, Signal.baseBandSignal);
                outFile.outFileHandler.Close();
                timeIntervalFile.outFileHandler.Close();
                Console.WriteLine("End of Test");
            }
            catch (Exception e)
            {
                Console.WriteLine(e);
            }
        }
    }
}
