using Microsoft.VisualStudio.TestPlatform.TestHost;
using System;
using Xunit;
using DSPLib;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;

namespace lasertest
{
    public class LaserProgramUnitTests
    {
        private Dsp Signal;
        
        public LaserProgramUnitTests()
        {
            Signal = new Dsp();
        }

        /// <summary>
        /// Converts real sequence to Complex sequence, imaginary part is zero
        /// </summary>
        /// <param name="sequence">= Real sequence</param>
        /// <returns></returns>
        private Complex[] ToComplex(double[] sequence)
        {
            Complex[] complexSequence = new Complex[sequence.Length];
            for (int i = 0; i < sequence.Length; i++)
            {
                Complex number = new Complex(sequence[i], 0);
                complexSequence[i] = number;
            }
            return complexSequence;
        }

        /// <summary>
        /// Extract the real part of a complex array
        /// </summary>
        /// <param name="sequence">= Complex array</param>
        /// <returns>= Real part of complex array</returns>
        private double[] RealPart(Complex[] sequence)
        {
            double[] realPart = new double[sequence.Length];
            for (int i=0; i<sequence.Length; i++)
            {
                realPart[i] = sequence[i].Real;
            }
            return realPart;
        }
        
        /// <summary>
        /// Extract the imaginary part of a complex array
        /// </summary>
        /// <param name="sequence">= Complex array</param>
        /// <returns>= Imaginary part of complex array</returns>
        private double[] ImaginaryPart(Complex[] sequence)
        {
            double[] imaginaryPart = new double[sequence.Length];
            for (int i=0; i<sequence.Length; i++)
            {
                imaginaryPart[i] = sequence[i].Imaginary;
            }
            return imaginaryPart;
        }

        [Fact]
        public void Test_MinBlockLength()
        {
            // Arrange
            double[] h = { 1, 1, 1 };

            // Act
            int hBlockLength = this.Signal.MinBlockLength(h.Length);

            // Assert
            Assert.Equal(8, hBlockLength);
        }

        [Fact]
        public void Test_AcquireSignal()
        {
            // Arrange
            double[] h = { 1, 1, 1 };
            double[] x = { 3, -1, 0, 1, 3, 2, 0, 1, 2, 1 };
            double[] expected = { 3, 2, 2, 0, 4, 6, 5, 3, 3, 4, 3, 1 };
            int M = h.Length;
            int N = this.Signal.MinBlockLength(h.Length);
            int L = N - M + 1;
            Complex[] appendSequence = new Complex[L - 1];
            Array.Clear(appendSequence, 0, L - 1);
            // FIR filter frequency response
            List<Complex> hComplex = this.ToComplex(h).ToList();
            Complex[] hFFT = Signal.ComputePaddedFFT(ref hComplex, N, appendSequence.ToList());
            double phase = 0;

            // Act
            foreach (double dataPoint in x)
            {
                Signal.AcquireSignal(dataPoint, phase, 0, 0, L, ref hFFT);
            }
            Signal.LastSamplesConvolution(ref hFFT, L);
            double[] xhConvRe = this.RealPart(Signal.baseBandSignal.ToArray());
            double[] xhConvImg = this.ImaginaryPart(Signal.baseBandSignal.ToArray());

            // Assert
            for (int i=0; i<xhConvRe.Length; i++)
            {
                Assert.True(Math.Abs(xhConvRe[i] - expected[i]) < 1e-12);
                Assert.True(Math.Abs(xhConvImg[i]) < 1e-12);
            }
        }

        [Fact]
        public void Test_PeakDetectionAlgorithm()
        {
            // Arrange
            var input = new List<double> {1.0, 1.0, 1.1, 1.0, 0.9, 1.0, 1.0, 1.1, 1.0, 0.9, 1.0,
                1.1, 1.0, 1.0, 0.9, 1.0, 1.0, 1.1, 1.0, 1.0, 1.0, 1.0, 1.1, 0.9, 1.0, 1.1, 1.0, 1.0, 0.9,
                1.0, 1.1, 1.0, 1.0, 1.1, 1.0, 0.8, 0.9, 1.0, 1.2, 0.9, 1.0, 1.0, 1.1, 1.2, 1.0, 1.5, 1.0,
                3.0, 2.0, 5.0, 3.0, 2.0, 1.0, 1.0, 1.0, 0.9, 1.0, 1.0, 3.0, 2.6, 4.0, 3.0, 3.2, 2.0, 1.0,
                1.0, 0.8, 4.0, 4.0, 2.0, 2.5, 1.0, 1.0, 1.0 };
            var timeStamps = new List<double>(new double[input.Count]);
            var expectedPeaks = new List<double> { 45, 49, 60, 68 };
            List<Complex> inputSignal = ToComplex(input.ToArray()).ToList();
            for (int i=0; i<timeStamps.Count; i++)
            {
                timeStamps[i] = i;
            }
            int lag = 30;
            double threshold = 5;
            double influence = 0;

            // Act
            PeakDetectionOutput detector = PeakDetection.PeakDetectionAlgorithm(inputSignal, timeStamps, lag, threshold, influence);

            // Assert
            Assert.Equal(expectedPeaks.Count, detector.peakIndexes.Count);
            for (int i=0; i<expectedPeaks.Count; i++)
            {
                Assert.Equal(expectedPeaks.ElementAt(i), detector.peakIndexes.ElementAt(i));
            }
        }
    }
}
