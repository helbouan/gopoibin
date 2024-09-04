package gopoibin

import (
	"math"

	"github.com/mjibson/go-dsp/fft"
)

// PoissonBinomialPMFSDFT calculates the distribution using the DFT method
func PoissonBinomialPMFSDFT(ps []float64) []float64 {
	psNum := len(ps)
	omega := 2 * math.Pi / float64(psNum+1)
	a := make([]float64, psNum+1)
	b := make([]float64, psNum+1)
	for i := 0; i < psNum+1; i++ {
		a[i] = 1.0
		b[i] = 0.0
	}
	// a[0] = 1.0
	for l := 1; l <= int(math.Ceil(float64(psNum)/2)); l++ {
		zlReal := make([]float64, psNum)
		zlImag := make([]float64, psNum)

		for i := 0; i < psNum; i++ {
			zlReal[i] = 1 - ps[i] + ps[i]*math.Cos(omega*float64(l))
			zlImag[i] = ps[i] * math.Sin(omega*float64(l))
		}

		var dl float64
		for i := 0; i < psNum; i++ {
			dl += math.Log(math.Sqrt(zlReal[i]*zlReal[i] + zlImag[i]*zlImag[i]))
		}
		dl = math.Exp(dl)

		var sumArgs float64
		for i := 0; i < psNum; i++ {
			sumArgs += math.Atan2(zlImag[i], zlReal[i])
		}

		a[l] = dl * math.Cos(sumArgs)
		b[l] = dl * math.Sin(sumArgs)
	}

	for l := int(math.Ceil(float64(psNum/2))) + 1; l <= psNum; l++ {
		a[l] = a[psNum+1-l]
		b[l] = -b[psNum+1-l]
	}

	// Combine real and imaginary parts
	x := make([]complex128, psNum+1)
	for i := 0; i <= psNum; i++ {
		x[i] = complex(a[i], b[i]) / complex(float64(psNum+1), 0)
	}

	// Perform FFT
	fftResult := fft.FFT(x)

	// Extract real parts
	result := make([]float64, psNum+1)
	for i := 0; i <= psNum; i++ {
		result[i] = real(fftResult[i])
	}

	return result
}
