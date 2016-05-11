
#include <vector>
#include <string>

#include <iostream>
#include <cstdlib>
#include <math.h>

// ROOT vector, matrix and XML libs
#include <TVectorT.h>
#include <TMatrixT.h>
#include <TDOMParser.h>
#include <TXMLNode.h>
#include <TXMLAttr.h>
#include <TList.h>

#ifndef MLP_NNReader_H
#define MLP_NNReader_H

#define F_TABLE_MAX 12.0F
#define F_TABLE_LENGTH 256

namespace nnet
{
	class NNReader;
}


/// Class for loading and running neural network stored in xml file.
class nnet::NNReader
{
public:
	enum ELayerType { kLayerBase = 0, kLayerInp = 1, kLayerLin = 2, kLayerSig = 3, kLayerTanh = 4, kLayerSoftmax = 5 };

	NNReader(const char* xmlFileName);
	virtual ~NNReader(void);

	unsigned int GetInputLength(void) const { return _nInps; }
	int GetOutputLength(void) const { return _nOuts; }
	ELayerType GetOutputType(void) { return _layers.back()->GetType(); }
	unsigned int NetworkSize(void) const;

	void Run(const TVectorT<float>& vInput);
	void Run(const std::vector<float>& vInput);
	void Run(const std::vector<double>& vInput);

	float GetOneOutput(int neuronIndex) const;
	TVectorT<float> GetAllOutputs() const;
	TVectorT<float> GetLinOutputs() const;

private:
	float OutputSigma2(unsigned int outIndex);
	TMatrixT<float> OutputCov(void);

	void HiddenGradient(unsigned int layerIndex);
	void ClearGradient(void);

	void F_dLin_dw(TMatrixT<float>& vDst, size_t row, unsigned int outIndex);
	void GradientToVector(TMatrixT<float>& vDst, size_t row) const;

	void F_dLin_dw(TVectorT<float>& vDst, unsigned int outIndex);
	void GradientToVector(TVectorT<float>& vDst) const;

	class Layer
	{
	public:
		virtual ~Layer();

		virtual float ActFunction(float xValue) const = 0;
		virtual float ActDerivative(unsigned int nIndex) const = 0;
		virtual TVectorT<float> ActGradient(unsigned int nIndex) const;

		virtual void Activate() = 0;
		TVectorT<float>& GetOutput() { return _output; }

		NNReader::ELayerType GetType(void) { return _type; }

		TVectorT<float>& Weights(unsigned int neuronIndex) { return _coeffs[neuronIndex]; }
		float& Weight(int neuronIndex, int weightIndex) { return _coeffs[neuronIndex][weightIndex]; }

		unsigned int NInputs(void) { return _nInputs; }
		unsigned int NLength(void) { return _nLength; }

		float& Delta(unsigned int neuronIndex) { return _deltas[neuronIndex]; }
		TVectorT<float>& Gradient(unsigned int neuronIndex) { return _gradient[neuronIndex]; }
		void ClearGradient(void);

	protected:
		Layer(int nNeurons, Layer* prevLayer);
		Layer();

		NNReader::ELayerType _type;

		int _nLength, _nInputs;
		TVectorT<float> _output;
		TVectorT<float>* _coeffs;
		TVectorT<float>* _gradient;
		TVectorT<float> _deltas;
		Layer* _previous;
	};

	class LayerSig : public Layer
	{
	public:
		LayerSig(int nNeurons, Layer* prevLayer);

		float ActFunction(float xValue) const { return _fsig(xValue); }
		float ActDerivative(unsigned int nIndex) const;

		void Activate();

	private:
		float _fsig(double dotprod) const;

		static float f_factor;
		static float f_f[F_TABLE_LENGTH], f_d[F_TABLE_LENGTH];
	};

	class LayerTanh : public Layer
	{
	public:
		LayerTanh(int nNeurons, Layer* prevLayer);

		float ActFunction(float xValue) const { return _ftanh(xValue); }
		float ActDerivative(unsigned int nIndex) const;

		void Activate();

	private:
		float _ftanh(double dotprod) const;

		static float f_factor;
		static float f_f[F_TABLE_LENGTH], f_d[F_TABLE_LENGTH];
	};

	class LayerSoftmax : public Layer
	{
	public:
		LayerSoftmax(int nNeurons, Layer* prevLayer);

		float ActFunction(float xValue) const;
		float ActDerivative(unsigned int nIndex) const { return 1.0F; }
		TVectorT<float> ActGradient(unsigned int nIndex) const;

		void Activate();
		void ModerateOutputs(const TMatrixT<float>& S, unsigned int mcSize = 256);
		void ModerateOutputs(const TMatrixT<float>& S, const TMatrixT<float>& N01);

		TVectorT<float>& GetLinOutput(void) { return _linout; }

	protected:
		TVectorT<float> _linout;
	};

	class LayerLin : public Layer
	{
	public:
		LayerLin(int nNeurons, Layer* prevLayer) : Layer(nNeurons, prevLayer) { }

		float ActFunction(float xValue) const { return xValue; }
		float ActDerivative(unsigned int nIndex) const { return 1.0F; }

		void Activate();
	};

	class LayerInp : public Layer
	{
	public:
		LayerInp(int nInps);
		LayerInp(int nInps, const TVectorT<float>& vMean, const TVectorT<float>& vStd);
		virtual ~LayerInp();

		float ActFunction(float xValue) const { return xValue; }
		float ActDerivative(unsigned int nIndex) const { return 1.0F; }

		void Activate();
		void Activate(const TVectorT<float>& vInput);
		void Activate(const std::vector<float>& vInput);
		void Activate(const std::vector<double>& vInput);

	private:
		void Normalize();
		bool _normalize;
		TVectorT<float> *_mean, *_istdev;
	};


	unsigned int _nInps; 
	int _nOuts;
	std::vector<Layer*> _layers;
	TMatrixT<float>* _invH;
	TMatrixT<float>* _N01;

	void SetupNetwork(const char* xmlfile);
};
#endif
