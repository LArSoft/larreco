#include "larreco/RecoAlg/ImagePatternAlgs/MLP/NNReader.h"

nnet::NNReader::NNReader(const char* xmlFileName) :
	_invH(0), _N01(0)
{
	try { SetupNetwork(xmlFileName); }
	catch (char const* e) { std::cout << e << std::endl; throw "NNetNotCreated"; }
}

void nnet::NNReader::SetupNetwork(const char* xmlFile)
{
	std::string type, node_name, attr_name;
	int size = 0; 

	TDOMParser *xml_parser = new TDOMParser();  // creates TDOMParser root object
	TXMLNode *node; 
	TXMLAttr *attr;
	TList *attr_list;

	xml_parser->SetValidate(0);  // switch validation off (no need for dtd file)
	int err_code = xml_parser->ParseFile(xmlFile);  // parses xmlFile and creates TXMLDocument root object
	if (err_code) throw "ErrXmlParser";
	node = xml_parser->GetXMLDocument()->GetRootNode()->GetChildren();  // node = pointer to the first child of RootNode (<model> in our case)

	std::string lay("layer");
	std::string typ("type");	
	std::string siz("size");
	std::string inp("Input");
	std::string sig("Sigmoid");
	std::string tnh("Tanh");
	std::string smx("Softmax");
	std::string lin("Linear");
	std::string vec("vector");
	std::string mat("matrix");
	std::string mxr("rows");
	std::string mxc("cols");
	std::string row("row");
	std::string rid("index");
	std::string hes("ihessian");
	std::string bay("bayes");

	bool no_layer = true;
	for (; node; node = node->GetNextNode())  // searches all children of <model>
	{   
		node_name = node->GetNodeName();
		if (node_name.compare(lay) == 0)
		{
			no_layer = false;
			bool no_size = true;
			bool no_type = true;

			if (node->HasAttributes())										
			{	
				attr_list = node->GetAttributes();
				for (attr = (TXMLAttr*)attr_list->First (); attr; attr = (TXMLAttr*)attr_list->After(attr))  // searches through all <layer> attributes
				{
					attr_name = attr->Key();
					if (attr_name.compare(typ) == 0) { type = attr->GetValue(); no_type = false; }
					if (attr_name.compare(siz) == 0) { size = atoi(attr->GetValue()); no_size= false; }
				}
			}
			if (no_size) { throw "ErrLayerSizeNotSpecified"; }
			if (no_type) { throw "ErrLayerTypeNotSpecified"; }

			if (type.compare(inp) == 0) 
			{ 
				if (_layers.size() == 0) { _layers.push_back(new LayerInp(size)); _nInps = size; }
			}
			else
			{
				if (type.compare(sig) == 0) _layers.push_back(new LayerSig(size, _layers.back()));
				else if (type.compare(tnh) == 0) _layers.push_back(new LayerTanh(size, _layers.back()));
				else if (type.compare(smx) == 0) _layers.push_back(new LayerSoftmax(size, _layers.back()));
				else if (type.compare(lin) == 0) _layers.push_back(new LayerLin(size, _layers.back()));
				else throw "ErrLayeTypeNotSupported";
				_nOuts = size;

				int i = 0;
				for (TXMLNode* vec_node = node->GetChildren(); vec_node; vec_node = vec_node->GetNextNode() )	// searches through all <layer> children
				{
					node_name = vec_node->GetNodeName();
					if (node_name.compare(vec) == 0)
					{
						std::string vec_coef;
						int vec_size = 0;
						if (vec_node->HasAttributes())
						{
							attr_list = vec_node->GetAttributes(); // creats root object TList of root objects TXMLAttr
							attr = (TXMLAttr*)attr_list->First();  // we assume the one and only attribute of <vector> is 'length'
							vec_size = atoi(attr->GetValue());
						}
						else { throw "ErrVectorSizeNotSpecified"; }
						if (vec_node -> HasChildren()) { vec_coef = vec_node->GetChildren()->GetContent(); } //gets the numbers specifying <vector> 	
						else { throw "ErrVectorCoefNotSpecified"; }

						int c, start;
						c = 0; start = 0;
						while (start != -1)	// changes string 'vec_coef' to floats. we assume numbers in vector are separated by ' ' 
						{				// this needs to be improved a bit, in fact numbers are separated by ';'
							(*_layers.back()).Weight(i, c) = (float)atof(vec_coef.substr(start, vec_coef.size() - start).c_str());
							start = (int)vec_coef.find_first_of(' ', start + 1);
							c++;
						}
						if (c != vec_size) { throw "ErrVectorSize"; } // checks if the number of above floats meets the <vector> length 
						i++;
					}
				}
				if (i!= size) { throw "ErrLayerSize"; } // checks if the number of <vector>s meets the size of the layer	
			}
		}
		else if (node_name.compare(bay) == 0)
		{
			for (TXMLNode* bay_node = node->GetChildren(); bay_node; bay_node = bay_node->GetNextNode()) // searches through all <bayes> children
			{
				node_name = bay_node->GetNodeName();
				if (node_name.compare(hes) == 0)
					for (TXMLNode* hes_node = bay_node->GetChildren(); hes_node; hes_node = hes_node->GetNextNode())
					{
						node_name = hes_node->GetNodeName();
						if ((node_name.compare(mat) == 0) && hes_node->HasAttributes())
						{
							unsigned int hrows = 0, hcols = 0;
							attr_list = hes_node->GetAttributes();
							for (attr = (TXMLAttr*)attr_list->First(); attr; attr = (TXMLAttr*)attr_list->After(attr))
							{
								attr_name = attr->Key();
								if (attr_name.compare(mxr) == 0) hrows = atoi(attr->GetValue());
								if (attr_name.compare(mxc) == 0) hcols = atoi(attr->GetValue());
							}
							if (hrows && hcols)
							{
								_invH = new TMatrixT<float>(hrows, hcols);
								for (TXMLNode* mat_node = hes_node->GetChildren(); mat_node; mat_node = mat_node->GetNextNode())
								{
									node_name = mat_node->GetNodeName();
									if ((node_name.compare(row) == 0) && (mat_node->HasChildren()))
									{
										attr_list = mat_node->GetAttributes();
										attr = (TXMLAttr*)attr_list->First();
										unsigned int rindex = atoi(attr->GetValue());

										std::string row_data = mat_node->GetChildren()->GetContent();

										int c, start;
										c = 0; start = 0;
										while ((start != -1) && (c < _invH->GetNcols()))
										{
											(*_invH)(rindex, c) = (float)atof(row_data.substr(start, row_data.size() - start).c_str());
											start = (int)row_data.find_first_of(' ', start + 1);
											c++;
										}
										if (c != _invH->GetNcols() ) { throw "ErrVectorSize"; }
									}
								}
							}
						}
					}
			}
		}
	}
	if (no_layer) { throw "ErrNoLayerInNetwork"; }

	if (_invH)
	{
		_N01 = new TMatrixT<float>(_nOuts, 1000);
		(*_N01) = 0.0F; // <------------------------ implement rnd values
		// _N01->RandomNormal();
	}
}

nnet::NNReader::~NNReader()
{
	while (_layers.size() > 0)
	{
		delete _layers.back();
		_layers.pop_back();
	}
	if (_invH) delete _invH;
	if (_N01) delete _N01;
}

unsigned int nnet::NNReader::NetworkSize(void) const
{
	unsigned int netSize = 0;
	for (size_t l = 0; l < _layers.size() - 1; l++)
		netSize += (_layers[l]->NLength() + 1) * _layers[l + 1]->NLength();
	return netSize;
}

void nnet::NNReader::Run(const TVectorT<float>& vInput)
{
	if ((size_t)vInput.GetNrows() == _nInps)
	{
		((LayerInp*)_layers[0])->Activate(vInput);
		for (unsigned int i = 1; i < _layers.size(); i++) _layers[i]->Activate();
		if (_invH && _N01)
		{
			LayerSoftmax* softmax = dynamic_cast< LayerSoftmax* >(_layers.back());
			if (softmax) softmax->ModerateOutputs(OutputCov(), *_N01);
		}
	}
	else throw "InpSizeErr";
}

void nnet::NNReader::Run(const std::vector<float>& vInput)
{
	if (vInput.size() == _nInps)
	{
		((LayerInp*)_layers[0])->Activate(vInput);
		for (unsigned int i = 1; i < _layers.size(); i++) _layers[i]->Activate();
		if (_invH && _N01)
		{
			LayerSoftmax* softmax = dynamic_cast< LayerSoftmax* >(_layers.back());
			if (softmax) softmax->ModerateOutputs(OutputCov(), *_N01);
		}
	}
	else throw "InpSizeErr";
}

void nnet::NNReader::Run(const std::vector<double>& vInput)
{
	if (vInput.size() == _nInps)
	{
		((LayerInp*)_layers[0])->Activate(vInput);
		for (unsigned int i = 1; i < _layers.size(); i++) _layers[i]->Activate();
		if (_invH && _N01)
		{
			LayerSoftmax* softmax = dynamic_cast< LayerSoftmax* >(_layers.back());
			if (softmax) softmax->ModerateOutputs(OutputCov(), *_N01);
		}
	}
	else throw "InpSizeErr";
}

float nnet::NNReader::GetOneOutput(int neuronIndex) const
{
	if (neuronIndex < _nOuts) return _layers.back()->GetOutput()[neuronIndex];
	else throw "OutSizeErr";
}

TVectorT<float> nnet::NNReader::GetAllOutputs() const
{
	TVectorT<float> vret(_nOuts);
	for (size_t i = 0; i < (size_t)_nOuts; ++i) vret[i] = _layers.back()->GetOutput()[i];
	//vret.Set(_layers.back()->GetOutput(), _nOuts);
	return vret;
}

TVectorT<float> nnet::NNReader::GetLinOutputs() const
{
	TVectorT<float> vret(_nOuts);
	LayerSoftmax* softmax = dynamic_cast< LayerSoftmax* >(_layers.back());
	if (softmax) vret = softmax->GetLinOutput();
	else
	{
		LayerLin* linear = dynamic_cast< LayerLin* >(_layers.back());
		if (linear)
		{
			for (size_t i = 0; i < (size_t)_nOuts; ++i) vret[i] = linear->GetOutput()[i];
			//vret.Set(linear->GetOutput(), _nOuts);
		}
	}
	return vret;
}

float nnet::NNReader::OutputSigma2(unsigned int outIndex)
{
	if (!_invH) return 0.0F;

	if (_nOuts == 1)
	{
		TVectorT<float> dLin_dw(NetworkSize());
		dLin_dw = 0.0F;

		F_dLin_dw(dLin_dw, outIndex);
		return Mult(dLin_dw, *_invH, dLin_dw);
	}
	else
	{
		TMatrixT<float> G(_nOuts, NetworkSize());
		G = 0.0F;

		for (int n = 0; n < _nOuts; n++)
		{
			F_dLin_dw(G, n, n);
			//F_dLin_dw(G.Row(n), n);
		}
		return (G * (*_invH) * G.T())(outIndex, outIndex);
	}
}

TMatrixT<float> nnet::NNReader::OutputCov(void)
{
	if (!_invH) return TMatrixT<float>(0, 0);

	if (_nOuts == 1)
	{
		TMatrixT<float> sigma2(1, 1);
		sigma2(0, 0) = OutputSigma2(0);
		return sigma2;
	}

	TMatrixT<float> G(_nOuts, NetworkSize());
	for (int n = 0; n < _nOuts; n++)
	{
		F_dLin_dw(G, n, n);
		//F_dLin_dw(G.Row(n), n);
	}

	return G * (*_invH) * G.T();
}

void nnet::NNReader::HiddenGradient(unsigned int layerIndex)
{
	unsigned int nThis = _layers[layerIndex]->NLength();
	unsigned int nNext = _layers[layerIndex + 1]->NLength();

	float delta;
	TVectorT<float> dsum(nThis + 1);
	dsum = 0.0F;

	for (unsigned int n = 0; n < nNext; n++) dsum += _layers[layerIndex + 1]->Delta(n) * _layers[layerIndex + 1]->Weights(n);
		// dsum.ScaleAdd(_layers[layerIndex + 1]->Weights(n), _layers[layerIndex + 1]->Delta(n));

	for (unsigned int n = 0; n < nThis; n++)
	{
		delta = dsum[n] * _layers[layerIndex]->ActDerivative(n);
		_layers[layerIndex]->Delta(n) = delta;
		_layers[layerIndex]->Gradient(n) += delta * _layers[layerIndex - 1]->GetOutput();
		//_layers[layerIndex]->Gradient(n).ScaleAdd(_layers[layerIndex - 1]->GetOutput(), delta);
	}
}

void nnet::NNReader::ClearGradient(void)
{
	for (unsigned int l = 1; l < _layers.size(); l++) _layers[l]->ClearGradient();
}

void nnet::NNReader::F_dLin_dw(TMatrixT<float>& vDst, size_t row, unsigned int outIndex)
{
	ClearGradient();

	//     last hidden layer's output       copy to  output layer's gradient of outIndex'th neuron
	for (int i = 0; i < _layers[_layers.size() - 2]->GetOutput().GetNrows(); ++i)
		_layers[_layers.size() - 1]->Gradient(outIndex)[i] = _layers[_layers.size() - 2]->GetOutput()[i];
	//_layers[_layers.size() - 2]->GetOutput().CopyTo(_layers[_layers.size() - 1]->Gradient(outIndex));
	_layers[_layers.size() - 1]->Delta(outIndex) = 1.0F;

	unsigned int l = _layers.size() - 1;
	while (l-- > 1) HiddenGradient(l);

	GradientToVector(vDst, row);
}

void nnet::NNReader::GradientToVector(TMatrixT<float>& vDst, size_t row) const
{
	for (size_t l = 1, pos = 0; l < _layers.size(); l++)
		for (size_t n = 0; n < _layers[l]->NLength(); n++)
			for (int c = 0; c < _layers[l]->Gradient(n).GetNrows(); c++)
			{
				vDst(row, pos) = _layers[l]->Gradient(n)[c];
				pos++;
			}
}

void nnet::NNReader::F_dLin_dw(TVectorT<float>& vDst, unsigned int outIndex)
{
	ClearGradient();

	//     last hidden layer's output       copy to  output layer's gradient of outIndex'th neuron
	for (int i = 0; i < _layers[_layers.size() - 2]->GetOutput().GetNrows(); ++i)
		_layers[_layers.size() - 1]->Gradient(outIndex)[i] = _layers[_layers.size() - 2]->GetOutput()[i];
	//_layers[_layers.size() - 2]->GetOutput().CopyTo(_layers[_layers.size() - 1]->Gradient(outIndex));
	_layers[_layers.size() - 1]->Delta(outIndex) = 1.0F;

	unsigned int l = _layers.size() - 1;
	while (l-- > 1) HiddenGradient(l);

	GradientToVector(vDst);
}

void nnet::NNReader::GradientToVector(TVectorT<float>& vDst) const
{
	for (size_t l = 1, pos = 0; l < _layers.size(); l++)
		for (size_t n = 0; n < _layers[l]->NLength(); n++)
			for (int c = 0; c < _layers[l]->Gradient(n).GetNrows(); c++)
			{
				vDst[pos] = _layers[l]->Gradient(n)[c];
				pos++;
			}
}


///   Layer Class (nested in NNReader)   ///
nnet::NNReader::Layer::Layer(int nNeurons, Layer* prevLayer) :
	_type(nnet::NNReader::kLayerBase),
	_nLength(nNeurons),
	_nInputs(prevLayer->_nLength),
	_output(nNeurons + 1),
	_deltas(nNeurons),
	_previous(prevLayer)
{
	_coeffs = new TVectorT<float>[nNeurons];
	_gradient = new TVectorT<float>[nNeurons];
	for (int i = 0; i < nNeurons; i++)
	{
		_coeffs[i].ResizeTo(_nInputs + 1);
		_coeffs[i] = 0.0F;

		_gradient[i].ResizeTo(_nInputs + 1);
		_gradient[i] = 0.0F;
	}
	_output[nNeurons] = 1.0F;
}

nnet::NNReader::Layer::Layer() :
	_nLength(0),
	_nInputs(0),
	_output(1),
	_coeffs(NULL),
	_gradient(NULL),
	_deltas(0),
	_previous(NULL)
{
}

nnet::NNReader::Layer::~Layer()
{
	if (_coeffs) delete [] _coeffs;
	if (_gradient) delete [] _gradient;
}

TVectorT<float> nnet::NNReader::Layer::ActGradient(unsigned int nIndex) const
{
	TVectorT<float> gradient(_nLength);
	gradient[nIndex] = ActDerivative(nIndex);
	return gradient;
}

void nnet::NNReader::Layer::ClearGradient(void)
{
	for (int n = 0; n < _nLength; n++) _gradient[n] = 0.0F;
	_deltas = 0.0F;
}


///   LayerSig Class (nested in NNReader)   ///
float nnet::NNReader::LayerSig::f_f[F_TABLE_LENGTH];
float nnet::NNReader::LayerSig::f_d[F_TABLE_LENGTH];
float nnet::NNReader::LayerSig::f_factor;

nnet::NNReader::LayerSig::LayerSig(int nNeurons, Layer* prevLayer) :
	Layer(nNeurons, prevLayer)
{
	_type = nnet::NNReader::kLayerSig;
	f_factor = (float)(F_TABLE_LENGTH - 1) / F_TABLE_MAX;
	for (unsigned int i = 0; i < F_TABLE_LENGTH; i++)
	{
		f_f[i] = (float)(1.0 / ( 1.0 + exp(-((double)i) / f_factor)));
		if (i > 0) f_d[i-1] = f_f[i] - f_f[i-1];
	}
}

float nnet::NNReader::LayerSig::ActDerivative(unsigned int nIndex) const
{
	float fValue = _output[nIndex];
	return fValue * (1.0F - fValue);
}

float nnet::NNReader::LayerSig::_fsig(double dotprod) const
{
	int i;
	float xd;

	if (dotprod >= 0.0)
	{
		xd = (float)(dotprod * f_factor);

		i = (int)xd;
		if (i < (F_TABLE_LENGTH - 1))
			return f_f[i] + f_d[i] * (xd - i);
		return f_f[F_TABLE_LENGTH - 1];
	}
	else
	{
		xd = (float)(-dotprod * f_factor);

		i = (int)xd;
		if (i < (F_TABLE_LENGTH - 1))
			return 1.0F - (f_f[i] + f_d[i] * (xd - i));
		return 1.0F - f_f[F_TABLE_LENGTH - 1];
	}
}

void nnet::NNReader::LayerSig::Activate()
{
	for (int i = 0; i < _nLength; i++)
		_output[i] = _fsig(_previous->GetOutput() * _coeffs[i]);
}

///   LayerTanh Class (nested in NNReader)   ///
float nnet::NNReader::LayerTanh::f_f[F_TABLE_LENGTH];
float nnet::NNReader::LayerTanh::f_d[F_TABLE_LENGTH];
float nnet::NNReader::LayerTanh::f_factor;

nnet::NNReader::LayerTanh::LayerTanh(int nNeurons, Layer* prevLayer) :
	Layer(nNeurons, prevLayer)
{
	_type = nnet::NNReader::kLayerTanh;
	f_factor = (float)(F_TABLE_LENGTH - 1) / F_TABLE_MAX;
	for (unsigned int i = 0; i < F_TABLE_LENGTH; i++)
	{
		f_f[i] = (float)tanh(((double)i) / f_factor);
		if (i > 0) f_d[i-1] = f_f[i] - f_f[i-1];
	}
}

float nnet::NNReader::LayerTanh::ActDerivative(unsigned int nIndex) const
{
	float fValue = _output[nIndex];
	return 1.0F - fValue * fValue;
}

float nnet::NNReader::LayerTanh::_ftanh(double dotprod) const
{
	int i;
	float xd;

	if (dotprod >= 0.0)
	{
		xd = (float)(dotprod * f_factor);

		i = (int)xd;
		if (i < (F_TABLE_LENGTH - 1))
			return f_f[i] + f_d[i] * (xd - i);
		return f_f[F_TABLE_LENGTH - 1];
	}
	else
	{
		xd = (float)(-dotprod * f_factor);

		i = (int)xd;
		if (i < (F_TABLE_LENGTH - 1))
			return -(f_f[i] + f_d[i] * (xd - i));
		return -f_f[F_TABLE_LENGTH - 1];
	}
}

void nnet::NNReader::LayerTanh::Activate()
{
	for (int i = 0; i < _nLength; i++)
		_output[i] = _ftanh(_previous->GetOutput() * _coeffs[i]);
}

/// LayerSoftmax Class (nested in NNReader) ///
nnet::NNReader::LayerSoftmax::LayerSoftmax(int nNeurons, Layer* prevLayer) :
	Layer(nNeurons, prevLayer),
	_linout(nNeurons)
{
	_type = nnet::NNReader::kLayerSoftmax;
}

float nnet::NNReader::LayerSoftmax::ActFunction(float xValue) const
{
	if (xValue < 30.0F) return (float)exp(xValue);
	else return (float)exp(30.0);
}

TVectorT<float> nnet::NNReader::LayerSoftmax::ActGradient(unsigned int nIndex) const
{
	float fValue = _output[nIndex];
	TVectorT<float> gradient(_nLength);
	for (int i = 0; i < gradient.GetNrows(); i++)
	{
		if (i == (int)nIndex) gradient[i] = fValue - fValue * fValue;
		else gradient[i] = -fValue * _output[i];
	}
	return gradient;
}

void nnet::NNReader::LayerSoftmax::Activate()
{
	double o, sum = 0.0;
	for (int i = 0; i < _nLength; i++)
	{
		o = _previous->GetOutput() * _coeffs[i];
		_linout[i] = o;

		if (o > 30.0) o = 30.0;

		o = exp(o); sum += o;
		_output[i] = (float)o;
	}

	if (sum > 0.0) _output *= (float)(1.0 / sum);
	else _output = 1.0F / _nLength;
	_output[_nLength] = 1.0F;
}

void nnet::NNReader::LayerSoftmax::ModerateOutputs(const TMatrixT<float>& S, unsigned int mcSize)
{
	TMatrixT<float> mcOutputs(_nLength, mcSize);
	mcOutputs = 0.0F; // <--------------------- implement rnd values
	//mcOutputs.RandomNormal(_linout, S);

	TVectorT<float> msum(_nLength);
	TVectorT<float> mout(_nLength);
	for (unsigned int i = 0; i < mcSize; i++)
	{
		for (int n = 0; n < _nLength; n++)
			mout[n] = ActFunction(mcOutputs(n, i));
		mout *= 1.0F / mout.Sum();
		msum += mout;
	}
	msum *= 1.0F / mcSize;
	for (int i = 0; i < _nLength; ++i) _output[i] = msum[i];
	_output[_nLength] = 1.0F;
}

void nnet::NNReader::LayerSoftmax::ModerateOutputs(const TMatrixT<float>& S, const TMatrixT<float>& N01)
{
	if (_nLength != N01.GetNrows()) throw "N01 rows must match the output length.";

	TMatrixT<float> mcOutputs(_nLength, N01.GetNcols());
	mcOutputs = 0.0F; // <--------------------- implement rnd values
	//mcOutputs.RandomNormal(_linout, S, N01);

	TVectorT<float> msum(_nLength);
	TVectorT<float> mout(_nLength);
	for (int i = 0; i < mcOutputs.GetNcols(); i++)
	{
		for (int n = 0; n < _nLength; n++)
			mout[n] = ActFunction(mcOutputs(n, i));
		mout *= 1.0F / mout.Sum();
		msum += mout;
	}
	msum *= 1.0F / mcOutputs.GetNcols();
	for (int i = 0; i < _nLength; ++i) _output[i] = msum[i];
	_output[_nLength] = 1.0F;
}

///   LayerLin Class (nested in NNReader)   ///
void nnet::NNReader::LayerLin::Activate()
{
	_type = nnet::NNReader::kLayerLin;
	for (int i = 0; i < _nLength; i++)
		_output[i] = (float)(_previous->GetOutput() * _coeffs[i]);
}

///   LayerInp Class (nested in NNReader)   ///
nnet::NNReader::LayerInp::LayerInp(int nInps) :
	_normalize(false)
{
	_type = nnet::NNReader::kLayerInp;
	_nLength = nInps;
	_nInputs = nInps;
	_output.ResizeTo(nInps + 1);
	_output = 0.0F; _output[nInps] = 1.0F;
}

nnet::NNReader::LayerInp::LayerInp(int nInps, const TVectorT<float>& vMean, const TVectorT<float>& vStd) :
	_normalize(true)
{
	if ((vMean.GetNrows() != nInps) || (vStd.GetNrows() != nInps)) throw "VectorSizeErr";

	_nLength = nInps;
	_nInputs = nInps;

	_output = TVectorT<float>(nInps + 1);
	_output = 0.0F; _output[nInps] = 1.0F;

	_mean = new TVectorT<float>(nInps + 1);
	for (int i = 0; i < nInps; i++) (*_mean)[i] = vMean[i];
	(*_mean)[nInps] = 0.0F;
	_istdev = new TVectorT<float>(nInps + 1);
	for (int i = 0; i < nInps; i++ ) (*_istdev)[i] = vStd[i];
	(*_istdev)[nInps] = 1.0F;
}

nnet::NNReader::LayerInp::~LayerInp()
{
	if (_normalize)
	{
		delete _mean;
		delete _istdev;
	}
}

void nnet::NNReader::LayerInp::Activate()
{
	_output = 0.0F;
}

void nnet::NNReader::LayerInp::Normalize()
{
	for (int i = 0; i < _output.GetNrows(); ++i)
	{
		_output[i] -= (*_mean)[i];
		_output[i] *= (*_istdev)[i];
	}
}

void nnet::NNReader::LayerInp::Activate(const TVectorT<float>& vInput)
{
	for (int i = 0; i < vInput.GetNrows(); ++i) _output[i] = vInput[i];
	if (_normalize) Normalize();
}

void nnet::NNReader::LayerInp::Activate(const std::vector<float>& vInput)
{
	for (unsigned int i = 0; i < vInput.size(); i++) _output[i] = vInput[i];
	if (_normalize) Normalize();
}

void nnet::NNReader::LayerInp::Activate(const std::vector<double>& vInput)
{
	for (unsigned int i = 0; i < vInput.size(); i++) _output[i] = (float)vInput[i];
	if (_normalize) Normalize();
}

