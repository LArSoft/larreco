////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       Graph
// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
//
// Iterface to run Tensorflow graph saved to a file. First attempts, not functional yet.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "tf_graph.h"

#include "tensorflow/core/public/session.h"
#include "tensorflow/core/platform/env.h"

// -------------------------------------------------------------------
tf::Graph::Graph(const char* graph_file_name, bool & success)
{
    success = false; // until all is done correctly

    //tensorflow::Scope root = tensorflow::Scope::NewRootScope();
    //fSession = new tensorflow::ClientSession(root);

    auto status = tensorflow::NewSession(tensorflow::SessionOptions(), &fSession);
    if (!status.ok())
    {
        std::cout << status.ToString() << std::endl;
        return;
    }

    tensorflow::GraphDef graph_def;
    status = tensorflow::ReadBinaryProto(tensorflow::Env::Default(), graph_file_name, &graph_def);
    if (!status.ok())
    {
        std::cout << status.ToString() << std::endl;
        return;
    }

    size_t ng = graph_def.node().size();
    fInputName = graph_def.node()[0].name();
    fOutputName = graph_def.node()[ng - 1].name();

    status = fSession->Create(graph_def);
    if (!status.ok())
    {
        std::cout << status.ToString() << std::endl;
        return;
    }

    success = true; // ok, graph oaded from the file
}

tf::Graph::~Graph()
{
    fSession->Close();
    delete fSession;
}
// -------------------------------------------------------------------

std::vector<float> tf::Graph::run(const std::vector< std::vector<float> > & x)
{
    if (x.empty() || x.front().empty()) return std::vector<float>();
    //size_t rows = x.size(), cols = x.front().size();
/*
    tensorflow::Tensor _x(tensorflow::DT_FLOAT, tensorflow::TensorShape({ rows, cols }));

    //_x.scalar<float>()() = x;

    auto input_map = _x.tensor<float, 2>();
    for (size_t r = 0; r < rows; ++r) {
        const auto & row = x[r];
        for (size_t c = 0; c < cols; ++c) { input_map(r, c) = row[c]; }
    }

    std::vector< std::pair<std::string, tensorflow::Tensor> > inputs = {
        { fInputName, _x }
    };

    std::cout << "run session" << std::endl;

    std::vector<tensorflow::Tensor> outputs;
    auto status = fSession->Run(inputs, { fOutputName }, {}, &outputs);

    if (status.ok())
    {
        std::cout << "get output" << std::endl;
        auto output_map = outputs.front().tensor<float, 2>();

        std::cout << "shape:";
        size_t nd = outputs.front().dims();
        for (size_t d = 0; d < nd; ++d) { std::cout << " " << outputs.front().dim_size(d); }
        std::cout << std::endl;

        std::vector< float > output(outputs.front().dim_size(nd-1));
        for (size_t i = 0; i < output.size(); ++i) { output[i] = output_map(0, i); }
        return output;
    }
    else
    {
        std::cout << status.ToString() << std::endl;
        return std::vector<float>();
    }
*/

    return std::vector<float>();

    // (There are similar methods for vectors and matrices here:
    // https://github.com/tensorflow/tensorflow/blob/master/tensorflow/core/public/tensor.h)
}
// -------------------------------------------------------------------


