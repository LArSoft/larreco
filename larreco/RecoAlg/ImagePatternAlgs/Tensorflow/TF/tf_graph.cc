////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       Graph
// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
//              P.Plonski,                      from DUNE, WUT, Sept. 2017
//
// Iterface to run Tensorflow graph saved to a file. First attempts, quite functional.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "tf_graph.h"

#include "tensorflow/core/public/session.h"
#include "tensorflow/core/platform/env.h"

#include "tensorflow/core/public/session_options.h"

// -------------------------------------------------------------------
tf::Graph::Graph(const char* graph_file_name, const std::vector<std::string> & outputs, bool & success)
{
    success = false; // until all is done correctly

    // Force tf to only use a single core so it doesn't eat batch farms
    tensorflow::SessionOptions options;
    tensorflow::ConfigProto &config = options.config;
    config.set_inter_op_parallelism_threads(1);
    config.set_intra_op_parallelism_threads(1);
    config.set_use_per_session_threads(false);

    auto status = tensorflow::NewSession(options, &fSession);
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

    // last node as output if no specific name provided
    if (outputs.empty()) { fOutputNames.push_back(graph_def.node()[ng - 1].name()); }
    else // or last nodes with names containing provided strings
    {
        std::string last, current, basename, name;
        for (size_t n = 0; n < ng; ++n)
        {
            name = graph_def.node()[n].name();
            auto pos = name.find("/");
            if (pos != std::string::npos) { basename = name.substr(0, pos); }
            else { continue; }

            bool found = false;
            for (const auto & s : outputs)
            {
                if (name.find(s) != std::string::npos) { found = true; break; }
            }
            if (found)
            {
                if (!last.empty() && (basename != current))
                {
                    fOutputNames.push_back(last);
                }
                current = basename;
                last = name;
            }
        }
        if (!last.empty()) { fOutputNames.push_back(last); }
    }
    if (fOutputNames.empty())
    {
        std::cout << "Output nodes not found in the graph." << std::endl;
        return;
    }

    status = fSession->Create(graph_def);
    if (!status.ok())
    {
        std::cout << status.ToString() << std::endl;
        return;
    }

    success = true; // ok, graph loaded from the file
}

tf::Graph::~Graph()
{
    fSession->Close();
    delete fSession;
}
// -------------------------------------------------------------------

std::vector<float> tf::Graph::run(const std::vector< std::vector<float> > & x)
{
    if (x.empty() || x.front().empty()) { return std::vector<float>(); }

    long long int rows = x.size(), cols = x.front().size();

    tensorflow::Tensor _x(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, rows, cols, 1 }));
    auto input_map = _x.tensor<float, 4>();

    for (long long int r = 0; r < rows; ++r) {
        const auto & row = x[r];
        for (long long int c = 0; c < cols; ++c) {
            input_map(0, r, c, 0) = row[c];
        }
    }

    auto result = run(_x);
    if (!result.empty()) { return result.front(); }
    else { return std::vector<float>(); }
}
// -------------------------------------------------------------------

std::vector< std::vector<float> > tf::Graph::run(
	const std::vector<  std::vector<  std::vector< std::vector<float> > > > & x,
	long long int samples)
{
    if ((samples == 0) || x.empty() || x.front().empty() || x.front().front().empty() || x.front().front().front().empty())
        return std::vector< std::vector<float> >();

    if ((samples == -1) || (samples > (long long int)x.size())) { samples = x.size(); }

    long long int
              rows = x.front().size(),
              cols = x.front().front().size(),
              depth = x.front().front().front().size();

    tensorflow::Tensor _x(tensorflow::DT_FLOAT, tensorflow::TensorShape({ samples, rows, cols, depth }));
    auto input_map = _x.tensor<float, 4>();
    for (long long int s = 0; s < samples; ++s) {
        const auto & sample = x[s];
        for (long long int r = 0; r < rows; ++r) {
            const auto & row = sample[r];
            for (long long int c = 0; c < cols; ++c) {
                const auto & col = row[c];
                for (long long int d = 0; d < depth; ++d) {
                    input_map(s, r, c, d) = col[d];
                }
            }
        }
    }

    return run(_x);
}
// -------------------------------------------------------------------

std::vector< std::vector< float > > tf::Graph::run(const tensorflow::Tensor & x)
{
    std::vector< std::pair<std::string, tensorflow::Tensor> > inputs = {
        { fInputName, x }
    };

    //std::cout << "run session" << std::endl;

    std::vector<tensorflow::Tensor> outputs;
    auto status = fSession->Run(inputs, fOutputNames, {}, &outputs);

    //std::cout << "out size " << outputs.size() << std::endl;

    if (status.ok())
    {
        size_t samples = 0, nouts = 0;
        for (size_t o = 0; o < outputs.size(); ++o)
        {
            if (o == 0) { samples = outputs[o].dim_size(0); }
            else if ((int)samples != outputs[o].dim_size(0))
            {
                throw std::string("TF outputs size inconsistent.");
            }
            nouts += outputs[o].dim_size(1);
        }
        //std::cout << "samples " << samples << " nouts " << nouts << std::endl;

        std::vector< std::vector< float > > result;
        result.resize(samples, std::vector< float >(nouts));

        size_t idx0 = 0;
        for (size_t o = 0; o < outputs.size(); ++o)
        {
            auto output_map = outputs[o].tensor<float, 2>();

            size_t n = outputs[o].dim_size(1);
            for (size_t s = 0; s < samples; ++s) {
                std::vector< float > & vs = result[s];
                for (size_t i = 0; i < n; ++i) {
                    vs[idx0 + i] = output_map(s, i);
                }
            }
            idx0 += n;
        }
        return result;
    }
    else
    {
        std::cout << status.ToString() << std::endl;
        return std::vector< std::vector< float > >();
    }
}
// -------------------------------------------------------------------

