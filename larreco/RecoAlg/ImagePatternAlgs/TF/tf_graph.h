////////////////////////////////////////////////////////////////////////////////////////////////////
//// Class:       Graph
//// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
////
//// Iterface to run Tensorflow graph saved to a file. First attempts, not functional yet.
////
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef Graph_h
#define Graph_h

#include <memory>
#include <vector>
#include <string>

namespace tensorflow
{
    class Session;
}

namespace tf
{

class Graph
{
public:
    static std::unique_ptr<Graph> create(const char* graph_file_name)
    {
        bool success;
        std::unique_ptr<Graph> ptr(new Graph(graph_file_name, success));
        if (success) { return ptr; }
        else { return nullptr; }
    }

    ~Graph();

    std::vector<float> run(const std::vector< std::vector<float> > & x);

private:
    /// Not-throwing constructor.
    Graph(const char* graph_file_name, bool & success);

    tensorflow::Session* fSession;
    std::string fInputName, fOutputName;
};

} // namespace tf

#endif
