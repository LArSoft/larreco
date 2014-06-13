#ifndef CMERGEMANAGER_CC
#define CMERGEMANAGER_CC

#include "CMergeManager.h"

namespace cluster {

  CMergeManager::CMergeManager(CMergePriority_t priority)
  {
    _fout = 0;
    _debug_mode = kNone;
    _merge_till_converge = false;
    _priority = priority;
    _merge_algo = nullptr;
    _separate_algo = nullptr;
    Reset();
  }

  void CMergeManager::Reset()
  {
    _in_clusters.clear();
    _out_clusters.clear();
    _book_keeper.Reset();
  }

  void CMergeManager::SetClusters(const std::vector<std::vector<util::PxHit> > &clusters)
  {

    // Clear & fill cluster info

    _in_clusters.clear();
    _out_clusters.clear();

    _in_clusters.reserve(clusters.size());

    ClusterParamsAlg tmp_alg;
    tmp_alg.SetMinNHits(0);
    tmp_alg.SetVerbose(false);

    for(auto const &c : clusters) {
      
      _in_clusters.push_back(tmp_alg);
      (*_in_clusters.rbegin()).Initialize();

      if((*_in_clusters.rbegin()).SetHits(c) < 1) continue;
      (*_in_clusters.rbegin()).DisableFANN();
      (*_in_clusters.rbegin()).FillParams(true,true,true,true,true,false);
      (*_in_clusters.rbegin()).FillPolygon();

    }    

    // Clear book keeper
    
    _book_keeper.Reset(clusters.size());

  }
  
  void CMergeManager::Process()
  {
    if(!_merge_algo) throw CRUException("No algorithm to run!");

    _merge_algo->SetAnaFile(_fout);
    if(_separate_algo) _merge_algo->SetAnaFile(_fout);
    
    _merge_algo->EventBegin(_in_clusters);
    if(_separate_algo) _separate_algo->EventBegin(_in_clusters);
    
    size_t ctr=0;
    bool keep_going=true;
    std::vector<CBookKeeper> book_keepers;
    std::vector<std::vector<unsigned short> > tmp_merged_indexes;
    std::vector<cluster::ClusterParamsAlg> tmp_merged_clusters;
    while(keep_going){
      
      // Configure input for RunMerge
      CBookKeeper bk;

      if(!ctr) tmp_merged_clusters = _in_clusters;
      else tmp_merged_clusters = _out_clusters;
      _out_clusters.clear();

      if(_debug_mode <= kPerIteration)
	
	std::cout 
	  << std::endl
	  << Form("  Merging iteration %zu: processing %zu clusters...",ctr,tmp_merged_clusters.size()) 
	  << std::endl;

      bk.Reset(tmp_merged_clusters.size());

      std::vector<bool> merge_switch(tmp_merged_clusters.size(),true);
      for(size_t i=0; i<tmp_merged_indexes.size(); ++i)
	
	if(tmp_merged_indexes.at(i).size()==1)
	  
	  merge_switch.at(i) = false;

      // Run separation algorithm
      if(_separate_algo)

	RunSeparate(tmp_merged_clusters, bk);

      // Run merging algorithm
      RunMerge(tmp_merged_clusters, merge_switch, bk);

      // Save output
      bk.PassResult(tmp_merged_indexes);

      if(bk.size() == tmp_merged_indexes.size())
	_out_clusters = tmp_merged_clusters;
      else {
	_out_clusters.reserve(tmp_merged_indexes.size());
	for(auto const& indexes_v : tmp_merged_indexes) {
	  
	  if(indexes_v.size()==1) {
	    _out_clusters.push_back(tmp_merged_clusters.at(indexes_v.at(0)));
	    continue;
	  }
	      
	  size_t tmp_hit_counts=0;
	  for(auto const& index : indexes_v) 
	    tmp_hit_counts += tmp_merged_clusters.at(index).GetHitVector().size();
	  std::vector<util::PxHit> tmp_hits;
	  tmp_hits.reserve(tmp_hit_counts);
	  
	  for(auto const& index : indexes_v) {
	    for(auto const& hit : tmp_merged_clusters.at(index).GetHitVector())
	      tmp_hits.push_back(hit);
	  }
	  _out_clusters.push_back(ClusterParamsAlg());
	  (*_out_clusters.rbegin()).SetVerbose(false);
	  (*_out_clusters.rbegin()).DisableFANN();

	  if((*_out_clusters.rbegin()).SetHits(tmp_hits) < 1) continue;
	  (*_out_clusters.rbegin()).FillParams(true,true,true,true,true,false);
	  (*_out_clusters.rbegin()).FillPolygon();
	}
	book_keepers.push_back(bk);
      }

      if(_debug_mode <= kPerIteration) {

	std::cout
	  << Form("  Input / Output cluster length: %zu/%zu",tmp_merged_clusters.size(),_out_clusters.size())
	  << std::endl;

	if(tmp_merged_clusters.size() == _out_clusters.size())

	  std::cout << "  Did not find any newly merged clusters..." <<std::endl;
	
	if(_out_clusters.size()==1)
	  
	  std::cout << "  Output cluster length is 1 (= no more merging needed)" << std::endl;

	if(!_merge_till_converge) 

	  std::cout << "  Non-iterative option: terminating..." << std::endl;

      }

      // Break if necessary
      if(!_merge_till_converge || tmp_merged_clusters.size() == _out_clusters.size() || _out_clusters.size()==1)

	break;

      ctr++;
    }
    
    // Gather the full book keeping result
    for(auto const& bk : book_keepers)

      _book_keeper.Combine(bk);

    _merge_algo->EventEnd();
    if(_separate_algo) _separate_algo->EventEnd();
    
  }

  void CMergeManager::RunMerge(const std::vector<cluster::ClusterParamsAlg> &in_clusters,
			       CBookKeeper &book_keeper) const
  {
    RunMerge(in_clusters,
	     std::vector<bool>(in_clusters.size(),true),
	     book_keeper);
  }

  void CMergeManager::RunMerge(const std::vector<cluster::ClusterParamsAlg> &in_clusters,
			       const std::vector<bool> &merge_flag,
			       CBookKeeper &book_keeper) const
  {
    if(merge_flag.size() != in_clusters.size())
      throw CRUException(Form("in_clusters (%zu) and merge_flag (%zu) vectors must be of same length!",
				   in_clusters.size(),
				   merge_flag.size()
				   )
			      );
    if(_debug_mode <= kPerIteration){
      
      std::cout
	<< Form("  Calling %s with %zu clusters...",__FUNCTION__,in_clusters.size())
	<<std::endl;

    }

    // Figure out ordering of clusters to process
    std::multimap<double,size_t> prioritized_index;
    for(size_t i=0; i<in_clusters.size(); ++i) {

      if(in_clusters.at(i).GetHitVector().size()<1) continue;

      switch(_priority) {
      case ::cluster::CMergeManager::kIndex:
	prioritized_index.insert(std::pair<double,size_t>(((double)(in_clusters.size()) - (double)i),i));
	break;
      case ::cluster::CMergeManager::kPolyArea:
	prioritized_index.insert(std::pair<double,size_t>(in_clusters.at(i).GetParams().PolyObject.Area(),i));
	break;
      case ::cluster::CMergeManager::kChargeSum:
	prioritized_index.insert(std::pair<double,size_t>(in_clusters.at(i).GetParams().sum_charge,i));
	break;
      case ::cluster::CMergeManager::kNHits:
	prioritized_index.insert(std::pair<double,size_t>((double)(in_clusters.at(i).GetNHits()),i));
	break;
      }
    }

    //
    // Merging
    //
    // Call prepare function for merging
    _merge_algo->IterationBegin(in_clusters);
    
    // Run over clusters and execute merging algorithms
    for(auto citer1 = prioritized_index.rbegin();
	citer1 != prioritized_index.rend();
	++citer1) {
      
      auto citer2 = citer1;

      UChar_t plane1 = in_clusters.at((*citer1).second).Plane();

      while(1) {
	citer2++;
	if(citer2 == prioritized_index.rend()) break;

	// Skip if not on the same plane
	UChar_t plane2 = in_clusters.at((*citer2).second).Plane();
	if(plane1 != plane2) continue;

	// Skip if this combination is not meant to be compared
	if(!(merge_flag.at((*citer2).second)) && !(merge_flag.at((*citer1).second)) ) continue;

	// Skip if this combination is not allowed to merge
	if(!(book_keeper.MergeAllowed((*citer1).second,(*citer2).second))) continue;

	if(_debug_mode <= kPerAlgoSet){
	  
	  std::cout
	    << Form("    \033[93mInspecting a pair (%zu, %zu) for merging... \033[00m",(*citer1).second, (*citer2).second)
	    << std::endl;
	}
	
	bool merge = _merge_algo->Bool(in_clusters.at((*citer1).second),in_clusters.at((*citer2).second));

	if(_debug_mode <= kPerAlgoSet) {
	  
	  if(merge) 
	    std::cout << "    \033[93mfound to be merged!\033[00m " 
		      << std::endl
		      << std::endl;
	  
	  else 
	    std::cout << "    \033[93mfound NOT to be merged...\033[00m" 
		      << std::endl
		      << std::endl;

	} // end looping over all sets of algorithms
	
	if(merge)

	  book_keeper.Merge((*citer1).second,(*citer2).second);

      } // end looping over all cluster pairs for citer1

    } // end looping over clusters

    _merge_algo->IterationEnd();
    
    if(_debug_mode <= kPerIteration && book_keeper.GetResult().size() != in_clusters.size()) {
      
      std::cout << "  Found following clusters to be merged..." << std::endl;
      for(auto const &indexes_v : book_keeper.GetResult()) {

	if(indexes_v.size()==1) continue;
	std::cout<<"    ";
	for(auto index : indexes_v)

	  std::cout<<index<<" ";
	std::cout<<" ... indexes to be merged!"<<std::endl;

      }
    }

  }

  void CMergeManager::RunSeparate(const std::vector<cluster::ClusterParamsAlg> &in_clusters,
				  CBookKeeper &book_keeper) const
  {
    /*
    if(separate_flag.size() != in_clusters.size())
      throw CRUException(Form("in_clusters (%zu) and separate_flag (%zu) vectors must be of same length!",
				   in_clusters.size(),
				   separate_flag.size()
				   )
			      );
    */
    if(_debug_mode <= kPerIteration){
      
      std::cout
	<< Form("  Calling %s with %zu clusters...",__FUNCTION__,in_clusters.size())
	<<std::endl;

    }

    // Figure out ordering of clusters to process
    std::multimap<double,size_t> prioritized_index;
    for(size_t i=0; i<in_clusters.size(); ++i) {

      if(in_clusters.at(i).GetHitVector().size()<1) continue;

      switch(_priority) {
      case ::cluster::CMergeManager::kIndex:
	prioritized_index.insert(std::pair<double,size_t>(((double)(in_clusters.size()) - (double)i),i));
	break;
      case ::cluster::CMergeManager::kPolyArea:
	prioritized_index.insert(std::pair<double,size_t>(in_clusters.at(i).GetParams().PolyObject.Area(),i));
	break;
      case ::cluster::CMergeManager::kChargeSum:
	prioritized_index.insert(std::pair<double,size_t>(in_clusters.at(i).GetParams().sum_charge,i));
	break;
      case ::cluster::CMergeManager::kNHits:
	prioritized_index.insert(std::pair<double,size_t>((double)(in_clusters.at(i).GetNHits()),i));
	break;
      }
    }

    //
    // Separation
    //

    // Call prepare functions for separation
    _separate_algo->IterationBegin(in_clusters);
    
    // Run over clusters and execute merging algorithms
    for(size_t cindex1 = 0; cindex1 < in_clusters.size(); ++cindex1) {
      
      UChar_t plane1 = in_clusters.at(cindex1).Plane();
      
      for(size_t cindex2 = cindex1+1; cindex2 < in_clusters.size(); ++cindex2) {
	
	// Skip if not on the same plane
	UChar_t plane2 = in_clusters.at(cindex2).Plane();
	if(plane1 != plane2) continue;
	
	// Skip if this combination is not meant to be compared
	//if(!(separate_flag.at(cindex2))) continue;
	
	if(_debug_mode <= kPerAlgoSet){
	  
	  std::cout
	    << Form("    \033[93mInspecting a pair (%zu, %zu) for separation... \033[00m",cindex1,cindex2)
	    << std::endl;
	}
	
	bool separate = _separate_algo->Bool(in_clusters.at(cindex1),in_clusters.at(cindex2));
	
	if(_debug_mode <= kPerAlgoSet) {
	  
	  if(separate) 
	    std::cout << "    \033[93mfound to be separated!\033[00m " 
		      << std::endl
		      << std::endl;
	  
	  else 
	    std::cout << "    \033[93mfound NOT to be separated...\033[00m" 
		      << std::endl
		      << std::endl;
	  
	} // end looping over all sets of algorithms
	
	if(separate)
	  
	  book_keeper.ProhibitMerge(cindex1,cindex2);
	  
      } // end looping over all cluster pairs for citer1
      
    } // end looping over clusters

    _separate_algo->IterationEnd();
  }

}

#endif
