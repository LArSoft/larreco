add_subdirectory(Segmentation2D)

simple_plugin(TrackShowerHits   "module"
			lardataobj_RawData
		 	lardataobj_RecoBase
			larreco_RecoAlg
			larreco_ClusterFinder_TrackShowerSplitter_Segmentation2D
		 	larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
               art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
         		art_root_io::tfile_support ROOT::Core
    			art_root_io::TFileService_service
 			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities
			canvas::canvas
			messagefacility::MF_MessageLogger

			cetlib::cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}

)

install_headers()
install_fhicl()
install_source()
install_scripts()
