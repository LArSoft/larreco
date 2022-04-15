cet_enable_asserts()

set( hitfinder_tool_lib_list
                        larreco_RecoAlg
                        larcorealg_Geometry
                        lardataobj_RecoBase
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        larvecutils::MarqFitAlg
                        nurandom::RandomUtils_NuRandomService_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art_root_io::tfile_support ROOT::Core
            art::Framework_Services_Optional_RandomNumberGenerator_service
                        art_root_io::TFileService_service
                        art::Framework_Services_System_FileCatalogMetadata_service
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        Boost::filesystem
    )

art_make(TOOL_LIBRARIES ${hitfinder_tool_lib_list} ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_source()
install_fhicl()
