
include { ST_LOAD_PREPROCESS_DATA;
          SC_LOAD_PREPROCESS_DATA;
        } from '../subworkflows/local/stLoadPreprocessData'
        
include { ST_SPATIALDE;
          ST_BAYESSPACE;
          ST_STDECONVOLVE;
          ST_SPOTLIGHT;
        } from '../modules/local/tools'

include { ST_POSTPROCESSING; 
        } from '../modules/local/post'

include { RMARKDOWNNOTEBOOK;
          EXPORT_PARAMETERS;
          EXPORT_SAMPLEINFO;
        } from '../modules/local/report'

workflow ST {

    take:
        samples

    main:
        EXPORT_PARAMETERS ()
        EXPORT_SAMPLEINFO ( samples )

        ST_LOAD_PREPROCESS_DATA ( samples )
        SC_LOAD_PREPROCESS_DATA ( samples )
          
        ST_LOAD_PREPROCESS_DATA.out
        .map( { it -> [ it[0], it[1] ] } )
        .set{ h5adst }

        ST_SPATIALDE( h5adst )
        
        ST_LOAD_PREPROCESS_DATA.out
        .map( { it -> [ it[0], it[2] ] } )
        .set{ csvst }

        ST_BAYESSPACE( csvst )
        
        ST_STDECONVOLVE ( samples
                          .join(csvst) )
                          
        SC_LOAD_PREPROCESS_DATA.out
        .map( { it -> [ it[0], it[2] ] } )
        .set{ csvsc }

        // ** SPOTlight *******************************************************
        // Samples should not have both scRNA-seq data and annotated scRNA-seq 
        // data specified. If both are specified in samplesheet, the annotated
        // scRNA-seq data will be used below. Only samples with matching
        // species scRNA-seq data are processed.
        csvs = csvst
               .concat(csvsc)
               .groupTuple()
               .map( { it -> [ it[0], it[1][1]==null ?
                                      it[1][0] : it[1][0] + it[1][1] ] } )

        String nameD = 'sc_annotation_data_dir'
        String nameC = 'sc_annotation_counts'
        String nameL = 'sc_annotation_labels'
        
        csvs_annotated = samples
                         .join(csvs)
                         .filter( { it[1][nameD]?.trim()!=null ?
                                   !it[1][nameD]?.trim().isEmpty() : false } )
        csvs_unannotated = samples
                           .join(csvs)
                           .filter( { it[2].size() > 3 } )
        
        //Add annotation files to the file set to get them staged
        csvs_annotated = csvs_annotated
                         .map( { it -> [it[0], it[1],
                                 it[2] + [it[1][nameD] + it[1][nameC],
                                          it[1][nameD] + it[1][nameL]]] } )

        input_spotlight = csvs_annotated.concat(csvs_unannotated)
            
        ST_SPOTLIGHT( input_spotlight )
        // ********************************************************************

        ST_SPATIALDE.out.csv.transpose()
        .concat(ST_STDECONVOLVE.out.csv.transpose())
        .concat(ST_SPOTLIGHT.out.csv.transpose())
        .concat(ST_BAYESSPACE.out.csv.transpose())
        .groupTuple()
        .set{ tools_csvs }
        
        ST_POSTPROCESSING ( h5adst
                            .join(tools_csvs) )
        
        // ** Report **********************************************************
        ST_LOAD_PREPROCESS_DATA.out
        .map( { it -> [ it[0], it[3] ] } )
        .set{ st_qc_pngs }
        
        SC_LOAD_PREPROCESS_DATA.out
        .map( { it -> [ it[0], it[3] ] } )
        .set{ sc_qc_pngs }
        
        st_qc_pngs.transpose()
        .concat(sc_qc_pngs.transpose())
        .concat(ST_POSTPROCESSING.out.png.transpose())
        .concat(ST_BAYESSPACE.out.png.transpose())
        .concat(ST_STDECONVOLVE.out.png.transpose())
        .concat(ST_SPOTLIGHT.out.png.transpose())
        .concat(ST_SPATIALDE.out.png.transpose())
        .groupTuple()
        .set{ pngs }

        RMARKDOWNNOTEBOOK ( samples
                            .join(pngs) )
        // ********************************************************************
}
