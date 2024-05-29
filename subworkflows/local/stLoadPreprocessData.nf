include { PREPARE_ST_DATA;
          PREPARE_SC_DATA;
          QC_FILTER_ST_DATA;
          QC_FILTER_SC_DATA;
        } from '../../modules/local/load'

include { NORMALIZE_CPM;
          RESAMPLE_OBS;
          NORMALIZE_SCT;
          UPDATE_ANNDATA;
        } from '../../modules/local/norm'

workflow ST_LOAD_PREPROCESS_DATA {
 
    take:
      samples
      
    main:
      PREPARE_ST_DATA ( samples )
      QC_FILTER_ST_DATA ( samples
                          .join( PREPARE_ST_DATA.out ) )

      prefix = 'st_data'
      
      if ( params.normalization == 'CPM' ) {
          NORMALIZE_CPM ( QC_FILTER_ST_DATA.out.h5ad, prefix)
          normalized_anndata = NORMALIZE_CPM.out.h5ad
          normalized_csv = NORMALIZE_CPM.out.csv
      }
      else if ( params.normalization == 'SCTransform' ) {
          RESAMPLE_OBS ( QC_FILTER_ST_DATA.out.h5ad, prefix )
          NORMALIZE_SCT ( RESAMPLE_OBS.out.csv, prefix )
          UPDATE_ANNDATA ( RESAMPLE_OBS.out.h5ad
                           .join(NORMALIZE_SCT.out.csv), prefix )
          normalized_anndata = UPDATE_ANNDATA.out.h5ad
          normalized_csv = UPDATE_ANNDATA.out.csv
      }

    emit:
      normalized_anndata
      .join(normalized_csv)
      .join(QC_FILTER_ST_DATA.out.png)
}


workflow SC_LOAD_PREPROCESS_DATA {
 
    take:
      samples
      
    main:
      samples_sc = samples.filter{ it[1]['sc_data_dir']!='' }

      PREPARE_SC_DATA ( samples_sc )
      QC_FILTER_SC_DATA ( samples_sc
                          .join( PREPARE_SC_DATA.out ) )
                          
      prefix = 'sc_data'
      
      if ( params.normalization == 'CPM' ) {
          NORMALIZE_CPM ( QC_FILTER_SC_DATA.out.h5ad, prefix)
          normalized_anndata = NORMALIZE_CPM.out.h5ad
          normalized_csv = NORMALIZE_CPM.out.csv
      }
      else if ( params.normalization == 'SCTransform' ) {
          RESAMPLE_OBS ( QC_FILTER_SC_DATA.out.h5ad, prefix )
          NORMALIZE_SCT ( RESAMPLE_OBS.out.csv, prefix )
          UPDATE_ANNDATA ( RESAMPLE_OBS.out.h5ad
                           .join(NORMALIZE_SCT.out.csv), prefix )
          normalized_anndata = UPDATE_ANNDATA.out.h5ad
          normalized_csv = UPDATE_ANNDATA.out.csv
      }

    emit:
      normalized_anndata
      .join(normalized_csv)
      .join(QC_FILTER_SC_DATA.out.png)
}
