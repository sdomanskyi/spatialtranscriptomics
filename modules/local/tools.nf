/* 
 * Spatially variable genes with SpatialDE
 */
process ST_SPATIALDE {

    tag "$sample_id"
    label 'python_process'
    maxRetries 0
    errorStrategy { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: '{*.csv,show/*.png}', saveAs: { "${file(it).getFileName()}" }, mode: 'copy', overwrite: false

    input:
    tuple val(sample_id), path(h5ad)
    
    output:
    tuple val(sample_id), file('show/*.png'), emit: png
    tuple val(sample_id), file('*.csv'), emit: csv
    
    script:
    """
    #!/usr/bin/env python
   
    import os
    os.environ["NUMBA_CACHE_DIR"] = "./tmp"

    import scanpy as sc
    import numpy as np
    import pandas as pd
    import SpatialDE

    # See more settings at:
    # https://scanpy.readthedocs.io/en/stable/generated/scanpy._settings.ScanpyConfig.html
    sc.settings.figdir = './'
    
    if not os.path.exists('./show/'):
        os.makedirs('./show/')
    
    st_adata = sc.read('${h5ad}')
    print(st_adata.shape)
    
    counts = pd.DataFrame(st_adata.X.todense(), columns=st_adata.var_names, index=st_adata.obs_names)
    coord = pd.DataFrame(st_adata.obsm['spatial'], columns=['x_coord', 'y_coord'], index=st_adata.obs_names)
    
    df_results = SpatialDE.run(coord, counts)
    
    df_results.index = df_results["g"]
    df_results = df_results.sort_values("qval", ascending=True)
    
    df_results.to_csv('stSpatialDE.csv')
    
    # Plotting top most-HVG
    keys = df_results.index.values[:${params.SpatialDE_plotTopHVG}]
    sc.pl.spatial(st_adata, img_key="hires", color=keys, alpha=0.7, save='/stSpatialDE.png', ncols=${params.SpatialDE_numberOfColumns})
    """
}

/* 
 * Resolution enhancement and spatial clustering with BayesSpace
 */
process ST_BAYESSPACE {

    tag "$sample_id"
    label 'r_process'
    maxRetries 0
    errorStrategy { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: false

    input:
    tuple val(sample_id), path(csv)

    output:
    tuple val(sample_id), path('*.csv'), emit: csv
    tuple val(sample_id), path('*.png'), emit: png

    """
    Rscript ${projectDir}/bin/characterization_BayesSpace.R \
    --filePath=./ \
    --nameX=st_data.matrix.csv.gz \
    --nameVar=st_data.var.csv.gz \
    --nameObs=st_data.obs.csv.gz \
    --numberHVG=$params.BayesSpace_numberHVG \
    --numberPCs=$params.BayesSpace_numberPCs \
    --minClusters=$params.BayesSpace_minClusters \
    --maxClusters=$params.BayesSpace_maxClusters \
    --optimalQ=$params.BayesSpace_optimalQ \
    --STplatform=$params.BayesSpace_STplatform \
    --gamma=$params.BayesSpace_gamma
    """
}

/* 
 * ST data deconvolution with STdeconvolve
 */
process ST_STDECONVOLVE {

    tag "$sample_id"
    label 'r_process'
    maxRetries 0
    errorStrategy { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: false

    input:
    tuple val(sample_id), val(meta), path(csv)

    output:
    tuple val(sample_id), path('*LDA.rds'), emit: model
    tuple val(sample_id), path('*.csv'), emit: csv
    tuple val(sample_id), path('*.png'), emit: png

    """
    Rscript ${projectDir}/bin/characterization_STdeconvolve.R \
    --filePath=./ \
    --nameX=st_data.matrix.csv.gz \
    --nameVar=st_data.var.csv.gz \
    --nameObs=st_data.obs.csv.gz \
    --outsPath=${meta.st_data_dir} \
    --mtxGeneColumn=$params.STdeconvolve_mtxGeneColumn \
    --countsFactor=$params.STdeconvolve_countsFactor \
    --corpusRemoveAbove=$params.STdeconvolve_corpusRemoveAbove \
    --corpusRemoveBelow=$params.STdeconvolve_corpusRemoveBelow \
    --LDAminTopics=$params.STdeconvolve_LDAminTopics \
    --LDAmaxTopics=$params.STdeconvolve_LDAmaxTopics \
    --STdeconvolveScatterpiesSize=$params.STdeconvolve_ScatterpiesSize \
    --STdeconvolveFeaturesSizeFactor=$params.STdeconvolve_FeaturesSizeFactor \
    --trainedLDA="${params.STdeconvolve_LDA}"
    """
}

/* 
 * ST data deconvolution with SPOTlight
 */
process ST_SPOTLIGHT {

    tag "$sample_id"
    label 'r_process'
    maxRetries 0
    errorStrategy { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: false

    input:
    tuple val(sample_id), val(meta), path(csv)

    output:
    tuple val(sample_id), path('*NMF.rds'), emit: model
    tuple val(sample_id), path('*.csv'), emit: csv
    tuple val(sample_id), path('*.png'), emit: png

    """
    Rscript $projectDir/bin/characterization_SPOTlight.R \
    --filePath="" \
    --nameX=st_data.matrix.csv.gz \
    --nameVar=st_data.var.csv.gz \
    --nameObs=st_data.obs.csv.gz \
    --SCnameX=sc_data.matrix.csv.gz \
    --SCnameVar=sc_data.var.csv.gz \
    --SCnameObs=sc_data.obs.csv.gz \
    --annoDataDir=${meta.sc_annotation_data_dir} \
    --annoFileCounts=${meta.sc_annotation_counts} \
    --annoFileCelltype=${meta.sc_annotation_labels} \
    --outsPath=${meta.st_data_dir} \
    --mtxGeneColumn=$params.SPOTlight_mtxGeneColumn \
    --countsFactor=$params.SPOTlight_countsFactor \
    --clusterResolution=$params.SPOTlight_clusterResolution \
    --numberHVG=$params.SPOTlight_numberHVG \
    --numberCellsPerCelltype=$params.SPOTlight_numberCellsPerCelltype \
    --SPOTlightScatterpiesSize=$params.SPOTlight_ScatterpiesSize
    """
}
