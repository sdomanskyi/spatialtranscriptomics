
process ST_POSTPROCESSING {

    tag "$sample_id"
    label 'python_process'
    maxRetries 0
    errorStrategy { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: '{show,umap,violin,heatmap,umap_density_clusters_}/*.png', saveAs: { "${file(it).getFileName()}" }, mode: 'copy', overwrite: false
    publishDir "${params.outdir}/${sample_id}", pattern: '{*processed.h5ad}', mode: 'copy', overwrite: false

    input:
    tuple val(sample_id), path(h5ad), path(csv)
    
    output:
    tuple val(sample_id), file('{show,umap,violin,heatmap,umap_density_clusters_}/*.png'), emit: png
    tuple val(sample_id), file('*processed.h5ad'), emit: h5ad

    """
    python ${projectDir}/bin/stClusteringWorkflow.py \
    --resolution=$params.Clustering_resolution \
    --nHVGs=$params.Clustering_numberHVG
    """
}
   