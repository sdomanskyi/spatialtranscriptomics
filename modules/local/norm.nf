
process NORMALIZE_CPM {

    tag "$sample_id"
    label 'python_process_low'
    maxRetries 0
    errorStrategy { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    
    input:
    tuple val(sample_id), path(h5ad)
    val(prefix)
    
    output:
    tuple val(sample_id), file('*.norm.h5ad'), emit: h5ad
    tuple val(sample_id), path('*.{matrix,obs,var}.csv.gz', arity: '3'), emit: csv
    
    script:
    """
    #!/usr/bin/env python

    import os
    os.environ["NUMBA_CACHE_DIR"] = "./tmp"
    import scanpy as sc
    import pandas as pd
    import numpy as np

    adata = sc.read('${h5ad}')

    sc.pp.normalize_total(adata, target_sum=10**6)
    sc.pp.log1p(adata)

    adata.write('${prefix}.norm.h5ad')  

    adata.var.to_csv('${prefix}.var.csv.gz')
    adata.obs.to_csv('${prefix}.obs.csv.gz')
    pd.DataFrame(np.array(adata.X.todense()).T).to_csv('${prefix}.matrix.csv.gz')                     
    """
}


process RESAMPLE_OBS {

    tag "$sample_id"
    label 'python_process_low'
    maxRetries 0
    errorStrategy { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    
    input:
    tuple val(sample_id), path(h5ad)
    val(prefix)
    
    output:
    tuple val(sample_id), file('*.resampled.h5ad'), emit: h5ad
    tuple val(sample_id), path('*.{matrix,obs,var}.csv.gz', arity: '3'), emit: csv
    
    script:
    """
    #!/usr/bin/env python

    import os
    os.environ["NUMBA_CACHE_DIR"] = "./tmp"
    
    import sys
    sys.path.append("${projectDir}/lib")
    from utils import resample_counts_inplace  

    import scanpy as sc
    import pandas as pd
    import numpy as np

    adata = sc.read('${h5ad}')

    resample_counts_inplace(adata, target_total=${params.target_total})

    adata.write('${prefix}.resampled.h5ad')  

    adata.var.to_csv('${prefix}.var.csv.gz')
    adata.obs.to_csv('${prefix}.obs.csv.gz')
    pd.DataFrame(np.array(adata.X.todense()).T).to_csv('${prefix}.matrix.csv.gz')                     
    """
}

  
process NORMALIZE_SCT {

    tag "$sample_id"
    label 'r_process'
    maxRetries 0
    errorStrategy { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    
    input:
    tuple val(sample_id), path(csv)
    val(prefix)
    
    output:
    tuple val(sample_id), path('*.{sct-matrix,sct-var,sct-obs}.csv.gz', arity: '3'), emit: csv
    
    script:
    """
    #!/usr/local/bin/Rscript

    library(Seurat)

    matrix_ct <- read.csv(paste0('${prefix}', '.matrix.csv.gz'), row.names=1)
    genes <- read.csv(paste0('${prefix}', '.var.csv.gz'), row.names=1)
    obs <- read.csv(paste0('${prefix}', '.obs.csv.gz'), row.names=1)

    #rownames(matrix_ct) <- rownames(genes)
    #colnames(matrix_ct) <- rownames(obs)
    print(dim(matrix_ct))

    se <- Seurat::CreateSeuratObject(counts=matrix_ct, min.cells=1, min.genes=1, project="data")
    print(dim(se))

    se <- Seurat::SCTransform(object=se, return.only.var.genes=FALSE, do.correct.umi=TRUE,
                              verbose=FALSE, variable.features.n=2000)
    write.csv(as.matrix(se@assays\$SCT@data), file=gzfile(paste0('${prefix}', '.sct-matrix.csv.gz')))
    write.csv(se@assays[["SCT"]]@data@Dimnames[[1]], file=gzfile(paste0('${prefix}', '.sct-var.csv.gz')))
    write.csv(se@assays[["SCT"]]@data@Dimnames[[2]], file=gzfile(paste0('${prefix}', '.sct-obs.csv.gz')))

    quit(status=0)
    """
}


process UPDATE_ANNDATA {

    tag "$sample_id"
    label 'python_process_low'
    maxRetries 0
    errorStrategy { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    
    input:
    tuple val(sample_id), path(h5ad), path(csv)
    val(prefix)
    
    output:
    tuple val(sample_id), file('*.norm.h5ad'), emit: h5ad
    tuple val(sample_id), path('*.{matrix,obs,var}.csv.gz', arity: '3'), emit: csv
    
    script:
    """
    #!/usr/bin/env python

    import os
    os.environ["NUMBA_CACHE_DIR"] = "./tmp"
    
    import sys
    sys.path.append("${projectDir}/lib")

    import scanpy as sc
    import pandas as pd
    import numpy as np
    from scipy.sparse import csr_matrix

    adata = sc.read('${h5ad}')

    ind_var = pd.read_csv('${prefix}.sct-var.csv.gz', index_col=0)['x'].values
    ind_obs = pd.read_csv('${prefix}.sct-obs.csv.gz', index_col=0)['x'].values
    ind_obs = np.array([v[1:] for v in ind_obs]).astype(int)
    
    adata = adata[adata.obs.index[ind_obs], adata.var.index[ind_var]]
    
    X = pd.read_csv('${prefix}.sct-matrix.csv.gz', index_col=0).values.T
    if type(adata.X) is csr_matrix:
        X = csr_matrix(X)
    adata.X = X

    adata.write('${prefix}.norm.h5ad')  

    adata.var.to_csv('${prefix}.var.csv.gz')
    adata.obs.to_csv('${prefix}.obs.csv.gz')
    pd.DataFrame(np.array(adata.X.todense()).T).to_csv('${prefix}.matrix.csv.gz')                     
    """
}
