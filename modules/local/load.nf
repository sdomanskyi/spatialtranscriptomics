
process PREPARE_ST_DATA {

    tag "$sample_id"
    label 'python_process_low'
    maxRetries 0
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: false

    input:
    tuple val(sample_id), val(meta)
    
    output:
    tuple val(sample_id), file("st_data.raw.h5ad")
    
    script:
    """
    #!/usr/bin/env python
    
    import os
    os.environ["NUMBA_CACHE_DIR"] = "./tmp"

    import sys
    sys.path.append("${projectDir}/lib")
    from utils import read_visium_mtx

    import scanpy as sc

    st_adata = read_visium_mtx("${meta.st_data_dir}")

    st_adata.var_names_make_unique()
    sc.pp.filter_cells(st_adata, min_counts=${params.STload_minCounts})
    sc.pp.filter_cells(st_adata, min_genes=${params.STload_minGenes})
    sc.pp.filter_genes(st_adata, min_cells=${params.STload_minCells})

    st_adata.write('st_data.raw.h5ad')                      
    """
}


process PREPARE_SC_DATA {

    tag "$sample_id"
    label 'python_process_low'
    maxRetries 0
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: false

    input:
    tuple val(sample_id), val(meta)
    
    output:
    tuple val(sample_id), file("sc_data.raw.h5ad")
    
    script:
    """
    #!/usr/bin/env python

    import os
    os.environ["NUMBA_CACHE_DIR"] = "./tmp"

    import scanpy as sc
    import numpy as np
    import pandas as pd
    
    outsPath = '${meta.sc_data_dir}filtered_feature_bc_matrix/'

    sc_adata = sc.read_mtx(outsPath + 'matrix.mtx.gz').T   
    print(sc_adata.shape)
    
    genes = pd.read_csv(outsPath + 'features.tsv.gz', header=None, sep='\t')  
    print(genes)
    
    if len(genes.columns)==1:
        genes[0] = genes[0].str.replace('_', '-')
        gs = genes[0]
    else:
        genes[1] = genes[1].str.replace('_', '-')
        gs = genes[1]
    
    sc_adata.var_names = gs.values    
    sc_adata.var['gene_symbols'] = gs.values   
    sc_adata.obs_names = pd.read_csv(outsPath + 'barcodes.tsv.gz', header=None)[0].values

    sc_adata.var_names_make_unique()
    sc.pp.filter_cells(sc_adata, min_genes=${params.SCload_minGenes})
    sc.pp.filter_cells(sc_adata, min_counts=${params.SCload_minCounts})
    sc.pp.filter_genes(sc_adata, min_cells=${params.SCload_minCells})

    sc_adata.write('sc_data.raw.h5ad')             
    """
}


process QC_FILTER_ST_DATA {

    tag "$sample_id"
    label 'python_process_low'
    maxRetries 0
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: '{*.png,show/*.png,violin/*.png}', saveAs: { "${file(it).getFileName()}" }, mode: 'copy', overwrite: false
    
    input:
    tuple val(sample_id), val(meta), path(h5ad)
    
    output:
    tuple val(sample_id), file('st_data.h5ad'), emit: h5ad
    tuple val(sample_id), file('{*.png,show/*.png,violin/*.png}'), emit: png
    
    script:
    """
    #!/usr/bin/env python

    import os
    os.environ["NUMBA_CACHE_DIR"] = "./tmp"
    
    import sys
    sys.path.append("${projectDir}/lib")
    from utils import histplotQC    

    sys.path.append("${projectDir}/")
    from importlib import resources
    import scanpy as sc
    import numpy as np
    import pandas as pd
    from matplotlib import pyplot as plt

    # See more settings at:
    # https://scanpy.readthedocs.io/en/stable/generated/scanpy._settings.ScanpyConfig.html
    sc.settings.figdir = './'
    
    pltFigSize = ${params.STpreprocess_pltFigSize}
    plt.rcParams["figure.figsize"] = (pltFigSize, pltFigSize)
    
    if not os.path.exists('./violin/'):
        os.makedirs('./violin/')
        
    if not os.path.exists('./show/'):
        os.makedirs('./show/')
    
    st_adata = sc.read('${h5ad}')
    print(st_adata.shape)
    
    print("${projectDir}/")
    species = '${meta.species}'
    with resources.path('assets', f'{species}.MitoCarta3.0.Symbol.csv.gz') as mitoPath:
        mito = pd.read_csv(mitoPath, index_col=0, header=None).index.values
        print(len(mito))
    
    hg = st_adata.var_names.str.startswith('GRCh38_')
    if (species == 'Human') and (hg.any()):
        st_adata.var_names = pd.Index([v[7:] for v in st_adata.var_names])
        st_adata = st_adata[:, hg]
  
    mm = st_adata.var_names.str.startswith('mm10___')
    if (species == 'Mouse') and (mm.any()):
        st_adata.var_names = pd.Index([v[7:] for v in st_adata.var_names])
        st_adata = st_adata[:, mm]

    st_adata.var_names_make_unique()

    st_adata.var["mt"] = st_adata.var_names.isin(mito)
    sc.pp.calculate_qc_metrics(st_adata, qc_vars=["mt"], inplace=True)
    
    print(st_adata.var["mt"].value_counts())

    st_adata.obs['in_tissue_cat'] = st_adata.obs['in_tissue'].replace({1: 'in', 0:'out'}).astype('category')
    
    sc.pl.violin(st_adata, ["pct_counts_mt", "total_counts", "n_genes_by_counts"],
                 jitter=0.4, groupby='in_tissue_cat', rotation=0, save='/stQCViolin.png')
                 
    in_tissue = st_adata.obs['in_tissue'].astype(int) == 1

    keys = ["in_tissue", "pct_counts_mt", "total_counts", "n_genes_by_counts"]
    sc.pl.spatial(st_adata[in_tissue], img_key="hires", color=keys, save='/st_QC_in.png')
    sc.pl.spatial(st_adata[~in_tissue], img_key="hires", color=keys, save='/st_QC_out.png')  
    
    bins = ${params.STpreprocess_histplotQCbins}
    tots = ${params.STpreprocess_histplotQCmaxTotalCounts}
    ming = ${params.STpreprocess_histplotQCminGeneCounts}
    
    fig, axs = plt.subplots(1, 5, figsize=(pltFigSize*5, pltFigSize))
    histplotQC(st_adata.obs["total_counts"], bins=bins, ax=axs[0])
    histplotQC(st_adata.obs["total_counts"][st_adata.obs["total_counts"] < tots], bins=bins, ax=axs[1])
    histplotQC(st_adata.obs["n_genes_by_counts"], bins=bins, ax=axs[2])
    histplotQC(st_adata.obs["n_genes_by_counts"][st_adata.obs["n_genes_by_counts"] < ming], bins=bins, ax=axs[3])
    histplotQC(st_adata.obs["pct_counts_mt"], bins=bins, ax=axs[4])
    fig.tight_layout()
    fig.savefig('st_histogrtam_all.png', facecolor='w')
    
    fig, axs = plt.subplots(1, 5, figsize=(pltFigSize*5, pltFigSize))
    histplotQC(st_adata[in_tissue].obs["total_counts"], bins=bins, ax=axs[0])
    histplotQC(st_adata[in_tissue].obs["total_counts"][st_adata[in_tissue].obs["total_counts"] < tots], bins=bins, ax=axs[1])
    histplotQC(st_adata[in_tissue].obs["n_genes_by_counts"], bins=bins, ax=axs[2])
    histplotQC(st_adata[in_tissue].obs["n_genes_by_counts"][st_adata[in_tissue].obs["n_genes_by_counts"] < ming], bins=bins, ax=axs[3])
    histplotQC(st_adata[in_tissue].obs["pct_counts_mt"], bins=bins, ax=axs[4])
    fig.tight_layout()
    fig.savefig('st_histogrtam_in.png', facecolor='w')
    plt.close(fig)
 
    st_adata = st_adata[in_tissue]
    sc.pp.filter_cells(st_adata, min_counts=${params.STpreprocess_minCounts})
    sc.pp.filter_cells(st_adata, min_genes=${params.STpreprocess_minGenes})
    sc.pp.filter_genes(st_adata, min_cells=${params.STpreprocess_minCells})

    st_adata.write('st_data.h5ad')                   
    """
}


process QC_FILTER_SC_DATA {

    tag "$sample_id"
    label 'python_process_low'
    maxRetries 0
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: '{*.png,show/*.png,violin/*.png}', saveAs: { "${file(it).getFileName()}" }, mode: 'copy', overwrite: false

    input:
    tuple val(sample_id), val(meta), path(h5ad)
    
    output:
    tuple val(sample_id), file('sc_data.h5ad'), emit: h5ad
    tuple val(sample_id), file('{*.png,show/*.png,violin/*.png}'), emit: png
    
    script:
    """
    #!/usr/bin/env python

    import os
    os.environ["NUMBA_CACHE_DIR"] = "./tmp"
    
    import sys
    sys.path.append("${projectDir}/lib")
    from utils import histplotQC    

    sys.path.append("${projectDir}/")
    from importlib import resources
    import scanpy as sc
    import numpy as np
    import pandas as pd
    from matplotlib import pyplot as plt

    # See more settings at:
    # https://scanpy.readthedocs.io/en/stable/generated/scanpy._settings.ScanpyConfig.html
    sc.settings.figdir = './'
    
    pltFigSize = ${params.STpreprocess_pltFigSize}
    plt.rcParams["figure.figsize"] = (pltFigSize, pltFigSize)
    
    if not os.path.exists('./violin/'):
        os.makedirs('./violin/')
        
    if not os.path.exists('./show/'):
        os.makedirs('./show/')
    
    sc_adata = sc.read('${h5ad}')
    print(sc_adata.shape)
    
    print("${projectDir}/")
    species = '${meta.species}'
    with resources.path('assets', f'{species}.MitoCarta3.0.Symbol.csv.gz') as mitoPath:
        mito = pd.read_csv(mitoPath, index_col=0, header=None).index.values
        print(len(mito))

    hg = sc_adata.var_names.str.startswith('GRCh38_')
    if (species == 'Human') and (hg.any()):
        sc_adata.var_names = pd.Index([v[7:] for v in sc_adata.var_names])
        sc_adata = sc_adata[:, hg]
        
    mm = sc_adata.var_names.str.startswith('mm10___')
    if (species == 'Mouse') and (mm.any()):
        sc_adata.var_names = pd.Index([v[7:] for v in sc_adata.var_names])
        sc_adata = sc_adata[:, mm]

    sc_adata.var_names_make_unique()

    sc_adata.var["mt"] = sc_adata.var_names.isin(mito)
    sc.pp.calculate_qc_metrics(sc_adata, qc_vars=["mt"], inplace=True)

    print(sc_adata.var.index, mito)
    print(sc_adata.var["mt"].value_counts())
    
    bins = ${params.SCpreprocess_histplotQCbins}
    tots = ${params.SCpreprocess_histplotQCmaxTotalCounts}
    ming = ${params.SCpreprocess_histplotQCminGeneCounts}
    
    fig, axs = plt.subplots(1, 5, figsize=(pltFigSize*5, pltFigSize))
    histplotQC(sc_adata.obs["total_counts"], bins=bins, ax=axs[0])
    histplotQC(sc_adata.obs["total_counts"][sc_adata.obs["total_counts"] < tots], bins=bins, ax=axs[1])
    histplotQC(sc_adata.obs["n_genes_by_counts"], bins=bins, ax=axs[2])
    histplotQC(sc_adata.obs["n_genes_by_counts"][sc_adata.obs["n_genes_by_counts"] < ming], bins=bins, ax=axs[3])
    histplotQC(sc_adata.obs["pct_counts_mt"], bins=bins, ax=axs[4])
    fig.tight_layout()
    fig.savefig('sc_histogrtam_all.png', facecolor='w')

    sc.pp.filter_cells(sc_adata, min_counts=${params.SCpreprocess_minCounts})
    sc.pp.filter_cells(sc_adata, min_genes=${params.SCpreprocess_minGenes})
    sc.pp.filter_genes(sc_adata, min_cells=${params.SCpreprocess_minCells})

    sc_adata.write('sc_data.h5ad')                     
    """
}




