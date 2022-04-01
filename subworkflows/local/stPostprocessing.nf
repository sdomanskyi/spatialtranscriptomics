nextflow.enable.dsl=2

/* 
 * Include requires tasks 
 */
include { ST_CLUSTERING               } from '../../modules/local/tasks'
include { ALL_REPORT                  } from '../../modules/local/tasks'
include { RMARKDOWNNOTEBOOK           } from '../../modules/nf-core/modules/rmarkdownnotebook/main'

def create_meta_map(item) {
    def meta = [:]
    meta.id = item[0]
    return [meta, projectDir + params.rmarkdown_template] 
}

   
/* 
 * Run postprocessing tools
 */
workflow ST_POSTPROCESSING {
 
    take:
      sample_ids
      outdir
      
    main:
      ST_CLUSTERING(  sample_ids, outdir)
      RMARKDOWNNOTEBOOK( ST_CLUSTERING.out.map{ create_meta_map(it) }, params, [])

    emit:
      ST_CLUSTERING.out
 }
