import org.yaml.snakeyaml.Yaml
import org.yaml.snakeyaml.DumperOptions
import groovy.json.JsonOutput

/**
 * Multiline code blocks need to have the same indentation level
 * as the `script:` section. This function re-indents code to the specified level.
 */
def indent_code_block(code, n_spaces) {
    def indent_str = " ".multiply(n_spaces)
    return code.stripIndent().split("\n").join("\n" + indent_str)
}

/**
 * Create a config YAML file from a groovy map
 *
 * @params task The process' `task` variable
 * @returns a line to be inserted in the bash script.
 */
def dump_params_yml(params) {
    DumperOptions options = new DumperOptions();
    options.setDefaultFlowStyle(DumperOptions.FlowStyle.BLOCK);
    def yaml = new Yaml(options)
    def yaml_str = yaml.dump(params)

    // Writing the .params.yml file directly as follows does not work.
    // It only works in 'exec:', but not if there is a `script:` section:
    // task.workDir.resolve('.params.yml').text = yaml_str

    // Therefore, we inject it into the bash script:
    return """\
        cat <<"END_PARAMS_SECTION" > ./.params.yml
        ${indent_code_block(yaml_str, 8)}
        END_PARAMS_SECTION
    """
}

process RMARKDOWNNOTEBOOK {
    tag "$sample_id"
    label 'process_low'

    //NB: You likely want to override this with a container containing all required
    //dependencies for your analysis. The container at least needs to contain the
    //yaml and rmarkdown R packages.
    conda (params.enable_conda ? "r-base=4.1.0 r-rmarkdown=2.9 r-yaml=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-31ad840d814d356e5f98030a4ee308a16db64ec5:0e852a1e4063fdcbe3f254ac2c7469747a60e361-0' :
        'quay.io/biocontainers/mulled-v2-31ad840d814d356e5f98030a4ee308a16db64ec5:0e852a1e4063fdcbe3f254ac2c7469747a60e361-0' }"
    publishDir "${params.outdir}/${sample_id}", pattern: '{*.html}', mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), val(meta), path(input_files)

    output:
    tuple val(meta), path("*.html")           , emit: report
    tuple val(meta), path ("artifacts/*")     , emit: artifacts, optional: true
    tuple val(meta), path ("session_info.log"), emit: session_info
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample_id}"
    def parametrize = (task.ext.parametrize == null) ?  true : task.ext.parametrize
    def implicit_params = (task.ext.implicit_params == null) ? true : task.ext.implicit_params
    def meta_params = (task.ext.meta_params == null) ? true : task.ext.meta_params

    // Dump parameters to yaml file.
    // Using a yaml file over using the CLI params because
    //  * no issue with escaping
    //  * allows to pass nested maps instead of just single values
    def params_cmd = ""
    def render_cmd = ""
    if (parametrize) {
        nb_params = [:]
        if (implicit_params) {
            nb_params["cpus"] = task.cpus
            nb_params["artifact_dir"] = "artifacts"
            nb_params["input_dir"] = "./"
        }
        if (meta_params) {
            nb_params["meta"] = meta
        }
        nb_params += params
        params_cmd = dump_params_yml(nb_params)
        render_cmd = """\
            params = yaml::read_yaml('.params.yml')
            rmarkdown::render('${prefix}.Rmd', params=params, envir=new.env())
        """
    } else {
        render_cmd = "rmarkdown::render('${prefix}.Rmd')"
    }

    """
    # Dump .params.yml heredoc (section will be empty if parametrization is disabled)
    ${indent_code_block(params_cmd, 4)}
    
    echo

    # Create output directory
    mkdir artifacts

    # Set parallelism for BLAS/MKL etc. to avoid over-booking of resources
    export MKL_NUM_THREADS="$task.cpus"
    export OPENBLAS_NUM_THREADS="$task.cpus"
    export OMP_NUM_THREADS="$task.cpus"

    # Work around  https://github.com/rstudio/rmarkdown/issues/1508
    # If the symbolic link is not replaced by a physical file
    # output- and temporary files will be written to the original directory.
    cp -L "${projectDir}/${params.rmarkdown_template}" "${prefix}.Rmd"

    # Render notebook
    Rscript - <<EOF
        ${indent_code_block(render_cmd, 8)}
        writeLines(capture.output(sessionInfo()), "session_info.log")
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmarkdown: \$(Rscript -e "cat(paste(packageVersion('rmarkdown'), collapse='.'))")
    END_VERSIONS
    """
}

process EXPORT_PARAMETERS {

    publishDir "${params.tracedir}", pattern: '{*.json}', mode: 'copy', overwrite: true

    output:
    path 'parameters.json'

    script:
    "echo '${JsonOutput.toJson(params)}' > parameters.json"
}

process EXPORT_SAMPLEINFO {

    publishDir "${params.outdir}/${sample_id}", pattern: '{*.json}', mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), val(meta)
    
    output:
    path 'info.json'

    script:
    "echo '${JsonOutput.toJson(meta)}' > info.json"
}
