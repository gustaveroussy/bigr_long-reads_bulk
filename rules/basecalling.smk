"""
##########################################################################
These rules launch basecalling
##########################################################################
"""

"""
This rule launches basecalling
"""

rule dorado_basecalling:
    input:
        pod5_folder = os.path.normpath(OUTPUT_DIR + "/tmp/{samples_name}/{batch_name}/")
    output:
        os.path.normpath(OUTPUT_DIR + "/tmp/calling/{samples_name}/{batch_name}.bam")
    params:
        basic_model = DORADO_MODEL[0],
        math_model = DORADO_MODEL[1],
        model_path = os.path.dirname(DORADO_MODEL[0]),
        input_paths = ",".join(INPUT_PATHS)
    threads: 2
    resources:
        mem_mb=30720,
        partition="gpgpuq",
        gres="gpu:a100:1",
        time_min=10079
    shell:
        """
        if [ ! -z {params.math_model} ]; then echo "Methylation model provided ..."; meth="--modified-bases-models {params.math_model}"; else meth=""; fi
        echo -e "Launching dorado...\n"

        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        singularity exec --no-home --nv -B {OUTPUT_DIR},{params.model_path},{params.input_paths},/usr -B ${{TMP_DIR}}:/tmp \
        {SING_ENV_DORADO} dorado basecaller {params.basic_model} {input.pod5_folder} $meth -x cuda:0 --emit-moves --no-trim > {output}
        
        echo "Finished"
        """

