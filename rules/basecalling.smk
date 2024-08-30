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
        os.path.normpath(OUTPUT_DIR + "/calling/{samples_name}/{batch_name}.bam")
    params:
        basic_model = DORADO_MODEL[0],
        math_model = DORADO_MODEL[1]
    threads: 2
    resources:
        mem_mb=30720,
        partition="gpgpuq",
        gres="gpu:a100:1",
        time_min=10079
    shell:
        """
        #EXPORT CUDA LOCATION
        export LD_LIBRARY_PATH=/usr/local/cuda-11.1/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
        export PATH=/usr/local/cuda/bin:$PATH

        if [ ! -z {params.math_model} ]; then echo "Methylation model provided ..."; meth="--modified-bases-models {params.math_model}"; else meth=""; fi
        echo -e "Launching dorado...\n"
        {TOOL_DORADO} basecaller {params.basic_model} {input.pod5_folder} $meth -x cuda:0 --emit-moves --no-trim -b 1728 > {output}
        echo "Finished"
        """

"""
#CHANGING LOCATION TO PROJECT PATH AND CREATING SUBDIRECTORY
cd ${WKDIR}
mkdir -m 771 -p ${WKDIR}/calling/${DATASET_FOLDER}
"""