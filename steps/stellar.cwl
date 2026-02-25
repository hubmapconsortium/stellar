class: CommandLineTool
cwlVersion: v1.1
baseCommand: ["python", "/stellar-main/STELLAR_run.py"]

requirements:
  DockerRequirement:
    dockerPull: hubmap/stellar:0.4.outofband.sprm-convert
  DockerGpuRequirement: {}
  EnvVarRequirement:
    envDef:
      CUDA_VISIBLE_DEVICES: "0"

inputs:
  h5ad_file:
    label: File containing cell data
    type: File
    inputBinding:
      position: 0

outputs:
  stellar_results_for_sprm:
    label: Directory containing outputs from STELLAR for SPRM
    type: Directory
    outputBinding:
      glob: stellar
