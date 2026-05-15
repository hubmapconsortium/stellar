class: CommandLineTool
cwlVersion: v1.1
baseCommand: ["python", "/stellar-main/STELLAR_run.py"]

requirements:
  DockerRequirement:
    dockerPull: hubmap/stellar:latest
  DockerGpuRequirement: {}
  EnvVarRequirement:
    envDef:
      CUDA_VISIBLE_DEVICES: "6"

inputs:
  h5ad_file:
    label: File containing cell data
    type: File
    inputBinding:
      position: 0
  tissue:
    label: tissue type string
    type: string
    inputBinding:
      position: 1

outputs:
  stellar_results_for_sprm:
    label: Directory containing outputs from STELLAR for SPRM
    type: Directory
    outputBinding:
      glob: stellar
