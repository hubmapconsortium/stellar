#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: hubmap/stellar-prep-convert:latest
  InlineJavascriptRequirement: {}
baseCommand: "/opt/convert_input_ometiffs.py"

inputs:
  directory:
    type: Directory
    inputBinding:
      position: 0
  tissue
    type: string
    inputBinding:
      position: 1

outputs:
  h5ad_file:
    type: File
    outputBinding:
      glob: "*.h5ad"
  spatialdata_zarrs:
    type: Directory[]
    outputBinding:
      glob: "*_spatialdata.zarr"
