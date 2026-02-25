#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.1
requirements:
  ScatterFeatureRequirement: {}

inputs:
  data_dir:
    type: Directory

  tissue:
    type: string
    label: "tissue to be annotated"

outputs:
  stellar_results_for_sprm:
    type: Directory
    outputSource: stellar/stellar_results_for_sprm
  spatialdata_zarrs:
    type: Directory[]
    outputSource: pre-convert/spatialdata_zarrs
  h5ad_file:
    type: File
    outputSource: pre-convert/h5ad_file

steps:
  pre-convert:
    run: steps/pre-convert.cwl
    in:
      directory: data_dir
    out:
      - h5ad_file
      - spatialdata_zarrs

  stellar:
    run: steps/stellar.cwl
    in:
      h5ad_file: pre-convert/h5ad_file
      tissue: tissue
    out:
      - stellar_results_for_sprm
