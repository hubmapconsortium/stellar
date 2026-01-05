#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.1
requirements:
  ScatterFeatureRequirement: {}

inputs:
  data_dir:
    type: Directory

outputs:
  stellar_results_for_sprm:
    type: Directory
    outputSource: stellar/stellar_results_for_sprm

steps:
  pre-convert:
    run: expr-mask-convert/pipeline.cwl
    in:
      directory: data_dir
    out:
      - h5ad_file
      - spatialdata_zarr

  stellar:
    run: steps/stellar.cwl
    in:
      h5ad_file: pre-convert/h5ad_file
    out:
      - stellar_results_for_sprm
