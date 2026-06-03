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
  provider:
    type: string?

outputs:
  stellar_results_for_sprm:
    type: Directory
    outputSource: stellar/stellar_results_for_sprm
  h5ad_file:
    type: File
    outputSource: pre-convert/h5ad_file
steps:
  pre-convert:
    run: steps/pre-convert.cwl
    in:
      directory: data_dir
      tissue: tissue
    out:
      - h5ad_file

  stellar:
    run: steps/stellar.cwl
    in:
      h5ad_file: pre-convert/h5ad_file
      tissue: tissue
    out:
      - stellar_results_for_sprm
