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
    type: string

outputs:
  stellar_results_for_sprm:
    type: Directory
    outputSource: stellar/stellar_results_for_sprm

steps:
  check-models:
    run: steps/check_models.cwl
    in:
      tissue: tissue
      provider: provider
      directory: data_dir
    out:
      - results

  pre-convert:
    run: steps/pre-convert.cwl
    in:
      directory: data_dir
      results: check-models/results
    out:
      - h5ad_file

  stellar:
    run: steps/stellar.cwl
    in:
      h5ad_file: pre-convert/h5ad_file
      results: check-models/results
    out:
      - stellar_results_for_sprm
