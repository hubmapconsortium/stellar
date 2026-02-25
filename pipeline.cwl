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

steps:
  pre-convert:
    run: steps/pre-convert.cwl
    in:
      directory: data_dir
    out:
      - h5ad_file

  stellar:
    run: steps/stellar.cwl
    in:
      h5ad_file: pre-convert/h5ad_file
      tissue: tissue
    out:
      - stellar_results_for_sprm
