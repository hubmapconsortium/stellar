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
  ome_tiff_normalize:
    in:
      data_dir:
        source: data_dir
    out:
      - output_dir
    run: ome-tiff-normalize/ome_tiff_normalize.cwl

  pre-convert:
    run: steps/pre-convert.cwl
    in:
      directory: ome_tiff_normalize/output_dir
    out:
      - h5ad_file

  stellar:
    run: steps/stellar.cwl
    in:
      h5ad_file: pre-convert/h5ad_file
    out:
      - stellar_results_for_sprm
