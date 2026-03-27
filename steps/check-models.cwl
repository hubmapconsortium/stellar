#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: hubmap/stellar:latest
  InlineJavascriptRequirement: { }
baseCommand: "/opt/check_models.py"

inputs:
  tissue:
    type: string
    inputBinding:
      position: 0

  provider:
    type: string
    inputBinding:
      position: 1

  directory:
    type: Directory
    inputBinding:
      position: 2

outputs:
  model_test_results:
    label: whether or not a model exists for stellar annotations
    type: File
    outputBinding:
      glob: results.txt

