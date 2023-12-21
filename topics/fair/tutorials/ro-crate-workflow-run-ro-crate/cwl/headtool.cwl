class: CommandLineTool
cwlVersion: v1.0

baseCommand: head

inputs:
  lines:
    type: int
    inputBinding:
      position: 1
      prefix: "--lines"
  input:
    type: File
    inputBinding:
      position: 2
outputs:
  output:
    type: File
    outputBinding:
      glob: head_out.txt
stdout: head_out.txt
