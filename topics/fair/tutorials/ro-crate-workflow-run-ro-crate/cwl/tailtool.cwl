class: CommandLineTool
cwlVersion: v1.0

baseCommand: tail

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
      glob: tail_out.txt
stdout: tail_out.txt
