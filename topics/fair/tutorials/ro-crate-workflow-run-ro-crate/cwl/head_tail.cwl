class: Workflow
cwlVersion: v1.0

inputs:
  input_file:
    type: File
  head_lines:
    type: int
  tail_lines:
    type: int
outputs:
  output_file:
    type: File
    outputSource: tail/output

steps:
  head:
    in:
      lines: head_lines
      input: input_file
    out: [output]
    run: headtool.cwl
  tail:
    in:
      lines: tail_lines
      input: head/output
    out: [output]
    run: tailtool.cwl
