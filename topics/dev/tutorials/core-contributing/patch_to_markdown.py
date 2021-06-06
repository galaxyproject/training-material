import os
import re

        
DIFF_START = re.compile(r"diff --git a/(.*) b/(.*)")


def main():
    with open("user_extensions.patch", "r") as f:
        contents = f.read()

    diff_block = ""
    in_diff = None

    def handle_diff():
        print("A diff block starting with %s for file %s" % (diff_block[0:25], in_diff))
        diff_of_filename = in_diff.group(1)
        diff_of_basename = os.path.basename(diff_of_filename)
        with open(f"{diff_of_basename}_diff.md", "w") as f:
            indented_block = "\n".join([f"> {line}" for line in diff_block.splitlines()])
            f.write(f"""
> ### {{% icon solution %}} ``{diff_of_filename}``
> 
> Possible changes to file ``{diff_of_filename}``:
> 
> ```diff
{indented_block}
> ```
{{: .solution }}
""")

    for line in contents.splitlines():
        match = DIFF_START.match(line)
        if match:
            if in_diff:
                handle_diff()
            in_diff = match
            diff_block = ""
        else:
            diff_block += f"{line}\n"
    if in_diff:
        handle_diff()


if __name__ == "__main__":
    main()
