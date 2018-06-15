
# Extract code blocks from an .Rmd file

import sys, textwrap

print("# This file is generated from the corresponding .Rmd file")
print()
print()

in_code = False
in_challenge = False
for line in sys.stdin:
    line = line.rstrip()
    if line.startswith("```"):
        print()
        in_code = not in_code
        assert in_code or line == "```", line
    elif in_code:
        print(line)
    elif line.startswith("#"):
        print("#" if in_challenge else "")
        in_challenge = "{.challenge}" in line
        if in_challenge:
            line = line.replace("{.challenge}","").rstrip()
        n = line.count("#")
        banner = "#"*n + " " + ("-" if n > 1 else "=") * (len(line)-n-1)
        print(banner)
        print(line)
        print(banner)
    elif in_challenge:
        for line2 in textwrap.wrap(line) or [""]:
            print("# " + line2)
