
# Extract code blocks from an .Rmd file

import sys

print "# This file is generated from the corresponding .Rmd file"
print
print

in_code = False
in_challenge = False
for line in sys.stdin:
    line = line.rstrip()
    if line.startswith("```"):
        print
        in_code = not in_code
        assert in_code or line == "```", line
    elif in_code:
        print line
    elif line.startswith("#"):
        n = line.count("#")
        banner = "#"*n + " " + ("-" if n > 1 else "=") * (len(line)-n-1)
        print
        print banner
        print line
        print banner
        print
        in_challenge = ".challenge" in line
    elif in_challenge:
        print line
