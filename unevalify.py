
import sys
for line in sys.stdin:
    if line.startswith("```{r"):
        line = "```{r eval=FALSE}\n"
    sys.stdout.write(line)
