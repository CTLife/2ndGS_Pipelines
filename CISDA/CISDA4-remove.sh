ls   *.bam*   |  grep  -P "^removed_" | wc -l
rename s/removed_//  *
ls   *.bam*  |  grep -P  "^removed_" | wc -l

