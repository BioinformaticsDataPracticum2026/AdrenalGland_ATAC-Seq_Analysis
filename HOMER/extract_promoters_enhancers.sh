mkdir peak_types


for dir in results/*/; do
    subdir=$(basename "$dir")
    input="${dir}annotated_peaks.txt"
    output="peak_types/${subdir}.txt"

    awk -F'\t' 'BEGIN {OFS="\t"}
    NR==1 {
        id=1
        for (i=2; i<=NF; i++) {
            if ($i=="Annotation") annot=i
            if ($i=="Distance to TSS") dist=i
        }
        print "PeakID", "Annotation"
        next
    }
    {
        new_annot = ""
        if (((tolower($annot) ~ /intergenic/ || tolower($annot) ~ /intron/) && $dist >= 2000) || (tolower($annot) ~ /enhancer/)) {
            new_annot = "enhancer"
        } else if (tolower($annot) ~ /promoter/) {
            new_annot  = "promoter"
        }

        if (new_annot != ""){
            print $id, new_annot
        }
        
    }' "$input" > "$output"

done