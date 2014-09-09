#!/bin/bash    
#
# Contributed by Tobias Baumann
# Mainz, sep. 2013
#


OUTDOC=doc.txt

function getinfo {
for FILE in $(ls *.c)
do
    echo "%----------------------------------------------------------------------" >> $OUTDOC
    echo "\cprotect\subsection{File: \verb-$FILE -}" >> $OUTDOC

    START=($(grep -nr "#START_DOC" $FILE | awk -F: '{print $1}'))
    END=($(grep -nr "#END_DOC" $FILE | awk -F: '{print $1}'))

	for ((i = 0 ; i < ${#START[@]} ; i++ ))
	do	
		STARTi=$((${START[$i]} + 1))
		ENDi=$((${END[$i]} - 1))
		sed "$STARTi ,$ENDi !d" $FILE >> $OUTDOC
	done
done
}


# Generate latex code for /src
getinfo
mv doc.txt ../doc/Manual/ 

# Generate latex code for /src/fdstag/
cd fdstag
getinfo
mv doc.txt ../../doc/Manual/doc2.txt

# Create *.tex file
cd ../../doc/Manual/
cat doc2.txt >> doc.txt
sed '/%INSERT_DOC/ r doc.txt' code.tex_org > code.tex

# Invoke latex and bibtex
$PDFLATEX DevelDoc.tex
$BIBTEX DevelDoc.aux
$PDFLATEX DevelDoc.tex

# Clean up
rm doc2.txt doc.txt
