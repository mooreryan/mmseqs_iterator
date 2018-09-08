BEGIN { OFS = "\t" }

{ print $1, $2, $11, annotation }
