# Rename all *..asc to *.asc in a SimDir
# Usage: ./clean_0D.sh <SimDir>
shopt -s globstar

for file in $1/**/*..asc; do
		echo "$file"
    mv -- "$file" "${file%..asc}.asc"
done

