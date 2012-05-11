
for file in *.txt; do
    echo ";; $file"
    timeout 120 ../../basil --arrangement-pivot < $file | tee $file.out|  grep "\(basis orbits\|time\)" | sed -e 's/.*orbits:/#(/' -e 's/.*time://' -e 's/ms/)/' 
done | tee data.rkt
